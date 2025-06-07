import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.callbacks import EarlyStopping,ModelCheckpoint, ReduceLROnPlateau
# from sklearn.model_selection import train_test_split
import numpy as np
# import uproot
import pickle
import os
import argparse
# import awkward as ak
import random
import submitit
import json

import omnifold
from ReweightMCDataLoading import *

# set gpu growth
def set_gpu_growth():
  gpus = tf.config.list_physical_devices('GPU')
  if gpus:
    try:
      # Currently, memory growth needs to be the same across GPUs
      for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
        try:
          import horovod.tensorflow as hvd
          hvd.init()
          tf.config.experimental.set_visible_devices(gpus[hvd.local_rank()], 'GPU')
        except:
          print("No horovod")
      logical_gpus = tf.config.list_logical_devices('GPU')
      print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
    except RuntimeError as e:
      # Memory growth must be set before GPUs have been initialized
      print(e)

# fitting function
def doFit(model, train_dataset, val_dataset, epochs, verbose, callbacks, name):
    
    hist =  model.fit(
        train_dataset,
        epochs = epochs,
        validation_data = val_dataset,
        verbose = verbose,
        callbacks = callbacks
    )

    print(f"Last val loss {hist.history['val_loss'][0]}")
    print("INFO: Dumping training history ...")
    with open(name.replace(".weights.h5",".pkl"),"wb") as f:
        pickle.dump(hist.history, f)

    # weights = reweight(evaluate_data, model, batch_size=10000)
    # print(weights.shape)
    # np.save(name.replace(".weights.h5",".reweight.npy"), weights)    

    return hist

def doEvaluate(model, data, data_mean, data_std, name, batch_size=10000):

    # copy data
    temp = data.copy()

    # standardize data
    if data_mean is not None and data_std is not None:
        standardized_data = (temp - data_mean) / data_std
        temp = standardized_data

    # expect padded values to be np.nan
    temp = np.nan_to_num(temp, nan=0)
    # get weights
    weights = reweight(temp, model, batch_size=batch_size)
    print(weights.shape)
    np.save(name.replace(".weights.h5",".reweight.npy"), weights)

def train(conf):

    # SLURM sets CUDA_VISIBLE_DEVICES, so only the allocated GPU is visible to the task
    set_gpu_growth()
    gpu_id = os.environ.get('CUDA_VISIBLE_DEVICES')
    print(f"Assigned GPU: {gpu_id}")

    print(conf)

    # update %j with actual job number
    output_directory = conf["output_directory"]
    try:
        job_env = submitit.JobEnvironment()
        job_id = str(job_env.job_id)
    except:
        job_id = "%08x" % random.randrange(16**8)
        
    output_directory = os.path.abspath(output_directory.replace("%j", job_id))
    os.makedirs(output_directory, exist_ok=True)
    print(output_directory)

    # make weights directory
    weights_folder_id = "%08x" % random.randrange(16**8)
    weights_folder = os.path.abspath(os.path.join(output_directory, f"./model_weights_{weights_folder_id}"))
    os.makedirs(weights_folder, exist_ok=True)

    # Particle level distribution reweighting 
    aleph_mc = loadDataParticles(
        filePath = mc_paths["ArchivedMC"]["path"],
        treeName = mc_paths["ArchivedMC"]["tree"],
        branches = mc_paths["ArchivedMC"]["branches"],
        maxNPart = conf["maxNPart"],
    )
    aleph_mc = convert_PxPyPz_to_EtaPhiPmag(aleph_mc)
    print(aleph_mc.shape)
    new_mc = loadDataParticles(
        filePath = mc_paths[conf["new_mc_name"]]["path"], 
        treeName = mc_paths[conf["new_mc_name"]]["tree"],
        branches = mc_paths[conf["new_mc_name"]]["branches"],
        maxNPart = conf["maxNPart"],
    )
    new_mc = convert_PxPyPz_to_EtaPhiPmag(new_mc)
    print(new_mc.shape)

    # prepare datasets
    train_dataset_step1, val_dataset_step1, data_mean_step1, data_std_step1 = create_train_val_datasets(data_0=aleph_mc, data_1=aleph_mc, test_size=conf["test_size"], batch_size=conf["batch_size"], normalize=True, standardize=True)
    train_dataset_step2, val_dataset_step2, data_mean_step2, data_std_step2 = create_train_val_datasets(data_0=aleph_mc, data_1=new_mc, test_size=conf["test_size"], batch_size=conf["batch_size"], normalize=True, standardize=True)
    
    # save the standardization parameters
    np.savez(
        os.path.abspath(os.path.join(weights_folder, "dataset_standardization.npz")),
        data_mean_step1 = data_mean_step1,
        data_std_step1 = data_std_step1,
        data_mean_step2 = data_mean_step2,
        data_std_step2 = data_std_step2
    )

    # prepare PET network
    model = omnifold.PET(
        num_feat = aleph_mc.shape[2], 
        num_evt = 0, 
        num_part = conf["maxNPart"], 
        num_heads = 1, 
        num_transformer = 1, 
        local = False,
        projection_dim = 64
    )

    # training settings
    opt = tf.keras.optimizers.Adam(learning_rate=conf["lr"])
    model.compile(opt, loss = omnifold.net.weighted_binary_crossentropy)
    callbacks = [
        ReduceLROnPlateau(patience=1000, min_lr=1e-7, verbose=conf["verbose"], monitor="val_loss"),
        EarlyStopping(patience=conf["early_stopping_patience"], restore_best_weights=True,  monitor="val_loss"),
    ]

    # Step 1 pretrain the model to get it to unity
    model_name = os.path.join(weights_folder, f'Reweight_Step1.weights.h5')
    callbacks.append(ModelCheckpoint(model_name, save_best_only=True, mode='auto', save_weights_only=True))
    print("Running Step 1 (pre-train reweight aleph to aleph) with model name: ", model_name)
    hist = doFit(model, train_dataset_step1, val_dataset_step1, conf["step1_epochs"], conf["verbose"], callbacks, model_name)
    doEvaluate(model, aleph_mc, data_mean_step1, data_std_step1, model_name)

    # Step 2 train the model for reweighting
    model_name = os.path.join(weights_folder, f'Reweight_Step2.weights.h5')
    callbacks.append(ModelCheckpoint(model_name, save_best_only=True, mode='auto', save_weights_only=True))
    print("Running Step 2 (reweight aleph to new mc) with model name: ", model_name)
    hist = doFit(model, train_dataset_step2, val_dataset_step2, conf["step2_epochs"], conf["verbose"], callbacks, model_name)
    doEvaluate(model, aleph_mc, data_mean_step2, data_std_step2, model_name)


if __name__ == "__main__":
    
    # set up command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--slurm", help="path to json file containing slurm configuration", default=None)
    parser.add_argument("--njobs", help="number of jobs to actually launch. default is all", default=-1, type=int)
    args = parser.parse_args()

    # settings
    top_dir = "/pscratch/sd/b/badea/aleph/unfold-ee-logtau/ReweightMC/results/"
    top_dir = os.path.abspath(os.path.join(top_dir, f'training-{"%08x" % random.randrange(16**8)}', "%j"))
    conf = {
        "output_directory": top_dir,
        "test_size": 0.2,
        "lr": 5e-4,
        "batch_size": 2048,
        "verbose": True,
        "maxNPart": 80,
        "new_mc_name": "Sherpa",
        "step1_epochs" : 2,
        "step2_epochs" : 200,
        "early_stopping_patience" : 20,
    }
    confs = [conf]

    # if no slurm config file provided then just launch job
    if args.slurm == None:
        print("No slurm config file provided. Running jobs locally.")
        for iC, conf in enumerate(confs):
            # only launch a single job
            if args.njobs != -1 and (iC+1) > args.njobs:
                continue
            train(conf)
    else:
        # read in query
        query_path = os.path.abspath(args.slurm)
        with open(query_path) as f:
            query = json.load(f)

        # submission
        executor = submitit.AutoExecutor(folder=top_dir)
        executor.update_parameters(**query.get("slurm", {}))
      
        # loop over configurations
        jobs = []
        with executor.batch():
            for iC, conf in enumerate(confs):
                if args.njobs != -1 and (iC+1) > args.njobs:
                    continue
                job = executor.submit(train, conf) # **conf
                jobs.append(job)

    # # training settings
    # test_size = 0.2
    # lr = 5e-4
    # batch_size = 512
    # verbose = True
    # top_dir = "/pscratch/sd/b/badea/aleph/unfold-ee-logtau/ReweightMC/results/"
    # weights_folder = os.path.join(top_dir, f'training-{"%08x" % random.randrange(16**8)}') # "./"
    # os.makedirs(weights_folder, exist_ok=True)
    # new_mc_name = "PYTHIA8"
    # maxNPart = 80

    # # Particle level distribution reweighting 
    # aleph_mc = loadDataParticles(
    #     filePath = mc_paths["ArchivedPYTHIA6"]["path"],
    #     treeName = mc_paths["ArchivedPYTHIA6"]["tree"],
    #     branches = mc_paths["ArchivedPYTHIA6"]["branches"],
    #     maxNPart = maxNPart
    # )
    # print(aleph_mc.shape)
    # new_mc = loadDataParticles(
    #     filePath = mc_paths[new_mc_name]["path"], 
    #     treeName = mc_paths[new_mc_name]["tree"],
    #     branches = mc_paths[new_mc_name]["branches"],
    #     maxNPart = maxNPart
    # )
    # print(new_mc.shape)

    # # prepare datasets
    # train_dataset_step1, val_dataset_step1 = create_train_val_datasets(data_0=aleph_mc, data_1=aleph_mc, test_size=test_size, batch_size=batch_size, normalize=True)
    # train_dataset_step2, val_dataset_step2 = create_train_val_datasets(data_0=aleph_mc, data_1=new_mc, test_size=test_size, batch_size=batch_size, normalize=True)

    # # prepare PET network
    # model = omnifold.PET(
    #     num_feat = aleph_mc.shape[2], 
    #     num_evt = 0, 
    #     num_part = maxNPart, 
    #     num_heads = 1, 
    #     num_transformer = 1, 
    #     local = False,
    #     projection_dim = 64
    # )

    # # training settings
    # opt = tf.keras.optimizers.Adam(learning_rate=lr)
    # model.compile(opt, loss = omnifold.net.weighted_binary_crossentropy)
    # callbacks = [
    #     ReduceLROnPlateau(patience=1000, min_lr=1e-7, verbose=verbose, monitor="val_loss"),
    #     EarlyStopping(patience=10, restore_best_weights=True,  monitor="val_loss"),
    # ]

    # # Step 1 pretrain the model to get it to unity
    # model_name = os.path.join(weights_folder, f'Reweight_Step1.weights.h5')
    # callbacks.append(ModelCheckpoint(model_name, save_best_only=True, mode='auto', save_weights_only=True))
    # print("Running Step 1 (pre-train reweight aleph to aleph) with model name: ", model_name)
    # epochs = 2
    # hist = doFitAndEvaluate(model, train_dataset_step1, val_dataset_step1, epochs, verbose, callbacks, model_name, aleph_mc) 

    # # Step 2 train the model for reweighting
    # model_name = os.path.join(weights_folder, f'Reweight_Step2.weights.h5')
    # callbacks.append(ModelCheckpoint(model_name, save_best_only=True, mode='auto', save_weights_only=True))
    # print("Running Step 2 (reweight aleph to new mc) with model name: ", model_name)
    # epochs = 50
    # hist = doFitAndEvaluate(model, train_dataset_step2, val_dataset_step2, epochs, verbose, callbacks, model_name, aleph_mc)

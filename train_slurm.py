'''
Author: Anthony Badea
Date: May 27, 2024
'''

# fix for keras v3.0 update
import os
# os.environ['TF_USE_LEGACY_KERAS'] = '1' # THIS NEEDS TO BE COMMENTED OUT IF ON NERSC

# python based
import tensorflow as tf
import tensorflow.keras.backend as K
import random
from pathlib import Path
import time
import argparse
import json
import submitit
import h5py
import numpy as np
import shutil
import subprocess

# custom code
import dataloader

# omnifold
import omnifold

# SLURM sets CUDA_VISIBLE_DEVICES, so only the allocated GPU is visible to the task
# gpu_id = os.environ.get('CUDA_VISIBLE_DEVICES')
# print(f"Assigned GPU: {gpu_id}")

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
        
def train(
    conf
):

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
        
    output_directory = Path(str(output_directory).replace("%j", job_id))    
    os.makedirs(output_directory, exist_ok=True)
    print(output_directory)
        
    # load the aleph reconstructed data, reconstructed mc, and generator mc
    reco_data, reco_mc, gen_mc, pass_reco, pass_gen = dataloader.DataLoader(conf)

    # create the event weights
    weights_mc = np.ones(gen_mc.shape[0], dtype=np.float32)
    weights_data = np.ones(reco_data.shape[0], dtype=np.float32)
    
    # make omnifold dataloaders ready for training
    data = omnifold.DataLoader(
      reco = reco_data,
      weight = weights_data,
      normalize = True,
      bootstrap = True if conf["job_type"] == "BootstrapData" else False
    )
    
    mc = omnifold.DataLoader(
      reco = reco_mc,
      pass_reco = pass_reco,
      gen = gen_mc,
      pass_gen = pass_gen,
      weight = weights_mc,
      normalize = True,
      bootstrap = True if conf["job_type"] == "BootstrapMC" else False
    )

    # make weights directory
    # weights_folder = Path(output_directory, "./model_weights").resolve() "%08x" % random.randrange(16**8)
    weights_folder_id = "%08x" % random.randrange(16**8)
    weights_folder = Path(output_directory, f"./model_weights_{weights_folder_id}").resolve()
    weights_folder.mkdir()
    weights_folder = str(weights_folder)

    # save the starting weights
    outFileName = Path(weights_folder, f"starting_weights.h5").resolve()
    with h5py.File(outFileName, 'w') as hf:
      hf.create_dataset("weights_data", data=data.weight)
      hf.create_dataset("weights_mc", data=mc.weight)
      
    # write conf to json in output directory for logging
    output_conf_name = Path(weights_folder, "./conf.json").resolve()
    with open(output_conf_name, 'w') as file:
      json.dump(conf, file, indent=4)  # indent=4 for pretty printing
    
    # prepare networks
    ndim = reco_data.shape[1] # Number of features we are going to create = thrust
    model1 = omnifold.MLP(ndim)
    model2 = omnifold.MLP(ndim)
    
    # prepare multifold
    mfold = omnifold.MultiFold(
      name = 'mfold_job{}'.format(job_id),
      model_reco = model1,
      model_gen = model2,
      data = data,
      mc = mc,
      weights_folder = weights_folder,
      log_folder = output_directory,
      batch_size = conf["batch_size"],
      epochs = conf["epochs"],
      lr = conf["lr"],
      niter = conf["niter"],
      verbose = conf["verbose"],
      early_stop = conf["early_stop"],
    )
    
    # launch training
    mfold.Unfold()
    
    # get weights
    omnifold_weights = mfold.reweight(gen_mc[pass_gen], mfold.model2, batch_size=1000)
    
    # save weights to h5 file
    outFileName = Path(weights_folder, f"omnifold_weights.h5").resolve()
    with h5py.File(outFileName, 'w') as hf:
      hf.create_dataset("weights", data=omnifold_weights)
      
if __name__ == "__main__":

    # set up command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--slurm", help="path to json file containing slurm configuration", default=None)
    parser.add_argument("--njobs", help="number of jobs to actually launch. default is all", default=-1, type=int)
    parser.add_argument('--verbose', action='store_true', default=False,help='Run the scripts with more verbose output')
    parser.add_argument('--run_systematics', action='store_true', default=False,help='Run the track and event selection systematic variations')
    parser.add_argument('--run_bootstrap_mc', action='store_true', default=False,help='Run the bootstrapping for MC')
    parser.add_argument('--run_bootstrap_data', action='store_true', default=False,help='Run the bootstrapping for data')
    parser.add_argument('--run_ensembling', action='store_true', default=False,help='Run the ensembling by retraining without changing the inputs')
    args = parser.parse_args()
        
    # read in query
    if Path(args.slurm).resolve().exists():
        query_path = Path(args.slurm).resolve()
    else:
        # throw
        raise ValueError(f"Could not locate {args.slurm} in query directory or as absolute path")
    with open(query_path) as f:
        query = json.load(f)

    # create top level output directory
    top_dir = Path("results", f'./training-{"%08x" % random.randrange(16**8)}', "%j").resolve()

    # shared training configuration
    training_conf = {
      'output_directory' : str(top_dir),
      'FILE_MC':'/global/homes/b/badea/aleph/data/processed/20220514/alephMCRecoAfterCutPaths_1994_ThrustReprocess.npz',
      'FILE_DATA':'/global/homes/b/badea/aleph/data/processed/20220514/LEP1Data1994_recons_aftercut-MERGED_ThrustReprocess.npz',
      'TrackVariation': 1, # nominal track selection as written here https://github.com/badeaa3/ALEPHOmnifold/blob/main/src/Thrust.cxx#L143
      'EvtVariation': 0, # nominal event selection
      'niter': 3,
      'lr': 1e-4,
      'batch_size': 128,
      'epochs': 50,
      'early_stop': 10,
      'verbose' : args.verbose,
      'job_type' : None,
      'i_ensemble_per_omnifold': 0, # i-th ensemble per omnifold
    }

    # number of repeated trainings per omnifold configuration
    '''
    If we need N ensemble per omnifold. Then we need:
    40*N trainings for ensembeling uncertainty
    40*N trainings for data bootstrap uncertainty
    40*N trainings for mc bootstrap uncertainty
    18*N trainings for systematic variations
    = 138*N
    '''
    
    n_training_per_node = 4 # number of trainings per node or per job launched
    # n_ensemble_per_omnifold = n_training_per_node # this will be scaled by 4 from the running 4 copies on each node
    
    # list of configurations to launch
    confs = []

    # add configurations for track and event selection systematic variations
    total_n_systematics = 12 # closest to 10 which divides by 4
    n_systematics = int(total_n_systematics / n_training_per_node)
    if args.run_systematics:
      for iN in range(n_systematics):

        # track selections defined https://github.com/badeaa3/ALEPHOmnifold/blob/main/src/Thrust.cxx#L141-L150
        for TrackVariation in range(2, 9):
          temp = training_conf.copy()
          temp["TrackVariation"] = TrackVariation
          temp["job_type"] = "Systematics"
          temp["i_ensemble_per_omnifold"] = iN
          confs.append(temp)

        # event selection defined https://github.com/badeaa3/ALEPHOmnifold/blob/main/src/Thrust.cxx#L192-L197
        for EvtVariation in [1]:
          temp = training_conf.copy()
          temp["EvtVariation"] = EvtVariation
          temp["job_type"] = "Systematics"
          temp["i_ensemble_per_omnifold"] = iN
          confs.append(temp)

    # bootstrap mc
    total_n_bootstraps_mc = 40
    n_bootstraps_mc = int(total_n_bootstraps_mc / n_training_per_node)
    if args.run_bootstrap_mc:
      for strapn in range(n_bootstraps_mc):
        temp = training_conf.copy()
        temp["job_type"] = "BootstrapMC"
        temp["i_ensemble_per_omnifold"] = strapn
        confs.append(temp)

    # bootstrap data
    total_n_bootstraps_data = 40
    n_bootstraps_data = int(total_n_bootstraps_data / n_training_per_node)
    if args.run_bootstrap_data:
      for strapn in range(n_bootstraps_data):
        temp = training_conf.copy()
        temp["job_type"] = "BootstrapData"
        temp["i_ensemble_per_omnifold"] = strapn
        confs.append(temp)
    
    # add configurations for ensembling
    total_n_ensembles = 400
    n_ensembles = int(total_n_ensembles / n_training_per_node)
    if args.run_ensembling:
      for i in range(n_ensembles):
        temp = training_conf.copy()
        temp["job_type"] = "Ensembling"
        temp["i_ensemble_per_omnifold"] = iN
        confs.append(temp)
        
    # if submitit false then just launch job
    if not query.get("submitit", False):
        for iC, conf in enumerate(confs):
            # only launch a single job
            if args.njobs != -1 and (iC+1) > args.njobs:
                continue
            train(conf)
        exit()
    

    # submission
    executor = submitit.AutoExecutor(folder=top_dir)
    executor.update_parameters(**query.get("slurm", {}))
    # the following line tells the scheduler to only run at most 2 jobs at once. By default, this is several hundreds
    # executor.update_parameters(slurm_array_parallelism=2)
    
    # loop over configurations
    jobs = []
    with executor.batch():
        for iC, conf in enumerate(confs):
            
            # only launch a single job
            if args.njobs != -1 and (iC+1) > args.njobs:
                continue
            
            # print(conf)

            job = executor.submit(train, conf) # **conf
            jobs.append(job)

'''
Author: Anthony Badea
Date: May 27, 2024
'''

# fix for keras v3.0 update
import os
os.environ['TF_USE_LEGACY_KERAS'] = '1' 

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

# custom code
import dataloader

# omnifold
import omnifold

# set gpu growth
gpus = tf.config.list_physical_devices('GPU')
if gpus:
  try:
    # Currently, memory growth needs to be the same across GPUs
    for gpu in gpus:
      tf.config.experimental.set_memory_growth(gpu, True)
    logical_gpus = tf.config.list_logical_devices('GPU')
    print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
  except RuntimeError as e:
    # Memory growth must be set before GPUs have been initialized
    print(e)


def train(
    conf
):

    print(conf)
    
    # update %j with actual job number
    output_directory = conf["output_directory"]
    try:
        job_env = submitit.JobEnvironment()
        job_id = job_env.job_id
        output_directory = Path(str(output_directory).replace("%j", str(job_id)))
    except:
        job_id = random.randrange(16**8)
        output_directory = Path(str(output_directory).replace("%j", "%08x" % job_id))
        
    os.makedirs(output_directory, exist_ok=True)
    print(output_directory)
    
    # load the aleph reconstructed data, reconstructed mc, and generator mc
    reco_data, reco_mc, gen_mc, pass_reco, pass_gen = dataloader.DataLoader(conf)

    # create the event weights
    weights_mc = np.random.poisson(1, gen_mc.shape[0]) if conf["poisson_weights"] == "mc" else np.ones(gen_mc.shape[0], dtype=np.float32)
    weights_data = np.random.poisson(1, data.shape[0]) if conf["poisson_weights"] == "data" else np.ones(data.shape[0], dtype=np.float32)

    # make omnifold dataloaders ready for training
    data = omnifold.DataLoader(
      reco = reco_data,
      weight = weights_data,
      normalize = True
    )
    
    mc = omnifold.DataLoader(
      reco = reco_mc,
      pass_reco = pass_reco,
      gen = gen_mc,
      pass_gen = pass_gen,
      weight = weights_mc,
      normalize=True
    )
    
    # make weights directory
    weights_folder = Path(output_directory, "./model_weights").resolve()
    weights_folder.mkdir()
    weights_folder = str(weights_folder)

    for itrial in range(conf['ntrial']):

      # clear previous session
      K.clear_session()
      
      # prepare networks
      ndim = reco_data.shape[1] # Number of features we are going to create = thrust
      model1 = omnifold.MLP(ndim)
      model2 = omnifold.MLP(ndim)

      # prepare multifold
      mfold = omnifold.MultiFold(
        name = 'mfold_trial{}_strapn{}'.format(itrial, conf["strapn"]),
        model_reco = model1,
        model_gen = model2,
        data = data,
        mc = mc,
        batch_size = conf["batch_size"],
        epochs = conf["epochs"],
        lr = conf["lr"],
        # nstrap = conf["strapn"],
        niter = conf["niter"],
        weights_folder = weights_folder,
        verbose = conf["verbose"],
        early_stop = conf["early_stop"],
      )

      # launch training
      mfold.Preprocessing()
      mfold.Unfold()
    
      # get weights
      omnifold_weights = mfold.reweight(gen_mc[pass_gen], mfold.model2)
    
      # save weights to h5 file
      outFileName = Path(output_directory, f"omnifold_weights_trial{itrial}.h5").resolve()
      with h5py.File(outFileName, 'w') as hf:
        hf.create_dataset("weights", data=omnifold_weights)

      # move log file
      mfold.log_file.close()
      shutil(mfold.log_file.name, output_directory)
      
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
      'FILE_MC':'/home/badea/e+e-/aleph/data/processed/20220514/alephMCRecoAfterCutPaths_1994_ThrustReprocess.npz',
      'FILE_DATA':'/home/badea/e+e-/aleph/data/processed/20220514/LEP1Data1994_recons_aftercut-MERGED_ThrustReprocess.npz',
      'TrackVariation': 0, # nominal track selection
      'EvtVariation': 0, # nominal event selection
      'ntrial': 2, # number of times to run the unfolding per training
      'niter': 1,
      'lr': 1e-4,
      'batch_size': 128,
      'epochs': 1,
      'early_stop': 10,
      'strapn' : 0,
      'verbose' : args.verbose,
      'poisson_weights' : None
    }
    
    # list of configurations to launch
    confs = []

    # add configurations for track and event selection systematic variations
    if args.run_systematics:
      for TrackVariation in range(0, 9):
        for EvtVariation in range(0, 2):
          temp = training_conf.copy() # copy overall
          temp["TrackVariation"] = TrackVariation
          temp["EvtVariation"] = EvtVariation
          confs.append(temp)

    # add configurations for bootstrap mc or data
    if args.run_bootstrap_mc or args.run_bootstrap_data:

      # double check they aren't both true
      if args.run_bootstrap_mc == True and args.run_bootstrap_data == True:
        print("Warning you set bootstrap mc and data to True. Exiting to make sure that's what you want")
        exit()

      # number of bootstraps
      nstraps = 40
      for strapn in range(nstraps):
        temp = training_conf.copy() # copy overall
        temp["strapn"] = strapn
        temp["poisson_weights"] = "mc" if args.run_bootstrap_mc else "data"
        confs.append(temp)

    # add configurations for ensembling
    if args.run_ensembling:
      nensemble = 40
      for i in range(nensemble):
        temp = training_conf.copy()
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
            
            print(conf)

            job = executor.submit(train, conf) # **conf
            jobs.append(job)

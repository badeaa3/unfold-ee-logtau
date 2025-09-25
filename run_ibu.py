'''
Author: Anthony Badea
Date: May 27, 2024
'''

# fix for keras v3.0 update
import os
import random
import time
import argparse
import json
import submitit
import h5py
import numpy as np
import shutil
import subprocess
import math
import itertools
import uproot

# custom code
import ibu

def getThrustBins(conf):
    if conf["obs"] == "tau":
        return np.linspace(0, 0.42, 43) # tau bins
    if conf["obs"] == "logtau":
        return None # to be implemented
    return None #

def unfold(
    conf
):
    
    # print conf
    print(conf)

    # print random seed information
    print("SLURM_JOBID:", os.environ.get("SLURM_JOBID"))
    print("SLURM_ARRAY_TASK_ID:", os.environ.get("SLURM_ARRAY_TASK_ID"))
    # Check Python's built-in RNG
    print("Python random sample:", random.random())
    # Check NumPy RNG
    print("NumPy random sample:", np.random.rand())
    
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

    # write conf to json in output directory for logging
    output_conf_name = os.path.abspath(os.path.join(output_directory, "conf.json"))
    with open(output_conf_name, 'w') as file:
      json.dump(conf, file, indent=4)  # indent=4 for pretty printing
    
    # load data files
    with uproot.open(os.path.join(conf["storage"], conf["data"])) as f:
      data = np.array(f["t/Thrust"]) # data
      data_mask = np.array(f["t/passEventSelection"])
      print(data.shape, data_mask.sum())

    # load mc reco
    with uproot.open(os.path.join(conf["storage"], conf["reco"])) as f:
      mc_reco = np.array(f["t/Thrust"])
      mc_reco_mask = np.array(f["t/passEventSelection"])
      print(mc_reco.shape, mc_reco_mask.sum())

    # load mc gen
    with uproot.open(os.path.join(conf["storage"], conf["gen"])) as f:
      mc_gen = np.array(f["tgen/Thrust"])
      print(mc_gen.shape)

    # load mc genBefore
    with uproot.open(os.path.join(conf["storage"], conf["gen"].replace("tgen", "tgenBefore"))) as f:
      mc_genBefore = np.array(f["tgenBefore/Thrust"])
      print(mc_genBefore.shape)
      
    # apply observable, default loading is thrust
    if conf["obs"] == "tau":
        data = 1-data
        mc_reco = 1-mc_reco
        mc_gen = 1-mc_gen
        mc_genBefore = 1-mc_genBefore
    elif conf["obs"] == "logtau":
        data = np.log(1-data)
        mc_reco = np.log(1-mc_reco)
        mc_gen = np.log(1-mc_gen)
        mc_genBefore = np.log(1-mc_genBefore)
    else:
        print("Defaulting to Thrust observable and binning")

    # calculate bin width
    bins = getThrustBins(conf)
    binwidth = bins[1] - bins[0]
    binwidth_det = binwidth
    binwidth_mc = binwidth

    # run a closure test where data is replaced by reco MC
    if "run_closure_test" in conf.keys() and conf["run_closure_test"]:
       print("Running a closure test where data is replaced by reco MC")
       data = mc_reco
       
    print(f"Number of data events (selected) {data.shape} ({data_mask.sum()})")
    print(f"Number of MC reco events (selected) {mc_reco.shape} ({mc_reco_mask.sum()})")
    print(f"Number of MC gen events {mc_gen.shape}")

    # create the event weights
    weights_mc = np.ones(mc_gen.shape[0], dtype=np.float32)
    weights_data = np.ones(data.shape[0], dtype=np.float32)

    # apply theory reweighting
    if "theory_variation_weights_path" in conf.keys():
        print("Using theory variation weights")
        theory_variation_weights = np.load(conf["theory_variation_weights_path"])
        weights_mc = theory_variation_weights

    # get the histograms for selected events
    gen_hist = np.histogram(mc_gen[mc_reco_mask], bins=bins, density=True, weights=weights_mc[mc_reco_mask])[0]
    data_hist = np.histogram(data[data_mask], bins=bins, density=True, weights=weights_data[data_mask])[0]

    # compute (and normalize) the response matrix between GEN and SIM for selected events
    response = np.histogram2d(mc_reco[mc_reco_mask], mc_gen[mc_reco_mask], bins=(bins, bins), weights=weights_mc[mc_reco_mask])[0]
    response /= (response.sum(axis=0) + 10**-50)

    # perform iterative bayesian unfolding
    ibu_phis = ibu.ibu(data_hist, response, gen_hist, binwidth_det, binwidth_mc, it=conf["niter"])
    ibu_phi_unc = ibu.ibu_unc(data_hist, response, mc_gen[mc_reco_mask], binwidth_det, bins, binwidth_mc, it=5, nresamples=20) # note bins_mc = bins here
    # ibu_phi_unc = ibu.ibu_unc(ob, it=itnum, nrespamples=50) # udpate to take in the actual values, this is bootstrapping. This relies on reweighting. Can also use this for the theory reweighting
    
    np.save(os.path.abspath(os.path.join(output_directory, "ibu_phis.npy")), ibu_phis)
    np.save(os.path.abspath(os.path.join(output_directory, "ibu_phi_unc.npy")), ibu_phi_unc)

    # compute hadronic event selection
    if conf["job_type"] == "Nominal":
        gen_hist = np.histogram(mc_gen, bins=bins, density=True)[0]
        gen_bhist = np.histogram(mc_genBefore, bins=bins, density=True)[0]
        corrs = np.ones(gen_bhist.shape)
        corrs = gen_bhist/(gen_hist + 10**-50)
        np.save(os.path.abspath(os.path.join(output_directory, "hadronic_event_sel_corr.npy")), corrs)
    
if __name__ == "__main__":

    # set up command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--slurm", help="path to json file containing slurm configuration", default=None)
    parser.add_argument("--njobs", help="number of jobs to actually launch. default is all", default=-1, type=int)
    parser.add_argument('--verbose', action='store_true', default=False, help='Run the scripts with more verbose output')
    parser.add_argument('--run_systematics', action='store_true', default=False, help='Run the track and event selection systematic variations')
    parser.add_argument('--run_bootstrap_mc', action='store_true', default=False, help='Run the bootstrapping for MC')
    parser.add_argument('--run_bootstrap_data', action='store_true', default=False, help='Run the bootstrapping for data')
    # parser.add_argument('--run_ensembling', action='store_true', default=False, help='Run the ensembling by retraining without changing the inputs')
    parser.add_argument('--run_closure_test', action='store_true', default=False, help="Run a closure test where data is replaced by reco MC.")
    # parser.add_argument('--run_hyperparameter_scan', action='store_true', default=False, help='Run the hyperparameter scan')
    parser.add_argument('--run_niter_scan', action='store_true', default=False, help='Run the number of iteration scan based on the optimized hyperparameters')
    parser.add_argument('--run_theory_uncert', action='store_true', default=False, help='Run the theory uncertainty scan')
    parser.add_argument('--top_dir', help="Top level directory for storing data. Default to nersc directory", default="/pscratch/sd/b/badea/aleph/unfold-ee-logtau/UniFold/results/")
    args = parser.parse_args()

    # create top level output directory
    top_dir = args.top_dir
    top_dir = os.path.abspath(os.path.join(top_dir, f'training-{"%08x" % random.randrange(16**8)}', "%j"))
    
    # load configuration
    with open("training_conf.json") as f:
      training_conf = json.load(f)
    training_conf["output_directory"] = top_dir
    training_conf["verbose"] = args.verbose
    print(training_conf)
    # update gen to be tgen rather than tgenBefore
    training_conf["gen"] = training_conf["gen"].replace("tgenBefore", "tgen") # binned unfolding and then apply the hadronic event selection correction after
    training_conf["niter"] = 5 # number of IBU iterations
    training_conf["obs"] = "tau" # tau or log(tau)
    
    # configurations
    confs = []

    # nominal
    temp = training_conf.copy()
    temp["job_type"] = "Nominal"
    confs.append(temp)
    
    # sysematic variations
    if args.run_systematics:
        SystematicVariationList = ["ntpc7", "pt04", "ech10", "no_neutrals", "with_met"]
        for SystematicVariation in SystematicVariationList:
            temp = training_conf.copy()
            temp["data"] = temp["data"].replace("nominal", SystematicVariation)
            temp["reco"] = temp["reco"].replace("nominal", SystematicVariation)
            temp["job_type"] = "Systematics"
            temp["i_ensemble_per_omnifold"] = i
            confs.append(temp)

    # add configurations for theory uncertainty scan
    if args.run_theory_uncert:
      # theory_variation_dir = "/pscratch/sd/b/badea/aleph/unfold-ee-logtau/ReweightMC/results/training-200471c7/"
      theory_variation_dir = "/home/badea/e+e-/aleph/UnfoldThrustResults/theory_reweighting/training-200471c7/"
      theory_variations = [
        ["Pythia8", os.path.join(theory_variation_dir, "39912440_0/model_weights_b7634c53/Reweight_Step2.reweight.npy")],
        ["Herwig", os.path.join(theory_variation_dir, "39912440_1/model_weights_cc44b19d/Reweight_Step2.reweight.npy")],
        ["Sherpa", os.path.join(theory_variation_dir, "39912440_2/model_weights_afd3a072/Reweight_Step2.reweight.npy")]
      ]
      for name, inFileName in theory_variations:
          temp = training_conf.copy()
          temp["job_type"] = f"TheoryUncertainty_{name}"
          temp["i_ensemble_per_omnifold"] = i
          temp["theory_variation_weights_path"] = inFileName
          confs.append(temp)

    # bootstrap mc
    n_bootstraps_mc = 1
    if args.run_bootstrap_mc:
      for i in range(n_bootstraps_mc):
        temp = training_conf.copy()
        temp["job_type"] = "BootstrapMC"
        temp["i_ensemble_per_omnifold"] = i
        confs.append(temp)

    # bootstrap data
    n_bootstraps_data = 1
    if args.run_bootstrap_data:
      for i in range(n_bootstraps_data):
        temp = training_conf.copy()
        temp["job_type"] = "BootstrapData"
        temp["i_ensemble_per_omnifold"] = i
        confs.append(temp)
    
    # add configuration for closure check with a single job
    if args.run_closure_test:
      temp = training_conf.copy()
      temp["run_closure_test"] = args.run_closure_test
      temp["job_type"] = "ClosureTest"
      confs.append(temp)

    # add configurations for niter scan
    if args.run_niter_scan:
      niter = list(range(1,7))
      for niter in niter:
        temp = training_conf.copy()
        temp["job_type"] = "NiterScan"
        temp["niter"] = niter
        confs.append(temp)

    # if no slurm config file provided then just launch job
    if args.slurm == None:
      
      print("No slurm config file provided. Running jobs locally.")
      for iC, conf in enumerate(confs):
          # only launch a single job
          if args.njobs != -1 and (iC+1) > args.njobs:
              continue
          unfold(conf)
    
    # if slurm config file provided then launch job on slurm
    else:
      
      # read in query
      query_path = os.path.abspath(args.slurm)
      if not os.path.exists(query_path):
        raise ValueError(f"Could not locate {args.slurm}")
      with open(query_path) as f:
        query = json.load(f)

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

              job = executor.submit(unfold, conf) # **conf
              jobs.append(job)

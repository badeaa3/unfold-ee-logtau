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
      mc_gen_uniqueID = np.array(f["tgen/uniqueID"])
      print(mc_gen.shape)

    # load mc genBefore
    with uproot.open(os.path.join(conf["storage"], conf["gen"].replace("tgen", "tgenBefore"))) as f:
      mc_genBefore = np.array(f["tgenBefore/Thrust"])
      mc_genBefore_uniqueID = np.array(f["tgenBefore/uniqueID"])
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
    np.save(os.path.abspath(os.path.join(output_directory, "bins.npy")), bins) # save bins by default
    binwidth = bins[1] - bins[0]
    binwidth_det = binwidth
    binwidth_mc = binwidth

    # run a closure test where data is replaced by reco MC
    if "run_closure_test" in conf.keys() and conf["run_closure_test"]:
       print("Running a closure test where data is replaced by reco MC")
       data = mc_reco
       data_mask = mc_reco_mask
       
    print(f"Number of data events (selected) {data.shape} ({data_mask.sum()})")
    print(f"Number of MC reco events (selected) {mc_reco.shape} ({mc_reco_mask.sum()})")
    print(f"Number of MC gen events {mc_gen.shape}")

    # create the event weights
    weights_mc = np.ones(mc_gen.shape[0], dtype=np.float32)
    weights_data = np.ones(data.shape[0], dtype=np.float32)

    # apply theory reweighting
    if "theory_variation_weights_path" in conf.keys():
        print("Using theory variation weights")
        theory_variation_weights = np.load(conf["theory_variation_weights_path"]) # need to update this since the reweighting is applied to the genBefore events so the weights_mc[mc_reco_mask] crashes
        # weights at tgenBefore level, need to pick up the correct ones for tgen
        intersect, ind_mc_gen, ind_mc_genBefore = np.intersect1d(mc_gen_uniqueID, mc_genBefore_uniqueID, return_indices=True)
        temp_w = np.zeros(len(ind_mc_gen))
        temp_w[ind_mc_gen] = theory_variation_weights[ind_mc_genBefore]
        weights_mc = temp_w

    ibu_phis = ibu.ibu_wrapper(data[data_mask], weights_data[data_mask], mc_reco[mc_reco_mask], mc_gen[mc_reco_mask], weights_mc[mc_reco_mask], binwidth, bins, it=conf["niter"], output_directory=output_directory)
    np.save(os.path.abspath(os.path.join(output_directory, "ibu_phis.npy")), ibu_phis)
    print(ibu_phis[-1])

    # perform bootstrapping
    if conf["job_type"] == "Nominal":
        
        ibu_boostrap_data = ibu.ibu_bootstrap(data[data_mask], weights_data[data_mask], mc_reco[mc_reco_mask], mc_gen[mc_reco_mask], weights_mc[mc_reco_mask], binwidth, bins, boot="data", it=conf["niter"], nresamples=conf["nresamples"])
        np.save(os.path.abspath(os.path.join(output_directory, "ibu_phi_bootstrap_data.npy")), ibu_boostrap_data)
        
        ibu_boostrap_mc = ibu.ibu_bootstrap(data[data_mask], weights_data[data_mask], mc_reco[mc_reco_mask], mc_gen[mc_reco_mask], weights_mc[mc_reco_mask], binwidth, bins, boot="mc", it=conf["niter"], nresamples=conf["nresamples"])
        np.save(os.path.abspath(os.path.join(output_directory, "ibu_phi_bootstrap_mc.npy")), ibu_boostrap_mc)
        
    # compute hadronic event selection
    # if conf["job_type"] == "Nominal":
    gen_hist = np.histogram(mc_gen[mc_reco_mask], bins=bins, density=True)[0]
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
    parser.add_argument('--run_nominal', action='store_true', default=False, help='Run the nominal unfolding')
    parser.add_argument('--run_systematics', action='store_true', default=False, help='Run the track and event selection systematic variations')
    # parser.add_argument('--run_bootstrap_mc', action='store_true', default=False, help='Run the bootstrapping for MC')
    # parser.add_argument('--run_bootstrap_data', action='store_true', default=False, help='Run the bootstrapping for data')
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
    training_conf["niter"] = 4 # number of IBU iterations
    training_conf["obs"] = "tau" # tau or log(tau)
    training_conf["storage"] = "/home/badea/e+e-/aleph/unfold-ee-logtau/DataProcessing/100725/7" # update the storage directory
    
    # configurations
    confs = []

    # nominal
    if args.run_nominal:
        temp = training_conf.copy()
        temp["job_type"] = "Nominal"
        temp["nresamples"] = 100
        confs.append(temp)
    
    # sysematic variations
    if args.run_systematics:

        SystematicVariationList = ["ntpc7", "pt04", "ech10", "no_neutrals", "with_met"] # all cut based systematics
        NeutralParticleMCVariations = ["nes_up", "nes_down", "ner"] # reco MC variations to account for mis-modeling of detector response/efficiency for neutral particles
        SystematicVariationList += NeutralParticleMCVariations
        
        for SystematicVariation in SystematicVariationList:
            temp = training_conf.copy()
            # only apply cut based variations to data
            if SystematicVariation not in NeutralParticleMCVariations:
                temp["data"] = temp["data"].replace("nominal", SystematicVariation)
            # apply all variations to reco mc
            temp["reco"] = temp["reco"].replace("nominal", SystematicVariation)
            temp["job_type"] = "Systematics"
            confs.append(temp)

    # add configurations for theory uncertainty scan
    if args.run_theory_uncert:

        # # boost results
        # # theory_variation_dir = "/pscratch/sd/b/badea/aleph/unfold-ee-logtau/ReweightMC/results/training-200471c7/"
        # theory_variation_dir = "/home/badea/e+e-/aleph/UnfoldThrustResults/theory_reweighting/training-200471c7/"
        # theory_variations = [
        #     ["Pythia8", os.path.join(theory_variation_dir, "39912440_0/model_weights_b7634c53/Reweight_Step2.reweight.npy")],
        #     ["Herwig", os.path.join(theory_variation_dir, "39912440_1/model_weights_cc44b19d/Reweight_Step2.reweight.npy")],
        #     ["Sherpa", os.path.join(theory_variation_dir, "39912440_2/model_weights_afd3a072/Reweight_Step2.reweight.npy")]
        # ]

        # ensembled 15 trainings on nersc
        theory_variation_dir = "/home/badea/e+e-/aleph/UnfoldThrustResults/theory_reweighting/training-bf3b5fc3/"
        theory_variations = [
            ["Pythia8", os.path.join(theory_variation_dir, "Reweight_Step2_Ensemble_Pythia8.npy")],
            ["Herwig", os.path.join(theory_variation_dir, "Reweight_Step2_Ensemble_Herwig.npy")],
            ["Sherpa", os.path.join(theory_variation_dir, "Reweight_Step2_Ensemble_Sherpa.npy")],
        ]
      
        for name, inFileName in theory_variations:
            temp = training_conf.copy()
            temp["job_type"] = f"TheoryUncertainty_{name}"
            temp["theory_variation_weights_path"] = inFileName
            confs.append(temp)

    # # bootstrap mc
    # n_bootstraps_mc = 1
    # if args.run_bootstrap_mc:
    #   for i in range(n_bootstraps_mc):
    #     temp = training_conf.copy()
    #     temp["job_type"] = "BootstrapMC"
    #     temp["i_ensemble_per_omnifold"] = i
    #     confs.append(temp)

    # # bootstrap data
    # n_bootstraps_data = 1
    # if args.run_bootstrap_data:
    #   for i in range(n_bootstraps_data):
    #     temp = training_conf.copy()
    #     temp["job_type"] = "BootstrapData"
    #     temp["i_ensemble_per_omnifold"] = i
    #     confs.append(temp)
    
    # add configuration for closure check with a single job
    if args.run_closure_test:
      temp = training_conf.copy()
      temp["run_closure_test"] = args.run_closure_test
      temp["job_type"] = "ClosureTest"
      confs.append(temp)

    # add configurations for niter scan
    if args.run_niter_scan:
        niter = list(range(1,5))
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

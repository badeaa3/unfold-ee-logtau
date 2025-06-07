import uproot
import numpy as np
import os


def DataLoader(config):

    # load data files
    with uproot.open(os.path.join(config["storage"], config["data"])) as f:
      data = np.array(f["t/Thrust"]) # data
      data_mask = np.array(f["t/passEventSelection"])
      print(data.shape, data_mask.sum())

    # load mc reco
    with uproot.open(os.path.join(config["storage"], config["reco"])) as f:
      mc_reco = np.array(f["t/Thrust"])
      mc_reco_mask = np.array(f["t/passEventSelection"])
      print(mc_reco.shape, mc_reco_mask.sum())

    # load mc gen
    with uproot.open(os.path.join(config["storage"], config["gen"])) as f:
      mc_gen = np.array(f["tgenBefore/Thrust"])
      mc_gen_mask = np.ones(f["tgenBefore"].num_entries) # no selection on gen level

      # pick up the correct event ID's to use directly tgenbefore
      with uproot.open(os.path.join(config["storage"], config["gen"].replace("tgenBefore", "tgen"))) as f_tgen:
        a = np.array(f_tgen["tgen/uniqueID"]) # generator level with hadronic event selection, events matched to reco
        b = np.array(f["tgenBefore/uniqueID"]) # generator level without hadronic event selection
        intersect, ind_a, ind_b = np.intersect1d(a, b, return_indices=True)

        # make temporary reco lists with same size as tgenBefore
        pass_reco = np.zeros(len(b))
        reco_vals = np.nan*np.ones(len(b))
        # now populate the list with the reco values according to the uniqueID order of tgenBefore
        pass_reco[ind_b] = mc_reco_mask[ind_a]
        reco_vals[ind_b] = mc_reco[ind_a]
        # now update reco values
        mc_reco = reco_vals
        mc_reco_mask = pass_reco

    # convert to log
    data = np.log(1-data)
    mc_reco = np.log(1-mc_reco)
    mc_gen = np.log(1-mc_gen)

    # reco could get a nan from populating above, so we need to replace it with -999
    mc_reco[np.isnan(mc_reco)] = -999.

    # final preparation
    data = np.expand_dims(data[data_mask],-1) # We only want data events passing a selection criteria
    mc_reco = np.expand_dims(mc_reco, -1) # expand dimension
    mc_reco_mask = mc_reco_mask == 1 # convert to boolean mask
    mc_gen = np.expand_dims(mc_gen, -1) # expand dimension
    mc_gen_mask = mc_gen_mask == 1 # convert to boolean mask
    
    return data, mc_reco, mc_gen, mc_reco_mask, mc_gen_mask

if __name__ == "__main__":

  training_conf = {
     'storage' : '/global/cfs/cdirs/m3246/aleph/derivations/20250527',
     'data' : 'LEP1Data1994_recons_aftercut-MERGED_thrust_nominal_t.root',
     'reco' : 'alephMCRecoAfterCutPaths_1994_thrust_nominal_t.root',
     'gen' : 'alephMCRecoAfterCutPaths_1994_thrust_no_event_sel_tgenBefore.root',
  }
  data, mc_reco, mc_gen, mc_reco_mask, mc_gen_mask = DataLoader(training_conf)
  print(data.shape, mc_reco.shape, mc_reco_mask.shape, mc_gen.shape, mc_gen_mask.shape)
  # print(data[:5], mc_reco[:5], mc_gen[:5], mc_reco_mask[:5], mc_gen_mask[:5])
  print("data: ", data[:10])
  print("reco: ", mc_reco[:10])
  print("reco mask: ", mc_reco_mask[:10])
  print("gen: ", mc_gen[:10])
  print("gen mask: ", mc_gen_mask[:10])
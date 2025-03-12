import uproot
import numpy as np

def DataLoader(config):

    # load files
    f_data = uproot.open(config["f_data"])
    f_mc = uproot.open(config["f_mc"])

    # pick up the correct entries
    i = config["SystematicVariation"]
    # data systematic variation
    data = np.array([x[i] for x in np.array(f_data["t/Thrust"])]) # data
    data_mask = np.array([x[i] for x in np.array(f_data["t/passEventSelection"])])
    # mc reco systematic variation
    mc_reco = np.array([x[i] for x in np.array(f_mc["t/Thrust"])]) # reco level
    mc_reco_mask = np.array([x[i] for x in np.array(f_mc["t/passEventSelection"])])
    # generator level only one element
    mc_gen = np.array([x[0] for x in np.array(f_mc["tgenBefore/Thrust"])]) # generator level without hadronic event selection
    # mc_gen_mask = np.array([x[0] for x in np.array(f_mc["tgenBefore/passEventSelection"])])
    mc_gen_mask = np.ones(f_mc['tgenBefore'].num_entries) # no selection on gen level
    
    # convert to log
    data = np.log(1-data)
    mc_reco = np.log(1-mc_reco)
    mc_gen = np.log(1-mc_gen)
    
    ###### pick up the correct event ID's to use directly tgenbefore
    a = np.array(f_mc["tgen/uniqueID"]) # generator level with hadronic event selection, events matched to reco
    b = np.array(f_mc["tgenBefore/uniqueID"]) # generator level without hadronic event selection
    intersect, ind_a, ind_b = np.intersect1d(a, b, return_indices=True)

    # make temporary reco lists with same size as tgenBefore
    pass_reco = np.zeros(len(b))
    reco_vals = -999.*np.ones(len(b))
    # now populate the list with the reco values according to the uniqueID order of tgenBefore
    pass_reco[ind_b] = mc_reco_mask[ind_a]
    reco_vals[ind_b] = mc_reco[ind_a]
    # now update reco values
    mc_reco = reco_vals
    mc_reco_mask = pass_reco

    # print(mc_reco.shape, mc_reco_mask.shape)
    # final preparation
    data = np.expand_dims(data[data_mask],-1) # We only want data events passing a selection criteria
    mc_reco = np.expand_dims(mc_reco, -1) # expand dimension
    mc_reco_mask = mc_reco_mask == 1 # convert to boolean mask
    mc_gen = np.expand_dims(mc_gen, -1) # expand dimension
    mc_gen_mask = mc_gen_mask == 1 # convert to boolean mask
    
    return data, mc_reco, mc_gen, mc_reco_mask, mc_gen_mask

if __name__ == "__main__":

  training_conf = {
    'f_data':'/global/homes/b/badea/aleph/data/ThrustDerivation/030725/LEP1Data1994_recons_aftercut-MERGED_thrust.root',
    'f_mc':'/global/homes/b/badea/aleph/data/ThrustDerivation/030725/alephMCRecoAfterCutPaths_1994_thrust.root',
    'SystematicVariation': 1, # chosen selection
  }
  data, mc_reco, mc_gen, mc_reco_mask, mc_gen_mask = DataLoader(training_conf)
  print(data.shape, mc_reco.shape, mc_reco_mask.shape, mc_gen.shape, mc_gen_mask.shape)
  print(data[:5], mc_reco[:5], mc_gen[:5], mc_reco_mask[:5], mc_gen_mask[:5])
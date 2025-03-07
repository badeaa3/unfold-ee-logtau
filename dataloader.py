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
    mc_gen_mask = np.array([x[0] for x in np.array(f_mc["tgenBefore/passEventSelection"])])

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

    # final preparation
    data = np.expand_dims(data[data_mask],-1) # We only want data events passing a selection criteria
    mc_reco = np.expand_dims(mc_reco, -1) # expand dimension
    mc_reco_mask = mc_reco_mask == 1 # convert to boolean mask
    mc_gen = np.expand_dims(mc_gen, -1) # expand dimension
    mc_gen_mask = mc_gen_mask == 1 # convert to boolean mask
    
    return data, mc_reco, mc_gen, mc_reco_mask, mc_gen_mask
    
# def DataLoaderOld(config):

#     # load everything
#     DATA = np.load(config['FILE_DATA'], allow_pickle=True)
#     MC = np.load(config['FILE_MC'], allow_pickle=True)
    
#     # pick out the correct entries
#     data = DATA["t_thrust"][:,config['TrackVariation']]
#     data_mask = DATA[f"t_passEventSelection_{config['EvtVariation']}"][:,config['TrackVariation']]
#     mc_reco = MC["t_thrust"][:,config['TrackVariation']]
#     mc_reco_mask = MC[f"t_passEventSelection_{config['EvtVariation']}"][:,config['TrackVariation']]
#     mc_gen = np.stack(MC["tgenBefore_thrust"]).flatten() # little hack because format is array([array([]), array([]), ....]) and we need array([,,,,])
#     mc_gen_mask = np.stack(MC[f"tgenBefore_passEventSelection"]).flatten()
#     print(mc_reco.shape, mc_gen.shape)

#     ###### pick up the correct event ID's to use directly tgenbefore
#     a = MC['tgen_uniqueID'] # WHES_id
#     b = MC['tgenBefore_uniqueID'] # WHOES_id
#     intersect, ind_a, ind_b = np.intersect1d(a, b, return_indices=True)
#     print(ind_a.shape, ind_b.shape,  a.shape, b.shape)

#     pass_reco = np.zeros(len(b))
#     pass_reco[ind_b] = mc_reco_mask[ind_a]
    
#     reco_vals = -999.*np.ones(len(b))
#     reco_vals[ind_b] = mc_reco[ind_a]

#     mc_reco = reco_vals
#     mc_reco_mask = pass_reco
#     # print(mc_reco.shape, mc_gen.shape)
#     ########
    
#     # append to reco to make the same length as tgenbefore
#     diff = mc_gen.shape[0] - mc_reco.shape[0]
#     # print(diff)
#     mc_reco = np.concatenate([mc_reco, np.ones(diff)])
#     mc_reco_mask = np.concatenate([mc_reco_mask, np.ones(diff).astype(bool)])
    
#     #We only want data events passing a selection criteria
#     data = np.expand_dims(data[data_mask],-1)
#     mc_reco = np.expand_dims(mc_reco, -1)
#     mc_reco_mask = mc_reco_mask == 1
#     mc_gen = np.expand_dims(mc_gen,-1)
#     mc_gen_mask = mc_gen_mask == 1
    
#     return data, mc_reco, mc_gen, mc_reco_mask, mc_gen_mask

if __name__ == "__main__":

    # training_conf = {
    #   'FILE_MC':'/global/homes/b/badea/aleph/data/processed/20220514/alephMCRecoAfterCutPaths_1994_ThrustReprocess.npz',
    #   'FILE_DATA':'/global/homes/b/badea/aleph/data/processed/20220514/LEP1Data1994_recons_aftercut-MERGED_ThrustReprocess.npz',
    #   'TrackVariation': 1, # nominal track selection as written here https://github.com/badeaa3/ALEPHOmnifold/blob/main/src/Thrust.cxx#L143
    #   'EvtVariation': 0, # nominal event selection
    # }
    # data, mc_reco, mc_gen, mc_reco_mask, mc_gen_mask = DataLoaderOld(training_conf)
    # print(data.shape, mc_reco.shape, mc_reco_mask.shape, mc_gen.shape, mc_gen_mask.shape)
    # print(data[:5], mc_reco[:5], mc_gen[:5], mc_reco_mask[:5], mc_gen_mask[:5])
    # print("new")

    training_conf = {
      'f_data':'/global/homes/b/badea/aleph/data/ThrustDerivation/030725/LEP1Data1994_recons_aftercut-MERGED_thrust.root',
      'f_mc':'/global/homes/b/badea/aleph/data/ThrustDerivation/030725/alephMCRecoAfterCutPaths_1994_thrust.root',
      'SystematicVariation': 1, # chosen selection
    }
    data, mc_reco, mc_gen, mc_reco_mask, mc_gen_mask = DataLoader(training_conf)
    print(data.shape, mc_reco.shape, mc_reco_mask.shape, mc_gen.shape, mc_gen_mask.shape)
    print(data[:5], mc_reco[:5], mc_gen[:5], mc_reco_mask[:5], mc_gen_mask[:5])
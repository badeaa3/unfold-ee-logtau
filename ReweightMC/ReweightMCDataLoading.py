import uproot
import numpy as np
import awkward as ak

# file paths
mc_paths = {
    "ArchivedPYTHIA6" : {
        "path" : "/pscratch/sd/b/badea/aleph/data/alephMCRecoAfterCutPaths_1994.root",
        "tree" : "tgenBefore" # t=reco, tgen = generator level after hadronic event selection, tgenBefore = generator level before hadronic event selection
    },
    "HERWIG7" : {
        "path" : "/pscratch/sd/b/badea/aleph/data/LEP1MCVariations/abaty/HERWIG7/2_10_2024_LEP1MC/LEP-Matchbox-S1000-1_0_0.root",
        "tree" : "t"
    },
    "SHERPA" : {
        "path" : "/pscratch/sd/b/badea/aleph/data/LEP1MCVariations/abaty/SHERPA/2_10_2024_LEP1MC/Sherpa_RNG100_0_0.root",
        "tree" : "t"
    },
    "PYTHIA8" : {
        "path" : "/pscratch/sd/b/badea/aleph/data/LEP1MCVariations/hannah/LEP1_pythia8_MC_withSphericity_v2.root",
        "tree" : "tgen"
    },
    "PYTHIA8_DIRE" : {
        "path": "/pscratch/sd/b/badea/aleph/data/LEP1MCVariations/hannah/LEP1_pythia8_MC_DIRE.root",
        "tree" : "tgen"
    },
    "PYTHIA8_VINCIA" : {
        "path" : "/pscratch/sd/b/badea/aleph/data/LEP1MCVariations/hannah/LEP1_pythia8_MC_VINCIA.root",
        "tree" : "tgen"
    }
}

def expit(x):
    return 1. / (1. + np.exp(-x))

def reweight(data, model, batch_size, verbose=True):
    f = expit(model.predict(data,batch_size=batch_size,verbose=verbose))
    weights = f / (1. - f)  # this is the crux of the reweight, approximates likelihood ratio
    weights = np.nan_to_num(weights[:,0],posinf=1)
    return weights

def loadBranchAndPad(branch, maxN, value=0):
    a = branch.array()
    a = ak.to_numpy(ak.fill_none(ak.pad_none(a, max(maxN,np.max(ak.num(a)))),value)) # must take maximum so that it pads to a uniform max and then cut down
    a = a[:,:maxN]
    return a

def loadDataParticles(filePath, treeName, branches, maxNPart, padValue = 0):
    # kinematics
    data = []
    with uproot.open(filePath) as rFile:
        tree = rFile[treeName]
        # print(tree.num_entries)
        # print(tree.keys())
        for key in branches:
            if key in tree.keys():
                temp = loadBranchAndPad(tree[key], maxNPart, value=padValue) # pad to get to maxNPart
                data.append(temp)
            else:
                print(f"Attempted to add key not in tree: {key} for file {filePath}")
    # stack and return
    data = np.stack(data, axis=-1)
    return data
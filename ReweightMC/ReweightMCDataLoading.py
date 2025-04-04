import uproot
import numpy as np
import awkward as ak
from sklearn.model_selection import train_test_split
import tensorflow as tf

# file paths
mc_paths = {
    "ArchivedPYTHIA6" : {
        "path" : "/pscratch/sd/b/badea/aleph/data/alephMCRecoAfterCutPaths_1994.root",
        "thrust_path" : "/global/homes/b/badea/aleph/data/ThrustDerivation/030725/alephMCRecoAfterCutPaths_1994_thrust.root",
        "tree" : "tgenBefore", # t=reco, tgen = generator level after hadronic event selection, tgenBefore = generator level before hadronic event selection
        "branches" : ["pmag", "eta", "phi"],
    },
    "HERWIG7" : {
        "path" : "/pscratch/sd/b/badea/aleph/data/LEP1MCVariations/abaty/HERWIG7/2_10_2024_LEP1MC/LEP-Matchbox-S1000-1_0_0.root",
        "tree" : "t",
        "branches" : ["pmag", "eta", "phi"],
    },
    "SHERPA" : {
        "path" : "/pscratch/sd/b/badea/aleph/data/LEP1MCVariations/abaty/SHERPA/2_10_2024_LEP1MC/Sherpa_RNG100_0_0.root",
        "tree" : "t",
        "branches" : ["pmag", "eta", "phi"],
    },
    "PYTHIA8" : {
        "path" : "/pscratch/sd/b/badea/aleph/data/LEP1MCVariations/hannah/LEP1_PYTHIA8_MC_TGenBefore.root",
        "tree" : "tgenBefore",
        "branches" : ["pmag", "eta", "phi"],
    },
    "PYTHIA8_DIRE" : {
        "path": "/pscratch/sd/b/badea/aleph/data/LEP1MCVariations/hannah/LEP1_pythia8_MC_DIRE.root",
        "tree" : "tgen",
        "branches" : ["p", "eta", "phi"],
    },
    "PYTHIA8_VINCIA" : {
        "path" : "/pscratch/sd/b/badea/aleph/data/LEP1MCVariations/hannah/LEP1_pythia8_MC_VINCIA.root",
        "tree" : "tgen",
        "branches" : ["p", "eta", "phi"],
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

def loadData(filePath, treeName, branches, SystematicVariation=0):
    i = SystematicVariation
    data = []
    with uproot.open(filePath) as rFile:
        for branch in branches:
            temp = np.array([x[i] for x in np.array(rFile[f"{treeName}/{branch}"])]) # reco level
            data.append(temp)
    # stack to form data
    data = np.stack(data,axis=1)
    return data

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

def create_train_val_datasets(data_0, data_1, test_size=0.2, batch_size=512, normalize=True):
    
    # create labels
    labels_data_0 = np.zeros(len(data_0),dtype=np.float32)
    labels_data_1 = np.ones(len(data_1),dtype=np.float32)

    # create weights
    weights_data_0 = np.ones(len(data_0),dtype=np.float32)
    weights_data_1 = np.ones(len(data_1),dtype=np.float32)  

    # normalize
    if normalize:
        weights_data_1 *= (np.sum(weights_data_0) / np.sum(weights_data_1))
        print(weights_data_0.sum(), weights_data_1.sum())

    # concatenate
    data = np.concatenate((data_0, data_1))
    # labels expected to be stack of [labels, weights]
    labels = np.concatenate((labels_data_0, labels_data_1))
    weights = np.concatenate((weights_data_0, weights_data_1))
    labels = np.stack([labels, weights], axis=1).astype(np.float32)
    # print(data.shape, labels.shape)

    # standardize data
    # data = (data - np.mean(data, axis=0)) / np.std(data, axis=0)

    # Split into training (80%) and validation (20%) sets randomly
    train_data, val_data, train_labels, val_labels = train_test_split(data, labels, test_size=test_size, stratify=labels)
    # print(train_data.shape, val_data.shape, train_labels.shape, val_labels.shape)
    
    # Create TensorFlow datasets
    train_dataset = tf.data.Dataset.from_tensor_slices((train_data, train_labels))
    val_dataset = tf.data.Dataset.from_tensor_slices((val_data, val_labels))

    # Batch and shuffle training data
    train_dataset = train_dataset.shuffle(buffer_size=len(train_data)).batch(batch_size).prefetch(tf.data.AUTOTUNE)
    val_dataset = val_dataset.batch(batch_size).prefetch(tf.data.AUTOTUNE)

    return train_dataset, val_dataset
import uproot
import numpy as np
import awkward as ak
from sklearn.model_selection import train_test_split
import tensorflow as tf

# file paths
mc_paths = {
        "ArchivedMC" : {"path": "/global/cfs/cdirs/m3246/aleph/derivations/20250527/alephMCRecoAfterCutPaths_1994_thrust_no_event_sel_tgenBefore.root", "tree" : "tgenBefore", "branches" : ["px", "py", "pz"]},
        "Pythia8" : {"path": "/global/cfs/cdirs/m3246/aleph/derivations/20250527/LEP1_PYTHIA8_MC_TGenBefore_NoISR_thrust_no_event_sel_tgenBefore.root", "tree" : "tgenBefore", "branches" : ["px", "py", "pz"]},
        "Herwig" : {"path": "/global/cfs/cdirs/m3246/aleph/derivations/20250527/Herwig_noISR_ALL_thrust_no_event_sel_tgenBefore.root", "tree" : "tgenBefore", "branches" : ["px", "py", "pz"]},
        "Sherpa" : {"path": "/global/cfs/cdirs/m3246/aleph/derivations/20250527/Sherpa_noISR_ALL_thrust_no_event_sel_tgenBefore.root", "tree" : "tgenBefore", "branches" : ["px", "py", "pz"]},
}

# # file paths
# mc_paths = {
#     "ArchivedPYTHIA6" : {
#         "path" : "/pscratch/sd/b/badea/aleph/data/ThrustDerivation/042825/alephMCRecoAfterCutPaths_1994_thrust_CleanNeutralAndConversion.root",
#         "tree" : "tgenBefore", # t=reco, tgen = generator level after hadronic event selection, tgenBefore = generator level before hadronic event selection
#         "branches" : ["px", "py", "pz"]
#     },
#     "PYTHIA8" : {
#         "path" : "/pscratch/sd/b/badea/aleph/data/LEP1MCVariations/hannah/ThrustDerivation/042925/LEP1_PYTHIA8_MC_TGenBefore_thrust.root",
#         "tree" : "tgenBefore",
#         "branches" : ["px", "py", "pz"]
#     },
# }

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

def padAwkward(a, maxN, value=0):
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

def loadDataParticles(filePath, treeName, branches, maxNPart, padValue = np.nan):
    # kinematics
    data = []
    with uproot.open(filePath) as rFile:
        tree = rFile[treeName]
        for key in branches:
            if key in tree.keys():
                # temp = loadBranchAndPad(tree[key], maxNPart, value=padValue) # pad to get to maxNPart
                temp = tree[key].array()
                temp = padAwkward(temp, maxNPart, value=padValue) # pad to get to maxNPart
                data.append(temp)
            else:
                print(f"Attempted to add key not in tree: {key} for file {filePath}")
    # stack and return
    data = np.stack(data, axis=-1)
    return data

def convert_PxPyPz_to_EtaPhiPmag(x):
    # x shape = [events, particles, 3], where 3 = [px,py,pz]
    # output = [events, particles, 3], where 3 = [eta, phi, pmag]

    px = x[:,:,0]
    py = x[:,:,1]
    pz = x[:,:,2]
    
    pmag = np.sqrt(px**2 + py**2 + pz**2)
    theta = np.arccos(pz / pmag)
    eta = -np.log(np.tan(theta / 2)) # used because the PET assumes eta provided for the pair wise distance. It doesn't handle theta where you have to wrap around for 2pi
    phi = np.arctan2(py, px)
    out = np.stack([eta, phi, pmag], axis=-1)

    return out

def create_train_val_datasets(data_0, data_1, test_size=0.2, batch_size=512, normalize=True, standardize=False):
    
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
    data_mean, data_std = None, None    
    if standardize:
        print("Number of nan entries pre-standardization: ", np.isnan(data).sum())
        
        # assume masked data is np.nan
        data_mean = np.nanmean(data, axis=(0, 1), keepdims=True)  # shape: (1, 1, F)
        data_std = np.nanstd(data, axis=(0, 1), keepdims=True)    # shape: (1, 1, F)
        print("Pre-standardization mean: ", data_mean.shape, data_mean)
        print("Pre-standardization std: ",data_std.shape, data_std)

        # Standardize: only apply to valid values
        standardized_data = (data - data_mean) / data_std
        standardized_data_mean = np.nanmean(standardized_data, axis=(0, 1), keepdims=True)  # shape: (1, 1, F)
        standardized_data_std = np.nanstd(standardized_data, axis=(0, 1), keepdims=True)    # shape: (1, 1, F)
        print("Post-standardization mean: ",standardized_data_mean.shape, standardized_data_mean)
        print("Post-standardization std: ", standardized_data_std.shape, standardized_data_std)

        # update data
        data = standardized_data
        print("Number of nan entries post-standardization: ", np.isnan(data).sum())

    # convert from np.nan to 0 for masked data
    data = np.nan_to_num(data, nan=0)  # convert np.nan to 0
    print("Number of nan entries after converting nan to 0: ", np.isnan(data).sum())

    # Split into training (80%) and validation (20%) sets randomly
    train_data, val_data, train_labels, val_labels = train_test_split(data, labels, test_size=test_size, stratify=labels)
    # print(train_data.shape, val_data.shape, train_labels.shape, val_labels.shape)
    
    # Create TensorFlow datasets
    train_dataset = tf.data.Dataset.from_tensor_slices((train_data, train_labels))
    val_dataset = tf.data.Dataset.from_tensor_slices((val_data, val_labels))

    # Batch and shuffle training data
    train_dataset = train_dataset.shuffle(buffer_size=len(train_data)).batch(batch_size).prefetch(tf.data.AUTOTUNE)
    val_dataset = val_dataset.batch(batch_size).prefetch(tf.data.AUTOTUNE)

    return train_dataset, val_dataset, data_mean, data_std


if __name__ == "__main__":
    dtype = "Sherpa"
    aleph_mc = loadDataParticles(
        filePath = mc_paths[dtype]["path"],
        treeName = mc_paths[dtype]["tree"],
        branches = mc_paths[dtype]["branches"],
        maxNPart = 80,
    )
    aleph_mc = convert_PxPyPz_to_EtaPhiPmag(aleph_mc)
    train_dataset_step1, val_dataset_step1, data_mean_step1, data_std_step1 = create_train_val_datasets(data_0=aleph_mc, data_1=aleph_mc, test_size=0.2, batch_size=512, normalize=True, standardize=True)
    print(data_mean_step1, data_std_step1)

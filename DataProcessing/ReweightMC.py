import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.callbacks import EarlyStopping,ModelCheckpoint, ReduceLROnPlateau
from sklearn.model_selection import train_test_split
import numpy as np
import uproot
import pickle
import os
import awkward as ak

import omnifold

def expit(x):
    return 1. / (1. + np.exp(-x))

def reweight(data, model, batch_size, verbose=True):
    f = expit(model.predict(data,batch_size=batch_size,verbose=verbose))
    weights = f / (1. - f)  # this is the crux of the reweight, approximates likelihood ratio
    weights = np.nan_to_num(weights[:,0],posinf=1)
    return weights

def loadData(filePath, treeName, vars, SystematicVariation=None):
    i = SystematicVariation
    data = []
    with uproot.open(filePath) as rFile:
        for var in vars:
            if SystematicVariation == None:
                temp = np.array(rFile[f"{treeName}/{var}"])
            else:
                temp = np.array([x[i] for x in np.array(rFile[f"{treeName}/{var}"])]) # reco level
            data.append(temp)
    # stack to form data
    data = np.stack(data,axis=1)
    return data

def loadBranchAndPad(branch, maxN, value=0):
    a = branch.array()
    a = ak.to_numpy(ak.fill_none(ak.pad_none(a, max(maxN,np.max(ak.num(a)))),value)) # must take maximum so that it pads to a uniform max and then cut down
    a = a[:,:maxN]
    return a

def loadDataParticles(filePath, treeName, branches, maxNPart):
    # kinematics
    data = []
    with uproot.open(filePath) as rFile:
        tree = rFile[treeName]
        # print(tree.num_entries)
        # print(tree.keys())
        for key in branches:
            if key in tree.keys():
                temp = loadBranchAndPad(tree[key], maxNPart) # pad to get to maxNPart
                data.append(temp)
            else:
                print(f"Attempted to add key not in tree: {key} for file {filePath}")
    # stack and return
    data = np.stack(data, axis=-1)
    return data

# copied from https://github.com/ViniciusMikuni/OmniLearn/blob/main/scripts/omnifold.py#L13C1-L20C32
# def weighted_binary_crossentropy(y_true, y_pred):
#     """Custom loss function with weighted binary cross-entropy."""

#     weights = tf.cast(tf.gather(y_true, [1], axis=1), tf.float32)  # Event weights
#     y_true = tf.cast(tf.gather(y_true, [0], axis=1), tf.float32)  # Actual labels

#     # Compute loss using TensorFlow's built-in function to handle numerical stability
#     loss = weights * tf.nn.sigmoid_cross_entropy_with_logits(labels=y_true, logits=y_pred)
#     return tf.reduce_mean(loss)

mc_paths = {
    "ArchivedPYTHIA6" : {
        "path" : "",
        "tree" : "t"
    },
    "HERWIG7" : {
        "path" : "/global/homes/b/badea/aleph/data/LEP1MCVariations/abaty/HERWIG7/2_10_2024_LEP1MC/LEP-Matchbox-S1000-1_0_0.root",
        "tree" : "t"
    },
    "SHERPA" : {
        "path" : "/global/homes/b/badea/aleph/data/LEP1MCVariations/abaty/SHERPA/2_10_2024_LEP1MC/Sherpa_RNG100_0_0.root",
        "tree" : "t"
    },
    "PYTHIA8" : {
        "path" : "/global/homes/b/badea/aleph/data/LEP1MCVariations/hannah/LEP1_pythia8_MC_withSphericity_v2.root",
        "tree" : "tgen"
    },
    "PYTHIA8_DIRE" : {
        "path": "/global/homes/b/badea/aleph/data/LEP1MCVariations/hannah/LEP1_pythia8_MC_DIRE.root",
        "tree" : "tgen"
    },
    "PYTHIA8_VINCIA" : {
        "path" : "/global/homes/b/badea/aleph/data/LEP1MCVariations/hannah/LEP1_pythia8_MC_VINCIA.root",
        "tree" : "tgen"
    }
}

# training settings
test_size = 0.2
lr = 5e-4
epochs = 1
batch_size = 512
verbose = True
weights_folder = "./"
new_mc_name = "HERWIG7"
model_name = os.path.join(weights_folder, f'Reweight_{new_mc_name}.weights.h5')
print(model_name)

# Event level distribution reweighting
# # questions
# # do we reweight before particle or event selections?
# aleph_mc = loadData(
#     filePath = "/global/homes/b/badea/aleph/data/ThrustDerivation/030725/alephMCRecoAfterCutPaths_1994_thrust.root", 
#     treeName = "t", # t=reco, tgen = generator level after hadronic event selection, tgenBefore = generator level before hadronic event selection
#     vars = ["Thrust"],
#     SystematicVariation = 0 # no event selection, use all reco events
# )
# new_mc = loadData(
#     filePath = new_mc_paths[new_mc_name]["path"], 
#     treeName = new_mc_paths[new_mc_name]["tree"],
#     vars = ["Thrust"],
#     SystematicVariation = None # 
# )
# print(aleph_mc.shape, new_mc.shape)


# Particle level distribution reweighting 
aleph_mc = loadDataParticles(
    filePath = mc_paths["PYTHIA8"]["path"], # "/global/homes/b/badea/aleph/data/ThrustDerivation/030725/alephMCRecoAfterCutPaths_1994_thrust.root", 
    treeName = mc_paths["PYTHIA8"]["tree"], # "t", # t=reco, tgen = generator level after hadronic event selection, tgenBefore = generator level before hadronic event selection
    branches = ["px", "py", "pz", "mass", "isCharged"],
    maxNPart = 80
)
print(aleph_mc.shape)
new_mc = loadDataParticles(
    filePath = mc_paths[new_mc_name]["path"], 
    treeName = mc_paths[new_mc_name]["tree"],
    branches = ["px", "py", "pz", "mass", "charge"],
    maxNPart = 80
)
print(new_mc.shape)

# create labels
labels_aleph_mc = np.zeros(len(aleph_mc),dtype=np.float32)
labels_new_mc = np.ones(len(new_mc),dtype=np.float32)

# concatenate
data = np.concatenate((aleph_mc, new_mc))
# labels expected to be stack of [labels, weights]
labels = np.concatenate((labels_aleph_mc, labels_new_mc))
weights = np.ones(len(labels)) # unit weights for all for now
labels = np.stack([labels, weights], axis=1).astype(np.float32)
print(data.shape, labels.shape)

# Split into training (80%) and validation (20%) sets randomly
train_data, val_data, train_labels, val_labels = train_test_split(data, labels, test_size=test_size, stratify=labels)
print(train_data.shape, val_data.shape, train_labels.shape, val_labels.shape)

# Create TensorFlow datasets
train_dataset = tf.data.Dataset.from_tensor_slices((train_data, train_labels))
val_dataset = tf.data.Dataset.from_tensor_slices((val_data, val_labels))

# Batch and shuffle training data
train_dataset = train_dataset.shuffle(buffer_size=len(train_data)).batch(batch_size).prefetch(tf.data.AUTOTUNE)
val_dataset = val_dataset.batch(batch_size).prefetch(tf.data.AUTOTUNE)
# exit()

# prepare MLP network
# ndim = train_data.shape[1]
# layer_sizes = [100, 100, 100]
# model = omnifold.MLP(ndim, layer_sizes = layer_sizes, activation="relu")

# prepare PET network
model = omnifold.PET(num_feat = train_data.shape[2], num_evt = 0, num_part = train_data.shape[1], num_heads = 1, num_transformer = 1, local = False, projection_dim = 64)

# training settings
opt = tf.keras.optimizers.Adam(learning_rate=lr)
model.compile(opt, loss = omnifold.net.weighted_binary_crossentropy)
callbacks = [
    ReduceLROnPlateau(patience=1000, min_lr=1e-7, verbose=verbose, monitor="val_loss"),
    EarlyStopping(patience=10, restore_best_weights=True,  monitor="val_loss"),
    ModelCheckpoint(model_name, save_best_only=True, mode='auto', save_weights_only=True),
]

hist =  model.fit(
    train_dataset,
    epochs = epochs,
    validation_data = val_dataset,
    verbose = verbose,
    callbacks = callbacks
)

print(f"Last val loss {hist.history['val_loss'][0]}")
print("INFO: Dumping training history ...")
with open(model_name.replace(".weights.h5",".pkl"),"wb") as f:
    pickle.dump(hist.history, f)

weights = reweight(new_mc, model, batch_size=10000)
print(weights.shape)
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.callbacks import EarlyStopping,ModelCheckpoint, ReduceLROnPlateau
from sklearn.model_selection import train_test_split
import numpy as np
import uproot
import pickle
import os
import awkward as ak
import random

import omnifold
from ReweightMCDataLoading import *

# training settings
test_size = 0.2
lr = 5e-4
epochs = 10
batch_size = 512
verbose = True
top_dir = "/pscratch/sd/b/badea/aleph/unfold-ee-logtau/ReweightMC/results/"
weights_folder = os.path.join(top_dir, f'training-{"%08x" % random.randrange(16**8)}') # "./"
os.makedirs(weights_folder, exist_ok=True)
new_mc_name = "PYTHIA8"
model_name = os.path.join(weights_folder, f'Reweight_{new_mc_name}.weights.h5')
print(model_name)
branches = ["px", "py", "pz", "mass"] # note: "charge" only exists for files processed with scan.cc (aleph mc, herwig, sherpa) not for the current pythia8 files

if __name__ == "__main__":
    
    # Particle level distribution reweighting 
    aleph_mc = loadDataParticles(
        filePath = mc_paths["ArchivedPYTHIA6"]["path"],
        treeName = mc_paths["ArchivedPYTHIA6"]["tree"],
        branches = branches,
        maxNPart = 80
    )
    print(aleph_mc.shape)
    new_mc = loadDataParticles(
        filePath = mc_paths[new_mc_name]["path"], 
        treeName = mc_paths[new_mc_name]["tree"],
        branches = branches,
        maxNPart = 80
    )
    print(new_mc.shape)

    # create labels
    labels_aleph_mc = np.ones(len(aleph_mc),dtype=np.float32)
    labels_new_mc = np.zeros(len(new_mc),dtype=np.float32)

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

    weights = reweight(aleph_mc, model, batch_size=10000)
    print(weights.shape)
    np.save(model_name.replace(".weights.h5",".reweight.npy"), weights)

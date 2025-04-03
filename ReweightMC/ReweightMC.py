import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.callbacks import EarlyStopping,ModelCheckpoint, ReduceLROnPlateau
# from sklearn.model_selection import train_test_split
import numpy as np
# import uproot
import pickle
import os
# import awkward as ak
import random

import omnifold
from ReweightMCDataLoading import *

# fitting function
def doFitAndEvaluate(model, train_dataset, val_dataset, epochs, verbose, callbacks, name, evaluate_data):
    
    hist =  model.fit(
        train_dataset,
        epochs = epochs,
        validation_data = val_dataset,
        verbose = verbose,
        callbacks = callbacks
    )

    print(f"Last val loss {hist.history['val_loss'][0]}")
    print("INFO: Dumping training history ...")
    with open(name.replace(".weights.h5",".pkl"),"wb") as f:
        pickle.dump(hist.history, f)

    weights = reweight(evaluate_data, model, batch_size=10000)
    print(weights.shape)
    np.save(name.replace(".weights.h5",".reweight.npy"), weights)    

    return hist

if __name__ == "__main__":
    
    # training settings
    test_size = 0.2
    lr = 5e-4
    epochs = 2
    batch_size = 512
    verbose = True
    top_dir = "/pscratch/sd/b/badea/aleph/unfold-ee-logtau/ReweightMC/results/"
    weights_folder = os.path.join(top_dir, f'training-{"%08x" % random.randrange(16**8)}') # "./"
    os.makedirs(weights_folder, exist_ok=True)
    new_mc_name = "PYTHIA8"
    maxNPart = 80

    # Particle level distribution reweighting 
    aleph_mc = loadDataParticles(
        filePath = mc_paths["ArchivedPYTHIA6"]["path"],
        treeName = mc_paths["ArchivedPYTHIA6"]["tree"],
        branches = mc_paths["ArchivedPYTHIA6"]["branches"],
        maxNPart = maxNPart
    )
    print(aleph_mc.shape)
    new_mc = loadDataParticles(
        filePath = mc_paths[new_mc_name]["path"], 
        treeName = mc_paths[new_mc_name]["tree"],
        branches = mc_paths[new_mc_name]["branches"],
        maxNPart = maxNPart
    )
    print(new_mc.shape)

    # prepare datasets
    train_dataset_step1, val_dataset_step1 = create_train_val_datasets(data_0=aleph_mc, data_1=aleph_mc, test_size=test_size, batch_size=batch_size, normalize=True)
    train_dataset_step2, val_dataset_step2 = create_train_val_datasets(data_0=new_mc, data_1=aleph_mc, test_size=test_size, batch_size=batch_size, normalize=True)

    # prepare PET network
    model = omnifold.PET(
        num_feat = aleph_mc.shape[2], 
        num_evt = 0, 
        num_part = maxNPart, 
        num_heads = 1, 
        num_transformer = 1, 
        local = False,
        projection_dim = 64
    )

    # training settings
    opt = tf.keras.optimizers.Adam(learning_rate=lr)
    model.compile(opt, loss = omnifold.net.weighted_binary_crossentropy)
    callbacks = [
        ReduceLROnPlateau(patience=1000, min_lr=1e-7, verbose=verbose, monitor="val_loss"),
        EarlyStopping(patience=10, restore_best_weights=True,  monitor="val_loss"),
    ]

    # Step 1 pretrain the model to get it to unity
    model_name = os.path.join(weights_folder, f'Reweight_Step1_ALEPHMC.weights.h5')
    callbacks.append(ModelCheckpoint(model_name, save_best_only=True, mode='auto', save_weights_only=True))
    print("Running Step 1 with model name: ", model_name)
    hist = doFitAndEvaluate(model, train_dataset_step1, val_dataset_step1, epochs, verbose, callbacks, model_name, aleph_mc) 

    # Step 2 train the model for reweighting
    model_name = os.path.join(weights_folder, f'Reweight_Step2_{new_mc_name}.weights.h5')
    callbacks.append(ModelCheckpoint(model_name, save_best_only=True, mode='auto', save_weights_only=True))
    print("Running Step 2 with model name: ", model_name)
    hist = doFitAndEvaluate(model, train_dataset_step2, val_dataset_step2, epochs, verbose, callbacks, model_name, aleph_mc)

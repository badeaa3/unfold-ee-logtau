'''
Author: Anthony Badea
Date: January 30, 2022

A python script that takes in the files and outputs a numpy array with the thrust value, event selections, and event weights.

File names:
data_fileName = 'LEP1Data{}_recons_aftercut-MERGED.root'.format(year), with year = 1992, 1993, 1994, 1995
mc_fileName = 'alephMCRecoAfterCutPaths_1994.root'

Tree definitions:
- t = data or reco
- tgen = generation level + hadronic event selection
- tgenBefore = generation level without hadronic event selection

Notes:
- see to understand event/track selections https://www.dropbox.com/scl/fi/gqe7qgm4ygnr7xontuke4/ALEPH-Omnifold-073020.pdf?dl=0

Systematic uncertainty variation from https://arxiv.org/pdf/1906.00489.pdf:
The required number of hits a track leaves in the ALEPH time projection chamber was varied from 4 to 7. 
From this variation, the tracking uncertainty is estimated to be 0.7% in the lab coordinate analysis and 0.3% in the thrust coordinate analysis. 
The hadronic event selection was studied by changing the required charged energy in an event to be 10 instead of 15 GeV. 
... 
An additional systematic of 0.2%-10% (0.1%-0.5%) in the lab (thrust) coordinate analysis is included to quantify the residual 
uncertainty in the reconstruction effect correction factor derived from the pythia 6.1 archived MC sample, which is mainly from the 
limited size of the archived MC sample.
'''

import uproot
import numpy as np
import argparse

def process(inFileName, outFileName):

    # load the data file and infer type (data or mc)
    f = uproot.open(inFileName)
    
    # list of relevenat event selections
    event_selections = [i for i in f['t'].keys() if 'passEventSelection' in i]

    if "LEP1Data" in inFileName: # data

        output = {f"t_{evsel}" : np.stack(np.array(f["t"][evsel])) for evsel in event_selections}
        output["t_thrust"] = np.stack(np.array(f["t"]["Thrust"]))
        output["t_uniqueID"] = np.array(f["t"]["uniqueID"])
        # return output

    else: # MC  

        output = {
            "t_thrust" : np.stack(np.array(f["t"]["Thrust"])), # reco (after detector simulation)
            "tgen_thrust" : np.array(f["tgen"]["Thrust"]), # with hadronic event selection
            "tgenBefore_thrust" : np.array(f["tgenBefore"]["Thrust"]), # without hadronic event selection
            "t_uniqueID" : np.array(f["t"]["uniqueID"]),
            "tgen_uniqueID" : np.array(f["t"]["uniqueID"]),
            "tgenBefore_uniqueID" : np.array(f["t"]["uniqueID"]),
        }
        # construct the final event selection
        for evsel in event_selections:
            output[f"t_{evsel}"] = np.stack(np.array(f["t"][evsel]))
        # no event selection on gen
        output[f"tgen_passEventSelection"] = np.ones(f['tgen'].num_entries)
        output[f"tgenBefore_passEventSelection"] = np.ones(f['tgenBefore'].num_entries)
            
        # return output

    # print summary
    for key, val in output.items():
        print(f"{key} of size {val.shape}")

    # save to file
    np.savez(outFileName, **output)

    return output

if __name__ == "__main__":

    # user options
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inFile", help="Input file")
    parser.add_argument("-o", "--outFile", help="Output file name", default="test.npz")
    ops = parser.parse_args()

    # perform reprocessing
    output = process(ops.inFile, ops.outFile)

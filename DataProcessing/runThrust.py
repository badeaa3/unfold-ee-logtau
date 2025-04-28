'''
Author: Anthony Badea
Date: April 30, 2022
Purpose: Wrapping job launcher for Thrust.cxx to multiprocess the event loop
'''

import os
import argparse
import multiprocessing as mp

def main():

    # user options
    ops = options()

    # make confs
    confs = []
    for thisdiv in range(ops.ndivs):
        confs.append({
                "i" : ops.inFile,
                "o" : ops.outFile.replace(".root", f"_{thisdiv}.root"),
                "divide" : ops.ndivs,
                "thisdiv" : thisdiv,
                "dryrun" : ops.dryrun,
                "debug" : ops.debug
            })

    # launch jobs
    if len(confs) == 1:
        runThrust(confs[0])
    else:
        results = mp.Pool(ops.ncpu).map(runThrust, confs)

    # hadd
    jobs = " ".join([d["o"] for d in confs])
    cmds = [f"hadd -f {ops.outFile} {jobs}", f"rm -rf {jobs}"]
    for cmd in cmds: 
        print("\n" + cmd)
    if not ops.dryrun:
        for cmd in cmds:
            os.system(cmd)

def options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inFile", help="Input file", required=True)
    parser.add_argument("-o", "--outFile", help="Output ROOT file", required=True)
    parser.add_argument("-n", "--ndivs", help="Number of divisions to partition the input file into", default=1, type=int)
    parser.add_argument("-j", "--ncpu", help="Number of cores to use for multiprocessing", default=1, type=int)
    parser.add_argument("--dryrun", help="Do not execute any commands", action="store_true")
    parser.add_argument("--debug", help="Enable debug mode", action="store_true")
    return parser.parse_args()

def runThrust(conf):
    cmd = f"./Thrust.exe -i {conf['i']} -o {conf['o']} --divide {conf['divide']} --thisdiv {conf['thisdiv']}"
    if conf["debug"]:
        cmd += " --debug"
    print(cmd)
    if not conf["dryrun"]:
        os.system(cmd)

if __name__ == "__main__":
    main()

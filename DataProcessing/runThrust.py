'''
Author: Anthony Badea
Date: April 30, 2022
Purpose: Wrapping job launcher for Thrust.cxx to multiprocess selection variations
'''

import os
import argparse
import multiprocessing as mp

# Global configuration for different file types
FILE_CONFIGS = {
    'LEP1Data': {
        'inFileType': 'LEP1Data',
        'treeNames': ['t']
    },
    'PYTHIA8': {
        'inFileType': 'PYTHIA8',
        'treeNames': ['tgenBefore']
    },
    'Herwig': {
        'inFileType': 'Herwig',
        'treeNames': ['tgenBefore']
    },
    'Sherpa': {
        'inFileType': 'Sherpa',
        'treeNames': ['tgenBefore']
    },
    'ALEPHMC': {  # Default
        'inFileType': 'ALEPHMC',
        'treeNames': ['t', 'tgen', 'tgenBefore']
    }
}

# Default selection parameters
DEFAULT_SELECTION = {
    'nTPCcut': 4,
    'chargedTracksAbsCosThCut': 0.94,
    'ptCut': 0.2,  # GeV
    'd0Cut': 2,    # cm
    'z0Cut': 10,   # cm
    'ECut': 0.4,   # GeV
    'neutralTracksAbsCosThCut': 0.98,
    'TotalTrkEnergyCut': 15,   # GeV
    'AbsCosSThetaCut': 0.82,
    'NTrkCut': 5,
    'NeuNchCut': 13,
    'EVisCut': 0,     # GeV
    'MissPCut': 9999, # GeV
    'keepChargedTracks': True,
    'keepNeutralTracks': True,
    'doMET': False
}

def config_to_string(config):
    """Convert configuration dictionary to command line string."""
    parts = []
    for key, value in config.items():
        if isinstance(value, bool):
            parts.append(f"{key}={'true' if value else 'false'}")
        elif isinstance(value, list):
            parts.append(f"{key}={','.join(map(str, value))}")
        else:
            parts.append(f"{key}={value}")
    return ";".join(parts)

def get_file_config(input_file):
    """Determine file type and tree names based on input file."""
    filename = os.path.basename(input_file)
    
    for file_type, config in FILE_CONFIGS.items():
        if file_type.lower() in filename.lower():
            return config.copy()
    
    # raise error if no known file type is found
    raise ValueError(f"Unknown file type for '{filename}'. Supported types: {list(FILE_CONFIGS.keys())}")

def create_selection_variations():
    """Create all selection variations to run."""
    variations = []
    
    # No event selections
    no_event_sel = DEFAULT_SELECTION.copy()
    no_event_sel.update({
        'TotalTrkEnergyCut': 0, 
        'AbsCosSThetaCut': 1, 
        'NTrkCut': 0, 
        'NeuNchCut': 0, 
        'EVisCut': 0
    })
    variations.append(('no_event_sel', no_event_sel))
    
    # Nominal selection
    variations.append(('nominal', DEFAULT_SELECTION.copy()))
    
    # nTPC variation: 4 -> 7
    ntpc7 = DEFAULT_SELECTION.copy()
    ntpc7['nTPCcut'] = 7
    variations.append(('ntpc7', ntpc7))
    
    # Charged tracks pT variation: 0.2 -> 0.4 GeV
    pt04 = DEFAULT_SELECTION.copy()
    pt04['ptCut'] = 0.4
    variations.append(('pt04', pt04))
    
    # Total charged energy variation: 15 -> 10 GeV
    ech10 = DEFAULT_SELECTION.copy()
    ech10['TotalTrkEnergyCut'] = 10
    variations.append(('ech10', ech10))
    
    # Thrust without neutral objects
    no_neutrals = DEFAULT_SELECTION.copy()
    no_neutrals['keepNeutralTracks'] = False
    variations.append(('no_neutrals', no_neutrals))
    
    # Thrust with missing momentum vector
    with_met = DEFAULT_SELECTION.copy()
    with_met['doMET'] = True
    variations.append(('with_met', with_met))
    
    return variations


def run_selection_variation(ops, file_config, sel_name, selection):
    """Run all divisions for a single selection variation across all trees and merge results."""
    tree_names = file_config['treeNames']
    
    # Process each tree separately
    for tree_name in tree_names:
        
        # Filter: generator-level trees (tgen, tgenBefore) should only run with no_event_sel
        if tree_name in ['tgen', 'tgenBefore'] and sel_name != 'no_event_sel':
            print(f"\n[{sel_name}] Skipping tree: {tree_name} (generator-level trees only run with no_event_sel)")
            continue
            
        print(f"\n[{sel_name}] Processing tree: {tree_name}")
        run_tree_jobs(ops, file_config, sel_name, selection, tree_name)

def run_tree_jobs(ops, file_config, sel_name, selection, tree_name):
    """Run all divisions for a single tree within a selection variation and merge results."""
    base_outfile = os.path.join(ops.outDir, os.path.splitext(os.path.basename(ops.inFile))[0])
    base_outfile = base_outfile + "_thrust"
    
    # Create job configurations for this tree
    job_configs = []
    for division in range(ops.ndivs):
        # Create output filename with tree name
        if ops.ndivs > 1:
            outfile = f"{base_outfile}_{sel_name}_{tree_name}_{division}.root"
        else:
            outfile = f"{base_outfile}_{sel_name}_{tree_name}.root"
        
        # Create config for single tree (not a list)
        tree_config = {
            'inFileType': file_config['inFileType'],
            'treeNames': tree_name  # Single tree name, not a list
        }
        
        # Combine tree config with selection parameters
        combined_config = {}
        combined_config.update(tree_config)
        combined_config.update(selection)
        
        job_configs.append({
            "input_file": ops.inFile,
            "output_file": outfile,
            "config_string": config_to_string(combined_config),
            "selection_name": f"{sel_name}_{tree_name}",
            "total_divisions": ops.ndivs,
            "current_division": division,
            "dry_run": ops.dryrun,
            "debug": ops.debug
        })
    
    print(f"  Running {ops.ndivs} divisions for tree {tree_name}...")
    
    # Execute all divisions for this tree using multiprocessing
    if len(job_configs) == 1:
        run_thrust_job(job_configs[0])
    else:
        with mp.Pool(ops.ncpu) as pool:
            results = pool.map(run_thrust_job, job_configs)
    
    # Merge divided files if necessary
    if ops.ndivs > 1:
        final_outfile = f"{base_outfile}_{sel_name}_{tree_name}.root"
        division_files = [f"{base_outfile}_{sel_name}_{tree_name}_{div}.root" for div in range(ops.ndivs)]
        
        # Merge files
        hadd_cmd = f"hadd -f {final_outfile} {' '.join(division_files)}"
        cleanup_cmd = f"rm -f {' '.join(division_files)}"
        
        print(f"  [{tree_name}] Merging files...")
        for cmd in [hadd_cmd, cleanup_cmd]:
            print(f"    {cmd}")
            if not ops.dryrun:
                os.system(cmd)
    
    print(f"  [{tree_name}] Completed!")

def run_thrust_job(job_config):
    """Execute a single thrust job."""
    cmd_parts = [
        "./Thrust.exe",
        f"-i {job_config['input_file']}",
        f"-o {job_config['output_file']}",
        f"-c \"{job_config['config_string']}\"",
        f"--divide {job_config['total_divisions']}",
        f"--thisdiv {job_config['current_division']}"
    ]
    
    if job_config["debug"]:
        cmd_parts.append("--debug")
    
    cmd = " ".join(cmd_parts)
    print(f"[{job_config['selection_name']}] {cmd}")
    
    if not job_config["dry_run"]:
        return os.system(cmd)
    return 0

def main():
    """Main execution function."""
    ops = parse_arguments()
    
    # Get file configuration and create variations
    file_config = get_file_config(ops.inFile)
    variations = create_selection_variations()
    
    total_trees = len(file_config['treeNames'])
    print(f"Running {len(variations)} selection variations with {total_trees} trees each, {ops.ndivs} divisions per tree")
    
    # Execute each selection variation sequentially, completing each one fully
    for sel_name, selection in variations:
        run_selection_variation(ops, file_config, sel_name, selection)

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Launch Thrust.cxx jobs with multiple selection variations")
    parser.add_argument("-i", "--inFile", required=True, help="Input ROOT file")
    parser.add_argument("-o", "--outDir", default="./", help="Output directory (default: ./)")
    parser.add_argument("-n", "--ndivs", type=int, default=1, help="Number of divisions to partition the input file into (default: 1)")
    parser.add_argument("-j", "--ncpu", type=int, default=1, help="Number of cores to use for multiprocessing (default: 1)")
    parser.add_argument("--dryrun", action="store_true", help="Print commands without executing them")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode in Thrust.exe")
    return parser.parse_args()

if __name__ == "__main__":
    main()

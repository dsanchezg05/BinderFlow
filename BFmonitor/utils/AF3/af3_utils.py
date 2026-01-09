'''
Script with utility functions for running AF3 in BFmonitor
'''

def extract_sequence_from_pdb(pdb_file):
    '''
    Extract the sequence from a PDB file using Biopython

    Input:
    - pdb_file: path to the PDB file

    Output:
    - sequence: extracted sequence as a string
    '''
    from Bio import SeqIO
    alphabetical = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    sequence = {}
    with open(pdb_file, 'r') as file:
        counter = 0
        for record in SeqIO.parse(file, 'pdb-atom'):
            sequence[alphabetical[counter]] = str(record.seq)
            counter += 1
    return sequence


def AF3_target_json_generator(directory, dir):
    '''
    Script to generate the JSON files required for AF3 for the target only
    '''
    import os 
    import sys
    from Bio import SeqIO
    from pathlib import Path
    import subprocess
    import json
    import random
    import warnings
    warnings.filterwarnings("ignore")

    seed = random.randint(1, 10000)
    # Select just one PDB file from the provided directory
    pdb_list = [str(file) for file in Path(directory).glob("*.pdb")]
    pdb_file = pdb_list[0] if pdb_list else None

    if pdb_file is not None:
        try:
            os.listdir(f"{dir}/AF3_run/AF3_json/")
        except:
            subprocess.run(["mkdir -p "+str(dir)+"/AF3_run/AF3_json"], shell=True)

    target_sequence = extract_sequence_from_pdb(pdb_file)['B']# Target sequence is just the B chain

     # Load the settings

    sequences_json1=[{
        "protein":{
        "id": 'B',
        "sequence": target_sequence
        }
    }]
    data = {
        "name": "msa_template",
        "modelSeeds": [seed],
        "sequences": sequences_json1,
        "dialect": "alphafold3",
        "version": 1
        }
    if not os.path.exists(f"{dir}/AF3_run/AF3_json/Target_data/Template_data.json"):
        subprocess.run(["mkdir -p "+str(dir)+"/AF3_run/AF3_json/Target_data"], shell=True)
        with open(f"{dir}/AF3_run/AF3_json/Target_data/Template_data.json", "w") as f:
            json.dump(data, f, indent=2)
    else:
        print("Target JSON already exists, skipping...")


def load_bash_vars(file_path):
    '''
    Script to load the bash variables from a .sh file
    Input:
    - file_path: path to the .sh file
    Output:
    - env: dictionary with the bash variables
    '''
    import subprocess
    # Run a bash shell that sources the file, then prints all environment vars
    result = subprocess.run(
        ['bash', '-c', f'set -a && source {file_path} && env'],
        capture_output=True, text=True, check=True
    )
    env = {}
    for line in result.stdout.splitlines():
        key, _, value = line.partition("=")
        env[key] = value
    return env

def AF3_binder_json_generator(directory, msa_dir,dir):
    '''
    Script to generate the JSON files required for AF3 starting from the PDB files

    Input:

    - pdbs_list: list of PDB files to process

    Output:

    - Template_data_X.json and Inference_data_X.json files in AF3_run/AF3_json/Structure_<pdb_name>/ for inference with AF3 of each PDB in the list

    In the code JSON1 stands for the data part (MSA of the target) while JSON2 stands for the inference part (binder+target) 
    '''
    import os 
    import sys
    from Bio import SeqIO
    from pathlib import Path
    import subprocess
    import json
    import random
    import warnings
    warnings.filterwarnings("ignore")
    # Get the seed
    seed = random.randint(1, 10000)
    # Extract the list of pdb files from the provided directory

    pdb_list = [str(file) for file in Path(directory).glob("*.pdb")]

    # If the pdb list is not empty, proceed to generate the JSON folder 
    if pdb_list != []:
        try:
            os.listdir(f"{dir}/AF3_run/AF3_json/")
        except:
            subprocess.run(["mkdir -p "+str(dir)+"/AF3_run/AF3_json"], shell=True)
        
        alphabetical = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";idx=0 
        pdb_dict={};run=0
        # Extract the sequences from the PDB files
        for file in pdb_list:
            pdb_name = file.split('/')[-1].split('.')[0]
            sequences_json2=[] ; num_seqs=0
            # Extract the sequences from the PDB file. Probably this can be done from the silent file, so there is no need to extract the PDBs in the first place
            sequences = extract_sequence_from_pdb(file)
            for chain, seq in sequences.items():
                # Create the JSON structure for AF3
                ## Json part of the binder
                if chain == 'A':
                    sequences_json2.append({
                        "protein":{
                            "id": chain,
                            "sequence": str(sequences[chain]),
                            "unpairedMsa": "",
                            "pairedMsa": "",
                        }
                    })
                elif chain == 'B':
                    sequences_json2.append({
                        "protein":{
                            "id": chain,
                            "sequence": str(sequences[chain]),
                            "unpairedMsaPath": str(f"{msa_dir}/Unpaired_msa.a3m"),
                            "pairedMsaPath": str(f"{msa_dir}/Paired_msa.a3m"),
                        }
                    })

                # General info for inference jsons
                data_json2 = {
                    "name": "inference_"+str(pdb_name),
                    "modelSeeds": [seed],
                    "sequences": sequences_json2,
                    "dialect": "alphafold3",
                    "version": 1
                }

            try:
                # Check if the PDB has already been processed
                os.listdir(f"{dir}/AF3_run/AF3_json/Structure_"+str(pdb_name))
                print("")
                print(str(pdb_name)+" json already generated, skipping...")
                print("")
            except:
                # If not, create the folder and write the JSON files
                loc_dir="mkdir "+str(dir)+"/AF3_run/AF3_json/Structure_"+str(pdb_name)
                subprocess.run([loc_dir], shell=True)
                with open(f"{dir}/AF3_run/AF3_json/Structure_"+str(pdb_name)+"/Inference_"+f"{pdb_name}.json", "w") as f2:
                    json.dump(data_json2, f2, indent=2)                



def AF3_check_done(selected_dir, num_inferences, data = True):
    '''
    Script to check if AF3 has finished running
    Input:
    - hits_dir: path to the results folder in which there should be the AF3 results
    - num_inferences: indicates the amount of models that have been run though AF3 in inference mode
    - data: boolean indicating if checking for data mode or inference mode
    Output:
    - True if AF3 has finished, False otherwise
    '''
    from pathlib import Path
    import os 

    path = Path(selected_dir)
    if not path.exists():
        return False
    inferences=0
    for entry in path.iterdir():
        if entry.is_dir():
            if data and entry.name == "msa_template":
                # Check for any JSON files inside the subdirectory
                json_files = list(entry.glob('*.json'))
                if not json_files:
                    return False
                else: return True
            elif not data and entry.name.startswith("inference_"):
                print('Detecting inference mode')
                summary_files = list(entry.glob('*summary_confidences.json'))
                if not summary_files:
                    return False
                else: inferences=inferences+1
        else: pass
    if num_inferences == inferences:
        return True
    return False



def AF3_summary_info(dir):

    '''
    Script to get summary information from AF3 predictions.
    Code to get the information from the AF3 JSON file and summarize it.
    It also gets the information from the binder design scoring step

    Output:

    - Scoring_Stats_AF3.csv file with the summary information for all the runs in the provided directory
    '''

    import pandas as pd
    import json
    import glob
    import os 
    import subprocess
    import sys
    SCRIPT_DIR = os.path.dirname(__file__)
    sys.path.append(f'{SCRIPT_DIR}/BFmonitor/utils/AF3/')
    from ipsae import main as ipsae
    import numpy as np
    # Create storing dictionary

    summary_data = {
        "description":[],
        "iPTM":[],
        "chain_pair_pae_min":[],
        "ranking_score":[],
        "plddt_binder_AF3":[],
        "ipSAE_AF3":[],
        "pDOCKQ_AF3":[],
        "pDOCKQ2_AF3":[],
        "LIS_AF3":[]
    }
    # Get all the directories from results
    
    results_dirs = str(f"{dir}/AF3_run/AF3_output")
    scoring_df=pd.read_csv(f"{dir}/Scoring_Stats.csv")
    for folder in glob.glob(f'{results_dirs}/inference*'):
        print(folder)
        description = str(os.path.basename(folder)).replace('inference_', '')
        # AF3 statitics
        summary_file_path= f"{results_dirs}/inference_{description}/inference_{description}_summary_confidences.json"
        if not os.path.exists(summary_file_path):
            print(f"JSON file {summary_file_path} does not exist. Skipping...")
            continue
        with open(summary_file_path, 'r') as file:
            data = json.load(file)
        # Load the confidences JSON file
        confidences_file_path= f"{results_dirs}/inference_{description}/inference_{description}_confidences.json"
        if not os.path.exists(confidences_file_path):
            print(f"JSON file {confidences_file_path} does not exist. Skipping...")
            continue
        with open(confidences_file_path, 'r') as f:
            readed = json.load(f)

        # AF3 plddt binder
        ## Get the number of residues in the target
        atoms=readed["atom_plddts"]
        atoms_chains_ids=readed["atom_chain_ids"]
        atoms_binders = [atoms[i] for i in range(len(atoms)) if atoms_chains_ids[i]=="A"]
        # pLDDT is computed per residue, not per atom
        command = f"""cat {dir}/AF3_run/AF3_output/inference_{description}/inference_{description}_model.cif | grep "ATOM*" | grep "A [0-9]" | awk '{{print $9}}'"""
        residues=subprocess.run([command], shell=True, capture_output=True).stdout
        residues = np.array([int(x) for x in residues.decode().split('\n') if x])
        # Residues extracts the atoms per residue for the chain A (binder) Proabbly Biopython can also be used
        residue_id = set(residues) # Getting unique residue ids
        plddt_per_residue = []
        for res_id in residue_id:
            plddt_per_residue.append(np.array(atoms_binders)[residues == res_id].mean()) # For each residue computing the plddt with its atoms
        plddt_af3 = np.array(plddt_per_residue).mean()
        print('pLDDT binder AF3:', plddt_af3)
        #af3 ipsae
        iPSAE, iPSAE_per_residue = ipsae(confidences_file_path, f"{results_dirs}/inference_{description}/inference_{description}_model.cif", pae_cutoff=15, dist_cutoff=10)
        ipsae_df = pd.DataFrame(iPSAE)
        filtered_ipsae = ipsae_df[ipsae_df['Type'] == 'max']
        ipSAE_value = filtered_ipsae['ipSAE'].values[0]
        pDOCKQ_value = filtered_ipsae['pDockQ'].values[0]
        pDOCKQ2_value = filtered_ipsae['pDockQ2'].values[0]
        LIS_value = filtered_ipsae['LIS'].values[0]

        # Extract the required information
        iptm = data.get('iptm', None)
        chain_pair_pae_min = max(max(data.get('chain_pair_pae_min', None)))
        ranking_score = data.get('ranking_score', None)

        # Append the information to the summary_data dictionary
        summary_data["description"].append(description)
        summary_data["iPTM"].append(iptm)
        summary_data["chain_pair_pae_min"].append(chain_pair_pae_min)
        summary_data["ranking_score"].append(ranking_score)
        summary_data["plddt_binder_AF3"].append(round(plddt_af3,2))
        summary_data["ipSAE_AF3"].append(round(ipSAE_value,2))
        summary_data["pDOCKQ_AF3"].append(round(pDOCKQ_value,2))
        summary_data["pDOCKQ2_AF3"].append(round(pDOCKQ2_value,2))
        summary_data["LIS_AF3"].append(round(LIS_value,2))


    # Create a DataFrame from the summary_data dictionary
    summary_df = pd.DataFrame(summary_data)

    # Merge with the existing CSV if it exists


    # Define the output CSV file path
    output_csv_path = f"{dir}/Scoring_Stats_AF3.csv"
    if os.path.exists(output_csv_path):
        print(f"Appending to existing summary file at {output_csv_path}")
        df=pd.read_csv(output_csv_path)
        summary_df=pd.concat([df, summary_df], ignore_index=True)
        merged_df = pd.merge(scoring_df, summary_df, left_on='description', right_on='description', how='inner')
        #remove duplicates
        rows=[]
        merged_fin=pd.DataFrame(columns=list(merged_df.columns))
        for r in range(0,merged_df.shape[0]):
            row=merged_df.iloc[r,:]
            if not row[0] in rows:
                rows.append(row[0])
                merged_fin.loc[len(merged_fin)] = row
        merged_fin.to_csv(output_csv_path, index=False)
    else:
        print(f"Creating new summary file at {output_csv_path}")
        merged_df = pd.merge(scoring_df, summary_df, left_on='description', right_on='description', how='inner')
        #remove duplicates
        rows=[]
        merged_fin=pd.DataFrame(columns=list(merged_df.columns))
        for r in range(0,merged_df.shape[0]):
            row=merged_df.iloc[r,:]
            if not row[0] in rows:
                rows.append(row[0])
                merged_fin.loc[len(merged_fin)] = row
        merged_fin.to_csv(output_csv_path, index=False)

    print(f"Summary information saved to {output_csv_path}")
    print("Preparing RMSD script...")



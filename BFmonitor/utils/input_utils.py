import os
import json
import base64
import io
import re
import subprocess
from dash_molstar.utils import molstar_helper
from dash_molstar.utils.representations import Representation
from Bio import PDB


###############################################
########## JSON MANAGEMENT FUNCTIONS ##########
###############################################


def json_uploader(contents,dir,process):
    '''
    Function to decode the uploaded json file and save it in the working directory
    -----------
    Inputs:
    - contents --> Json file uploaded by the user
    - dir --> Working directory where the file will be saved
    -----------
    Outputs:
    - json_string --> String with the content of the json file
    ----------- 
    '''
    _, content_string = contents[0].split(',')
    decoded = base64.b64decode(content_string)
    if process == "Protein_Design":
        data = json.load(io.StringIO(decoded.decode('utf-8')))
        json_string=json.dumps(data, indent=2)
        with open(str(dir)+"/input_binderflow_Design.json","w") as f:
            f.write(json_string)
        return json_string
    elif process == "Partial_Diffusion":
        data = json.load(io.StringIO(decoded.decode('utf-8')))
        json_string=json.dumps(data, indent=2)
        with open(str(dir)+"/input_binderflow_PD.json","w") as f:
            f.write(json_string)  
        return json_string

    data = json.load(io.StringIO(decoded.decode('utf-8')))
    json_string=json.dumps(data, indent=2)
    with open(str(dir)+"/input_binderflow.json","w") as f:
        f.write(json_string)
    return json_string

def get_ckp_paths():
    '''
    Function to get the checkpoint paths from the config file
    -----------
    Outputs:
    - checkpoint_files --> List of checkpoint file names
    - ckp_abspaths --> Dict that links ckp files and its absolute paths 
    -----------
    '''
        # Get the checkpoint path from the config file
    SCRIPT_DIR = os.path.dirname(__file__).split("BFmonitor/utils")[0]
    config_file = os.path.join(SCRIPT_DIR, 'config.sh')
    checkpoint_dir = ''
    with open(config_file, 'r') as f:
        for line in f:
            if line.startswith('RFD_PATH'):
                checkpoint_dir = line.split('=')[1].strip().strip('"') + '/models/'
                break
    checkpoint_files = [f for f in os.listdir(checkpoint_dir) if f.endswith('.pt')]
    ckp_abspaths = {}
    for ckp in checkpoint_files:
        ckp_abspaths[ckp] = os.path.join(checkpoint_dir, ckp)
    return checkpoint_files, ckp_abspaths

def generate_json_from_flags(process, size, dir, template, max_threads, residues_molstar, rfd_designs, pmpnn_seqs, pmpnn_cycles, noise_steps, noise_scale, ckp, hits, SCRIPT_DIR):
    '''
    Function to generate the json file for BFmonitor from the input flags
    
    -----------
    Inputs:
    - process --> Generation process (initial_generation, partial_diffusion, sequence_diversity)
    - size --> Size range for the binder
    - dir --> Working directory where the temp pdb file is located
    - template --> Template type for RFD (antibody, protein, custom)
    - max_threads --> Maximum number of threads to use
    - residues_molstar --> Residues selected as hotspots in molstar
    - rfd_designs --> Number of designs to generate
    - pmpnn_seqs --> Number of sequences to generate with pMPNN
    - pmpnn_cycles --> Number of relaxation cycles with pMPNN
    - noise_steps --> Number of noise steps for RFD (Partial Diffusion only)
    - noise_scale --> Noise scale for RFD (Partial Diffusion only)
    - ckp --> Checkpoint file name
    - hits --> Number of hits to generate
    - SCRIPT_DIR --> Directory where the BinderFlow code is located
    -----------
    Outputs:
    - json_file --> JSON file content as a string
    -----------    
    '''
    # Get the contigs map and pdb location
    rfd_ctgs, pdb_loc = get_contigs_ref_pdbLOC(process, size, dir)

    # Get the checkpoint absolute paths 
    _, ckp_abspaths = get_ckp_paths()
    ckp = ckp_abspaths[ckp] if ckp != None else ckp_abspaths['Complex_base_ckpt.pt']
    

    # Get the common JSON features
    json_dict = {}
    json_dict["input"] = pdb_loc
    json_dict["template"] = template if template != 'None' else pdb_loc 
    json_dict["max_threads"] = max_threads if max_threads != None else 1
    json_dict["pmp_nseqs"] = pmpnn_seqs if pmpnn_seqs != None else 1
    json_dict["pmp_relax_cycles"] = pmpnn_cycles if pmpnn_cycles != None else 0
    json_dict["noise_steps"] = noise_steps if noise_steps != None else 20
    json_dict["noise_scale"] = noise_scale if noise_scale != None else 1.0
    json_dict["checkpoint"] = ckp
    json_dict["node"] = '' 
    json_dict["residues"] = "None"
    json_dict["hits_number"] = hits if hits != None else 96
    json_dict["rfd_ndesigns"] = rfd_designs if rfd_designs != None else 10
    json_dict["rfd_contigs"] = str(rfd_ctgs)
    # Get specific JSON features
    if process == "Protein_Design":

        # Check if there are residues selected as hotspots
        if str(residues_molstar).startswith("[B"):
            used_res=str(residues_molstar).strip()
        else:
            used_res= " "
        json_dict["rfd_hotspots"] = used_res
        json_dict["partial_diff"] = "False"
        json_dict["sequence_diversity"] = "False"

    elif process == "Partial_Diffusion":
        if str(rfd_ctgs) == "":
            return ("No valid PDB for Partial Diffusion",None)
        json_dict["rfd_hotspots"] = ''
        json_dict["partial_diff"] = "True"
        json_dict["sequence_diversity"] = "False"
        json_file=json.dumps(json_dict,indent=2)

    elif process == "Sequence_Diversity":
        json_dict["rfd_hotspots"] = ''
        json_dict["partial_diff"] = "False"
        json_dict["sequence_diversity"] = "True"
        json_file=json.dumps(json_dict,indent=2)

    json_file=json.dumps(json_dict,indent=2)
    with open(str(dir)+"/bfmonitor_input.json","w") as f:
        json.dump(json_dict, f, indent=2)
    return json_file


def get_contigs_ref_pdbLOC(process, size, dir):
    '''
    Function to get the contigs map for the generation process and the location of the temp pdb file 
    -----------
    Inputs:

    - process --> Generation process (initial_generation, partial_diffusion, sequence_diversity)
    - size --> Size range for the binder
    - dir --> Working directory where the temp pdb file is located
    -----------
    Outputs:
    - rfd_ctgs --> Contigs map for the RFD input
    - pdb_loc --> Location of the temp pdb file
    -----------'''
    # The input pdb file must be in the input folder, get absolute path
    pdb_path = os.path.abspath(f"{dir}/input/{process}_input.pdb")
    SCRIPT_DIR=os.path.dirname(os.path.abspath(__file__)).split("BFmonitor/utils")[0]

    #print(pdb_input)
    if process == "Protein_Design":
        command = 'python3 ' + f'{SCRIPT_DIR}/binderflow/scripts/contigs_map_getter.py -i ' + pdb_path
        contigs=subprocess.run(command, shell=True, capture_output=True, text=True).stdout
        rfd_contigs= f"[ {size[0]}-{size[1]}/0 {contigs[2:-3]} ]"

    elif process == "Partial_Diffusion":
        command = 'python3 ' + f'{SCRIPT_DIR}/binderflow/scripts/contigs_map_getter.py -i ' + pdb_path + ' --partial_diff True'
        rfd_c = subprocess.run(command, shell=True, capture_output=True, text=True).stdout
        rfd_contigs=rfd_c[:-1]
    
    else:
        rfd_contigs = ""
    return rfd_contigs, pdb_path

###############################################
########## INPUT STRUCTURE FUNCTIONS ##########
###############################################


def corp_PDB(residues, dir):
    '''
    Function to corp the uploaded pdb file in Protein_Design process
    -----------
    Inputs:
    - dir --> Working directory where the file will be saved
    - residues --> The sequence between those residues will be eliminated
    -----------
    Outputs:
    - results --> Corped pdb file content
    -----------
    '''

    pdb_file=f"{dir}/input/Protein_Design_input.pdb"
    #print(pdb_file)
    with open(pdb_file,"r") as f:
        readed=f.read()
    with open(pdb_file,"r") as lin:
        readed_l=lin.readlines()

    lines=[]
    res=residues[1:-1]

    lines.append(re.findall(fr'ATOM.*{res.split(",")[0]}.*',readed)[0])
    lines.append(re.findall(fr'ATOM.*{res.split(",")[-1]}.*',readed)[-1])
    print(lines)
    new_doc=""; ini=1; fin=0
    for l in readed_l:
        #print(str(l))
        if str(l)[:-1] != lines[0] and ini==1:
            new_doc = new_doc + str(l)
        else:
            ini=0
        if str(l)[:-1] == lines[1]:
            fin=1
        if fin == 1 and str(l)[:-1] != lines[1]:
            new_doc = new_doc + str(l)
        
    #print(new_doc)
    with open(pdb_file,"w") as n:
        n.write(new_doc)

    return 0


def decoded_pdbFILE(pdb_file):
    '''
    Function to decode the uploaded pdb file. This is needed since Dash uploads the file in base64 format
    -----------
    Inputs:
    - dir --> Working directory where the file will be saved
    - process --> Generation process (initial_generation, partial_diffusion, sequence_diversity)
    - pdb_file --> PDB file uploaded by the user
    -----------
    Outputs:
    - decoded --> Decoded pdb file content
    -----------
    '''
    _, content_string = pdb_file[0].split(',')
    decoded = base64.b64decode(content_string).decode("utf-8")
    return decoded

def load_input_pdb(filename, process, pdb_file, dir,surface):
        '''
        Load input PDB file for the specified process.
        -----------
        Inputs:
        - filename --> Name of the uploaded PDB file, as Dash provides it
        - process --> Generation process (initial_generation, partial_diffusion, sequence_diversity)
        - pdb_file --> PDB file uploaded by the user
        - dir --> Working directory where the file will be saved
        - surface --> Representation of the surface for molstar visualization
        -----------
        Outputs:
        - data --> Parsed molecule data for visualization
        - error --> Error message if any issues occur
        -----------
        '''
        try:
            if str(filename[0]).endswith(".pdb"):
                decoded = decoded_pdbFILE(pdb_file)

                # Create the input directory if it doesn't exist
                if not os.path.exists(f"{dir}/input"):
                    os.makedirs(f"{dir}/input")
                # Write down the pdb
                pdb_path = f'{dir}/input/{process}_input.pdb'
                with open(pdb_path,"w") as tmp:
                    tmp.write(decoded)
                            

                #Edit pdb if Protein Design
                if process == "Protein_Design":

                    process_pdb_B_1000(pdb_path, dir, process, filename)
                    
                    # Parse the molecule using molstar_helper
                    chainB = molstar_helper.get_targets(chain="B")
                    chainB_representation = Representation(type=surface, color="hydrophobicity")
                    component_B = molstar_helper.create_component("Chain B", chainB, chainB_representation)
                    components=[component_B]

                    # Load the molecule data
                    data = molstar_helper.parse_molecule(pdb_path, fmt='pdb',component=components)

                    return [data]
                elif process == "Partial_Diffusion" or process == "Sequence_Diversity":
                    partial_diff_numeration(pdb_path)
                    components=partial_diff_components(pdb_path, surface)

                    data = molstar_helper.parse_molecule(pdb_path, fmt='pdb',component=components)
                    return [data]
        except:
            return 'There was an error loading the PDB file. Please, check that the file is correct.'

def load_template_pdb(process, pdb_file, dir,):
    '''
    Function to load the template pdb file for the specified process.
    -----------
    Inputs:
    - process --> Generation process (initial_generation, partial_diffusion, sequence_diversity)
    - dir --> Working directory 
    - pdb_file --> PDB file uploaded by the user
    -----------
    Outputs:
    - pdb_path --> Path to the saved template pdb file
    - Saves the template pdb file in the working directory
    -----------
    '''
    decoded = decoded_pdbFILE(pdb_file)

    # Create the input directory if it doesn't exist
    if not os.path.exists(f"{dir}/input"):
        os.makedirs(f"{dir}/input")
    # Write down the pdb
    pdb_path = f'{dir}/input/{process}_template.pdb'
    with open(pdb_path,"w") as tmp:
        tmp.write(decoded)
    return pdb_path

def partial_diff_components(pdb_file, surface):
    '''
    Function to create the components for the partial diffusion and sequence diversity input visualization
    -----------
    Inputs:
    - pdb_file --> Path to the pdb file
    - surface --> Representation of the surface for molstar visualization
    -----------
    Outputs:
    - components --> List of components for the molstar visualization
    -----------
    '''
    components=[]
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('input_structure', pdb_file)
    chains = ""
    for model in structure:
        for ch in model:
            ch = ch.id
            if ch == "B":
                representation = Representation(type=surface, color="hydrophobicity")
            else:
                representation = Representation(type="cartoon", color="hydrophobicity")

            component = molstar_helper.create_component(f"Chain {ch}", molstar_helper.get_targets(chain=ch), representation)
            components.append(component)
    return components

###########################################
########## NUMERATION FUNCTIONS ###########
###########################################

def process_pdb_B_1000(path, dir, process, filename, template=False):
    '''
    Function to process the pdb file to change chain B residue numbers +1000 using Biopython
    It also removes HETATM lines
    -----------
    Inputs:
    - path --> Path to the pdb file
    - dir --> Working directory where the processed file will be saved
    - process --> Generation process (initial_generation, partial_diffusion, sequence_diversity)
    - filename --> Name of the uploaded PDB file, as Dash provides it
    -----------
    Outputs:
    
    - Saves a new pdb file in the working directory with chain B residue numbers +1000
    -----------
    '''

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('input_structure', path)
    for model in structure:
        for chain in model:
            # Create a list to hold residues to be removed (HETATM)
            residues_to_remove = []
            if len(chain.id) == 1:
                if chain.id !='B':
                    # Change the id of the chain A to chain B
                    chain.id = 'B'
            else: 
                pass
            if chain.id == 'B':
                for residue in chain:
                    if residue.id[0] == " ":
                        residue.id = (residue.id[0], residue.id[1] + 1000, residue.id[2])
                    else:
                        residues_to_remove.append(residue)
                # Remove HETATM residues
                for residue in residues_to_remove:
                    chain.detach_child(residue.id)
    io = PDB.PDBIO()
    io.set_structure(structure)
    # Create the input directory if it doesn't exist
    if template == False:
        if not os.path.exists(f"{dir}/input"):
            os.makedirs(f"{dir}/input")
        io.save(f"{dir}/input/{process}_input.pdb")
    else:
        if not os.path.exists(f"{dir}/input"):
            os.makedirs(f"{dir}/input")
        io.save(f"{dir}/input/{process}_template.pdb")


def partial_diff_numeration(pdb_file):
    '''
    Function to renumber chain B residues to continue from chain A in a partial diffusion pdb file
    -----------
    Inputs:
    pdb_file --> Path to the pdb file
    -----------
    Outputs:
    Saves a new pdb file in the working directory with chain B residue numbers continuing from chain A
    '''
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    chain_a = structure[0]['A']
    chain_b = structure[0]['B']

    # Get the sequence of chain A
    seq_a = ''.join([residue.get_resname() for residue in chain_a])
    
    # Get the length of chain A
    length_a = len(seq_a)/3 #Dividing by 3 to get the number of residues (each residue is represented by 3 characters)

    # Renumber chain B to match the sequence of chain A
    for i, residue in enumerate(chain_b):
        residue.id = (' ', i + length_a + 1, ' ')

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

##############################################
########## RUN BINDERFLOW FUNCTIONS ##########
##############################################

def run_binderflow(SCRIPT_DIR, dir, pd_dict_str, process):
    '''
    Function to run BinderFlow in the background using nohup
    -----------
    Inputs:
    - SCRIPT_DIR --> Directory where the BinderFlow code is located
    - dir --> Working directory where the json file is located
    -----------
    Outputs:
    - Executes BinderFlow in the background
    -----------
    '''
    # get the realpath of dir
    dir = os.path.realpath(dir)
    if process != "Partial_Diffusion":
        cmd = (
        f"nohup {SCRIPT_DIR}/binderflow.sh --json {dir}/bfmonitor_input.json --metrics No_data"
        f"> {dir}/project_binderflow_Design.log 2>&1 &"
        )
    if process == "Partial_Diffusion":
        pd_dict_string=str(pd_dict_str).replace(" ","")
        cmd = (
        f"nohup {SCRIPT_DIR}/binderflow.sh --json {dir}/bfmonitor_input.json --metrics {pd_dict_string}"
        f"> {dir}/project_binderflow_Design.log 2>&1 &"
        )
    subprocess.run(cmd, shell=True, cwd=dir)
    return f"Binderflow Executed\nResults available in {dir}"


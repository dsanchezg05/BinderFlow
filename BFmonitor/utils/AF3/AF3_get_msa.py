"""Extracts protein/RNA MSA from the AlphaFold 3 JSON and saves them as individual a3m files."""
import os
from alphafold3.common import folding_input
import argparse


parser=argparse.ArgumentParser(description="Extracts MSA from AlphaFold 3 JSON input.")
parser.add_argument('--input_json_path', type=str, required=True, help='Path to the input JSON file containing AlphaFold 3 data.')
parser.add_argument('--output_dir', type=str, required=True, help='Directory to save the extracted MSA files.')
args = parser.parse_args()



INPUT_JSON_PATH = args.input_json_path

with open(INPUT_JSON_PATH, 'rt') as f:
  af_json = f.read()

af_input = folding_input.Input.from_json(af_json)

# Get the input directory name
OUTPUT_DIR = args.output_dir
for protein_chain in af_input.protein_chains:
  if protein_chain.unpaired_msa is not None:
    with open(os.path.join(OUTPUT_DIR, f'Unpaired_msa.a3m'), 'wt') as f:
      f.write(protein_chain.unpaired_msa)

  if protein_chain.paired_msa is not None:
    with open(os.path.join(OUTPUT_DIR, f'Paired_msa.a3m'), 'wt') as f:
      f.write(protein_chain.paired_msa)

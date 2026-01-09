#!/usr/bin/env python3 


import pandas as pd 
import argparse
import subprocess

'''
Read the csv with the scorings and the sc files, get all the info, and end the run once enough binders are generated

Input:

--number: Max number of hits that want to be generated

Output:

campaign_done: An empty file that marks that all wanted binders have been generated and stops the run
'''

# Define number of designs desired
parser=argparse.ArgumentParser()
parser.add_argument('--number', '-n' , type=int, default=100, help='Maximum number of successful binder desired')
parser.add_argument('--partial_diff', '-pd' , type=str, help='Partial Diffusion')
parser.add_argument('--plddt', '-pt' , type=float, help='Original model PLDDT')
parser.add_argument('--pae', '-pi' , type=float, help='Original model PAE')
parser.add_argument('--ipsae', '-ip' , type=float, help='Original model ipSAE')
args=parser.parse_args()

number=args.number
prt_diff=args.partial_diff
pt=args.plddt
pi=args.pae
ip=args.ipsae

print(prt_diff)

# define Thresholds
if str(prt_diff) == "True":
	print("using modified thersholds")
	pae_interaction_thres=pi
	plddt_binder_thres=pt
	ipSAE=ip
else: 
	pae_interaction_thres=15
	plddt_binder_thres=80
	ipSAE=0.6

print(pae_interaction_thres)
print(plddt_binder_thres)
print(ipSAE)
# CUTRE_thres=10
# unsat_hbonds_thres=4
# hbond_thres=3
# binder_surface_hyd_thres=0.35
# shape_complementarity_thres=0.55
# dSASA_thres=1000 # is a percentage better ?

# Read the csv file
try:
	df=pd.read_csv('Scoring_Stats.csv')

	# Load the variable 

	hits_number=0

	# Check the conditions

	hits_number = (
	# 	(df['CUTRE'] <= CUTRE_thres) &
	# 	(df['interface_unsat_hbonds'] <= unsat_hbonds_thres) &
	# 	(df['interface_hbonds'] >= hbond_thres) &
	# 	(df['binder_surf_hyd'] <= binder_surface_hyd_thres) &
	# 	(df['Shape_complementarity'] >= shape_complementarity_thres) &
	# 	(df['dSASA'] >= dSASA_thres)
		(df["ipSAE"] >= ipSAE) &
		(df['pae_interaction'] <= pae_interaction_thres) & 
		(df['plddt_binder'] >= plddt_binder_thres)).sum()
	
	command='touch campaign_done'
	print("Hits: {0}".format(hits_number))
	if hits_number >= number:
		print('done')
		subprocess.run(command, shell=True)
	else:
		print('continue') 

except FileNotFoundError:
	print('No Run has been completed yet')

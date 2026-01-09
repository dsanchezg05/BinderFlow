import os
import sys
import pandas as pd
import numpy
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--directory", "-d", help='Working directory')
args = parser.parse_args()

dir = args.directory

merged_df=pd.read_csv(f"{dir}/hits/hits_clustered.csv")

#cols= args.col
val="cluster,plddt_binder_AF3".split(",")
#print(val)
sorted_df=merged_df.sort_values(by=val,ascending=False).reset_index(drop=True)
try:
    sorted_df.drop(["Unnamed: 0"],axis=1, inplace=True)
except:
    pass
unique_cluster=sorted_df["cluster"].unique()
print(unique_cluster)
R=[]
for cl in unique_cluster:
    #print(sorted_df[sorted_df["cluster"] == cl])
    #print("\n")
    for row in range(0,sorted_df[sorted_df["cluster"] == cl].shape[0]):
        if sorted_df[sorted_df["cluster"] == cl].iloc[row,:]["ipSAE_AF3"] > 0.6:
            R.append(sorted_df[sorted_df["cluster"] == cl].iloc[row,:])
        break

rep_df=pd.DataFrame(columns=list(sorted_df.columns))
for r in R:
    rep_df.loc[len(rep_df)] = r

rep_df.sort_values(by=["cluster"]).reset_index(drop=True).to_csv(f'{dir}/hits/Representative_models.csv')
print(rep_df.sort_values(by=["cluster"]).reset_index(drop=True)["description"])
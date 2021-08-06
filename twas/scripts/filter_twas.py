import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--twas_folder", help="Folder with TWAS results")
parser.add_argument("--models_folder", help="Folder where the prediction models are saved")
parser.add_argument("--models",help="Comma separated list of models used for TWAS analysis")
parser.add_argument("--out",help="Output file")
args = parser.parse_args()


###Load TWAS results###
data=[]
for file in os.listdir(args.twas_folder):
    temp=pd.read_csv(os.path.join(args.twas_folder,file))
    temp["model"]=file.replace(".asso.csv","")
    data.append(temp)

data=pd.concat(data)

###Get number of gene-model pairs##
models=args.models.split(",")
npairs=0
for model in models:
    npairs+=pd.read_csv(f"{args.models_folder}/{model}/output/{model}.cov",sep="\t",usecols=["GENE"]).drop_duplicates().shape[0]

###Filter results###
threshold=0.05/npairs
data=data[data.pvalue<=threshold]

###Write results
outdir=os.path.abspath(os.path.join(args.out, os.pardir))
if not os.path.exists(outdir):
    os.makedirs(outdir)

data.to_csv(args.out,index=False)


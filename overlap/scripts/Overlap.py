import datatable as dt
import pandas as pd
import argparse as argp 

parser = argp.ArgumentParser()
parser.add_argument("--gwas", help="GWAS summary statistics file")
parser.add_argument("--ld",help="plink output file")
parser.add_argument("--eqtls",help="File with eQTL information")
parser.add_argument("--out",help="output file")
args = parser.parse_args()

##Read Files 
gwas=pd.read_table(args.gwas,usecols=["rsid","pval","beta","se","ref","alt"])
ld=dt.fread(cmd=f"sed -r 's/\\s+/\\t/g' {args.ld} | cut -f5,6,7").to_pandas()
eqtls=pd.read_table(args.eqtls)

significant=gwas[gwas.pval<5*10**-8].rsid.tolist()
ld.columns=["chr","pos","rsid"]
ld.chr="chr"+ld.chr.astype(str).str[:]
ld.drop_duplicates(inplace=True)
ld["Category"]=["Lead" if snp in significant else "LD" for snp in ld.rsid]

##extract gwas info 
ld=pd.merge(ld,gwas,how="left",on="rsid")

##Overlap
overlap=pd.merge(ld,eqtls,on=["chr","pos"]).drop_duplicates()
overlap.to_csv(args.out,index=False)

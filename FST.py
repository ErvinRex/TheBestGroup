#install & import all necessary packages
pip install plotnine

conda install -c conda-forge scikit-allel

import numpy as np
import pandas as pd
import allel; print('scikit-allel', allel.__version__)
from plotnine import ggplot, aes, geom_line, geom_bar

#Reading the vcf file and creating a dateframe from it
df = allel.vcf_to_dataframe("Chr22_Phased.vcf.gz", fields=None, exclude_fields=None, types=None, numbers=None,
                       alt_number=3, fills=None, region=None,
                       tabix='tabix', transformers=None, buffer_size=16384, chunk_length=65536, log=None)

#Reading the vcf file and creating an array
callset = allel.read_vcf("Chr22_Phased.vcf.gz")

#Extrapolating SNPs positional array
pos = callset["variants/POS"]

#convert samples to df and only keep sample name, population code, and superpop code 
samples = pd.read_excel("Samples.xlsx",  usecols=['Sample name', 'Population code', 'Superpopulation code']) #convert excel sheet to df
samples.rename(columns={'Sample name':'sample_name'}, inplace=True) #rename Sample name to sample_name
samples.rename(columns={'Population code':'Pop_code'}, inplace=True) #rename Population code to Pop_code
samples.rename(columns={'Superpopulation code':'Superpop_code'}, inplace=True) #rename Superpopulation code to Superpop_code

#get dataframe with samples of 5 population of interest
samples_5pop = (samples.loc[samples['Pop_code'].isin(["LWK", "CLM", "CHS", "TSI", "STU"])])

#Get list of sample names of interest 
Sample_names = samples_5pop['sample_name'].tolist()

#get data for only samples of interest
callset_5pop = allel.read_vcf('Chr22_Phased.vcf.gz', samples=Sample_names)

#get genotype data into array
gt_5pop = allel.GenotypeArray(callset_5pop['calldata/GT'])

#Identifying Population_code to work with
Africa = 'LWK'
Americas = 'CLM'
East_Asian = 'CHS'
Europe = 'TSI'
South_Asian = 'STU'

n_samples_Africa = np.count_nonzero(df_samples[df_samples['Population code'] == Africa])
n_samples_Americas = np.count_nonzero(df_samples[df_samples['Population code'] == Americas])
n_samples_East_Asian = np.count_nonzero(df_samples[df_samples['Population code'] == East_Asian])
n_samples_Europe = np.count_nonzero(df_samples[df_samples['Population code'] == Europe])
n_samples_South_Asian = np.count_nonzero(df_samples[df_samples['Population code'] == South_Asian])

#dictionary mapping population names to sample indices
subpops = {
    Africa: df_samples[df_samples['Population code'] == Africa].index,
    Americas: df_samples[df_samples['Population code'] == Americas].index,
    East_Asian: df_samples[df_samples['Population code'] == East_Asian].index,
    Europe: df_samples[df_samples['Population code'] == Europe].index,
    South_Asian: df_samples[df_samples['Population code'] == South_Asian].index,
}

# Producing allele frequency arrays for populations
def Allele_Frequncy(filename, vcf, Code):
    samples = pd.read_excel(filename,  usecols=['Sample name', 'Population code', 'Superpopulation code']) #convert excel sheet to df
    samples_Code = samples[samples['Population code'] == Code] #get a dataframe with only samples from pop code
    names = samples_Code['Sample name'].tolist() #get list of samples from population of interest
    callset = allel.read_vcf(vcf, samples= names) #get the callset for samples of interest
    gt = allel.GenotypeArray(callset['calldata/GT']) #get the genotype array for samples of interest 
    ac = gt.count_alleles() #get the allele count 
    ac_df = pd.DataFrame(ac) #convert allele count to dataframe
    ac_df.columns = ("Ref_AC", "Alt_AC") #change the column names 
    poslist = callset['variants/POS'] #get the list of positons 
    return ac

acs = Allele_Frequncy('Samples.xlsx', 'Chr22_Phased.vcf.gz', 'LWK')
acs2= Allele_Frequncy('Samples.xlsx', 'Chr22_Phased.vcf.gz', 'CLM')

#Creating positional list for plotting
pos1 = pos.tolist()

#Creating windows for plotting
windows1 = []
for x in range(0, (len(pos1)-10000)):
    a = (pos1[x], pos1[x+10000])
    windows1.append(a)

#Hudson Fst with sliding window withs specified 
D1, windows, counts = allel.windowed_hudson_fst(pos, acs, acs2, windows = windows1)

#get the x values  which is the mean of each window
x1 = windows.mean(axis=1)

#make a dataframe of the positons and FST values 
HFs= pd.DataFrame(D1)
HFs["position"] = x1
HFs.columns = ["FST", "Position"]
HFs['row_ID'] = HFs.reset_index().index

#make empty rows 0 
HFs['FST'] = HFs['FST'].replace(np.nan, 0) #replace NAN to 0

#plotting the dataframe against poistion 
(
    ggplot(HFs) #what dataframe to use 
     + aes(x="Position", y = "FST") + geom_line()  # specify columns that are the x and y axis and a line graph 
)
#plottinng against the index 
HFs['row_ID'] = HFs.reset_index().index
(
    ggplot(HFs) #what dataframe to use 
     + aes(x="row_ID", y = "FST") + geom_line()  # specify columns that are the x and y axis and a line graph 



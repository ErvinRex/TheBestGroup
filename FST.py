#import all necessary packages
import matplotlib
import numpy as np
import scipy
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('notebook')
import pandas
import allel; print('scikit-allel', allel.__version__)

conda install -c conda-forge scikit-allel
pip install openpyxl

#Reading the vcf file and creating a dateframe from it
df = allel.vcf_to_dataframe("Chr22_Phased.vcf.gz", fields=None, exclude_fields=None, types=None, numbers=None,
                       alt_number=3, fills=None, region=None,
                       tabix='tabix', transformers=None, buffer_size=16384, chunk_length=65536, log=None)

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

#Taking a look at the sample data
df_samples = pandas.read_csv('Samples2.csv', sep=',')

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

# allele counts - this is where I have problems
acs = gt_5pop.count_alleles_subpops(subpops)
acs

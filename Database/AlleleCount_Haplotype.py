''''Code to get alllel count database from phased vcf
Link to Phased VCF (519930289) :http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/ '''

#import statements 
import allel
import pandas as pd

#function to get the allele count dataframe for a particular population 
def Allele_Count(filename, vcf, Code):
    samples = pd.read_excel(filename,  usecols=['Sample name', 'Population code', 'Superpopulation code']) #convert excel sheet to df
    samples_Code = samples[samples['Population code'] == Code] #get a dataframe with only samples from pop code
    names = samples_Code['Sample name'].tolist() #get list of samples from population of interest
    callset = allel.read_vcf(vcf, samples= names) #get the callset for samples of interest
    gt = allel.GenotypeArray(callset['calldata/GT']) #get the genotype array for samples of interest 
    ac = gt.count_alleles() #get the allele count 
    ac_df = pd.DataFrame(ac)
    ac_df.columns = ["Ref_AC", "Alt_AC"]
    return ac_df

#get the dataframe for each population allele count
AC_LWK = Allele_Count("Samples.xlsx", "Chr22_Phased.vcf.gz", "LWK")
AC_CLM = Allele_Count("Samples.xlsx", "Chr22_Phased.vcf.gz", "CLM")
AC_CHS = Allele_Count("Samples.xlsx", "Chr22_Phased.vcf.gz", "CHS")
AC_TSI = Allele_Count("Samples.xlsx", "Chr22_Phased.vcf.gz", "TSI")
AC_STU = Allele_Count("Samples.xlsx", "Chr22_Phased.vcf.gz", "STU")

#code to convert  VCF to dataframe using allel.vcf_to_dataframe()
df = allel.vcf_to_dataframe('Chr22_Phased.vcf.gz', fields='*', numbers={'ALT': 1}) #* extracts all the info, SNP have only 1 alt

#create dataframe with columns of interest
df_AC= df[['POS']] #position
df_AC['is_snp']= df[['is_snp']] # Column to determine if SNP

#Get the Allele counts for all the populations into the dataframe 
df_AC["Ref_AC_LWK"] = AC_LWK[["Ref_AC"]]
df_AC["Alt_AC_LWK"] = AC_LWK[["Alt_AC"]]

df_AC["Ref_AC_CLM"] = AC_CLM[["Ref_AC"]]
df_AC["Alt_AC_CLM"] = AC_CLM[["Alt_AC"]]

df_AC["Ref_AC_CHS"] = AC_CHS[["Ref_AC"]]
df_AC["Alt_AC_CHS"] = AC_CHS[["Alt_AC"]]

df_AC["Ref_AC_TSI"] = AC_TSI[["Ref_AC"]]
df_AC["Alt_AC_TSI"] = AC_TSI[["Alt_AC"]]

df_AC["Ref_AC_STU"] = AC_STU[["Ref_AC"]]
df_AC["Alt_AC_STU"] = AC_STU[["Alt_AC"]]

#isolate SNP to make sure it is true 
df_ACSNP = df_AC[df_AC["is_snp"]== True] #only keeps rows where is_SNP  value is True

#remove is_snp column 
df_ACSNP.drop(columns=["is_snp"])

df_ACSNP['row_ID'] = df_ACSNP.reset_index().index

#Save it as zipped file
df_ACSNP.to_csv("AllelCount.csv.gz", 
           index=False, 
           compression="gzip")



'''Code to get the Haplotype array as csv from the phased VCF'''

#function to get haplotype for a population code and output it as a csv file 

def Haplotype_Dataframe(code, filename):
    samples = pd.read_excel("Samples.xlsx",  usecols=['Sample name', 'Population code', 'Superpopulation code']) #convert excel sheet to df
    samples = samples[samples['Population code'] == code] #get df for specific population code
    names = samples['Sample name'].tolist() #get list of sample names to list 
    callset = allel.read_vcf('Chr22_Phased.vcf.gz', samples=names) #get data from vcf file of samples of interest
    gt = allel.GenotypeArray(callset['calldata/GT']) #get the genotype array
    Hap = gt.to_haplotypes() #get the haplotype array 
    df = pd.DataFrame(Hap) #convert haplotype array to dataframe 
    pos = callset['variants/POS'] #get position array from the vcf file 
    pos = pos.tolist() #convert position array to list 
    df['POS'] = pos #add position column to dataframe 
    df1 = pd.read_csv("AllelCount.csv.gz") #convert allele count document to dataframe 
    poslist = df1["POS"].to_list() #get position list from allele count dataframe 
    boolean_series = df.POS.isin(poslist) #create boolean series of positons in positon list
    filter_df = df[boolean_series] #filter dataframe to only have positions of interest
    filter_df['row_ID'] = filter_df.reset_index().index #add row ID column so that it can be the primary key
    sample_name = callset['samples'] #get array of sample names from vcf
    sample_name = sample_name.tolist() #change array to list 
    col_names= [] #create empty list 
    for x in sample_name: #loop to add sample names in repeat as each sample has 2 columns 
        col_names.append(x)
        col_names.append(x+'_2')
    col_names.append("POS") #add postion column name to ist 
    col_names.append("row_ID") #add row ID column name to list 
    filter_df.columns = col_names #change dataframe column names to sample names 
    filter_df.to_csv(filename,  #convert to csv file 
           index=False, 
           compression="gzip") #ensure file is zipped 
    

#Haplotype to dataframe for each population 
Haplotype_Dataframe("LWK", "Haplotype_LWK.csv.gz")
Haplotype_Dataframe("CLM", "Haplotype_CLM.csv.gz")
Haplotype_Dataframe("CHS", "Haplotype_CHS.csv.gz")
Haplotype_Dataframe("TSI", "Haplotype_TSI.csv.gz")
Haplotype_Dataframe("STU", "Haplotype_STU.csv.gz")


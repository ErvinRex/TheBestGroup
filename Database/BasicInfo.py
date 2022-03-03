''' The code to obtain the basic info table from the phased vcf file and sample information sheet
    Link to Phased VCF (519930289) :http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
    Link to anotated text (1047074970): http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/
    Link to vcf with the rsID (120868377): http://ftp.ensembl.org/pub/release-105/variation/vcf/homo_sapiens/
    Sample info download: https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
    Missing Gene names from Biomart: https://www.ensembl.org/biomart/martview/a87f1bf2f45e1bdbe3511cbe50e82514
 '''

#import statements 
import allel
import pandas as pd
import numpy as np

#code to convert  VCF to dataframe using allel.vcf_to_dataframe()
df = allel.vcf_to_dataframe('Chr22_Phased.vcf.gz', fields='*', numbers={'ALT': 1})
    #fields '*' allows it to extract all the information 
    #SNP only have one alternate 

#create dataframe with columns of interest
df_Final= df[['CHROM']] #Chromosome number
df_Final['Position']= df[['POS']] #Position
df_Final['REF']= df[['REF']] #Reference Allele
df_Final['ALT']= df[['ALT']] #Alternate Allele
df_Final['is_snp']= df[['is_snp']] # Column to determine if  variant is a SNP

#Function to get the Allele frequency for a population 
def Allele_Frequncy(filename, Code): #sample information filename  and population code
    samples = pd.read_excel(filename,  usecols=['Sample name', 'Population code', 'Superpopulation code']) #convert excel sheet to df
    samples_Code = samples[samples['Population code'] == Code] #get a dataframe with only samples from pop code
    names = samples_Code['Sample name'].tolist() #get list of samples from population of interest
    callset = allel.read_vcf('Chr22_Phased.vcf.gz', samples= names) #get the callset for samples of interest
    gt = allel.GenotypeArray(callset['calldata/GT']) #get the genotype array for samples of interest 
    ac = gt.count_alleles() #get the allele count 
    ac_df = pd.DataFrame(ac) #convert allele count to dataframe
    ac_df.columns = ("Ref_AC", "Alt_AC") #change the column names 
    poslist = callset['variants/POS'] #get the list of positons 
    ac_df["Positions"] = poslist # add a column with the list of positons 
    ac_df["AF"] = ac_df["Alt_AC"]/(ac_df["Ref_AC"]+ac_df["Alt_AC"]) #get the allel frequncy for the alternate allele 
    return ac_df #return dataframe

#Allele frequncy table for each of the 5 populations 
AF_LWK = Allele_Frequncy("Samples.xlsx", "LWK")
AF_CLM = Allele_Frequncy("Samples.xlsx", "CLM")
AF_CHS = Allele_Frequncy("Samples.xlsx", "CHS")
AF_TSI = Allele_Frequncy("Samples.xlsx", "TSI")
AF_STU = Allele_Frequncy("Samples.xlsx", "STU")

#add the allele frequncy columns to final dataframe
df_Final['AF_LWK']= AF_LWK[['AF']]
df_Final['AF_CLM']= AF_CLM[['AF']]
df_Final['AF_CHS']= AF_CHS[['AF']]
df_Final['AF_TSI']= AF_TSI[['AF']]
df_Final['AF_STU']= AF_STU[['AF']]

#Function to get the Allele frequency for a population 
def Genotype_Frequncy(filename, Code): #sample information filename  and population code
    samples = pd.read_excel(filename,  usecols=['Sample name', 'Population code', 'Superpopulation code']) #convert excel sheet to df
    samples_Code = samples[samples['Population code'] == Code] #get a dataframe with only samples from pop code
    names = samples_Code['Sample name'].tolist() #get list of samples from population of interest
    callset = allel.read_vcf('Chr22_Phased.vcf.gz', samples= names) #get the callset for samples of interest
    gt = allel.GenotypeArray(callset['calldata/GT']) #get the genotype array for samples of interest
    i=0 # starts with first position
    #empty list for hommozygous ref, alternate and heterozygous GT count
    hom_ref = []  
    hom_alt = [] 
    het = [] 
    while i < len(gt): #while it less than or equal to length of the array 
        n = gt[i] #start with row i
        #get the counts for each genotype
        hr = n.count_hom_ref()
        ha = n.count_hom_alt() 
        ht = n.count_het()  
        #add counts to list
        hom_ref.append(hr)  
        hom_alt.append(ha)
        het.append(ht)
        i = i+1 #moves to next row 
    #creates dataframe with counts
    df = pd.DataFrame({"hom_ref": hom_ref}) 
    df["hom_alt"] = hom_alt
    df["het"] = het
    #get the hommozygous ref/alt  and heterozyous genotype frequency and add to new dataframe
    df["homAlt_gt"] = df["hom_alt"]/(df["hom_alt"]+df["hom_ref"]+df["het"]) 
    df_fin = df[["homAlt_gt"]]
    df_fin["homRef_gt"] = df["hom_ref"]/(df["hom_alt"]+df["hom_ref"]+df["het"])
    df_fin["het_gt"] = df["het"]/(df["hom_alt"]+df["hom_ref"]+df["het"])
    #get positions and add to final dataframe 
    poslist = callset['variants/POS'] 
    df_fin["Position"]= poslist 
    return df_fin

#Get the genotypes for all 5 populations 
GT_LWK = Genotype_Frequncy("Samples.xlsx", "LWK")
GT_CLM = Genotype_Frequncy("Samples.xlsx", "CLM")
GT_CHS = Genotype_Frequncy("Samples.xlsx", "CHS")
GT_TSI = Genotype_Frequncy("Samples.xlsx", "TSI")
GT_STU = Genotype_Frequncy("Samples.xlsx", "STU")

#add the Genotype frequncy columns to final dataframe
df_Final['GTF_LWK_HomRef']= GT_LWK[['homRef_gt']]
df_Final['GTF_LWK_HomAlt']= GT_LWK[['homAlt_gt']]
df_Final['GTF_LWK_Het']= GT_LWK[['het_gt']]

df_Final['GTF_CLM_HomRef']= GT_CLM[['homRef_gt']]
df_Final['GTF_CLM_HomAlt']= GT_CLM[['homAlt_gt']]
df_Final['GTF_CLM_Het']= GT_CLM[['het_gt']]

df_Final['GTF_CHS_HomRef']= GT_CHS[['homRef_gt']]
df_Final['GTF_CHS_HomAlt']= GT_CHS[['homAlt_gt']]
df_Final['GTF_CHS_Het']= GT_CHS[['het_gt']]

df_Final['GTF_TSI_HomRef']= GT_TSI[['homRef_gt']]
df_Final['GTF_TSI_HomAlt']= GT_TSI[['homAlt_gt']]
df_Final['GTF_TSI_Het']= GT_TSI[['het_gt']]

df_Final['GTF_STU_HomRef']= GT_STU[['homRef_gt']]
df_Final['GTF_STU_HomAlt']= GT_STU[['homAlt_gt']]
df_Final['GTF_STU_Het']= GT_STU[['het_gt']]


#filter dataframe to only contain SNP 
df_FinalSNP = df_Final[df_Final["is_snp"]== True]


#try to save it as zipped file
df_FinalSNP.to_csv("BasicInfo.csv.gz", 
           index=False, 
           compression="gzip")



''' Get the Ancestral Allele and derrived allele frequency'''
#code to convert  VCF with rS and AA to dataframe using allel.vcf_to_dataframe()
df_rs = allel.vcf_to_dataframe('rsID_chr22.vcf.gz', fields='*') 

#get position list from the df_FinalSNP with only SNPs 
poslist = df_FinalSNP["Position"].to_list()

#filtering rsID table to only have the positions with SNPs 
boolean_series = df_rs.POS.isin(poslist)
df_rsSNP = df_rs[boolean_series]

#make new dataframe with only AA, rsID, and position 
df_ID = df_rsSNP[["POS"]]
df_ID["rsID"] = df_rsSNP["ID"]
df_ID["AA"] = df_rsSNP["AA"]

#merge the info from the basic info and the rsID table based on poitison 
BasicInfo = pd.merge(left=df_FinalSNP, right=df_ID, how='left', left_on='Position', right_on='POS')



'''Get Gene names from annotation table '''
 #change csv to dataframe
df_Annot = pd.read_csv("Annotation.csv.gz")

#filtering Annotation table to have only position from list 
boolean_series2 = df_Annot.POS.isin(poslist)
df_AnnotSNP = df_Annot[boolean_series2]

#make a new dataframe with columns of interest 
df_Ann = df_AnnotSNP[["POS"]]
df_Ann["Gene_Name"] = df_AnnotSNP["GENE"]
df_Ann["GeneID"] = df_AnnotSNP['GENEID']

#merge the info from the basic info and the rsID table based on poitison 
BasicInfo = pd.merge(left=BasicInfo, right=df_Ann, how='left', left_on='Position', right_on='POS')

#remove excess columns
BasicInfo = BasicInfo.drop(columns=["is_snp", "POS_x",])

#round Allele Frequency and Genotype frequency columns to 3 decimal places 
BasicInfo = BasicInfo.round({'AF_LWK': 3, 'AF_CLM': 3,'AF_CHS': 3,'AF_TSI': 3,'AF_STU': 3,
'GTF_LWK_HomRef': 3,'GTF_LWK_HomAlt': 3, 'GTF_LWK_Het': 3, 
'GTF_CLM_HomRef': 3,'GTF_CLM_HomAlt': 3, 'GTF_CLM_Het': 3,
'GTF_CHS_HomRef': 3,'GTF_CHS_HomAlt': 3, 'GTF_CHS_Het': 3,
'GTF_TSI_HomRef': 3,'GTF_TSI_HomAlt': 3, 'GTF_TSI_Het': 3,
'GTF_STU_HomRef': 3,'GTF_STU_HomAlt': 3, 'GTF_STU_Het': 3,})



''' Code to calulcate the derrived allele frequency'''


#Get the reference, alternate and AA allele columns into list 
ref = BasicInfo['REF'].tolist()
alt = BasicInfo['ALT'].tolist() 
AA = BasicINfo['AA'].tolist()

#function to convert get the derrived allel frequency from AA and allele frequency
def Derrived_AF(colname):
    AF = df[colname].tolist() #change allel frequency column to list 
    der_AF = [] #empty list 
    for (a, b, c, d) in zip(AA, ref, alt, AF): 
        if a == b: #if the reference and AA are the same then return the alt AF
            der_AF.append(d)
        elif a == c: #if the alt and the AA then return the reference AF
            der_AF.append(1-d)
        else:
            der_AF.append('NaN') #if there is no match then return nothing 
    return der_AF

#get the derrived allele frequency for each population 
DAF_LWK = Derrived_AF('AF_LWK')
DAF_CLM = Derrived_AF('AF_CLM')
DAF_CHS = Derrived_AF('AF_CHS')
DAF_TSI = Derrived_AF('AF_TSI')
DAF_STU = Derrived_AF('AF_STU')

#add derrived allele frequency basic info table 
BasicInfo["DAF_LWK"] = DAF_LWK
BasicInfo["DAF_CLM"] = DAF_CLM
BasicInfo["DAF_CHS"] = DAF_CHS
BasicInfo["DAF_TSI"] = DAF_TSI
BasicInfo["DAF_STU"] = DAF_STU

#remove duplicated rows
BasicInfo = BasicInfo.drop_duplicates()



'''Get missing gene names from biomart and gene alias'''

#convert csv to dataframe
df_mart= pd.read_csv("mart_export.csv") 

#remove positions that don't have a gene name
df_mart = df_mart.dropna(subset=['Gene name']) 
df_mart['row_ID'] = df_mart.reset_index().index #add a row ID column 

#empty list for position, genename, and gene ID
pos = [] 
gname= []
gID = []

#while loop to get the positons, gene name and gene ID into list
#one row for each gene, therefore want to convert to have a row for each position
i=0
while i < len(df_mart): 
    row = df_mart.loc[df_mart['row_ID']==i] #get the first row 
    w = list(range(int(row['Gene start (bp)']), (int(row['Gene end (bp)'])))) #convert start and stop positions to list
    pos.append(w) #add to position list 
    x = list(list(row['Gene name']) * len(w)) #make list of name as long as the position list 
    gname.append(x) #add to gene name list 
    y = list(list(row['Gene stable ID']) * len(w)) #make list of ID as long as the positon list 
    gID.append(y) #add to gene ID list 
    i=i+1 #move to next row 

#flatten lists to get one list
position = sum(pos, [])  
Gene_name = sum(gname, []) 
Gene_ID = sum(gID, [])

#get lists into one dataframe
Gene_info = pd.DataFrame(position, columns = ["position"])
Gene_info['Gene name'] = Gene_name
Gene_info['Gene ID'] = Gene_ID

#get a gene info dataframe with positions from basic info table 
POS = BasicInfo["Position"].to_list() #get positions as list 
boolean_series = Gene_info.position.isin(POS) #create boolean list of position in position list
GeneInfo_df = Gene_info[boolean_series] 

#make GeneInfo_df with BasicInfo
BasicInfo = pd.merge(left=BasicInfo, right=GeneInfo_df, how='left', left_on='Position', right_on='position')

#combine gene annotation names with basic Info names 
#replace nan to 0
BasicInfo['Gene ID'] = BasicInfo['Gene ID'].replace(np.nan, 0)
BasicInfo['GeneID'] = BasicInfo['GeneID'].replace(np.nan, 0)

#Merge 2 ensemble ID into one list 
Bio_gene = BasicInfo['Gene name'].to_list() #get ID from annotation col into one list
Ann_gene = BasicInfo['Gene_Name'].to_list() #get ID from biomart into one list 

#for loop to merge the two lsits 
Genes = []
for (a,b) in zip(Bio_gene, Ann_gene):
    if a != 0:
        ID.append(a)
    elif b != 0:
        ID.append(b)
    else:
        ID.append(0)

BasicInfo['Gene Name'] = Genes


##repeat for gene ID
BasicInfo['Gene ID'] = BasicInfo['Gene ID'].replace(np.nan, 0)
BasicInfo['GeneID'] = BasicInfo['GeneID'].replace(np.nan, 0)

#Merge 2 ensemble ID into one list 
Bio_ID = BasicInfo['Gene ID'].to_list() #get ID from annotation col into one list
Ann_ID = BasicInfo['GeneID'].to_list() #get ID from biomart into one list 

#for loop to merge the two lsits 
ID = []
for (a,b) in zip(Bio_ID, Ann_ID):
    if a != 0:  
        ID.append(a) 
    elif b != 0:
        ID.append(b)
    else:
        ID.append(0)

BasicInfo['Gene_ID '] = ID

#drop the old columns
BasicInfo = BasicInfo.drop(columns=["Gene_Name", 'GeneID', 'Gene ID', 'Gene name'])

#replace any 0 with nan
BasicInfo['Gene_ID '] = BasicInfo['Gene_ID '].replace(0, np.nan)
BasicInfo['Gene Name'] = BasicInfo['Gene Name'].replace(0, np.nan)

#rename columns 
BasicInfo = BasicInfo.rename(columns={'Gene_ID ':'GeneID', 'Gene Name':'Gene_Name'}, inplace=True)


'''Getting Alias'''


#get gene synonym dataframe 
alias = df_mart.dropna(subset=['Gene Synonym']) #remove psotions that do have a gene synonym
alias ['row_ID'] = alias.reset_index().index #add a row ID column 

#function to get positions and gene synonyms
pos_al = []
al = []
i=0
while i < len(alias):
    row = alias.loc[alias['row_ID']==i]
    w = list(range(int(row['Gene start (bp)']), (int(row['Gene end (bp)']))))
    pos_al.append(w)
    x = list(list(row['Gene Synonym']) * len(w))
    al.append(x)
    i=i+1

#flatten list of lists into one list 
position_alias= sum(pos_al, [])
Gene_Alias = sum(al, [])

#create dataframe with position and gene alias 
Alias = pd.DataFrame(position_alias, columns = ["position"])
Alias['Gene Alias'] = Gene_Alias

#get a gene info dataframe with positions from basic info table 
POS = BasicInfo["Position"].to_list()
boolean_series1 = Alias.position.isin(POS) #create boolean list of position in position list
Alias_df = Alias[boolean_series1] #filter list based on boolean series

#remove duplicates
Alias_df = Alias_df.drop_duplicates() 

#get the duplicated positions into one row and multiple alias into different columns
mask = Alias_df['position'].duplicated(keep=False) #mask all the duplicaed positions 
Alias_df.loc[mask, 'dup'] = Alias_df.groupby('position').cumcount().add(1).astype(str) #make new column with the duplicate number 
Alias_df1 = Alias_df.set_index(['position', 'dup'])['Gene Alias'].unstack() #make duplicated columns based on duplicaed columns
GeneAlias_df = Alias_df1.reset_index().rename_axis(None).rename_axis(None, axis=1) #flatten dataframe

#rename  and reoder columns
GeneAlias_df.columns = ['Position', 'NaN', 'Alias1', 'Alias10', 'Alias11', 'Alias12', 'Alias13', 'Alias14', 'Alias2', 'Alias3', 'Alias4', 'Alias5', 'Alias6', 'Alias7', 'Alias8', 'Alias9']
GeneAlias_df = GeneAlias_df[['Position', 'NaN', 'Alias1', 'Alias2', 'Alias3', 'Alias4', 'Alias5', 'Alias6', 'Alias7', 'Alias8', 'Alias9', 'Alias10', 'Alias11', 'Alias12', 'Alias13', 'Alias14']]

#drop empty column 
GeneAlias_df = GeneAlias_df.drop(columns=['NaN'])

#merge dataframe
BasicInfo = pd.merge(left=BasicInfo, right=GeneAlias_df, how='left', left_on='Position', right_on='Position')

#merge the alias columns into one seperated with -
cols = ['Gene_Name', 'Alias1', 'Alias2', 'Alias3', 'Alias4', 'Alias5', 'Alias6', 'Alias7', 'Alias8', 'Alias9', 'Alias10', 'Alias11', 'Alias12', 'Alias13', 'Alias14'] #list of columns to merge
BasicInfo["Alias"] = BasicInfo[cols].apply(lambda x: "-".join(x.dropna()), axis= 1) #merge columsn and drop empty values 

#drop columns that are not needed 
BasicInfo = BasicInfodrop(columns=['Alias1', 'Alias2', 'Alias3', 'Alias4', 'Alias5', 'Alias6', 'Alias7', 'Alias8', 'Alias9', 'Alias10', 'Alias11', 'Alias12', 'Alias13', 'Alias14', 'POS_y', 'position', 'row_ID'])

#drop duplicates 
BasicInfo = BasicInfo.drop_duplicates() 

#add row ID = will be primary key
BasicInfo = BasicInfo.reset_index().index

#Save it as zipped file
BasicInfo.to_csv("BasicInfo.csv.gz", 
           index=False, 
           compression="gzip")

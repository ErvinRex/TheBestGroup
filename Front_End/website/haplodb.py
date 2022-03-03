import pandas as pd 

class Haplo_df():
    Haplo_df_LWK_read = pd.read_csv('website/database/Haplotype_LWK.csv.gz')
    Haplo_df_LWK = Haplo_df_LWK_read
    Haplo_df_CLM_read = pd.read_csv('website/database/Haplotype_CLM.csv.gz')
    Haplo_df_CLM = Haplo_df_CLM_read
    Haplo_df_CHS_read = pd.read_csv('website/database/Haplotype_CHS.csv.gz')
    Haplo_df_CHS = Haplo_df_CHS_read
    Haplo_df_TSI_read = pd.read_csv('website/database/Haplotype_TSI.csv.gz')
    Haplo_df_TSI = Haplo_df_TSI_read
    Haplo_df_STU_read = pd.read_csv('website/database/Haplotype_STU.csv.gz')
    Haplo_df_STU = Haplo_df_STU_read

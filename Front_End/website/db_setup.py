import sqlite3
import pandas as pd

######################## ADD A if TABLE exists clause for when dockerising ? ############
def create_the_database():
    #connect to db
    con = sqlite3.connect('website/data.db')

    #cursor object
    cur = con.cursor()

    #create Allele Count table
    cur.execute('''CREATE TABLE AlleleCount(
                            POS, is_SNP, Ref_AC_LWK, Alt_AC_LWK, Ref_AC_CLM, Alt_AC_CLM,
                            Ref_AC_CHS, Alt_AC_CHS, Ref_AC_TSI, Alt_AC_TSI,
                            Ref_AC_STU, Alt_AC_STU, row_ID PRIMARY KEY)''')


    # load the data into a Pandas DataFrame -> AlleleCount
    AlleleCount = pd.read_csv('website/database/AlleleCount.csv.gz')

    # write the data to a sqlite table
    AlleleCount.to_sql('AlleleCount', con, if_exists='append', index = False)

    con.commit()

    #create basics info table
    cur.execute('''CREATE TABLE BasicInfo(
                            CHROM, Position, REF, ALT, is_snp, 
                            AF_LWK, AF_CLM, AF_CHS, AF_TSI, AF_STU, 
                            GTF_LWK_HomRef, GTF_LWK_HomAlt, GTF_LWK_Het,
                            GTF_CLM_HomRef, GTF_CLM_HomAlt, GTF_CLM_Het,
                            GTF_CHS_HomRef, GTF_CHS_HomAlt, GTF_CHS_Het,
                            GTF_TSI_HomRef, GTF_TSI_HomAlt, GTF_TSI_Het,
                            GTF_STU_HomRef, GTF_STU_HomAlt, GTF_STU_Het, 
                            rsID, AA, POS_y, Gene_Name, GeneID, Alias, DAF_LWK, DAF_CLM, DAF_CHS, DAF_TSI, DAF_STU, row_ID PRIMARY KEY)''')

    # load the data into a Pandas DataFrame -> basicsInfo
    BasicInfo = pd.read_csv('website/database/BasicInfo.csv.gz')

    # write the data to a sqlite table
    BasicInfo.to_sql('BasicInfo', con, if_exists='append', index = False)

    con.commit()


#Import, db=SQLalchemy
from . import db

#SQLAlchemy models for BasicInfo table
class BasicInfo(db.Model):
    __tablename__ = 'BasicInfo'

    CHROM = db.Column(db.String)
    Position = db.Column(db.Integer)
    REF = db.Column(db.String)
    ALT = db.Column(db.String)
    AF_LWK = db.Column(db.Float)
    AF_CLM = db.Column(db.Float)
    AF_CHS = db.Column(db.Float)
    AF_TSI = db.Column(db.Float)
    AF_STU = db.Column(db.Float)
    GTF_LWK_HomRef = db.Column(db.Float)
    GTF_LWK_HomAlt = db.Column(db.Float)
    GTF_LWK_Het = db.Column(db.Float)
    GTF_CLM_HomRef = db.Column(db.Float)
    GTF_CLM_HomAlt = db.Column(db.Float)
    GTF_CLM_Het = db.Column(db.Float)
    GTF_CHS_HomRef = db.Column(db.Float)
    GTF_CHS_HomAlt = db.Column(db.Float)
    GTF_CHS_Het = db.Column(db.Float)
    GTF_TSI_HomRef = db.Column(db.Float)
    GTF_TSI_HomAlt = db.Column(db.Float)
    GTF_TSI_Het = db.Column(db.Float)
    GTF_STU_HomRef = db.Column(db.Float)
    GTF_STU_HomAlt = db.Column(db.Float)
    GTF_STU_Het = db.Column(db.Float)
    rsID = db.Column(db.String)
    AA = db.Column(db.String)
    POS_y = db.Column(db.Integer)
    Gene_Name = db.Column(db.String)
    GeneID = db.Column(db.String)
    Alias = db.Column(db.String)
    DAF_LWK = db.Column(db.Float)
    DAF_CLM = db.Column(db.Float)
    DAF_CHS = db.Column(db.Float)
    DAF_TSI = db.Column(db.Float)
    DAF_STU = db.Column(db.Float)
    row_ID = db.Column(db.Integer, primary_key=True)

#SQLAlchemy models for Allelecount table
class AlleleCount(db.Model):
    __tablename__ = 'AlleleCount'

    POS = db.Column(db.Integer)
    is_SNP = db.Column(db.Boolean)
    Ref_AC_LWK = db.Column(db.Integer)
    Alt_AC_LWK = db.Column(db.Integer)
    Ref_AC_CLM = db.Column(db.Integer)
    Alt_AC_CLM = db.Column(db.Integer)
    Ref_AC_CHS = db.Column(db.Integer)
    Alt_AC_CHS = db.Column(db.Integer)
    Ref_AC_TSI = db.Column(db.Integer)
    Alt_AC_TSI = db.Column(db.Integer)
    Ref_AC_STU = db.Column(db.Integer)
    Alt_AC_STU = db.Column(db.Integer)
    row_ID = db.Column(db.Integer, primary_key=True)

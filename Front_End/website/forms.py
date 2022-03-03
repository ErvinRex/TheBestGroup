#Import statments
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, IntegerField, SelectField, SelectMultipleField
from wtforms.validators import DataRequired

#Forms for home page
class rsIDForm(FlaskForm):
    rsID = StringField('Enter rsID:', validators=[DataRequired()])
    submit = SubmitField('Submit')

class geneForm(FlaskForm):
    gene = StringField('Enter Gene:', validators=[DataRequired()])
    submit = SubmitField('Submit')

class geneidForm(FlaskForm):
    geneid = StringField('Enter GeneID:', validators=[DataRequired()])
    submit = SubmitField('Submit')

class posForm(FlaskForm):
    pos_start = StringField('Enter Position  Start:', validators=[DataRequired()])
    pos_end = StringField('End:', validators=[DataRequired()])
    submit = SubmitField('Submit')

#Forms for statistical tests and popultions selection page
class popsForm(FlaskForm):
    pops = SelectMultipleField('Select Population (use cmd for multiple)', choices=['LWK', 'CLM', 'CHS', 'TSI', 'STU'], validators=[DataRequired()])
    submit = SubmitField('Submit')

class statsForm(FlaskForm):
    stats = SelectMultipleField('Select one Statistic from each box', choices=["Tajima's D", "Wattersons estimator", "Haplotype diversity", "Fst"], validators=[DataRequired()])
    submit = SubmitField('Submit')


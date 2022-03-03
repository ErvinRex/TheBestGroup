#Import packages and modules
from flask import Blueprint, render_template, redirect, url_for, send_file
from .models import BasicInfo, AlleleCount
from .forms import rsIDForm, geneForm, posForm, geneidForm, popsForm, statsForm
import pandas as pd
import numpy as np
from .stats import Tajima_D, Watterson_Theta, fst_windowed, Haplotype_Diversity
from plotnine import ggplot, aes, geom_line, geom_bar
import itertools
import sqlite3
import plotly
import plotly.express as px
from .haplodb import Haplo_df
from plotly.subplots import make_subplots
import plotly.graph_objects as go 

#Blueprint to be registered with the application factory 
bp = Blueprint('home', __name__)

#Global variables to pass data objects between routes
class Data():
    rsID = None
    Gene = None
    GeneID = None
    reflist_all_pops= None
    altlist_all_pops= None
    stats = None
    positions_stats = None
    pos_start = None
    pos_end = None
    search_data = None
    pops = None
    Fst_plot = None
    TjD_plot = None
    WeT_plot = None
    HaD_plot = None
    TjD_dwn = None
    WeT_dwn  = None
    HaD_dwn = None
    Fst_dwn = None

data = Data()

#Home route - contains search forms for rsID, Gene, GeneID, Chromosome Position
@bp.route('/', methods=['GET', 'POST'])
def home():

    #Reset global variables when arriving at/returning to home page
    data.rsID = None
    data.Gene = None
    data.GeneID = None
    data.reflist_all_pops= None
    data.altlist_all_pops= None
    data.stats = None
    data.positions_stats = None
    data.pos_start = None
    data.pos_end = None
    data.search_data = None
    data.pops = None
    #data.error = None
    data.WeT_plot = None
    data.HaD_plot = None
    data.Fst_plot = None
    data.TjD_plot = None
    data.TjD_dwn = None
    data.WeT_dwn = None
    data.HaD_dwn = None
    data.Fst_dwn = None


    #SEARCH RSID
    form_rsID = rsIDForm()
    rsID = None
    
    #Check for form submission 
    if form_rsID.validate_on_submit():
        #Collect rsID from form and store in global variable data.rsID
        rsID = form_rsID.rsID.data
        data.rsID = rsID

        return redirect(url_for('home.thread', page_num=1)) #Return to results page 1 
    

    #SEARCH GENE
    form_gene = geneForm()
    Gene = None

    #Check for form submission 
    if form_gene.validate_on_submit():
        #Collect gene name from form and store in global variable data.Gene
        Gene = form_gene.gene.data
        data.Gene = Gene
        loader = 'loading'
        return redirect(url_for('home.thread', page_num=1)) #Return to results page 1
   

   #SEARCH GENEID
    form_geneid = geneidForm()
    geneID = None
    
    #Check for form submission 
    if form_geneid.validate_on_submit():
        #Collect gene ID from form and store in global variable data.GeneID
        geneID = form_geneid.geneid.data
        data.GeneID = geneID

        return redirect(url_for('home.thread', page_num=1)) #Return to results page 1


    #SEARCH POSITION
    form_pos = posForm()
    pos_start = None
    pos_end = None

    #Check for form submission 
    if form_pos.validate_on_submit():

        #Collect start position from form
        pos_start = form_pos.pos_start.data
        pos_start = int(pos_start)
        data.pos_start = pos_start

        #Collect end position from form
        pos_end = form_pos.pos_end.data
        pos_end = int(pos_end)
        data.pos_end = pos_end

        return redirect(url_for('home.thread', page_num=1)) #Return to results page 1


    #Return home template to display the user search options
    return render_template('home.html', rsID = rsID, form_rsID = form_rsID,
                                        form_gene = form_gene, Gene = Gene,
                                        form_geneid = form_geneid, geneID = geneID,
                                        form_pos = form_pos, pos_start = pos_start, pos_end = pos_end
                                        )



####===--> ABOUT ROUTE <--===####
#Route for the application about page ***
@bp.route('/about')
def about():
    return render_template('about.html')

####===--> DOCUMENTATION ROUTE <--===####
#Route for the application about page ***
@bp.route('/how_to')
def how_to():
    return render_template('how_to.html')

####===--> THREAD ROUTE <--===####
#Route for basic information search results
@bp.route('/<int:page_num>', methods=['GET', 'POST'])
def thread(page_num):
    error = None
    #Return rsID basic information if rsID queried 
    if data.rsID:
        threads = BasicInfo.query.filter(BasicInfo.rsID == data.rsID).paginate(per_page=50, page=page_num) #Query BasicInfo table & Return first 50 results on page 1
        thready = BasicInfo.query.filter(BasicInfo.rsID == data.rsID).first()
        
        #Return message if no matches to BasicInfo table
        if thready == None:
            error = ('Sorry, no results for' + ' ' + str(data.rsID))
        data.search_data = data.rsID

    #Return Gene basic information if Gene queried 
    if data.Gene:
        #threads = BasicInfo.query.filter(BasicInfo.Gene_Name == data.Gene).paginate(per_page=50, page=page_num) #Query BasicInfo table & Return first 50 results on page 1
        data.search_data = data.Gene
        thready = BasicInfo.query.filter(BasicInfo.Alias.like('%{}%'.format(data.Gene))).first()

        #threads = BasicInfo.query.filter(or_(BasicInfo.Gene_Name == data.Gene, BasicInfo.GeneID == data.Gene)).paginate(per_page=50, page=page_num)
        threads = BasicInfo.query.filter(BasicInfo.Alias.like('%{}%'.format(data.Gene))).paginate(per_page=50, page=page_num)
        #Return message if no matches to BasicInfo table
        if thready == None:
            error = ('Sorry, no results for' + ' ' + str(data.Gene))
    
    #Return Gene ID basic information if Gene ID queried 
    if data.GeneID:
        threads = BasicInfo.query.filter(BasicInfo.GeneID == data.GeneID).paginate(per_page=50, page=page_num) #Query BasicInfo table & Return first 50 results on page 1
        data.search_data = data.GeneID
        thready = BasicInfo.query.filter(BasicInfo.GeneID == data.GeneID).first()

        #Return message if no matches to BasicInfo table
        if thready == None:
            error = ('Sorry, no results for' + ' ' + str(data.GeneID))
                                            
    #Return basic information for position range query 
    if data.pos_start:
        threads = BasicInfo.query.filter(BasicInfo.Position >= data.pos_start, BasicInfo.Position <= data.pos_end).paginate(per_page=50, page=page_num) #Return first 50 results on page 1
        #if pos not in db, modify pos till matches with db
        data.search_data = 'positions' + ' ' + str(data.pos_start) + ' ' + '-' + ' ' + str(data.pos_end)
        thready = BasicInfo.query.filter(BasicInfo.Position >= data.pos_start, BasicInfo.Position <= data.pos_end).first()

        #Return message if no matches to BasicInfo table
        if thready == None:
            error = 'Sorry, no results for those positions'
        
    #Return the basic search template with the variables for the basic information table
    return render_template('basic_search.html', threads=threads, search_data = data.search_data, error=error, rsID = data.rsID)



####===--> STATS + POPS SELECTION ROUTE <--===####
@bp.route('/stats', methods=['GET', 'POST'])
def stats():

    positions = None

    #If gene queried, get associated postions list
    if data.Gene:
        #Get positions for queried Gene
        pos_for_stats = BasicInfo.query.filter(BasicInfo.Alias.like('%{}%'.format(data.Gene))).all() #Query BasicInfo table
        positions = []
        #Place postitions into list for slicing 
        for x in pos_for_stats:
            positions.append(x.Position)
    
    #If gene ID queried, get associated postions list
    if data.GeneID:
        #Get positions for queried Gene ID
        pos_for_stats = BasicInfo.query.filter(BasicInfo.GeneID == data.GeneID).all() #Query BasicInfo table
        positions = []
        print(len(positions))
        #Place postitions into list for slicing 
        for x in pos_for_stats:
            positions.append(x.Position)
    
    #If chromosome positions queried, get associated postions list
    if data.pos_start:
        #Get positions for queried positions
        positions = []
        #Place postitions into list for slicing
        positions.append(data.pos_start)
        positions.append(data.pos_end)
    
    #Get allele count data for user query (if not rsID/returns multiple SNPs)
    if data.rsID == None:
        
        #Slice poitions list to query first and last postitions in AlleleCount table
        s=positions[0]
        e=positions[-1]
        allele_query = AlleleCount.query.filter(AlleleCount.POS >= s, AlleleCount.POS <= e).all() #Query AlleleCount table using postions associated with user query e.g Gene

        #Get positions list from AlleleCount table to ensure data list lengths are equal for Allelecount and Haplotype diversity stats
        poslist=[]
        for x in allele_query:
                poslist.append(x.POS)
        
        #Store allele count positions data in global variable data.positions_stats
        data.positions_stats = poslist
    
    #Population selection form
    form_pops = popsForm()

    #Check for user submission and get data from Population selection form, store in pops variable
    if form_pops.validate_on_submit():
        
        pops = form_pops.pops.data
        data.pops = pops
        form_pops.pops.data = None

        #Lists of Allelecount data (both reference allele counts and alternate allele counts) 
        #For every population the user has selected for
        reflist_all_pops=[]
        altlist_all_pops=[]
        for pop in pops:

            reflist=[]
            altlist=[]

            #Get allele count data for population in list, using positions from allele_query
            for x in allele_query:
                if pop == 'LWK':
                    reflist.append(x.Ref_AC_LWK)
                    altlist.append(x.Alt_AC_LWK)

                if pop == 'CLM':
                    reflist.append(x.Ref_AC_CLM)
                    altlist.append(x.Alt_AC_CLM)
                
                if pop == 'CHS':
                    reflist.append(x.Ref_AC_CHS)
                    altlist.append(x.Alt_AC_CHS)
                
                if pop == 'TSI':
                    reflist.append(x.Ref_AC_TSI)
                    altlist.append(x.Alt_AC_TSI)
                
                if pop == 'STU':
                    reflist.append(x.Ref_AC_STU)
                    altlist.append(x.Alt_AC_STU)

            #Create list containing allele count data for all populations selected
            reflist_all_pops.append(reflist)
            altlist_all_pops.append(altlist)

            #Store allele count data for all populations in global variables reflist_all_pops & altlist_all_pops
            data.reflist_all_pops = reflist_all_pops
            data.altlist_all_pops = altlist_all_pops

    #Statistical tests selection form
    form_stats = statsForm()

    #Check for user submission of stats form and store in stats variable 
    if form_stats.validate_on_submit():
        stats = form_stats.stats.data
        form_stats.stats.data = None
        data.stats = stats
    
    #Return template that displays the options to select the users choice of populations and statistical tests
    return render_template('pops_stats_selection.html', form_pops=form_pops, form_stats=form_stats)


####===--> ROUTE FOR STATISTICAL TESTS LOGIC <--===####
@bp.route('/stats_2', methods=['GET', 'POST'])
def stats_2():
    
    overall_HaD = None
    overall_WeT = None
    overall_TjD = None
    overall_Fst = None

    #If only one population is selected, display message for user
    #Message tells user that Fst requires two populations selected 
    Jinja_path_block = None
    Fst_path = 0
    if len(data.pops) == 1:
        Fst_path += 1
    else:
        pass

    if 'Fst' in data.stats:
        Fst_path += 1
    else:
        pass
    
    if Fst_path == 2:
        Jinja_path_block = 'block'
    else:
        pass
    
    #### Stats logic ####
    #If Tajima's D selected, run the corresponding logic 
    if "Tajima's D" in data.stats:
        
        TjD_val = [] #Variable for displaying the overall statistical value for the queried positions
        df_TjD_dwn = pd.DataFrame(data.positions_stats, columns = ['positions']) #Dataframe for plots data download
        fig_TjD = make_subplots(rows=1, cols=1) #Empty plot for user selected population data to be added
        
        #If multiple population are selected
        if len(data.reflist_all_pops) > 1:
            
            #Loop through populations and the corresponding AlleleCount reference and alternate data lists
            for i in range(0,len(data.reflist_all_pops)):
                df_TjD, TjD = Tajima_D(data.positions_stats, data.reflist_all_pops[i], data.altlist_all_pops[i]) #Run stats function, returns dataframe and overall stat value
                df_TjD_dwn['TajimaD' + str(data.pops[i])] = df_TjD['TajimaD'] #Add to dataframe for download 
                fig_TjD.add_trace(go.Scatter(x=df_TjD["Position"], y=df_TjD["TajimaD"], name = str(data.pops[i]))) #Add plot data to graph for population in loop
                TjD_val.append(TjD) #Append overall stat value of search query to be displayed on stats output page
        
        #If a single population is selected
        else:

            #Using AlleleCount reference and alternate data list single item
            df_TjD, TjD = Tajima_D(data.positions_stats, data.reflist_all_pops[0], data.altlist_all_pops[0]) #Run stats function, returns dataframe and overall stat value
            df_TjD_dwn['TajimaD' + str(data.pops[0])] = df_TjD['TajimaD'] #Add to dataframe for download 
            TjD_val.append(TjD) #Append overall stat value of search query to be displayed on stats output page
            fig_TjD.add_trace(go.Scatter(x=df_TjD["Position"], y=df_TjD["TajimaD"], name = str(data.pops[0]), showlegend=True)) #Add plot data to graph for population 
        
        fig_TjD.update_layout(xaxis_title='Chromosome position', yaxis_title="Tajima's D")
        data.TjD_plot = fig_TjD.to_html(full_html=False) #Global variable used for displaying plot
        data.TjD_dwn = df_TjD_dwn #Global variable for downloading plots data
        overall_TjD = list(zip(TjD_val, data.pops)) #Zip overall stats value with population codes
        
    #If Wattersons estimator selected, run the corresponding logic 
    if "Wattersons estimator" in data.stats:

        df_WeT_dwn = pd.DataFrame(data.positions_stats, columns = ['positions'])  #Dataframe for plots data download
        WeT_val = [] #Variable for displaying the overall statistical value for the queried positions
        fig_WeT = make_subplots(rows=1, cols=1) #Empty plot for user selected population data to be added

        #If multiple population are selected
        if len(data.reflist_all_pops) > 1:
            
            #Loop through populations and the corresponding AlleleCount reference and alternate data lists
            for i in range(0,len(data.reflist_all_pops)):
                df_WeT, WeT = Watterson_Theta(data.positions_stats, data.reflist_all_pops[i], data.altlist_all_pops[i]) #Run stats function, returns dataframe and overall stat value
                df_WeT_dwn['Watterson_Theta' + str(data.pops[i])] = df_WeT['Watterson_Theta'] #Add to dataframe for download
                fig_WeT.add_trace(go.Scatter(x=df_WeT["Position"], y=df_WeT["Watterson_Theta"], name = str(data.pops[i]))) #Add plot data to graph for population in loop
                WeT_val.append(WeT) #Append overall stat value of search query to be displayed on stats output page
            
        #If a single population is selected        
        else:

            #Using AlleleCount reference and alternate data list single item
            df_WeT, WeT = Watterson_Theta(data.positions_stats, data.reflist_all_pops[0], data.altlist_all_pops[0]) #Run stats function, returns dataframe and overall stat value
            df_WeT_dwn['Watterson_Theta' + str(data.pops[0])] = df_WeT['Watterson_Theta'] #Add to dataframe for download
            WeT_val.append(WeT) #Append overall stat value of search query to be displayed on stats output page
            fig_WeT.add_trace(go.Scatter(x=df_WeT["Position"], y=df_WeT["Watterson_Theta"], name = str(data.pops[0]), showlegend=True)) #Add plot data to graph for population 
        
        fig_WeT.update_layout(xaxis_title='Chromosome position', yaxis_title="Wattersons estimator Theta")
        data.WeT_plot = fig_WeT.to_html(full_html=False) #Global variable used for displaying plot
        data.WeT_dwn = df_WeT_dwn #Global variable for downloading plots data
        overall_WeT = list(zip(WeT_val, data.pops)) #Zip overall stats value with population codes
    
    #If Haplotype diversity selected, run the corresponding logic 
    if "Haplotype diversity" in data.stats:
        
        #Get Haplotype dataframes from Haplo_df module
        hapdf_LWK = Haplo_df.Haplo_df_LWK
        hapdf_CLM = Haplo_df.Haplo_df_CLM
        hapdf_CHS = Haplo_df.Haplo_df_CHS
        hapdf_TSI = Haplo_df.Haplo_df_TSI
        hapdf_STU = Haplo_df.Haplo_df_STU

        df_HaD_dwn = pd.DataFrame(data.positions_stats, columns = ['positions']) #Dataframe for plots data download
        fig_HaD = make_subplots(rows=1, cols=1) #Empty plot for user selected population data to be added
        HaD_val = [] #Variable for displaying the overall statistical value for the queried positions

        #For each conditional statement:
        #Run stats function, returns dataframe and overall stat value
        #Add plot data to graph for population
        #Add to dataframe for download
        #Append overall stat value of search query to be displayed on stats output page
        for pop in data.pops:
            if pop == 'LWK':                 
                df_LWK, HaD = Haplotype_Diversity(hapdf_LWK, data.positions_stats)
                fig_HaD.add_trace(go.Scatter(x=df_LWK["Position"], y=df_LWK["HapDiv"], name = 'LWK'))
                df_HaD_dwn['Haplotype_diversity_LWK'] = df_LWK['HapDiv']
                HaD_val.append(HaD)

            elif pop == 'CHS':                 
                df_CHS, HaD = Haplotype_Diversity(hapdf_CHS, data.positions_stats)
                fig_HaD.add_trace(go.Scatter(x=df_CHS["Position"], y=df_CHS["HapDiv"], name = 'CHS'))
                df_HaD_dwn['Haplotype_diversity_CHS'] = df_CHS['HapDiv']
                HaD_val.append(HaD)

            elif pop == 'CLM':                 
                df_CLM, HaD = Haplotype_Diversity(hapdf_CLM, data.positions_stats)
                fig_HaD.add_trace(go.Scatter(x=df_CLM["Position"], y=df_CLM["HapDiv"], name = 'CLM'))
                df_HaD_dwn['Haplotype_diversity_CLM'] = df_CLM['HapDiv']
                HaD_val.append(HaD)

            elif pop == 'TSI':                 
                df_TSI, HaD = Haplotype_Diversity(hapdf_TSI, data.positions_stats)
                fig_HaD.add_trace(go.Scatter(x=df_TSI["Position"], y=df_TSI["HapDiv"], name = 'TSI'))
                df_HaD_dwn['Haplotype_diversity_TSI'] = df_TSI['HapDiv']
                HaD_val.append(HaD)

            else:                 
                df_STU, HaD = Haplotype_Diversity(hapdf_STU, data.positions_stats)
                fig_HaD.add_trace(go.Scatter(x=df_STU["Position"], y=df_STU["HapDiv"], name = 'STU'))
                df_HaD_dwn['Haplotype_diversity_STU'] = df_STU['HapDiv']
                HaD_val.append(HaD)
        
        fig_HaD.update_layout(xaxis_title='Chromosome position', yaxis_title="Haplotype diversity")
        fig_HaD.update_traces(showlegend=True) #Display legend
        data.HaD_plot = fig_HaD.to_html(full_html=False) #Global variable used for displaying plot
        data.HaD_dwn = df_HaD_dwn #Global variable for downloading plots data
        overall_HaD = list(zip(HaD_val, data.pops)) #Zip overall stats value with population codes
    
    #If Fst selected, run the corresponding logic 
    message_for_user = None
    if "Fst" in data.stats:
        #If only one population selected, display message to user
        if Jinja_path_block:
            message_for_user = 'Please select at least 2 popuations for Fst'
        
        #If multiple populations selected
        else:
        
            fst_pairs_ref = list(itertools.combinations(data.reflist_all_pops, 2)) #get reference list combos for fst dataframes
            fst_pairs_alt = list(itertools.combinations(data.altlist_all_pops, 2)) #get alternate list combos for fst dataframes
            key_name_combo = list(itertools.combinations(data.pops, 2)) #get list of population code combinations 

            zipped_fst = list(zip(fst_pairs_ref, fst_pairs_alt)) #zip reference and alternate list combos

            #get dataframes for each population combination
            fig_Fst = make_subplots(rows=1, cols=1) #empty plot
            Fst_val = []
            key_index = 0
            names_overall = []
            df_Fst_dwn = pd.DataFrame(data.positions_stats, columns = ['positions']) #empty dataframe for plot data download
            for i in zipped_fst:
                
                df_ref = pd.DataFrame(i[0])
                xref = df_ref.transpose()
                df_alt = pd.DataFrame(i[1])
                xalt = df_alt.transpose()
                
                #Generate legend for plot
                key_name = key_name_combo[key_index]
                name = str(key_name[0] + ' ' + '-' + ' ' + key_name[1])
                dwn_col_name = str('Fst-' + name)

                #Store legend names for overall values list
                names_overall.append(name)
            
                HFs, AHF = fst_windowed(data.positions_stats, xref, xalt) #Run stats function, returns dataframe
                Fst_val.append(AHF)
                fig_Fst.add_trace(go.Scatter(x=HFs["Position"], y=HFs["FST"], name = name, showlegend=True)) #Add plot data for each population combination 
                df_Fst_dwn[dwn_col_name] = HFs['FST'] #Add data to dataframe for download 
                key_index += 1

            fig_Fst.update_layout(xaxis_title='Chromosome position', yaxis_title="Fst")
            data.Fst_plot = fig_Fst.to_html(full_html=False) #Global variable used for displaying plot
            data.Fst_dwn = df_Fst_dwn #Global variable for downloading plots data
            overall_Fst = list(zip(Fst_val, names_overall)) #Zip overall stats value with population codes

    return render_template('stats_output.html', Fst_plot=data.Fst_plot, TjD_plot=data.TjD_plot, WeT_plot=data.WeT_plot, HaD_plot=data.HaD_plot,
                                            overall_TjD = overall_TjD, overall_WeT = overall_WeT, overall_HaD = overall_HaD, message_for_user = message_for_user, overall_Fst = overall_Fst)
    

####===--> ROUTE FOR STATISTICAL OUTPUTS <--===####
@bp.route('/stats_3', methods=['GET', 'POST'])
def stats_3():
    return render_template('stats_output.html')

####===--> ROUTES FOR STATISTICAL DOWNLOADS <--===####
#Route for Tajima's D plot data download
@bp.route('/getCSV_TjD', methods=['GET', 'POST'])
def getCSV_TjD():
    dfTjD = data.TjD_dwn
    dfTjD.to_csv('website/TjD-plot-data.csv', index=False)
    
    return send_file('TjD-plot-data.csv')

#Route for Waterson's estimator theta plot data download
@bp.route('/getCSV_WeT', methods=['GET', 'POST'])
def getCSV_WeT():
    dfWeT = data.WeT_dwn
    dfWeT.to_csv('website/WeT-plot-data.csv', index=False)
   
    return send_file('WeT-plot-data.csv')

#Route for Haplotype Diversity plot data download
@bp.route('/getCSV_HaD', methods=['GET', 'POST'])
def getCSV_HaD():
    dfHaD = data.HaD_dwn
    dfHaD.to_csv('website/HaD-plot-data.csv', index=False)
    
    return send_file('HaD-plot-data.csv')

#Route for Fst plot data download
@bp.route('/getCSV_Fst', methods=['GET', 'POST'])
def getCSV_Fst():
    dfFst = data.Fst_dwn
    dfFst.to_csv('website/Fst-plot-data.csv', index=False)
   
    return send_file('Fst-plot-data.csv')
#import statments 
import allel
import pandas as pd
import numpy as np
from plotnine import ggplot, aes, geom_line
import plotly.express as px

#function to run tajima D and get a sliding window as a dataframe and overall value

def Tajima_D(pos, ref, alt):
    ac = np.array(list(zip(ref, alt))) #convert 2 lists into aray
    if len(pos) < 20: #get the window size 
        size = 2
    else:
        size = len(pos)/10
    size = int(size)
    window_ac = [] 
    for x in range(0, (len(pos)-size)): #loop to get start and end positoons for windows 
        a = (pos[x], pos[x+size])
        window_ac.append(a)
    D, windows, counts = allel.windowed_tajima_d(pos = pos, ac = ac, windows = window_ac) #tajima D test with sliding windows 
    TjD = allel.tajima_d(ac=ac, pos=pos) #tajima D value for whole region
    x = windows.mean(axis=1) #get the mean for each window
    df = pd.DataFrame(D) #make a dataframe from tajima D resuls 
    df["position"] = x #make the positon to x xis 
    df.columns = ["TajimaD", "Position"] #change the name of the positions 
    df['TajimaD'] = df['TajimaD'].replace(np.nan, 0) #replace NAN to 0
    
    return (df, TjD)

#function to get the watterson theta values 
def Watterson_Theta(pos, ref, alt):
    ac = np.array(list(zip(ref, alt))) #convert 2 lists into aray
    if len(pos) < 20: #get the window size 
        size = 2
    elif len(pos) <100000:
        size = len(pos)/10
    else:
        size = len(pos)/100
    size = int(size)
    window_ac = [] 
    for x in range(0, (len(pos)-size)): #loop to get start and end positoons for windows 
        a = (pos[x], pos[x+size])
        window_ac.append(a)
    theta_hat_w, windows, n_bases, counts = allel.windowed_watterson_theta(pos, ac, windows = window_ac) #watterson theta value with sliding windows
    wt = allel.watterson_theta(ac=ac, pos=pos) #watteron theta value without sliding widow
    x = windows.mean(axis=1) #get the mean for each window
    df= pd.DataFrame(theta_hat_w) #make a dataframe from tajima D resuls 
    df["position"] = x #make the positon to x xis 
    df.columns = ["Watterson_Theta", "Position"] #change the name of the positions 
    
    return (df, wt)

    
#get the Haplotype diversity based on the haplotype array 
def Haplotype_Diversity(df, pos):
    boolean_series = df.POS.isin(pos) #create boolean list of position in position list
    filter_df = df[boolean_series] #filter dataframe to only have positions of interest
    poslist = filter_df["POS"].to_list()
    filter_df = filter_df.drop(columns=["POS", "row_ID"]) #removes columns that are not needed in array
    ht = filter_df.to_numpy() #convert dataframe to array
    if len(poslist) < 20: #get the window size 
        size = 2
    elif len(poslist) <100000:
        size = len(poslist)/10
    else:
        size = len(poslist)/100
    size = int(size) #get size of windows (1/100 of array length) as integer
    hdiv = allel.moving_haplotype_diversity(ht, size = size, step = 1, stop = (len(ht)-1)) #get the haplotype diversity 
    h_div = allel.haplotype_diversity(ht) #haplotype diversity without sliding window
    window = [] #empty window as list 
    i = 0  #start with the first position 
    while i < (len(poslist)-(size)): #while loop to get the start and end position for each window in function
        a = (poslist[i], poslist[i+size])
        window.append(a)
        i = i + 1 #next position in the step value
    window = np.array(window) #convert window list to array
    x = window.mean(axis=1) #get the mean value for each position 
    df_hdiv = pd.DataFrame(hdiv) #convert haplotype diversity value to dataframe 
    df_hdiv["position"] = x #add position column to dataframe
    df_hdiv.columns = ["HapDiv", "Position"] #rename columns 
    return (df_hdiv, h_div)

def fst_windowed(pos, ref, alt):
    '''
    Finds, using the hudson method, a windowed FST value for all SNPs in respective populations and plots this onto an interactive graph
    Parameters:
    pos (pos.array): Positional array of SNPs
    ref (ref.allel_count): Reference dataframe allele counts of both populations
    alt (alt.allel_count): Alternate dataframe allele counts of both populations
    Returns:
    fig.show(): Interactive graph for FST between both populations
    '''

    #Convert allel counts into arrays (here I have to split the two lists)
    #names = ref.columns.str.split('_').str[2]
    ref1 = ref.iloc[:,[0]].to_numpy()
    ref2 = ref.iloc[:,[1]].to_numpy()
    alt1 = alt.iloc[:,[0]].to_numpy()
    alt2 = alt.iloc[:,[1]].to_numpy()
    acs = np.column_stack((ref1 , alt1))
    acs2 = np.column_stack((ref2 , alt2))
    #Convert position into positional array
    pos=np.array(pos)
    #Create windows for plotting
    if len(pos) < 20: #get the window size 
        size = 2
    elif len(pos) <100:
        size = 5
    else:
        size = len(pos)/100
    size = int(size)
    windows1 = [] 
    for x in range(0, (len(pos)-size)): #loop to get start and end positons for windows 
        a = (pos[x], pos[x+size])
        windows1.append(a)
    
    #Average hudson fst
    AHF, se, vb, vj = allel.average_hudson_fst(ac1= acs, ac2 = acs2, blen = len(pos))

    #Hudson Fst with sliding windows specified 
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
    return HFs, AHF
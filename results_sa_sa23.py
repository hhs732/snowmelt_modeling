###       /bin/bash runTestCases_docker.sh
#nSnow - time-varying number of snow layers
#nSoil - time-varying number of soil layers
#nLayers - time-varying number of total layers
#We also we need to know where along the [mid|ifc][Snow|Soil|Toto]andTime dimension the information for a specific time step is stored. SUMMA stores this information in the variables
#
#ifcSnowStartIndex (may not be present if there is no snow)
#ifcSoilStartIndex
#ifcTotoStartIndex
#midSnowStartIndex (may not be present if there is no snow)
#midSoilStartIndex
#midTotoStartIndex
#Note that these indices are 1-based, so if you are using python, C, or any other 0-based language, you'll need to substract 1 when you use these indices to index a variable arranged along a [mid|ifc][Snow|Soil|Toto]andTime dimension.
#
#This information can now be used to extract the information you need. This may best be demonstrated with an example. For example, if you save the variable mLayerVolFracWat in the history file, then you'll find that this variable is arranged along the midTotoAndTime and hru dimensions in the NetCDF file
#
#double mLayerVolFracWat(midTotoAndTime, hru)
#If you want to extract the information for the first timestep then in pseudo-code you want to do something like this (assuming a 0-based language like python)
#
## Note that the following is not real code. Typing this into python will not work as-is
#hru = 0 # python indices are 0-based, so this is the first hru
#timestep = 0 # python indices are 0-based, so this is the first time step
#layers = nLayers[timestep, hru] # extract the number of layers associated with the first timestep
#startIndex = midTotoStartIndex[timestep, hru] - 1 # - 1 since the SUMMA indices are 1-based and python indices are 0-based
#endIndex = startIndex + layers
#moistureProfile = mLayerVolFracWat[startIndex:endIndex, hru] # moisture profile at the first timestep
#To get the actual vertical locations (rather than just the vertical indices), you should also include mLayerHeight in the history file. Since this is arranged along the same dimension, you can extract the actual vertical locations (in m) as mLayerHeight[startIndex:endIndex, hru]
#%%
import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats
from sklearn.metrics import mean_squared_error
import itertools
import csv
#%% SWE observation data 
date_swe = ['2006-11-01 11:10','2006-11-30 12:30','2007-01-01 11:10','2007-01-30 10:35','2007-03-05 14:30','2007-03-12 14:00', 
            '2007-03-19 12:30','2007-03-26 12:30','2007-04-02 12:30','2007-04-18 08:35','2007-04-23 10:30','2007-05-02 08:40', 
            '2007-05-09 08:50','2007-05-16 09:00','2007-05-23 08:30','2007-05-30 09:00','2007-06-06 08:15', 
            
            '2007-12-03 10:45','2008-01-01 11:30','2008-01-31 12:00','2008-03-03 14:30','2008-03-24 09:10','2008-04-01 09:55', 
            '2008-04-14 14:45','2008-04-22 12:30','2008-04-28 12:30','2008-05-06 09:15','2008-05-12 12:45','2008-05-19 10:40',
            '2008-05-26 08:45','2008-06-02 12:45','2008-06-08 08:45'] 
            
swe_mm = [58,  169, 267, 315, 499, 523, 503, 549, 611, 678, 654, 660, 711, 550, 443, 309, 84, 
          141, 300, 501, 737, 781, 837, 977, 950, 873, 894, 872, 851, 739, 538, 381]  

#obs_swe_date = pd.DataFrame (np.column_stack([date_swe,swe_mm]), columns=['date_swe','swe_mm'])
obs_swe = pd.DataFrame (swe_mm, columns=['swe_mm'])
obs_swe.set_index(pd.DatetimeIndex(date_swe),inplace=True)

max_swe_obs = max(obs_swe['swe_mm'])
max_swe_date_obs = obs_swe[obs_swe ['swe_mm']== max_swe_obs].index.tolist()    
#%% Snow depth observation data
with open("snowDepth_2006_2008.csv") as safd1:
    reader1 = csv.reader(safd1)
    raw_snowdepth1 = [r for r in reader1]
sa_snowdepth_column1 = []
for csv_counter1 in range (len (raw_snowdepth1)):
    for csv_counter2 in range (2):
        sa_snowdepth_column1.append(raw_snowdepth1[csv_counter1][csv_counter2])
sa_snowdepth=np.reshape(sa_snowdepth_column1,(len (raw_snowdepth1),2))
sa_sd_obs=[np.array(val) for val in sa_snowdepth[1:len(raw_snowdepth1)-1,1:]]
#sa_sd_obs = [float(value) for value in sa_snowdepth1]
sa_sd_obs_date = pd.DatetimeIndex(sa_snowdepth[1:len(raw_snowdepth1)-1,0])

snowdepth_obs_df = pd.DataFrame(sa_sd_obs, columns = ['observed_snowdepth']) 
snowdepth_obs_df.set_index(sa_sd_obs_date,inplace=True)

#%%
p1 = [0.1] #LAIMIN
p2 = [1] #LAIMAX
p3 = [0.1] #winterSAI
p4 = [0.9] #summerLAI
p5 = [0.5] #rootingDepth
p6 = [0.5] #heightCanopyTop
p7 = [0.01] #heightCanopyBottom
p8 = [0.89] #throughfallScaleSnow
p9 = [55] #newSnowDenMin 
p10 = [500000, 1000000, 1300000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p11 = [0.8, 0.9, 0.94] #albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p12 = [0.6, 0.68, 0.74] #albedoMinVisible |       0.7500 |       0.5000 |       0.7500
p13 = [0.55, 0.65, 0.7] #albedoMaxNearIR |       0.6500 |       0.5000 |       0.7500
p14 = [0.2, 0.3, 0.4] #albedoMinNearIR  |       0.3000 |       0.1500 |       0.4500
p15 = [0.002] #[0.001, 0.002] #z0Snow
p16 = [6]# 1, 3, 6] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
p17 = [4] #2, 3, 4] #mw_exp exponent for meltwater flow
p18 = [0.6] #0.2, 0.4 , 0.6] #fixedThermalCond_snow

def hru_ix_ID(p10, p11, p12, p13, p14):
    ix10 = np.arange(1,len(p10)+1)
    ix11 = np.arange(1,len(p11)+1)
    ix12 = np.arange(1,len(p12)+1)
    ix13 = np.arange(1,len(p13)+1)
    ix14 = np.arange(1,len(p14)+1)

    c = list(itertools.product(ix10,ix11,ix12,ix13,ix14))
    ix_numlist=[]
    for tup in c:
        ix_numlist.append(''.join(map(str, tup)))
    new_list = [float(i) for i in ix_numlist]

    return(new_list)  

hruidx = hru_ix_ID(p10, p11, p12, p13, p14)

hruidxID = []
for index in hruidx:
    hruidxID.append(int(index))
    

hru_num = np.size(hruidxID)
years = ['2006','2007']
out_names = ['lj1110',
#              'lj1120','lj1130','lj1210','lj1220','lj1230','lj1310','lj1320','lj1330','lj2110','lj2120','lj2130','lj2210','lj2220','lj2230','lj2310','lj2320','lj2330',
#              'mj1110','mj1120','mj1130','mj1210','mj1220','mj1230','mj1310','mj1320','mj1330','mj2110','mj2120','mj2130','mj2210','mj2220','mj2230','mj2310','mj2320','mj2330',
#              'sj1110','sj1120','sj1130','sj1210','sj1220','sj1230','sj1310','sj1320','sj1330',
#              #'sj2110',
#              'sj2120','sj2130','sj2210','sj2220','sj2230','sj2310','sj2320','sj2330',
#              'lc1111','lc1112','lc1113','lc1121','lc1122','lc1123','lc1131','lc1132','lc1133','lc1211','lc1212','lc1213','lc1221','lc1222','lc1223','lc1231','lc1232','lc1233','lc1311','lc1312','lc1313','lc1321','lc1322','lc1323','lc1331','lc1332','lc1333','lc2111','lc2112','lc2113','lc2121','lc2122','lc2123','lc2131','lc2132','lc2133','lc2211',
#              'sc1111','sc1112','sc1113','sc1121','sc1122','sc1123','sc1131','sc1132','sc1133','sc1211','sc1212','sc1213','sc1221','sc1222','sc1223','sc1231','sc1232',
              'sc1233']

paramModel = (np.size(out_names))*(hru_num)
hru_names =[]
for i in out_names:
    hru_names.append(['{}{}'.format(j, i) for j in hruidxID])
hru_names1 = np.reshape(hru_names,(paramModel,1))
hru_names_df = pd.DataFrame (hru_names1)
#%% reading output_swe files
av_ncfiles = ["SA2/sa_sa2_lj1110_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj1120_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj1130_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj1210_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj1220_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj1230_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj1310_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj1320_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj1330_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj2110_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj2120_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj2130_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj2210_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj2220_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj2230_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj2310_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj2320_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lj2330_2006-2007_senatorVariableDecayRate_1.nc",
#              
#              "SA2/sa_sa2_mj1110_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj1120_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj1130_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj1210_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj1220_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj1230_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj1310_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj1320_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj1330_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj2110_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj2120_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj2130_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj2210_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj2220_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj2230_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj2310_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj2320_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_mj2330_2006-2007_senatorVariableDecayRate_1.nc",
#              
#              "SA2/sa_sa2_sj1110_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj1120_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj1130_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj1210_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj1220_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj1230_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj1310_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj1320_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj1330_2006-2007_senatorVariableDecayRate_1.nc",
##              "SA2/sa_sa2_sj2110_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj2120_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj2130_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj2210_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj2220_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj2230_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj2310_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj2320_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sj2330_2006-2007_senatorVariableDecayRate_1.nc",
#              
#              "SA2/sa_sa2_lc1111_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1112_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1113_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1121_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1122_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1123_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1131_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1132_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1133_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1211_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1212_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1213_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1221_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1222_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1223_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1231_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1232_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1233_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1311_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1312_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1313_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1321_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1322_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1323_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1331_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1332_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc1333_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc2111_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc2112_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc2113_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc2121_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc2122_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc2123_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc2131_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc2132_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc2133_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_lc2211_2006-2007_senatorVariableDecayRate_1.nc",
#              
#              "SA2/sa_sa2_sc1111_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1112_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1113_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1121_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1122_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1123_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1131_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1132_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1133_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1211_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1212_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1213_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1221_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1222_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1223_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1231_2006-2007_senatorVariableDecayRate_1.nc",
#              "SA2/sa_sa2_sc1232_2006-2007_senatorVariableDecayRate_1.nc",
              "SA2/sa_sa2_sc1233_2006-2007_senatorVariableDecayRate_1.nc"
              ]

av_all = []
for ncfiles in av_ncfiles:
    av_all.append(Dataset(ncfiles))

av_sd = []
for dfssd in av_all:
    av_sd.append(pd.DataFrame(dfssd['scalarSnowDepth'][:]))
av_sd_df = pd.concat (av_sd, axis=1)
av_sd_df.columns =  hru_names_df[0]

# extract all swe results
av_swe = []
for dfsswe in av_all:
    av_swe.append(pd.DataFrame(dfsswe['scalarSWE'][:]))
av_swe_df = pd.concat (av_swe, axis=1)
av_swe_df.columns =  hru_names_df[0]

#merge two years of swe data in a dataframe 'scalarSWE'
#av_swe_lst = []
#numfile = 0
#while numfile < (len(av_ncfiles)):
#    av_swe_lst.append(pd.concat ([av_swe[numfile],av_swe[numfile+1]], axis=0, ignore_index=True))
#    numfile = numfile+2
#av_swe_df = []
#for numfil in range (len(av_swe)):
#    for hrus in range (hru_num):
#        av_swe_df.append(np.array(av_swe[numfil][hrus]))
#av_swe_df = pd.DataFrame(np.array(av_swe_df).T)
#av_swe_df.columns =  hru_names_df[0]
#%%
for varname in av_all[0].variables.keys():
    var = av_all[0].variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)
#%% output time step
TimeSa2006 = av_all[0].variables['time'][:] # get values
#TimeSa2007 = av_all[1].variables['time'][:] # get values
#TimeSa = np.concatenate((TimeSa2006, TimeSa2007), axis=0)
t_unitSa = av_all[0].variables['time'].units # get unit  "days since 1950-01-01T00:00:00Z"

try :

    t_cal = av_all[0].variables['time'].calendar

except AttributeError : # Attribute doesn't exist

    t_cal = u"gregorian" # or standard

#tvalueSa = num2date(TimeSa, units=t_unitSa, calendar=t_cal)
tvalueSa = num2date(TimeSa2006, units=t_unitSa, calendar=t_cal)
DateSa = [i.strftime("%Y-%m-%d %H:%M") for i in tvalueSa] # -%d %H:%M to display dates as string #i.strftime("%Y-%m-%d %H:%M")  

#%% defining counter
counter = pd.DataFrame(np.arange(0,np.size(av_sd_df['11111lj1110'])),columns=['counter'])
counter.set_index(pd.DatetimeIndex(DateSa),inplace=True)
#%%
av_swe_df.set_index(pd.DatetimeIndex(DateSa),inplace=True)
av_swe_df2 = pd.concat([counter, av_swe_df], axis=1)

av_sd_df.set_index(pd.DatetimeIndex(DateSa),inplace=True)
av_sd_df2 = pd.concat([counter, av_sd_df], axis=1)
#%% day of snow disappearance (based on snowdepth)-final output
av_sd_df5000 = av_sd_df2[:][5000:8737]

zerosnowdate = []
for val in hru_names_df[0]:
    zerosnowdate.append(np.where(av_sd_df5000[val]==0))
zerosnowdate_omg = [item[0] for item in zerosnowdate] #change tuple to array
for i,item in enumerate(zerosnowdate_omg):
    if len(item) == 0:
        zerosnowdate_omg[i] = 3737
for i,item in enumerate(zerosnowdate_omg):
    zerosnowdate_omg[i] = zerosnowdate_omg[i]+5000
        
first_zerosnowdate =[]
for i,item in enumerate(zerosnowdate_omg):
    if np.size(item)>1:
        #print np.size(item)
        first_zerosnowdate.append(item[0])
    if np.size(item)==1:
        first_zerosnowdate.append(item)
    
#first_zerosnowdate_df = pd.DataFrame(np.reshape(first_zerosnowdate, ((np.size(hru_names1)),0)).T, columns=out_names)
dosd_df = pd.DataFrame(np.array(first_zerosnowdate)).T
dosd_df.columns = hru_names_df[0]
#first_zerosnowdate_df_obs = pd.DataFrame(np.array([[5985],[6200]]).T,columns=out_names)
dosd_obs = pd.DataFrame(np.array([5985]),columns=['2006'])

dosd_residual=[]
for hru in dosd_df.columns:
    dosd_residual.append((dosd_df[hru][0]-dosd_obs['2006'])/24)

dosd_residual_df = pd.DataFrame(np.reshape(np.array(dosd_residual),(np.size(out_names),hru_num)).T, columns=out_names)

#%%
#plt.xticks(x, hru[::3], rotation=25)
#for namefile in out_names:
#    x = list(np.arange(1,244))
#    fig = plt.figure(figsize=(20,15))
#    plt.bar(x,zerosnowdate_residual_df[namefile])
#    plt.title(namefile, fontsize=42)
#    plt.xlabel('hrus',fontsize=30)
#    plt.ylabel('residual dosd (day)', fontsize=30)
#    #vax.yaxis.set_label_coords(0.5, -0.1) 
#    plt.savefig('SA2/'+namefile)

#%%out_names = ['AvInitial','Avas2l','Avas3m','Avtc2s','Avns3p','Avns4c']# 'Avwp2e','Avas2l','Avas3m','Avtc3t','Avtc4m','Avns3p','Avce2s','Avtc2s','Avtc3t','Avtc4m','Avns2a','Avns3p','Avns4c']
swe_obs2006 = pd.DataFrame (obs_swe['swe_mm']['2006-11-01':'2007-06-06'], columns=['swe_mm'])
swe_obs2007 = pd.DataFrame (obs_swe['swe_mm']['2007-12-03':'2008-06-08'], columns=['swe_mm'])
date_swe2006 = ['2006-11-01 11:10','2006-11-30 12:30','2007-01-01 11:10','2007-01-30 10:35','2007-03-05 14:30', 
                '2007-03-12 14:00','2007-03-19 12:30','2007-03-26 12:30','2007-04-02 12:30','2007-04-18 08:35', 
                '2007-04-23 10:30','2007-05-02 08:40','2007-05-09 08:50','2007-05-16 09:00','2007-05-23 08:30',
                '2007-05-30 09:00','2007-06-06 08:15']
swe_obs2006.set_index(pd.DatetimeIndex(date_swe2006),inplace=True)

#%%
#DateSa2 = [i.strftime("%Y-%m-%d") for i in tvalueSa]
#sax = np.arange(0,np.size(DateSa2))
#sa_xticks = DateSa2
#safig, saax = plt.subplots(1,1, figsize=(20,15))
#plt.xticks(sax, sa_xticks[::1000], rotation=25, fontsize=20)
#saax.xaxis.set_major_locator(ticker.AutoLocator())
#plt.yticks(fontsize=20)
#for hru in hru_names_df[0]:
#    plt.plot(av_swe_df[hru])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]
#
#plt.plot(swe_obs2006, 'ok', markersize=10)
#
#plt.title('rainbow_SWE', position=(0.04, 0.88), ha='left', fontsize=40)
#plt.xlabel('Time 2005-2006', fontsize=30)
#plt.ylabel('SWE(mm)', fontsize=30)
##plt.show()
#plt.savefig('SA2/swelj.png')
#%%**************************************************************************************************
# *********************** finding max corespondance swe for '2007-05-09 08:50'***********************
#'2007-04-18' 4776: 4800, '2007-04-23' 4896:4920, '2007-05-02' 5112:5136
#Group1: '2007-03-12 14:00' (3902),'2007-03-19 12:30 (4068)','2007-03-26 12:30 (4236)','2007-04-02 12:30'(4404),
#Group2: '2007-04-18 08:35' (4784),'2007-04-23 10:30 (4907)','2007-05-02 08:40'(5121), 
maxSWE = []
for names in hru_names_df[0]:
    maxSWE.append(av_swe_df[names][5289])

#%%**************************************************************************************************
# ********************** calculating snowmelt rate based on SWE *************************************
minSWE = []
minSWEdate = []
for names in hru_names_df[0]:
    if av_swe_df[names][5960]>0:
        minSWE.append(av_swe_df[names][5960])
        minSWEdate.append(5960)
    else: 
        minSWE.append(float(av_swe_df[names][dosd_df[names]-1]))
        minSWEdate.append(float(dosd_df[names]-1))

mdeltaday = []
mdeltaSWE = []
meltingrate = [] #cm/day
for counterhd in range (np.size(minSWE)):
    mdeltaday.append(float(minSWEdate[counterhd]-5289))
    mdeltaSWE.append(float(maxSWE[counterhd]-minSWE[counterhd]))
    meltingrate.append(float(0.1*24*mdeltaSWE[counterhd]/mdeltaday[counterhd]))
#%% **************************************************************************************************
# ************************** calculating cold content ************************************************
nlayerTemp = []
for dfstemp in av_all:
    nlayerTemp.append(pd.DataFrame(dfstemp['mLayerTemp'][:][:]))
        
#%% group1 number of snowlayer 
nsnow0312 = []
nsnow0319 = []
nsnow0326 = []
nsnow0402 = []

nsnow0418 = []
nsnow0423 = []
nsnow0502 = []
for dfsg1 in av_all:
    nsnow0312.append(pd.DataFrame(dfsg1['nSnow'][:][3902]))
    nsnow0319.append(pd.DataFrame(dfsg1['nSnow'][:][4068]))
    nsnow0326.append(pd.DataFrame(dfsg1['nSnow'][:][4236]))
    nsnow0402.append(pd.DataFrame(dfsg1['nSnow'][:][4404]))
    #group2 number of snowlayer 
    nsnow0418.append(pd.DataFrame(dfsg1['nSnow'][:][4784]))
    nsnow0423.append(pd.DataFrame(dfsg1['nSnow'][:][4907]))
    nsnow0502.append(pd.DataFrame(dfsg1['nSnow'][:][5121]))
    
sumlayer0312 = []
sumlayer0319 = []
sumlayer0326 = []
sumlayer0402 = []
sumlayer0418 = []
sumlayer0423 = []
sumlayer0502 = []
for dfsg1s in av_all:
    sumlayer0312.append(pd.DataFrame(sum(dfsg1s['nLayers'][:][0:3902])))
    sumlayer0319.append(pd.DataFrame(sum(dfsg1s['nLayers'][:][0:4068])))
    sumlayer0326.append(pd.DataFrame(sum(dfsg1s['nLayers'][:][0:4236])))
    sumlayer0402.append(pd.DataFrame(sum(dfsg1s['nLayers'][:][0:4404])))
    sumlayer0418.append(pd.DataFrame(sum(dfsg1s['nLayers'][:][0:4784])))
    sumlayer0423.append(pd.DataFrame(sum(dfsg1s['nLayers'][:][0:4907])))
    sumlayer0502.append(pd.DataFrame(sum(dfsg1s['nLayers'][:][0:5121])))      
#%% #finding snow layers temperature
#group1
snowlayertemp0312 = []
snowlayertemp0319 = []
snowlayertemp0326 = []
snowlayertemp0402 = []
#group2
snowlayertemp0418 = []
snowlayertemp0423 = []
snowlayertemp0502 = []
for nfslt in range (np.size(av_all)):
    for hruslt in range (hru_num):
        snowlayertemp0312.append(list(nlayerTemp[nfslt][hruslt][sumlayer0312[nfslt][0][hruslt]:sumlayer0312[nfslt][0][hruslt]+nsnow0312[nfslt][0][hruslt]]))
        snowlayertemp0319.append(list(nlayerTemp[nfslt][hruslt][sumlayer0319[nfslt][0][hruslt]:sumlayer0319[nfslt][0][hruslt]+nsnow0319[nfslt][0][hruslt]]))
        snowlayertemp0326.append(list(nlayerTemp[nfslt][hruslt][sumlayer0326[nfslt][0][hruslt]:sumlayer0326[nfslt][0][hruslt]+nsnow0326[nfslt][0][hruslt]]))
        snowlayertemp0402.append(list(nlayerTemp[nfslt][hruslt][sumlayer0402[nfslt][0][hruslt]:sumlayer0402[nfslt][0][hruslt]+nsnow0402[nfslt][0][hruslt]]))

        snowlayertemp0418.append(list(nlayerTemp[nfslt][hruslt][sumlayer0418[nfslt][0][hruslt]:sumlayer0418[nfslt][0][hruslt]+nsnow0418[nfslt][0][hruslt]]))
        snowlayertemp0423.append(list(nlayerTemp[nfslt][hruslt][sumlayer0423[nfslt][0][hruslt]:sumlayer0423[nfslt][0][hruslt]+nsnow0423[nfslt][0][hruslt]]))
        snowlayertemp0502.append(list(nlayerTemp[nfslt][hruslt][sumlayer0502[nfslt][0][hruslt]:sumlayer0502[nfslt][0][hruslt]+nsnow0502[nfslt][0][hruslt]]))

#%% finding ice and liquid fraction of snow layers and height of each snow layer
#nporespace = []
nvolfracIce = []
nvolfracliq = []
nheight = []
for df in av_all:
    nvolfracIce.append(pd.DataFrame(df['mLayerVolFracIce'][:])) 
    nvolfracliq.append(pd.DataFrame(df['mLayerVolFracLiq'][:]))
    nheight.append(pd.DataFrame(df['mLayerHeight'][:]))
    #nporespace.append(pd.DataFrame(df['mLayerPoreSpace'][:])) 

#volumetric fraction of ice in snow layers
#group1
volfracIce0312 = []
volfracIce0319 = []
volfracIce0326 = []
volfracIce0402 = []
#group2
volfracIce0418 = []
volfracIce0423 = []
volfracIce0502 = []
for nfslt in range (np.size(av_all)):
    for hruslt in range (hru_num):
        volfracIce0312.append(list(nvolfracIce[nfslt][hruslt][sumlayer0312[nfslt][0][hruslt]:sumlayer0312[nfslt][0][hruslt]+nsnow0312[nfslt][0][hruslt]]))
        volfracIce0319.append(list(nvolfracIce[nfslt][hruslt][sumlayer0319[nfslt][0][hruslt]:sumlayer0319[nfslt][0][hruslt]+nsnow0319[nfslt][0][hruslt]]))
        volfracIce0326.append(list(nvolfracIce[nfslt][hruslt][sumlayer0326[nfslt][0][hruslt]:sumlayer0326[nfslt][0][hruslt]+nsnow0326[nfslt][0][hruslt]]))
        volfracIce0402.append(list(nvolfracIce[nfslt][hruslt][sumlayer0402[nfslt][0][hruslt]:sumlayer0402[nfslt][0][hruslt]+nsnow0402[nfslt][0][hruslt]]))
        volfracIce0418.append(list(nvolfracIce[nfslt][hruslt][sumlayer0418[nfslt][0][hruslt]:sumlayer0418[nfslt][0][hruslt]+nsnow0418[nfslt][0][hruslt]]))
        volfracIce0423.append(list(nvolfracIce[nfslt][hruslt][sumlayer0423[nfslt][0][hruslt]:sumlayer0423[nfslt][0][hruslt]+nsnow0423[nfslt][0][hruslt]]))
        volfracIce0502.append(list(nvolfracIce[nfslt][hruslt][sumlayer0502[nfslt][0][hruslt]:sumlayer0502[nfslt][0][hruslt]+nsnow0502[nfslt][0][hruslt]]))
#volumetric fraction of liquid in snow layers
#group1
volfracliq0312 = []
volfracliq0319 = []
volfracliq0326 = []
volfracliq0402 = []
#group2
volfracliq0418 = []
volfracliq0423 = []
volfracliq0502 = []
for nfslt in range (np.size(av_all)):
    for hruslt in range (hru_num):
        volfracliq0312.append(list(nvolfracliq[nfslt][hruslt][sumlayer0312[nfslt][0][hruslt]:sumlayer0312[nfslt][0][hruslt]+nsnow0312[nfslt][0][hruslt]]))
        volfracliq0319.append(list(nvolfracliq[nfslt][hruslt][sumlayer0319[nfslt][0][hruslt]:sumlayer0319[nfslt][0][hruslt]+nsnow0319[nfslt][0][hruslt]]))
        volfracliq0326.append(list(nvolfracliq[nfslt][hruslt][sumlayer0326[nfslt][0][hruslt]:sumlayer0326[nfslt][0][hruslt]+nsnow0326[nfslt][0][hruslt]]))
        volfracliq0402.append(list(nvolfracliq[nfslt][hruslt][sumlayer0402[nfslt][0][hruslt]:sumlayer0402[nfslt][0][hruslt]+nsnow0402[nfslt][0][hruslt]]))
        volfracliq0418.append(list(nvolfracliq[nfslt][hruslt][sumlayer0418[nfslt][0][hruslt]:sumlayer0418[nfslt][0][hruslt]+nsnow0418[nfslt][0][hruslt]]))
        volfracliq0423.append(list(nvolfracliq[nfslt][hruslt][sumlayer0423[nfslt][0][hruslt]:sumlayer0423[nfslt][0][hruslt]+nsnow0423[nfslt][0][hruslt]]))
        volfracliq0502.append(list(nvolfracliq[nfslt][hruslt][sumlayer0502[nfslt][0][hruslt]:sumlayer0502[nfslt][0][hruslt]+nsnow0502[nfslt][0][hruslt]]))
# height of each snow layer
#group1
height0312 = []
height0319 = []
height0326 = []
height0402 = []
#group2
height0418 = []
height0423 = []
height0502 = []
for nfslt in range (np.size(av_all)):
    for hruslt in range (hru_num):
        height0312.append(list(nheight[nfslt][hruslt][sumlayer0312[nfslt][0][hruslt]:sumlayer0312[nfslt][0][hruslt]+nsnow0312[nfslt][0][hruslt]]))
        height0319.append(list(nheight[nfslt][hruslt][sumlayer0319[nfslt][0][hruslt]:sumlayer0319[nfslt][0][hruslt]+nsnow0319[nfslt][0][hruslt]]))
        height0326.append(list(nheight[nfslt][hruslt][sumlayer0326[nfslt][0][hruslt]:sumlayer0326[nfslt][0][hruslt]+nsnow0326[nfslt][0][hruslt]]))
        height0402.append(list(nheight[nfslt][hruslt][sumlayer0402[nfslt][0][hruslt]:sumlayer0402[nfslt][0][hruslt]+nsnow0402[nfslt][0][hruslt]]))
        height0418.append(list(nheight[nfslt][hruslt][sumlayer0418[nfslt][0][hruslt]:sumlayer0418[nfslt][0][hruslt]+nsnow0418[nfslt][0][hruslt]]))
        height0423.append(list(nheight[nfslt][hruslt][sumlayer0423[nfslt][0][hruslt]:sumlayer0423[nfslt][0][hruslt]+nsnow0423[nfslt][0][hruslt]]))
        height0502.append(list(nheight[nfslt][hruslt][sumlayer0502[nfslt][0][hruslt]:sumlayer0502[nfslt][0][hruslt]+nsnow0502[nfslt][0][hruslt]]))
    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
    
    
    
    
    
    
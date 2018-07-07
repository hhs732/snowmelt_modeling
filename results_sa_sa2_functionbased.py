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
from allNcFiles import av_ncfiles
#%% output file

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

#%% defining hrus_name
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
#%% defining functions
def readAllNcfilesAsDataset(allNcfiles):
    allNcfilesDataset = []
    for ncfiles in allNcfiles:
        allNcfilesDataset.append(Dataset(ncfiles))
    return allNcfilesDataset
#hruname = hru_names_df[0]
def readVariablefromNcfilesDatasetasDF(NcfilesDataset,variable,hruname):
    variableNameList = []
    for datasets in NcfilesDataset:
        variableNameList.append(pd.DataFrame(datasets[variable][:][:]))
    variableNameDF = pd.concat (variableNameList, axis=1)
    variableNameDF.columns = hruname
    counter = pd.DataFrame(np.arange(0,np.size(variableNameDF[hru_names_df[0][0]])),columns=['counter'])
    counter.set_index(variableNameDF.index,inplace=True)
    variableNameDF = pd.concat([counter, variableNameDF], axis=1)
    return variableNameDF

def date(allNcfilesDataset,formatDate):
    Time = allNcfilesDataset[0].variables['time'][:] 
    #TimeSa = np.concatenate((TimeSa2006, TimeSa2007), axis=0)
    t_unit = allNcfilesDataset[0].variables['time'].units 
    
    try :
        t_cal = allNcfilesDataset[0].variables['time'].calendar
    
    except AttributeError : # Attribute doesn't exist error
        t_cal = u"gregorian" # or standard
    #tvalueSa = num2date(TimeSa, units=t_unitSa, calendar=t_cal)
    tvalue = num2date(Time, units=t_unit, calendar=t_cal)
    Date = [i.strftime(formatDate) for i in tvalue] # "%Y-%m-%d %H:%M"
    return Date

def readSpecificDatafromAllHRUs(variablename,hruname,day):
    dayData = []
    for names in hruname:
        dayData.append(variablename[names][day])
    return dayData

def sumBeforeSpecificDatafromAllHRUs(variablename,hruname,day):
    sumData = []
    for names in hruname:
        sumData.append(sum(variablename[names][0:day]))
    return sumData   

def snowLayerAttributeforSpecificDate(layerattributefile,hruname,sumlayer,snowlayer): #like snowlayetemp, volFracosIce, ....
    snowlayerattribute = []
    for names in hruname:
        snowlayerattribute.append(list(layerattributefile[names][sumlayer[names][0]:sumlayer[names][0]+snowlayer[names][0]]))
    return snowlayerattribute
#%%  reading output files
av_all = readAllNcfilesAsDataset(av_ncfiles)
DateSa = date(av_all,"%Y-%m-%d %H:%M")

av_sd_df = readVariablefromNcfilesDatasetasDF(av_all,'scalarSnowDepth',hru_names_df[0])
av_sd_df.set_index(pd.DatetimeIndex(DateSa),inplace=True)
av_swe_df =  readVariablefromNcfilesDatasetasDF(av_all,'scalarSWE',hru_names_df[0])
av_swe_df.set_index(pd.DatetimeIndex(DateSa),inplace=True)

nvolfracIce = readVariablefromNcfilesDatasetasDF(av_all,'mLayerVolFracIce',hru_names_df[0])
nvolfracliq = readVariablefromNcfilesDatasetasDF(av_all,'mLayerVolFracLiq',hru_names_df[0])
nheight = readVariablefromNcfilesDatasetasDF(av_all,'mLayerHeight',hru_names_df[0])

#%%
for varname in av_all[0].variables.keys():
    var = av_all[0].variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)
#%% day of snow disappearance (based on snowdepth)-final output
av_sd_df5000 = av_sd_df[:][5000:8737]

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

maxSWE = readSpecificDatafromAllHRUs(av_swe_df,hru_names_df[0],5289)

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
nlayerTemp =  readVariablefromNcfilesDatasetasDF(av_all,'mLayerTemp',hru_names_df[0])
nsnow =  readVariablefromNcfilesDatasetasDF(av_all,'nSnow',hru_names_df[0])
nlayer = readVariablefromNcfilesDatasetasDF(av_all,'nLayers',hru_names_df[0])

#%% number of snowlayer 
nsnow0312 = pd.DataFrame(readSpecificDatafromAllHRUs(nsnow,hru_names_df[0],3902)).T; nsnow0312.columns = hru_names_df[0]
nsnow0319 = pd.DataFrame(readSpecificDatafromAllHRUs(nsnow,hru_names_df[0],4068)).T; nsnow0319.columns = hru_names_df[0]
nsnow0326 = pd.DataFrame(readSpecificDatafromAllHRUs(nsnow,hru_names_df[0],4236)).T; nsnow0326.columns = hru_names_df[0]
nsnow0402 = pd.DataFrame(readSpecificDatafromAllHRUs(nsnow,hru_names_df[0],4404)).T; nsnow0402.columns = hru_names_df[0]

nsnow0418 = pd.DataFrame(readSpecificDatafromAllHRUs(nsnow,hru_names_df[0],4784)).T; nsnow0418.columns = hru_names_df[0]
nsnow0423 = pd.DataFrame(readSpecificDatafromAllHRUs(nsnow,hru_names_df[0],4907)).T; nsnow0423.columns = hru_names_df[0]
nsnow0502 = pd.DataFrame(readSpecificDatafromAllHRUs(nsnow,hru_names_df[0],5121)).T; nsnow0502.columns = hru_names_df[0]

# sum of all layers befor target layer
sumlayer0312 = pd.DataFrame(sumBeforeSpecificDatafromAllHRUs(nlayer,hru_names_df[0],3902)).T; sumlayer0312.columns = hru_names_df[0]
sumlayer0319 = pd.DataFrame(sumBeforeSpecificDatafromAllHRUs(nlayer,hru_names_df[0],4068)).T; sumlayer0319.columns = hru_names_df[0]
sumlayer0326 = pd.DataFrame(sumBeforeSpecificDatafromAllHRUs(nlayer,hru_names_df[0],4236)).T; sumlayer0326.columns = hru_names_df[0]
sumlayer0402 = pd.DataFrame(sumBeforeSpecificDatafromAllHRUs(nlayer,hru_names_df[0],4404)).T; sumlayer0402.columns = hru_names_df[0]
sumlayer0418 = pd.DataFrame(sumBeforeSpecificDatafromAllHRUs(nlayer,hru_names_df[0],4784)).T; sumlayer0418.columns = hru_names_df[0]
sumlayer0423 = pd.DataFrame(sumBeforeSpecificDatafromAllHRUs(nlayer,hru_names_df[0],4907)).T; sumlayer0423.columns = hru_names_df[0]
sumlayer0502 = pd.DataFrame(sumBeforeSpecificDatafromAllHRUs(nlayer,hru_names_df[0],5121)).T; sumlayer0502.columns = hru_names_df[0]

#%%snow layer temperature
##group1
snowlayertemp0312 = snowLayerAttributeforSpecificDate(nlayerTemp,hru_names_df[0],sumlayer0312,nsnow0312)
snowlayertemp0319 = snowLayerAttributeforSpecificDate(nlayerTemp,hru_names_df[0],sumlayer0319,nsnow0319)
snowlayertemp0326 = snowLayerAttributeforSpecificDate(nlayerTemp,hru_names_df[0],sumlayer0326,nsnow0326)
snowlayertemp0402 = snowLayerAttributeforSpecificDate(nlayerTemp,hru_names_df[0],sumlayer0402,nsnow0402)
#group2
snowlayertemp0418 = snowLayerAttributeforSpecificDate(nlayerTemp,hru_names_df[0],sumlayer0418,nsnow0418)
snowlayertemp0423 = snowLayerAttributeforSpecificDate(nlayerTemp,hru_names_df[0],sumlayer0423,nsnow0423)
snowlayertemp0502 = snowLayerAttributeforSpecificDate(nlayerTemp,hru_names_df[0],sumlayer0502,nsnow0502)

#%% volumetric fraction of ice in snow layers
#group1
volfracIce0312 = snowLayerAttributeforSpecificDate(nvolfracIce,hru_names_df[0],sumlayer0312,nsnow0312)
volfracIce0319 = snowLayerAttributeforSpecificDate(nvolfracIce,hru_names_df[0],sumlayer0319,nsnow0319)
volfracIce0326 = snowLayerAttributeforSpecificDate(nvolfracIce,hru_names_df[0],sumlayer0326,nsnow0326)
volfracIce0402 = snowLayerAttributeforSpecificDate(nvolfracIce,hru_names_df[0],sumlayer0402,nsnow0402)
#group2
volfracIce0418 = snowLayerAttributeforSpecificDate(nvolfracIce,hru_names_df[0],sumlayer0418,nsnow0418)
volfracIce0423 = snowLayerAttributeforSpecificDate(nvolfracIce,hru_names_df[0],sumlayer0423,nsnow0423)
volfracIce0502 = snowLayerAttributeforSpecificDate(nvolfracIce,hru_names_df[0],sumlayer0502,nsnow0502)

#volumetric fraction of liquid in snow layers
#group1
volfracLiq0312 = snowLayerAttributeforSpecificDate(nvolfracliq,hru_names_df[0],sumlayer0312,nsnow0312)
volfracLiq0319 = snowLayerAttributeforSpecificDate(nvolfracliq,hru_names_df[0],sumlayer0319,nsnow0319)
volfracLiq0326 = snowLayerAttributeforSpecificDate(nvolfracliq,hru_names_df[0],sumlayer0326,nsnow0326)
volfracLiq0402 = snowLayerAttributeforSpecificDate(nvolfracliq,hru_names_df[0],sumlayer0402,nsnow0402)
#group2
volfracLiq0418 = snowLayerAttributeforSpecificDate(nvolfracliq,hru_names_df[0],sumlayer0418,nsnow0418)
volfracLiq0423 = snowLayerAttributeforSpecificDate(nvolfracliq,hru_names_df[0],sumlayer0423,nsnow0423)
volfracLiq0502 = snowLayerAttributeforSpecificDate(nvolfracliq,hru_names_df[0],sumlayer0502,nsnow0502)

# height of each snow layer
#group1
height0312 = snowLayerAttributeforSpecificDate(nheight,hru_names_df[0],sumlayer0312,nsnow0312)
height0319 = snowLayerAttributeforSpecificDate(nheight,hru_names_df[0],sumlayer0319,nsnow0319)
height0326 = snowLayerAttributeforSpecificDate(nheight,hru_names_df[0],sumlayer0326,nsnow0326)
height0402 = snowLayerAttributeforSpecificDate(nheight,hru_names_df[0],sumlayer0402,nsnow0402)
#group2
height0418 = snowLayerAttributeforSpecificDate(nheight,hru_names_df[0],sumlayer0418,nsnow0418)
height0423 = snowLayerAttributeforSpecificDate(nheight,hru_names_df[0],sumlayer0423,nsnow0423)
height0502 = snowLayerAttributeforSpecificDate(nheight,hru_names_df[0],sumlayer0502,nsnow0502)
#%%
def depthOfLayers(heightfile):
    finalHeight = []
    sumSofar=0
    for lsts in heightfile:
        lsts.reverse()
        height_ls = []
        for values in lsts:
            height=2*(abs(values)-sumSofar)
            height_ls.append(height)
            sumSofar+=height
        finalHeight.append(height_ls)
    return finalHeight
    
#%%
height0312rt=depthOfLayers(height0312)
 
 
    
    
    
    
    
    
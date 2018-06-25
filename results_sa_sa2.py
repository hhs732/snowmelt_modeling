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
date_swe = ['2006-11-01', '2006-11-30', '2007-01-01', '2007-01-30', '2007-03-05', '2007-03-12', 
            '2007-03-19', '2007-03-26', '2007-04-02', '2007-04-18', '2007-04-23', '2007-05-02', 
            '2007-05-09', '2007-05-16', '2007-05-23', '2007-05-30', '2007-06-06', 
            
            '2007-12-03', '2008-01-01', '2008-01-31', '2008-03-03', '2008-03-24', '2008-04-01', 
            '2008-04-14', '2008-04-22', '2008-04-28', '2008-05-06', '2008-05-12', '2008-05-19',
            '2008-05-26', '2008-06-02', '2008-06-08'] 
            
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
#             'lj1120','lj1130','lj1210','lj1220','lj1230','lj1310','lj1320','lj1330','lj2110','lj2120','lj2130','lj2210','lj2220','lj2230','lj2310','lj2320','lj2330', 
#             'mj1110','mj1120','mj1130','mj1210','mj1220','mj1230','mj1310','mj1320','mj1330','mj2110','mj2120','mj2130','mj2210','mj2220','mj2230','mj2310','mj2320','mj2330',
#             'sj1110','sj1120','sj1130','sj1210',
             'sj1220']

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
              "SA2/sa_sa2_sj1220_2006-2007_senatorVariableDecayRate_1.nc"
#              "SA2/sa_sa2_sj1230_2006-2007_senatorVariableDecayRate_1.nc",
              ] 
av_all = []
for ncfiles in av_ncfiles:
    av_all.append(Dataset(ncfiles))

av_sd = []
for dfs in av_all:
    av_sd.append(pd.DataFrame(dfs['scalarSnowDepth'][:]))
av_sd_df = pd.concat (av_sd, axis=1)
av_sd_df.columns =  hru_names_df[0]
# extract all swe results
av_swe = []
for dfs in av_all:
    av_swe.append(pd.DataFrame(dfs['scalarSWE'][:]))
#merge two years of swe data in a dataframe 
#av_swe_lst = []
#numfile = 0
#while numfile < (len(av_ncfiles)):
#    av_swe_lst.append(pd.concat ([av_swe[numfile],av_swe[numfile+1]], axis=0, ignore_index=True))
#    numfile = numfile+2

av_swe_df = []
for numfil in range (len(av_swe)):
    for hrus in range (hru_num):
        av_swe_df.append(av_swe[numfil][hrus])
av_swe_df = pd.DataFrame(np.array(av_swe_df).T, columns =  hru_names_df)
#%%
#for varname in av_all[10].variables.keys():
#    var = av_all[10].variables[varname]
#    print (varname, var.dtype, var.dimensions, var.shape)
    
#nsnow = []
#for dfs in av_all:
#    nsnow.append(pd.DataFrame(dfs['nSnow'][:]))
#    
#nlayer = []
#for dfs in av_all:
#    nlayer.append(pd.DataFrame(dfs['nLayers'][:]))

#'2007-04-18' 4776: 4800, '2007-04-23' 4896:4920, '2007-05-02' 5112:5136
nsnow = []
for dfs in av_all:
    nsnow.append(pd.DataFrame(dfs['nSnow'][:][4776:4800]))
    nsnow.append(pd.DataFrame(dfs['nSnow'][:][4896:4920]))
    nsnow.append(pd.DataFrame(dfs['nSnow'][:][5112:5136]))

nlayer = []
for dfs in av_all:
    nlayer.append(pd.DataFrame(dfs['nLayers'][:][4776:4800]))
    nlayer.append(pd.DataFrame(dfs['nLayers'][:][4896:4920]))
    nlayer.append(pd.DataFrame(dfs['nLayers'][:][5112:5136]))

sumlayer1 = []
sumlayer2 = []
for dfs in av_all:
    sumlayer1.append(pd.DataFrame(sum(dfs['nLayers'][:][0:4776])))
    sumlayer1.append(pd.DataFrame(sum(dfs['nLayers'][:][0:4896])))
    sumlayer1.append(pd.DataFrame(sum(dfs['nLayers'][:][0:5112])))
    sumlayer2.append(pd.DataFrame(sum(dfs['nLayers'][:][0:4800])))
    sumlayer2.append(pd.DataFrame(sum(dfs['nLayers'][:][0:4920])))
    sumlayer2.append(pd.DataFrame(sum(dfs['nLayers'][:][0:5136])))      

nlayerTemp = []
for dfs in av_all:
    nlayerTemp.append(pd.DataFrame(dfs['mLayerTemp'][:][:]))
#%%finding all layers temperature
nlayertemp0418 = []
nlayertemp0423 = []
nlayertemp0502 = []
nf = 0
while nf < (3*np.size(av_ncfiles)):
    for hru in range (hru_num):
        nlayertemp0418.append(nlayerTemp[nf/3][hru][sumlayer1[nf][0][hru]:sumlayer2[nf][0][hru]])
        nlayertemp0423.append(nlayerTemp[nf/3][hru][sumlayer1[nf+1][0][hru]:sumlayer2[nf+1][0][hru]])
        nlayertemp0502.append(nlayerTemp[nf/3][hru][sumlayer1[nf+2][0][hru]:sumlayer2[nf+2][0][hru]])
    nf = nf+3
#%% finding snow layer temperature
nlayertemp0418_rsh = np.reshape(nlayertemp0418,(np.size(av_ncfiles),hru_num))
nlayertemp0423_rsh = np.reshape(nlayertemp0423,(np.size(av_ncfiles),hru_num))
nlayertemp0502_rsh = np.reshape(nlayertemp0502,(np.size(av_ncfiles),hru_num))
#nltcheck = nlayertemp0418_rsh[1][242]
#%%
nsnowlayertemp0418 = []
nsnowlayertemp0423 = []
nsnowlayertemp0502 = []
for nf in range (np.size(av_ncfiles)):
    for hru in range (hru_num):
        st = 0
        for tm in range (np.size(nsnow[0][0])):
            nsnowlayertemp0418.append(nlayertemp0418_rsh[nf][hru][st:st+nsnow[nf][hru][tm]])
            nsnowlayertemp0423.append(nlayertemp0423_rsh[nf][hru][st:st+nsnow[nf][hru][tm]])
            nsnowlayertemp0502.append(nlayertemp0502_rsh[nf][hru][st:st+nsnow[nf][hru][tm]])
            st = st+nlayer[nf][hru][tm]

#%%
#snowlayertemp_df = pd.DataFrame(np.reshape(snowlayertemp, ((np.size(hru_names1)),np.size(nsnow[0][0]))).T, columns=hru_names_df)
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

av_swe_df.set_index(pd.DatetimeIndex(DateSa),inplace=True)
#%%
#TimeSa2006 = av_all[0].variables['time'][:] # get values
#tvalueSa2006 = num2date(TimeSa2006, units=t_unitSa, calendar=t_cal)
#DateSa2006 = [i.strftime("%Y-%m") for i in tvalueSa2006]
#
#TimeSa2007 = av_all[1].variables['time'][:] # get values
#tvalueSa2007 = num2date(TimeSa2007, units=t_unitSa, calendar=t_cal)
#DateSa2007 = [i.strftime("%Y-%m") for i in tvalueSa2007]
#%% day of snow disappearance-final output
av_sd_df.set_index(pd.DatetimeIndex(DateSa),inplace=True)
counter = pd.DataFrame(np.arange(0,np.size(av_sd_df['11111lj1110'])),columns=['counter'])
counter.set_index(av_sd_df.index,inplace=True)
av_sd_df2 = pd.concat([counter, av_sd_df], axis=1)
   
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
first_zerosnowdate_df = pd.DataFrame(np.array(first_zerosnowdate)).T
first_zerosnowdate_df.columns = hru_names_df
#first_zerosnowdate_df_obs = pd.DataFrame(np.array([[5985],[6200]]).T,columns=out_names)
first_zerosnowdate_df_obs = pd.DataFrame(np.array([5985]),columns=['2006'])

zerosnowdate_residual=[]
for hru in first_zerosnowdate_df.columns:
    zerosnowdate_residual.append((first_zerosnowdate_df[hru][0]-first_zerosnowdate_df_obs['2006'])/24)

zerosnowdate_residual_df = pd.DataFrame(np.reshape(np.array(zerosnowdate_residual),(np.size(out_names),hru_num)).T, columns=out_names)

#%%
#plt.xticks(x, hru[::3], rotation=25)
for namefile in out_names:
    x = list(np.arange(1,244))
    fig = plt.figure(figsize=(20,15))
    plt.bar(x,zerosnowdate_residual_df[namefile])
    plt.title(namefile, fontsize=42)
    plt.xlabel('hrus',fontsize=30)
    plt.ylabel('residual dosd (day)', fontsize=30)
    #vax.yaxis.set_label_coords(0.5, -0.1) 
    plt.savefig('SA2/'+namefile)

#%%out_names = ['AvInitial','Avas2l','Avas3m','Avtc2s','Avns3p','Avns4c']# 'Avwp2e','Avas2l','Avas3m','Avtc3t','Avtc4m','Avns3p','Avce2s','Avtc2s','Avtc3t','Avtc4m','Avns2a','Avns3p','Avns4c']
param_nam_list = ['LAIMIN','LAIMAX','winterSAI','summerLAI','rootingDepth','heightCanopyTop','heightCanopyBottom','throughfallScaleSnow','newSnowDenMin',
                  'albedoDecayRate', 'albedoMaxVisible', 'albedoMinVisible', 'albedoMaxNearIR', 'albedoMinNearIR',  
                  'z0Snow', 'albedoRefresh', 'mw_exp', 'fixedThermalCond_snow'] 

swe_obs2006 = pd.DataFrame (obs_swe['swe_mm']['2006-11-01':'2007-06-06'], columns=['swe_mm'])
swe_obs2007 = pd.DataFrame (obs_swe['swe_mm']['2007-12-03':'2008-06-08'], columns=['swe_mm'])
date_swe2006 = ['2006-11-01', '2006-11-30', '2007-01-01', '2007-01-30', '2007-03-05', '2007-03-12', '2007-03-19', '2007-03-26', '2007-04-02', '2007-04-18', '2007-04-23', '2007-05-02','2007-05-09', '2007-05-16', '2007-05-23', '2007-05-30', '2007-06-06']
swe_obs2006.set_index(pd.DatetimeIndex(date_swe2006),inplace=True)

#%%
DateSa2 = [i.strftime("%Y-%m-%d") for i in tvalueSa]
sax = np.arange(0,np.size(DateSa2))
sa_xticks = DateSa2
safig, saax = plt.subplots(1,1, figsize=(20,15))
plt.xticks(sax, sa_xticks[::1000], rotation=25, fontsize=20)
saax.xaxis.set_major_locator(ticker.AutoLocator())
plt.yticks(fontsize=20)
for hru in av_swe_df.columns:
    plt.plot(av_swe_df[hru])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]

plt.plot(swe_obs2006, 'ok', markersize=10)

plt.title('rainbow_SWE', position=(0.04, 0.88), ha='left', fontsize=40)
plt.xlabel('Time 2005-2006', fontsize=30)
plt.ylabel('SWE(mm)', fontsize=30)
#plt.show()
plt.savefig('SA2/swelj.png')
#%%finding max snowdepth and swe
#max_swe=[]
#for hrus in hru_names_df[0]:
#    max_swe.append(av_swe_df2[hrus].max())
#max_residual_SWE = max_swe - max_swe_obs
#
#swe_corsp_max2date = []
#for hrus in hru_names_df[0]:
#    swe_corsp_max2date.append(av_swe_df2[hrus][max_swe_date_obs])
#max_residual_swe_corsp = pd.DataFrame((swe_corsp_max2date - max_swe_obs), columns=['resCorspMaxSWE'])
#max_residual_swe_corsp.set_index(hru_names_df[0],inplace=True)
#residual_df = pd.concat([zerosnowdate_dif_obs,max_residual_swe_corsp], axis=1)
#residual_df_finale = residual_df.drop(['T0AvASsSWcTCj11112.0', 'T0AvASsSWcTCj11113.0','T2AvASsSWcTCj11112.0', 'T2AvASsSWcTCj11113.0','T4AvASsSWcTCj11112.0', 'T4AvASsSWcTCj11113.0',
#                                       'H2AvASsSWcTCj11112.0', 'H2AvASsSWcTCj11113.0','H4AvASsSWcTCj11112.0', 'H4AvASsSWcTCj11113.0'])
#    
#av_sd_df_finale = av_sd_df2.drop(['T0AvASsSWcTCj11112.0', 'T0AvASsSWcTCj11113.0','T2AvASsSWcTCj11112.0', 'T2AvASsSWcTCj11113.0','T4AvASsSWcTCj11112.0', 'T4AvASsSWcTCj11113.0',
#                                  'H2AvASsSWcTCj11112.0', 'H2AvASsSWcTCj11113.0','H4AvASsSWcTCj11112.0', 'H4AvASsSWcTCj11113.0'], axis=1)
#av_swe_df_finale = av_swe_df2.drop(['T0AvASsSWcTCj11112.0', 'T0AvASsSWcTCj11113.0','T2AvASsSWcTCj11112.0', 'T2AvASsSWcTCj11113.0','T4AvASsSWcTCj11112.0', 'T4AvASsSWcTCj11113.0',
#                                    'H2AvASsSWcTCj11112.0', 'H2AvASsSWcTCj11113.0','H4AvASsSWcTCj11112.0', 'H4AvASsSWcTCj11113.0'], axis=1)


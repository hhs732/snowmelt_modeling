###       /bin/bash runTestCases_docker.sh
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
#%% SWE data
date_swe = ['2006-11-01 08:00', '2006-11-30 08:00', '2007-01-01 08:00', '2007-01-30 08:00', '2007-03-05 08:00', '2007-03-12 08:00', 
            '2007-03-19 08:00', '2007-03-26 08:00', '2007-04-02 08:00', '2007-04-18 08:00', '2007-04-23 08:00', '2007-05-02 08:00', 
            '2007-05-09 08:00', '2007-05-16 08:00', '2007-05-23 08:00', '2007-05-30 08:00', '2007-06-06 08:00', 
            
            '2007-12-03 08:00', '2008-01-01 08:00', '2008-01-31 08:00', '2008-03-03 08:00', '2008-03-24 08:00', '2008-04-01 08:00', 
            '2008-04-14 08:00', '2008-04-22 08:00', '2008-04-28 08:00', '2008-05-06 08:00', '2008-05-12 08:00', '2008-05-19 08:00',
            '2008-05-26 08:00', '2008-06-02 08:00', '2008-06-08 08:00'] 
            
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
hruidxID = list(np.arange(101,105))
hru_num = np.size(hruidxID)
out_names = ['2006p','2007p']
paramModel = (np.size(out_names))*(hru_num)
hru_names =[]
for i in out_names:
    hru_names.append(['{}{}'.format(i, j) for j in hruidxID])
hru_names1 = np.reshape(hru_names,(paramModel,1))
hru_names_df = pd.DataFrame (hru_names1)
#%% reading output_swe files
av_ncfiles = ["Calib/param2_calibration2_swampAngel_2006-2007_senatorVariableDecayRate_1.nc", 
              "Calib/param2_calibration2_swampAngel_2007-2008_senatorVariableDecayRate_1.nc"] 
av_all = []
for ncfiles in av_ncfiles:
    av_all.append(Dataset(ncfiles))

for varname in av_all[0].variables.keys():
    var = av_all[0].variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)

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
av_swe_lst = []
numyear = 0
while numyear < (len(av_ncfiles)):
    av_swe_lst.append(pd.concat ([av_swe[numyear],av_swe[numyear+1]], axis=0, ignore_index=True))
    numyear = numyear+2

av_swe_df = []
for numfil in range (len(av_swe_lst)):
    for hrus in range (hru_num):
        av_swe_df.append(av_swe_lst[numfil][hrus])
av_swe_df = pd.DataFrame(np.array(av_swe_df).T, columns =  ['101p','102p','103p','104p'])

#%% output time step
TimeSa2006 = av_all[0].variables['time'][:] # get values
TimeSa2007 = av_all[1].variables['time'][:] # get values
TimeSa = np.concatenate((TimeSa2006, TimeSa2007), axis=0)
t_unitSa = av_all[1].variables['time'].units # get unit  "days since 1950-01-01T00:00:00Z"

try :

    t_cal = av_all[1].variables['time'].calendar

except AttributeError : # Attribute doesn't exist

    t_cal = u"gregorian" # or standard

tvalueSa = num2date(TimeSa, units=t_unitSa, calendar=t_cal)
DateSa = [i.strftime("%Y-%m-%d %H:%M") for i in tvalueSa] # -%d %H:%M to display dates as string #i.strftime("%Y-%m-%d %H:%M")  

av_swe_df.set_index(pd.DatetimeIndex(DateSa),inplace=True)
#%%
TimeSa2006 = av_all[0].variables['time'][:] # get values
tvalueSa2006 = num2date(TimeSa2006, units=t_unitSa, calendar=t_cal)
DateSa2006 = [i.strftime("%Y-%m") for i in tvalueSa2006]

TimeSa2007 = av_all[1].variables['time'][:] # get values
tvalueSa2007 = num2date(TimeSa2007, units=t_unitSa, calendar=t_cal)
DateSa2007 = [i.strftime("%Y-%m") for i in tvalueSa2007]
#%% day of snow disappearance-final output
av_sd_df.set_index(pd.DatetimeIndex(DateSa2007),inplace=True)
counter = pd.DataFrame(np.arange(0,np.size(av_sd_df['2006p101'])),columns=['counter'])
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
    
first_zerosnowdate_df = pd.DataFrame(np.reshape(first_zerosnowdate, ((np.size(out_names)),hru_num)).T, columns=out_names)
first_zerosnowdate_df_obs = pd.DataFrame(np.array([[5985],[6200]]).T,columns=out_names)
zerosnowdate_residual=[]
for year in out_names:
    for count in range (np.size(hruidxID)):
        zerosnowdate_residual.append((first_zerosnowdate_df[year][count]-first_zerosnowdate_df_obs[year])/24)

zerosnowdate_residual_df = pd.DataFrame(np.reshape(np.array(zerosnowdate_residual),(2,4)).T, columns=out_names)
#%%
#plt.xticks(x, hru[::3], rotation=25)
for modnames in out_names:
    x = list(np.arange(1,5))
    fig = plt.figure(figsize=(20,15))
    plt.bar(x,zerosnowdate_residual_df[modnames])
    plt.title(modnames, fontsize=42)
    plt.xlabel('parameters',fontsize=30)
    plt.ylabel('residual dosd', fontsize=30)
    #vax.yaxis.set_label_coords(0.5, -0.1) 
    plt.savefig(modnames)
#%%
sbx = np.arange(0,np.size(DateSa))

sb_xticks = DateSa
sbfig, sbax = plt.subplots(1,1)
plt.xticks(sbx, sb_xticks[::1000], rotation=25)
sbax.xaxis.set_major_locator(ticker.AutoLocator())
#%%out_names = ['AvInitial','Avas2l','Avas3m','Avtc2s','Avns3p','Avns4c']# 'Avwp2e','Avas2l','Avas3m','Avtc3t','Avtc4m','Avns3p','Avce2s','Avtc2s','Avtc3t','Avtc4m','Avns2a','Avns3p','Avns4c']
param_nam_list = ['winterSAI', 'summerLAI', 'rootingDepth', 'heightCanopyTop', 'heightCanopyBottom', 'throughfallScaleSnow', 
                  'albedoDecayRate', 'albedoMaxVisible', 'albedoMinVisible', 'albedoMaxNearIR', 'albedoMinNearIR', 'albedoRefresh',
                  'fixedThermalCond_snow', 'Fcapil', 'k_snow', 'mw_exp', 'z0Snow', 'critRichNumber', 'Louis79_bparam', 
                  'Louis79_cStar', 'Mahrt87_eScale', 'newSnowDenMin', 'newSnowDenMult', 'newSnowDenScal', 'constSnowDen', 
                  'newSnowDenAdd', 'newSnowDenMultTemp', 'newSnowDenMultWind', 'newSnowDenMultAnd'] 

swe_obs2006 = pd.DataFrame (obs_swe['swe_mm']['2006-11-01 08:00':'2007-06-06 08:00'])
swe_obs2007 = pd.DataFrame (obs_swe['swe_mm']['2007-12-03 08:00':'2008-06-08 08:00'])

#%%
fig = plt.figure(figsize=(20,15))
plt.plot(av_swe_df['101p'])
plt.plot(av_swe_df['102p'])
plt.plot(av_swe_df['103p'])
plt.plot(av_swe_df['104p'])
plt.plot(obs_swe, 'ok')#, linewidth=1)
plt.legend(['101pns','102pns','103pns','104pns','observed_swe']) 
plt.title('swe kg/m2',fontsize=50)#, position=(0.04, 0.88), ha='left', fontsize=12)
plt.xlabel('Time 2006-2008',fontsize=30)
plt.ylabel('swe (mm)',fontsize=30)
#plt.show()
plt.savefig('sweCalib4pns_param2.png')
#%%

#%% finding max snowdepth and swe
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


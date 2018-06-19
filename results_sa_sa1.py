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
            '2008-05-26 08:00', '2008-06-02 08:00', '2008-06-08 08:00', 
            
            '2008-12-02 08:00', '2009-01-01 08:00', '2009-02-01 08:00', '2009-02-28 08:00', '2009-03-09 08:00', '2009-03-16 08:00',
            '2009-03-24 08:00', '2009-03-30 08:00', '2009-04-07 08:00', '2009-04-15 08:00', '2009-04-22 08:00', '2009-04-29 08:00', 
            '2009-05-06 08:00', '2009-05-13 08:00', 
            
            '2009-11-27 08:00', '2009-12-31 08:00', '2010-01-31 08:00', '2010-03-02 08:00', '2010-03-21 08:00', '2010-04-05 08:00',
            '2010-04-12 08:00', '2010-04-20 08:00', '2010-04-26 08:00', '2010-05-03 08:00', '2010-05-11 08:00', '2010-05-17 08:00',
            '2010-05-24 08:00', 
            
            '2010-11-02 08:00', '2010-12-04 08:00', '2011-01-02 08:00', '2011-02-03 08:00', '2011-03-01 08:00', '2011-03-29 08:00',
            '2011-04-06 08:00', '2011-04-11 08:00', '2011-04-22 08:00', '2011-05-03 08:00', '2011-05-12 08:00', '2011-05-23 08:00', 
            '2011-06-01 08:00', '2011-06-08 08:00', '2011-06-14 08:00', 
            
            '2011-11-29 08:00', '2012-01-02 08:00', '2012-02-01 08:00', '2012-03-05 08:00', '2012-03-26 08:00', '2012-04-02 08:00', 
            '2012-04-08 08:00', '2012-04-16 08:00', '2012-04-23 08:00', '2012-04-30 08:00', '2012-05-07 08:00'] 
            
            #'2012-12-01 08:00', '2013-01-02 08:00', '2013-02-04 08:00', '2013-03-01 08:00', '2013-03-22 08:00', '2013-04-01 08:00',
            #'2013-04-10 08:00', '2013-04-16 08:00', '2013-04-22 08:00', '2013-04-30 08:00', '2013-05-06 08:00', '2013-05-13 08:00']

swe_mm = [58,  169, 267, 315, 499, 523, 503, 549, 611, 678, 654, 660, 711, 550, 443, 309, 84, 
          141, 300, 501, 737, 781, 837, 977, 950, 873, 894, 872, 851, 739, 538, 381, 
          133, 380, 456, 564, 512, 568, 626, 627, 715, 772, 764, 699, 698, 389, 
          89,  255, 347, 481, 608, 646, 682, 585, 553, 608, 520, 440, 302,  
          50,  165, 361, 454, 611, 704, 717, 774, 867, 951, 984, 999, 915, 699, 450, 
          130, 188, 290, 494, 542, 425, 433, 453, 413, 283, 150] 
          #55,  182, 305, 419, 481, 489, 508, 569, 624, 528, 405, 325]  

#obs_swe_date = pd.DataFrame (np.column_stack([date_swe,swe_mm]), columns=['date_swe','swe_mm'])
obs_swe = pd.DataFrame (swe_mm, columns=['swe_mm'])
obs_swe.set_index(pd.DatetimeIndex(date_swe),inplace=True)

max_swe_obs = max(obs_swe['swe_mm'])
max_swe_date_obs = obs_swe[obs_swe ['swe_mm']== max_swe_obs].index.tolist()    
#%% Snow depth data
with open("snowDepth_2007_2008.csv") as safd:
    reader = csv.reader(safd)
    raw_snowdepth = [r for r in reader]
sa_snowdepth_column = []
for csv_counter1 in range (len (raw_snowdepth)):
    for csv_counter2 in range (2):
        sa_snowdepth_column.append(raw_snowdepth[csv_counter1][csv_counter2])
sa_snowdepth=np.reshape(sa_snowdepth_column,(len (raw_snowdepth),2))
sa_snowdepth = sa_snowdepth[1:]
sa_sd_obs_date = pd.DatetimeIndex(sa_snowdepth[:,0])
sa_sd_obs = [float(value) for value in sa_snowdepth[:,1]]
snowdepth_obs_df = pd.DataFrame(sa_sd_obs, columns = ['observed snowdepth']) 
snowdepth_obs_df.set_index(sa_sd_obs_date,inplace=True)
#%%
hruidxID = list(np.arange(101,170))
hru_num = np.size(hruidxID)
out_names = ['AvInitial','Avas2l','Avas3m','Avtc2s', 'Avtc3t','Avtc4m','Avns3p','Avns4c']# 'Avwp2e','Avas2l','Avas3m','Avce2s','Avtc2s','Avtc3t','Avtc4m','Avns2a','Avns3p','Avns4c']
paramModel = (np.size(out_names))*(hru_num)
hru_names =[]
for i in out_names:
    hru_names.append(['{}{}'.format(i, j) for j in hruidxID])
hru_names1 = np.reshape(hru_names,(paramModel,1))
hru_names_df = pd.DataFrame (hru_names1)
#%% reading output_swe files
av_ncfiles = ["SA1/AvInitial_swampAngel_2007-2008_senatorVariableDecayRate_1.nc", 
              "SA1/Avas2l_swampAngel_2007-2008_senatorVariableDecayRate_1.nc", 
              "SA1/Avas3m_swampAngel_2007-2008_senatorVariableDecayRate_1.nc", 
              "SA1/Avtc2s_swampAngel_2007-2008_senatorVariableDecayRate_1.nc", 
              "SA1/Avtc3t_swampAngel_2007-2008_senatorVariableDecayRate_1.nc", 
              "SA1/Avtc4m_swampAngel_2007-2008_senatorVariableDecayRate_1.nc", 
              "SA1/Avns3p_swampAngel_2007-2008_senatorVariableDecayRate_1.nc", 
              "SA1/Avns4c_swampAngel_2007-2008_senatorVariableDecayRate_1.nc"]
              #"Avns2a_swampAngel_2007-2008_senatorVariableDecayRate_1.nc",
              #"Avce2s_swampAngel_2007-2008_senatorVariableDecayRate_1.nc",
              #"Avwp2e_swampAngel_2007-2008_senatorVariableDecayRate_1.nc", 
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

av_swe = []
for dfs in av_all:
    av_swe.append(pd.DataFrame(dfs['scalarSWE'][:]))
av_swe_df = pd.concat (av_swe, axis=1)
av_swe_df.columns = hru_names_df[0]

#%% output time step
TimeSa = av_all[0].variables['time'][:] # get values
t_unitSa = av_all[0].variables['time'].units # get unit  "days since 1950-01-01T00:00:00Z"

try :

    t_cal = av_all[0].variables['time'].calendar

except AttributeError : # Attribute doesn't exist

    t_cal = u"gregorian" # or standard

tvalueSa = num2date(TimeSa, units=t_unitSa, calendar=t_cal)
DateSa = [i.strftime("%Y-%m-%d %H:%M") for i in tvalueSa] # -%d %H:%M to display dates as string #i.strftime("%Y-%m-%d %H:%M")        
#%% day of snow disappearance-final output
av_sd_df.set_index(pd.DatetimeIndex(DateSa),inplace=True)
counter = pd.DataFrame(np.arange(0,np.size(av_sd_df['AvInitial101'])),columns=['counter'])
counter.set_index(av_sd_df.index,inplace=True)
av_sd_df2 = pd.concat([counter, av_sd_df], axis=1)

av_swe_df.set_index(pd.DatetimeIndex(DateSa),inplace=True)
counter = pd.DataFrame(np.arange(0,np.size(av_swe_df['AvInitial101'])),columns=['counter'])
counter.set_index(av_swe_df.index,inplace=True)
av_swe_df2 = pd.concat([counter, av_swe_df], axis=1)
#%%   
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
    
dayofsnowdisappearance_indicator_test = first_zerosnowdate[0]
zerosnowdate_residual_test = pd.DataFrame((first_zerosnowdate - dayofsnowdisappearance_indicator_test)/24, columns=['resSnowDisDate'])
zerosnowdate_residual_test.set_index(hru_names_df[0],inplace=True)

first_zerosnowdate_df = pd.DataFrame(np.reshape(first_zerosnowdate, ((np.size(out_names)),hru_num)).T, columns=out_names)
dayofsnowdisappearance_indicator = []
zerosnowdate_residual = []
out_names_df = pd.DataFrame(out_names)
for q in range (np.size(out_names)):
    dayofsnowdisappearance_indicator.append(first_zerosnowdate_df[out_names_df[0][q]][0])
    zerosnowdate_residual.append((first_zerosnowdate_df[out_names_df[0][q]] - dayofsnowdisappearance_indicator[q])/24)

zerosnowdate_residual_df = pd.DataFrame(np.array(zerosnowdate_residual).T, columns=out_names)
#%%
#plt.xticks(x, hru[::3], rotation=25)
for modnames in out_names:
    x = list(np.arange(1,70))
    fig = plt.figure(figsize=(20,15))
    plt.bar(x,zerosnowdate_residual_df[modnames])
    plt.title(modnames, fontsize=42)
    plt.xlabel('parameters')
    plt.ylabel('day of snow disappreance')
    #vax.yaxis.set_label_coords(0.5, -0.1) 
    plt.savefig(modnames)
#%%
sd_obs2008 = snowdepth_obs_df['observed snowdepth'][:]
DateSa2 = [i.strftime("%Y-%m") for i in tvalueSa]
sbx = np.arange(0,np.size(DateSa2))

sb_xticks = DateSa2
sbfig, sbax = plt.subplots(1,1)
plt.xticks(sbx, sb_xticks[::1000], rotation=25)
sbax.xaxis.set_major_locator(ticker.AutoLocator())
#%%out_names = ['AvInitial','Avas2l','Avas3m','Avtc2s','Avns3p','Avns4c']# 'Avwp2e','Avas2l','Avas3m','Avtc3t','Avtc4m','Avns3p','Avce2s','Avtc2s','Avtc3t','Avtc4m','Avns2a','Avns3p','Avns4c']
param_nam_list = ['winterSAI', 'summerLAI', 'rootingDepth', 'heightCanopyTop', 'heightCanopyBottom', 'throughfallScaleSnow', 
                  'albedoDecayRate', 'albedoMaxVisible', 'albedoMinVisible', 'albedoMaxNearIR', 'albedoMinNearIR', 'albedoRefresh',
                  'fixedThermalCond_snow', 'Fcapil', 'k_snow', 'mw_exp', 'z0Snow', 'critRichNumber', 'Louis79_bparam', 
                  'Louis79_cStar', 'Mahrt87_eScale', 'newSnowDenMin', 'newSnowDenMult', 'newSnowDenScal', 'constSnowDen', 
                  'newSnowDenAdd', 'newSnowDenMultTemp', 'newSnowDenMultWind', 'newSnowDenMultAnd'] 

swe_obs2008 = obs_swe['swe_mm']['2007-12-03 08:00':'2008-06-08 08:00']

date_swem = ['2007-12-03 08:00', '2008-01-01 08:00', '2008-01-31 08:00', '2008-03-03 08:00', '2008-03-24 08:00', '2008-04-01 08:00', 
            '2008-04-14 08:00', '2008-04-22 08:00', '2008-04-28 08:00', '2008-05-06 08:00', '2008-05-12 08:00', '2008-05-19 08:00',
            '2008-05-26 08:00', '2008-06-02 08:00', '2008-06-08 08:00']
swe_m = [141*0.001, 300*0.001, 501*0.001, 737*0.001, 781*0.001, 837*0.001, 977*0.001, 950*0.001, 873*0.001, 894*0.001, 872*0.001, 851*0.001, 739*0.001, 538*0.001, 381*0.001] 
obs_swem = pd.DataFrame (swe_m, columns=['swe_m'])
obs_swem.set_index(pd.DatetimeIndex(date_swem),inplace=True)

fig = plt.figure(figsize=(20,15))
plt.plot(av_swe_df2['Avtc2s101'])
plt.plot(av_swe_df2['Avtc2s102'])
plt.plot(av_swe_df2['Avtc2s103'])
plt.plot(swe_obs2008, 'ok')#, linewidth=1)
plt.legend(['swe1','swe2','swe3','observed_swe']) 
plt.title('swe kg/m2')#, position=(0.04, 0.88), ha='left', fontsize=12)
plt.xlabel('Time 2010-2011')
plt.ylabel('sd (m)')
#plt.show()
plt.savefig('swe.png')
#%%
fig = plt.figure(figsize=(20,15))

plt.plot(av_sd_df2['Avtc2s101'])
plt.plot(av_sd_df2['Avtc2s102'])
plt.plot(av_sd_df2['Avtc2s103'])
plt.plot(snowdepth_obs_df ['observed snowdepth'][:], 'k')
#plt.plot(obs_swem, 'ok')#, linewidth=1)

plt.legend(['sd1','sd2','sd3','observed_sd'])  #'swe1','swe2','swe3','observed_swe', 
plt.title('snowdepth (m)')#, position=(0.04, 0.88), ha='left', fontsize=12)
plt.xlabel('Time 2010-2011')
plt.ylabel('sd (m)')
#plt.show()
plt.savefig('Sd.png')
#%%
#snow_density = []
#for sd in range (len (av_swe_df2['AvInitial101'])):
#    snow_density.append(av_swe_df2['AvInitial101'][sd]/av_sd_df2['AvInitial101'][sd])
#
#fig = plt.figure(figsize=(20,15))
#plt.plot(snow_density) 
#plt.savefig('Density.png')                      


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


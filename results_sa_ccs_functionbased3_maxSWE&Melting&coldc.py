###       /bin/bash runTestCases_docker.sh
import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
from datetime import datetime
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats
from sklearn.metrics import mean_squared_error
import itertools
import csv
from allNcFilesCCS import av_ncfiles

#%% defining functions  13132sj2230

def mySubtract(myList,num):
    return list(np.subtract(myList,num))
def myMultiply(myList,num):
    return list(np.multiply(myList,num))
def sum2lists (list1,list2):
    return list(np.add(list1,list2))

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
    counter = pd.DataFrame(np.arange(0,np.size(variableNameDF[hruname[0]])),columns=['counter'])
    counter.set_index(variableNameDF.index,inplace=True)
    variableNameDF = pd.concat([counter, variableNameDF], axis=1)
    return variableNameDF

def readSomePartofVariableDatafromNcfilesDatasetasDF(NcfilesDataset,variable,hruname,paretoNameDF):
    variableNameList = []
    for datasets in NcfilesDataset:
        variableNameList.append(pd.DataFrame(datasets[variable][:]))
    variableDF = pd.concat (variableNameList, axis=1)
    variableDF.columns = hruname
    desiredDataframe = []
    for paretos in paretoNameDF:
        desiredDataframe.append(variableDF[paretos])
    return desiredDataframe

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

def snowLayerAttributeforSpecificDate(layerattributefile,hruname,sumlayer,snowlayer): #like snowlayertemp, volFracosIce, ....
    snowlayerattribute = []
    for names in hruname:
        snowlayerattribute.append(list(layerattributefile[names][sumlayer[names][0]:sumlayer[names][0]+snowlayer[names][0]]))
    return snowlayerattribute

def depthOfLayers(heightfile):
    finalHeight = []
    for lsts in heightfile:
        sumSofar=0
        lstscopy= lsts[:]
        lstscopy.reverse()
        height_ls = []
        for values in lstscopy:
            height=2*(abs(values)-sumSofar)
            height_ls.append(height)
            sumSofar+=height
        #print ("original:", height_ls) 
        height_ls.reverse()
        #print ("after reverse:", height_ls)
        finalHeight.append(height_ls)
    return finalHeight

def coldContentFunc(hruname,volfracLiq,volfracIce,snowlayertemp,layerHeight):
    densityofWater = 997. #kg/m³
    densityofIce = 917. #kg/m3
    heatCapacityIce2 = -2102./1000000. #Mj kg-1 m3-1
    coldcontent = []
    for nlst in range (np.size(hruname)):
        swe = np.array(sum2lists(myMultiply(volfracLiq[nlst],densityofWater/1000.),myMultiply(volfracIce[nlst],densityofIce/1000.))) #[swe] = m
        temp = np.array(mySubtract(snowlayertemp[nlst],273.2))
        HCItHS = np.array(myMultiply(heatCapacityIce2,layerHeight[nlst]))
        cct = sum(list(swe*temp*HCItHS))
        coldcontent.append(cct)
    return coldcontent

def SWEandSWEDateforSpecificDate(hruname,hour,swe_df,dosd):
    SWE = []
    SWEdate = []
    for names in hruname:
        if swe_df[names][hour]>0:
            SWE.append(swe_df[names][hour])
            SWEdate.append(hour)
        else: 
            SWE.append(float(swe_df[names][dosd[names]-1]))
            SWEdate.append(float(dosd[names]-1))
    return SWE,SWEdate

def meltingRateBetween2days(swe1,swe2,sweDate1,sweDate2):
    mdeltaday = []; mdeltaSWE = []; meltingrate = []; #cm/day
    for counterhd in range (np.size(swe1)):
        mdeltaday.append(float(sweDate2[counterhd]-sweDate1[counterhd]))
        mdeltaSWE.append(float(swe1[counterhd]-swe2[counterhd]))
        if mdeltaday[counterhd]==0:
            meltingrate.append(float(0))
        else: meltingrate.append(float(0.1*24*mdeltaSWE[counterhd]/mdeltaday[counterhd]))
    return meltingrate
    
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

swe_obs2006 = pd.DataFrame (obs_swe['swe_mm']['2006-11-01':'2007-06-06'], columns=['swe_mm'])
swe_obs2007 = pd.DataFrame (obs_swe['swe_mm']['2007-12-03':'2008-06-08'], columns=['swe_mm'])
date_swe2006 = ['2006-11-01 11:10','2006-11-30 12:30','2007-01-01 11:10','2007-01-30 10:35','2007-03-05 14:30', 
                '2007-03-12 14:00','2007-03-19 12:30','2007-03-26 12:30','2007-04-02 12:30','2007-04-18 08:35', 
                '2007-04-23 10:30','2007-05-02 08:40','2007-05-09 08:50','2007-05-16 09:00','2007-05-23 08:30',
                '2007-05-30 09:00','2007-06-06 08:15']
swe_obs2006.set_index(pd.DatetimeIndex(date_swe2006),inplace=True)  
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

    
years = ['2006','2007']
out_names = ['ccs','hs']
 
bestHRUs =['13132sj2230','21121sc1233','22111sc1233','21121sc2233','22111sc2233','22121sc2331','13321lj1230','13231lj2130',
            '12231lc1123','12331lc1131','23111lc1232','13231lc2122','12231lc2123','13222lc2231','22121lc2231','21121lc2233',
            '13222mj1130','12332mj1330','22221mj1330','21231mj2130','23131mj2210','22231mj2320','22122mj2330','33211mj2330',
            '12232mc1122','31131mc1131','31121mc1132','32111mc1132','32121mc1221','13312mc1233','22231mc1321','31131mc1323',
            '13223mc1331','31122mc1332','32112mc1332','12232mc2113','13123mc2122','21231mc2122','12123mc2123','12332mc2222',
            '13322mc2222','23112mc2223','21232mc2231','33112mc2231','11233mc2232','31221mc2232','32211mc2232','12213mc2233',
            '23221mc2312','22122mc2322']

hru_num = np.size(bestHRUs)
paramModel = (np.size(out_names))*(hru_num)

hru_names =[]
for i in out_names:
    hru_names.append(['{}{}'.format(i, j) for j in bestHRUs])

hru_names1 = np.reshape(hru_names,(paramModel,1))
hru_names_df = pd.DataFrame (hru_names1)
#%%  reading output files
av_all = readAllNcfilesAsDataset(av_ncfiles)
DateSa = date(av_all,"%Y-%m-%d %H:%M")
#test0=av_all[0]['scalarSWE'][:,0]

#test1=av_all[0]['scalarSWE'][:,1]
av_sd_df = readVariablefromNcfilesDatasetasDF(av_all,'scalarSnowDepth',hru_names_df[0])
av_sd_df.set_index(pd.DatetimeIndex(DateSa),inplace=True)
av_swe_df =  readVariablefromNcfilesDatasetasDF(av_all,'scalarSWE',hru_names_df[0])
#test00 = av_swe_df['11111lj1110']
#test11 = av_swe_df['11112lj1110']

av_swe_df.set_index(pd.DatetimeIndex(DateSa),inplace=True)
#Testsj = av_swe_df['22122mc2322']
#%% ploting annual swe curves
DateSa2 = date(av_all,"%Y-%m-%d")
sax = np.arange(0,np.size(DateSa2))
sa_xticks = DateSa2
safig, saax = plt.subplots(1,1, figsize=(20,15))
plt.xticks(sax, sa_xticks[::1000], rotation=25, fontsize=20)
saax.xaxis.set_major_locator(ticker.AutoLocator())
plt.yticks(fontsize=20)
for hru in hru_names_df[0]:
    plt.plot(av_swe_df[hru])#, sbx, swe_obs2006, 'k--', linewidth=0.5)#, label='wwe', color='maroon') param_nam_list[q] color_list[q]

plt.plot(swe_obs2006, 'ok', markersize=10)

plt.title('rainbow_SWE', position=(0.04, 0.88), ha='left', fontsize=40)
plt.xlabel('Time 2006-2007', fontsize=40)
plt.ylabel('SWE(mm)', fontsize=40)
#plt.legend()
#plt.show()
plt.savefig('ccs/'+'sweCCsandHS.png')
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
dosd_obs = pd.DataFrame(np.array([5943]),columns=['2006'])

dosd_residual=[]
for hru in dosd_df.columns:
    dosd_residual.append((dosd_df[hru][0]-dosd_obs['2006'])/24)

dosd_residual_df = pd.DataFrame(np.reshape(np.array(dosd_residual),(np.size(out_names),hru_num)).T, columns=out_names)

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
#%%**************************************************************************************************
# *********************** finding max corespondance swe for '2007-05-09 08:50'***********************
#'2007-04-18' 4776: 4800, '2007-04-23' 4896:4920, '2007-05-02' 5112:5136
#Group1: '2007-03-12 14:00' (3902),'2007-03-19 12:30 (4068)','2007-03-26 12:30 (4236)','2007-04-02 12:30'(4404),
#Group2: '2007-04-18 08:35' (4784),'2007-04-23 10:30 (4907)','2007-05-02 08:40'(5121), 

maxSWE = readSpecificDatafromAllHRUs(av_swe_df,hru_names_df[0],5289)
maxSWE_obs = [711]  

av_swe_df2 = av_swe_df.copy(); av_swe_df2.set_index(av_swe_df['counter'],inplace=True)
realMaxSWE = av_swe_df2.max()
realMaxSWE_date = av_swe_df2.idxmax()
#%%**************************************************************************************************
# ********************** calculating snowmelt rate based on SWE *************************************
sweM1,SWE1date = SWEandSWEDateforSpecificDate(hru_names_df[0],5289,av_swe_df,dosd_df)
sweM2,SWE2date = SWEandSWEDateforSpecificDate(hru_names_df[0],5457,av_swe_df,dosd_df)
sweM3,SWE3date = SWEandSWEDateforSpecificDate(hru_names_df[0],5793,av_swe_df,dosd_df)
sweM4,SWE4date = SWEandSWEDateforSpecificDate(hru_names_df[0],5960,av_swe_df,dosd_df)
#%%
meltingrate1 = meltingRateBetween2days(sweM1,sweM2,SWE1date,SWE2date)
meltingrate2 = meltingRateBetween2days(sweM2,sweM3,SWE2date,SWE3date)
meltingrate3 = meltingRateBetween2days(sweM3,sweM4,SWE3date,SWE4date)

meltingrateAvg_mod = []
for countermr in range (np.size(meltingrate1)):
    meltingrateAvg_mod.append((meltingrate1[countermr]+meltingrate2[countermr]+meltingrate3[countermr])/3)
#%%
sweMR = [711, 550, 309, 84]
mrDate = ['2007-05-09 08:50 5289','2007-05-16 09:00 5457','2007-05-30 09:00 5793','2007-06-06 08:15 5960']  
meltingrate1_obs = np.array([0.1*24*(711-550.)/(5457.-5289)])
meltingrate2_obs = np.array([0.1*24*(550.-309)/(5793.-5457)])
meltingrate3_obs = np.array([0.1*24*(309-84.)/(5960-5793.)])
meltingrateAvg_obs = (meltingrate1_obs+meltingrate2_obs+meltingrate3_obs)/3.
#'2007-05-09 08:50':5289, to '2007-06-06 08:15': 5960, 
#swe_mm = [711, 84]
meltingRate_obs = [0.1*24*(711-84.)/(5960-5289.)] #cm/day
#%% new criteria-swe before max 
swe2bfrmax = readSpecificDatafromAllHRUs(av_swe_df,hru_names_df[0],4906)
swe3bfrmax = readSpecificDatafromAllHRUs(av_swe_df,hru_names_df[0],4785)

swe2bmax_obs = [654]
swe3bmax_obs = [678]
#%% defining criteria
meltingRateCrit = [abs(values) for values in mySubtract(meltingrateAvg_mod,meltingrateAvg_obs)]
maxSWEcrit = [abs(values) for values in mySubtract(maxSWE,maxSWE_obs)]
swe2bmaxCrit = [abs(values) for values in mySubtract(swe2bfrmax,swe2bmax_obs)]
swe3bmaxCrit = [abs(values) for values in mySubtract(swe3bfrmax,swe3bmax_obs)]
#fig = plt.figure(figsize=(20,15))
#xs = meltingRateCrit
#ys = maxSWEcrit
#plt.scatter(xs, ys)
#plt.title('criteria for best combos')
#plt.xlabel('delta_maxSWE (mm)',fontsize=40)
#plt.ylabel('delta_meltingRate (cm/day)',fontsize=40)
#plt.savefig('SA2/'+'maxswe_meltinRateAvg')
#%%
meltingRateCrit_df = pd.DataFrame(meltingRateCrit, columns=['meltingRate'])
maxSWECrit_df = pd.DataFrame(maxSWEcrit, columns=['maxSWE'])
swe2bmaxCrit_df = pd.DataFrame(swe2bmaxCrit, columns=['swe2bmaxCrit'])
swe3bmaxCrit_df = pd.DataFrame(swe3bmaxCrit, columns=['swe3bmaxCrit'])

criteria_df = pd.concat([meltingRateCrit_df, maxSWECrit_df, swe2bmaxCrit_df, swe3bmaxCrit_df], axis=1) #coldcontentcrit_df, 
criteria_df.set_index(hru_names_df[0],inplace=True)
Apareto_model_param = pd.DataFrame(criteria_df.index[((criteria_df['maxSWE']) <= 5) & ((criteria_df['meltingRate'])<=0.03)].tolist()) # & ((criteria_df['coldContent'])<=7)

#Apareto_model_param = pd.DataFrame(criteria_df.index[((criteria_df['maxSWE']) <= 45) & ((criteria_df['swe2bmaxCrit']) <= 67) & ((criteria_df['swe3bmaxCrit']) <= 67) & ((criteria_df['meltingRate'])<=0.1)].tolist()) # & ((criteria_df['coldContent'])<=7)

#for hru in Apareto_model_param[0]:
#    plt.plot(av_swe_df[hru])
#plt.plot(swe_obs2006, 'ok', markersize=5)
##plt.legend()
#plt.savefig('best_combo')
#%% **************************************************************************************************
## ************************** calculating cold content ************************************************
#observed cold content in each day
heatCapacityIce1 = -2102. #J kg-1 K-1
swe0305 = [0.83,0.94,0.65,0.91,0.68,0.61,0.37];  T0305 = [-0.45,-1.45,-2.45,-3.4,-5.35,-8,-4.7] #14:30 #3734
swe0312 = [0.81,0.91,0.85,0.76,0.69,0.61,0.38,0.22]; T0312 = [-0.8,-1.8,-2.4,-3.2,-3.8,-4.3,-4.6,-1] #14:00 #3902
swe0319 = [0.83,0.84,0.95,0.75,0.69,0.36,0.27,0.34]; T0319 = [-0.4,-1,-1.4,-1.5,-1.5,-1.5,-0.7,0] # 12:30 #4068
swe0326 = [0.93,0.87,0.75,0.81,0.70,0.80,0.63]; T0326 = [-0.06,-0.2,-0.2,-0.2,-0.3,-0.7,-0.3] #12;30 #4236

cc0305 = [sum(np.multiply((heatCapacityIce1/1000000.),np.multiply(swe0305,T0305)))]
cc0312 = [sum(np.multiply((heatCapacityIce1/1000000.),np.multiply(swe0312,T0312)))]
cc0319 = [sum(np.multiply((heatCapacityIce1/1000000.),np.multiply(swe0319,T0319)))]
cc0326 = [sum(np.multiply((heatCapacityIce1/1000000),np.multiply(swe0326,T0326)))] #after this date, snowpack got isothermal

coldContenAvg_obs = np.mean([cc0305,cc0312,cc0319])
#%% calculating modeled cold content in each day for 2 sets of best combos
nvolfracIce0 = pd.concat(readSomePartofVariableDatafromNcfilesDatasetasDF(av_all,'mLayerVolFracIce',hru_names_df[0],hru_names_df[0]), axis=1)
nvolfracliq0 = pd.concat(readSomePartofVariableDatafromNcfilesDatasetasDF(av_all,'mLayerVolFracLiq',hru_names_df[0],hru_names_df[0]), axis=1)
nheight0 = pd.concat(readSomePartofVariableDatafromNcfilesDatasetasDF(av_all,'mLayerHeight',hru_names_df[0],hru_names_df[0]), axis=1)
nlayerTemp0 =  pd.concat(readSomePartofVariableDatafromNcfilesDatasetasDF(av_all,'mLayerTemp',hru_names_df[0],hru_names_df[0]), axis=1)
nsnow0 =  pd.concat(readSomePartofVariableDatafromNcfilesDatasetasDF(av_all,'nSnow',hru_names_df[0],hru_names_df[0]), axis=1)
nlayer0 = pd.concat(readSomePartofVariableDatafromNcfilesDatasetasDF(av_all,'nLayers',hru_names_df[0],hru_names_df[0]), axis=1)

nvolfracWat = av_all[0].variables['mLayerVolFracWat'][:]
#%%
#nvolfracIce1 = pd.concat(readSomePartofVariableDatafromNcfilesDatasetasDF(av_all,'mLayerVolFracIce',hru_names_df[0],Apareto_model_param1[0]), axis=1)
#nvolfracliq1 = pd.concat(readSomePartofVariableDatafromNcfilesDatasetasDF(av_all,'mLayerVolFracLiq',hru_names_df[0],Apareto_model_param1[0]), axis=1)
#nheight1 = pd.concat(readSomePartofVariableDatafromNcfilesDatasetasDF(av_all,'mLayerHeight',hru_names_df[0],Apareto_model_param1[0]), axis=1)
#nlayerTemp1 =  pd.concat(readSomePartofVariableDatafromNcfilesDatasetasDF(av_all,'mLayerTemp',hru_names_df[0],Apareto_model_param1[0]), axis=1)
#nsnow1 =  pd.concat(readSomePartofVariableDatafromNcfilesDatasetasDF(av_all,'nSnow',hru_names_df[0],Apareto_model_param1[0]), axis=1)
#nlayer1 = pd.concat(readSomePartofVariableDatafromNcfilesDatasetasDF(av_all,'nLayers',hru_names_df[0],Apareto_model_param1[0]), axis=1)

#%% number of snowlayer
nvolfracliq0 = pd.concat(readSomePartofVariableDatafromNcfilesDatasetasDF(av_all,'mLayerVolFracWat',hru_names_df[0],hru_names_df[0]), axis=1)

nsnow03260 = pd.DataFrame(readSpecificDatafromAllHRUs(nsnow0,hru_names_df[0],4236)).T; nsnow03260.columns = hru_names_df[0]
#nsnow03261 = pd.DataFrame(readSpecificDatafromAllHRUs(nsnow1,Apareto_model_param1[0],4236)).T; nsnow03261.columns = Apareto_model_param1[0]

# sum of all layers befor target layer
sumlayer03260 = pd.DataFrame(sumBeforeSpecificDatafromAllHRUs(nlayer0,hru_names_df[0],4236)).T; sumlayer03260.columns = hru_names_df[0]
#sumlayer03261 = pd.DataFrame(sumBeforeSpecificDatafromAllHRUs(nlayer1,Apareto_model_param1[0],4236)).T; sumlayer03261.columns = Apareto_model_param1[0]

#%%snow layer temperature
snowlayertemp0326 = snowLayerAttributeforSpecificDate(nlayerTemp0,hru_names_df[0],sumlayer03260,nsnow03260)

#%% volumetric fraction of ice in snow layers
volfracIce0326 = snowLayerAttributeforSpecificDate(nvolfracIce0,hru_names_df[0],sumlayer03260,nsnow03260)

#%% volumetric fraction of liquid in snow layers
volfracLiq0326 = snowLayerAttributeforSpecificDate(nvolfracliq0,hru_names_df[0],sumlayer03260,nsnow03260)

#%% height of each snow layer
height0326 = snowLayerAttributeforSpecificDate(nheight0,hru_names_df[0],sumlayer03260,nsnow03260)

height0326layer = depthOfLayers(height0326)

#%% cold content in each day
coldcontent0326 = coldContentFunc(hru_names_df[0],volfracLiq0326,volfracIce0326,snowlayertemp0326,height0326layer)

#%%
coldcontent0326_df = pd.DataFrame(np.reshape(coldcontent0326,[2,50]).T,columns=['cc_ccs','cc_hhs']) #cc0326
meltingrateAvg_mod_df = pd.DataFrame(np.reshape(meltingrateAvg_mod,[2,50]).T,columns=['mr_ccs','mr_hhs']) #meltingrateAvg_obs
maxSWE_df = pd.DataFrame(np.reshape(maxSWE,[2,50]).T,columns=['mxs_ccs','mxs_hhs']) #maxSWE_obs


#d1 = [coldcontent0326_df['cc_hhs'],coldcontent0326_df['cc_ccs']]
#fig, axs = plt.subplots(1, 2)
#bp1 = axs[0].boxplot(d1)
#axs[0].set_title('cold content (MJ kg-1 K-1)')
#bp1['boxes'][0].set(color='seagreen', linewidth=2, facecolor = 'skyblue', hatch = '/')
#bp1['boxes'][1].set(color='orange', linewidth=2, facecolor = 'olive', hatch = '/')

#d2 = [meltingrateAvg_mod_df['mr_hhs'],meltingrateAvg_mod_df['mr_ccs']]
#bp2 = axs[1].boxplot(d2)
#axs[1].set_title('melting rate (cm/day)')
#bp2['boxes'][0].set(color='seagreen', linewidth=2, facecolor = 'skyblue', hatch = '/')
#bp2['boxes'][1].set(color='orange', linewidth=2, facecolor = 'olive', hatch = '/')
#boxprops = dict(linestyle='--', linewidth=3, color='darkgoldenrod')
#flierprops = dict(marker='o', markerfacecolor='green', markersize=12, linestyle='none')
#medianprops = dict(linestyle='-.', linewidth=2.5, color='firebrick')
#meanpointprops = dict(marker='D', markeredgecolor='black', markerfacecolor='firebrick')
#meanlineprops = dict(linestyle='--', linewidth=2.5, color='purple')

boxprops = dict(linewidth=15, color='darkgreen')
flierprops = dict(markerfacecolor='darkgreen', markersize=12, linestyle='none')
medianprops = dict(linewidth=15, color='darkgreen')
meanpointprops = dict( markeredgecolor='darkgreen', markerfacecolor='darkgreen')
meanlineprops = dict(linewidth=5, color='darkgreen')


d1 = [maxSWE_df['mxs_ccs'],maxSWE_df['mxs_hhs']]
d2 = [coldcontent0326_df['cc_ccs'],coldcontent0326_df['cc_hhs']]
d3 = [meltingrateAvg_mod_df['mr_ccs'],meltingrateAvg_mod_df['mr_hhs']]

fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20,15), sharey=True)
axes[0].boxplot(d1, widths=(0.6, 0.6), boxprops=boxprops)
axes[0].set_ylabel('max SWE (mm)', fontsize=32)

axes[1].boxplot(d2, widths=(0.6, 0.6), boxprops=boxprops)
axes[1].set_ylabel('cold content (mj m-2 s-1)', fontsize=32)

axes[2].boxplot(d3, widths=(0.6, 0.6), boxprops=boxprops)
axes[2].set_ylabel('melting rate (cm/day)', fontsize=32)

for ax in axes.flatten():
    ax.set_xticklabels(['hst', 'ccs'], fontsize=32)
#axes.set_yticks(fontsize=22)

#bp3 = plt.boxplot(, patch_artist=True)
#bp3['boxes'][0].set(color='seagreen', linewidth=5, facecolor = 'skyblue', hatch = '/')
#bp3['boxes'][1].set(color='orange', linewidth=5, facecolor = 'olive', hatch = '/')
#plt.title('max swe (mm)', fontsize=32)


plt.savefig('ccs/'+'boxplot.png')
#plt.close()
#%% #3d Plot
#fig = plt.figure(figsize=(20,15))
#ax = fig.add_subplot(111, projection='3d')
#
#xs = coldcontentcrit
#ys = meltingRateCrit
#zs = maxSWEcrit
#ax.scatter(xs, ys, zs)
#
#ax.set_xlabel('cold content (Mj/m2)',fontsize=20)
#ax.set_ylabel('melting rate (cm/day)',fontsize=20)
#ax.set_zlabel('maxSWE (mm)',fontsize=20)
#plt.savefig('ccs/'+'+cc')
##plt.show()















    
    
    
    
    
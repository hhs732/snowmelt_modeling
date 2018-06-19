#%matplotlib inline    /bin/bash runTestCases_docker.sh
import numpy as np
import matplotlib.pyplot as plt 
from netCDF4 import Dataset,netcdftime,num2date
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import csv
#%% hru names
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

p15 = [0.001] #[0.001, 0.002] #z0Snow

p16 = [1]# 1, 3, 6] #albedoRefresh |       1.0000 |       1.0000 |      10.0000

p17 = [2] #2, 3, 4] #mw_exp exponent for meltwater flow

p18 = [0.4] #0.2, 0.4 , 0.6] #fixedThermalCond_snow

#p21 = [0.700, 1.000, 1.500] #Mahrt87_eScale  
#p14 = [0.040, 0.060, 0.080] #Fcapil
#p15 = [0.100, 0.015, 0.350] #k_snow
#p18 = [0.150, 0.200, 0.400] #critRichNumber  
#p19 = [9.300, 9.400, 9.500] #Louis79_bparam  
#p20 = [5.200, 5.300, 5.400] #Louis79_cStar 
#p15 = [105] #51, 70, 105] #constSnowDen 70.00, 100.0, 170.0    55 56 57 

def hru_ix_ID(p10, p11, p12, p13, p14):#, p15, p16, p17, p18):
    ix10 = np.arange(1,len(p10)+1)
    ix11 = np.arange(1,len(p11)+1)
    ix12 = np.arange(1,len(p12)+1)
    ix13 = np.arange(1,len(p13)+1)
    ix14 = np.arange(1,len(p14)+1)
#    ix15 = np.arange(1,len(p15)+1)
#    ix16 = np.arange(1,len(p16)+1)
#    ix17 = np.arange(1,len(p17)+1)
#    ix18 = np.arange(1,len(p18)+1)
    c = list(itertools.product(ix10,ix11,ix12,ix13,ix14))#,ix15,ix16,ix17,ix18))
    ix_numlist=[]
    for tup in c:
        ix_numlist.append(''.join(map(str, tup)))
    new_list = [float(i) for i in ix_numlist]

    return(new_list)  

hruidxID = hru_ix_ID(p10, p11, p12, p13, p14)#, p15, p16, p17, p18)
#
hru_num = np.size(hruidxID)
#%% #Swamp Angel forcing data
#swampangel_forcing = open('swamp_angel_forcingdata2_corrected.csv', 'rb')
#sa_forcing = csv.reader(swampangel_forcing)#, delimiter=',')
with open("swamp_angel_forcingdata2_corrected_precipCalib_final.csv") as safd:
    reader = csv.reader(safd)
    data_forcing = [r for r in reader]
data_forcing2 = data_forcing[1:]
sa_fd_column = []
for csv_counter1 in range (len (data_forcing2)):
    for csv_counter2 in range (26):
        sa_fd_column.append(float(data_forcing2[csv_counter1][csv_counter2]))
sa_forcing=np.reshape(sa_fd_column,(len (data_forcing2),26))
#%% #Senator Beck basin forcing data-Air pressure
#fr = Dataset('testCases_data/inputData/fieldData/reynolds/forcing_above_aspen.nc')
sbFD = Dataset('SenatorBeck_forcing.nc')
#ap_sb = sbFD.variables['airpres'][:]
TimeSb = sbFD.variables['time'][:] # get values
t_unitS = sbFD.variables['time'].units # get unit  "days since 1950-01-01T00:00:00Z"

try :

    t_cal = sbFD.variables['time'].calendar

except AttributeError : # Attribute doesn't exist

    t_cal = u"gregorian" # or standard

tvalueSb = num2date(TimeSb, units=t_unitS, calendar=t_cal)
sbDate = [i.strftime("%Y-%m-%d %H:%M") for i in tvalueSb]

#%% # ***** How to get time value into Python DateTIme Objects *****
# time step for swamp Angel
first_day = 5752.0000
#last_day = 7943.04166667
interval = 0.04166667
time_lentgh = 61369
saTime = []
for counter1 in range (time_lentgh):
    saTime.append(first_day+interval*counter1)
#saTime = np.concatenate([TimeSb,extended_year])

tvalueSa = num2date(saTime, units=t_unitS, calendar=t_cal)
saDate = [i.strftime("%Y-%m-%d %H:%M") for i in tvalueSa] # to display dates as string #i.strftime("%Y-%m-%d %H:%M")
#%% swamp angel forcing data dataframe
sa_forcing2 = sa_forcing[:,[3,6,7,8,9,10,11,12,13,14]]   #[:,1:]  #np.delete(sa_forcing, 0, axis=1)
sa_df = pd.DataFrame (sa_forcing2, columns=['day','pptrate','SWRadAtm','LWRadAtm','airtemp','windspd','airpres',
                                            'spechum','relative_humidity[%]','dewpoint_temperature[K]'])
sa_df.set_index(pd.DatetimeIndex(saDate),inplace=True)
#%% swamp angel Temp and ppt average to select the year
temp_data=pd.Series(np.array(sa_df['airtemp']),index=pd.DatetimeIndex(saDate))
temp_meanyr=temp_data.resample("A").mean()

ppt_data=pd.Series(np.array(sa_df['pptrate']),index=pd.DatetimeIndex(saDate))
ppt_meanyr=ppt_data.resample("A").sum()
#%% # I go through the a forcing file from the test cases to see what our netcdfs need to look like
print sbFD.file_format
# read out variables, data types, and dimensions of original forcing netcdf
for varname in sbFD.variables.keys():
    var = sbFD.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)

for dimname in sbFD.dimensions.keys():
    dim = sbFD.dimensions[dimname]
    #print(dimname, len(dim), dim.isunlimited())

atemp = sbFD.variables['airtemp'][:]
#%% make new nc file
new_fc_sa = Dataset("SwampAngel_forcing_sa2.nc",'w',format='NETCDF3_CLASSIC')
# define dimensions 
hru = new_fc_sa.createDimension('hru', hru_num)
time = new_fc_sa.createDimension('time', None)
# define variables
hruid = new_fc_sa.createVariable('hruId', np.int32,('hru',))
lat = new_fc_sa.createVariable('latitude', np.float64,('hru',))
lon = new_fc_sa.createVariable('longitude', np.float64,('hru',))
ds = new_fc_sa.createVariable('data_step', np.float64)
times = new_fc_sa.createVariable('time', np.float64,('time',))
lwrad = new_fc_sa.createVariable('LWRadAtm', np.float64,('time','hru'), fill_value = -999.0)
swrad = new_fc_sa.createVariable('SWRadAtm', np.float64,('time','hru'), fill_value = -999.0)
airpres = new_fc_sa.createVariable('airpres', np.float64,('time','hru'), fill_value = -999.0)
airtemp = new_fc_sa.createVariable('airtemp', np.float64,('time','hru'), fill_value = -999.0)
pptrate = new_fc_sa.createVariable('pptrate', np.float64,('time','hru'), fill_value = -999.0)
spechum = new_fc_sa.createVariable('spechum', np.float64,('time','hru'), fill_value = -999.0)
windspd = new_fc_sa.createVariable('windspd', np.float64,('time','hru'), fill_value = -999.0)
# give variables units
times.units = 'days since 1990-01-01 00:00:00'
ds.units = 'seconds'
lwrad.units = 'W m-2'
swrad.units = 'W m-2'
airpres.units = 'Pa'
airtemp.units = 'K'
pptrate.units = 'kg m-2 s-1'
spechum.units = 'g g-1'
windspd.units = 'm s-1'
# give variables value type
lwrad.vtype = 'scalarv'
swrad.vtype = 'scalarv'
airpres.vtype = 'scalarv'
airtemp.vtype = 'scalarv'
pptrate.vtype = 'scalarv'
spechum.vtype = 'scalarv'
windspd.vtype = 'scalarv'
# read out to compare with original
#for varname in new_fc_sa.variables.keys():
#    var = new_fc_sa.variables[varname]
#    print (varname, var.dtype, var.dimensions, var.shape)
#%% define hru id, time step (1hr), lat and lon
step = np.array([3600])

lat_sa = np.array([37.906914133])#33333) np.array(sbFD.variables['latitude'][:])
len_lat = np.repeat(lat_sa[:,np.newaxis], hru_num, axis=1); len_lat=len_lat.reshape(hru_num,)

long_sa = np.array([360 - 107.711322011])#11111)np.array(sbFD.variables['longitude'][:])
len_lon= np.repeat(long_sa[:,np.newaxis], hru_num, axis=1); len_lon=len_lon.reshape(hru_num,)
#lat_sb = np.array(sbFD.variables['latitude'][:])
#len_lat = np.repeat(lat_sb[:,np.newaxis], hru_num, axis=1); len_lat=len_lat.reshape(hru_num,)
#long_sb =np.array(sbFD.variables['longitude'][:])
#len_lon= np.repeat(long_sb[:,np.newaxis], hru_num, axis=1); len_lon=len_lon.reshape(hru_num,)
#%%
# assign newly created variables with lists of values from NLDAS and Sagehen data
hruid[:] = hruidxID 
lat[:] = len_lat
lon[:] = len_lon
ds[:] = step

new_ix = np.array(saTime) #new_ix = sbFD.variables['time'][:]
times[:] = new_ix

lwr_sa = np.array(sa_df['LWRadAtm'])  #lwr_sb = sbFD.variables['LWRadAtm'][:]
lwr_sa_hru = np.repeat(lwr_sa[:,np.newaxis], hru_num, axis=1)
lwrad[:] = lwr_sa_hru

swr_sa = np.array(sa_df['SWRadAtm'])
swr_sa_hru = np.repeat(swr_sa[:,np.newaxis], hru_num, axis=1)
swrad[:] = swr_sa_hru

ap_sa = np.array(sa_df['airpres'])
ap_sa_hru = np.repeat(ap_sa[:,np.newaxis], hru_num, axis=1)
airpres[:] = ap_sa_hru

at_sa = np.array(sa_df['airtemp'])
at_sa_hru = np.repeat(at_sa[:,np.newaxis], hru_num, axis=1) 
airtemp[:] = at_sa_hru

ws_sa = np.array(sa_df['windspd'])
ws_sa_hru = np.repeat(ws_sa[:,np.newaxis], hru_num, axis=1) 
windspd[:] = ws_sa_hru

sh_sa = np.array(sa_df['spechum'])
sh_sa_hru = np.repeat(sh_sa[:,np.newaxis], hru_num, axis=1) 
spechum[:] = sh_sa_hru

ppt_sa = np.array(sa_df['pptrate'])
ppt_sa_hru = np.repeat(ppt_sa[:,np.newaxis], hru_num, axis=1) 
pptrate[:] = ppt_sa_hru

#ppt_sa0 = np.array(sa_df['pptrate'])
#ppt_sa1 = []
#for cpi in range (np.size(ppt_sa0)):
#    if at_sa [cpi] <= 274 :
#        ppt_sa1.append(ppt_sa0[cpi]+ppt_sa0[cpi]*(1-((np.exp(4.61-0.04*((ws_sa[cpi])**1.75)))/100)))
#    else: ppt_sa1.append(ppt_sa0[cpi])
#ppt_sa = np.array(ppt_sa1)
#ppt_sa_hru = np.repeat(ppt_sa[:,np.newaxis], hru_num, axis=1) 
#pptrate[:] = ppt_sa_hru

#ppt_sa0 = np.array(sa_df['pptrate'])
#ppt_sa1 = []
#for cpi in range (np.size(ppt_sa0)):
#    if at_sa [cpi] <= 274 :
#        ppt_sa1.append(ppt_sa0[cpi]+ppt_sa0[cpi]*(1-((np.exp(4.61-0.16*(ws_sa[cpi]**1.28)))/100)))
#    else: ppt_sa1.append(ppt_sa0[cpi])
#ppt_sa = np.array(ppt_sa1)
#ppt_sa_hru = np.repeat(ppt_sa[:,np.newaxis], hru_num, axis=1) 
#pptrate[:] = ppt_sa_hru

#%%specific humididty calculations***********************************************
#at0_sb = sbFD.variables['airtemp'][:]
#
#e_t = (ap_sb * sh_sb)/0.622
#p_da = ap_sb - e_t
#e_star_t = 611*(np.exp((17.27*(at0_sb-273.15))/(at0_sb-273.15+237.3)))
#rh = e_t/e_star_t
#
#e_star_t2 = 611*(np.exp((17.27*(at0_sb+2-273.15))/(at0_sb+2-273.15+237.3)))
#e_star_t4 = 611*(np.exp((17.27*(at0_sb+4-273.15))/(at0_sb+4-273.15+237.3)))
#
#e_t2 = rh * e_star_t2
#e_t4 = rh * e_star_t4
#
#p_t2 = p_da + e_t2
#p_t4 = p_da + e_t4
#
#sh_t2 = 0.622 * e_t2 / p_t2
#sh_t4 = 0.622 * e_t4 / p_t4
#%%******************************************************************************
test = new_fc_sa.variables['pptrate'][:]

# close the file to write it
new_fc_sa.close()
#%%
testfd = Dataset("SwampAngel_forcing_sa2.nc")
#print testfd.file_format
# read out variables, data types, and dimensions of original forcing netcdf
for varname in testfd.variables.keys():
    var = testfd.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)









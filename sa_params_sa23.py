###       /bin/bash runTestCases_dockerD.sh   
# 2007 - 2008 as wet year for sensirivity analysis 1st step
import numpy as np
from netCDF4 import Dataset
import itertools
#%%  all parameters
#20,1, 'SHDFAC NROOT   RS      RGL      HS      SNUP  MAXALB   LAIMIN  LAIMAX   EMISSMIN EMISSMAX ALBEDOMIN ALBEDOMAX   Z0MIN    Z0MAX'
#7,      .80,   3,     40.,   100.,   36.35,   0.04,    70.,    0.52,   2.90,   .920,    .960,     .19,      .23,      .10,     .12,     'Grassland' 
p1 = [0.1] #LAIMIN
p2 = [1] #LAIMAX
p3 = [0.1] #winterSAI
p4 = [0.9] #summerLAI
p5 = [0.5] #rootingDepth
p6 = [0.5] #heightCanopyTop
p7 = [0.1] #heightCanopyBottom
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

#p18 = [0.6] #0.2, 0.4 , 0.6] #fixedThermalCond_snow

paramfile = Dataset("summa_zParamTrial_variableDecayRate_sa_sa2_j2330.nc",'w',format='NETCDF3_CLASSIC') #create new paramtrail.nc file

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
#%%
# function to create lists of each parameter, this will iterate through to make sure all combinations are covered
def param_fill(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17):#, p18): 
    b = list(itertools.product(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17))#, p18)) 
    p1l =[]; p2l =[]; p3l =[]; p4l=[]; p5l = []; p6l =[]; p7l =[]; p8l =[]; p9l =[]; p10l=[]; p11l =[]; p12l =[]; p13l =[]; p14l=[]; p15l = []; p16l =[]; p17l=[]#; p18l = []
    for tup in b:
        p1l.append(tup[0]); p2l.append(tup[1]); p3l.append(tup[2]); p4l.append(tup[3]); p5l.append(tup[4]); p6l.append(tup[5]); p7l.append(tup[6]); p8l.append(tup[7]); 
        p9l.append(tup[8]); p10l.append(tup[9]); p11l.append(tup[10]); p12l.append(tup[11]); p13l.append(tup[12]); p14l.append(tup[13]); p15l.append(tup[14]); p16l.append(tup[15]); p17l.append(tup[16])#; p18l.append(tup[17])
    return(p1l, p2l, p3l, p4l, p5l, p6l, p7l, p8l, p9l, p10l, p11l, p12l, p13l, p14l, p15l, p16l, p17l)#, p18l)  

# call the function on the parameters
valst1 = param_fill(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17)#, p18)  

#%% #create new paramtrail.nc file and adding vaiables to it --- summa_zParamTrial_variableDecayRate_test
hru = paramfile.createDimension('hru', None)
hidx = paramfile.createVariable('hruIndex', np.float64,('hru',)) # add hruIndex variable

param_nam_list = ['LAIMIN','LAIMAX','winterSAI','summerLAI','rootingDepth','heightCanopyTop','heightCanopyBottom','throughfallScaleSnow','newSnowDenMin',
                  'albedoDecayRate', 'albedoMaxVisible', 'albedoMinVisible', 'albedoMaxNearIR', 'albedoMinNearIR',  
                  'z0Snow', 'albedoRefresh', 'mw_exp']#, 'fixedThermalCond_snow'] 

for param in param_nam_list:
    paramfile.createVariable(param, np.float64,('hru',))

constant_params = ['frozenPrecipMultip','rootDistExp','theta_sat','theta_res','vGn_alpha','vGn_n','k_soil','critSoilWilting','critSoilTranspire']
for params in constant_params:
    paramfile.createVariable(params, np.float64,('hru',))
#paramfile.close()
#%% parameterTrial, Local attributes and initial conditions for senatore beck
pt = Dataset('summa_zParamTrial_variableDecayRate.nc')
la = Dataset('summa_zLocalAttributes_senatorSheltered.nc') #('settings/wrrPaperTestCases/figure07/summa_zLocalAttributes_riparianAspen.nc')
ic = Dataset('summa_zInitialCond.nc') #('settings/wrrPaperTestCases/figure07/summa_zInitialCond.nc')
for j in pt.variables:
    print j
#%% # add values for the constant variables in HRUs for parameter Trail file
for varname in pt.variables.keys():
    var = pt.variables[varname][0]
    c = np.full((hru_num,),var)
    try :
        paramfile.variables[varname][:]=c
    except IndexError: # size of data array does not conform to slice
        pass
#%% creating changing variables and adding values
j = 0 
for var in param_nam_list:
    paramfile.variables[var][:]=valst1[j]
    j=j+1
# don't forget the HRU Index!!
paramfile.variables['hruIndex'][:]=hruidxID

for varname in paramfile.variables.keys():
    var = paramfile.variables[varname]
    #print varname, var.dtype, var.dimensions, var.shape

#print paramfile.variables['hruIndex'][:]
paramfile.close()
#%% 
varcheck = Dataset ("C:\Users\HHS\summaTestCases_2.x\settings\wrrPaperTestCases\swampAngel\paramtrial\summa_zParamTrial_variableDecayRate_sa_sa2_j1110.nc")
print varcheck.variables['heightCanopyTop'][:]
#check2 =  varcheck.variables['albedoMaxNearIR'][:]
#I checked it in Check.py code
#%% # local attributes file
# create a new localAtribute file ---- summa_zLocalAttributes_swampAngel_vtest
local_atrbt = Dataset("summa_zLocalAttributes_swampAngel_sa2.nc",'w',format='NETCDF3_CLASSIC')
# define dimensions 
hru = local_atrbt.createDimension('hru', hru_num) 
time = local_atrbt.createDimension('gru', 1)
# define variables
h2gid = local_atrbt.createVariable('hru2gruId', np.int32,('hru',))
dhruindx = local_atrbt.createVariable('downHRUindex', np.int32,('hru',))
slopeindx = local_atrbt.createVariable('slopeTypeIndex', np.int32,('hru',))
soilindx = local_atrbt.createVariable('soilTypeIndex', np.int32,('hru',))
vegindx = local_atrbt.createVariable('vegTypeIndex', np.int32,('hru',))
mh = local_atrbt.createVariable('mHeight', np.float64,('hru',))
cl = local_atrbt.createVariable('contourLength', np.float64,('hru',))
tanslope = local_atrbt.createVariable('tan_slope', np.float64,('hru',))
elev = local_atrbt.createVariable('elevation', np.float64,('hru',))
lon = local_atrbt.createVariable('longitude', np.float64,('hru',))
lat = local_atrbt.createVariable('latitude', np.float64,('hru',))
hruarea = local_atrbt.createVariable('HRUarea', np.float64,('hru',))
hruid = local_atrbt.createVariable('hruId', np.int32,('hru',))
gruid = local_atrbt.createVariable('gruId', np.int32,('gru',))
# give variables units
mh.units = 'm'
cl.units = 'm'
tanslope.units = 'm m-1'
elev.units = 'm'
lat.units = 'decimal degree north'
lon.units = 'decimal degree east'
hruarea.units = 'm^2'
#%% # add values for the constant variables in HRUs for local atribute file
for varname in la.variables.keys():
    var = la.variables[varname][0]
    #print var
    c2 = np.full((hru_num,),var)
    #print c2
    try :
        local_atrbt.variables[varname][:]=c2
    except IndexError: # size of data array does not conform to slice
        pass
    #local_atrbt.variables[varname][:]=c2

#%% # get the hru, gru, and hru2gru in local_atribute file
newgru = np.array([11111])
local_atrbt.variables['gruId'][:] = newgru

c3 = np.repeat(newgru[:,np.newaxis], hru_num, axis=1); newlad = c3.reshape(hru_num,)
local_atrbt.variables['hru2gruId'][:] = c3

local_atrbt.variables['hruId'][:] = hruidxID

newMH = np.array([3.5])
newLat = np.array([37.906914133])#33333) np.array(sbFD.variables['latitude'][:])
newLon = np.array([360 - 107.711322011])#11111)np.array(sbFD.variables['longitude'][:])
newElev = np.array([3371])

len_lat = np.repeat(newLat[:,np.newaxis], hru_num, axis=1); len_lat=len_lat.reshape(hru_num,)
len_lon= np.repeat(newLon[:,np.newaxis], hru_num, axis=1); len_lon=len_lon.reshape(hru_num,)
len_MH= np.repeat(newMH[:,np.newaxis], hru_num, axis=1); len_MH=len_MH.reshape(hru_num,)
len_Elev= np.repeat(newElev[:,np.newaxis], hru_num, axis=1); len_Elev=len_Elev.reshape(hru_num,)

lat[:] = len_lat
lon[:] = len_lon
mh[:] = len_MH
elev[:] = len_Elev

local_atrbt.close()
#%%
lacheck = Dataset('summa_zLocalAttributes_swampAngel_sa2.nc')

print lacheck.variables['mHeight'][:]
#for j in laCheck.variables:
#    print j
for varname in lacheck.variables.keys():
    var = lacheck.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)    
#lacheck.close()
#%% # initial conditions file. summa_zInitialCond_vtest

in_condi = Dataset("summa_zInitialCond_sa2.nc",'w',format='NETCDF3_CLASSIC')
#print ic.variables.keys()

# define dimensions 
midtoto = in_condi.createDimension('midToto',8)
midsoil = in_condi.createDimension('midSoil',8)
idctoto = in_condi.createDimension('ifcToto',9)
scalarv = in_condi.createDimension('scalarv', 1)
# this is the number you will change to the number of HRU's from your param trial file
hrud = in_condi.createDimension('hru', hru_num)
# define variables
mlvfi = in_condi.createVariable('mLayerVolFracIce', np.float64, ('midToto', 'hru'))
scat = in_condi.createVariable('scalarCanairTemp', np.float64, ('scalarv', 'hru'))
nsnow = in_condi.createVariable('nSnow', np.int32, ('scalarv', 'hru'))
ilh = in_condi.createVariable('iLayerHeight', np.float64, ('ifcToto', 'hru'))
mlmh = in_condi.createVariable('mLayerMatricHead', np.float64, ('midSoil', 'hru'))
ssa = in_condi.createVariable('scalarSnowAlbedo', np.float64, ('scalarv', 'hru'))
dti = in_condi.createVariable('dt_init', np.float64, ('scalarv', 'hru'))
mlt = in_condi.createVariable('mLayerTemp', np.float64, ('midToto', 'hru'))
ssmp = in_condi.createVariable('scalarSfcMeltPond', np.float64, ('scalarv', 'hru'))
sct = in_condi.createVariable('scalarCanopyTemp', np.float64, ('scalarv', 'hru'))
ssd = in_condi.createVariable('scalarSnowDepth', np.float64, ('scalarv', 'hru'))
nsoil = in_condi.createVariable('nSoil', np.int32, ('scalarv', 'hru'))
sswe = in_condi.createVariable('scalarSWE', np.float64, ('scalarv', 'hru'))
scl = in_condi.createVariable('scalarCanopyLiq', np.float64, ('scalarv', 'hru'))
mlvf = in_condi.createVariable('mLayerVolFracLiq', np.float64, ('midToto', 'hru'))
mld = in_condi.createVariable('mLayerDepth', np.float64, ('midToto', 'hru'))
sci = in_condi.createVariable('scalarCanopyIce', np.float64, ('scalarv', 'hru'))
sas = in_condi.createVariable('scalarAquiferStorage', np.float64, ('scalarv', 'hru'))
#%% # add values for the intial condition variables in HRUs

for varname in ic.variables.keys():
    infovar = ic.variables[varname]
    var = ic.variables[varname][:]
    cic = np.repeat(var[:,np.newaxis], hru_num, axis=1); newic = cic.reshape(infovar.shape[0],hru_num)
    in_condi.variables[varname][:]=newic

print in_condi.variables['scalarSnowAlbedo'][:]

in_condi.close()
#%%
iccheck = Dataset("summa_zInitialCond_sa2.nc")
#for varname in iccheck.variables.keys():
#    var = iccheck.variables[varname]
#    print (varname, var.dtype, var.dimensions, var.shape)
print iccheck.variables['nSoil'][:]
#%%





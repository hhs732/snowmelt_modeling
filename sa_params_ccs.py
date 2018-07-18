###       /bin/bash runTestCases_docker_ccs.sh   snwDensity snwDensity
# 2007 - 2008 as wet year for sensirivity analysis 1st step
import numpy as np
from netCDF4 import Dataset
import itertools
#%%  all parameters
#p1 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1] #LAIMIN
#p2 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] #LAIMAX
#p3 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1] #winterSAI
#p4 = [0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9] #summerLAI
#p5 = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5] #rootingDepth
#p6 = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5] #heightCanopyTop
#p7 = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01] #heightCanopyBottom
#p8 = [0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89] #throughfallScaleSnow
#p9 = [55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55,55] #newSnowDenMin 
#
#p10 =[500000,1300000,1300000,1300000,1300000,500000,1000000,1300000,500000,1300000,1300000,500000,500000,1000000,500000,500000,
#      500000,1000000,1000000,1300000,500000,1300000,1300000,500000,1000000,1000000] #[500000, 1000000, 1300000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
#p11 =[0.9,0.8,0.8,0.9,0.9,0.94,0.9,0.8,0.94,0.8,0.9,0.9,0.94,0.8,0.9,0.9,0.94,0.94,0.8,0.94,0.8,0.8,0.9,0.9,0.94,0.9] #[0.8, 0.9, 0.94] #albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
#p12 =[0.68,0.6,0.6,0.6,0.6,0.74,0.68,0.6,0.68,0.6,0.6,0.68,0.6,0.68,0.6,0.74,0.74,0.6,0.68,0.6,0.68,0.68,0.68,0.68,0.68,0.6] #[0.6, 0.68, 0.74] #albedoMinVisible |       0.7500 |       0.5000 |       0.7500
#p13 =[0.7,0.7,0.65,0.55,0.65,0.55,0.7,0.7,0.65,0.65,0.55,0.7,0.65,0.7,0.65,0.7,0.65,0.55,0.7,0.55,0.7,0.65,0.55,0.55,0.65,0.65] #[0.55, 0.65, 0.7] #albedoMaxNearIR |       0.6500 |       0.5000 |       0.7500
#p14 =[0.3,0.2,0.2,0.2,0.2,0.3,0.2,0.2,0.4,0.3,0.3,0.3,0.4,0.2,0.4,0.3,0.3,0.3,0.3,0.3,0.4,0.2,0.2,0.4,0.2,0.3] #[0.2, 0.3, 0.4] #albedoMinNearIR  |       0.3000 |       0.1500 |       0.4500
#
#p15 =[0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002,0.002] #[0.002] #[0.001, 0.002] #z0Snow
#p16 =[1,1,1,1,3,3,6,6,6,6,6,1,1,1,1,3,3,3,3,3,3,3,3,3,6,6] #[6]# 1, 3, 6] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
#p17 =[3,4,4,4,3,4,3,3,4,4,4,2,3,3,3,3,3,3,4,4,4,4,4,4,2,3] #[4] #2, 3, 4] #mw_exp exponent for meltwater flow
#p18 =[0.4,0.2,0.4,0.4,0.2,0.6,0.2,0.6,0.2,0.4,0.4,0.6,0.4,0.4,0.6,0.4,0.4,0.6,0.2,0.2,0.4,0.4,0.4,0.6,0.4,0.4] #[0.6] #0.2, 0.4 , 0.6] #fixedThermalCond_snow
#21121sc1233
#22111sc1233
#21121sc2233
#22111sc2233
#22121sc2331
p1 = [0.1,0.1,0.1,0.1,0.1] #LAIMIN
p2 = [1,1,1,1,1] #LAIMAX
p3 = [0.1,0.1,0.1,0.1,0.1] #winterSAI
p4 = [0.9,0.9,0.9,0.9,0.9] #summerLAI
p5 = [0.5,0.5,0.5,0.5,0.5] #rootingDepth
p6 = [0.5,0.5,0.5,0.5,0.5] #heightCanopyTop
p7 = [0.01,0.01,0.01,0.01,0.01] #heightCanopyBottom
p8 = [0.89,0.89,0.89,0.89,0.89] #throughfallScaleSnow
p9 = [55,55,55,55,55] #newSnowDenMin 

p10 =[1000000,1000000,1000000,1000000,1000000] #[500000, 1000000, 1300000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p11 =[0.8,0.9,0.8,0.9,0.9] #[0.8, 0.9, 0.94] #albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p12 =[0.6,0.6,0.6,0.6,0.6] #[0.6, 0.68, 0.74] #albedoMinVisible |       0.7500 |       0.5000 |       0.7500
p13 =[0.65,0.55,0.65,0.55,0.65] #[0.55, 0.65, 0.7] #albedoMaxNearIR |       0.6500 |       0.5000 |       0.7500
p14 =[0.2,0.2,0.2,0.2,0.2] #[0.2, 0.3, 0.4] #albedoMinNearIR  |       0.3000 |       0.1500 |       0.4500

p15 =[0.001,0.001,0.002,0.002,0.002] #[0.002] #[0.001, 0.002] #z0Snow
p16 =[3,3,3,3,6] #[6]# 1, 3, 6] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
p17 =[4,4,4,4,4] #[4] #2, 3, 4] #mw_exp exponent for meltwater flow
p18 =[0.6,0.6,0.6,0.6,0.2]
paramfile = Dataset("summa_zParamTrial_variableDecayRate_sa_ccs_sc.nc",'w',format='NETCDF3_CLASSIC') #create new paramtrail.nc file

hruidxID = list(np.arange(101,106))
hru_num = np.size(hruidxID)
#%% #create new paramtrail.nc file and adding vaiables to it --- summa_zParamTrial_variableDecayRate_test
hru = paramfile.createDimension('hru', None)
hidx = paramfile.createVariable('hruIndex', np.float64,('hru',)) # add hruIndex variable

param_nam_list = ['LAIMIN','LAIMAX','winterSAI','summerLAI','rootingDepth','heightCanopyTop','heightCanopyBottom','throughfallScaleSnow','newSnowDenMin',
                  'albedoDecayRate', 'albedoMaxVisible', 'albedoMinVisible', 'albedoMaxNearIR', 'albedoMinNearIR',  
                  'z0Snow', 'albedoRefresh', 'mw_exp', 'fixedThermalCond_snow'] 

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
#for j in pt.variables:
#    print j
#%% # add values for the constant variables in HRUs for parameter Trail file
for varname in pt.variables.keys():
    var = pt.variables[varname][0]
    c = np.full((hru_num,),var)
    try :
        paramfile.variables[varname][:]=c
    except IndexError: # size of data array does not conform to slice
        pass
#%% creating changing variables and adding values
# add values for the changing variables in HRUs

paramfile.variables['LAIMIN'][:]=p1
paramfile.variables['LAIMAX'][:]=p2
paramfile.variables['winterSAI'][:]=p3
paramfile.variables['summerLAI'][:]=p4
paramfile.variables['rootingDepth'][:]=p5
paramfile.variables['heightCanopyTop'][:]=p6

paramfile.variables['heightCanopyBottom'][:]=p7
paramfile.variables['throughfallScaleSnow'][:]=p8
paramfile.variables['newSnowDenMin'][:]=p9
paramfile.variables['albedoDecayRate'][:]=p10
paramfile.variables['albedoMaxVisible'][:]=p11
paramfile.variables['albedoMinVisible'][:]=p12

paramfile.variables['albedoMaxNearIR'][:]=p13
paramfile.variables['albedoMinNearIR'][:]=p14
paramfile.variables['z0Snow'][:]=p15
paramfile.variables['albedoRefresh'][:]=p16

paramfile.variables['mw_exp'][:]=p17
paramfile.variables['fixedThermalCond_snow'][:]=p18

paramfile.variables['hruIndex'][:]=hruidxID

for varname in paramfile.variables.keys():
    var = paramfile.variables[varname]
    print varname, var.dtype, var.dimensions, var.shape

print paramfile.variables['hruIndex'][:]
paramfile.close()
#%% 
#varcheck = Dataset ('summa_zParamTrial_variableDecayRate_sa_sa2_s2333.nc')
#print varcheck.variables['theta_res'][:]
#check2 =  varcheck.variables['albedoMaxNearIR'][:]
#I checked it in Check.py code
#%% # local attributes file
# create a new localAtribute file ---- summa_zLocalAttributes_swampAngel_vtest
local_atrbt = Dataset("summa_zLocalAttributes_swampAngel_ccs_sc.nc",'w',format='NETCDF3_CLASSIC')
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
newgru = np.array([111])
local_atrbt.variables['gruId'][:] = newgru

c3 = np.repeat(newgru[:,np.newaxis], hru_num, axis=1); newlad = c3.reshape(hru_num,)
local_atrbt.variables['hru2gruId'][:] = c3

local_atrbt.variables['hruId'][:] = hruidxID

#print local_atrbt.variables['hruId'][:]

local_atrbt.close()
#%%
lacheck = Dataset('summa_zLocalAttributes_swampAngel_ccs_sc.nc')

print lacheck.variables['vegTypeIndex'][:]
#for j in laCheck.variables:
#    print j
for varname in lacheck.variables.keys():
    var = lacheck.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)    
lacheck.close()
#%% # initial conditions file. summa_zInitialCond_vtest

in_condi = Dataset("summa_zInitialCond_ccs_sc.nc",'w',format='NETCDF3_CLASSIC')
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
iccheck = Dataset("summa_zInitialCond_ccs_sc.nc")
#for varname in iccheck.variables.keys():
#    var = iccheck.variables[varname]
#    print (varname, var.dtype, var.dimensions, var.shape)
print iccheck.variables['nSoil'][:]
#%%
#1 2 2 3 2 mc1 1 2 2
#3 1 1 3 1 mc1 1 3 1
#3 1 1 2 1 mc1 1 3 2
#3 2 1 1 1 mc1 1 3 2
#3 2 1 2 1 mc1 2 2 1
#1 3 3 1 2 mc1 2 3 3
#2 2 2 3 1 mc1 3 2 1
#3 1 1 3 1 mc1 3 2 3
#1 3 2 2 3 mc1 3 3 1
#3 1 1 2 2 mc1 3 3 2
#3 2 1 1 2 mc1 3 3 2
#1 2 2 3 2 mc2 1 1 3
#1 3 1 2 3 mc2 1 2 2
#2 1 2 3 1 mc2 1 2 2
#1 2 1 2 3 mc2 1 2 3
#1 2 3 3 2 mc2 2 2 2
#1 3 3 2 2 mc2 2 2 2
#2 3 1 1 2 mc2 2 2 3
#2 1 2 3 2 mc2 2 3 1
#3 3 1 1 2 mc2 2 3 1
#1 1 2 3 3 mc2 2 3 2
#3 1 2 2 1 mc2 2 3 2
#3 2 2 1 1 mc2 2 3 2
#1 2 2 1 3 mc2 2 3 3
#2 3 2 2 1 mc2 3 1 2
#2 2 1 2 2 mc2 3 2 2          




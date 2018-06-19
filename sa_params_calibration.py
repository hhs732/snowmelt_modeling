###       /bin/bash runTestCases_docker.sh 
# 2007 - 2008 as wet year for sensirivity analysis 1st step
import numpy as np
from netCDF4 import Dataset
import itertools
#%% Vegetation Cover
#20,1, 'SHDFAC NROOT   RS      RGL      HS      SNUP  MAXALB   LAIMIN  LAIMAX   EMISSMIN EMISSMAX ALBEDOMIN ALBEDOMAX   Z0MIN    Z0MAX'
#7,      .80,   3,     40.,   100.,   36.35,   0.04,    70.,    0.52,   2.90,   .920,    .960,     .19,      .23,      .10,     .12,     'Grassland' 
p1 = [0.1, 0.1, 0.1, 0.1] #LAIMIN
p2 = [1, 1, 1, 1] #LAIMAX
#EMISSMIN = 0.92
#EMISSMAX = 0.96
#ALBEDOMIN = 0.19
#ALBEDOMAX = 0.23
#Z0MIN = 0.10
#Z0MAX = 0.12

p3 = [0.1, 0.1, 0.1, 0.1] #winterSAI
p4 = [1, 1, 1, 1]  #summerLAI
p5 = [0.5, 0.5, 0.5, 0.5] #rootingDepth
p6 = [0.5, 0.5, 0.5, 0.5] #heightCanopyTop
p7 = [0.01, 0.01] #heightCanopyBottom
p8 = [0.89, 0.89, 0.89, 0.89] #throughfallScaleSnow

p9 = [0.4, 0.4, 0.4, 0.4] #0.2, 0.35 , 0.6] #fixedThermalCond_snow

p10 = [1000000,1000000, 1000000,1000000] #albedoDecayRate ??????????????????????
p11 = [0.85, 0.85, 0.85, 0.85] #albedoMaxVisible 
p12 = [0.6, 0.6, 0.6, 0.6] #albedoMinVisible  
p13 = [0.55, 0.55, 0.55, 0.55] #albedoMaxNearIR
p14 = [0.2, 0.2, 0.2, 0.2] #albedoMinNearIR
p15 = [1, 7, 1, 7] #albedoRefresh

p16 = [2, 2, 4, 4] #2, 3, 4] #mw_exp exponent for meltwater flow

p17 = [57, 57, 57, 57] #51, 60, 70] #(m)newSnowDenMin 60 80 100     46 47 48
#p17 = [88, 104, 120] #51, 70, 105] #(c)constSnowDen 70.00, 100.0, 170.0    55 56 57 
p18 = [0.001, 0.001, 0.001, 0.001] #z0Snow

#p21 = [0.700, 1.000, 1.500] #Mahrt87_eScale  

hruidxID = list(np.arange(101,105))
hru_num = np.size(hruidxID)

#%% #create new paramtrail.nc file and adding vaiables to it --- summa_zParamTrial_variableDecayRate_test
paramfile = Dataset("summa_zParamTrial_variableDecayRate_sa_calibration.nc",'w',format='NETCDF3_CLASSIC') #create new paramtrail.nc file

hru = paramfile.createDimension('hru', None)
hidx = paramfile.createVariable('hruIndex', np.float64,('hru',)) # add hruIndex variable

# add 2 new variables to paramfile that we are going to change (use the param list to make the variables for the netCDF)
param_nam_list = ['LAIMIN', 'LAIMAX', 'winterSAI', 'summerLAI', 'rootingDepth', 'heightCanopyTop', 'heightCanopyBottom', 'throughfallScaleSnow', 
                  'fixedThermalCond_snow', 
                  'albedoDecayRate', 'albedoMaxVisible', 'albedoMinVisible', 'albedoMaxNearIR', 'albedoMinNearIR', 'albedoRefresh',
                  'mw_exp', 'newSnowDenMin', 'z0Snow'] 
for param in param_nam_list:
    paramfile.createVariable(param, np.float64,('hru',))

# add any variable to paramfile that we are NOT going to change. Any variable that you are not going to change, it should include in this list
# Ava case in constant_parameter : 'windReductionParam'
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
    # add values for the constant variables in HRUs
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

paramfile.variables['fixedThermalCond_snow'][:]=p9

paramfile.variables['albedoDecayRate'][:]=p10
paramfile.variables['albedoMaxVisible'][:]=p11
paramfile.variables['albedoMinVisible'][:]=p12
paramfile.variables['albedoMaxNearIR'][:]=p13
paramfile.variables['albedoMinNearIR'][:]=p14
paramfile.variables['albedoRefresh'][:]=p15

paramfile.variables['mw_exp'][:]=p16
paramfile.variables['newSnowDenMin'][:]=p17
paramfile.variables['z0Snow'][:]=p18

paramfile.variables['hruIndex'][:]=hruidxID

for varname in paramfile.variables.keys():
    var = paramfile.variables[varname]
    print varname, var.dtype, var.dimensions, var.shape

print paramfile.variables['hruIndex'][:]
paramfile.close()
#%% 8         0.039       0.387      -2.667       1.449 4.43084e-06        4.90       0.065        1.34       0.435       0.383       0.218     3.47E-5    0.805E-5       0.114        0.60    'SANDY LOAM'
# 'theta_res   theta_sat   vGn_alpha       vGn_n      k_soil          BB      DRYSMC          HC      MAXSMC      REFSMC      SATPSI       SATDK       SATDW      WLTSMC 
varcheck = Dataset ('summa_zParamTrial_variableDecayRate_sa_calibration.nc')
print varcheck.variables['theta_res'][:]
#check2 =  varcheck.variables['albedoMaxNearIR'][:]
#check3 =  varcheck.variables['albedoMinVisible'][:]
#I checked it in Check.py code
#%% # local attributes file
# create a new localAtribute file ---- summa_zLocalAttributes_swampAngel_vtest
local_atrbt = Dataset("summa_zLocalAttributes_swampAngel_calibration.nc",'w',format='NETCDF3_CLASSIC')
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
lacheck = Dataset('summa_zLocalAttributes_swampAngel_calibration.nc')

print lacheck.variables['vegTypeIndex'][:]
#for j in laCheck.variables:
#    print j
for varname in lacheck.variables.keys():
    var = lacheck.variables[varname]
    print (varname, var.dtype, var.dimensions, var.shape)    
lacheck.close()
#%% # initial conditions file. summa_zInitialCond_vtest

in_condi = Dataset("summa_zInitialCond_calibration.nc",'w',format='NETCDF3_CLASSIC')
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
iccheck = Dataset("summa_zInitialCond_calibration.nc")
#for varname in iccheck.variables.keys():
#    var = iccheck.variables[varname]
#    print (varname, var.dtype, var.dimensions, var.shape)
print iccheck.variables['scalarSWE'][:]




#%%





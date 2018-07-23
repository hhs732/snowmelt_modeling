#%%
#SC
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
#%%
#Sj
p1 = [0.1] #LAIMIN
p2 = [1] #LAIMAX
p3 = [0.1] #winterSAI
p4 = [0.9] #summerLAI
p5 = [0.5] #rootingDepth
p6 = [0.5] #heightCanopyTop
p7 = [0.01] #heightCanopyBottom
p8 = [0.89] #throughfallScaleSnow
p9 = [55] #newSnowDenMin 


p10 =[500000] #[500000, 1000000, 1300000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p11 =[0.94] #[0.8, 0.9, 0.94] #albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p12 =[0.6] #[0.6, 0.68, 0.74] #albedoMinVisible |       0.7500 |       0.5000 |       0.7500
p13 =[0.7] #[0.55, 0.65, 0.7] #albedoMaxNearIR |       0.6500 |       0.5000 |       0.7500
p14 =[0.3] #[0.2, 0.3, 0.4] #albedoMinNearIR  |       0.3000 |       0.1500 |       0.4500

p15 =[0.002] #[0.002] #[0.001, 0.002] #z0Snow
p16 =[3] #[6]# 1, 3, 6] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
p17 =[4] #[4] #2, 3, 4] #mw_exp exponent for meltwater flow
#p18 =[0.6,0.6,0.6,0.6,0.2]
paramfile = Dataset("summa_zParamTrial_variableDecayRate_sa_ccs_sj.nc",'w',format='NETCDF3_CLASSIC') #create new paramtrail.nc file

hruidxID = list(np.arange(101,102))
hru_num = np.size(hruidxID)
#%%
#lj
p1 = [0.1,0.1] #LAIMIN
p2 = [1,1] #LAIMAX
p3 = [0.1,0.1] #winterSAI
p4 = [0.9,0.9] #summerLAI
p5 = [0.5,0.5] #rootingDepth
p6 = [0.5,0.55] #heightCanopyTop
p7 = [0.01,0.01] #heightCanopyBottom
p8 = [0.89,0.89] #throughfallScaleSnow
p9 = [55,55] #newSnowDenMin

p10 =[500000,500000] #[500000, 1000000, 1300000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p11 =[0.94,0.94] #[0.8, 0.9, 0.94] #albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p12 =[0.74,0.68] #[0.6, 0.68, 0.74] #albedoMinVisible |       0.7500 |       0.5000 |       0.7500
p13 =[0.65, 0.7] #[0.55, 0.65, 0.7] #albedoMaxNearIR |       0.6500 |       0.5000 |       0.7500
p14 =[0.2,0.2] #[0.2, 0.3, 0.4] #albedoMinNearIR  |       0.3000 |       0.1500 |       0.4500

p15 =[0.001,0.002] #[0.002] #[0.001, 0.002] #z0Snow
p16 =[3,1] #[6]# 1, 3, 6] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
p17 =[4,4] #[4] #2, 3, 4] #mw_exp exponent for meltwater flow
#p18 =[0.6,0.6,0.6,0.6,0.2]
paramfile = Dataset("summa_zParamTrial_variableDecayRate_sa_ccs_lj.nc",'w',format='NETCDF3_CLASSIC') #create new paramtrail.nc file

hruidxID = list(np.arange(101,103))
hru_num = np.size(hruidxID)
#%%
p1 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1] #LAIMIN
p2 = [1,1,1,1,1,1,1,1] #LAIMAX
p3 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1] #winterSAI
p4 = [0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9] #summerLAI
p5 = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5] #rootingDepth
p6 = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5] #heightCanopyTop
p7 = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01] #heightCanopyBottom
p8 = [0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89] #throughfallScaleSnow
p9 = [55,55,55,55,55,55,55,55] #newSnowDenMin 
#13 2 2 2mj11 30
#12 3 3 2mj13 30
#22 2 2 1mj13 30
#21 2 3 1mj21 30
#23 1 3 1mj22 10
#22 2 3 1mj23 20
#22 1 2 2mj23 30
#33 2 1 1mj23 30

p10 =[500000,500000,1000000,1000000,1000000,1000000,1000000,1300000] #[500000, 1000000, 1300000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p11 =[0.94,0.9,0.9,0.8,0.94,0.9,0.9,0.94] #[0.8, 0.9, 0.94] #albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p12 =[0.68,0.74,0.68,0.68,0.6,0.68,0.6,0.68] #[0.6, 0.68, 0.74] #albedoMinVisible |       0.7500 |       0.5000 |       0.7500
p13 =[0.65,0.7,0.65,0.7,0.7,0.7,0.65,0.55] #[0.55, 0.65, 0.7] #albedoMaxNearIR |       0.6500 |       0.5000 |       0.7500
p14 =[0.3,0.3,0.2,0.2,0.2,0.2,0.3,0.2] #[0.2, 0.3, 0.4] #albedoMinNearIR  |       0.3000 |       0.1500 |       0.4500

p15 =[0.001,0.001,0.001,0.002,0.002,0.002,0.002,0.002] #[0.002] #[0.001, 0.002] #z0Snow
p16 =[1,6,6,1,3,6,6,6] #[6]# 1, 3, 6] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
p17 =[4,4,4,4,2,3,4,4] #[4] #2, 3, 4] #mw_exp exponent for meltwater flow
#p18 =[0.6,0.6,0.6,0.6,0.2]
paramfile = Dataset("summa_zParamTrial_variableDecayRate_sa_ccs_mj.nc",'w',format='NETCDF3_CLASSIC') #create new paramtrail.nc file

hruidxID = list(np.arange(101,109))
hru_num = np.size(hruidxID)
#%%
#lc
p1 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1] #LAIMIN
p2 = [1,1,1,1,1,1,1,1] #LAIMAX
p3 = [0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1] #winterSAI
p4 = [0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9] #summerLAI
p5 = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5] #rootingDepth
p6 = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5] #heightCanopyTop
p7 = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01] #heightCanopyBottom
p8 = [0.89,0.89,0.89,0.89,0.89,0.89,0.89,0.89] #throughfallScaleSnow
p9 = [55,55,55,55,55,55,55,55] #newSnowDenMin 
#
#12 2 3 1lc1 1 23
#12 3 3 1lc1 1 31
#23 1 1 1lc1 2 32
#13 2 3 1lc2 1 22
#12 2 3 1lc2 1 23
#13 2 2 2lc2 2 31
#22 1 2 1lc2 2 31
#21 1 2 1lc2 2 33


p10 =[500000,500000,1000000,500000,500000,500000,1000000,1000000] #[500000, 1000000, 1300000] ##albedoDecayRate |       1.0d+6 |       0.1d+6 |       5.0d+6 
p11 =[0.9,0.9,0.94,0.94,0.9,0.94,0.9,0.8] #[0.8, 0.9, 0.94] #albedoMaxVisible |       0.9500 |       0.7000 |       0.9500
p12 =[0.68,0.74,0.6,0.68,0.68,0.68,0.6,0.6] #[0.6, 0.68, 0.74] #albedoMinVisible |       0.7500 |       0.5000 |       0.7500

p13 =[0.7,0.7,0.55,0.7,0.7,0.65,0.65,0.65] #[0.55, 0.65, 0.7] #albedoMaxNearIR |       0.6500 |       0.5000 |       0.7500
p14 =[0.2,0.2,0.2,0.2,0.2,0.3,0.2,0.2] #[0.2, 0.3, 0.4] #albedoMinNearIR  |       0.3000 |       0.1500 |       0.4500

p15 =[0.001,0.001,0.001,0.002,0.002,0.002,0.002,0.002] #[0.002] #[0.001, 0.002] #z0Snow
p16 =[1,1,3,1,1,3,3,3] #[6]# 1, 3, 6] #albedoRefresh |       1.0000 |       1.0000 |      10.0000
p17 =[3,4,4,3,3,4,4,4] #[4] #2, 3, 4] #mw_exp exponent for meltwater flow
p18 =[0.6,0.2,0.4,0.4,0.6,0.2,0.2,0.6]
paramfile = Dataset("summa_zParamTrial_variableDecayRate_sa_ccs_lc.nc",'w',format='NETCDF3_CLASSIC') #create new paramtrail.nc file

hruidxID = list(np.arange(101,109))
hru_num = np.size(hruidxID)
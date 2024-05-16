#!/usr/bin/env python
# coding: utf-8

# # Spin Test

# In[1]:


import numpy as np
from netneurotools.freesurfer import find_parcel_centroids
from netneurotools.stats import gen_spinsamples
from scipy.stats import pearsonr, spearmanr
import scipy.io as sio

# In[2]:


lh_annot= r"D:\Nuova cartella\DOCs\results\dynamicINTrepertoire\spintest_Mass\lh.HCP-MMP1.annot" #substitute with your annot files
rh_annot= r"D:\Nuova cartella\DOCs\results\dynamicINTrepertoire\spintest_Mass\rh.HCP-MMP1.annot"

N_spins=10000
spin_seed=1234 #for consistency


#%% compute nulls


# Compute centroids
centroids, hemiid = find_parcel_centroids(lhannot =lh_annot, rhannot=rh_annot, method='surface', version='fsaverage')

centroids = np.delete(centroids, [0, 181], axis=0)
hemiid = np.delete(hemiid, [0, 181])

 # Generate spin-test nulls (method can be 'original', 'hungarian' or 'vasa') - may take a while if nparc > 500
nulls=gen_spinsamples(centroids, hemiid, n_rotate = N_spins, method = 'vasa', seed = spin_seed, return_cost = False )


#%% define function


def spin_test(x,y,nulls,N_spins):
    
    r_emp = pearsonr(x,y)[0]
    r_null = np.empty(N_spins)
    
    # Compute p-val
    for i in range(N_spins):
        r_null[i] = pearsonr(x,y[nulls[:,i]])[0]
    
    p_perm = (1 + sum(abs(r_null - np.nanmean(r_null)) > abs(r_emp - np.nanmean(r_null)))) / (N_spins + 1)
    
    return r_emp, p_perm 



#%% apply test

    n_maps = 10 # substitute here the number of maps
mapsnmyelination = sio.loadmat() #load myelination map and MEG maps array
Cs = mapsnmyelination["C"]
x = Cs
y = mapsnmyelination["myelinationIndices"] 

for i in range(n_maps):
    [r_emp,p_perm] = spin_test(np.squeeze(x[i,:]),np.squeeze(y),nulls,N_spins)
    print(p_perm, r_emp)

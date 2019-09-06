# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 15:06:06 2017

@author: zehaojin
"""

##########only allows band=2,tau=1e11.0

import pickle,pyatomdb, numpy as np, matplotlib.pyplot as plt
from constants import version

###################
rmf="response/west.rmf"
arf="response/west.arf"
###################
'''
#################To see the energy bins in emf file
ret=pyatomdb.spectrum.get_response_ebins(rmf)
print ret
################One of the results: energy bin size=0.01 keV
'''

dust_perc=input('enter the initial dust percentage, e.g., 95:\n')
print 'percentage of dust=',dust_perc


band=2
print 'band=',band

precision=5.0
print 'precision=',precision

te=input('enter temperature(currently 1.0 or 2.0 or 4.0 or 10.0):\n')
te=float(te)
print 'temperature %s keV' %(str(te))

print 'initializing...'


filename='spec_data_v%s_p%s_b%s_te%s_d%s_tau11.0' %(str(version),str(precision),str(band),str(te),str(dust_perc))
raw_spec=pickle.load(open('data_package/'+filename,'rb'))



#####create energy bins(0.01keV)and array of its mid value
energybins=np.linspace(6.35,6.74,40)
energybins_midvalue=energybins+0.005

###########plot raw input data
plt.figure('%skeV_Chandra_%sdust_Raw_input_data:Energy_vs_Emissivity' %(str(te),str(dust_perc)))
plt.title('Energy vs Emissivity(energy bin=0.01keV) Te=%skeV %sdust 11.0logtau\n' %(str(te),str(dust_perc)))
plt.bar(energybins,raw_spec,0.01, color='b',alpha=0.7,edgecolor='y')
plt.xlabel('Energy (keV)')
plt.ylabel('Emissivity $({\mathrm{erg} \mathrm{cm}^3 \mathrm{s}^{-1}} \mathrm{keV}^{-1})$')
plt.savefig('results/response/%skeV_Chandra_%sdust_Raw_input_data:Energy_vs_Emissivity.eps' %(str(te),str(dust_perc)))

#######cancel the erg in raw_spec
for raw_index in range(len(raw_spec)):
    raw_spec[raw_index]=raw_spec[raw_index]*6.242e8/energybins_midvalue[raw_index]
###plot   
plt.figure('%skeV_Chandra_%sdust_Processed_input_data:Energy_vs_Counts' %(str(te),str(dust_perc)))
plt.title('Energy vs Counts(energy bin=0.01keV) Te=%skeV %sdust 11.0logtau\n' %(str(te),str(dust_perc)))
plt.xlabel('Energy (keV)')
plt.ylabel('Spectrum $({counts* \mathrm{cm}^3 \mathrm{s}^{-1}} \mathrm{keV}^{-1})$')
plt.bar(energybins,raw_spec,0.01, color='g',alpha=0.7,edgecolor='y')
plt.savefig('results/response/%skeV_Chandra_%sdust_Processed_input_data:Energy_vs_Counts.eps' %(str(te),str(dust_perc)))

###create the same size of spectrum as the rmf file, and set them all to zero
spec=np.zeros(1070)

###fill the spectrum bins with our data at the right place
for raw_spec_index,spec_index in enumerate(range(605,645)):
    spec[spec_index]=raw_spec[raw_spec_index]

####apply response
energy_grid,folded_spectrum=pyatomdb.spectrum.apply_response(spec,rmf,arf)

######calculate width of each bin. width is an array
width=[]
for i in range(435,462):
    width+=[energy_grid[i+1]-energy_grid[i]]

####plot reponse image
plt.figure('%skeV_Chandra_%sdust_after_responsing_by_west.rmf_west.arf' %(str(te),str(dust_perc)))
plt.title('Chandra model  Te=%skeV %sdust 11.0logtau\n' %(str(te),str(dust_perc)))
plt.bar(energy_grid[435:462],folded_spectrum[435:462],width,color='r',edgecolor='y')
plt.xlim(6.35,6.75)
plt.xlabel('Energy(keV)')
plt.ylabel('Spectrum $({counts* \mathrm{cm}^3 \mathrm{s}^{-1}} \mathrm{keV}^{-1})$')
plt.savefig('results/response/%skeV_Chandra_%sdust_after_responsing_by_west.rmf_west.arf.eps' %(str(te),str(dust_perc)))
plt.show()

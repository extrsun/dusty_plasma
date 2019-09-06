# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 16:19:05 2017

@author: zehaojin
"""



##########only allows band=2,tau=1e11.0,te=1/2/4/10

import pickle,pyatomdb, numpy as np, matplotlib.pyplot as plt
from constants import version

###################
rmf="response/ah_sxs_7ev_basefilt_20090216.rmf"
arf="response/sxt-s_100208_ts02um_flatsky-vig7kev_intallpxl.arf"
###################
'''
#################To see the energy bins in emf file
ret=pyatomdb.spectrum.get_response_ebins(rmf)
print ret
###############One of the results: energy bin size=0.001 keV
'''

dust_perc=input('enter the initial dust percentage, e.g., 95:\n')
print 'percentage of dust=',dust_perc

band=2
print 'band=',band

precision=0.5
print 'precision=',precision

te=input('enter temperature(currently 1.0 or 2.0 or 4.0 or 10.0):\n')
te=float(te)
print 'temperature %s keV' %(str(te))

print 'initializing...'



filename='spec_data_v%s_p%s_b%s_te%s_d%s_tau11.0' %(str(version),str(precision),str(band),str(te),str(dust_perc))
raw_spec=pickle.load(open('data_package/'+filename,'rb'))



#####create energy bins(0.001keV)and array of its mid value
energybins=np.linspace(6.35,6.74,400)
energybins_midvalue=energybins+0.0005

###########plot raw input data
plt.figure('%skeV_Astro-H_%sdust_Raw_input_data:Energy_vs_Emissivity' %(str(te),str(dust_perc)))
plt.title('Energy vs Emissivity(energy bin=0.001 keV) Te=%skeV %sdust 11.0logtau\n' %(str(te),str(dust_perc)))
plt.xlabel('Energy (keV)')
plt.ylabel('Emissivity $({\mathrm{erg} \mathrm{cm}^3 \mathrm{s}^{-1}} \mathrm{keV}^{-1})$')
plt.bar(energybins,raw_spec,0.001, color='b',alpha=0.7,edgecolor='b')
plt.savefig('results/response/v%s_%skeV_Astro-H_%sdust_Raw_input_data:Energy_vs_Emissivity.eps' %(str(version),str(te),str(dust_perc)))


#######cancel the erg in raw_spec
for raw_index in range(len(raw_spec)):
    raw_spec[raw_index]=raw_spec[raw_index]*6.242e8/energybins_midvalue[raw_index]
###plot
plt.figure('%skeV_Astro-H_%sdust_Processed_input_data:Energy_vs_Counts' %(str(te),str(dust_perc)))
plt.title('Energy vs Counts(energy bin=0.001keV) Te=%skeV %sdust 11.0logtau\n' %(str(te),str(dust_perc)))
plt.xlabel('Energy (keV)')
plt.ylabel('Spectrum $({counts* \mathrm{cm}^3 \mathrm{s}^{-1}} \mathrm{keV}^{-1})$')
plt.bar(energybins,raw_spec,0.001, color='g',alpha=0.7,edgecolor='g')
plt.savefig('results/response/v%s_%skeV_Astro-H_%sdust_Processed_input_data:Energy_vs_Counts.eps' %(str(version),str(te),str(dust_perc)))


###create the same size of spectrum as the rmf file, and set them all to zero
spec=np.zeros(16384)

###fill the spectrum bins with our data at the right place
for raw_spec_index,spec_index in enumerate(range(6350,6750)):
    spec[spec_index]=raw_spec[raw_spec_index]

####apply response
energy_grid,folded_spectrum=pyatomdb.spectrum.apply_response(spec,rmf,arf)

'''
Astro-H resulted energy bin is equally spaced, so it doesn't need this part
width=[]
for i in range(6350,6750):
    width+=[energy_grid[i+1]-energy_grid[i]]
'''
###plot response image
plt.figure('%skeV_Astro-H_%sdust_after_responsing_by_Astro-H_SXS' %(str(te),str(dust_perc)))
plt.title('Astro-H SXS model Te=%skeV %sdust 11.0logtau\n' %(str(te),str(dust_perc)))
plt.bar(energy_grid[6350:6750],folded_spectrum[6350:6750],0.001,color='r',edgecolor='r')
plt.xlim(6.35,6.75)
plt.xlabel('Energy(keV) [bin=0.001keV]')
plt.ylabel('Spectrum $({counts* \mathrm{cm}^3 \mathrm{s}^{-1}} \mathrm{keV}^{-1})$')
plt.savefig('results/response/v%s_%skeV_Astro-H_%sdust_after_responsing_by_Astro-H_SXS.eps' %(str(version),str(te),str(dust_perc)))
plt.show()

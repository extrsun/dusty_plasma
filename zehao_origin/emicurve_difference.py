# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 18:16:37 2017

@author: zehaojin
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np

print 'To modify default initial conditions, generate associate data package by pltemicurve_data.py, then Edit this code'

######Edit this variable to change which dust percentage to include. default 0 and 95
dust_perc_to_include=(0,95)
######Edit this variable to change energy band to include. default 2
band=2
######Edit this variable to change which two version to compare. default 3.0.8 vs 3.0.3
version_to_compare=('3.0.9','3.0.3')
######Edit this variable to change precision. default 0.2
precision=0.2
######temperature default 1,2,4,10; so index=4
te=(1.0,2.0,4.0,10.0)
temperature_index_counts=4
######usually leave it default
taulist=np.linspace(1,1e13,100001)


##print initial conditions
print 'dust percentage to include=',dust_perc_to_include
print 'band=',band
print 'version to compare=',version_to_compare
print 'precision=',precision
print ': gas   -- dust   - total'

##loop through temperaures
for dust_perc in dust_perc_to_include:
    
    ##get filename done for two versions
    version=version_to_compare[0]
    filename0='emicurve_data_v%s_p%s_b%s_d%s' %(str(version),str(precision),str(band),str(dust_perc))
    version=version_to_compare[1]
    filename1='emicurve_data_v%s_p%s_b%s_d%s' %(str(version),str(precision),str(band),str(dust_perc))
    
    ##load data
    data0=pickle.load(open('data_package/'+filename0,'rb'))
    data1=pickle.load(open('data_package/'+filename1,'rb'))

    ###create a figure~
    plt.figure('Diff_d%s_v%s_v%s' %(str(dust_perc),version_to_compare[0],version_to_compare[1]))
    plt.title('Difference in emicurve v%s-v%s %spercentage dust' %(version_to_compare[0],version_to_compare[1],str(dust_perc)))
    plt.xlabel('tau $({\mathrm{cm}^{-3}\mathrm{s}})$')
    plt.ylabel('percentage')
    ##loop through temperatures
    for temperature_index in range(temperature_index_counts):
        gasflux_difference=(data0['gasflux'][temperature_index,:]-data1['gasflux'][temperature_index,:])/data1['gasflux'][temperature_index,:]
        dustflux_difference=(data0['dustflux'][temperature_index,:]-data1['dustflux'][temperature_index,:])/data1['gasflux'][temperature_index,:]
        totalflux_difference=(data0['totalflux'][temperature_index,:]-data1['totalflux'][temperature_index,:])/data1['gasflux'][temperature_index,:]
        plt.plot(taulist,gasflux_difference,':')
        plt.plot(taulist,dustflux_difference,'--')
        plt.plot(taulist,totalflux_difference,label=str(te[temperature_index]))
        plt.xscale('log')
        #plt.yscale('symlog')
    ###legend 
    plt.legend(loc=0)


    ###x lim
    full=''
    plt.xlim(1e9,1e13)
    ###########
    #plt.xlim(0,1e13)
    #full='_full'

    ###save it as eps
    plt.savefig('results/emicurve_diff/Diff_d%s_v%s_v%s%s.eps' %(str(dust_perc),version_to_compare[0],version_to_compare[1],full))


###show it
plt.show()



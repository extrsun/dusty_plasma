# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 17:01:22 2017

@author: zehaojin
"""

import pickle,matplotlib.pyplot as plt
import constants

#####get input
dust_init=input('enter initial dust percentage, e.g., 95:\n')
print('percentage of dust=',dust_init)
###precision stuff
precision=input('precision(just press enter->default=0.2):\n')   
#eV, this number should be <0.5 in order to get high-reso spectra and right centroid (avoid taking conti- spec into account 
#    because criteria line < 0.01*maxflux is used)
                ###Zehao's note:0.2 for pltmicurve.py,pltgausscenter.py,plot_ionfrac.py
                ###             1.5 for pltspec.py dust percent=95
                ###             2.5 for pltspec.py dust percent=0
if precision=="":
    precision=0.2
precision=float(precision)
print('precision=',precision)

print('loading data...')

####load data package
dat=pickle.load(open('data_package/ionfrac_%s_v%s.pkl' %(str(dust_init),constants.version),'rb'))
####get data in the package
taulist=dat['taulist']
telist=dat['telist']
##telist=[1.0]#just for test
ionfrac=dat['ionfrac']

####do the same thing to 0 dust
dat0=pickle.load(open('data_package/ionfrac_0_v%s.pkl' %constants.version,'rb'))
#taulist0=dat0['taulist']#omit cuz it will be the same as taulist
#telist0=dat0['telist']#omit cuz it will be the same as telist
ionfrac0=dat0['ionfrac']


print('plotting......')
####loop through 4 temperatures
for ite,te in enumerate(telist):
    print('plotting... Temperature=%skeV' %te)
    ####create a plot for this temperature    
    plt.figure(ite)
    ####draw the dust's condition first####
    ###create subplot for dust's condition
    plt.subplot(121)
    ####loop through ions    
    for ion in range(0,28):
        ####build an empty list to store ionfrac data value
        ionfrac_value=[]
        ###get the data value by loop through tau
        for itau, tau in enumerate(taulist):
            ionfrac_value.append(ionfrac[ite,itau,ion])
        ###draw the plot for this ion
        plt.plot(taulist,ionfrac_value,label='Fe%s' %str(ion))
    ####set formats for the superposition of the 28 plots above
    plt.loglog()
    plt.xlim(1e9,1e13)
    plt.ylim(1e-6,1)
    plt.title('ionic fraction %s dust T=%skeV' %(str(float(dust_init)/100.0),str(te)))
    plt.xlabel('tau $({\mathrm{cm}^{-3}\mathrm{s}})$')
    plt.ylabel('ionic fraction')
    ####put the legend out of the graph
    plt.legend(loc='center left', bbox_to_anchor=(-0.3, 0.5))
    
    ####non-dust's condition below, similiar as the one above
    plt.subplot(122)
    for ion in range(0,28):
        ionfrac_value0=[]
        for itau, tau in enumerate(taulist):
            ionfrac_value0.append(ionfrac0[ite,itau,ion])
        plt.plot(taulist,ionfrac_value0,label='Fe%s' %str(ion))
    plt.loglog()
    plt.xlim(1e9,1e13)
    plt.ylim(1e-6,1)
    plt.title('ionic fraction 0 dust T=%skeV' %str(te))
    plt.xlabel('tau $({\mathrm{cm}^{-3}\mathrm{s}})$')
    plt.ylabel('ionic fraction')
    plt.legend(loc='center right', bbox_to_anchor=(1.2, 0.5))
#    manager = plt.get_current_fig_manager()
#    manager.resize(*manager.window.maxsize())

###okay cool
print('opening graphs...\nemm.... save them as .eps yourself if needed')
plt.show()

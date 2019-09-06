# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 16:46:47 2017

@author: zehaojin
"""

import numpy,pickle,dust
from constants import *


dust_init=input('enter initial dust percentage, e.g., 95:\n')
print 'percentage of dust=',dust_init


band=raw_input('what band? 1. 6.35--7.15(kalpha including lyman); 2. 6.35--6.75(excluding lyman)(press enter->default=2):\n')
if band=="":
	band=2
band=int(band)
print 'band=',band

precision=raw_input('precision(just press enter->default=0.2):\n')   
#eV, this number should be <0.5 in order to get high-reso spectra and right centroid (avoid taking conti- spec into account 
#    because criteria line < 0.01*maxflux is used)
                ###Zehao's note:0.2 for pltmicurve.py,pltgausscenter.py,plot_ionfrac.py
                ###             1.5 for pltspec.py dust percent=95
                ###             2.5 for pltspec.py dust percent=0
if precision=="":
    precision=0.2
precision=float(precision)
print 'precision=',precision
print 'start calculating...'

dat=pickle.load(open('data_package/ionfrac_%s_v%s.pkl' %(str(dust_init),version),'rb'))
taulist=dat['taulist']

if band==1:
	dat2=pickle.load(open('data_package/eperionperbin_6.35-7.15_v%s_p%s.pkl' %(version,precision),'rb'))
if band==2:
	dat2=pickle.load(open('data_package/eperionperbin_6.35-6.75_v%s_p%s.pkl' %(version,precision),'rb'))




gasflux=numpy.zeros([len(telist),len(taulist)])
dustflux=numpy.zeros([len(telist),len(taulist)])
totalflux=numpy.zeros([len(telist),len(taulist)])


eperion=sum(dat2['eperionperbin'],axis=2) #sum up energy bins to get emissivity per ion
print 'summed emissivity per ion=',eperion

for ite,te in enumerate(telist):
	print 'calculating Temperature=',te
	for itau, tau in enumerate(taulist):
		etmp=eperion[ite,:]*dat['ionfrac'][ite,itau,:]
		gasflux[ite,itau]=sum(etmp)
		dustflux[ite,itau]=dat['ionfrac'][ite,itau,0]*dust.calc_dust_kalpha_emi(te)
		totalflux[ite,itau]=gasflux[ite,itau]+dustflux[ite,itau]
  
print 'saving data...'
###build a dictionary
data={}
data['gasflux']=gasflux
data['dustflux']=dustflux
data['totalflux']=totalflux
###save the data
filename='emicurve_data_v%s_p%s_b%s_d%s' %(str(version),str(precision),str(band),str(dust_init))
pickle.dump(data,open('data_package/'+filename,'wb'))



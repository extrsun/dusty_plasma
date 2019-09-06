import pyatomdb, numpy, pickle, dust, matplotlib.pyplot as plt
from constants import *


dust_init=input('enter initial dust percentage, e.g., 95:\n')
print('percentage of dust=',dust_init)


band=input('what band? 1. 6.35--7.15(kalpha including lyman); 2. 6.35--6.75(excluding lyman)(press enter->default=2):\n')
if band=="":
	band=2
band=int(band)
print('band=',band)

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
print('start plotting...')

dat=pickle.load(open('data_package/ionfrac_%s_v%s.pkl' %(str(dust_init),version),'rb'))
taulist=dat['taulist']

if band==1:
	dat2=pickle.load(open('data_package/eperionperbin_6.35-7.15_v%s_p%s.pkl' %(version,precision),'rb'))
if band==2:
    # dat2=pickle.load(open('data_package/eperionperbin_6.35-6.75_v%s_p%s.pkl' %(version,precision),'rb'))
    dat2=pickle.load(open('data_package/eperionperbin_7.0-7.9_v%s_p%s.pkl' %(version,precision),'rb'))





gasflux=numpy.zeros([len(telist),len(taulist)])
dustflux=numpy.zeros([len(telist),len(taulist)])
totalflux=numpy.zeros([len(telist),len(taulist)])


eperion=sum(dat2['eperionperbin'],axis=2) #sum up energy bins to get emissivity per ion
print('summed emissivity per ion=',eperion)

for ite,te in enumerate(telist):
	print('plotting Temperature=',te)
	for itau, tau in enumerate(taulist):
		etmp=eperion[ite,:]*dat['ionfrac'][ite,itau,:]
		gasflux[ite,itau]=sum(etmp)
		dustflux[ite,itau]=dat['ionfrac'][ite,itau,0]*dust.calc_dust_kalpha_emi(te)
		totalflux[ite,itau]=gasflux[ite,itau]+dustflux[ite,itau]


for ite,te in enumerate(telist):
	plt.loglog(taulist,gasflux[ite,:],':')
	plt.loglog(taulist,dustflux[ite,:],'--')
	plt.loglog(taulist,totalflux[ite,:],label=repr(te))
plt.title(str(dust_init)+'% of dust')
plt.xlabel('tau $({\mathrm{cm}^{-3}\mathrm{s}})$')
plt.ylabel('emissivity $({\mathrm{erg} \mathrm{cm}^3 \mathrm{s}^{-1}})$')
########
plt.ylim(1e-28,1e-23)

#######
full=''
plt.xlim(1e9,1e13)
#########
#full='_full'
plt.legend(loc=0,fancybox=True,framealpha=0.5)
plt.savefig('results/emicurve/pltemicurve_%sdust_v%s%s.eps' %(str(dust_init),version,full))



plt.show()

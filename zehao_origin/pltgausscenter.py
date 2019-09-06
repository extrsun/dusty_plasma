import pyatomdb, numpy, pickle,pylab,dust
import matplotlib.pyplot as plt
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
print('initializing...')

if band==1:
	dat2=pickle.load(open('data_package/eperionperbin_6.35-7.15_v%s_p%s.pkl' %(version,precision),'rb'))
if band==2:
    dat2=pickle.load(open('data_package/eperionperbin_6.35-6.75_v%s_p%s.pkl' %(version,precision),'rb'))
    # dat2=pickle.load(open('data_package/eperionperbin_7.0-7.9_v%s_p%s.pkl' %(version,precision),'rb'))
    
    


dat=pickle.load(open('data_package/ionfrac_%s_v%s.pkl' %(str(dust_init),version),'rb'))
dat0=pickle.load(open('data_package/ionfrac_0_v%s.pkl' %version,'rb'))


telist=dat['telist']

taulist=dat['taulist']

print('full tau list=\n',taulist)

tauind=numpy.logspace(1,5,100)    #there are 10^5 tau data points in ionfrac.pkl,choose 100, log spaced, just to reduce calculation time 
tauind=numpy.array(list([1,2,3,4,5,6,7,8,9,10])+list(tauind),dtype=int)
print('selected tau index to calculate=\n',tauind)

shorttaulist=taulist[tauind]
print('selected tau to calculate=\n',shorttaulist,'\ntotal count of selected tau=',len(shorttaulist))
print('selected temperature(keV) to calculate=',telist)

input('press enter to start plotting...')


dust_dblline_center=E_keV 
dust_dblline_yield=yie  
dust_center=(E_keV[1]*yie[1]+E_keV[2]*yie[2])/(yie[1]+yie[2])


fig=pylab.figure(1)

#kcent stores dust model result, kcent0 stores pure-gas model result
kcent=numpy.zeros([len(telist),len(tauind)])
kcent0=numpy.zeros([len(telist),len(tauind)])
diff=numpy.zeros([len(telist),len(tauind)])

for ite,te in enumerate(telist):
	for iiindex,index in enumerate(tauind):

			print('calculating...  temperature list=',telist)
			print('temperature=',te,'keV','index of tau=',index)
			#numerator stores dust model result, numer pure-gas.
			numerator=0.
			denominator=0.
			numer=0.
			denom=0.
			numerator_array=numpy.zeros(len(dat2['ecent']),dtype=float)
			denominator_array=numpy.zeros(len(dat2['ecent']),dtype=float)
			numer_array=numpy.zeros(len(dat2['ecent']),dtype=float)
			denom_array=numpy.zeros(len(dat2['ecent']),dtype=float)
			for iebin,ebin in enumerate(dat2['ecent']):
				#dust model
				numerator_array[iebin]=ebin*sum(dat2['eperionperbin'][ite,:,iebin]*dat['ionfrac'][ite,index,:])
				denominator_array[iebin]=sum(dat2['eperionperbin'][ite,:,iebin]*dat['ionfrac'][ite,index,:])
				#pure-gas model
				numer_array[iebin]=ebin*sum(dat2['eperionperbin'][ite,:,iebin]*dat0['ionfrac'][ite,index,:])
				denom_array[iebin]=sum(dat2['eperionperbin'][ite,:,iebin]*dat0['ionfrac'][ite,index,:])
			for iebin,ebin in enumerate(dat2['ecent']):
				#choose strong lines to calc centroid
				if denominator_array[iebin]> (0.01*max(denominator_array)):
					numerator+=numerator_array[iebin]
					denominator+=denominator_array[iebin]
				if denom_array[iebin]>(0.01*max(denom_array)):
					numer+=numer_array[iebin]
					denom+=denom_array[iebin]
			#dust
			dust_emi=dust.calc_dust_kalpha_emi(te)
			dust_flux=dat['ionfrac'][ite,index,0]*dust_emi
			numerator+=(dust_flux*dust_center)
			denominator+=dust_flux


			kcent[ite,iiindex]=numerator/denominator
			kcent0[ite,iiindex]=numer/denom
			diff[ite,iiindex]=kcent0[ite,iiindex]-kcent[ite,iiindex]

ax1=fig.add_subplot(111)

print('plotting...\n')#,len(shorttaulist),len(diff)
for ite,te in enumerate(telist):
	plt.plot(shorttaulist,diff[ite,:]*1000,label=repr(te))
	plt.semilogx()
#ax.set_ylim([-20,200])
plt.xlabel('tau $({\mathrm{cm}^{-3}\mathrm{s}})$')
plt.ylabel('centroid shift(eV)')
plt.xlim(1e9,1e13)
plt.ylim(0,800)

ax1.legend(loc=0)
plt.savefig('results/Gauss_center/centeroid_%s.eps' %(version))

fig2=pylab.figure(2)
ax=fig2.add_subplot(121)
for ite,te in enumerate(telist):
    plt.plot(shorttaulist,kcent0[ite,:]*1000,label=repr(te))
    plt.semilogx()
print('kcent with 0 dust=\n',kcent0[0,:])
# print('kcent with %s dust=\n' %dust_init/100.0,kcent[0,:])
print('difference=\n',diff[0,:])
plt.xlabel('tau $({\mathrm{cm}^{-3}\mathrm{s}})$')
plt.ylabel('centroid (eV)')
plt.xlim(1e9,1e13)
plt.ylim(6900,7950)
plt.legend(loc=0)
ax=fig2.add_subplot(122)
for ite,te in enumerate(telist):
    plt.plot(shorttaulist,kcent[ite,:]*1000,label=repr(te))
    plt.semilogx()

plt.xlabel('tau $({\mathrm{cm}^{-3}\mathrm{s}})$')
plt.ylabel('centroid (eV)')
plt.xlim(1e9,1e13)
plt.ylim(6900,7950)
plt.legend(loc=0)
plt.tight_layout()
plt.savefig('results/Gauss_center/centeroid_shift_%s.eps' %(version),bbox_inches='tight')
plt.show()

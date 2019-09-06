from constants import *
import numpy,pickle,dust
import matplotlib.pyplot as plt

dust_perc=input('enter the initial dust percentage, e.g., 95:\n')
print 'percentage of dust=',dust_perc

band=raw_input('what band? 1. 6.35--7.15(kalpha including lyman); 2. 6.35--6.75(excluding lyman)(press enter->default=2):\n')
if band=="":
	band=2
band=int(band)
print 'band=',band,'\n\n'

print "Note: To fit Huanqin's graph, Zehao suggest:\nprecision=1.5 for dust percent=95;\nprecision=2.5 for dust percent=0\nBoth stop at tau=10^11\n"

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

print 'initializing...'


if band==1:
	dat3=pickle.load(open('data_package/eperionperbin_6.35-7.15_v%s_p%s.pkl' %(version,precision),'rb'))
if band==2:
	dat3=pickle.load(open('data_package/eperionperbin_6.35-6.75_v%s_p%s.pkl' %(version,precision),'rb'))
 
dat=pickle.load(open('data_package/ionfrac_%s_v%s.pkl' %(str(dust_perc),version),'rb'))


telist=dat['telist']
taulist=dat['taulist']

width=dat3['binwidth']
ecent=dat3['ecent']

tauind=[5,10,100,300,1000,2000,4000,10000,100000]
print "energy bin's width=",width,'keV'
print 'list of centers of each energy bin=\n',ecent 

raw_input('press enter to start drawing')

eperionperbin=dat3['eperionperbin']

dust_dblline_center=E_keV #
dust_dblline_yield=yie    #
print 'dust_dblline_center=',dust_dblline_center[1]


print 'energy:\n',ecent-width/2.
for ite,te in enumerate(telist):
	for index in tauind:	
		gas_spec=numpy.zeros(len(ecent))	
		dust_spec=numpy.zeros(len(ecent))
		for ibin in range(len(ecent)):
			gas_spec[ibin]=sum(eperionperbin[ite,:,ibin]*dat['ionfrac'][ite,index,:])/width
		dust_total_emi=dust.calc_dust_kalpha_emi(te)
		print '\nionfrac rate=',dat['ionfrac'][ite,index,0],'  @Temperature=%s, tau=%s' %(str(te),str(taulist[index]))
		dust_one=dat['ionfrac'][ite,index,0]*dust_total_emi*E_keV[1]*yie[1]/(E_keV[1]*yie[1]+E_keV[2]*yie[2])
		dust_two=dat['ionfrac'][ite,index,0]*dust_total_emi*E_keV[2]*yie[2]/(E_keV[1]*yie[1]+E_keV[2]*yie[2])
		print 'dust_total_emi=',dust_total_emi
		print 'dust_one=',dust_one
		for ii in range(len(ecent)):
			if abs(ecent[ii]-dust_dblline_center[1]) < (width/2.):
				dust_spec[ii]+=dust_one/width
				#print 'find it!',ii
				break
			#else:
				#print 'oh noooooo'
		for ii in range(len(ecent)):
			if abs(ecent[ii]-dust_dblline_center[2]) <= (width/2.):
				dust_spec[ii]+=dust_two/width
				#print '2nd one!',ii
				break
			#else:	
				#print 'what?? where is the second one!'
		plt.figure(ite)
		plt.xlabel('energy (keV)')
		plt.ylabel('emissivity $({\mathrm{erg} \mathrm{cm}^3 \mathrm{s}^{-1}} \mathrm{keV}^{-1})$')
		plt.bar(ecent-width/2., gas_spec, width, color='b',alpha=0.7,label='gas_spec',edgecolor='b')
		plt.bar(ecent-width/2., dust_spec,width, color='r',bottom=gas_spec,alpha=0.7,label='dust_spec',edgecolor='r')
		plt.legend(loc=0)
		
		plt.xlim(6.35,6.75)
		#if dust_perc==0:
			#plt.xlim(6.35,6.75)
			#plt.ylim(0,6*1e-25)
		#if dust_perc==95:
			#plt.xlim(6.3,6.8)
			#plt.ylim(0,1e-25)
		
		plt.title('tau:10^'+str(round(log10(taulist[index]),1))+' Te:'+str(te)+'keV '+str(dust_perc)+'% dust')
		plt.savefig('results/spec/'+str(te)+'keV_'+str(round(log10(taulist[index]),1))+'logtau_p%s_v%s_' %(str(precision),str(version))+str(dust_perc)+'dust.eps')
		print "energy bin's width=",width,'keV'
		print 'close the graph to continue'
  		plt.show()
		next=raw_input('PRESS ENTER(to draw next tau) \nOR type anything then press enter(to skip to next temperature)\n')
		if next!='':
			break




'''
dat=pickle.load(open('ionfrac_dust.pkl','rb'))


for ite,te in enumerate(telist):
	for itau, tau in enumerate(taulist):
		etmp=dat3['eperionperbin'][ite,:]*dat['ionfrac'][ite,itau,:]
		gasflux[ite,itau]=sum(etmp)
		dustflux[ite,itau]=dat['ionfrac'][ite,itau,0]*dust.calc_dust_kalpha_emi(te)
		totalflux[ite,itau]=gasflux[ite,itau]+dustflux[ite,itau]

ax=fig.add_subplot(122)

for ite,te in enumerate(telist):
	ax.loglog(taulist,gasflux[ite,:],':',label=repr(te))
	ax.loglog(taulist,dustflux[ite,:],'--')
	ax.loglog(taulist,totalflux[ite,:],)

ax.legend(loc=0)
'''
#plt.show()


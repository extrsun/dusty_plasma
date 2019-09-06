#using function pyatomdb.spectrum.get_index to get emissivity for each ion each energy bin
#Divide Kalpha band into many bins. For each ion, calculate its emissivity in each bin. Write data into eperionperbin_wave_range.pkl
#result(out) stored in eperionperbin_energy_range.pkl
#out is a dictionary: out['ecent']: middle energy of every bin, out['binwidth']: width of each bin
#				out['eperionperbin']: 3-D list emissivity for each ion each bin [Te,ion, bins]


from constants import *

band=input('what band? 1. 6.35--7.15(kalpha including lyman); 2. 6.35--6.75(excluding lyman)(press enter->default=2):\n')
if band=="":
	band=2
band=int(band)
print('band=',band)

#######################
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
print('working on it...')
########################

if band==1:
	number=int(400/precision)
if band==2:
	number=int(200/precision)

#get bin grids
if band==1:
	ebins=numpy.linspace(6.35,7.15,number*2+1)
	filename='eperionperbin_6.35-7.15_v%s_p%s.pkl' %(version,precision)
if band==2:
    ebins=numpy.linspace(6.35,6.75,number+1)
    filename='eperionperbin_6.35-6.75_v%s_p%s.pkl' %(version,precision)
    # ebins=numpy.linspace(7.05,7.95,number+1)
    # filename='eperionperbin_7.0-7.9_v%s_p%s.pkl' %(version,precision)


ecent=(ebins[1:]+ebins[:-1])/2
emiss={}

linefile = '$ATOMDB/apec_nei_line.fits'
cocofile = '$ATOMDB/apec_nei_comp.fits'


eperionperbin=numpy.zeros([len(telist),28,len(ecent)],dtype=float)

for ite,te in enumerate(telist):
	ihdu=pyatomdb.spectrum.get_index(te,linefile)
	for z1 in range(1,28):
		spec=pyatomdb.spectrum.make_ion_spectrum(ebins,\
												ihdu,26,z1,\
												linefile=linefile,\
												cocofile=cocofile)
		spec=spec*ecent*pyatomdb.const.ERG_KEV
		eperionperbin[ite,z1,:]=spec

out={}
out['eperionperbin']=eperionperbin
out['ecent']=ecent
out['binwidth']=ecent[1]-ecent[0]
print('eperionperbin=\n',eperionperbin)
print('ebins=\n',ebins)
print('binwidth=',out['binwidth'])
pickle.dump(out,open('data_package/'+filename,'wb'))


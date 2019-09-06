import pyatomdb, numpy, pickle,pylab,dust
import matplotlib.pyplot as plt
from constants import *


dat1=pickle.load(open('data_package/ionfrac_95_v3.0.8.pkl','rb'))
taulist=dat1['taulist']
telist=dat1['telist']
dat2=pickle.load(open('data_package/ionfrac_0_v3.0.8.pkl','rb'))
y1=[]
#tau=linspace(1,1e13,1001)

plt.figure(1)

plt.subplot(122)
for z in range(0,28):

	for itau, tau in enumerate(taulist):
		a=dat1['ionfrac'][3,itau,:]

	for k in a:
		y1.append(k[z])
	plt.plot(taulist,y1)
	'''
	if z<18:
		plt.plot(taulist,y)
	else:
		plt.plot(taulist,y,label='Fe%s'%z)
		pass
	'''
plt.ylim(1e-6,1e0)
plt.xlim(1e9,1e13)
plt.loglog()
plt.title('ionic fraction 0.95 dust T=%s kev'%telist[3])
plt.xlabel('tau $({\mathrm{cm}^{-3}\mathrm{s}})$')
plt.ylabel('ionic fraction')
#plt.legend(loc=0)

'''
plt.subplot(121)
for z in range(0,28):
	h=[]
	for itau, tau in enumerate(taulist):
		m=dat2['ionfrac'][3,itau,:]
		h.append(m)
	plt.plot(taulist,h)
plt.ylim(1e-6,1e0)
plt.xlim(1e9,1e13)
plt.loglog()
plt.title('ionic fraction 0 dust T=%s kev'%telist[3])
plt.xlabel('tau $({\mathrm{cm}^{-3}\mathrm{s}})$')
plt.ylabel('ionic fraction')
'''
#plt.legend(loc=0)

plt.show()
##############
'''
import pyatomdb, numpy, pickle,pylab,dust
import matplotlib.pyplot as plt
from constants import *


dat1=pickle.load(open('ionfrac_95.pkl','rb'))
taulist=dat1['taulist']
telist=dat1['telist']
dat2=pickle.load(open('ionfrac_0.pkl','rb'))

#tau=linspace(1,1e13,1001)

plt.figure(1)

plt.subplot(122)
for z in range(0,28):
	a=[]
	y=[]
	#for itau in range(len(taulist)):
	for itau, tau in enumerate(taulist):
		x=dat1['ionfrac'][3,itau,:]
		y.append(x)
	for k in y:
		a.append(k[z])
	plt.plot(taulist,a)
plt.ylim(1e-6,1e0)
plt.xlim(1e9,1e13)
plt.loglog()
plt.title('ionic fraction 0.95 dust T=%s kev'%telist[3])
plt.xlabel('tau $({\mathrm{cm}^{-3}\mathrm{s}})$')
plt.ylabel('ionic fraction')
#plt.legend(loc=0)

#plt.title('ionic fraction 0.95 dust T=',telist[3])
'''
'''
plt.subplot(121)
for z in range(0,28):
	b=[]
	h=[]
	#for itau in range(len(taulist)):
	for itau, tau in enumerate(taulist):
		m=dat2['ionfrac'][3,itau,:]
		h.append(m)
	for k in h:
		b.append(k[z])
	plt.plot(taulist,b)
plt.ylim(1e-6,1e0)
plt.xlim(1e9,1e13)
plt.loglog()
plt.title('ionic fraction 0 dust T=%s kev'%telist[3])
plt.legend(loc=0)
'''
'''
plt.show()

TK=1e7
tau=linspace(1,1e13,1000001)
ele=zeros(28,dtype=float)
ele[0]=0
ele[1]=1.0
y=calc_nei(26,tau,ele,TK)
fig1=plt.figure('fig1')
for z in range(0,28):
	a=[]
	for k in y:
		a.append(k[z])
	plt.plot(tau,a)
plt.loglog()
plt.show()
'''
############
'''
import pyatomdb, numpy, pickle,pylab,dust
import matplotlib.pyplot as plt
from constants import *


dat1=pickle.load(open('ionfrac_95.pkl','rb'))
taulist=dat1['taulist']
telist=dat1['telist']
dat2=pickle.load(open('ionfrac_0.pkl','rb'))

#tau=linspace(1,1e13,1001)

plt.figure(1)

plt.subplot(122)
for z in range(0,28):
	a=[]
	y=[]
	for itau, tau in enumerate(taulist):
		x=dat1['ionfrac'][3,itau,:]
		y.append(x)
	for k in y:
		a.append(k[z])
	plt.plot(taulist,a)
plt.ylim(1e-6,1e0)
plt.xlim(1e9,1e13)
plt.loglog()
plt.title('ionic fraction 0.95 dust T=%s kev'%telist[3])
plt.xlabel('tau $({\mathrm{cm}^{-3}\mathrm{s}})$')
plt.ylabel('ionic fraction')
#plt.legend(loc=0)

plt.subplot(121)
for z in range(0,28):
	b=[]
	h=[]
	#for itau in range(len(taulist)):
	for itau, tau in enumerate(taulist):
		m=dat2['ionfrac'][3,itau,:]
		h.append(m)
	for k in h:
		b.append(k[z])
	if z<18:
		plt.plot(taulist,y)
	else:
		plt.plot(taulist,y,label='Fe%s'%z)
		pass
	#plt.plot(taulist,b)
plt.ylim(1e-6,1e0)
plt.xlim(1e9,1e13)
plt.loglog()
plt.title('ionic fraction 0 dust T=%s kev'%telist[3])
plt.legend(loc=0)

plt.show()
'''

'''
TK=1e7
tau=linspace(1,1e13,1000001)
ele=zeros(28,dtype=float)
ele[0]=0
ele[1]=1.0
y=calc_nei(26,tau,ele,TK)
fig1=plt.figure('fig1')
for z in range(0,28):
	a=[]
	for k in y:
		a.append(k[z])
	plt.plot(tau,a)
plt.loglog()
plt.show()
'''

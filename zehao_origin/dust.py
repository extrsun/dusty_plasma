from scipy import *
from numpy import *
import matplotlib.pyplot as plt

ryd=0.0136 
#keV；1 ryd = 0.0136 keV
Ip=7.112 
#keV，I为铁的 K 层电离势阱 7.112 keV。
keV2erg=1.602e-9
abund=4.67735e-05
#？

E_keV=[0,6.4077,6.3915]
#For isolated Fe atoms, the ﬂuorescent Fe Ka line consists of two components,K α 1 =6.404keV 和 K α 2 =6.391keV
yie=[0,0.2026,0.1013] 
#florescent yields;K-shell ionization rates must be multiplied by the ﬂuorescent yields for production of Ka photons (equal to 0.30 for Fe)

def calc_dust_kalpha_emi(Te_keV):
	
	beta=Ip/Te_keV
	S=keV2erg*1.13e-7*(ryd/Ip)**1.5*beta**0.5/exp(beta)/(beta+0.4)
	#此为k-shell collisional ionization rate coefficient by electron with a maxwell velocity distribution.
	emiKalpha=S*(E_keV[1]*yie[1]+E_keV[2]*yie[2])*abund
	#why *abund
	return emiKalpha
'''
这个函数后面的plot gausscenter没有用到，直接写的下面的公式
def calc_dust_kalpha_center():
	center=(E_keV[1]*yie[1]+E_keV[2]*yie[2])/(yie[1]+yie[2])
	return center
'''	
'''
plt.plot(tau_13,dustKalpha)
plt.loglog()
plt.ylim([1e-30,1e-22])
plt.show()
print dustKalpha
'''

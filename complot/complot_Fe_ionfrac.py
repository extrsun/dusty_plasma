"""
Purpose: make a comparison plot: Zehao's code v.s. my approach

In this file, create the ionic fraction plot

Written by sw, Jun 20, 2019
"""

# Import used modules
import pyatomdb
import pickle, os
import numpy as np
import astropy.io.fits as pyfits
import matplotlib as mpl
import matplotlib.pyplot as plt

def complot_Fe_ionfrac_module(log=False):
  #system parameters
  rootpath = os.getcwd()+'/'
  kTlist  = np.array([1.,2.,4.,10.]) # in keV
  n_kT = len(kTlist)

  Z = 26
  elsymb = pyatomdb.atomic.Ztoelsymb(Z)

  if log:
    subdirec = 'log'
  else:
    subdirec = 'linear'
  modname = ['Zehao', 'Mine']
  ionfracfile_name = ['zehao_mimic/data_package/ionfrac_95_v3.0.9.pkl', \
    'my_approach/'+subdirec+'/dustsnr_ionfrac_Fe.pkl']
  nmod = 2

  ionfrac1 = pickle.load(open(rootpath + ionfracfile_name[0], 'rb'))
  ionfrac2 = pickle.load(open(rootpath + ionfracfile_name[1], 'rb'))

  begin_ion = 17
  ion_range = range(begin_ion+1,27)

  for M in range(0,n_kT):
    plt.plot(ionfrac1['taulist'], ionfrac1['ionfrac'][M,:,0], \
      linewidth=1.0, label='Neutral', color='C0')
    for iZ in ion_range:
      plt.plot(ionfrac1['taulist'], ionfrac1['ionfrac'][M,:,iZ+1], \
        color='C'+"%s"%(iZ-begin_ion), linewidth=1.0, \
        label=elsymb+' '+pyatomdb.atomic.int2roman(iZ+1))
    plt.plot(ionfrac2['taulist'], ionfrac2['ionfrac'][M,0,:], \
      linewidth=2.0, linestyle='dashed', color='C0')
    for iZ in ion_range:
      plt.plot(ionfrac2['taulist'], ionfrac2['ionfrac'][M,iZ+1,:], \
        color='C'+"%s"%(iZ-begin_ion), linewidth=2.0, linestyle='dashed')
    plt.plot(ionfrac1['taulist'], np.zeros(len(ionfrac1['taulist'])), \
      linewidth=1.0, label=modname[0])
    plt.plot(ionfrac2['taulist'], np.zeros(len(ionfrac2['taulist'])), \
      linewidth=2.0, linestyle='dashed',label=modname[1])
    plt.xlabel('Ionization Timescale (cm$^{-3}$ s)')
    plt.ylabel('Ionic Fraction')
    plt.title("Fe for dust of 95 percentage, with $kT$=%04.1f keV" % kTlist[M])
    plt.xlim(1e10,1e13)
    plt.ylim(1e-6,1.0)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc=0)
    plt.savefig(rootpath+'figures/'+subdirec+'/comp_ionfrac_%04.1f.eps' % kTlist[M])
    plt.close()
   
    n_tau   = len(ionfrac2['taulist'])
    pltfrac = np.zeros([10, n_tau], dtype=float)
    pltfrac[0,:] = np.sum(ionfrac2['ionfrac'][M,1:20,:], axis=0)
    for iZ in range(20,28):
      pltfrac[iZ-19,:] = ionfrac2['ionfrac'][M,iZ,:]
    pltfrac[9,:] = ionfrac2['ionfrac'][M,0,:]
    
    labels = ["<Fe XX", "Fe XX", "Fe XXI", "Fe XXII", "Fe XXIII", \
              "Fe XXIV", "Fe XXV", "Fe XXVI", "Fe XXVII", 'Neutral']

    plt.stackplot(ionfrac2['taulist'], pltfrac[0,:], pltfrac[1,:], \
         pltfrac[2,:], pltfrac[3,:], pltfrac[4,:], pltfrac[5,:], \
         pltfrac[6,:], pltfrac[7,:], pltfrac[8,:], pltfrac[9,:], \
         labels=labels)
    plt.legend(loc='upper left')
    plt.xlim(1e10,1e13)
    plt.ylim(0.0,1.0)
    plt.xscale('log')
    plt.title("Fe for dust of 95 percentage, with $kT$=%04.1f keV" % kTlist[M])
    plt.xlabel('Ionization Timescale (cm$^{-3}$ s)')
    plt.ylabel('Ionic Fraction')
    plt.savefig(rootpath+'figures/'+subdirec+'/comp_accumfrac_%04.1f.eps' % kTlist[M])
    plt.close()


def main():
    complot_Fe_ionfrac_module(log=True)
    complot_Fe_ionfrac_module(log=False)

if __name__ == '__main__':
    main()

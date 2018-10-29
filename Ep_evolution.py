from __future__ import division
#~ import os
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import uncertainties as unc  
import uncertainties.unumpy as unumpy  
from uncertainties import ufloat

###########################################################################################################################################################################################################################

font = {
        'weight' : 'bold',
        'size'   : 15}

mpl.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 1.25
mpl.rcParams['axes.labelweight']='bold'
mpl.rcParams['axes.labelsize']='large'
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['ytick.major.width'] = 2

lstyle = 'None'

#Section to take the input parameters

#Band
band_ifl = ascii.read("results_band.txt", header_start = None, comment = '#')
band_alpha   = np.array(band_ifl['col4'])
band_alpha_n = np.array(band_ifl['col5'])
band_alpha_p = np.array(band_ifl['col6'])
band_beta    = np.array(band_ifl['col7'])    
band_beta_n  = np.array(band_ifl['col8'])
band_beta_p  = np.array(band_ifl['col9'])
band_Ep = np.array(band_ifl['col10'])
band_Ep_n = np.array(band_ifl['col11'])
band_Ep_p = np.array(band_ifl['col12'])


#BandC
bandC_ifl = ascii.read("results_bandC_freeze.txt", header_start = None, comment = '#')
bandC_alpha   = np.array(bandC_ifl['col4'])
bandC_alpha_n = np.array(bandC_ifl['col5'])
bandC_alpha_p = np.array(bandC_ifl['col6'])
bandC_beta    = np.array(bandC_ifl['col7'])    
bandC_beta_n  = np.array(bandC_ifl['col8'])
bandC_beta_p  = np.array(bandC_ifl['col9'])
bandC_Ep = np.array(bandC_ifl['col10'])
bandC_Ep_n = np.array(bandC_ifl['col11'])
bandC_Ep_p = np.array(bandC_ifl['col12'])
bandC_cut = np.array(bandC_ifl['col22'])
bandC_cut = np.array(bandC_ifl['col23'])
bandC_cut = np.array(bandC_ifl['col24'])

#BB+BandC
bb_bandC_ifl = ascii.read("bb_bandC_freeze_results.txt", header_start = None, comment = '#')
bb_bandC_kT =  np.array(bb_bandC_ifl['col4'])
bb_bandC_kT_n =  np.array(bb_bandC_ifl['col5'])
bb_bandC_kT_p =  np.array(bb_bandC_ifl['col6'])
bb_bandC_alpha   = np.array(bb_bandC_ifl['col10'])
bb_bandC_alpha_n = np.array(bb_bandC_ifl['col11'])
bb_bandC_alpha_p = np.array(bb_bandC_ifl['col12'])
bb_bandC_beta    = np.array(bb_bandC_ifl['col13'])    
bb_bandC_beta_n  = np.array(bb_bandC_ifl['col14'])
bb_bandC_beta_p  = np.array(bb_bandC_ifl['col15'])
bb_bandC_Ep = np.array(bb_bandC_ifl['col16'])
bb_bandC_Ep_n = np.array(bb_bandC_ifl['col17'])
bb_bandC_Ep_p = np.array(bb_bandC_ifl['col18'])
bb_bandC_cut = np.array(bb_bandC_ifl['col28'])
bb_bandC_cut = np.array(bb_bandC_ifl['col29'])
bb_bandC_cut = np.array(bb_bandC_ifl['col30'])

#BB+Band
bb_band_ifl = ascii.read("results_bb_band.txt", header_start = None, comment = '#')
bb_band_alpha   = np.array(bb_band_ifl['col4'])
bb_band_alpha_n = np.array(bb_band_ifl['col5'])
bb_band_alpha_p = np.array(bb_band_ifl['col6'])
bb_band_beta    = np.array(bb_band_ifl['col7'])    
bb_band_beta_n  = np.array(bb_band_ifl['col8'])
bb_band_beta_p  = np.array(bb_band_ifl['col9'])
bb_band_Ep = np.array(bb_band_ifl['col10'])
bb_band_Ep_n = np.array(bb_band_ifl['col11'])
bb_band_Ep_p = np.array(bb_band_ifl['col12'])
bb_band_kT =  np.array(bb_band_ifl['col19'])
bb_band_kT_n =  np.array(bb_band_ifl['col20'])
bb_band_kT_p =  np.array(bb_band_ifl['col21'])

#Cutoff Powerlaw
cutoffpl_ifl = ascii.read("pre_results_cutoffpl.txt", header_start = None, comment = '#')
cutoffpl_alpha   = np.array(cutoffpl_ifl['col4'])
cutoffpl_alpha_n = np.array(cutoffpl_ifl['col5'])
cutoffpl_alpha_p = np.array(cutoffpl_ifl['col6'])
cutoffpl_Ep = np.array(cutoffpl_ifl['col7'])		#PLot This
cutoffpl_Ep_n = np.array(cutoffpl_ifl['col8'])
cutoffpl_Ep_p = np.array(cutoffpl_ifl['col9'])

#BB+PL
bbpl_pre_ifl = ascii.read("pre_results_bbpl.txt", header_start = None, comment = '#')
bbpl_pre_kT   = np.array(bbpl_pre_ifl['col4'])
bbpl_pre_kT_n = np.array(bbpl_pre_ifl['col5'])
bbpl_pre_kT_p = np.array(bbpl_pre_ifl['col6'])
bbpl_pre_alpha = np.array(bbpl_pre_ifl['col10'])
bbpl_pre_alpha_n = np.array(bbpl_pre_ifl['col11'])
bbpl_pre_alpha_p = np.array(bbpl_pre_ifl['col12'])

#Bkn2pow
bkn2powC_parameters = ascii.read("GRB160509A_bkn2powC_par.txt", header_start = None, comment = '#')
bkn2powC_a1 = np.array(bkn2powC_parameters['col4'])
bkn2powC_a1_n = np.array(bkn2powC_parameters['col5'])
bkn2powC_a1_p = np.array(bkn2powC_parameters['col6'])
bkn2powC_E1 = np.array(bkn2powC_parameters['col7'])		#include this in plot 1
bkn2powC_E1_n = np.array(bkn2powC_parameters['col8'])
bkn2powC_E1_p = np.array(bkn2powC_parameters['col9'])
bkn2powC_a2 = np.array(bkn2powC_parameters['col10'])
bkn2powC_a2_n = np.array(bkn2powC_parameters['col11'])
bkn2powC_a2_p = np.array(bkn2powC_parameters['col12'])
bkn2powC_E2 = np.array(bkn2powC_parameters['col13'])		#include this in plot 1
bkn2powC_E2_n = np.array(bkn2powC_parameters['col14'])
bkn2powC_E2_p = np.array(bkn2powC_parameters['col15'])
bkn2powC_a3 = np.array(bkn2powC_parameters['col16'])
bkn2powC_a3_n = np.array(bkn2powC_parameters['col17'])
bkn2powC_a3_p = np.array(bkn2powC_parameters['col18'])

#3 Comp
comp3_parameters = ascii.read("GRB160509AA_three_comp_parameters.txt", header_start = None, comment = '#')
comp3_kT = np.array(comp3_parameters['col4'])
comp3_kT_n = np.array(comp3_parameters['col5'])
comp3_kT_p = np.array(comp3_parameters['col6'])
comp3_pi1 = np.array(comp3_parameters['col10'])
comp3_pi1_n = np.array(comp3_parameters['col11'])
comp3_pi1_p = np.array(comp3_parameters['col12'])
comp3_E1 = np.array(comp3_parameters['col13'])
comp3_E1_n = np.array(comp3_parameters['col14'])
comp3_E1_p = np.array(comp3_parameters['col15'])
comp3_pi2 = np.array(comp3_parameters['col19'])
comp3_pi2_n = np.array(comp3_parameters['col20'])
comp3_pi2_p = np.array(comp3_parameters['col21'])
comp3_E2 = np.array(comp3_parameters['col22'])
comp3_E2_n = np.array(comp3_parameters['col23'])
comp3_E2_p = np.array(comp3_parameters['col24'])

############################################################################################################################################################################################################################

#print time_main_plot

fig = plt.figure()
ax0 = plt.subplot(1, 1, 1)
axes = plt.gca()

plt.ylabel(r'E$_{p}$',fontsize= 11 )
ax0.get_yaxis().set_label_coords(-0.1,0.5)
plt.xlabel(r'Time since GBM trigger$\/$(s)',fontsize=11)
mask = np.array([False, False, True, True, True, True, False, True, False, False, False, False, False, False, False, False])
plt.ylabel(r'E$_{p}$',fontsize= 11 )
ax0.get_yaxis().set_label_coords(-0.1,0.5)
plt.xlabel(r'Time since GBM trigger$\/$(s)',fontsize=18)

#############################################################################################################################################################################################################################1

# Plotting first plot with error
Ktoffset = 3
E1offset = 1
Epoffset = 1
time_main = np.arange(8, 28, 1)
time_main_plot = time_main[0:len(time_main)-1] + 0.5
time_pre_plot = (np.array([-1.472, 1.778, 2.814, 4.651]) + np.array([1.778, 2.814, 4.651, 8.0]))/2
#Band
plt.errorbar(time_main_plot,band_Ep, yerr=[-1*band_Ep_n, band_Ep_p], fmt = 'o', ls = lstyle, capsize=0, markerfacecolor='b',color='b', ms=7, label='B (Ep)', markeredgewidth=0.5, markeredgecolor='b')
#BandC
plt.errorbar(time_main_plot[0:11],bandC_Ep[0:11], yerr=[-1*bandC_Ep_n[0:11], bandC_Ep_p[0:11]], fmt = 'o', ls = lstyle, capsize=0, markerfacecolor='None',color='b', ms=7, label='BC (Ep)', markeredgewidth=0.5, markeredgecolor='b')
#BB+BandC
plt.errorbar(time_main_plot[0:7],bb_bandC_Ep[0:11], yerr=[-1*bb_bandC_Ep_n[0:11], bb_bandC_Ep_p[0:11]], fmt = 'o', ls = lstyle, capsize=0, markerfacecolor='None',color='firebrick', ms=7, label='BB+BC (Ep)', markeredgewidth=0.5, markeredgecolor='b')
plt.errorbar(time_main_plot[0:7],Ktoffset*bb_bandC_kT[0:11], yerr=[-Ktoffset*bb_bandC_kT_n[0:11], Ktoffset*bb_bandC_kT_p[0:11]], fmt = 's', ls = lstyle, capsize=0, markerfacecolor='None',color='firebrick', ms=7, label='BB+BC (kT)', markeredgewidth=0.5, markeredgecolor='b')
#BB+Band
plt.errorbar(time_main_plot[0:16][mask],bb_band_Ep[mask], yerr=[-1*bb_band_Ep_n[mask], bb_band_Ep_p[mask]], fmt = 'o', ls = lstyle, capsize=0, markerfacecolor='firebrick',color='firebrick', ms=7, label='BB+B (Ep)', markeredgewidth=0.5, markeredgecolor='b')
plt.errorbar(time_main_plot[0:16][mask],Ktoffset*bb_band_kT[mask], yerr=[-Ktoffset*bb_band_kT_n[mask], Ktoffset*bb_band_kT_p[mask]], fmt = 's', ls = lstyle, capsize=0, markerfacecolor='firebrick',color='firebrick', ms=7, label='BB+B (kT)', markeredgewidth=0.5, markeredgecolor='b')
#bkn2pow
plt.errorbar(time_main_plot[0:13], E1offset*bkn2powC_E1[0:13], yerr=[-E1offset*bkn2powC_E1_n[0:13], E1offset*bkn2powC_E1_p[0:13]], fmt = 'rx', ls = lstyle, capsize=0, markerfacecolor='None',color='red', ms=7, label='Bkn2powC (Ep1)', markeredgewidth=0.5, markeredgecolor='b')
plt.errorbar(time_main_plot[0:13],bkn2powC_E2[0:13], yerr=[-1*bkn2powC_E2_n[0:13], bkn2powC_E2_p[0:13]], fmt = 'gx', ls = lstyle, capsize=0, markerfacecolor='None',color='green', ms=7, label='Bkn2powC (Ep2)', markeredgewidth=0.5, markeredgecolor='b')
#3comp
plt.errorbar(time_main_plot[0:13],Ktoffset*comp3_kT[0:13], yerr=[-Ktoffset*comp3_kT_n[0:13], Ktoffset*comp3_kT_p[0:13]], fmt = 'go', ls = lstyle, capsize=0, markerfacecolor='Green',color='green', ms=7, label='3Comp (kT)', markeredgewidth=0.5, markeredgecolor='green')
plt.errorbar(time_main_plot[0:13],comp3_E1[0:13], yerr=[-1*comp3_E1_n[0:13], comp3_E1_p[0:13]], fmt = 'ro', ls = lstyle, capsize=0, markerfacecolor='Yellow',color='firebrick', ms=7, label='3Comp (Ep1)', markeredgewidth=0.5, markeredgecolor='b')
#plt.errorbar(time_main_plot[0:13],comp3_E2[0:13], yerr=[-1*comp3_E2_n[0:13], comp3_E2_p[0:13]], fmt = 'go', ls = lstyle, capsize=0, markerfacecolor='Pink',color='firebrick', ms=7, label='3Comp (Ep2)', markeredgewidth=0.5, markeredgecolor='b')
#Cutoff PL
plt.errorbar(time_pre_plot,Epoffset*cutoffpl_Ep, yerr=[-Epoffset*cutoffpl_Ep_n, Epoffset*cutoffpl_Ep_p], fmt = '*', ls = lstyle, capsize=0, markerfacecolor='Red',color='Cyan', ms=7, label='Cpl (Ep)', markeredgewidth=0.5, markeredgecolor='green')
#Preplot
plt.errorbar(time_pre_plot,Ktoffset*bbpl_pre_kT, yerr=[-Ktoffset*bbpl_pre_kT_n, Ktoffset*bbpl_pre_kT_p], fmt = '*', ls = lstyle, capsize=0, markerfacecolor='Blue',color='magenta', ms=7, label='BBPL (kT)', markeredgewidth=0.5, markeredgecolor='b')
#plt.errorbar(time_pre_plot,cutoffpl_Ep, yerr=[-1*cutoffpl_Ep_n, cutoffpl_Ep_p], fmt = 'P', ls = lstyle, capsize=0, markerfacecolor='b',color='b', ms=7, label='CPL', markeredgewidth=0.5, markeredgecolor='b')


# Plotting the vertical lines on the first plot.
plt.legend(numpoints=1,prop={'size':8})
plt.legend(numpoints=1,prop={'size':11}, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=7, mode="expand", borderaxespad=0.)
plt.ylabel(r'Parameters in keV',fontsize= 20 )
ax0.get_yaxis().set_label_coords(-0.1,0.5)
plt.xlabel(r'Time since GBM trigger$\/$(s)',fontsize=22)
c_time_stamps = np.array([-1.47, 8, 15, 27])		#these are the time bin start points during the burst
ax0.axvline(x=c_time_stamps[0], ymin=0.0, linewidth=1.0, color='r',lw = 1.3, ls='--')
ax0.axvline(x=c_time_stamps[1], ymin=0.0, linewidth=1.0, color='r',lw = 1.3, ls='--')
ax0.axvline(x=c_time_stamps[2], ymin=0.0, linewidth=1.0, color='r',lw = 1.3, ls='--')
ax0.axvline(x=c_time_stamps[3], ymin=0.0, linewidth=1.0, color='r',lw = 1.3, ls='--')
plt.tight_layout()
#~ plt.xscale('log')
plt.yscale('log')
ax0.set_ylim([10,1400])
ax0.yaxis.set_ticks_position('both')
ax0.xaxis.set_ticks_position('both')
#plt.savefig('Ep_and_kT_Evolution_fine_bins.pdf')
plt.tight_layout()
plt.tick_params(labeltop=False, labelright=True)
plt.show()

##############################################################################################################################################################################################################################

#Plotting second plot with error
xmin = bkn2powC_E2_n
xmax = bkn2powC_E2_p
ymin = bkn2powC_E1_n
ymax = bkn2powC_E1_p
x = bkn2powC_E2
y = bkn2powC_E1
rat=x/y
s_x = np.mean([abs(xmin), xmax], axis = 0)
s_y = np.mean([abs(ymin), ymax], axis = 0)
sigma_y = abs(rat)*np.sqrt((s_x**2)/(x**2) + (s_y**2)/(y**2))

uncer = unumpy.uarray(rat,sigma_y)
print uncer

fig = plt.figure(figsize=(6,6))
gs = gridspec.GridSpec(3, 1, height_ratios=[2, 1, 2])

#ax1.set_ylabel(r'Energy [keV]',fontsize=18)
#ax2.set_yscale('log')
#ax6 = plt.subplot(gs[4], sharex = ax1)
#plt.errorbar(time_main_plot, bkn2powC_E2, yerr=[-1*bkn2powC_E2_n, bkn2powC_E2_p], fmt = '*', ls = lstyle, capsize=0, markerfacecolor='firebrick',color='firebrick', ms=7, label='BBPL', markeredgewidth=0.5, markeredgecolor='b')
#ax5.axhline(y=1.0/2,xmin=0,xmax=3,ls='-', c="firebrick",linewidth=2.0, label='jitter 1/2')

ax1 = plt.subplot(gs[0])
plt.errorbar(time_main_plot, bkn2powC_E1, yerr=[-1*bkn2powC_E1_n, bkn2powC_E1_p], fmt = 's', ls = lstyle, capsize=0, markerfacecolor='r',color='cyan', ms=7, label='E$_c$1', markeredgewidth=0.5, markeredgecolor='r')
plt.errorbar(time_main_plot, bkn2powC_E2, yerr=[-1*bkn2powC_E2_n, bkn2powC_E2_p], fmt = 'o', ls = lstyle, capsize=0, markerfacecolor='b',color='crimson', ms=7, label='E$_c$2', markeredgewidth=0.5, markeredgecolor='b')
#ax1.get_yaxis().set_label_coords(-0.08,0.5)
ax1.set_yscale('log')
#ax1.set_ylim([5,200])
ax1.set_ylabel("E$_c$")
ax1.set_ylim([5,1.4e3])

ax2 = plt.subplot(gs[1], sharex = ax1)
ax2.errorbar(time_main_plot, rat, yerr=sigma_y, fmt = '*', ls = lstyle, capsize=0, markerfacecolor='firebrick',color='firebrick', ms=7, label='Ratio', markeredgewidth=0.5, markeredgecolor='b')
ax2.set_ylim([0,13])
ax2.set_ylabel("E$_c$2 / E$_c$1")
ax2.axhline(y=2.0/3,xmin=0,xmax=3, ls='dashdot', c="blue",linewidth=2.0, label='SCS -2/3')

ax3 = plt.subplot(gs[2], sharex = ax1)
plt.errorbar(time_main_plot, bkn2powC_a1, yerr=[-1*bkn2powC_a1_n, bkn2powC_a1_p], fmt = 's', ls = lstyle, capsize=0, markerfacecolor='none',color='red', ms=7, label='Bkn2powC a1', markeredgewidth=0.5, markeredgecolor='r')
plt.errorbar(time_main_plot, bkn2powC_a2, yerr=[-1*bkn2powC_a2_n, bkn2powC_a2_p], fmt = '*', ls = lstyle, capsize=0, markerfacecolor='none',color='green', ms=7, label='Bkn2powC a2', markeredgewidth=0.5, markeredgecolor='g')
plt.errorbar(time_main_plot, bkn2powC_a3, yerr=[-1*bkn2powC_a3_n, bkn2powC_a3_p], fmt = 'o', ls = lstyle, capsize=0, markerfacecolor='none',color='blue', ms=7, label='Bkn2powC a3', markeredgewidth=0.5, markeredgecolor='b')
ax3.set_ylabel("Slopes")
ax3.set_ylim([-2,7])
ax3.axhline(y=3.0/2,xmin=0,xmax=3, ls=':', c="firebrick",linewidth=2.0, label='FCS -3/2')

##############################################################################################################################################################################################################################

ax1.legend(numpoints=1,prop={'size':8})
ax1.legend(numpoints=1,prop={'size':11}, loc=3, ncol=7, borderaxespad=0.)
ax2.legend(numpoints=1,prop={'size':8})
ax2.legend(numpoints=1,prop={'size':11}, loc=9, ncol=7, borderaxespad=0.)
ax3.legend(numpoints=1,prop={'size':8})
ax3.legend(numpoints=1,prop={'size':11}, loc=2, ncol=7, borderaxespad=0.)

ax1.yaxis.set_ticks_position('both')
ax1.tick_params(labeltop=True, labelright=True)
ax2.yaxis.set_ticks_position('both')
ax2.tick_params(labeltop=False, labelright=True)
ax3.yaxis.set_ticks_position('both')
ax3.tick_params(labeltop=False, labelright=True)
#ax4.yaxis.set_ticks_position('both')
#ax4.tick_params(labeltop=False, labelright=True)
ax1.xaxis.set_ticks_position('both')
plt.xlabel(r'Time since GBM trigger$\/$(s)',fontsize=18)
#~ xticklabels =ax7.get_xticklabels()+ax6.get_xticklabels() + ax1.get_xticklabels()+ ax2.get_xticklabels() + ax3.get_xticklabels() + ax5.get_xticklabels() + ax4.get_xticklabels() # +  ax8.get_xticklabels()
#~ yticklabels = ax2.get_yticklabels()  + ax3.get_yticklabels() + ax4.get_yticklabels() + ax5.get_yticklabels() + ax6.get_yticklabels() + ax7.get_yticklabels()+ ax8.get_yticklabels()
#~ plt.setp(yticklabels, visible=False)
#~ plt.setp(xticklabels, visible=False)
plt.subplots_adjust(hspace=0.00)
#~ fig.text(0.03, 0.5, 'Counts/sec', rotation="vertical", va="center", fontsize=18)
#plt.tight_layout()
plt.show()


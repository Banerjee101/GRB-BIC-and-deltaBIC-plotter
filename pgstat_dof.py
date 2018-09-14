import os
import math
import scipy
import numpy as np
import pandas as pd
from astropy.io import ascii
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches	
from matplotlib import rcParams, cycler 
from scipy import stats
import matplotlib as mpl

########################################

font = {
        'weight' : 'bold',
        'size'   : 15}

mpl.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['axes.labelweight']='bold'
mpl.rcParams['axes.labelsize']='large'
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['ytick.major.width'] = 2
lstyle = 'None'


#For reading the data from ASCII
bkn2power_ifl_pgstat, bkn2power_ifl_dof = np.loadtxt("GRB160509A_bkn2powC_stat.txt", unpack=True, usecols=[1,2])
band_ifl_pgstat, band_ifl_dof = np.loadtxt("result_band_stat.txt", unpack=True, usecols=[1,2])
three_comp_ifl_pgstat, three_comp_ifl_dof = np.loadtxt("GRB160509AA_three_comp_stat.txt", unpack=True, usecols=[1,2])
bb_bandC_ifl_pgstat, bb_bandC_ifl_dof = np.loadtxt("results_bandC_stat_freeze.txt", unpack=True, usecols=[1,2])
LC1_rate, LC1_time=np.loadtxt("LC_LLE20000100000.qdp", unpack=True, usecols=[1, 0])
LC2_rate, LC2_time=np.loadtxt("LC_n3_1s830.qdp", unpack=True, usecols=[1, 0])

#print t_plot
#print len(t_plot), len(bkn2power_pgstat)

pha_bins = 476				#Set the PHA bins from data (generalize this)
Nlog = np.log(pha_bins)
range1=8.5				#set the reference time for the plot
x_lim_min=7.5				#set the x limit minimum
x_lim_max=28.0				#set the x limit maximum

#This is for creating BIC 
bkn2power_BIC = bkn2power_ifl_pgstat + (pha_bins-bkn2power_ifl_dof)*Nlog
band_BIC = band_ifl_pgstat + (pha_bins-band_ifl_dof)*Nlog
three_comp_BIC = three_comp_ifl_pgstat + (pha_bins-three_comp_ifl_dof)*Nlog
bb_bandC_BIC = bb_bandC_ifl_pgstat + (pha_bins-bb_bandC_ifl_dof)*Nlog



#calculating the differences
#bkn2power_band=abs(bkn2power_BIC-band_BIC)			#C1
#bkn2power_three_comp=abs(bkn2power_BIC-three_comp_BIC)		#C2
#bkn2power_bb_bandC=abs(bkn2power_BIC-bb_bandC_BIC)		#C3
#band_three_comp=abs(band_BIC-three_comp_BIC)			#C4
#band_bb_bandC=abs(band_BIC-bb_bandC_BIC)			#C5
#three_comp_bb_bandC=abs(three_comp_BIC-bb_bandC_BIC)		#C6

# for loop to compare the values and arrange in ascending order

#Subplots
fig, (ay1, ay2) = plt.subplots(2, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})

plt.xlim(xmin=x_lim_min)
plt.xlim(xmax=x_lim_max)
ay1.plot(np.arange(range1, range1+len(bkn2power_BIC)), bkn2power_BIC, "v-", lw = 1.0, markersize=7, label='bkn2pow' )
ay1.plot(np.arange(range1, range1+len(band_BIC)), band_BIC, "^-", lw = 1.0, markersize=7, label='Band')
ay1.plot(np.arange(range1, range1+len(three_comp_BIC)), three_comp_BIC, "<-", lw = 1.0, markersize=7, label='BB + CPL + CPL')
ay1.plot(np.arange(range1, range1+len(bb_bandC_BIC)), bb_bandC_BIC, ">-", lw = 1.0, markersize=7, label='BB + BandC')
ay11=ay1.twinx()
ay11.axis('off')
ay11.plot(LC1_time, LC1_rate, lw=1.0, ls='steps-post', color='red', label='LC1')
ay12=ay1.twinx()
ay12.axis('off')
ay12.plot(LC2_time, LC2_rate, lw=1.0, ls='steps-post', label='LC2', color ='blue')
ay1.legend(numpoints=1,prop={'size':10},loc="upper right")

#combining the arrays into a text file
col_stack = np.column_stack([bkn2power_BIC, band_BIC, three_comp_BIC, bb_bandC_BIC])
col_stackT = col_stack.transpose()
np.savetxt('BIC_Stack.txt', col_stackT, delimiter='  ')
steps = 0
for h in range(0,np.size(col_stackT,1)):
	list1 = np.loadtxt("BIC_Stack.txt", unpack=True, usecols=[h])
	list2 = np.array(['a','b','c','d'])
	tups = zip(list1, list2); tups.sort(); zip(*tups)
	#print tups
	d1=abs(tups[0][0]-tups[1][0])
	d2=abs(tups[0][0]-tups[2][0])
	d3=abs(tups[0][0]-tups[3][0])
	plt.xlim(xmin=x_lim_min)
	plt.xlim(xmax=x_lim_max)
	ay2.plot((range1+steps),d1,"-*", lw = 1.0, markersize=5, color= 'blue')
	ay2.plot((range1+steps),d2,"-s", lw = 1.0, markersize=5, color= 'red')
	ay2.plot((range1+steps),d3,"-8", lw = 1.0, markersize=5, color= 'green')
	m1,m2,m3,m4=tups[0][1],tups[1][1],tups[2][1],tups[3][1]
	ay2.text((range1+steps),d1, m1+m2, fontsize = 10)
	ay2.text((range1+steps),d2, m1+m3, fontsize = 10)
	ay2.text((range1+steps),d3, m1+m4, fontsize = 10)
	red_patch = mpatches.Patch(color='red', label='a=bkn2pow')
	blue_patch = mpatches.Patch(color='blue', label='b=Band')
	yellow_patch = mpatches.Patch(color='yellow', label='c=BB + CPL + CPL')
	green_patch = mpatches.Patch(color='green', label='d=BB + BandC')
	ay2.legend(handles=[red_patch, blue_patch, yellow_patch, green_patch], fontsize = 9)	
	steps=steps+1
#ay2.plot(np.arange(8.5, 8.5+len(bkn2power_BIC)),d1,"*-", lw = 1.0, markersize=5, label='bkn2pow - Band')
#ay2.plot(np.arange(8.5, 8.5+len(bkn2power_three_comp)),bkn2power_three_comp,"s-", lw = 1.0, markersize=5, label='bkn2pow - (BB + CPL + CPL)')
#ay2.plot(np.arange(8.5, 8.5+len(bkn2power_bb_bandC)),bkn2power_bb_bandC,"8-", lw = 1.0, markersize=5, label='bkn2pow - (BB + Band)')
#ay2.plot(np.arange(8.5, 8.5+len(band_three_comp)),band_three_comp,"o-", lw = 1.0, markersize=5, label='Band - (BB + CPL + CPL)')
#ay2.plot(np.arange(8.5, 8.5+len(band_bb_bandC)),band_bb_bandC,"d-", lw = 1.0, markersize=5, label='Band - (BB + Band)')
#ay2.plot(np.arange(8.5, 8.5+len(three_comp_bb_bandC)),three_comp_bb_bandC,"h-", lw = 1.0, markersize=5, label='(BB + CPL + CPL) - (BB + Band)')
#ay2.legend(numpoints=1,prop={'size':10},loc="upper right")

ay2.axhline(y=2, color='black', lw=0.7)
ay2.axhline(y=6, color='black', lw=0.7)
ay2.axhline(y=10, color='black', lw=0.7)
ay2.set_xlabel('time (s)')
ay1.set_ylabel('BIC')
ay2.set_ylabel('delta(BIC)')
#~ plt.xscale('log')
#plt.tight_layout()
plt.show()

import os
import math
import scipy
import texttable as tt
import itertools as itt
import numpy as np
import pandas as pd
from astropy.io import ascii
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches	
from matplotlib import rcParams, cycler 
from scipy import stats
import matplotlib as mpl

#####################################################################################################################################################

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
mo1_ifl_pgstat, mo1_ifl_dof = np.loadtxt("GRB160509A_bkn2powC_stat.txt", unpack=True, usecols=[1,2])
mo2_ifl_pgstat, mo2_ifl_dof = np.loadtxt("result_band_stat.txt", unpack=True, usecols=[1,2])
mo3_ifl_pgstat, mo3_ifl_dof = np.loadtxt("GRB160509AA_three_comp_stat.txt", unpack=True, usecols=[1,2])
mo4_ifl_pgstat, mo4_ifl_dof = np.loadtxt("results_bandC_stat_freeze.txt", unpack=True, usecols=[1,2])

LC1_rate, LC1_time=np.loadtxt("LC_LLE20000100000.qdp", unpack=True, usecols=[1, 0])
LC2_rate, LC2_time=np.loadtxt("LC_n3_1s830.qdp", unpack=True, usecols=[1, 0])

#####################################################################################################################################################

#setting the input parameters
pha_bins = 476									#Set the PHA bins from data (generalize this)
range1=8.5									#set the reference time for the plot
x_lim_min=7.5									#set the x limit minimum
x_lim_max=28.0									#set the x limit maximum
corr_lim=6.0									#Set the correlation limit here
mo1,mo2,mo3,mo4='bkn2pow', 'Band', 'BB + CPL + CPL', 'BB + BandC'		#Set the 4 model names here, should be the same as input order

#####################################################################################################################################################

#This is for creating BIC 
Nlog = np.log(pha_bins)
mo1_BIC = mo1_ifl_pgstat + (pha_bins-mo1_ifl_dof)*Nlog
mo2_BIC = mo2_ifl_pgstat + (pha_bins-mo2_ifl_dof)*Nlog
mo3_BIC = mo3_ifl_pgstat + (pha_bins-mo3_ifl_dof)*Nlog
mo4_BIC = mo4_ifl_pgstat + (pha_bins-mo4_ifl_dof)*Nlog

#####################################################################################################################################################

#Subplot section that makes the 
fig, (ay1, ay2) = plt.subplots(2, sharex=True, gridspec_kw = {'height_ratios':[2, 1]})

plt.xlim(xmin=x_lim_min)
plt.xlim(xmax=x_lim_max)
ay1.plot(np.arange(range1, range1+len(mo1_BIC)), mo1_BIC, "v-", lw = 1.0, markersize=7, label=mo1)
ay1.plot(np.arange(range1, range1+len(mo2_BIC)), mo2_BIC, "^-", lw = 1.0, markersize=7, label=mo2)
ay1.plot(np.arange(range1, range1+len(mo3_BIC)), mo3_BIC, "<-", lw = 1.0, markersize=7, label=mo3)
ay1.plot(np.arange(range1, range1+len(mo4_BIC)), mo4_BIC, ">-", lw = 1.0, markersize=7, label=mo4)
ay11=ay1.twinx()
ay11.axis('off')
ay11.plot(LC1_time, LC1_rate, lw=1.0, ls='steps-post', color='red', label='LC1')
ay12=ay1.twinx()
ay12.axis('off')
ay12.plot(LC2_time, LC2_rate, lw=1.0, ls='steps-post', label='LC2', color ='blue')
ay1.legend(numpoints=1,prop={'size':10},loc="upper right")

#combining the arrays into a text file
col_stack = np.column_stack([mo1_BIC, mo2_BIC, mo3_BIC, mo4_BIC])
col_stackT = col_stack.transpose()
np.savetxt('BIC_Stack.txt', col_stackT, delimiter='  ')
steps = 0
table=tt.Texttable(max_width=0)
headings=['Time Stamps', 'Strongest Correlation Model', 'Delta BIC 1', 'Delta BIC 2', 'Delta BIC 3', 'Models with ascending order of BICs' ]
table.header(headings)
for h in range(0,np.size(col_stackT,1)):
	list1 = np.loadtxt("BIC_Stack.txt", unpack=True, usecols=[h])
	list2 = np.array(['a','b','c','d'])
	list3 = np.array([mo1, mo2, mo3, mo4])
	tups = zip(list1, list2, list3); tups.sort(); zip(*tups)
	d1=abs(tups[0][0]-tups[1][0])
	d2=abs(tups[0][0]-tups[2][0])
	d3=abs(tups[0][0]-tups[3][0])
	time_stamps=range1+steps
	plt.xlim(xmin=x_lim_min)
	plt.xlim(xmax=x_lim_max)
	ay2.plot(time_stamps,d1,"-*", lw = 1.0, markersize=5, color= 'blue')
	ay2.plot(time_stamps,d2,"-s", lw = 1.0, markersize=5, color= 'red')
	ay2.plot(time_stamps,d3,"-8", lw = 1.0, markersize=5, color= 'green')
	m1,m2,m3,m4=tups[0][1][0],tups[1][1][0],tups[2][1][0],tups[3][1][0]
	ay2.text(time_stamps,d1, m1+m2, fontsize = 10)
	ay2.text(time_stamps,d2, m1+m3, fontsize = 10)
	ay2.text(time_stamps,d3, m1+m4, fontsize = 10)
	red_patch = mpatches.Patch(color='red', label='a='+mo1)
	blue_patch = mpatches.Patch(color='blue', label='b='+mo2)
	yellow_patch = mpatches.Patch(color='yellow', label='c='+mo3)
	green_patch = mpatches.Patch(color='green', label='d='+mo4)
	ay2.legend(handles=[red_patch, blue_patch, yellow_patch, green_patch], fontsize = 9)
	steps=steps+1
	
	#Section to sort the strongest correlation model
	m11,m22,m33,m44=tups[0][2],tups[1][2],tups[2][2],tups[3][2]
	model_conv=''
	if d1<=corr_lim:
		model_conv=m11+' and '+m22
	elif d1>=corr_lim:
		model_conv=' No strong correlation '
	for row in zip(itt.repeat(time_stamps, np.size(col_stackT,1)), [model_conv], [d1], [d2], [d3], [m11+' < '+m22+' < '+m33+' < '+m44]):
		table.add_row(row)
s = table.draw()
print s
#ascii.write(s, 'TestTable1.txt', format='latex')	

#####################################################################################################################################################	

#This is the plotting section
ay2.axhline(y=2, color='black', lw=0.7)
ay2.axhline(y=6, color='black', lw=0.7)
ay2.axhline(y=10, color='black', lw=0.7)
ay2.set_xlabel('time (s)')
ay1.set_ylabel('BIC')
ay2.set_ylabel('delta(BIC)')
#~ plt.xscale('log')
#plt.tight_layout()
plt.show()

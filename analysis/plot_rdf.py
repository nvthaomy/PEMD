""" overlay density profiles from multiple AA MDs of different Uext
need to generate data file with python ~/bin/PEMD/analysis/1d-histogram.py trajectory.dcd topology -axis  -atoms 
 """

from matplotlib import cm, ticker
import numpy as np
import os, sys, glob
import matplotlib
import re
showPlots = True
try:
  os.environ["DISPLAY"] #Detects if display is available
except KeyError:
  showPlots = False
  matplotlib.use('Agg') #Need to set this so doesn't try (and fail) to open interactive graphics window
import matplotlib.pyplot as plt
#plt.style.use('seaborn-dark')
matplotlib.rc('font', size=7)
matplotlib.rc('axes', titlesize=7)
colors = ['#6495ED','r','#6da81b','#483D8B','#FF8C00','#2E8B57','#800080','#008B8B','#1E90FF'] 
################################

plotTitle = 'rdf vs NaCl xp0.1 N12 f1 LJPME 298K'
ylabel = 'g(r)'
xlabel = 'r (nm)'
legends = [0,100,500,1000,2000,4000]
Dirs = [
'../xp0.1_N12_f1_V157_LJPME_298K/f1/w0.13/NPT/run2',
'../xp0.1_N12_f1_9NaCl_V157_LJPME_298K/f1/w0.13',
'../xp0.1_N12_f1_47NaCl_V157_LJPME_298K/f1/w0.13',
'../xp0.1_N12_f1_94NaCl_V157_LJPME_298K/f1/w0.13',
'../xp0.1_N12_f1_189NaCl_V157_LJPME_298K/f1/w0.13',
'../xp0.1_N12_f1_379NaCl_V157_LJPME_298K/f1/w0.13']


#x = [0,0.17, 0.25, 0.33, 0.5, 0.67, 0.75, 0.83, 0.92, 1.]
#Dirs = [
#'../xp0.1_N12_f0_V157_LJPME_298K/f0/w0.13/NPT/run2',
#'../xp0.1_N12_f0.17_V157_LJPME_298K/f0.17/w0.13/NPT/run2',
#'../xp0.1_N12_f0.25_V157_LJPME_298K/f0.25/w0.13/NPT',
#'../xp0.1_N12_f0.33_V157_LJPME_298K/f0.33/w0.13/run1',
#'../xp0.1_N12_f0.5_V157_LJPME_298K/f0.5/w0.13/NPT/run2',
#'../xp0.1_N12_f0.67_V157_LJPME_298K/f0.67/w0.13/NPT',
#'../xp0.1_N12_f0.75_V157_LJPME_298K/f0.75/w0.13/NPT',
#'../xp0.1_N12_f0.83_V157_LJPME_298K/f0.83/w0.13/NPT',
#'../xp0.1_N12_f0.92_V157_LJPME_298K/f0.92/w0.13/NPT',
#'../xp0.1_N12_f1_V157_LJPME_298K/f1/w0.13/NPT/run2']

dataFName = 'rdf_OD1_Na+.dat'
rColId = 0
rdfColId = 1
##################
#check
if len(legends) != len(Dirs):
    Exception('Mismatch in sizes of legends and Dirs')
rdfs = []
rs = []
cwd = os.getcwd()
for i, dir in enumerate(Dirs): 
    f = open(os.path.join(dir+'/',dataFName), 'r')
    data = np.loadtxt(f)
    r = data[:,rColId]
    rdf = data[:,rdfColId] 
    rdfs.append(rdf)
    rs.append(r)

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
ax.set_prop_cycle('color', colors)
for i,legend in enumerate(legends):
    ax.plot(rs[i], rdfs[i], marker='None', ls='-', lw=1, ms=4,label = legend)
plt.ylabel(ylabel)
plt.xlabel(xlabel)
plt.legend(loc='best')
title = plotTitle
plt.title(title,loc='center')
plt.savefig('_'.join(title.split())+'.png',dpi=500,bbox_inches='tight')
plt.show()
    

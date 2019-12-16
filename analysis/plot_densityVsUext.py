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

plotExt = 'Uext on NaCl'
uexts = [0.5,1,3]
uextDirs = ['uext0.5_nacl','uext1_nacl','uext3_nacl']
histDict = {'NaCl':['zhist_Na+_Cl-.dat',1,'Number density'], 'All':['zhist_all.dat',1,'Count']}
dataDict = {}
cwd = os.getcwd()

for atom, [hist,col,ylabel] in histDict.items():
    for i, dir in enumerate(uextDirs):
        filehandler = open(os.path.join(cwd,dir+'/',hist),'r')
        filehandler.seek(0)
        AllData=np.loadtxt(filehandler)
        Data=AllData[:,col]
        x = AllData[:,0]
        #write bin values and histogram into dataDict
        dataDict.update({(atom, uexts[i]):[x,Data]})        

for key in sorted(histDict.keys()):
    fig,ax = plt.subplots(nrows=1, ncols=1, figsize=[3,2])
    ax.set_prop_cycle('color', colors)
    for (atom, uext), [x,hist] in sorted(dataDict.items()):
        if atom == key:
            ax.plot(x, hist ,marker='o',ls='-',lw=1 ,mfc="None",ms=0, label = '$U_{ext}=$'+'{}'.format(uext))
    ax.legend(loc='best',prop={'size': 5})
    plt.ylabel(histDict[key][2])
    plt.xlabel('x')
    title = 'hist {} {}'.format(key,plotExt)
    plt.title(title,loc='center')
    plt.savefig('_'.join(title.split())+'.png',dpi=500,bbox_inches='tight')
plt.show()
    

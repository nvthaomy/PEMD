#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 17:56:07 2019

@author: nvthaomy
"""
import random
import numpy as np
import math
import os
from subprocess import call

""" write out tleap input file to build a single polyelectrolyte chain with specified
    tacticity and various fraction of deprotonation
    input: N: degree of polymerization
           f: degree of ionization
           Pm: fraction of meso diads"""
class Molecule:
    def __init__(self,N,f,Pm,pattern,headName,tailName):
        self.DOP = N
        self.f = f
        self.Pm = Pm
        self.pattern = pattern
        self.headName = headName
        self.tailName = tailName
    def GetIndices(self,num1,num2,N):
        """get indices of ionized monomers if f<=0.5 or neutral monomers if f>0.5
           N: degree of polymerization
           num2: number of ionized monomers (for f<=0.5) or neutral monomer (for f>0.5)
           num1: list of upper and lower bounds of number of neutral monomers if f  <= 0.5 or 
                number of ionized monomers if f > 0.5
           output: list of indices"""
        i=[] #list of indices of deprotonated monomers
        num2 = (round(num2))
        if len(num1) > 1: #for case when num is not an integer or divisible by 5
            i.append(random.randint(0,min(num1)))
            print ("place the first deprot (or prot if f>0.5) monomer at position %i" %(i[0]+1))
            counter = 0 
            old_j = 0
            for j in range(1,int(num2)):
                spacing = num1[random.randint(0,len(num1)-1)] #draw random spacing from the list of upper and lower bounds of number of protonated(deprotonated) monomers 
                i.append(i[j-1] + spacing + 1)
                while i[j] > N-1:
                    i[j] -= 1        
                    if j != old_j: #only count if the next monomer index exceed the monomer chain length
                        counter += 1 
                        if counter > 1:
                            print ("\nWarning: Index of monomer is out of range!")
                    old_j = j 
        else:
            i.append(0)
            counter = 0
            old_j = 0
            print ("place the first deprot/prot monomer at position %i" %(i[0]+1))
            for j in range(1,int(num2)):
                i.append(i[j-1] + num1[0] + 1)
                while i[j] > N-1:
                        i[j] -= 1        
                        if j != old_j: #only count if the next monomer index exceed the monomer chain length
                            counter += 1 
                            if counter > 1:
                                print ("\nWarning: Index of monomer is out of range!")
                        old_j = j 
        i = [int(x) for x in i] #convert indices to integer
        return i
            
    def UniformChargePattern(self):
        print ("\nGenerating evenly distributed charge pattern for f = %3.2f" %self.f) 
        N=self.DOP
        f = self.f
        if f <= 0.5:
            nIonized = (f * N)
            nNeutral = [(N/nIonized - 1)] 
            print ("\nNumber of neutral monomers between each ionized monomer is: %3.3f" %nNeutral[0])
            if not abs(nNeutral[0]-round(nNeutral[0])) < 0.0001: #not an integer
                if nNeutral[0]*10 % 5 == 0: #screen out charge fraction that result in nNeutral= x.5
                   nNeutral = [math.ceil(nNeutral[0])]
                else:
                    nMin = math.floor(nNeutral[0])
                    nMax = math.ceil(nNeutral[0])
                    mean = nNeutral[0] #average number of protonated moners in between two deprotonated ones
                    nNeutral = [nMin]
                    calc_mean = sum(nNeutral)/len(nNeutral)
                    #generate list of possible number of neutral monomers between ionized monomers 
                    while abs(calc_mean - mean) > 0.001:
                        if calc_mean > mean:
                            nNeutral.append(nMin)
                            calc_mean = sum(nNeutral)/len(nNeutral)
                        else:
                            nNeutral.append(nMax)
                            calc_mean = sum(nNeutral)/len(nNeutral)
                        print ('\nNumber of neutral monomers between ionized monomers will be drawn from this list {}'.format(nNeutral))
                
            iIonized = self.GetIndices(nNeutral,nIonized,N)
            print ("\n Indices of ionized monomer for f = {:3.2f} are {}".format(f,iIonized))
            if abs(len(iIonized)-nIonized) > 0.001:
                raise Exception("Number of ionized monomers does not match the value of charge fraction!")
            chargePattern = N*['n']
            for i in iIonized:
                chargePattern[i] = 'i'
            
        elif f > 0.5:
            nNeutral = ((1 - f) * N)
            nIonized = [(N/nNeutral - 1)]
            print ("\nNumber of ionized monomers between each neutral monomer is: %3.3f" %nIonized[0])
            if not abs(nIonized[0]-round(nIonized[0])) < 0.0001: #not an integer
                if nIonized[0]*10 % 5 == 0: #screen out charge fraction that result in nIonized= x.5
                    nIonized = [math.ceil(nIonized[0])]
                else:
                    nMin = math.floor(nIonized[0])
                    nMax = math.ceil(nIonized[0])
                    mean = nIonized[0] #average number of ionized moners in between two protonated ones
                    nIonized = [nMin]
                    calc_mean = sum(nIonized)/len(nIonized)
                    while abs(calc_mean - mean) > 0.001:
                        if calc_mean > mean:
                            nIonized.append(nMin)
                            calc_mean = sum(nIonized)/len(nIonized)
                        else:
                            nIonized.append(nMax)
                            calc_mean = sum(nIonized)/len(nIonized)
                    print ('\nNunmber of ionized monomers will be drawn from this list: {}'.format(nIonized))
            iNeutral = self.GetIndices(nIonized,nNeutral,N)
            print ("\nIndices of protonated monomer for f = {:3.2f} are {}".format(f,iNeutral))
            if abs(len(iNeutral)-nNeutral) > 0.001:
                raise Exception("Number of neutral monomers does not match the value of charge fraction!\n")
            chargePattern = N*['i']
            for i in iNeutral:
                chargePattern[i] = 'n' 
        return chargePattern
    
    def ChargePattern(self):
        """Evaluate if charge pattern is random or evenly distributed
        and enerate charge pattern"""
        N=self.DOP
        f = self.f
        pattern = self.pattern
        f0= N*['n']
        f1= N*['i']
        if f == 0.:
            chargePattern = f0
        elif f == 1.:
            chargePattern = f1
        else:
            if pattern == 'random':
                print ("\nGenerating randomly distributed charge pattern for f = %3.2f" %f) 
                nIonized = N*f
                chargePattern = f0
                iIonized = random.sample(range(N),int(nIonized))
                for i in iIonized:
                    chargePattern[i] = 'i' 
            else:
                chargePattern = self.UniformChargePattern()
        print('Charge pattern {}'.format(chargePattern))
        return chargePattern           
                
        
    def GetTactCharge(self):
        """Append tacticity with charge pattern"""
        tact_charge=[]
        #get tacticity and charge pattern
        tacticity = self.Tacticity()
        chargePattern = self.ChargePattern()
        for i in range(0,self.DOP):
            a = tacticity[i] + chargePattern[i]
            tact_charge.append(a)
        tact_charge[0]=tact_charge[0].upper()
        tact_charge[-1]=tact_charge[-1].upper()
        return tact_charge
        
    def Tacticity(self):
            N=self.DOP
            Pm = self.Pm
            dyad=[]
            tact=[]
            #follow Bernoullian statistics (common with free radical polymerization)
            #stereo of the next monomer is independent of the stereochemistry of growing chain
            for i in range(0,N-3): # -1 due to get number of dyads, -2 due to the end monomers do not have stereochemistry
                rand = np.random.random()
                if rand <= Pm:
                        dyad.append('m')
                else:
                        dyad.append('r')
            #count triad fractions
            dyad_str = ''.join(dyad)
            f_m =  dyad_str.count('m')/len(dyad)
            f_mm = dyad_str.count('mm')/(len(dyad)-1)
            f_rr = dyad_str.count('rr')/(len(dyad)-1)
            f_mr = 1 - f_mm - f_rr
            print('Diad sequence (excluding end monomers) from meso diad fraction of {}:\n{}'.format(Pm,dyad))
            print('Actual f_m {}'.format(f_m))
            print('f_mm {}'.format(f_mm))
            print('f_rr {}'.format(f_rr))
            print('f_mr + f_rm {}'.format(f_mr))
            tact.append(self.headName) #head group is achiral
            tact.append('u') #arbitrarily pick the second monomer to be "up"
            for i in range(0,len(dyad)):
                if dyad[i] == 'm':
                    tact.append(tact[-1])
                else:
                    if tact[-1] == 'u':
                        tact.append('d')
                    else:
                        tact.append('u')
            tact.append(self.tailName)
            return tact
def BuildPoly(f,N,Pm,pattern, monName, neutralMons, ionizedMons,headName,tailName,libs,ff,runTleap=True):
    with open("build_AA"+str(N)+".log", 'w') as log:     
        Poly = Molecule(N,f,Pm,pattern,headName,tailName)
        molName = monName + str(N) 
        log.write("\nCalculating number of deprotonated monomers and"\
                  "modifying input fraction of deprotonation:")        
        if f !=0:
            nIonized = float(f * N)
            if not nIonized.is_integer(): 
                log.write("\nNeed to modify f = %3.2f to " %f)
                nIonized = round(nIonized)
                f = nIonized/N #new charge fraction
                log.write("f = %3.2f" %f)
        log.write("\nNew DOI is : %3.2f"%f)
        log.write('\nBuilding {} of with N = {}, f = {:3.2f}, charge pattern = {},'\
                  ' meso fraction  = {}'.format(molName,N,f,pattern,Pm))

        #get list consists of strings defining ionization and tacticity of monomers
        tact_charge = Poly.GetTactCharge()
        
        #fix name of head and tail monomers
        if monName == 'AA': #PAA
            if tact_charge[0][-1] == 'I':
                tact_charge[0] = tact_charge[0][0:-1] + 'D'
            elif tact_charge[0][-1] == 'N':
                tact_charge[0] = tact_charge[0][0:-1] + 'P'
            if tact_charge[-1][-1] == 'I':
                tact_charge[-1] = tact_charge[-1][0:-1] + 'D'
            elif tact_charge[-1][-1] == 'N':
                tact_charge[-1] = tact_charge[-1][0:-1] + 'P'

        elif monName == 'AH': #PAH
            if tact_charge[0][-1] == 'I':
                tact_charge[0] = tact_charge[0][0:-1] + 'P'
            elif tact_charge[0][-1] == 'N':
                tact_charge[0] = tact_charge[0][0:-1] + 'D'
            if tact_charge[-1][-1] == 'I':
                tact_charge[-1] = tact_charge[-1][0:-1] + 'P'
            elif tact_charge[-1][-1] == 'N':
                tact_charge[-1] = tact_charge[-1][0:-1] + 'D'
        else:
            raise Exception('Do not recognize monomer name {}'.format(monName))
        sequence = ' '.join(tact_charge)    
        pdbOut = '{}_f{:3.2f}.pdb'.format(molName,f)  

        buildFile = ''
        if runTleap:
            log.write("\nWriting tleap input file to build a single polymer with different deprotonated fraction")
            libstr = ""
            for lib in libs:
                libstr += 'loadOFF {}\n'.format(lib)
            s = """source leaprc.{ff}
{libstr}un =loadpdb {neu1}
ui =loadpdb {ion1}
dn =loadpdb {neu2}
di =loadpdb {ion2}
set un head un.1.1
set un tail un.1.3
set ui head ui.1.1
set ui tail ui.1.3
set dn head dn.1.1
set dn tail dn.1.3
set di head di.1.1
set di tail di.1.3

#f = {f}
x = sequence{{{seq}}}
savepdb x {pdbOut}
quit""".format(libstr=libstr, neu1 = neutralMons[0], ion1 = ionizedMons[0], neu2 = neutralMons[1], ion2 = ionizedMons[1], f = f, seq = sequence, pdbOut = pdbOut, ff=ff)

            buildFile = "build_"+str(molName)+".in"
            file = open(buildFile,"w")
            file.write(s)
            log.write("\nDone writing tleap input file")
            cmd = 'tleap -f {} > build_{}.out'.format(buildFile,molName)
        log.close()
    return buildFile, pdbOut, f, sequence
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-ff", type=str, default = 'gaff2', help="gaff2 or ff99")
    parser.add_argument("-f", type=float, required=True, help="degree of ionization")
    parser.add_argument("-N",type=int, required=True,
                        help="degree of polymerization")
    parser.add_argument("-Pm","--Pm",type=float,default = 1, 
                        help="Meso diad fraction, isotactic if = 1, syndiotactic if = 0")
    parser.add_argument("-r","--random", action = "store_true",
                        help="Random deprotonation, default pattern is random")
    parser.add_argument("-e","--evendist", action = "store_true",
                        help="Evenly distributed deprotonation")
    parser.add_argument("-n", help='monomer name')
    parser.add_argument("-pdbN", nargs=2, default=[], help = 'pdbs of two neutral chiral monomers')
    parser.add_argument("-pdbI", nargs=2, default=[], help = 'pdbs of two ionized chiral monomers')
    parser.add_argument("-hn", help = 'name of head monomer (AH, NH)')
    parser.add_argument("-tn", help = 'name of tail monomer (AT,NT)')
    parser.add_argument("-l",nargs ="+", default = None,
                        help="path to tleap libraries for monomers")
    args = parser.parse_args()    
    f = args.f
    N = args.N
    Pm = args.Pm
    monName = args.n
    neutralMons = args.pdbN
    ionizedMons = args.pdbI
    headName = args.hn
    tailName = args.tn
    libs =  args.l
    pattern = 'random' #random deprotonation
    ff = args.ff
    if args.evendist:
        pattern = 'even'
    BuildPoly(f,N,Pm,pattern, monName, neutralMons, ionizedMons, headName,tailName,libs,ff)


     




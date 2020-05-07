#!/bin/bash

#polyelectrolytes
monNames=(AA AH)
fs=(0.7 0.7)
Ns=(10 10)
Pms=(1 1)
pattern=(-e -e)
monNames=(AA AH)
pdbNeus=(AA_prot1.pdb AA_prot2.pdb AH_deprot1.pdb AH_deprot2.pdb)
pdbIons=(AA_deprot1.pdb AA_deprot2.pdb AH_prot1.pdb AH_prot2.pdb)
headNames=(AH NH)
tailNames=(AT NT)
libs=(PAAP9N9c5.lib PAAD9N9c4.lib PAHP1N11c5.lib PAHD1N11c6.lib) 
ffs=(gaff2 ff99)
water=opc

#packmol
PAA=AA10_f0.70.pdb
PAH=AH10_f0.70.pdb
hoh=${water}.pdb
na=na.pdb
cl=cl.pdb

pdbs=($PAA $PAH $hoh $na $cl)
nMols=(1 1 100 20 20)
mixName='mixture'
x=2
y=2
z=2

#loadFF
mixPdb=${mixName}.pdb

#write sim
top=${mixName}.parm7
crd=${mixName}.crd
jobName=testSim
steps=100
coolingsteps=100
meltingsteps=100
anisoP= #-a
stride=10
checkstride=10
temp=298.15
melttemp=400
pressure= #empty for NVT, -p {pressure in atm} for NPT

nPE=${#monNames[@]}
for ((i=0;i<$nPE;i++)); do
    buildFile=build_${monNames[$i]}${Ns[$i]}.in
    buildOut=build_${monNames[$i]}${Ns[$i]}.out
    k1=$(expr ${i}*2| bc)
    k2=$(expr ${i}*2+1| bc)
    echo === building ${monNames[$i]}${Ns[$i]} ===
    python singleChain.py -f ${fs[$i]} -N ${Ns[$i]} -Pm ${Pms[$i]} ${pattern[$i]} -n ${monNames[$i]} -pdbN ${pdbNeus[$k1]} ${pdbNeus[$k2]} -pdbI ${pdbIons[$k1]} ${pdbIons[$k2]} -hn ${headNames[$i]} -tn ${tailNames[$i]} -l ${libs[$k1]} ${libs[$k2]} -ff ${ffs[$i]}

    while [ ! -f $buildFile ]; do sleep 1; done
    tleap -s -f $buildFile > $buildOut
    while [ ! -f $buildOut ]; do sleep 1; done
done

echo === packing molecules ===
python packmol.py -n ${nMols[@]} -pdb ${pdbs[@]} -xyz $x $y $z -o $mixName
packmol < ${mixName}.inp
while [ ! -f $mixPdb ]; do sleep 1; done

echo === loading forcefield ===
python loadFF.py -pdb $mixPdb -w $water -l ${libs[@]} -ff ${ffs[@]}
while [ ! -f loadFF.in ]; do sleep 1; done
tleap -s -f loadFF.in > loadFF.out
while [ ! -f ${mixName}.parm7 ]; do sleep 1; done

echo === writing simulation script ===
python writeSim.py $top $crd -jn $jobName -s $steps -ms $meltingsteps -cs $coolingsteps -xyz $x $y $z $anisoP -stride $stride -cstride $checkstride -temp $temp -mtemp $melttemp $pressure 

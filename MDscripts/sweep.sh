#!/bin/bash

ffDir=/home/mnguyen/bin/PEMD/MDscripts/lib/
mainDir=/home/mnguyen/bin/PEMD/MDscripts/
#polyelectrolytes
monNames=(AA AH)
fs=(1.0 1.0)
Ns=(24 24)
Pms=(1 1)
pattern=(-e -e)
pdbNeus=(${ffDir}AA_prot1.pdb ${ffDir}AA_prot2.pdb ${ffDir}AH_deprot1.pdb ${ffDir}AH_deprot2.pdb)
pdbIons=(${ffDir}AA_deprot1.pdb ${ffDir}AA_deprot2.pdb ${ffDir}AH_prot1.pdb ${ffDir}AH_prot2.pdb)
headNames=(AH NH)
tailNames=(AT NT)
libs=(${ffDir}PAAP9N9c5.lib ${ffDir}PAAD9N9c4.lib ${ffDir}PAHP1N11c5.lib ${ffDir}PAHD1N11c6.lib) 
ffs=(gaff2 ff99)
water=opc

#packmol
PAA=AA24_f1.00.pdb
PAH=AH24_f1.00.pdb
hoh=${ffDir}${water}.pdb
na=${ffDir}na.pdb
cl=${ffDir}cl.pdb

pdbs=($PAA $PAH $hoh $na $cl)
nMols=(8 8 9849 0 0)
mixName='xp0.1_8AA24f1_8AH24f1_9849hoh'
x=5.4
y=5.4
z=30.

#loadFF
mixPdb=${mixName}.pdb

#write sim
top=${mixName}.parm7
crd=${mixName}.crd
jobName=testSim
steps=200000000 #2e8
coolingsteps=4000000 #4e6
meltingsteps=2000000 #2e6
anisoP=-a
stride=5000
checkstride=1000
temp=298.15
melttemp=400
pressure='-p 1.0' #empty for NVT, -p {pressure in atm} for NPT

nPE=${#monNames[@]}

mkdir build
cd build
for ((i=0;i<$nPE;i++)); do
    buildFile=build_${monNames[$i]}${Ns[$i]}.in
    buildOut=build_${monNames[$i]}${Ns[$i]}.out
    k1=$(expr ${i}*2| bc)
    k2=$(expr ${i}*2+1| bc)
    echo === building ${monNames[$i]}${Ns[$i]} ===
    python ${mainDir}singleChain.py -f ${fs[$i]} -N ${Ns[$i]} -Pm ${Pms[$i]} ${pattern[$i]} -n ${monNames[$i]} -pdbN ${pdbNeus[$k1]} ${pdbNeus[$k2]} -pdbI ${pdbIons[$k1]} ${pdbIons[$k2]} -hn ${headNames[$i]} -tn ${tailNames[$i]} -l ${libs[$k1]} ${libs[$k2]} -ff ${ffs[$i]}

    while [ ! -f $buildFile ]; do sleep 1; done
    echo "tleap -s -f $buildFile > $buildOut"
    tleap -s -f $buildFile > $buildOut
    while [ ! -f $buildOut ]; do sleep 1; done
done

echo === packing molecules ===
python ${mainDir}packmol.py -n ${nMols[@]} -pdb ${pdbs[@]} -xyz $x $y $z -o $mixName
echo "packmol < ${mixName}.inp"
packmol < ${mixName}.inp
while [ ! -f $mixPdb ]; do sleep 1; done

echo === loading forcefield ===
python ${mainDir}loadFF.py -pdb $mixPdb -w $water -l ${libs[@]} -ff ${ffs[@]}
while [ ! -f loadFF.in ]; do sleep 1; done
echo "tleap -s -f loadFF.in > loadFF.out"
tleap -s -f loadFF.in > loadFF.out
while [ ! -f ${mixName}.parm7 ]; do sleep 1; done

echo === writing simulation script ===
python ${mainDir}writeSim.py $top $crd -jn $jobName -s $steps -ms $meltingsteps -cs $coolingsteps -xyz $x $y $z $anisoP -stride $stride -cstride $checkstride -temp $temp -mtemp $melttemp $pressure 

cd ../; mkdir sim; cd sim
cp ../build/$top .; cp ../build/$crd .; cp ../build/sim.py .; cp ../build/run.sh .

#!/bin/bash
#SBATCH --ignore-pbs
#SBATCH --nodes=1 --partition=gpu --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=72:00:00
#SBATCH --job-name=3M_6
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=my@ucsb.edu
/bin/hostname
srun --gres=gpu:1 /usr/bin/nvidia-smi
cd $SLURM_SUBMIT_DIR

export PATH="/home/mnguyen/miniconda3/bin:$PATH"

maindir=$(pwd)
libDir=/home/mnguyen/bin/PEMD/MDscripts
NPTdir=sim/run1
top=/home/mnguyen/PE/PAA_PAH/AA/xp0.05_7AA12f1_7AH12f1_592nacl/modifiedO-Na+/6.0.263sigma_0.635eps/build/xp0.05_7AA12f1_7AH12f1_592nacl_10396hoh.parm7
crd=/home/mnguyen/PE/PAA_PAH/AA/xp0.05_7AA12f1_7AH12f1_592nacl/modifiedO-Na+/6.0.263sigma_0.635eps/build/xp0.05_7AA12f1_7AH12f1_592nacl_10396hoh.crd

temp=298.15
platformName=CUDA #CUDA, CPU, OpenCL
x=5.4 #in nm
y=5.4 #in nm
trajStride=1
checkStride=1
check="'$maindir/$NPTdir/checkpnt.chk'"  #enclose in "''" if not boolean  
state=False  #enclose in "''" if not boolean                              
eqSteps=5 #250000 #0.5ns
steps=5 #200000000 #400ns

echo 'getting average volume and energy'
cd $NPTdir; python ~/bin/scripts/log2txt.py -i log.txt ; cd $maindir

mkdir simNVE; cd simNVE
python $libDir/avg.py $maindir/$NPTdir/data.txt $libDir/sim_NVE_template.py -e 'Total_Energy_(kJ/mole)' -v 'Box_Volume_(nm^3)' -x $x -y $y

sed -i "s#__top__#${top}#g" sim.py
sed -i "s#__crd__#${crd}#g" sim.py
sed -i "s#__check__#${check}#g" sim.py                                                    
sed -i "s#__state__#${state}#g" sim.py                                                    
sed -i "s/__temp__/${temp}/g" sim.py
sed -i "s/__eqSteps__/${eqSteps}/g" sim.py                                                
sed -i "s/__steps__/${steps}/g" sim.py 
sed -i "s/__platformName__/${platformName}/g" sim.py
sed -i "s/__checkStride__/${checkStride}/g" sim.py
sed -i "s/__trajStride__/${trajStride}/g" sim.py
python sim.py


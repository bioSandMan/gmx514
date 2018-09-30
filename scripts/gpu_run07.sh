#!/bin/bash
#SBATCH -N 1                      # number of nodes
#SBATCH -n 8                      # number of cores
#SBATCH --mem 1024                # memory pool for all cores
module load gcc/4.9.1
echo $(hostname)
source /usr/modules/init/bash
gmx="/mnt/cmci/cmira/bin/gmx514/bin/gmx" # PBS script will need to load module
mdrun="/mnt/cmci/cmira/bin/gmx514/bin/gmx mdrun"

mdrun="$mdrun -ntomp 8 -pin on -pinoffset 16 -pinstride 1"
run=07
echo "Running "
$gmx grompp -f expanded.mdp -c methane.gro -p methane.top -o out$run.tpr -maxwarn 4
$mdrun -deffnm out$run >> aimgpu$run.log 2>&1
echo "Script complete."
exit 0

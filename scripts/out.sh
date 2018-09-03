#!/bin/bash
## Job Name
#PBS -N chrism_lambdaJobs
## Specifies number of nodes and cores per node
#PBS -l nodes=1:ppn=8
## Specifies que name; default, gpu, hughmem
#PBS -q gpu


NODEFILE="nodefile.txt"
cat $PBS_NODEFILE>$NODEFILE
cd $PBS_O_WORKDIR

# Set GROMACS environment variables
gro=/home/chrism/bin/gromacs-5.1.4/bin
gmx=$gro/gmx-514
mdrun="$gmx mdrun"

output="prod"
in="out"
log="run.log"
mdrun="$mdrun -nt 8 -pin on -pinoffset 0"

echo "Running "
nlambda=`grep coul-lambdas $in.mdp | wc -w`
((nlambda=nlambda-2))
lastgro="$in.gro"
echo "  performing alchemical sims with $nlambda lambda values"
for ((i=0; i<$nlambda; i++)); do
    echo "    running lambda column $i"
    ibase="${output}.${i}"
    sed 's/ILS/'$i'/' $in.mdp > ${ibase}.mdp
    $gmx grompp -maxwarn 3 -f $ibase -po $ibase.mdout -c $lastgro\
            -p $in -o $ibase >> $log 2>&1
    $mdrun -s $ibase -o $ibase -x $ibase -c $ibase.last -e $ibase\
            -g $ibase -cpo $ibase -dhdl $ibase.dhdl >> $log 2>&1
    lastgro="$ibase.last.gro"
done

echo "Script complete. Can remove generated files with -c option."

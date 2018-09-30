#!/bin/bash
module load gcc/6.1.0
gmx="/mnt/ceph/cmira/gmx514cpu/bin/gmx" # PBS script will need to load module
mdrun="/mnt/ceph/cmira/gmx514cpu/bin/gmx mdrun"

output="prod"
in="tripep_A2V"
log="run.log"
mdrun="$mdrun -ntmpi 1 -ntomp 8 -pin on -pinoffset 0 -pinstride 1"

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

echo "Script complete."

#!/bin/bash
# find first xvg file in the current directory
pattern="*.xvg"
files=( $pattern )
xvgfile="${files[0]}"

# find number of lambdas in xvg file
number=`grep "^@" $xvgfile | tail -1 | awk '{ print $2 }' | sed 's/[^0-9]*//g'`
n_lambdas=$(expr $number - 2)
prefix="prod"
suffix="dhdl.xvg"

out="A2VTI100ps.out"

for ((j=0; j<=n_lambdas; j++)); do
    cat $prefix.$j.$suffix | sed '/^#/ d' | sed '/^@/ d' | awk '{for (i=2;i<=3;i++){a[i]+=$i;}} END {for (i=2;i<=3;i++){printf "%.5f", a[i]/NR; printf "\t"};printf "\n"}' >> $out
done

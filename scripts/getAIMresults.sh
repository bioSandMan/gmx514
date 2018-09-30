#!/bin/bash
lambdas=11
tail_lines=115
out="aim100ps11ls0"
for i in {1..8}; do tail -$tail_lines out0$i.log | grep -A $lambdas "1  0.000  0.000" | tr -d "<<" >> $out$i.out; done

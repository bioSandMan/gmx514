#!/bin/bash
cd ~/tripep_A2V_TI/moreLambdas/run01
sbatch -p volatile --nodelist=n018 alchem.sh
cd ~/tripep_A2V_TI/moreLambdas/run02
sbatch -p volatile --nodelist=n044 alchem.sh
cd ~/tripep_A2V_TI/moreLambdas/run03
sbatch -p volatile --nodelist=n045 alchem.sh
cd ~/tripep_A2V_TI/moreLambdas/run04
sbatch -p volatile --nodelist=n046 alchem.sh
cd ~/tripep_A2V_TI/moreLambdas/run05
sbatch -p volatile --nodelist=n047 alchem.sh
cd ~/tripep_A2V_TI/moreLambdas/run06
sbatch -p volatile --nodelist=n048 alchem.sh
cd ~/tripep_A2V_TI/moreLambdas/run07
sbatch -p volatile --nodelist=n049 alchem.sh
cd ~/tripep_A2V_TI/moreLambdas/run08
sbatch -p volatile --nodelist=n050 alchem.sh

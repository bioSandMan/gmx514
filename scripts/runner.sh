#!/bin/bash
cd ~/tripep_A2V_TI/moreLambdas/run01
sbatch -p volatile --nodelist=n018 cpu_run.sh
cd ~/tripep_A2V_TI/moreLambdas/run02
sbatch -p volatile --nodelist=n044 cpu_run.sh
cd ~/tripep_A2V_TI/moreLambdas/run03
sbatch -p volatile --nodelist=n045 cpu_run.sh
cd ~/tripep_A2V_TI/moreLambdas/run04
sbatch -p volatile --nodelist=n046 cpu_run.sh
cd ~/tripep_A2V_TI/moreLambdas/run05
sbatch -p volatile --nodelist=n047 cpu_run.sh
cd ~/tripep_A2V_TI/moreLambdas/run06
sbatch -p volatile --nodelist=n048 cpu_run.sh
cd ~/tripep_A2V_TI/moreLambdas/run07
sbatch -p volatile --nodelist=n049 cpu_run.sh
cd ~/tripep_A2V_TI/moreLambdas/run08
sbatch -p volatile --nodelist=n050 cpu_run.sh

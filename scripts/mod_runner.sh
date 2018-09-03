#!/bin/bash
cd ~/tripep_A2V_TI/run01
sbatch -p volatile --nodelist=n018 mod_cpu_run.sh
cd ~/tripep_A2V_TI/run02
sbatch -p volatile --nodelist=n044 mod_cpu_run.sh
cd ~/tripep_A2V_TI/run03
sbatch -p volatile --nodelist=n045 mod_cpu_run.sh
cd ~/tripep_A2V_TI/run04
sbatch -p volatile --nodelist=n046 mod_cpu_run.sh
cd ~/tripep_A2V_TI/run05
sbatch -p volatile --nodelist=n047 mod_cpu_run.sh
cd ~/tripep_A2V_TI/run06
sbatch -p volatile --nodelist=n048 mod_cpu_run.sh
cd ~/tripep_A2V_TI/run07
sbatch -p volatile --nodelist=n049 mod_cpu_run.sh
cd ~/tripep_A2V_TI/run08
sbatch -p volatile --nodelist=n050 mod_cpu_run.sh

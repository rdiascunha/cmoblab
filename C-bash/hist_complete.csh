#!/bin/csh

#cat ../cMD/cluster_04/rg_0-3us.dat ../cMD/cluster_06/rg_0-4us.dat ../cMD/cluster_08/rg_0-1us.dat > cMD_rg_8us.dat
# From Rg cpptraj files to create a Histogram (delete the first line and the frames column)
#python histogram.py -f cMD_rg_complete.dat -o cMD_rg_hist.dat --min 10 --max 70 -b 240
#python histogram.py -f GaMD_rg_complete.dat -o GaMD_rg_hist.dat --min 10 --max 70 -b 240
#python histogram.py -f SIRAH_rg_complete.dat -o SIRAH_rg_hist.dat --min 10 --max 70 -b 240
# Normalize by area the Histogram
#python norm_histogram.py -i cMD_rg_hist.dat -o norm_cMD_rg_hist.dat
#python norm_histogram.py -i GaMD_rg_hist.dat -o norm_GaMD_rg_hist.dat
#python norm_histogram.py -i SIRAH_rg_hist.dat -o norm_SIRAH_rg_hist.dat
#python norm_histogram.py -i ../cluster_04_12LJ/RG/HIST/SIRAH_rg_complete_hist.dat -o norm_SIRAH_12LJ_rg_hist.dat
# Plot the normalized histogram
#python plot_histogram.py -f norm_cMD_rg_hist.dat norm_GaMD_rg_hist.dat norm_SIRAH_rg_hist.dat norm_SIRAH_12LJ_rg_hist.dat -l cMD GaMD SIRAH SIRAH-scaled 
python plot_histogram.py -f norm_SIRAH_rg_hist.dat norm_SIRAH_12LJ_rg_hist.dat -l SIRAH SIRAH-scaled 

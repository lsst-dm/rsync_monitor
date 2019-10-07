#!/bin/bash
if [ "$#" -ne 1 ]; then 
    echo "auto-check-rsync.bash DIR"
    echo "DIR is the directory where monitor_rsync.py lives"
    echo "This script is desinged to be a one-liner that inserts into your crontab."
    exit 1
fi 
setup lsst_distrib
today=`date +"%Y-%m-%d"`
threedaysago=`date -d '-3day' +"%Y-%m-%d"`
$1/monitor_rsync.py -c $1/config/auxTel_L1Archiver_monitor.yaml $threedaysago $today
$1/monitor_rsync.py -c $1/config/auxTel_L1Archiver_monitor.yaml $threedaysago $today
$1/monitor_rsync.py -c $1/config/auxTel_L1Archiver_monitor.yaml $threedaysago $today

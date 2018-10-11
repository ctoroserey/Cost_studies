#!bin/bash/

# Move data from lab_tasks to where I'm doing analyses
dataDir='/Users/ctoro/git_clones/lab_tasks/Claudio/Replication_2018'
analDir='/Users/ctoro/git_clones/R/MA681/MA681_Final/data'

# Cognitive effort
cp ${dataDir}/cog/data/*_log.csv ${analDir}/cog
rm -f ${analDir}/cog/attention*

# Easy effort
cp ${dataDir}/pheasycal/data/*_log.csv ${analDir}/easy
rm -f ${analDir}/easy/attention*
rm -f ${analDir}/easy/*_calibration_log.csv

# Physical effort
cp ${dataDir}/phys/data/*_log.csv ${analDir}/phys
rm -f ${analDir}/phys/attention*
rm -f ${analDir}/phys/*_calibration_log.csv

# Passive delay
cp ${dataDir}/wait/data/*_log.csv ${analDir}/wait
rm -f ${analDir}/wait/attention*


#$ -N rootsolve
#$ -cwd
#$ -S /bin/bash
#$ -o job_reports/
#$ -e job_reports/
#$ -t 1-1000
#$ -l h_vmem=5G

if [ -n "${1}" ]; then
        SGE_TASK_ID=${1}
fi
i=${SGE_TASK_ID}

outdir="/clusterdata/uqgheman/jummy/qualitative-modeling/results/"

R --no-save --args ${outdir} ${i} < /clusterdata/uqgheman/jummy/qualitative-modeling/R/run_solve_simulations.R

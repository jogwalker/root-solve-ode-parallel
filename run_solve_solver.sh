
#$ -N runODE
#$ -cwd
#$ -S /bin/bash
#$ -o job_reports/
#$ -e job_reports/
#$ -t 1-153
#$ -l h_vmem=2G

if [ -n "${1}" ]; then
        SGE_TASK_ID=${1}
fi
i=${SGE_TASK_ID}

outdir="/clusterdata/uqgheman/jummy/root-solve-ode-parallel/results/runODE/"

R --no-save --args ${outdir} ${i} < /clusterdata/uqgheman/jummy/root-solve-ode-parallel/R/run_solve_solver.R

HOME=/Users/zhaomingyin/gitlocal/optkit/
RESULT_HOME=${HOME}data/results/disttime/
PY_HOME=${HOME}python/disdcj/
PY_EXE=process_dcjdis_time.py
R_HOME=${HOME}R_code/disdcj/
PLOT_HOME=${HOME}data/plots/disdcj/


#process sim data
echo "processing simulated data"
ARG1=${RESULT_HOME}result_optkit_sim_exem
ARG2=${RESULT_HOME}result_optkit_sim_matc
ARG3=${RESULT_HOME}result_gredo_sim
ARG4=${RESULT_HOME}result_sim.csv
ARG5=sim
python ${PY_HOME}${PY_EXE} ${ARG1} ${ARG2} ${ARG3} ${ARG4} ${ARG5}

#process real data
echo "processing real data"
ARG1=${RESULT_HOME}result_optkit_real_exem
ARG2=${RESULT_HOME}result_optkit_real_matc
ARG3=${RESULT_HOME}result_gredo_real
ARG4=${RESULT_HOME}result_real.csv
ARG5=real
python ${PY_HOME}${PY_EXE} ${ARG1} ${ARG2} ${ARG3} ${ARG4} ${ARG5}

###plot using R goes here

#plot time result
Rscript ${R_HOME}plot_distdcj_time.R ${RESULT_HOME}result_sim.csv  ${PLOT_HOME}sim_time.png
Rscript ${R_HOME}plot_distdcj_time.R ${RESULT_HOME}result_real.csv ${PLOT_HOME}real_time.png

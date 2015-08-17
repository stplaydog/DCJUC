#source /opt/intel/Compiler/11.1/059/bin/iccvars.sh intel64
source /opt/intel/Compiler/11.1/064/bin/intel64/iccvars_intel64.csh
source  /opt/intel/impi/3.2.2.006/bin64/mpivars.sh
unset LIB

num_nodes=8
NNODES=$(sort -u $PBS_NODEFILE |grep -c ^)
RANKS=$(grep -c ^ $PBS_NODEFILE)
mpdboot -r ssh -n ${num_nodes} -f $PBS_NODEFILE
chmod 600 $HOME/.mpd.conf



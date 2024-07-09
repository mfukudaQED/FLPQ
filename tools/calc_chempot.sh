#!/bin/bash
#SBATCH -p i8cpu
#SBATCH -N 4
#SBATCH -t 0:30:0

#cd ${SLURM_SUBMIT_DIR}
NUM_CPU_CORES_PER_NODE=128 # ohtaka:128 cores/node

#################################
### Parallel calculation info ###
#################################

TOTAL_NUM_NODES=${SLURM_JOB_NUM_NODES}
TOTAL_NUM_CORES=$((${TOTAL_NUM_NODES}*${NUM_CPU_CORES_PER_NODE}))

NUM_NODES=4
NUM_PROCS_PER_NODES=32
NUM_THREADS=4
MEMORY_PER_CPU_CORE=1840

NUM_PROCS=$((${NUM_PROCS_PER_NODES}*${NUM_NODES}))
TOTAL_CORES_PER_TASK=$((${NUM_PROCS}*${NUM_THREADS}))
NUM_PARALLEL_JOBS=$((${TOTAL_NUM_CORES}/${TOTAL_CORES_PER_TASK}))
export OMP_NUM_THREADS=${NUM_THREADS}

RUN_O="srun --exclusive --mem-per-cpu=${MEMORY_PER_CPU_CORE} -n ${NUM_PROCS} -c ${NUM_THREADS} -N ${NUM_NODES}"
RUN_F="srun --exclusive --mem-per-cpu=${MEMORY_PER_CPU_CORE} -n 64 -c 4 -N 2"


################################
### Apps, input, output info ###
################################

PATH_OPENMX="/home/i0017/i001700/program/flpq/openmx_flpq/OpenMX_for_FLPQ/source"
PATH_FLPQ="/home/i0017/i001700/program/flpq/flpq/FLPQ/source"
EXE_O="${PATH_OPENMX}/openmx"
EXE_F="${PATH_FLPQ}/flpq"
INP_O=INPUT.dat
INP_F=flpq.inp
LOG_O=RESULT.std
LOG_F=flpq.log
SYSTEM_NAME="INPUT"

FLAG_LOG_JOBINOFO=true
LOG_JOBINFO="jobinfo_${SLURM_JOBID}.log"

FLAG_GNU_PARALLEL=false
LOG_GNU_PARALLEL="runtask.log"
FILENAMELIST="filenamelist.txt"

#############################
### Setup of modules etc. ###
#############################

ulimit -s unlimited
module purge
#module load oneapi_compiler/2023.0.0 oneapi_mkl/2023.0.0 openmpi/4.1.5-oneapi-2023.0.0-classic
module load oneapi_compiler/2023.0.0 oneapi_mkl/2023.0.0 oneapi_mpi/2023.0.0

#export MKL_DEBUG_CPU_TYPE=5 # This option is not available in new MKL.
export I_MPI_COLL_EXTERNAL=no                                                                                                                                   
export I_MPI_FABRICS=shm:ofi  
#export FI_PROVIDER=psm3
#unset I_MPI_PMI_LIBRARY
export UCX_TLS='self,sm,ud'

### For Spglib ###
export SPGLIB_ROOT=/home/i0017/i001700/program/alamode/etc/spglib
export LD_LIBRARY_PATH=$SPGLIB_ROOT/lib64:$LD_LIBRARY_PATH

### For HDF5 ###
export PATH=$PATH:/home/local/ap/hdf5/1.10.6
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/local/ap/hdf5/1.10.6/lib

### For gnuplot ###                                                                                                                                             
PATH_GNUPLOT="~/program/gnuplot-5.4.3/installed_files/bin/gnuplot" 

################
### Workflow ###
################

NAMEDIR=$PWD

task+=("")
#task+=("cp ${NAMEDIR}/${INP} ${NAMEDIR}/{}")
#task+=("$RUN $EXE $INP -nt ${OMP_NUM_THREADS} > $LOG")

task_all="echo"                                                                                                                                                     
for itask in "${task[@]}"                                                                                                                                           
do                                                                                                                                                                  
    task_all="$task_all && $itask"                                                                                                                                  
done  


#task_all="$(IFS=" "; echo "${task[*]}")"
#echo $task_all



################################
## Write info to LOG_JOBINFO ###
################################

if ${FLAG_LOG_JOBINOFO}; then
    JOB_SCRIPT_NAME=$(basename "$0")
    date > ${LOG_JOBINFO}
    cat << EOF >> ${LOG_JOBINFO}

#########################
### Job resource info ###
#########################
JOB_SUB_DIR=${SLURM_SUBMIT_DIR}
SLURM_JOB_NAME=${SLURM_JOB_NAME}

SLURM_JOB_ID=${SLURM_JOB_ID}
SLURM_JOB_PARTITION=${SLURM_JOB_PARTITION}

TOTAL_NUM_NODES=${TOTAL_NUM_NODES}
TOTAL_NUM_CORES=${TOTAL_NUM_CORES}

SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST}
#########################


#######################################
### Parallel calc info for each job ###
#######################################
NUM_PARALLEL_JOBS=${NUM_PARALLEL_JOBS}
NUM_NODES=${NUM_NODES}
NUM_PROCESS_PER_NODE=${NUM_PROCS_PER_NODES}
NUM_PROCS=${NUM_PROCS}
NUM_THREADS=${NUM_THREADS}
#######################################


###################
### Module info ###
###################

EOF

    module list &>> ${LOG_JOBINFO}
    #module list -t | tail -n +2 | parallel -j 1 "module show {}"    >> ${LOG_JOBINFO}
    #
    cat << EOF >> ${LOG_JOBINFO}


#####################
### Workflow info ###
#####################
${task_all}

EOF
fi


#################################################
### Handling for massive parallel calculation ###
#################################################

# The filehandle limit of GNU parallel is around 250 jobs.
# To avoid the limitation, you just spawn more GNU parallels.
if ${FLAG_GNU_PARALLEL}; then
    if [ ${NUM_PARALLEL_JOBS} -le 250 ]; then
        FLAG_LARGE_PARALLEL=false
        PARALLEL="parallel --delay 0.1 -j ${NUM_PARALLEL_JOBS} --joblog ${LOG_GNU_PARALLEL} --resume"
        #PARALLEL="parallel --delay 0.1 -j ${NUM_PARALLEL_JOBS} --joblog ${LOG_GNU_PARALLEL} --resume-failed"
    else
        FLAG_LARGE_PARALLEL=true
        SIZE_SPAWN2=250
        SIZE_SPAWN1=$(( $((1+${NUM_PARALLEL_JOBS})) /${SIZE_SPAWN2} ))
        PARALLEL1="parallel --pipe -L ${SIZE_SPAWN2} --round-robin -j ${SIZE_SPAWN1} cat"
        PARALLEL2="parallel -a - -j ${SIZE_SPAWN2} --joblog ${LOG_GNU_PARALLEL} --resume"
    fi
fi


############################
### Perform calculations ###
############################

#if ${FLAG_GNU_PARALLEL}; then
#    if ${FLAG_LARGE_PARALLEL}; then
#        cat ${FILENAMELIST} | ${PARALLEL1} | ${PARALLEL2} "${task_all}" &
#    else
#        cat ${FILENAMELIST} | ${PARALLEL} "${task_all}" &
#    fi
#else
#    #eval ${task_all} &
#fi

### For flpq ### 
cd ${NAMEDIR}
mkdir ${NAMEDIR}/flpq
cp ${INP_O}\# ${NAMEDIR}/flpq/${INP_O}
cd ${NAMEDIR}/flpq

cat << 'EOF' >> ${INP_O}
DM.export   on
EOF

$RUN_O $EXE_O $INP_O -nt ${OMP_NUM_THREADS} > $LOG_O

cd ${NAMEDIR}
cp ${INP_F} ${NAMEDIR}/flpq/${INP_F}
cd ${NAMEDIR}/flpq
$RUN_F $EXE_F $INP_F > $LOG_F

### End flpq ###

### For flpq_ew ###
cd ${NAMEDIR}
mkdir ${NAMEDIR}/flpq_ew
cp ${INP_O}\# ${NAMEDIR}/flpq_ew/${INP_O}
cd ${NAMEDIR}/flpq_ew

cat << 'EOF' >> ${INP_O}
DM.export   window
DM.specify.energy.range   on
DM.energy.range  -2.0  0.0
EOF

$RUN_O $EXE_O $INP_O -nt ${OMP_NUM_THREADS} > $LOG_O

cd ${NAMEDIR}
cp ${INP_F} ${NAMEDIR}/flpq_ew/${INP_F}
cd ${NAMEDIR}/flpq_ew
$RUN_F $EXE_F $INP_F > $LOG_F

wait


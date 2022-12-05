#!/bin/bash

## SLURM options
# Options to sbatch start with '#SBATCH'. To disable an option, change
# the prefix, e.g. '#xSBATCH'

# Time limit after which the job will be killed. The default time limit
# is 60 minutes. Specified as HH:MM:SS or D-HH:MM
#SBATCH --time=00:05:00

# number of CPUs to utilize (this job is one "task")
#SBATCH --cpus-per-task=4

# Memory per node. Job will crash if this limit is exceeded.
# Default is 1000M per allocated core. Use values that will permit multiple jobs
# to run simultaneously when possible, e.g. a memory limit of 6000M will allow
# 4 jobs to run on a node with 24000M of RAM, while specifying 6G (=6144M) would
# only allow 3 jobs to run.
#SBATCH --mem=4000M

# Prevent multithreading of the cores
#SBATCH --hint=nomultithread

# number of nodes to utilize. Can be set to >1 for MPI jobs. Otherwise, 
# leave this set to 1 to make sure all allocated CPUs are on the same node.
#SBATCH -N 1

# Partition (queue) for the job
# - "normal" has a time limit of 30 days. Per-user resource limits apply.
# - "debug" has a time limit of 1 hour, and low resource limits,
#   but will pre-empt other jobs (memory permitting).
# - "idle" has no time limit or per-user limit, but jobs can be stopped
# (requeued) by jobs in the normal or debug queues.
#SBATCH --partition=normal
# or, for short, '-p normal'

# Infiniband. If a job requires Infiniband (either for MPI
# communication, or for faster access to d01-d12, the following option
# should be specified. Note that Infiniband is not available on
# c025-c034 (the 24-core nodes)
#xSBATCH --constraint=IB

# Output redirection. Default output file is 'slurm-%j.out' where %j is the
# job number. Setting "--output" alone will redirect both stdout and stderr
# to the same file.
#xSBATCH --output=somefile-%j.txt
#xSBATCH --error=otherfile-%j.txt

# Send email when the job begins and ends
#xSBATCH --mail-user=your.address@example.com
#xSBATCH --mail-type=BEGIN,END
# other options: FAIL, REQUEUE, TIME_LIMIT

# Set a name for the job. This will appear in the output of 'squeue'. The
# default is the name of the job script. This option is probably most useful
# as a command-line argument to sbatch.
#xSBATCH -J some_name

# To submit the job, run:
#
#     sbatch <scriptname.sh> <runfile.jl>
#

julia --project=../. $1
# Report additional information if job was killed because of memory limit
oom_check $?

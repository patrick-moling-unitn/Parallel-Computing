# Parallel Computing project - First Deliverable - SpmV

### This repository contains
+ **main.cpp** which needs to be compiled in order to create the main one
+ **mmio.c** and **mmio.h** that are used during the compilation for adding matrix reading support
+ **schedulable_job.pbs** which is an example of a job that can be scheduled to the cluster
+ **685_bus.mtx**, **Ragusa16.mtx**, **bcsstk19.mtx**, **bcsstm24.mtx**, **onetone2.mtx** that are the matrixes used for this project
+ **RawOutputs.zip** containing all the logs given by the cluster when the .pbs files in **ExecutedClusterJobs.zip** where executed
  
### How to reproduce the results
Download the files from git and copy the folder into the cluster
> scp -r ./YourGitFolderName/ yourname.yoursurname@hpc.unitn.it:~/YourDestinationFolder

Inside your cluster, in YourDestinationFolder, compile the **main.cpp**
> g++ -std=c++11 -fopenmp -o main main.cpp mmio.c

In case of errors due to Windows' *^M* execute the following commands and repeat the compilation attempt
> dos2unix schedulable_job.pbs
> 
> dos2unix main.cpp

Run the job you'd like to
> qsub schedulable_job.pbs

Check if the job has finished
> qstat -u yourname.yoursurname

Verify the output once done
> cat execution.out

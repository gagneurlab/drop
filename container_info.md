# Using the containerised version of DROP
The prebuilt container is available here: https://quay.io/repository/sarahbeecroft9/drop

In the example, let's use the Singularity container engine.

Getting the container from the repo onto your machine:
`singularity pull docker pull quay.io/sarahbeecroft9/drop`
This will create a file called drop_latest.sif

## some examples
Here, the my_drop_analysis_dir is the directory you wish to read/write your drop results into/out of. You can also specify the full path as opposed to the relative path used in these examples. It is important to use the `-B my_drop_analysis_dir:/home/ubuntu` flag, which will bind-mount your local analysis directory to the home directory in the container. 

### Example 1A: Run `drop demo` in `my_drop_analysis_dir`
`singularity exec -B my_drop_analysis_dir:/home/ubuntu drop.sif drop demo`

### Example 1B: Run snakemake pipeline in `my_drop_analysis_dir`
`singularity exec -B my_drop_analysis_dir:/home/ubuntu drop.sif snakemake -c 10`

### Example 2A: Run `drop init` in `my_drop_analysis_dir`
`singularity exec -B my_drop_analysis_dir:/home/ubuntu drop.sif drop init`

### Example 2B: Run snakemake pipeline in `my_drop_analysis_dir`
`singularity exec -B my_drop_analysis_dir:/home/ubuntu drop.sif snakemake -c 10`

Troubleshooting
===============

In case you have any issues, please check the previous issues on `git <https://github.com/gagneurlab/drop>`_ and the Troubleshoot section of the DROP paper. If you have a new issue, please create a new one on git.

Common errors
-------------

The ``MAE:mae_allelicCounts`` step is susceptible to fail if:

1. the chromosomes styles of the reference genome and the BAM files do not match

Solution: Identify the chromosomes style of the BAM file. Obtain an appropriate reference genome file and specify it in the config file.

2. the BAM file does not have the correct ``Read Groups`` documentation both the header and reads. You can often identify if the BAM file has any problems by using the command ``gatk ValidateSamFile -I path/to/bam_file.bam``

Solution:
To fix this is often dependent on the individual case, but some combination of the following tools is quite helpful:  

* `samtools reheader <http://www.htslib.org/doc/samtools-reheader.html>`_
* `gatk AddOrReplaceReadGroups <https://gatk.broadinstitute.org/hc/en-us/articles/5358911906459-AddOrReplaceReadGroups-Picard->`_

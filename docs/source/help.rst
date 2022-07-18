Troubleshooting
===============

In case you have any issues, please open an issue on `git <https://github.com/gagneurlab/drop>`_.

A common problem is that during the ``MAE:mae_allelicCounts`` step, if the BAM file does not have the correct ``Read Groups`` documentation both the header and reads.  
You can often identify if the BAM file is the problem by using the command ``gatk ValidateSamFile -I path/to/bam_file.bam``

To fix this is often dependent on the individual case, but some combination of the following tools is quite helpful:  

* `samtools reheader <http://www.htslib.org/doc/samtools-reheader.html>`_
* `gatk AddOrReplaceReadGroups <https://gatk.broadinstitute.org/hc/en-us/articles/5358911906459-AddOrReplaceReadGroups-Picard->`_

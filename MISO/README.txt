MISO needs python 2.7, so it's better to use a new conda environment. Create one using:
conda create --name my_new_env

Then, activate it using
conda activate my_new_env

Then, install MISO using
conda install misopy

Copy the BAM files of interest in a folder called bam-data. Also add the corresponding index (.bai) files.

Edit the ./settings/new_splice_site_events.gff by giving appropriate genomic coordinates.
Edit the ./settings/plot_timmdc1_new_exon.conf file with the samples and files of interest.

Edit the last line of the run_miso.sh script to indicate the event and config file.

Execute the run_miso.sh script. Takes around 10s.
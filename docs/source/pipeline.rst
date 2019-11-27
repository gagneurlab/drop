Executing the Pipeline
======================

DROP is executed in the same way as `Snakemake
<https://snakemake.readthedocs.io/en/stable/executing/cli.html>`_. 

Open a terminal in your project repository. Execute 

``snakemake -n``. 


That will perform a *dry-run*, which means it will display all the steps (or rules) that need to be executed. To also display the reason why those rules need to be exeucted, run ``snakemake -nr``. Finally, a simplified dry-run can be achieved by executing ``snakemake -nq``. 

Calling 

``snakemake``

without any parameters will execute the whole workflow. DROP's steps are computationally heavy, therefore it is a good idea to run them in parallel. Snakemake automatically determines the steps that can be parallelized. The user simply needs to specify the maximum number of cores that Snakemake can take, e.g.:

``snakemake --cores 10``

While running, Snakemake *locks* the directory. If, for a whatever reason, the pipeline was interrupted, the directory might be kept locked. Therefore, run ``snakemake --unlock`` to unlock it.
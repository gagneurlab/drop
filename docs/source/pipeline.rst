Executing the Pipeline
======================

DROP is `Snakemake <https://snakemake.readthedocs.io/en/stable/executing/cli.html>`_ pipeline, so it is called with the ``snakemake`` command.

Dry run
-------

Open a terminal in your project repository. Execute 

.. code-block:: bash
    
    snakemake -n 

This will perform a *dry-run*, which means it will display all the steps (or rules) that need to be executed. To also display the reason why those rules need to be exeucted, run 

.. code-block:: bash

    snakemake -nr

Finally, a simplified dry-run can be achieved by executing

.. code-block:: bash

    snakemake -nq
    
Calling ``snakemake`` without any parameters will execute the whole workflow. 


Parallelizing jobs
------------------

DROP's steps are computationally heavy, therefore it is a good idea to run them in parallel. Snakemake automatically determines the steps that can be parallelized. The user simply needs to specify the maximum number of cores that Snakemake can take, e.g. for 10 cores:

.. code-block:: bash

    snakemake --cores 10

If the ``--cores`` flag is not specified, snakemake will use a single core by default.


Executing parts of the pipeline
-------------------------------

Every single module can be called independently.

.. code-block:: bash

    snakemake <subworkflow>
    
========================  =======================================================================
<subworkflow>                Description                                                       
========================  =======================================================================
``aberrantExpression``     Aberrant expression pipeline
``aberrantSplicig``        Aberrant splicing pipeline
``mae``                    Monoalleic expression pipeline
========================  =======================================================================


Unlocking the pipeline
----------------------

While running, Snakemake *locks* the directory. If, for a whatever reason, the pipeline was interrupted, the directory might be kept locked. Therefore, run 

.. code-block:: bash

    snakemake --unlock

to unlock it.


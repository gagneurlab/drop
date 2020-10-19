Pipeline Commands
=================

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


Executing subworkflows
----------------------

Every single module can be called independently.

.. code-block:: bash

    snakemake <subworkflow>
    
========================  =======================================================================
Subworkflow                Description                                                       
========================  =======================================================================
``aberrantExpression``     Aberrant expression pipeline
``aberrantSplicing``       Aberrant splicing pipeline
``mae``                    Monoalleic expression pipeline
========================  =======================================================================

An example for calling the aberrant expression pipeline with 10 cores would be 

.. code-block:: bash

    snakemake aberrantExpression --cores 10

Rerunning the pipeline
----------------------

When DROP is updated or jobs fail, the following commands can be used to rerun and troubleshoot.


Unlocking the pipeline
++++++++++++++++++++++

While running, Snakemake *locks* the directory. If, for a whatever reason, the pipeline was interrupted, the directory might be kept locked. Therefore, call 

.. code-block:: bash

    snakemake unlock

to unlock it. This will call snakemake's ``unlock`` command for every module

.. _dropUpdate:

Updating DROP
+++++++++++++
Every time a project is initialized, a temporary folder ``.drop`` will be created in the project folder.
If a new version of drop is installed, the ``.drop`` folder has to be updated for each project that has been
initialized using an older version.
To do this run:

.. code-block:: bash

    drop update

Skipping recomputation of files
+++++++++++++++++++++++++++++++

If snakemake is interrupted and restarted, it will continue with the last unsuccessful job in the job graph. If a script is updated with minor change, e.g. when calling ``drop update``, all jobs of the modified script and its downstream steps will be rerun. However, in some cases one might want to keep the intermediate files instead and continue with the missing files. In order to do so, first execute

.. code-block:: bash
   
   snakemake <rule> --touch

for whichever rule or module you want to continue the computation. The ``--touch`` command touches all output files required by the pipeline that have already been computed. Omitting the rule will lead to accessing the complete pipeline. Afterwards, use 

.. code-block:: bash

    snakemake unlock
    
to unlock the submodules, so that the jobs that need to be computed can be identified.


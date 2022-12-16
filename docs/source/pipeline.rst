Pipeline Commands
=================

DROP is a `Snakemake <https://snakemake.readthedocs.io/en/stable/executing/cli.html>`_ pipeline, so it is called with the ``snakemake`` command.

Dry run
-------

Open a terminal in your project repository. Execute

.. code-block:: bash

    snakemake --cores 1 -n

This will perform a *dry-run*, which means it will display all the steps (or rules) that need to be executed. To also display the reason why those rules need to be executed, run

.. code-block:: bash

    snakemake --cores 1 -nr

A simplified dry-run can be achieved using the ``-q`` parameter.

.. code-block:: bash

    snakemake --cores 1 -nq

Calling ``snakemake --cores 1`` without any additional parameters will execute the whole workflow. Snakemake requires you to designate the number of cores.


Parallelizing jobs
------------------

DROP's steps are computationally heavy, therefore it is highly recommended to run them in parallel. Snakemake automatically determines the steps that can be parallelized. The user simply needs to specify the maximum number of cores that Snakemake can take, e.g. for 10 cores:

.. code-block:: bash

    snakemake --cores 10


Executing subworkflows
----------------------

Every module can be called independently.

.. code-block:: bash

    snakemake <subworkflow>

========================  =======================================================================
Subworkflow                Description
========================  =======================================================================
``aberrantExpression``     Aberrant expression pipeline
``aberrantSplicing``       Aberrant splicing pipeline
``mae``                    Monoallelic expression pipeline
``rnaVariantCalling``      RNA Variant Calling pipeline
========================  =======================================================================

For example, to run the aberrant expression pipeline with 10 cores, execute the following

.. code-block:: bash

    snakemake aberrantExpression --cores 10

Rerunning the pipeline
----------------------

When DROP is updated or jobs fail, the following commands can be used to rerun and troubleshoot.


Unlocking the pipeline
++++++++++++++++++++++

While running, Snakemake *locks* the directory. If, for a whatever reason, the pipeline was interrupted, the directory might be kept locked. If this is the case, call

.. code-block:: bash

    snakemake unlock

.. _dropUpdate:

Updating DROP
+++++++++++++
The developers of DROP are active in making DROP a better tool. As a result there are often bug fixes
or improvements that are implemented and released in new versions. You can check them out in the *What's new* section of the
`README. <https://github.com/gagneurlab/drop#whats-new>`_ 

When updating DROP we recommend using the conda/mamba functions to maintain any dependencies that could be related.

.. code-block:: bash

    mamba update drop

If you were working with a pip installation of DROP then you would need to reinstall using pip directly from github.

.. code-block:: bash

    pip install git+https://github.com/gagneurlab/drop.git

Once you have successfully bumped the DROP version to the latest, you will still need to update your project folder.
`drop update` will reset the local project's `Scripts/` directory to match the installed version,
so be sure to save any additional scripts or analyses in another location.

To complete your update, you must run the following to get your local directory to match the version:

.. code-block:: bash

    drop update

Skipping recomputation of files
+++++++++++++++++++++++++++++++

If the pipeline is interrupted and restarted, it will continue with the last unsuccessful job in the job graph. If a script is updated with minor change, e.g. when calling ``drop update``, all the jobs of the modified script and its downstream steps will be rerun. However, in some cases one might want to keep the intermediate files instead and continue with the missing files. In order to do so, first execute

.. code-block:: bash

   snakemake <rule> --touch

for whichever rule or module you want to continue the computation. The ``--touch`` command touches all output files required by the pipeline that have already been computed. Omitting the rule will lead to accessing the complete pipeline. Afterwards, run

.. code-block:: bash

    snakemake unlock

Overall, we recommend reading the snakemake documentation for further fine-tuning of the execution.

Installation
============

DROP is available on `bioconda <https://anaconda.org/bioconda/drop>`_ for python 3.6 and above.
We recommend using a dedicated conda environment.

.. code-block:: bash

    # create environment
    conda create -n drop_env python=3.6
    conda activate drop_env

    # install drop
    conda install -c bioconda drop

Installation time: ~ 10min

Test whether the pipeline runs through by setting up the demo dataset in an empty directory (e.g. ``~/drop_demo``).

.. code-block:: bash

    mkdir ~/drop_demo
    cd ~/drop_demo

    # demo will download the necessary data and pipeline files
    drop demo

The pipeline can be run using ``snakemake`` commands

.. code-block:: bash

    snakemake -n # dryrun
    snakemake

Initialize a project
++++++++++++++++++++
The demo project can be modified to be used for a new project.
Alternatively, a new DROP project can be set up using ``drop init``.

.. code-block:: bash

    cd <path-to-project>
    drop init

This will create an empty ``config.yaml`` file that needs to be filled according to the project data.
You also need to prepare a sample annotation file.
Go to :ref:`prepare` for more details.

.. _otherversions:

Other DROP versions
-------------------

The developer version of DROP can be found in the `repository <https://github.com/gagneurlab/drop>`_ under the branch
``dev``.
Make sure that the [dependencies](#dependencies) are installed.

.. code-block:: bash

    # activate your python environment if you are using one, e.g. drop_env
    conda activate drop_env

Then install DROP from github using ``pip``.
For this recursively clone the repository with all its submodules and then install from directory.

.. code-block:: bash

    git clone -b dev https://github.com/gagneurlab/drop.git --recurse-submodules
    pip install ./drop

Alternatively, you can also install it directly without cloning

.. code-block:: bash

    pip install git+https://github.com/gagneurlab/drop.git@dev

Dependencies
------------
The easiest way to ensure that all dependencies are installed is to install the
`bioconda package <https://anaconda.org/bioconda/drop>`_ into a conda environment.
Other versions of drop can be installed after the bioconda package has been installed.

Installation without conda
++++++++++++++++++++++++++
Alternatively, DROP can be installed without ``conda``.
In this case the following dependencies must be met:

+ python >= 3.6
     + pip >= 19.1
+ `samtools <https://www.htslib.org/download/>`_ >= 1.7
+ `bcftools <https://github.com/samtools/bcftools>`_ >= 1.7
+ `tabix <https://www.htslib.org/download/>`_
+ `GATK <https://software.broadinstitute.org/gatk/>`_
+ `graphviz <https://www.graphviz.org/>`_
+ `pandoc <https://pandoc.org/>`_
+ `R <https://www.r-project.org/>`_ >= 3.5 and corresponding `bioconductor <https://bioconductor.org/install/>`_ version

If you are using an already existing R installation, make sure that the R and ``bioconductor`` versions match.
Otherwise, use the newest versions of R and bioconductor.
The necessary R packages will be installed with the first pipeline call.
As this is a lengthy process, it might be desirable to install them in advance, if a local copy of the repository exists.

.. code-block:: bash

    # optional
    Rscript <path-to-drop-repo>/drop/installRPackages.R drop/requirementsR.txt


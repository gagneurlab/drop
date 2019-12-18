Installation
============

Here is the `repository <https://github.com/gagneurlab/drop>`_.

Dependencies
------------
Before installing DROP, make sure that the below listed dependencies are present and runnable.
Programming languages:

+ python >= 3.6.7
     + pip >= 19.1
     + we recommend using a virtual environment e.g. anaconda
+ R >= 3.5 (https://www.r-project.org/)


R packages
++++++++++

Bioconductor and base R packages need to be installed. The packages are listed in ``drop/requirementsR.txt``. A script for installing these packages is provided. From the repository root just execute:

.. code-block:: bash
    
    Rscript drop/installRPackages.R drop/requirements.R

    
Other packages
++++++++++++++

+ samtools >= 1.7 (https://www.htslib.org/download/)
+ bcftools (newest) (https://github.com/samtools/bcftools)
+ tabix (https://www.htslib.org/download/)
+ GATK (https://software.broadinstitute.org/gatk/)
+ graphviz (https://www.graphviz.org/)
+ pandoc (https://pandoc.org/)


Install DROP
------------

Make sure that all of the above listed [dependencies](#dependencies) are installed.
Then install DROP from github using ``pip``. For this you can either recursively clone the repository with all its submodules and then install from directory.

.. code-block:: bash

    git clone https://github.com/gagneurlab/drop.git --recurse-submodules
    # activate your python environment if you are using one
    # conda activate drop_env
    cd drop
    pip install .


Alternatively, you can also install it directly without cloning

.. code-block:: bash
    
    pip install git+https://github.com/gagneurlab/drop.git

Installation time (including all dependencies): ~ 1h


Initialize a project
++++++++++++++++++++

DROP projects are initialized in a separate directory dedictated to the analysis project. Calling the initialization command creates the necessary files.

.. code-block:: bash
    
    cd <project/path>
    drop init

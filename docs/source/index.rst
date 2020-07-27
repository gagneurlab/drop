DROP - Detection of RNA Outliers Pipeline
==========================================

DROP is intended to help researchers use RNA-Seq data in order to detect genes with aberrant expression,
aberrant splicing and mono-allelic expression. It consists of three independent modules for each of those strategies.
After installing DROP, the user needs to fill in the config file and sample annotation table (:ref:`prepare`).
Then, DROP can be executed in multiple ways (:ref:`pipeline`).

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   prepare
   pipeline
   license
   help

Quickstart
-----------

DROP is available on `bioconda <https://anaconda.org/bioconda/drop>`_ for python 3.6 and above.
We recommend using a dedicated conda environment.

.. code-block:: bash

    conda install -c bioconda drop

Initialize project

.. code-block:: bash

    cd <path-to-project>
    drop demo

Call the pipeline

.. code-block:: bash

    snakemake


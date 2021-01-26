Installation
============

DROP is available on `bioconda <https://anaconda.org/bioconda/drop>`_ .
In case the conda channel priority is set to ``strict``, it should be reset to ``flexible``:

.. code-block::

    conda config --set channel_priority true

We recommend using a dedicated conda environment (here: ``drop_env``) for installing drop.

.. code-block:: bash

    conda create -n drop_env -c conda-forge -c bioconda drop

Installation time: ~ 10min

Test whether the pipeline runs through by setting up the demo dataset in an empty directory (e.g. ``~/drop_demo``).

.. code-block:: bash

    mkdir ~/drop_demo
    cd ~/drop_demo

    # demo will download the necessary data and pipeline files
    drop demo

The pipeline can be run using `snakemake <snakemake.readthedocs.io/>`_ commands

.. code-block:: bash

    snakemake -n # dryrun
    snakemake --cores 1

Initialize a project
--------------------
The demo project can be modified to be used for a new project.
Alternatively, a new DROP project can be set up using ``drop init``.

.. code-block:: bash

    cd <path/to/project>
    drop init

This will create an empty ``config.yaml`` file that needs to be filled according to the project data.
You also need to prepare a sample annotation file.
Go to :doc:`prepare` for more details.


.. _otherversions:

Other DROP versions
-------------------

The developer version of DROP can be found in the `repository <https://github.com/gagneurlab/drop>`_ under the branch
``dev``.
Make sure that the :any:`prerequisites` are installed, preferably in a conda environment.
Then install DROP from github using ``pip``.

.. code-block:: bash

    pip install git+https://github.com/gagneurlab/drop.git@dev


Alternatively, you can clone the desired branch of the repository and install from directory.

.. code-block:: bash

    git clone -b dev https://github.com/gagneurlab/drop.git
    pip install ./drop

If the package needs to be updated frequently, it is more useful to use the ``-e` option of ``pip``.
Any new update pulled from the repository will be available without reinstall.
Note, that this requires an explicit call to update any existing project (:any:`dropUpdate`).

.. code-block::

    pip install -e ./drop

    # update project directory
    cd <path/to/project>
    drop update


.. _prerequisites:

Prerequisites
-------------

The easiest way to ensure that all dependencies are installed is to install the bioconda package, as described above.
Once the environment is set up and installation was successful, other versions of drop can be installed with ``pip``,
overwriting the conda version of ``DROP`` (see :any:`otherversions`).


Installation without conda
++++++++++++++++++++++++++
Alternatively, DROP can be installed without ``conda``. In this case the following dependencies must be met:

* Programming languages:

  * `python <https://www.python.org/>`_ >= 3.6 and `pip <https://pip.pypa.io/en/stable/installing/>`_ >= 19.1

  * `R <https://www.r-project.org/>`_ >= 3.6, <=4.0.2 and corresponding `bioconductor <https://bioconductor.org/install/>`_ version

* Commandline tools:

    * `GNU bc <https://www.gnu.org/software/bc/>`_

    * `GNU wget <https://www.gnu.org/software/wget/>`_

    * `tabix <https://www.htslib.org/download/>`_

    * `samtools <https://www.htslib.org/download/>`_ >= 1.7

    * `bcftools <https://github.com/samtools/bcftools>`_ >= 1.7

    * `GATK <https://software.broadinstitute.org/gatk/>`_ >= 4.1.8

    * `graphviz <https://www.graphviz.org/>`_

    * `pandoc <https://pandoc.org/>`_


.. note::

    If you are using an already existing R installation, make sure that the R and bioconductor versions match.
    Otherwise, use the newest versions of R and bioconductor.

At first invocation, all necessary R packages will be installed with the first pipeline call.
As this is a lengthy process, it might be desirable to install them in advance, if a local copy of the repository exists.

.. code-block:: bash

    # optional
    Rscript <path/to/drop/repo>/drop/installRPackages.R drop/requirementsR.txt


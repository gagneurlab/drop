# Use Mambaforge as base - it's optimized for conda-forge and bioconda
FROM condaforge/mambaforge:latest

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        build-essential \
        git \
        wget \
		unzip \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove -y \
    && apt-get clean

# Prevent conda from activating the base environment
ENV CONDA_AUTO_ACTIVATE_BASE=false

# Install DROP using mamba
# Combining commands to reduce layers and cleanup to reduce image size
RUN mamba create -n drop_env -c conda-forge -c bioconda drop \
    bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg38 \
    bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg19 \
    --override-channels -y && \
    mamba clean -afy && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.pyc' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete

# Set up the environment activation script
RUN echo '. /opt/conda/etc/profile.d/conda.sh && conda activate drop_env' >> /.bashrc && \
    echo '. /opt/conda/etc/profile.d/conda.sh && conda activate drop_env' >> /.profile

# Make sure the conda environment is activated for any command
ENV PATH /opt/conda/envs/drop_env/bin:$PATH

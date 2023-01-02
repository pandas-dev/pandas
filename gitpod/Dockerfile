#
# Dockerfile for pandas development
#
# Usage:
# -------
#
# To make a local build of the container, from the 'Docker-dev' directory:
# docker build  --rm -f "Dockerfile" -t <build-tag> "."
#
# To use the container use the following command. It assumes that you are in
# the root folder of the pandas git repository, making it available as
# /home/pandas in the container. Whatever changes you make to that directory
# are visible in the host and container.
# The docker image is retrieved from the pandas dockerhub repository
#
# docker run --rm -it -v $(pwd):/home/pandas pythonpandas/pandas-dev:<image-tag>
#
# By default the container will activate the conda environment pandas-dev
# which contains all the dependencies needed for pandas development
#
# To build and install pandas run:
#               python setup.py build_ext -j 4
#               python -m pip install -e . --no-build-isolation
#
# This image is based on: Ubuntu 20.04 (focal)
# https://hub.docker.com/_/ubuntu/?tab=tags&name=focal
# OS/ARCH: linux/amd64
FROM gitpod/workspace-base:latest

ARG MAMBAFORGE_VERSION="22.9.0-1"
ARG CONDA_ENV=pandas-dev
ARG PANDAS_HOME="/home/pandas"


# ---- Configure environment ----
ENV CONDA_DIR=/home/gitpod/mambaforge3 \
    SHELL=/bin/bash
ENV PATH=${CONDA_DIR}/bin:$PATH \
    WORKSPACE=/workspace/pandas

# -----------------------------------------------------------------------------
# ---- Creating as root - note: make sure to change to gitpod in the end ----
USER root

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

# Configure apt and install packages
RUN apt-get update \
    && apt-get -y install --no-install-recommends apt-utils dialog 2>&1 \
    #
    # Install tzdata and configure timezone (fix for tests which try to read from "/etc/localtime")
    && apt-get -y install tzdata \
    && ln -fs /usr/share/zoneinfo/Etc/UTC /etc/localtime \
    && dpkg-reconfigure -f noninteractive tzdata \
    #
    # Verify git, process tools, lsb-release (common in install instructions for CLIs) installed
    && apt-get -y install git iproute2 procps iproute2 lsb-release \
    #
    # cleanup
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog

# Allows this Dockerfile to activate conda environments
SHELL ["/bin/bash", "--login", "-o", "pipefail", "-c"]

# -----------------------------------------------------------------------------
# ---- Installing mamba  ----
RUN wget -q -O mambaforge3.sh \
    "https://github.com/conda-forge/miniforge/releases/download/$MAMBAFORGE_VERSION/Mambaforge-$MAMBAFORGE_VERSION-Linux-x86_64.sh" && \
    bash mambaforge3.sh -p ${CONDA_DIR} -b && \
    rm mambaforge3.sh

# -----------------------------------------------------------------------------
# ---- Copy needed files ----
# basic workspace configurations
COPY ./gitpod/workspace_config /usr/local/bin/workspace_config

RUN chmod a+rx /usr/local/bin/workspace_config && \
    workspace_config

# the container to create a conda environment from it
COPY environment.yml /tmp/environment.yml

# ---- Create conda environment ----
RUN mamba env create -f /tmp/environment.yml && \
    conda activate $CONDA_ENV && \
    mamba install ccache -y && \
    # needed for docs rendering later on
    python -m pip install --no-cache-dir sphinx-autobuild && \
    conda clean --all -f -y && \
    rm -rf /tmp/*

# -----------------------------------------------------------------------------
# Always make sure we are not root
USER gitpod

# Doing a local shallow clone - keeps the container secure
# and much slimmer than using COPY directly or making a
# remote clone
#ARG BASE_CONTAINER="pandas/pandas-dev:latest"
FROM gitpod/workspace-base:latest as clone

USER root
# Install base utilities
RUN apt-get update \
    && apt-get install -y build-essential \
    && apt-get install -y wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# init conda
RUN conda init

# Use the speedy libmamba solver
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba

# # -----------------------------------------------------------------------------
# # Always return to non privileged user
RUN chown -R gitpod:gitpod /opt/conda/pkgs
RUN chown -R gitpod:gitpod /opt/conda/envs
RUN chown -R gitpod:gitpod /home/gitpod/.conda
USER gitpod

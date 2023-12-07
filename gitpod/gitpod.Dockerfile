FROM gitpod/workspace-base:2023-11-24-15-04-57

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

# Init conda and use the speedy libmamba solver
RUN conda init && conda install -n base conda-libmamba-solver && conda config --set solver libmamba

# Install dependencies
RUN conda env update --file https://raw.githubusercontent.com/pandas-dev/pandas/main/environment.yml --prune

# # -----------------------------------------------------------------------------
# # Always return to non privileged user
RUN chown -R gitpod:gitpod /opt/conda/pkgs
RUN chown -R gitpod:gitpod /opt/conda/envs
RUN chown -R gitpod:gitpod /home/gitpod/.conda
RUN chown -R gitpod:gitpod /home/gitpod/.cache
USER gitpod

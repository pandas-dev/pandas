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

# Use the speedy libmamba solver
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba

# # -----------------------------------------------------------------------------
# # Using the pandas-dev Docker image as a base
# # This way, we ensure we have all the needed compilers and dependencies
# # while reducing the build time
# #FROM ${BASE_CONTAINER} as build

# # -----------------------------------------------------------------------------
# USER root

# # -----------------------------------------------------------------------------
# # ---- ENV variables ----
# # ---- Directories needed ----
# ENV WORKSPACE=/workspace/pandas/ \
#     CONDA_ENV=pandas-dev

# # Allows this micromamba.Dockerfile to activate conda environments
# SHELL ["/bin/bash", "--login", "-o", "pipefail", "-c"]

# # Copy over the shallow clone
# COPY --from=clone --chown=gitpod /tmp/pandas ${WORKSPACE}

# # Everything happens in the /workspace/pandas directory
# WORKDIR ${WORKSPACE}

# # Build pandas to populate the cache used by ccache
# RUN git config --global --add safe.directory /workspace/pandas
# RUN conda activate ${CONDA_ENV} && \
#     python -m pip install -ve . --no-build-isolation --config-settings=editable-verbose=true && \
#     ccache -s

# # Gitpod will load the repository into /workspace/pandas. We remove the
# # directory from the image to prevent conflicts
# RUN rm -rf ${WORKSPACE}

# # -----------------------------------------------------------------------------
# # Always return to non privileged user
# RUN chown -R gitpod:gitpod /home/gitpod/.cache/
USER gitpod

# Doing a local shallow clone - keeps the container secure
# and much slimmer than using COPY directly or making a
# remote clone
ARG BASE_CONTAINER="pythonpandas/pandas-dev:latest"
FROM gitpod/workspace-base:latest as clone

# the clone should be deep enough for versioneer to work
RUN git clone https://github.com/pandas-dev/pandas --depth 12 /tmp/pandas

# -----------------------------------------------------------------------------
# Using the pandas-dev Docker image as a base
# This way, we ensure we have all the needed compilers and dependencies
# while reducing the build time
FROM ${BASE_CONTAINER} as build

# -----------------------------------------------------------------------------
USER root

# -----------------------------------------------------------------------------
# ---- ENV variables ----
# ---- Directories needed ----
ENV WORKSPACE=/workspace/pandas/ \
    CONDA_ENV=pandas-dev

# Allows this micromamba.Dockerfile to activate conda environments
SHELL ["/bin/bash", "--login", "-o", "pipefail", "-c"]

# Copy over the shallow clone
COPY --from=clone --chown=gitpod /tmp/pandas ${WORKSPACE}

# Everything happens in the /workspace/pandas directory
WORKDIR ${WORKSPACE}

# Build pandas to populate the cache used by ccache
RUN git config --global --add safe.directory /workspace/pandas
RUN conda activate ${CONDA_ENV} && \
    python setup.py build_ext --inplace && \
    ccache -s

# Gitpod will load the repository into /workspace/pandas. We remove the
# directory from the image to prevent conflicts
RUN rm -rf ${WORKSPACE}

# -----------------------------------------------------------------------------
# Always return to non privileged user
USER gitpod

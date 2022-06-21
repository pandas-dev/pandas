FROM quay.io/condaforge/miniforge3

# if you forked pandas, you can pass in your own GitHub username to use your fork
# i.e. gh_username=myname
ARG gh_username=WillAyd
ARG pandas_home="/home/pandas"

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
    && apt-get -y install build-essential git iproute2 procps iproute2 lsb-release \
    #
    # cleanup
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog

# Without this versioneer will not be able to determine the pandas version.
# This is because of a security update to git that blocks it from reading the config folder if
# it is not owned by the current user. We hit this since the "mounted" folder is not hit by the
# Docker container.
# xref https://github.com/pypa/manylinux/issues/1309
RUN git config --global --add safe.directory "$pandas_home"

# Clone pandas repo
RUN mkdir "$pandas_home" \
    && git clone "https://github.com/$gh_username/pandas.git" "$pandas_home" \
    && cd "$pandas_home" \
    && git checkout cmake-build \
    && git remote add upstream "https://github.com/pandas-dev/pandas.git" \
    && git pull upstream main

# Because it is surprisingly difficult to activate a conda environment inside a DockerFile
# (from personal experience and per https://github.com/ContinuumIO/docker-images/issues/89),
# we just update the base/root one from the 'environment.yml' file instead of creating a new one.
#
# Set up environment
RUN conda install -y mamba
RUN mamba env update -n base -f "$pandas_home/environment.yml"

# Build C extensions and pandas
SHELL ["/bin/bash", "-c"]
RUN . /opt/conda/etc/profile.d/conda.sh \
    && conda activate base \
    && cd "$pandas_home" \
    && export \
    && cmake . \
    && cmake --build . --parallel \
    && python -m pip install --no-build-isolation -e .

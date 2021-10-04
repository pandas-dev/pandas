FROM quay.io/condaforge/miniforge3

# if you forked pandas, you can pass in your own GitHub username to use your fork
# i.e. gh_username=myname
ARG gh_username=pandas-dev
ARG pandas_home="/home/pandas"

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

# Configure apt and install packages
RUN apt-get update \
    && apt-get -y install --no-install-recommends apt-utils dialog 2>&1 \
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

# Clone pandas repo
RUN mkdir "$pandas_home" \
    && git clone "https://github.com/$gh_username/pandas.git" "$pandas_home" \
    && cd "$pandas_home" \
    && git remote add upstream "https://github.com/pandas-dev/pandas.git" \
    && git pull upstream master

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
    && python setup.py build_ext -j 4 \
    && python -m pip install -e .

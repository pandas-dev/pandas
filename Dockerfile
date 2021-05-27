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



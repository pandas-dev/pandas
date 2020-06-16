FROM gitpod/workspace-full:latest

USER root

RUN apt-get update \
    && apt-get -y install --no-install-recommends apt-utils dialog 2>&1 \
    #
    # Verify git, process tools, lsb-release (common in install instructions for CLIs) installed
    && apt-get -y install git iproute2 procps iproute2 lsb-release \
    #
    # Install C compilers (gcc not enough, so just went with build-essential which admittedly might be overkill),
    # needed to build pandas C extensions
    && apt-get -y install build-essential \
    #
    # cleanup
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

USER gitpod

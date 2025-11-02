FROM python:3.11.13
WORKDIR /home/pandas

# https://docs.docker.com/reference/dockerfile/#automatic-platform-args-in-the-global-scope
ARG TARGETPLATFORM

RUN apt-get update && \
    apt-get --no-install-recommends -y upgrade && \
    apt-get --no-install-recommends -y install \
    build-essential \
    bash-completion \
    # hdf5 needed for pytables installation
    libhdf5-dev \
    # libgles2-mesa needed for pytest-qt
    libgles2-mesa-dev && \
    rm -rf /var/lib/apt/lists/*

COPY requirements-dev.txt /tmp

RUN case "$TARGETPLATFORM" in \
    linux/arm*) \
        # Drop PyQt5 for ARM GH#61037
        sed -i "/^pyqt5/Id" /tmp/requirements-dev.txt \
        ;; \
    esac && \
    python -m pip install --no-cache-dir --upgrade pip && \
    python -m pip install --no-cache-dir -r /tmp/requirements-dev.txt
RUN git config --global --add safe.directory /home/pandas

ENV SHELL="/bin/bash"
CMD ["/bin/bash"]

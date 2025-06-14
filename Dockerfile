FROM python:3.10.8
WORKDIR /home/pandas

RUN apt-get update && \
    apt-get --no-install-recommends -y upgrade && \
    apt-get --no-install-recommends -y install \
    build-essential \
    bash-completion \
    # Install Qt5 dependencies for pytest-qt, only for m chip Macs
    #-y qt5-qmake qtbase5-dev\
    # hdf5 needed for pytables installation
    libhdf5-dev \
    # libgles2-mesa needed for pytest-qt
    libgles2-mesa-dev && \
    rm -rf /var/lib/apt/lists/*

COPY requirements-dev.txt /tmp
RUN python -m pip install --no-cache-dir --upgrade pip && \
    python -m pip install --no-cache-dir -r /tmp/requirements-dev.txt
RUN git config --global --add safe.directory /home/pandas

ENV SHELL="/bin/bash"
CMD ["/bin/bash"]

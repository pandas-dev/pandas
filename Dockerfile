FROM python:3.10.8
WORKDIR /home/pandas

RUN apt-get update && \
    apt-get -y upgrade && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y \
    build-essential \
    bash-completion \
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

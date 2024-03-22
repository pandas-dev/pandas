FROM python:3.10.8
WORKDIR /home/pandas

RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y build-essential

# hdf5 needed for pytables installation
# libgles2-mesa needed for pytest-qt
RUN apt-get install -y libhdf5-dev libgles2-mesa-dev

RUN python -m pip install --upgrade pip
COPY requirements-dev.txt /tmp
RUN python -m pip install -r /tmp/requirements-dev.txt
RUN git config --global --add safe.directory /home/pandas
CMD ["/bin/bash"]

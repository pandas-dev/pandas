FROM python:3.13.0a3
WORKDIR /home/pandas

RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y build-essential

# hdf5 needed for pytables installation
# libgles2-mesa needed for pytest-qt
RUN apt-get install -y libhdf5-dev libgles2-mesa-dev

RUN python -m pip install --upgrade pip
RUN python -m pip install \
    -r https://raw.githubusercontent.com/pandas-dev/pandas/main/requirements-dev.txt
CMD ["/bin/bash"]

FROM python:3.10.8
WORKDIR /home/pandas

RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y build-essential

# hdf5 needed for pytables installation
RUN apt-get install -y libhdf5-dev

RUN python -m pip install --upgrade pip
RUN python -m pip install \
    -r https://raw.githubusercontent.com/pandas-dev/pandas/main/requirements-dev.txt

ARG USERNAME=pandas-dev
ARG USER_UID=1000
ARG USER_GID=$USER_UID
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME -s /bin/bash \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

CMD ["/bin/bash"]

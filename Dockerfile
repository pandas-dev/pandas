FROM python:3.10.8
WORKDIR /home/pandas

RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y build-essential

# hdf5 needed for pytables installation
RUN apt-get install -y libhdf5-dev

RUN apt-get install -y vim less

ARG USERNAME=pandas-dev
ARG USER_UID=1000
ARG USER_GID=$USER_UID
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME -s /bin/bash \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

RUN apt-get install -y curl bzip2 && sudo -u $USERNAME bash -c \
    "mkdir ~/bin && \
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | \
    tar -C ~ -xj bin/micromamba \
    && mv ~/bin ~/micromamba-bin"

CMD ["/bin/bash"]

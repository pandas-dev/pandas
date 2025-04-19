#!/bin/bash
# Basic configurations for the workspace

set -e

# gitpod/workspace-base needs at least one file here
touch /home/gitpod/.bashrc.d/empty

# Add git aliases
git config --global alias.co checkout
git config --global alias.ci commit
git config --global alias.st status
git config --global alias.br branch
git config --global alias.hist "log --pretty=format:'%h %ad | %s%d [%an]' --graph --date=short"
git config --global alias.type 'cat-file -t'
git config --global alias.dump 'cat-file -p'

# Enable basic vim defaults in ~/.vimrc
echo "filetype plugin indent on" >>~/.vimrc
echo "set colorcolumn=80" >>~/.vimrc
echo "set number" >>~/.vimrc
echo "syntax enable" >>~/.vimrc

# Vanity custom bash prompt - makes it more legible
echo "PS1='\[\e]0;\u \w\a\]\[\033[01;36m\]\u\[\033[m\] > \[\033[38;5;141m\]\w\[\033[m\] \\$ '" >>~/.bashrc

# Enable prompt color in the skeleton .bashrc
# hadolint ignore=SC2016
sed -i 's/^#force_color_prompt=yes/force_color_prompt=yes/' /etc/skel/.bashrc

# .gitpod.yml is configured to install pandas from /workspace/pandas
echo "export PYTHONPATH=${WORKSPACE}" >>~/.bashrc

# make conda activate command available from /bin/bash (login and interactive)
if [[ ! -f "/etc/profile.d/conda.sh" ]]; then
    ln -s ${CONDA_DIR}/etc/profile.d/conda.sh /etc/profile.d/conda.sh
fi
echo ". ${CONDA_DIR}/etc/profile.d/conda.sh" >>~/.bashrc
echo "conda activate pandas-dev" >>~/.bashrc

# Enable prompt color in the skeleton .bashrc
# hadolint ignore=SC2016
sed -i 's/^#force_color_prompt=yes/force_color_prompt=yes/' /etc/skel/.bashrc

# .gitpod.yml is configured to install pandas from /workspace/pandas
echo "export PYTHONPATH=/workspace/pandas" >>~/.bashrc

# Set up ccache for compilers for this Dockerfile
# REF: https://github.com/conda-forge/compilers-feedstock/issues/31
echo "conda activate pandas-dev" >>~/.startuprc
echo "export CC=\"ccache \$CC\"" >>~/.startuprc
echo "export CXX=\"ccache \$CXX\"" >>~/.startuprc
echo "source ~/.startuprc" >>~/.profile
echo "source ~/.startuprc" >>~/.bashrc

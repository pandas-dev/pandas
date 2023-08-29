#! /usr/bin/env bash

set -e

~/micromamba-bin/micromamba \
    create \
    -y \
    -r ~/micromamba \
    -f environment.yml \
    -n test \
    --rc-file ci/condarc.yml

mv ~/.bashrc ~/.bashrc.orig
~/micromamba-bin/micromamba shell init --shell bash --root-prefix=~/micromamba
mv ~/.bashrc ~/.mamba-bashrc
cp ~/.bashrc.orig ~/.bashrc
echo micromamba activate test >> ~/.mamba-bashrc
cat ~/.mamba-bashrc >> ~/.bashrc
cat ~/.mamba-bashrc >> ~/.bash_profile
source ~/.mamba-bashrc


if [[ $(git tag | wc -c) -eq 0 ]]
then
    echo -e '\n\nERROR:'
    echo Your repository has been cloned without tags.
    echo To fix this, add the original pandas GitHub repository as an upstream
    echo reference, e.g., via the following:
    echo git remote add upstream https://github.com/pandas-dev/pandas.git
    echo Then, run
    echo git fetch --tags upstream
    echo to get the main repository\'s tags, and optionally run
    echo git push --tags
    echo to upload the tags to your GitHub fork.
    echo
    echo Finally, rebuild the devcontainer.
    exit 1
fi

pip install black

pip uninstall -y pandas
pip install -e . --no-build-isolation -v

pre-commit install

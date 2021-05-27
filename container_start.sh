# Switch back to dialog for any ad-hoc use of apt-get
export DEBIAN_FRONTEND=dialog

# Set up environment
conda install -y mamba
mamba env update -p .venv -f "environment.yml"

source activate /workspaces/pandas/.venv

# build + install pandas
python setup.py build_ext -j 4
python -m pip install -e .
FROM condaforge/miniforge3:23.3.1-1

# Init conda and use the speedy libmamba solver
RUN conda init && conda install -n base conda-libmamba-solver && conda config --set solver libmamba

# Install dependencies
RUN conda env update --file https://raw.githubusercontent.com/pandas-dev/pandas/main/environment.yml --prune


#!/bin/bash

source activate pandas

echo "[install 3.6 downstream deps]"

conda install -n pandas -c conda-forge pandas-datareader xarray geopandas seaborn statsmodels scikit-learn dask

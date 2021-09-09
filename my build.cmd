call conda activate pandas-dev
python setup.py build_ext -j 4
python -m pip install -e . --no-build-isolation --no-use-pep517
pause

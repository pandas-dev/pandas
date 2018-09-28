@rem https://github.com/numba/numba/blob/master/buildscripts/incremental/build.cmd
call activate %CONDA_ENV%

@rem Build numba extensions without silencing compile errors
python setup.py build_ext -q --inplace

@rem Install pandas locally
python -m pip install -e .

if %errorlevel% neq 0 exit /b %errorlevel%

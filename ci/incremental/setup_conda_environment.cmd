@rem https://github.com/numba/numba/blob/master/buildscripts/incremental/setup_conda_environment.cmd
@rem The cmd /C hack circumvents a regression where conda installs a conda.bat
@rem script in non-root environments.
set CONDA_INSTALL=cmd /C conda install -q -y
set PIP_INSTALL=pip install -q

@echo on

@rem Deactivate any environment
call deactivate
@rem Display root environment (for debugging)
conda list
@rem Clean up any left-over from a previous build
conda remove --all -q -y -n %CONDA_ENV%
@rem Scipy, CFFI, jinja2 and IPython are optional dependencies, but exercised in the test suite
conda env create -n %CONDA_ENV% --file=ci\azure-windows-%CONDA_PY%.yml

call activate %CONDA_ENV%
conda list

if %errorlevel% neq 0 exit /b %errorlevel%

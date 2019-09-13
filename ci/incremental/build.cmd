@rem https://github.com/numba/numba/blob/master/buildscripts/incremental/build.cmd

@rem Build extensions
python setup.py build_ext -q -i

@rem Install pandas
python -m pip install --no-build-isolation -e .

if %errorlevel% neq 0 exit /b %errorlevel%

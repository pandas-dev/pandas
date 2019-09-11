@rem https://github.com/numba/numba/blob/master/buildscripts/incremental/build.cmd

@rem Build wheel
python setup.py bdist_wheel -d wheelhouse

@rem Install wheel
python -m pip install --no-deps wheelhouse/*.whl

if %errorlevel% neq 0 exit /b %errorlevel%

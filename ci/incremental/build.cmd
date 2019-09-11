@rem https://github.com/numba/numba/blob/master/buildscripts/incremental/build.cmd

@rem Build wheel
python setup.py bdist_wheel -d wheelhouse

@rem Install wheel
FOR /F "Tokens=*" %A IN ('DIR "wheelhouse\*.whl"') DO @( python -m pip install --no-deps %~nxA )

if %errorlevel% neq 0 exit /b %errorlevel%

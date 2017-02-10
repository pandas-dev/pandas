# nosetests -A "not slow and not network" pandas --with-id $*
pytest pandas --skip-slow

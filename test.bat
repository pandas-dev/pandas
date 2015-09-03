:: test on windows

nosetests --exe -A "not slow and not network and not disabled" pandas %*

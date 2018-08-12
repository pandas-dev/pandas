:: test on windows

pytest -W error:ResourceWarning --skip-slow --skip-network pandas -n 2 -r sxX --strict %*

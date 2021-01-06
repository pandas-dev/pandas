:: test on windows
set PYTHONHASHSEED=314159265
pytest --skip-slow --skip-network --skip-db -m "not single" -n 8 -r sXX --strict pandas

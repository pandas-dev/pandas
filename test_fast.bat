:: test on windows
set PYTHONHASHSEED=314159265
pytest --skip-slow --skip-network --skip-db -m "not single" -n 4 -r sXX pandas

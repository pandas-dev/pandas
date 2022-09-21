set test_command=import pandas as pd; print(pd.__version__); ^
pd.test(extra_args=['-m not clipboard and not single_cpu', '--skip-slow', '--skip-network', '--skip-db', '-n=2']); ^
pd.test(extra_args=['-m not clipboard and single_cpu', '--skip-slow', '--skip-network', '--skip-db'])

python --version
pip install pytz six numpy python-dateutil
pip install hypothesis==6.52.1 pytest>=6.2.5 pytest-xdist pytest-asyncio>=0.17
pip install --find-links=pandas/dist --no-index pandas
python -c "%test_command%"

set test_command=import pandas as pd; print(pd.__version__); ^
pd.test(extra_args=['-m not clipboard and not single_cpu and not slow and not network and not db', '-n 2']); ^
pd.test(extra_args=['-m not clipboard and single_cpu and not slow and not network and not db'])

python --version
pip install pytz six numpy python-dateutil
pip install hypothesis>=6.34.2 pytest>=7.0.0 pytest-xdist>=2.2.0 pytest-asyncio>=0.17
pip install --find-links=pandas/dist --no-index pandas
python -c "%test_command%"

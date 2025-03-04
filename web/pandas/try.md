# Try pandas in your browser (experimental)

Try our experimental [JupyterLite](https://jupyterlite.readthedocs.io/en/stable/) live shell with `pandas`, powered by [Pyodide](https://pyodide.org/en/stable/).

**Please note it can take a while (>30 seconds) before the shell is initialized and ready to run commands.**

**Running it requires a reasonable amount of bandwidth and resources (>70 MiB on the first load), so it may not work properly on all devices or networks.**

<iframe
  src="./lite/repl/index.html?toolbar=1&kernel=python&execute=0&code=import%20pandas%20as%20pd%0Adf%20%3D%20pd.DataFrame%28%7B%22num_legs%22%3A%20%5B2%2C%204%5D%2C%20%22num_wings%22%3A%20%5B2%2C%200%5D%7D%2C%20index%3D%5B%22falcon%22%2C%20%22dog%22%5D%29%0Adf"
  style="width: 100%; max-width: 650px; height: 600px; border: 1px solid #130753;"
></iframe>

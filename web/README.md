Directory containing the pandas website (hosted at https://pandas.pydata.org).

The website sources are in `web/pandas/`, which also include a `config.yml` file
containing the settings to build the website. The website is generated with the
command `./pandas_web.py pandas`. See `./pandas_web.py --help` and the header of
the script for more information and options.

**🚀 Quick Start for Contributors**

1. Ensure Python 3.8+ and pip are installed
2. Install dependencies: `pip install -r requirements.txt`
3. **Build the site:** `./pandas_web.py pandas`
4. **Preview locally:** `cd web/build && python -m http.server 8000`
5. Open http://localhost:8000 in your browser

After building the website, a local http server is needed to access the built website (it is not possible to open the local files with the browser, since the links and
the image sources are absolute to where they are served from). The easiest way
to run an http server locally is to run `python -m http.server` from the
`web/build/` directory.

See the [Good First Issue label](https://github.com/pandas-dev/pandas/labels/Good%20first%20issue) 
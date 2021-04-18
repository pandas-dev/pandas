Directory containing the pandas website (hosted at https://pandas.pydata.org).

The website sources are in `web/pandas/`, which also include a `config.yml` file
containing the settings to build the website. The website is generated with the
command `./pandas_web.py pandas`. See `./pandas_web.py --help` and the header of
the script for more information and options.

After building the website, to navigate it, it is needed to access the web using
an http server (a not open the local files with the browser, since the links and
the image sources are absolute to where they are served from). The easiest way
to run an http server locally is to run `python -m http.server` from the
`web/build/` directory.

# The interactive `pandas` REPL

An interactive REPL to easily try `pandas` in the browser, powered by JupyterLite.

![image](https://user-images.githubusercontent.com/591645/175000291-e8c69f6f-5f2c-48d7-817c-cff05ab2cde9.png)

## Build

The interactive REPL can be as a part of the documentation build process
with Sphinx using the `jupyterlite-sphinx` extension. Please refer to the
`doc/make.py` file and the [pandas docs development workflow](https://pandas.pydata.org/docs/development/contributing_documentation.html#how-to-build-the-pandas-documentation) for more information.

## Configuration

The `doc/source/` folder contains shared configuration files for the interactive terminal powered by JupyterLite:

- `jupyter_lite_config.json`: build time configuration, used when building the assets with the `jupyter lite build` command
- `jupyter-lite.json` run time configuration applied when launching the application in the browser

This interactive `pandas` JupyterLite deployment enables optimizations by removing unused shared packages
and disabling source maps, which can make the assets smaller and faster to load, at the cost of debugging
capabilities.

To learn more about it, check out the JupyterLite documentation:

- Optimizations: https://jupyterlite.readthedocs.io/en/latest/howto/configure/advanced/optimizations.html
- JupyterLite schema: https://jupyterlite.readthedocs.io/en/latest/reference/schema-v0.html
- `jupyterlite-sphinx` extension: https://jupyterlite-sphinx.readthedocs.io/en/stable/
- `jupyter lite` CLI reference: https://jupyterlite.readthedocs.io/en/latest/reference/cli.html

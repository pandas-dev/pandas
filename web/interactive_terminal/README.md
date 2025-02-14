# The interactive `pandas` REPL

An interactive REPL to easily try `pandas` in the browser, powered by JupyterLite.

![image](https://user-images.githubusercontent.com/591645/175000291-e8c69f6f-5f2c-48d7-817c-cff05ab2cde9.png)

## Build

The interactive REPL is built with the `jupyter lite` CLI.

First make sure `jupyterlite` and a kernel are installed:

```bash
python -m pip install jupyterlite-core
python -m pip install jupyterlite-pyodide-kernel
```

Then in `web/interactive_terminal`, run the following command:

```bash
jupyter lite build
```

## Configuration

This folder contains configuration files for the interactive terminal powered by JupyterLite:

- `jupyter_lite_config.json`: build time configuration, used when building the assets with the `jupyter lite build` command
- `jupyter-lite.json` run time configuration applied when launching the application in the browser

This interactive `pandas` JupyterLite deployment enables a couple of optimizations to only include the `repl` app in the generated static assets, and disables source maps, which can make the assets smaller and faster to load, at the cost of
debugging capabilities.

To learn more about it, check out the JupyterLite documentation:

- Optimizations: https://jupyterlite.readthedocs.io/en/latest/howto/configure/advanced/optimizations.html
- JupyterLite schema: https://jupyterlite.readthedocs.io/en/latest/reference/schema-v0.html
- CLI reference: https://jupyterlite.readthedocs.io/en/latest/reference/cli.html

# The interactive `pandas` terminal

An interactive terminal to easily try `pandas` in the browser, powered by JupyterLite.

![image](https://user-images.githubusercontent.com/591645/175000291-e8c69f6f-5f2c-48d7-817c-cff05ab2cde9.png)

## Build

The interactive terminal is built with the `jupyterlite` CLI.

First make sure `jupyterlite` is installed:

```bash
python -m pip install jupyterlite
```

Then in `web/interactive_terminal`, run the following command:

```bash
jupyter lite build
```

## Configuration

This folder contains configuration files for the interactive terminal powered by JupyterLite:

- `jupyter_lite_config.json`: build time configuration, used when building the assets with the `jupyter lite build` command
- `jupyter-lite.json` run time configuration applied when launching the application in the browser

The interactive `pandas` terminal application enables a couple of optimizations to only include the `repl` app in the generated static assets.
To learn more about it, check out the JupyterLite documentation:

- Optimizations: https://jupyterlite.readthedocs.io/en/latest/howto/configure/advanced/optimizations.html
- JupyterLite schema: https://jupyterlite.readthedocs.io/en/latest/reference/schema-v0.html
- CLI reference: https://jupyterlite.readthedocs.io/en/latest/reference/cli.html

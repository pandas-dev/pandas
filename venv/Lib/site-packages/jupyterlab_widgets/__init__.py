# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from ._version import __version__


def _jupyter_labextension_paths():
    import sys
    from pathlib import Path

    labext_name = '@jupyter-widgets/jupyterlab-manager'
    here = Path(__file__).parent.resolve()
    src_prefix = here.parent / 'labextension'

    if not src_prefix.exists():
        src_prefix = Path(sys.prefix) / f'share/jupyter/labextensions/{labext_name}'

    return [{'src': str(src_prefix), 'dest': labext_name}]

__all__ = ['_jupyter_labextension_paths', '__version__']

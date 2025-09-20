"""Interactive widgets for the Jupyter notebook.

Provide simple interactive controls in the notebook.
Each widget corresponds to an object in Python and Javascript,
with controls on the page.

You can display widgets with IPython's display machinery::

    from ipywidgets import IntSlider
    from IPython.display import display
    slider = IntSlider(min=1, max=10)
    display(slider)

Moving the slider will change the value. Most widgets have a current value,
accessible as a `value` attribute.
"""
from ._version import __version__
from warnings import warn

def _jupyter_nbextension_paths():
    return [{
        'section': 'notebook',
        'src': 'static',
        'dest': 'jupyter-js-widgets',
        'require': 'jupyter-js-widgets/extension'
    }]

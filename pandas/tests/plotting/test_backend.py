import pytest

import pandas


def _mocked_import_module(name):
    """
    Mock of ``importlib.import_module``. Depending of the name of the module
    received will return:

    - 'sample_backend': A mock of a valid plotting backend
    - 'module_not_a_backend': A backend (object) with none of the backend
        methods
    - 'backend_missing_area_plot': A backend (object) with all backend
        attributes except ``AreaPlot``
    """
    class PlottingBackendModuleMock:
        LinePlot = None
        BarPlot = None
        BarhPlot = None
        HistPlot = None
        BoxPlot = None
        KdePlot = None
        AreaPlot = None
        PiePlot = None
        ScatterPlot = None
        HexBinPlot = None
        hist_series = None
        hist_frame = None
        boxplot = None
        boxplot_frame = None
        boxplot_frame_groupby = None

    if name == 'correct_backend':
        return PlottingBackendModuleMock
    elif name == 'module_not_a_backend':
        return object()
    elif name == 'backend_missing_area_plot':
        mod = PlottingBackendModuleMock
        del mod.AreaPlot
        return mod

    raise ValueError('Unknown mocked backend: {}'.format(name))


def test_matplotlib_backend_error():
    msg = ('matplotlib is required for plotting when the default backend '
           '"matplotlib" is selected.')
    try:
        import matplotlib  # noqa
    except ImportError:
        with pytest.raises(ImportError, match=msg):
            pandas.set_option('plotting.backend', 'matplotlib')


def test_backend_is_not_module():
    msg = ('"not_an_existing_module" does not seem to be an installed module. '
           'A pandas plotting backend must be a module that can be imported')
    with pytest.raises(ValueError, match=msg):
        pandas.set_option('plotting.backend', 'not_an_existing_module')


def test_backend_not_a_backend_module(monkeypatch):
    required_objs = ['LinePlot', 'BarPlot', 'BarhPlot', 'HistPlot',
                     'BoxPlot', 'KdePlot', 'AreaPlot', 'PiePlot',
                     'ScatterPlot', 'HexBinPlot', 'hist_series',
                     'hist_frame', 'boxplot', 'boxplot_frame',
                     'boxplot_frame_groupby']
    msg = ('"module_not_a_backend" does not seem to be a valid backend. '
           'Valid backends are modules that implement the next '
           'objects:\n{}'.format('\n'.join(required_objs)))
    monkeypatch.setattr('pandas.core.config_init.importlib.import_module',
                        _mocked_import_module)
    with pytest.raises(ValueError, match=msg):
        pandas.set_option('plotting.backend', 'module_not_a_backend')


def test_backend_has_missing_objects(monkeypatch):
    msg = ('"backend_missing_area_plot" does not seem to be a complete '
           'backend. Valid backends must implement the next objects:\n'
           'AreaPlot')
    monkeypatch.setattr('pandas.core.config_init.importlib.import_module',
                        _mocked_import_module)
    with pytest.raises(ValueError, match=msg):
        pandas.set_option('plotting.backend', 'backend_missing_area_plot')


def test_backend_is_correct(monkeypatch):
    monkeypatch.setattr('pandas.core.config_init.importlib.import_module',
                        _mocked_import_module)
    pandas.set_option('plotting.backend', 'correct_backend')
    assert pandas.get_option('plotting.backend') == 'correct_backend'

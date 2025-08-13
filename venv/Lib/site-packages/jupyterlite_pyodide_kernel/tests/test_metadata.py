import jupyterlite_pyodide_kernel


def test_labextension_meta():
    paths = jupyterlite_pyodide_kernel._jupyter_labextension_paths()
    assert len(paths) == 1

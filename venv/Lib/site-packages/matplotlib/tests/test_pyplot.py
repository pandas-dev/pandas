import difflib
import inspect

import numpy as np
import sys
from pathlib import Path

import pytest

import matplotlib as mpl
from matplotlib.testing import subprocess_run_for_testing
from matplotlib import pyplot as plt


def test_pyplot_up_to_date(tmp_path):
    pytest.importorskip("black", minversion="24.1")

    gen_script = Path(mpl.__file__).parents[2] / "tools/boilerplate.py"
    if not gen_script.exists():
        pytest.skip("boilerplate.py not found")
    orig_contents = Path(plt.__file__).read_text()
    plt_file = tmp_path / 'pyplot.py'
    plt_file.write_text(orig_contents, 'utf-8')

    subprocess_run_for_testing(
        [sys.executable, str(gen_script), str(plt_file)],
        check=True)
    new_contents = plt_file.read_text('utf-8')

    if orig_contents != new_contents:
        diff_msg = '\n'.join(
            difflib.unified_diff(
                orig_contents.split('\n'), new_contents.split('\n'),
                fromfile='found pyplot.py',
                tofile='expected pyplot.py',
                n=0, lineterm=''))
        pytest.fail(
            "pyplot.py is not up-to-date. Please run "
            "'python tools/boilerplate.py' to update pyplot.py. "
            "This needs to be done from an environment where your "
            "current working copy is installed (e.g. 'pip install -e'd). "
            "Here is a diff of unexpected differences:\n%s" % diff_msg
        )


def test_copy_docstring_and_deprecators(recwarn):
    @mpl._api.rename_parameter(mpl.__version__, "old", "new")
    @mpl._api.make_keyword_only(mpl.__version__, "kwo")
    def func(new, kwo=None):
        pass

    @plt._copy_docstring_and_deprecators(func)
    def wrapper_func(new, kwo=None):
        pass

    wrapper_func(None)
    wrapper_func(new=None)
    wrapper_func(None, kwo=None)
    wrapper_func(new=None, kwo=None)
    assert not recwarn
    with pytest.warns(mpl.MatplotlibDeprecationWarning):
        wrapper_func(old=None)
    with pytest.warns(mpl.MatplotlibDeprecationWarning):
        wrapper_func(None, None)


def test_pyplot_box():
    fig, ax = plt.subplots()
    plt.box(False)
    assert not ax.get_frame_on()
    plt.box(True)
    assert ax.get_frame_on()
    plt.box()
    assert not ax.get_frame_on()
    plt.box()
    assert ax.get_frame_on()


def test_stackplot_smoke():
    # Small smoke test for stackplot (see #12405)
    plt.stackplot([1, 2, 3], [1, 2, 3])


def test_nrows_error():
    with pytest.raises(TypeError):
        plt.subplot(nrows=1)
    with pytest.raises(TypeError):
        plt.subplot(ncols=1)


def test_ioff():
    plt.ion()
    assert mpl.is_interactive()
    with plt.ioff():
        assert not mpl.is_interactive()
    assert mpl.is_interactive()

    plt.ioff()
    assert not mpl.is_interactive()
    with plt.ioff():
        assert not mpl.is_interactive()
    assert not mpl.is_interactive()


def test_ion():
    plt.ioff()
    assert not mpl.is_interactive()
    with plt.ion():
        assert mpl.is_interactive()
    assert not mpl.is_interactive()

    plt.ion()
    assert mpl.is_interactive()
    with plt.ion():
        assert mpl.is_interactive()
    assert mpl.is_interactive()


def test_nested_ion_ioff():
    # initial state is interactive
    plt.ion()

    # mixed ioff/ion
    with plt.ioff():
        assert not mpl.is_interactive()
        with plt.ion():
            assert mpl.is_interactive()
        assert not mpl.is_interactive()
    assert mpl.is_interactive()

    # redundant contexts
    with plt.ioff():
        with plt.ioff():
            assert not mpl.is_interactive()
    assert mpl.is_interactive()

    with plt.ion():
        plt.ioff()
    assert mpl.is_interactive()

    # initial state is not interactive
    plt.ioff()

    # mixed ioff/ion
    with plt.ion():
        assert mpl.is_interactive()
        with plt.ioff():
            assert not mpl.is_interactive()
        assert mpl.is_interactive()
    assert not mpl.is_interactive()

    # redundant contexts
    with plt.ion():
        with plt.ion():
            assert mpl.is_interactive()
    assert not mpl.is_interactive()

    with plt.ioff():
        plt.ion()
    assert not mpl.is_interactive()


def test_close():
    try:
        plt.close(1.1)
    except TypeError as e:
        assert str(e) == (
            "'fig' must be an instance of matplotlib.figure.Figure, int, str "
            "or None, not a float")


def test_subplot_reuse():
    ax1 = plt.subplot(121)
    assert ax1 is plt.gca()
    ax2 = plt.subplot(122)
    assert ax2 is plt.gca()
    ax3 = plt.subplot(121)
    assert ax1 is plt.gca()
    assert ax1 is ax3


def test_axes_kwargs():
    # plt.axes() always creates new axes, even if axes kwargs differ.
    plt.figure()
    ax = plt.axes()
    ax1 = plt.axes()
    assert ax is not None
    assert ax1 is not ax
    plt.close()

    plt.figure()
    ax = plt.axes(projection='polar')
    ax1 = plt.axes(projection='polar')
    assert ax is not None
    assert ax1 is not ax
    plt.close()

    plt.figure()
    ax = plt.axes(projection='polar')
    ax1 = plt.axes()
    assert ax is not None
    assert ax1.name == 'rectilinear'
    assert ax1 is not ax
    plt.close()


def test_subplot_replace_projection():
    # plt.subplot() searches for axes with the same subplot spec, and if one
    # exists, and the kwargs match returns it, create a new one if they do not
    fig = plt.figure()
    ax = plt.subplot(1, 2, 1)
    ax1 = plt.subplot(1, 2, 1)
    ax2 = plt.subplot(1, 2, 2)
    ax3 = plt.subplot(1, 2, 1, projection='polar')
    ax4 = plt.subplot(1, 2, 1, projection='polar')
    assert ax is not None
    assert ax1 is ax
    assert ax2 is not ax
    assert ax3 is not ax
    assert ax3 is ax4

    assert ax in fig.axes
    assert ax2 in fig.axes
    assert ax3 in fig.axes

    assert ax.name == 'rectilinear'
    assert ax2.name == 'rectilinear'
    assert ax3.name == 'polar'


def test_subplot_kwarg_collision():
    ax1 = plt.subplot(projection='polar', theta_offset=0)
    ax2 = plt.subplot(projection='polar', theta_offset=0)
    assert ax1 is ax2
    ax1.remove()
    ax3 = plt.subplot(projection='polar', theta_offset=1)
    assert ax1 is not ax3
    assert ax1 not in plt.gcf().axes


def test_gca():
    # plt.gca() returns an existing axes, unless there were no axes.
    plt.figure()
    ax = plt.gca()
    ax1 = plt.gca()
    assert ax is not None
    assert ax1 is ax
    plt.close()


def test_subplot_projection_reuse():
    # create an Axes
    ax1 = plt.subplot(111)
    # check that it is current
    assert ax1 is plt.gca()
    # make sure we get it back if we ask again
    assert ax1 is plt.subplot(111)
    # remove it
    ax1.remove()
    # create a polar plot
    ax2 = plt.subplot(111, projection='polar')
    assert ax2 is plt.gca()
    # this should have deleted the first axes
    assert ax1 not in plt.gcf().axes
    # assert we get it back if no extra parameters passed
    assert ax2 is plt.subplot(111)
    ax2.remove()
    # now check explicitly setting the projection to rectilinear
    # makes a new axes
    ax3 = plt.subplot(111, projection='rectilinear')
    assert ax3 is plt.gca()
    assert ax3 is not ax2
    assert ax2 not in plt.gcf().axes


def test_subplot_polar_normalization():
    ax1 = plt.subplot(111, projection='polar')
    ax2 = plt.subplot(111, polar=True)
    ax3 = plt.subplot(111, polar=True, projection='polar')
    assert ax1 is ax2
    assert ax1 is ax3

    with pytest.raises(ValueError,
                       match="polar=True, yet projection='3d'"):
        ax2 = plt.subplot(111, polar=True, projection='3d')


def test_subplot_change_projection():
    created_axes = set()
    ax = plt.subplot()
    created_axes.add(ax)
    projections = ('aitoff', 'hammer', 'lambert', 'mollweide',
                   'polar', 'rectilinear', '3d')
    for proj in projections:
        ax.remove()
        ax = plt.subplot(projection=proj)
        assert ax is plt.subplot()
        assert ax.name == proj
        created_axes.add(ax)
    # Check that each call created a new Axes.
    assert len(created_axes) == 1 + len(projections)


def test_polar_second_call():
    # the first call creates the axes with polar projection
    ln1, = plt.polar(0., 1., 'ro')
    assert isinstance(ln1, mpl.lines.Line2D)
    # the second call should reuse the existing axes
    ln2, = plt.polar(1.57, .5, 'bo')
    assert isinstance(ln2, mpl.lines.Line2D)
    assert ln1.axes is ln2.axes


def test_fallback_position():
    # check that position kwarg works if rect not supplied
    axref = plt.axes([0.2, 0.2, 0.5, 0.5])
    axtest = plt.axes(position=[0.2, 0.2, 0.5, 0.5])
    np.testing.assert_allclose(axtest.bbox.get_points(),
                               axref.bbox.get_points())

    # check that position kwarg ignored if rect is supplied
    axref = plt.axes([0.2, 0.2, 0.5, 0.5])
    axtest = plt.axes([0.2, 0.2, 0.5, 0.5], position=[0.1, 0.1, 0.8, 0.8])
    np.testing.assert_allclose(axtest.bbox.get_points(),
                               axref.bbox.get_points())


def test_set_current_figure_via_subfigure():
    fig1 = plt.figure()
    subfigs = fig1.subfigures(2)

    plt.figure()
    assert plt.gcf() != fig1

    current = plt.figure(subfigs[1])
    assert plt.gcf() == fig1
    assert current == fig1


def test_set_current_axes_on_subfigure():
    fig = plt.figure()
    subfigs = fig.subfigures(2)

    ax = subfigs[0].subplots(1, squeeze=True)
    subfigs[1].subplots(1, squeeze=True)

    assert plt.gca() != ax
    plt.sca(ax)
    assert plt.gca() == ax


def test_pylab_integration():
    IPython = pytest.importorskip("IPython")
    mpl.testing.subprocess_run_helper(
        IPython.start_ipython,
        "--pylab",
        "-c",
        ";".join((
            "import matplotlib.pyplot as plt",
            "assert plt._REPL_DISPLAYHOOK == plt._ReplDisplayHook.IPYTHON",
        )),
        timeout=60,
    )


def test_doc_pyplot_summary():
    """Test that pyplot_summary lists all the plot functions."""
    pyplot_docs = Path(__file__).parent / '../../../doc/api/pyplot_summary.rst'
    if not pyplot_docs.exists():
        pytest.skip("Documentation sources not available")

    def extract_documented_functions(lines):
        """
        Return a list of all the functions that are mentioned in the
        autosummary blocks contained in *lines*.

        An autosummary block looks like this::

            .. autosummary::
               :toctree: _as_gen
               :template: autosummary.rst
               :nosignatures:

               plot
               errorbar

        """
        functions = []
        in_autosummary = False
        for line in lines:
            if not in_autosummary:
                if line.startswith(".. autosummary::"):
                    in_autosummary = True
            else:
                if not line or line.startswith("   :"):
                    # empty line or autosummary parameter
                    continue
                if not line[0].isspace():
                    # no more indentation: end of autosummary block
                    in_autosummary = False
                    continue
                functions.append(line.strip())
        return functions

    lines = pyplot_docs.read_text().split("\n")
    doc_functions = set(extract_documented_functions(lines))
    plot_commands = set(plt._get_pyplot_commands())
    missing = plot_commands.difference(doc_functions)
    if missing:
        raise AssertionError(
            f"The following pyplot functions are not listed in the "
            f"documentation. Please add them to doc/api/pyplot_summary.rst: "
            f"{missing!r}")
    extra = doc_functions.difference(plot_commands)
    if extra:
        raise AssertionError(
            f"The following functions are listed in the pyplot documentation, "
            f"but they do not exist in pyplot. "
            f"Please remove them from doc/api/pyplot_summary.rst: {extra!r}")


def test_minor_ticks():
    plt.figure()
    plt.plot(np.arange(1, 10))
    tick_pos, tick_labels = plt.xticks(minor=True)
    assert np.all(tick_labels == np.array([], dtype=np.float64))
    assert tick_labels == []

    plt.yticks(ticks=[3.5, 6.5], labels=["a", "b"], minor=True)
    ax = plt.gca()
    tick_pos = ax.get_yticks(minor=True)
    tick_labels = ax.get_yticklabels(minor=True)
    assert np.all(tick_pos == np.array([3.5, 6.5]))
    assert [l.get_text() for l in tick_labels] == ['a', 'b']


def test_switch_backend_no_close():
    plt.switch_backend('agg')
    fig = plt.figure()
    fig = plt.figure()
    assert len(plt.get_fignums()) == 2
    plt.switch_backend('agg')
    assert len(plt.get_fignums()) == 2
    plt.switch_backend('svg')
    assert len(plt.get_fignums()) == 2


def figure_hook_example(figure):
    figure._test_was_here = True


def test_figure_hook():

    test_rc = {
        'figure.hooks': ['matplotlib.tests.test_pyplot:figure_hook_example']
    }
    with mpl.rc_context(test_rc):
        fig = plt.figure()

    assert fig._test_was_here


def test_multiple_same_figure_calls():
    fig = plt.figure(1, figsize=(1, 2))
    with pytest.warns(UserWarning, match="Ignoring specified arguments in this call"):
        fig2 = plt.figure(1, figsize=np.array([3, 4]))
    with pytest.warns(UserWarning, match="Ignoring specified arguments in this call"):
        plt.figure(fig, figsize=np.array([5, 6]))
    assert fig is fig2
    fig3 = plt.figure(1)  # Checks for false warnings
    assert fig is fig3


def test_register_existing_figure_with_pyplot():
    from matplotlib.figure import Figure
    # start with a standalone figure
    fig = Figure()
    assert fig.canvas.manager is None
    with pytest.raises(AttributeError):
        # Heads-up: This will change to returning None in the future
        # See docstring for the Figure.number property
        fig.number
    # register the Figure with pyplot
    plt.figure(fig)
    assert fig.number == 1
    # the figure can now be used in pyplot
    plt.suptitle("my title")
    assert fig.get_suptitle() == "my title"
    # it also has a manager that is properly wired up in the pyplot state
    assert plt._pylab_helpers.Gcf.get_fig_manager(fig.number) is fig.canvas.manager
    # and we can regularly switch the pyplot state
    fig2 = plt.figure()
    assert fig2.number == 2
    assert plt.figure(1) is fig
    assert plt.gcf() is fig


def test_close_all_warning():
    fig1 = plt.figure()

    # Check that the warning is issued when 'all' is passed to plt.figure
    with pytest.warns(UserWarning, match="closes all existing figures"):
        fig2 = plt.figure("all")


def test_matshow():
    fig = plt.figure()
    arr = [[0, 1], [1, 2]]

    # Smoke test that matshow does not ask for a new figsize on the existing figure
    plt.matshow(arr, fignum=fig.number)


def assert_same_signature(func1, func2):
    """
    Assert that `func1` and `func2` have the same arguments,
    i.e. same parameter count, names and kinds.

    :param func1: First function to check
    :param func2: Second function to check
    """
    params1 = inspect.signature(func1).parameters
    params2 = inspect.signature(func2).parameters

    assert len(params1) == len(params2)
    assert all([
        params1[p].name == params2[p].name and
        params1[p].kind == params2[p].kind
        for p in params1
    ])


def test_setloglevel_signature():
    assert_same_signature(plt.set_loglevel, mpl.set_loglevel)


def test_subplots_reuse_existing_figure_error():
    """Test interaction of plt.subplots(num=...) with existing figures."""
    # Create a figure with a specific number first.
    fig = plt.figure(1)

    # Case 1: Reusing without clear=True should raise ValueError
    with pytest.raises(ValueError, match="already exists"):
        plt.subplots(num=1)

    # Case 2: Reusing WITH clear=True should work fine (no error)
    fig_new, axs = plt.subplots(num=1, clear=True)
    assert fig_new is fig

    # Case 3: Test passing the actual Figure object (The "Narrow Check")
    with pytest.raises(ValueError, match="cannot be a FigureBase instance"):
        plt.subplots(num=fig)

    plt.close(1)


def test_subplot_mosaic_reuse_existing_figure_error():
    """Test that plt.subplot_mosaic raises ValueError when reusing a figure."""
    fig = plt.figure(2)

    # 1. Test passing the existing figure number
    with pytest.raises(ValueError, match="already exists"):
        plt.subplot_mosaic([['A']], num=2)

    # 2. Test passing the actual Figure object
    with pytest.raises(ValueError, match="cannot be a FigureBase instance"):
        plt.subplot_mosaic([['A']], num=fig)

    # 3. Test that clear=True allows reuse without error
    fig_new, axd = plt.subplot_mosaic([['A']], num=2, clear=True)
    assert fig_new is fig

    plt.close(2)

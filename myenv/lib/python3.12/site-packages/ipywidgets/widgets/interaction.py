# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""Interact with functions using widgets."""

from collections.abc import Iterable, Mapping
from inspect import signature, Parameter
from inspect import getcallargs
from inspect import getfullargspec as check_argspec
import sys

from IPython import get_ipython
from . import (Widget, ValueWidget, Text,
    FloatSlider, IntSlider, Checkbox, Dropdown,
    VBox, Button, DOMWidget, Output)
from IPython.display import display, clear_output
from traitlets import HasTraits, Any, Unicode, observe
from numbers import Real, Integral
from warnings import warn



empty = Parameter.empty


def show_inline_matplotlib_plots():
    """Show matplotlib plots immediately if using the inline backend.

    With ipywidgets 6.0, matplotlib plots don't work well with interact when
    using the inline backend that comes with ipykernel. Basically, the inline
    backend only shows the plot after the entire cell executes, which does not
    play well with drawing plots inside of an interact function. See
    https://github.com/jupyter-widgets/ipywidgets/issues/1181/ and
    https://github.com/ipython/ipython/issues/10376 for more details. This
    function displays any matplotlib plots if the backend is the inline backend.
    """
    if 'matplotlib' not in sys.modules:
        # matplotlib hasn't been imported, nothing to do.
        return

    try:
        import matplotlib as mpl
        from ipykernel.pylab.backend_inline import flush_figures
    except ImportError:
        return

    if (mpl.get_backend() == 'module://ipykernel.pylab.backend_inline' or
        mpl.get_backend() == 'module://matplotlib_inline.backend_inline'):
        flush_figures()


def interactive_output(f, controls):
    """Connect widget controls to a function.

    This function does not generate a user interface for the widgets (unlike `interact`).
    This enables customisation of the widget user interface layout.
    The user interface layout must be defined and displayed manually.
    """

    out = Output()
    def observer(change):
        kwargs = {k:v.value for k,v in controls.items()}
        show_inline_matplotlib_plots()
        with out:
            clear_output(wait=True)
            f(**kwargs)
            show_inline_matplotlib_plots()
    for k,w in controls.items():
        w.observe(observer, 'value')
    show_inline_matplotlib_plots()
    observer(None)
    return out


def _matches(o, pattern):
    """Match a pattern of types in a sequence."""
    if not len(o) == len(pattern):
        return False
    comps = zip(o,pattern)
    return all(isinstance(obj,kind) for obj,kind in comps)


def _get_min_max_value(min, max, value=None, step=None):
    """Return min, max, value given input values with possible None."""
    # Either min and max need to be given, or value needs to be given
    if value is None:
        if min is None or max is None:
            raise ValueError('unable to infer range, value from: ({}, {}, {})'.format(min, max, value))
        diff = max - min
        value = min + (diff / 2)
        # Ensure that value has the same type as diff
        if not isinstance(value, type(diff)):
            value = min + (diff // 2)
    else:  # value is not None
        if not isinstance(value, Real):
            raise TypeError('expected a real number, got: %r' % value)
        # Infer min/max from value
        if value == 0:
            # This gives (0, 1) of the correct type
            vrange = (value, value + 1)
        elif value > 0:
            vrange = (-value, 3*value)
        else:
            vrange = (3*value, -value)
        if min is None:
            min = vrange[0]
        if max is None:
            max = vrange[1]
    if step is not None:
        # ensure value is on a step
        tick = int((value - min) / step)
        value = min + tick * step
    if not min <= value <= max:
        raise ValueError('value must be between min and max (min={}, value={}, max={})'.format(min, value, max))
    return min, max, value

def _yield_abbreviations_for_parameter(param, kwargs):
    """Get an abbreviation for a function parameter."""
    name = param.name
    kind = param.kind
    default = param.default
    not_found = (name, empty, empty)
    if kind in (Parameter.POSITIONAL_OR_KEYWORD, Parameter.KEYWORD_ONLY):
        if name in kwargs:
            value = kwargs.pop(name)
        elif default is not empty:
            value = default
        else:
            yield not_found
        yield (name, value, default)
    elif kind == Parameter.VAR_KEYWORD:
        # In this case name=kwargs and we yield the items in kwargs with their keys.
        for k, v in kwargs.copy().items():
            kwargs.pop(k)
            yield k, v, empty


class interactive(VBox):
    """
    A VBox container containing a group of interactive widgets tied to a
    function.

    Parameters
    ----------
    __interact_f : function
        The function to which the interactive widgets are tied. The `**kwargs`
        should match the function signature.
    __options : dict
        A dict of options. Currently, the only supported keys are
        ``"manual"`` (defaults to ``False``), ``"manual_name"`` (defaults
        to ``"Run Interact"``) and ``"auto_display"`` (defaults to ``False``).
    **kwargs : various, optional
        An interactive widget is created for each keyword argument that is a
        valid widget abbreviation.

    Note that the first two parameters intentionally start with a double
    underscore to avoid being mixed up with keyword arguments passed by
    ``**kwargs``.
    """
    def __init__(self, __interact_f, __options={}, **kwargs):
        VBox.__init__(self, _dom_classes=['widget-interact'])
        self.result = None
        self.args = []
        self.kwargs = {}

        self.f = f = __interact_f
        self.clear_output = kwargs.pop('clear_output', True)
        self.manual = __options.get("manual", False)
        self.manual_name = __options.get("manual_name", "Run Interact")
        self.auto_display = __options.get("auto_display", False)

        new_kwargs = self.find_abbreviations(kwargs)
        # Before we proceed, let's make sure that the user has passed a set of args+kwargs
        # that will lead to a valid call of the function. This protects against unspecified
        # and doubly-specified arguments.
        try:
            check_argspec(f)
        except TypeError:
            # if we can't inspect, we can't validate
            pass
        else:
            getcallargs(f, **{n:v for n,v,_ in new_kwargs})
        # Now build the widgets from the abbreviations.
        self.kwargs_widgets = self.widgets_from_abbreviations(new_kwargs)

        # This has to be done as an assignment, not using self.children.append,
        # so that traitlets notices the update. We skip any objects (such as fixed) that
        # are not DOMWidgets.
        c = [w for w in self.kwargs_widgets if isinstance(w, DOMWidget)]

        # If we are only to run the function on demand, add a button to request this.
        if self.manual:
            self.manual_button = Button(description=self.manual_name)
            c.append(self.manual_button)

        self.out = Output()
        c.append(self.out)
        self.children = c

        # Wire up the widgets
        # If we are doing manual running, the callback is only triggered by the button
        # Otherwise, it is triggered for every trait change received
        # On-demand running also suppresses running the function with the initial parameters
        if self.manual:
            self.manual_button.on_click(self.update)

            # Also register input handlers on text areas, so the user can hit return to
            # invoke execution.
            for w in self.kwargs_widgets:
                if isinstance(w, Text):
                    w.continuous_update = False
                    w.observe(self.update, names='value')
        else:
            for widget in self.kwargs_widgets:
                widget.observe(self.update, names='value')
            self.update()

    # Callback function
    def update(self, *args):
        """
        Call the interact function and update the output widget with
        the result of the function call.

        Parameters
        ----------
        *args : ignored
            Required for this method to be used as traitlets callback.
        """
        self.kwargs = {}
        if self.manual:
            self.manual_button.disabled = True
        try:
            show_inline_matplotlib_plots()
            with self.out:
                if self.clear_output:
                    clear_output(wait=True)
                for widget in self.kwargs_widgets:
                    value = widget.get_interact_value()
                    self.kwargs[widget._kwarg] = value
                self.result = self.f(**self.kwargs)
                show_inline_matplotlib_plots()
                if self.auto_display and self.result is not None:
                    display(self.result)
        except Exception as e:
            ip = get_ipython()
            if ip is None:
                self.log.warning("Exception in interact callback: %s", e, exc_info=True)
            else:
                ip.showtraceback()
        finally:
            if self.manual:
                self.manual_button.disabled = False

    # Find abbreviations
    def signature(self):
        return signature(self.f)

    def find_abbreviations(self, kwargs):
        """Find the abbreviations for the given function and kwargs.
        Return (name, abbrev, default) tuples.
        """
        new_kwargs = []
        try:
            sig = self.signature()
        except (ValueError, TypeError):
            # can't inspect, no info from function; only use kwargs
            return [ (key, value, value) for key, value in kwargs.items() ]

        for param in sig.parameters.values():
            for name, value, default in _yield_abbreviations_for_parameter(param, kwargs):
                if value is empty:
                    raise ValueError('cannot find widget or abbreviation for argument: {!r}'.format(name))
                new_kwargs.append((name, value, default))
        return new_kwargs

    # Abbreviations to widgets
    def widgets_from_abbreviations(self, seq):
        """Given a sequence of (name, abbrev, default) tuples, return a sequence of Widgets."""
        result = []
        for name, abbrev, default in seq:
            if isinstance(abbrev, Widget) and (not isinstance(abbrev, ValueWidget)):
                raise TypeError("{!r} is not a ValueWidget".format(abbrev))
            widget = self.widget_from_abbrev(abbrev, default)
            if widget is None:
                raise ValueError("{!r} cannot be transformed to a widget".format(abbrev))
            if not hasattr(widget, "description") or not widget.description:
                widget.description = name
            widget._kwarg = name
            result.append(widget)
        return result

    @classmethod
    def widget_from_abbrev(cls, abbrev, default=empty):
        """Build a ValueWidget instance given an abbreviation or Widget."""
        if isinstance(abbrev, ValueWidget) or isinstance(abbrev, fixed):
            return abbrev

        if isinstance(abbrev, tuple):
            widget = cls.widget_from_tuple(abbrev)
            if default is not empty:
                try:
                    widget.value = default
                except Exception:
                    # ignore failure to set default
                    pass
            return widget

        # Try single value
        widget = cls.widget_from_single_value(abbrev)
        if widget is not None:
            return widget

        # Something iterable (list, dict, generator, ...). Note that str and
        # tuple should be handled before, that is why we check this case last.
        if isinstance(abbrev, Iterable):
            widget = cls.widget_from_iterable(abbrev)
            if default is not empty:
                try:
                    widget.value = default
                except Exception:
                    # ignore failure to set default
                    pass
            return widget

        # No idea...
        return None

    @staticmethod
    def widget_from_single_value(o):
        """Make widgets from single values, which can be used as parameter defaults."""
        if isinstance(o, str):
            return Text(value=str(o))
        elif isinstance(o, bool):
            return Checkbox(value=o)
        elif isinstance(o, Integral):
            min, max, value = _get_min_max_value(None, None, o)
            return IntSlider(value=o, min=min, max=max)
        elif isinstance(o, Real):
            min, max, value = _get_min_max_value(None, None, o)
            return FloatSlider(value=o, min=min, max=max)
        else:
            return None

    @staticmethod
    def widget_from_tuple(o):
        """Make widgets from a tuple abbreviation."""
        if _matches(o, (Real, Real)):
            min, max, value = _get_min_max_value(o[0], o[1])
            if all(isinstance(_, Integral) for _ in o):
                cls = IntSlider
            else:
                cls = FloatSlider
            return cls(value=value, min=min, max=max)
        elif _matches(o, (Real, Real, Real)):
            step = o[2]
            if step <= 0:
                raise ValueError("step must be >= 0, not %r" % step)
            min, max, value = _get_min_max_value(o[0], o[1], step=step)
            if all(isinstance(_, Integral) for _ in o):
                cls = IntSlider
            else:
                cls = FloatSlider
            return cls(value=value, min=min, max=max, step=step)

    @staticmethod
    def widget_from_iterable(o):
        """Make widgets from an iterable. This should not be done for
        a string or tuple."""
        # Dropdown expects a dict or list, so we convert an arbitrary
        # iterable to either of those.
        if isinstance(o, (list, dict)):
            return Dropdown(options=o)
        elif isinstance(o, Mapping):
            return Dropdown(options=list(o.items()))
        else:
            return Dropdown(options=list(o))

    # Return a factory for interactive functions
    @classmethod
    def factory(cls):
        options = dict(manual=False, auto_display=True, manual_name="Run Interact")
        return _InteractFactory(cls, options)


class _InteractFactory:
    """
    Factory for instances of :class:`interactive`.

    This class is needed to support options like::

        >>> @interact.options(manual=True)
        ... def greeting(text="World"):
        ...     print("Hello {}".format(text))

    Parameters
    ----------
    cls : class
        The subclass of :class:`interactive` to construct.
    options : dict
        A dict of options used to construct the interactive
        function. By default, this is returned by
        ``cls.default_options()``.
    kwargs : dict
        A dict of **kwargs to use for widgets.
    """
    def __init__(self, cls, options, kwargs={}):
        self.cls = cls
        self.opts = options
        self.kwargs = kwargs

    def widget(self, f):
        """
        Return an interactive function widget for the given function.

        The widget is only constructed, not displayed nor attached to
        the function.

        Returns
        -------
        An instance of ``self.cls`` (typically :class:`interactive`).

        Parameters
        ----------
        f : function
            The function to which the interactive widgets are tied.
        """
        return self.cls(f, self.opts, **self.kwargs)

    def __call__(self, __interact_f=None, **kwargs):
        """
        Make the given function interactive by adding and displaying
        the corresponding :class:`interactive` widget.

        Expects the first argument to be a function. Parameters to this
        function are widget abbreviations passed in as keyword arguments
        (``**kwargs``). Can be used as a decorator (see examples).

        Returns
        -------
        f : __interact_f with interactive widget attached to it.

        Parameters
        ----------
        __interact_f : function
            The function to which the interactive widgets are tied. The `**kwargs`
            should match the function signature. Passed to :func:`interactive()`
        **kwargs : various, optional
            An interactive widget is created for each keyword argument that is a
            valid widget abbreviation. Passed to :func:`interactive()`

        Examples
        --------
        Render an interactive text field that shows the greeting with the passed in
        text::

            # 1. Using interact as a function
            def greeting(text="World"):
                print("Hello {}".format(text))
            interact(greeting, text="Jupyter Widgets")

            # 2. Using interact as a decorator
            @interact
            def greeting(text="World"):
                print("Hello {}".format(text))

            # 3. Using interact as a decorator with named parameters
            @interact(text="Jupyter Widgets")
            def greeting(text="World"):
                print("Hello {}".format(text))

        Render an interactive slider widget and prints square of number::

            # 1. Using interact as a function
            def square(num=1):
                print("{} squared is {}".format(num, num*num))
            interact(square, num=5)

            # 2. Using interact as a decorator
            @interact
            def square(num=2):
                print("{} squared is {}".format(num, num*num))

            # 3. Using interact as a decorator with named parameters
            @interact(num=5)
            def square(num=2):
                print("{} squared is {}".format(num, num*num))
        """
        # If kwargs are given, replace self by a new
        # _InteractFactory with the updated kwargs
        if kwargs:
            kw = dict(self.kwargs)
            kw.update(kwargs)
            self = type(self)(self.cls, self.opts, kw)

        f = __interact_f
        if f is None:
            # This branch handles the case 3
            # @interact(a=30, b=40)
            # def f(*args, **kwargs):
            #     ...
            #
            # Simply return the new factory
            return self

        # positional arg support in: https://gist.github.com/8851331
        # Handle the cases 1 and 2
        # 1. interact(f, **kwargs)
        # 2. @interact
        #    def f(*args, **kwargs):
        #        ...
        w = self.widget(f)
        try:
            f.widget = w
        except AttributeError:
            # some things (instancemethods) can't have attributes attached,
            # so wrap in a lambda
            f = lambda *args, **kwargs: __interact_f(*args, **kwargs)
            f.widget = w
        show_inline_matplotlib_plots()
        display(w)
        return f

    def options(self, **kwds):
        """
        Change options for interactive functions.

        Returns
        -------
        A new :class:`_InteractFactory` which will apply the
        options when called.
        """
        opts = dict(self.opts)
        for k in kwds:
            try:
                # Ensure that the key exists because we want to change
                # existing options, not add new ones.
                _ = opts[k]
            except KeyError:
                raise ValueError("invalid option {!r}".format(k))
            opts[k] = kwds[k]
        return type(self)(self.cls, opts, self.kwargs)


interact = interactive.factory()
interact_manual = interact.options(manual=True, manual_name="Run Interact")


class fixed(HasTraits):
    """A pseudo-widget whose value is fixed and never synced to the client."""
    value = Any(help="Any Python object")
    description = Unicode('', help="Any Python object")
    def __init__(self, value, **kwargs):
        super().__init__(value=value, **kwargs)
    def get_interact_value(self):
        """Return the value for this widget which should be passed to
        interactive functions. Custom widgets can change this method
        to process the raw value ``self.value``.
        """
        return self.value

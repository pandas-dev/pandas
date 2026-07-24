from matplotlib import cbook
from matplotlib.artist import Artist


class Container(tuple):
    """
    Base class for containers.

    Containers are classes that collect semantically related Artists such as
    the bars of a bar plot.
    """

    def __repr__(self):
        return f"<{type(self).__name__} object of {len(self)} artists>"

    def __new__(cls, *args, **kwargs):
        return tuple.__new__(cls, args[0])

    def __init__(self, kl, label=None):
        self._callbacks = cbook.CallbackRegistry(signals=["pchanged"])
        self._remove_method = None
        self._label = str(label) if label is not None else None

    def remove(self):
        for c in cbook.flatten(
                self, scalarp=lambda x: isinstance(x, Artist)):
            if c is not None:
                c.remove()
        if self._remove_method:
            self._remove_method(self)

    def get_children(self):
        return [child for child in cbook.flatten(self) if child is not None]

    get_label = Artist.get_label
    set_label = Artist.set_label
    add_callback = Artist.add_callback
    remove_callback = Artist.remove_callback
    pchanged = Artist.pchanged


class BarContainer(Container):
    """
    Container for the artists of bar plots (e.g. created by `.Axes.bar`).

    The container can be treated as a tuple of the *patches* themselves.
    Additionally, you can access these and further parameters by the
    attributes.

    Attributes
    ----------
    patches : list of :class:`~matplotlib.patches.Rectangle`
        The artists of the bars.

    errorbar : None or :class:`~matplotlib.container.ErrorbarContainer`
        A container for the error bar artists if error bars are present.
        *None* otherwise.

    datavalues : None or array-like
        The underlying data values corresponding to the bars.

    orientation : {'vertical', 'horizontal'}, default: None
        If 'vertical', the bars are assumed to be vertical.
        If 'horizontal', the bars are assumed to be horizontal.

    """

    def __init__(self, patches, errorbar=None, *, datavalues=None,
                 orientation=None, **kwargs):
        self.patches = patches
        self.errorbar = errorbar
        self.datavalues = datavalues
        self.orientation = orientation
        super().__init__(patches, **kwargs)

    @property
    def bottoms(self):
        """
        Return the values at the lower end of the bars.

        .. versionadded:: 3.11
        """
        if self.orientation == 'vertical':
            return [p.get_y() for p in self.patches]
        elif self.orientation == 'horizontal':
            return [p.get_x() for p in self.patches]
        else:
            raise ValueError("orientation must be 'vertical' or 'horizontal'.")

    @property
    def tops(self):
        """
        Return the values at the upper end of the bars.

        .. versionadded:: 3.11
        """
        if self.orientation == 'vertical':
            return [p.get_y() + p.get_height() for p in self.patches]
        elif self.orientation == 'horizontal':
            return [p.get_x() + p.get_width() for p in self.patches]
        else:
            raise ValueError("orientation must be 'vertical' or 'horizontal'.")

    @property
    def position_centers(self):
        """
        Return the centers of bar positions.

        .. versionadded:: 3.11
        """
        if self.orientation == 'vertical':
            return [p.get_x() + p.get_width() / 2 for p in self.patches]
        elif self.orientation == 'horizontal':
            return [p.get_y() + p.get_height() / 2 for p in self.patches]
        else:
            raise ValueError("orientation must be 'vertical' or 'horizontal'.")


class ErrorbarContainer(Container):
    """
    Container for the artists of error bars (e.g. created by `.Axes.errorbar`).

    The container can be treated as the *lines* tuple itself.
    Additionally, you can access these and further parameters by the
    attributes.

    Attributes
    ----------
    lines : tuple
        Tuple of ``(data_line, caplines, barlinecols)``.

        - data_line : A `~matplotlib.lines.Line2D` instance of x, y plot markers
          and/or line.
        - caplines : A tuple of `~matplotlib.lines.Line2D` instances of the error
          bar caps.
        - barlinecols : A tuple of `~matplotlib.collections.LineCollection` with the
          horizontal and vertical error ranges.

    has_xerr, has_yerr : bool
        ``True`` if the errorbar has x/y errors.

    """

    def __init__(self, lines, has_xerr=False, has_yerr=False, **kwargs):
        self.lines = lines
        self.has_xerr = has_xerr
        self.has_yerr = has_yerr
        super().__init__(lines, **kwargs)


class PieContainer:
    """
    Container for the artists of pie charts (e.g. created by `.Axes.pie`).

    .. versionadded:: 3.11

    .. warning::
        The class name ``PieContainer`` name is provisional and may change in future
        to reflect development of its functionality.

    You can access the wedge patches and further parameters by the attributes.

    Attributes
    ----------
    wedges : list of `~matplotlib.patches.Wedge`
        The artists of the pie wedges.

    values : `numpy.ndarray`
        The data that the pie is based on.

    fracs : `numpy.ndarray`
        The fraction of the pie that each wedge represents.

    texts : list of list of `~matplotlib.text.Text`
        The artists of any labels on the pie wedges.  Each inner list has one
        text label per wedge.

    """
    def __init__(self, wedges, values, normalize):
        self.wedges = wedges
        self._texts = []
        self._values = values
        self._normalize = normalize

    @property
    def texts(self):
        # Only return non-empty sublists.  An empty sublist may have been added
        # for backwards compatibility of the Axes.pie return value (see __getitem__).
        return [t_list for t_list in self._texts if t_list]

    @property
    def values(self):
        result = self._values.copy()
        result.flags.writeable = False
        return result

    @property
    def fracs(self):
        if self._normalize:
            result = self._values / self._values.sum()
        else:
            result = self._values

        result.flags.writeable = False
        return result

    def add_texts(self, texts):
        """Add a list of `~matplotlib.text.Text` objects to the container."""
        self._texts.append(texts)

    def remove(self):
        """Remove all wedges and texts from the axes"""
        for artist_list in self.wedges, self._texts:
            for artist in cbook.flatten(artist_list):
                artist.remove()

    def __getitem__(self, key):
        # needed to support unpacking into a tuple for backward compatibility of the
        # Axes.pie return value
        return (self.wedges, *self._texts)[key]


class StemContainer(Container):
    """
    Container for the artists created in a :meth:`.Axes.stem` plot.

    The container can be treated like a namedtuple ``(markerline, stemlines,
    baseline)``.

    Attributes
    ----------
    markerline : `~matplotlib.lines.Line2D`
        The artist of the markers at the stem heads.

    stemlines : `~matplotlib.collections.LineCollection`
        The artists of the vertical lines for all stems.

    baseline : `~matplotlib.lines.Line2D`
        The artist of the horizontal baseline.
    """
    def __init__(self, markerline_stemlines_baseline, **kwargs):
        """
        Parameters
        ----------
        markerline_stemlines_baseline : tuple
            Tuple of ``(markerline, stemlines, baseline)``.
            ``markerline`` contains the `.Line2D` of the markers,
            ``stemlines`` is a `.LineCollection` of the main lines,
            ``baseline`` is the `.Line2D` of the baseline.
        """
        markerline, stemlines, baseline = markerline_stemlines_baseline
        self.markerline = markerline
        self.stemlines = stemlines
        self.baseline = baseline
        super().__init__(markerline_stemlines_baseline, **kwargs)

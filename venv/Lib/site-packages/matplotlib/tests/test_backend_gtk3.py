import os
from matplotlib import pyplot as plt

import pytest
from unittest import mock


@pytest.mark.backend("gtk3agg", skip_on_importerror=True)
def test_save_figure_return():
    from gi.repository import Gtk
    fig, ax = plt.subplots()
    ax.imshow([[1]])
    with mock.patch("gi.repository.Gtk.FileFilter") as fileFilter:
        filt = fileFilter.return_value
        filt.get_name.return_value = "Portable Network Graphics"
        with mock.patch("gi.repository.Gtk.FileChooserDialog") as dialogChooser:
            dialog = dialogChooser.return_value
            dialog.get_filter.return_value = filt
            dialog.get_filename.return_value = "foobar.png"
            dialog.run.return_value = Gtk.ResponseType.OK
            fname = fig.canvas.manager.toolbar.save_figure()
            os.remove("foobar.png")
            assert fname == "foobar.png"

            with mock.patch("gi.repository.Gtk.MessageDialog"):
                dialog.get_filename.return_value = None
                dialog.run.return_value = Gtk.ResponseType.OK
                fname = fig.canvas.manager.toolbar.save_figure()
                assert fname is None

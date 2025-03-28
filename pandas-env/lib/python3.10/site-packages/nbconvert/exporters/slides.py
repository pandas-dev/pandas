"""HTML slide show Exporter class"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from copy import deepcopy
from warnings import warn

from traitlets import Bool, Unicode, default

from nbconvert.preprocessors.base import Preprocessor

from .html import HTMLExporter


class _RevealMetadataPreprocessor(Preprocessor):
    # A custom preprocessor adding convenience metadata to cells

    def preprocess(self, nb, resources=None):
        nb = deepcopy(nb)

        for cell in nb.cells:
            # Make sure every cell has a slide_type
            try:
                slide_type = cell.metadata.get("slideshow", {}).get("slide_type", "-")
            except AttributeError:
                slide_type = "-"
            cell.metadata.slide_type = slide_type

        # Find the first visible cell
        for index, cell in enumerate(nb.cells):
            if cell.metadata.slide_type not in {"notes", "skip"}:
                cell.metadata.slide_type = "slide"
                cell.metadata.slide_start = True
                cell.metadata.subslide_start = True
                first_slide_ix = index
                break
        else:
            msg = "All cells are hidden, cannot create slideshow"
            raise ValueError(msg)

        in_fragment = False

        for index, cell in enumerate(nb.cells[first_slide_ix + 1 :], start=(first_slide_ix + 1)):
            previous_cell = nb.cells[index - 1]

            # Slides are <section> elements in the HTML, subslides (the vertically
            # stacked slides) are also <section> elements inside the slides,
            # and fragments are <div>s within subslides. Subslide and fragment
            # elements can contain content:
            # <section>
            #   <section>
            #     (content)
            #     <div class="fragment">(content)</div>
            #   </section>
            # </section>

            # Get the slide type. If type is subslide or slide,
            # end the last slide/subslide/fragment as applicable.
            if cell.metadata.slide_type == "slide":
                previous_cell.metadata.slide_end = True
                cell.metadata.slide_start = True
            if cell.metadata.slide_type in {"subslide", "slide"}:
                previous_cell.metadata.fragment_end = in_fragment
                previous_cell.metadata.subslide_end = True
                cell.metadata.subslide_start = True
                in_fragment = False

            elif cell.metadata.slide_type == "fragment":
                cell.metadata.fragment_start = True
                if in_fragment:
                    previous_cell.metadata.fragment_end = True
                else:
                    in_fragment = True

        # The last cell will always be the end of a slide
        nb.cells[-1].metadata.fragment_end = in_fragment
        nb.cells[-1].metadata.subslide_end = True
        nb.cells[-1].metadata.slide_end = True

        return nb, resources


class SlidesExporter(HTMLExporter):
    """Exports HTML slides with reveal.js"""

    # Overrides from HTMLExporter
    #################################
    export_from_notebook = "Reveal.js slides"

    @default("template_name")
    def _template_name_default(self):
        return "reveal"

    @default("file_extension")
    def _file_extension_default(self):
        return ".slides.html"

    @default("template_extension")
    def _template_extension_default(self):
        return ".html.j2"

    # Extra resources
    #################################
    reveal_url_prefix = Unicode(
        help="""The URL prefix for reveal.js (version 3.x).
        This defaults to the reveal CDN, but can be any url pointing to a copy
        of reveal.js.

        For speaker notes to work, this must be a relative path to a local
        copy of reveal.js: e.g., "reveal.js".

        If a relative path is given, it must be a subdirectory of the
        current directory (from which the server is run).

        See the usage documentation
        (https://nbconvert.readthedocs.io/en/latest/usage.html#reveal-js-html-slideshow)
        for more details.
        """
    ).tag(config=True)

    @default("reveal_url_prefix")
    def _reveal_url_prefix_default(self):
        if "RevealHelpPreprocessor.url_prefix" in self.config:
            warn(
                "Please update RevealHelpPreprocessor.url_prefix to "
                "SlidesExporter.reveal_url_prefix in config files.",
                stacklevel=2,
            )
            return self.config.RevealHelpPreprocessor.url_prefix
        return "https://unpkg.com/reveal.js@4.0.2"

    reveal_theme = Unicode(
        "simple",
        help="""
        Name of the reveal.js theme to use.

        We look for a file with this name under
        ``reveal_url_prefix``/css/theme/``reveal_theme``.css.

        https://github.com/hakimel/reveal.js/tree/master/css/theme has
        list of themes that ship by default with reveal.js.
        """,
    ).tag(config=True)

    reveal_transition = Unicode(
        "slide",
        help="""
        Name of the reveal.js transition to use.

        The list of transitions that ships by default with reveal.js are:
        none, fade, slide, convex, concave and zoom.
        """,
    ).tag(config=True)

    reveal_scroll = Bool(
        False,
        help="""
        If True, enable scrolling within each slide
        """,
    ).tag(config=True)

    reveal_number = Unicode(
        "",
        help="""
        slide number format (e.g. 'c/t'). Choose from:
        'c': current, 't': total, 'h': horizontal, 'v': vertical
        """,
    ).tag(config=True)

    reveal_width = Unicode(
        "",
        help="""
        width used to determine the aspect ratio of your presentation.
        Use the horizontal pixels available on your intended presentation
        equipment.
        """,
    ).tag(config=True)

    reveal_height = Unicode(
        "",
        help="""
        height used to determine the aspect ratio of your presentation.
        Use the horizontal pixels available on your intended presentation
        equipment.
        """,
    ).tag(config=True)

    font_awesome_url = Unicode(
        "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.css",
        help="""
        URL to load font awesome from.

        Defaults to loading from cdnjs.
        """,
    ).tag(config=True)

    def _init_resources(self, resources):
        resources = super()._init_resources(resources)
        if "reveal" not in resources:
            resources["reveal"] = {}
        resources["reveal"]["url_prefix"] = self.reveal_url_prefix
        resources["reveal"]["theme"] = self.reveal_theme
        resources["reveal"]["transition"] = self.reveal_transition
        resources["reveal"]["scroll"] = self.reveal_scroll
        resources["reveal"]["number"] = self.reveal_number
        resources["reveal"]["height"] = self.reveal_height
        resources["reveal"]["width"] = self.reveal_width
        return resources

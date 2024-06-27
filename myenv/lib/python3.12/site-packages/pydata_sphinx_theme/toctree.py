"""Methods to build the toctree used in the html pages."""

from functools import lru_cache
from typing import List, Union
from urllib.parse import urlparse

import sphinx
from bs4 import BeautifulSoup
from docutils import nodes
from docutils.nodes import Node
from sphinx import addnodes
from sphinx.addnodes import toctree as toctree_node
from sphinx.application import Sphinx
from sphinx.environment.adapters.toctree import TocTree
from sphinx.locale import _

from .utils import traverse_or_findall


def add_inline_math(node: Node) -> str:
    """Render a node with HTML tags that activate MathJax processing.

    This is meant for use with rendering section titles with math in them, because
    math outputs are ignored by pydata-sphinx-theme's header.

    related to the behaviour of a normal math node from:
    https://github.com/sphinx-doc/sphinx/blob/master/sphinx/ext/mathjax.py#L28
    """
    return (
        '<span class="math notranslate nohighlight">' rf"\({node.astext()}\)" "</span>"
    )


def add_toctree_functions(
    app: Sphinx, pagename: str, templatename: str, context, doctree
) -> None:
    """Add functions so Jinja templates can add toctree objects."""

    @lru_cache(maxsize=None)
    def generate_header_nav_html(
        n_links_before_dropdown: int = 5, dropdown_text: str = "More"
    ) -> str:
        """Generate top-level links that are meant for the header navigation.

        We use this function instead of the TocTree-based one used for the
        sidebar because this one is much faster for generating the links and
        we don't need the complexity of the full Sphinx TocTree.

        This includes two kinds of links:

        - Links to pages described listed in the root_doc TocTrees
        - External links defined in theme configuration

        Additionally it will create a dropdown list for several links after
        a cutoff.

        Parameters:
            n_links_before_dropdown:The number of links to show before nesting the remaining links in a Dropdown element.
            dropdown_text:Text of the dropdown element button.
        """
        try:
            n_links_before_dropdown = int(n_links_before_dropdown)
        except Exception:
            raise ValueError(
                f"n_links_before_dropdown is not an int: {n_links_before_dropdown}"
            )
        toctree = TocTree(app.env)

        # Find the active header navigation item so we decide whether to highlight
        # Will be empty if there is no active page (root_doc, or genindex etc)
        if sphinx.version_info[:2] >= (7, 2):
            from sphinx.environment.adapters.toctree import _get_toctree_ancestors

            active_header_page = [
                *_get_toctree_ancestors(app.env.toctree_includes, pagename)
            ]
        else:
            active_header_page = toctree.get_toctree_ancestors(pagename)
        if active_header_page:
            # The final list item will be the top-most ancestor
            active_header_page = active_header_page[-1]

        # Find the root document because it lists our top-level toctree pages
        root = app.env.tocs[app.config.root_doc]

        # Iterate through each toctree node in the root document
        # Grab the toctree pages and find the relative link + title.
        links_html = []
        # TODO: use `root.findall(toctree_node)` once docutils min version >=0.18.1
        for toc in traverse_or_findall(root, toctree_node):
            for title, page in toc.attributes["entries"]:
                # if the page is using "self" use the correct link
                page = toc.attributes["parent"] if page == "self" else page

                # If this is the active ancestor page, add a class so we highlight it
                current = " current active" if page == active_header_page else ""

                # sanitize page title for use in the html output if needed
                if title is None:
                    title = ""
                    for node in app.env.titles[page].children:
                        if isinstance(node, nodes.math):
                            title += add_inline_math(node)
                        else:
                            title += node.astext()

                # set up the status of the link and the path
                # if the path is relative then we use the context for the path
                # resolution and the internal class.
                # If it's an absolute one then we use the external class and
                # the complete url.
                is_absolute = bool(urlparse(page).netloc)
                link_status = "external" if is_absolute else "internal"
                link_href = page if is_absolute else context["pathto"](page)

                # create the html output
                links_html.append(
                    f"""
                    <li class="nav-item{current}">
                      <a class="nav-link nav-{link_status}" href="{link_href}">
                        {title}
                      </a>
                    </li>
                """
                )

        # Add external links defined in configuration as sibling list items
        for external_link in context["theme_external_links"]:
            links_html.append(
                f"""
                <li class="nav-item">
                  <a class="nav-link nav-external" href="{ external_link["url"] }">
                    { external_link["name"] }
                  </a>
                </li>
                """
            )

        # The first links will always be visible
        links_solo = links_html[:n_links_before_dropdown]
        out = "\n".join(links_solo)

        # Wrap the final few header items in a "more" dropdown
        links_dropdown = [
            # 🐲 brittle code, relies on the assumption that the code above
            # gives each link in the nav a `nav-link` CSS class
            html.replace("nav-link", "nav-link dropdown-item")
            for html in links_html[n_links_before_dropdown:]
        ]

        if links_dropdown:
            links_dropdown_html = "\n".join(links_dropdown)
            out += f"""
            <li class="nav-item dropdown">
                <button class="btn dropdown-toggle nav-item" type="button" data-bs-toggle="dropdown" aria-expanded="false" aria-controls="pst-header-nav-more-links">
                    {_(dropdown_text)}
                </button>
                <ul id="pst-header-nav-more-links" class="dropdown-menu">
                    {links_dropdown_html}
                </ul>
            </li>
            """

        return out

    # Cache this function because it is expensive to run, and because Sphinx
    # somehow runs this twice in some circumstances in unpredictable ways.
    @lru_cache(maxsize=None)
    def generate_toctree_html(
        kind: str, startdepth: int = 1, show_nav_level: int = 1, **kwargs
    ) -> Union[BeautifulSoup, str]:
        """Return the navigation link structure in HTML.

        This is similar to Sphinx's own default TocTree generation, but it is modified
        to generate TocTrees for *second*-level pages and below (not supported
        by default in Sphinx).
        This is used for our sidebar, which starts at the second-level page.

        It also modifies the generated TocTree slightly for Bootstrap classes
        and structure (via BeautifulSoup).

        Arguments are passed to Sphinx "toctree" function (context["toctree"] below).

        ref: https://www.sphinx-doc.org/en/master/templating.html#toctree

        Parameters:
            kind : "sidebar" or "raw". Whether to generate HTML meant for sidebar navigation ("sidebar") or to return the raw BeautifulSoup object ("raw").
            startdepth : The level of the toctree at which to start. By default, for the navbar uses the normal toctree (`startdepth=0`), and for the sidebar starts from the second level (`startdepth=1`).
            show_nav_level : The level of the navigation bar to toggle as visible on page load. By default, this level is 1, and only top-level pages are shown, with drop-boxes to reveal children. Increasing `show_nav_level` will show child levels as well.
            kwargs: passed to the Sphinx `toctree` template function.

        Returns:
            HTML string (if kind == "sidebar") OR BeautifulSoup object (if kind == "raw")
        """
        if startdepth == 0:
            toc_sphinx = context["toctree"](**kwargs)
        else:
            # select the "active" subset of the navigation tree for the sidebar
            toc_sphinx = index_toctree(app, pagename, startdepth, **kwargs)

        soup = BeautifulSoup(toc_sphinx, "html.parser")

        # pair "current" with "active" since that's what we use w/ bootstrap
        for li in soup("li", {"class": "current"}):
            li["class"].append("active")

        # Remove sidebar links to sub-headers on the page
        for li in soup.select("li"):
            # Remove
            if li.find("a"):
                href = li.find("a")["href"]
                if "#" in href and href != "#":
                    li.decompose()

        if kind == "sidebar":
            # Add bootstrap classes for first `ul` items
            for ul in soup("ul", recursive=False):
                ul.attrs["class"] = ul.attrs.get("class", []) + ["nav", "bd-sidenav"]

            # Add collapse boxes for parts/captions.
            # Wraps the TOC part in an extra <ul> to behave like chapters with toggles
            # show_nav_level: 0 means make parts collapsible.
            if show_nav_level == 0:
                partcaptions = soup.find_all("p", attrs={"class": "caption"})
                if len(partcaptions):
                    new_soup = BeautifulSoup(
                        "<ul class='list-caption'></ul>", "html.parser"
                    )
                    for caption in partcaptions:
                        # Assume that the next <ul> element is the TOC list
                        # for this part
                        for sibling in caption.next_siblings:
                            if sibling.name == "ul":
                                toclist = sibling
                                break
                        li = soup.new_tag("li", attrs={"class": "toctree-l0"})
                        li.extend([caption, toclist])
                        new_soup.ul.append(li)
                    soup = new_soup

            # Add icons and labels for collapsible nested sections
            add_collapse_checkboxes(soup)

            # Open the sidebar navigation to the proper depth
            for ii in range(int(show_nav_level)):
                for checkbox in soup.select(
                    f"li.toctree-l{ii} > input.toctree-checkbox"
                ):
                    checkbox.attrs["checked"] = None

        return soup

    @lru_cache(maxsize=None)
    def generate_toc_html(kind: str = "html") -> BeautifulSoup:
        """Return the within-page TOC links in HTML."""
        if "toc" not in context:
            return ""

        soup = BeautifulSoup(context["toc"], "html.parser")

        # Add toc-hN + visible classes
        def add_header_level_recursive(ul, level):
            if ul is None:
                return
            if level <= (context["theme_show_toc_level"] + 1):
                ul["class"] = ul.get("class", []) + ["visible"]
            for li in ul("li", recursive=False):
                li["class"] = li.get("class", []) + [f"toc-h{level}"]
                add_header_level_recursive(li.find("ul", recursive=False), level + 1)

        add_header_level_recursive(soup.find("ul"), 1)

        # Add in CSS classes for bootstrap
        for ul in soup("ul"):
            ul["class"] = ul.get("class", []) + ["nav", "section-nav", "flex-column"]

        for li in soup("li"):
            li["class"] = li.get("class", []) + ["nav-item", "toc-entry"]
            if li.find("a"):
                a = li.find("a")
                a["class"] = a.get("class", []) + ["nav-link"]

        # If we only have one h1 header, assume it's a title
        h1_headers = soup.select(".toc-h1")
        if len(h1_headers) == 1:
            title = h1_headers[0]
            # If we have no sub-headers of a title then we won't have a TOC
            if not title.select(".toc-h2"):
                out = ""
            else:
                out = title.find("ul")
        # Else treat the h1 headers as sections
        else:
            out = soup

        # Return the toctree object
        if kind == "html":
            return out
        else:
            return soup

    def navbar_align_class() -> List[str]:
        """Return the class that aligns the navbar based on config."""
        align = context.get("theme_navbar_align", "content")
        align_options = {
            "content": ("col-lg-3", "col-lg-9", "me-auto"),
            "left": ("", "", "me-auto"),
            "right": ("", "", "ms-auto"),
        }
        if align not in align_options:
            raise ValueError(
                "Theme option navbar_align must be one of"
                f"{align_options.keys()}, got: {align}"
            )
        return align_options[align]

    context["generate_header_nav_html"] = generate_header_nav_html
    context["generate_toctree_html"] = generate_toctree_html
    context["generate_toc_html"] = generate_toc_html
    context["navbar_align_class"] = navbar_align_class


def add_collapse_checkboxes(soup: BeautifulSoup) -> None:
    """Add checkboxes to collapse children in a toctree."""
    # based on https://github.com/pradyunsg/furo

    toctree_checkbox_count = 0

    for element in soup.find_all("li", recursive=True):
        # We check all "li" elements, to add a "current-page" to the correct li.
        classes = element.get("class", [])

        # expanding the parent part explicitly, if present
        if "current" in classes:
            parentli = element.find_parent("li", class_="toctree-l0")
            if parentli:
                parentli.select("p.caption ~ input")[0].attrs["checked"] = ""

        # Nothing more to do, unless this has "children"
        if not element.find("ul"):
            continue

        # Add a class to indicate that this has children.
        element["class"] = classes + ["has-children"]

        # We're gonna add a checkbox.
        toctree_checkbox_count += 1
        checkbox_name = f"toctree-checkbox-{toctree_checkbox_count}"

        # Add the "label" for the checkbox which will get filled.
        if soup.new_tag is None:
            continue

        label = soup.new_tag(
            "label", attrs={"for": checkbox_name, "class": "toctree-toggle"}
        )
        label.append(soup.new_tag("i", attrs={"class": "fa-solid fa-chevron-down"}))
        if "toctree-l0" in classes:
            # making label cover the whole caption text with css
            label["class"] = "label-parts"
        element.insert(1, label)

        # Add the checkbox that's used to store expanded/collapsed state.
        checkbox = soup.new_tag(
            "input",
            attrs={
                "type": "checkbox",
                "class": ["toctree-checkbox"],
                "id": checkbox_name,
                "name": checkbox_name,
            },
        )

        # if this has a "current" class, be expanded by default
        # (by checking the checkbox)
        if "current" in classes:
            checkbox.attrs["checked"] = ""

        element.insert(1, checkbox)


def get_local_toctree_for(
    self: TocTree, indexname: str, docname: str, builder, collapse: bool, **kwargs
) -> List[BeautifulSoup]:
    """Return the "local" TOC nodetree (relative to `indexname`)."""
    # this is a copy of `TocTree.get_toctree_for`, but where the sphinx version
    # always uses the "root" doctree:
    #     doctree = self.env.get_doctree(self.env.config.root_doc)
    # we here use the `indexname` additional argument to be able to use a subset
    # of the doctree (e.g. starting at a second level for the sidebar):
    #     doctree = app.env.tocs[indexname].deepcopy()

    doctree = self.env.tocs[indexname].deepcopy()

    toctrees = []
    if "includehidden" not in kwargs:
        kwargs["includehidden"] = True
    if "maxdepth" not in kwargs or not kwargs["maxdepth"]:
        kwargs["maxdepth"] = 0
    else:
        kwargs["maxdepth"] = int(kwargs["maxdepth"])
    kwargs["collapse"] = collapse

    # TODO: use `doctree.findall(addnodes.toctree)` once docutils min version >=0.18.1
    for toctreenode in traverse_or_findall(doctree, addnodes.toctree):
        toctree = self.resolve(docname, builder, toctreenode, prune=True, **kwargs)
        if toctree:
            toctrees.append(toctree)
    if not toctrees:
        return None
    result = toctrees[0]
    for toctree in toctrees[1:]:
        result.extend(toctree.children)
    return result


def index_toctree(
    app: Sphinx, pagename: str, startdepth: int, collapse: bool = True, **kwargs
):
    """Returns the "local" (starting at `startdepth`) TOC tree containing the current page, rendered as HTML bullet lists.

    This is the equivalent of `context["toctree"](**kwargs)` in sphinx
    templating, but using the startdepth-local instead of global TOC tree.
    """
    # this is a variant of the function stored in `context["toctree"]`, which is
    # defined as `lambda **kwargs: self._get_local_toctree(pagename, **kwargs)`
    # with `self` being the HMTLBuilder and the `_get_local_toctree` basically
    # returning:
    #     return self.render_partial(TocTree(self.env).get_toctree_for(
    #         pagename, self, collapse, **kwargs))['fragment']

    if "includehidden" not in kwargs:
        kwargs["includehidden"] = False
    if kwargs.get("maxdepth") == "":
        kwargs.pop("maxdepth")

    toctree = TocTree(app.env)
    if sphinx.version_info[:2] >= (7, 2):
        from sphinx.environment.adapters.toctree import _get_toctree_ancestors

        ancestors = [*_get_toctree_ancestors(app.env.toctree_includes, pagename)]
    else:
        ancestors = toctree.get_toctree_ancestors(pagename)
    try:
        indexname = ancestors[-startdepth]
    except IndexError:
        # eg for index.rst, but also special pages such as genindex, py-modindex, search
        # those pages don't have a "current" element in the toctree, so we can
        # directly return an empty string instead of using the default sphinx
        # toctree.get_toctree_for(pagename, app.builder, collapse, **kwargs)
        return ""

    toctree_element = get_local_toctree_for(
        toctree, indexname, pagename, app.builder, collapse, **kwargs
    )
    return app.builder.render_partial(toctree_element)["fragment"]

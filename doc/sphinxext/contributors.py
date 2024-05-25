"""Sphinx extension for listing code contributors to a release.

Usage::

   .. contributors:: v0.23.0..v0.23.1

This will be replaced with a message indicating the number of
code contributors and commits, and then list each contributor
individually. For development versions (before a tag is available)
use::

    .. contributors:: v0.23.0..v0.23.1|HEAD

While the v0.23.1 tag does not exist, that will use the HEAD of the
branch as the end of the revision range.
"""

from announce import build_components
from docutils import nodes
from docutils.parsers.rst import Directive
import git


class ContributorsDirective(Directive):
    required_arguments = 1
    name = "contributors"

    def run(self):
        """
        Execute the ContributorsDirective directive.

        This function retrieves information about code contributors and commits within
        a specified version range and generates a list of contributors to include in the
        documentation.

        Returns
        -------
        List
            List containing author message and contributors.

        Examples
        --------
        >>> directive = ContributorsDirective()
        >>> directive.arguments = ['v0.23.0..v0.23.1']
        >>> nodes = directive.run()
        >>> nodes
        [<nodes.paragraph>, <nodes.bullet_list>]
        """
        range_ = self.arguments[0]
        if range_.endswith("x..HEAD"):
            return [nodes.paragraph(), nodes.bullet_list()]
        try:
            components = build_components(range_)
        except git.GitCommandError as exc:
            return [
                self.state.document.reporter.warning(
                    f"Cannot find contributors for range {repr(range_)}: {exc}",
                    line=self.lineno,
                )
            ]
        else:
            message = nodes.paragraph()
            message += nodes.Text(components["author_message"])

            listnode = nodes.bullet_list()

            for author in components["authors"]:
                para = nodes.paragraph()
                para += nodes.Text(author)
                listnode += nodes.list_item("", para)

        return [message, listnode]


def setup(app):
    """
    Setup function for the Sphinx extension.

    This function initializes the Sphinx extension, adding the 'contributors' directive
    to the app.

    Parameters
    ----------
    app : sphinx.application.Sphinx
        The Sphinx application object.
    
    Returns
    -------
    dict
        A dictionary containing version information and safety flags for the extension.

    Examples
    --------
    >>> app = App()
    >>> setup_info = setup(app)
    >>> setup_info
    {'version': '0.1', 'parallel_read_safe': True, 'parallel_write_safe': True}
    """
    app.add_directive("contributors", ContributorsDirective)

    return {"version": "0.1", "parallel_read_safe": True, "parallel_write_safe": True}

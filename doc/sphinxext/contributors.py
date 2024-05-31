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
    """
    A custom Sphinx directive to generate a contributors section for a specified revision range.

    This directive processes a given revision range to build and display a list of contributors 
    in a Sphinx document. If the revision range ends with 'x..HEAD', an empty paragraph and 
    bullet list are returned. If an error occurs while fetching the contributors, a warning is 
    issued.

    Parameters
    ----------
    required_arguments : int
        Number of required arguments for the directive, set to 1.
    name : str
        The name of the directive, set to "contributors".

    Returns
    -------
    list
        A list of Sphinx nodes containing the contributors section.

    Examples
    --------
    .. contributors:: v1.0.0..v2.0.0

    This directive can be used in Sphinx documentation to automatically generate and include a 
    contributors section based on a specified revision range.
    """
    required_arguments = 1
    name = "contributors"

    def run(self):
        """
        Process the directive and generate the contributors section.

        This method retrieves the specified revision range from the directive's arguments. 
        If the range ends with 'x..HEAD', it returns an empty paragraph and bullet list. 
        Otherwise, it attempts to build the components for the contributors section using 
        the `build_components` function. If an error occurs during this process, it issues 
        a warning. On success, it creates a paragraph node with the author message and a 
        bullet list node with the list of authors.

        Returns
        -------
        list
            A list of Sphinx nodes representing the contributors section. This includes 
            a paragraph node for the author message and a bullet list node with the authors.

        Raises
        ------
        git.GitCommandError
            If an error occurs while retrieving the contributors, a warning node is returned 
            with the error message.

        Examples
        --------
        This method is called automatically by the Sphinx framework when processing the 
        `contributors` directive. For example, in a Sphinx document:

        .. contributors:: v1.0.0..v2.0.
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
    Setup function to register the ContributorsDirective with a Sphinx application.

    This function adds the 'contributors' directive to the Sphinx application, enabling the 
    use of the custom directive to generate contributors sections in the documentation.

    Parameters
    ----------
    app : Sphinx
        The Sphinx application object.

    Returns
    -------
    dict
        A dictionary with metadata about the extension, indicating version and parallel 
        read/write safety.

    Examples
    --------
    In your Sphinx `conf.py` file, include the following line to register the extension:
    
    >>> def setup(app):
    >>>     app.add_directive("contributors", ContributorsDirective)

    """
    app.add_directive("contributors", ContributorsDirective)

    return {"version": "0.1", "parallel_read_safe": True, "parallel_write_safe": True}
#!/usr/bin/env python3
"""
Script to generate contributor and pull request lists

This script generates contributor and pull request lists for release
announcements using GitHub v3 protocol. Use requires an authentication token in
order to have sufficient bandwidth, you can get one following the directions at
`<https://help.github.com/articles/creating-an-access-token-for-command-line-use/>_
Don't add any scope, as the default is read access to public information. The
token may be stored in an environment variable as you only get one chance to
see it.

Usage::

    $ ./scripts/announce.py <token> <revision range>

The output is utf8 rst.

Dependencies
------------

- gitpython
- pygithub

Some code was copied from scipy `tools/gh_lists.py` and `tools/authors.py`.

Examples
--------

From the bash command line with $GITHUB token.

    $ ./scripts/announce.py $GITHUB v1.11.0..v1.11.1 > announce.rst

"""

import codecs
import os
import re
import textwrap

from git import Repo

# Contributors to be renamed.
CONTRIBUTOR_MAPPING = {"znkjnffrezna": "znetbgcubravk"}

UTF8Writer = codecs.getwriter("utf8")
this_repo = Repo(os.path.join(os.path.dirname(__file__), "..", ".."))

author_msg = """\
A total of %d people contributed patches to this release.  People with a
"+" by their names contributed a patch for the first time.
"""

pull_request_msg = """\
A total of %d pull requests were merged for this release.
"""


def get_authors(revision_range):
    """
    Get authors from git log.

    This function retrieves the authors from the git log within the specified
    revision range. It compares the authors from the current release to those
    from the previous release to determine the authors for the current release.

    Parameters
    ----------
    revision_range : str
        The revision range to get the authors from, specified in the format
        'previous_revision..current_revision'.

    Returns
    -------
    authors : list
        A list of authors sorted alphabetically, with new authors (not present
        in the previous release) marked with a '+'.

    See Also
    --------
    CONTRIBUTOR_MAPPING : dict
        A mapping of contributors' names for renaming purposes.

    Notes
    -----
    - The function assumes that the `this_repo` object is an instance of a Git
      repository and that it has methods like `git.log` and `git.shortlog` to
      fetch commit logs and short logs, respectively.
    - The `revision_range` parameter should be a valid range in the format
      'previous_revision..current_revision'. If the current revision includes
      `HEAD`, it should be specified as 'previous_revision..current_revision|HEAD'.
    - The function handles `Co-authored-by` commits separately to account for
      backported changes, which are typically merged by bots.
    - The function discards contributions by automated merge bots (e.g., 'Homu')
      to provide a cleaner list of human contributors.
    - Contributor names are updated according to the `CONTRIBUTOR_MAPPING` to
      ensure consistent naming.

    Examples
    --------
    >>> get_authors('v1.0.0..v1.0.1')
    ['Joris Van den Bossche', 'Tom Augspurger', 'Jeff Reback', 'Philip Cloud', 'Stephan Hoyer']
    >>> get_authors('v1.0.1..HEAD')
    ['Joris Van den Bossche', 'Tom Augspurger', 'Jeff Reback', 'Philip Cloud', 'Stephan Hoyer', 'Simon Hawkins +']
    """
    pat = "^.*\\t(.*)$"
    lst_release, cur_release = (r.strip() for r in revision_range.split(".."))

    if "|" in cur_release:
        # e.g. v1.0.1|HEAD
        maybe_tag, head = cur_release.split("|")
        assert head == "HEAD"
        if maybe_tag in this_repo.tags:
            cur_release = maybe_tag
        else:
            cur_release = head
        revision_range = f"{lst_release}..{cur_release}"

    # authors, in current release and previous to current release.
    # We need two passes over the log for cur and prev, one to get the
    # "Co-authored by" commits, which come from backports by the bot,
    # and one for regular commits.
    xpr = re.compile(r"Co-authored-by: (?P<name>[^<]+) ")
    cur = set(
        xpr.findall(
            this_repo.git.log("--grep=Co-authored", "--pretty=%b", revision_range)
        )
    )
    cur |= set(re.findall(pat, this_repo.git.shortlog("-s", revision_range), re.M))

    pre = set(
        xpr.findall(this_repo.git.log("--grep=Co-authored", "--pretty=%b", lst_release))
    )
    pre |= set(re.findall(pat, this_repo.git.shortlog("-s", lst_release), re.M))

    # Homu is the author of auto merges, clean him out.
    cur.discard("Homu")
    pre.discard("Homu")

    # Rename contributors according to mapping.
    for old_name, new_name in CONTRIBUTOR_MAPPING.items():
        old_name_decoded = codecs.decode(old_name, "rot13")
        new_name_decoded = codecs.decode(new_name, "rot13")
        if old_name_decoded in pre:
            pre.discard(old_name_decoded)
            pre.add(new_name_decoded)
        if old_name_decoded in cur:
            cur.discard(old_name_decoded)
            cur.add(new_name_decoded)

    # Append '+' to new authors.
    authors = [s + " +" for s in cur - pre] + list(cur & pre)
    authors.sort()
    return authors


def get_pull_requests(repo, revision_range):
    """
    Retrieve a list of pull requests merged in a given revision range of a repository.

    This function identifies pull requests from regular merges, Homu auto merges,
    and fast forward squash-merges within the specified revision range, and returns
    a list of pull request objects from the GitHub repository.

    Parameters
    ----------
    repo : object
        A GitHub repository object that provides the method `get_pull` to retrieve pull request data.
    revision_range : str
        The revision range to search for merged pull requests, specified in the format
        'start_revision..end_revision'.

    Returns
    -------
    list
        A list of pull request objects corresponding to the pull request numbers found in the specified
        revision range.

    Notes
    -----
    - The function assumes that the `this_repo` object is an instance of a Git repository
      and that it has a method `git.log` to fetch commit logs.
    - Pull request numbers are identified from commit messages that follow specific patterns:
      - "Merge pull request #<number>"
      - "Auto merge of #<number>"
      - "fast forward squash-merge (#<number>)"
    - The function sorts the pull request numbers before retrieving the pull request objects.

    Examples
    --------
    >>> get_pull_requests(repo, 'v1.0.0..v1.1.0')
    [<PullRequest #123>, <PullRequest #124>, <PullRequest #125>]
    """
    prnums = []

    # From regular merges
    merges = this_repo.git.log("--oneline", "--merges", revision_range)
    issues = re.findall("Merge pull request \\#(\\d*)", merges)
    prnums.extend(int(s) for s in issues)

    # From Homu merges (Auto merges)
    issues = re.findall("Auto merge of \\#(\\d*)", merges)
    prnums.extend(int(s) for s in issues)

    # From fast forward squash-merges
    commits = this_repo.git.log(
        "--oneline", "--no-merges", "--first-parent", revision_range
    )
    issues = re.findall("^.*\\(\\#(\\d+)\\)$", commits, re.M)
    prnums.extend(int(s) for s in issues)

    # get PR data from GitHub repo
    prnums.sort()
    prs = [repo.get_pull(n) for n in prnums]
    return prs


def build_components(revision_range, heading="Contributors"):
    """
    Build components for the contributors section based on a revision range.

    This function generates the components needed for constructing a formatted
    contributors section. It extracts the list of authors who contributed within
    the specified revision range and prepares a dictionary with the heading, a
    message about the authors, and the list of authors.

    Parameters
    ----------
    revision_range : str
        The revision range to get the contributors from, specified in the format
        'start_revision..end_revision'.
    heading : str, optional
        The heading for the contributors section. Default is "Contributors".

    Returns
    -------
    dict
        A dictionary with the following keys:
        - "heading" : str
            The heading for the contributors section.
        - "author_message" : str
            A message about the number of authors.
        - "authors" : list
            A list of authors who contributed within the specified revision range.

    Notes
    -----
    - The function assumes that the `get_authors` function is available and returns a list of authors
      for the given revision range.
    - The `author_msg` string is assumed to be a template string that takes a single integer, representing
      the number of authors.

    Examples
    --------
    >>> build_components('v1.0.0..v1.1.0')
    {'heading': 'Contributors', 'author_message': 'There are 5 contributors:', 'authors': ['Author One', 'Author Two']}
    """
    lst_release, cur_release = (r.strip() for r in revision_range.split(".."))
    authors = get_authors(revision_range)

    return {
        "heading": heading,
        "author_message": author_msg % len(authors),
        "authors": authors,
    }


def build_string(revision_range, heading="Contributors"):
    """
    Build a formatted string of contributors for a given revision range.

    This function constructs a formatted string listing the contributors within the specified
    revision range, including a heading and a list of authors. The heading is underlined for
    emphasis, and the authors are listed with bullet points.

    Parameters
    ----------
    revision_range : str
        The revision range to get the contributors from, specified in the format
        'start_revision..end_revision'.
    heading : str, optional
        The heading for the contributors section. Default is "Contributors".

    Returns
    -------
    str
        A formatted string with the heading, underlined heading, a message, and a list of authors.

    Notes
    -----
    - The function assumes that the `build_components` function is available and returns a dictionary
      with keys "heading", "author_message", and "authors".
    - The length of the underline is dynamically set to match the length of the heading.
    - Authors are listed with bullet points, each preceded by an asterisk (*) and a newline character.
    - The `textwrap.dedent` method is used to maintain proper formatting of the template string.

    Examples
    --------
    >>> build_string('v1.0.0..v1.1.0')
    'Contributors\n============\n\nList of authors who contributed:\n* Author One\n* Author Two'
    """
    components = build_components(revision_range, heading=heading)
    components["uline"] = "=" * len(components["heading"])
    components["authors"] = "* " + "\n* ".join(components["authors"])

    # Don't change this to an fstring. It breaks the formatting.
    tpl = textwrap.dedent(
        """\
    {heading}
    {uline}

    {author_message}
    {authors}"""
    ).format(**components)
    return tpl


def main(revision_range):
    """
    Generate and print the contributors list for a given revision range.

    This function builds a formatted string of contributors who have made changes
    within the specified revision range and prints it to the standard output.

    Parameters
    ----------
    revision_range : str
        The revision range to get the contributors from, specified in the format
        'start_revision..end_revision'.

    Returns
    -------
    None

    Notes
    -----
    - The function relies on `build_string` to generate the formatted contributors list.
    - The contributors list includes a heading, an underlined heading, a message about the authors,
      and a list of authors with new contributors marked.

    Examples
    --------
    >>> main('v1.0.0..v1.1.0')

    Authors
    ============
    - TomAugspurger
    - gfyoung
    - datapythonista
    - jreback
    - jschendel
    - ...
    """
    text = build_string(revision_range)
    print(text)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Generate author lists for release")
    parser.add_argument("revision_range", help="<revision>..<revision>")
    args = parser.parse_args()
    main(args.revision_range)

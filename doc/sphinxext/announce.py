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
    Extracts a list of code authors within a specified revision range.

    This function analyzes the git log for the given revision range and identifies
    all the contributors who authored or co-authored commits. It considers both
    regular commits and "Co-authored-by" commits used for backports.

    Parameters
    ----------
    revision_range : str 
        A string representing the revision range for analysis. The Expected format is: "<start-revision>..<end-revision>".

    Returns
    -------
    list
        A alphabetically sorted list of author names (strings) who contributed within the revision range. The new authors are marked with a '+'.

    Examples
    --------    
    >>> get_authors("v1.0.0..v1.0.1")
    ['Author1', 'Author2', 'Author3', 'Author4']

    >>> get_authors('v1.0.0..HEAD')
    ['Author1', 'Author2', 'Author3', 'Author4', 'Author5 +'] 
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
    Retrieve pull requests from a Git repository within a specified revision range.

    This function extracts pull request numbers from various types of merge commits (regular merges, Homu auto merges, and fast-forward squash-merges) in a specified revision range of a repository. It then retrieves detailed pull request data from the GitHub repository based on these numbers.

    Parameters
    ----------
    repo : Repository object
        The GitHub repository object from which to retrieve pull request data.
    revision_range : str
        The range of revisions to search for pull request merges, specified in a 
        format recognized by `git log`.

    Returns
    -------
    list 
        A sorted list of pull request objects corresponding to the pull request numbers found in the specified revision range.

    See Also
    --------
    Repository.get_pull : Retrieve a pull request by number from the repository.

    Examples
    --------
    >>> prs = get_pull_requests(repo, "v1.0.0...v2.0.0")
    >>> prs
    [<Pull Request 24>, <Pull Request 25>, <Pull Request 26>]
    >>> for pr in prs:
    ...     print(pr.number, pr.title)
    1 Fix bug in feature X
    2 Improve documentation
    3 Add new feature Y
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
    Build the components for a contributors section based on a revision range.

    This function extracts the list of authors who contributed within a specified revision range and constructs the components needed for a contributors section, including a heading and an author message.

    Parameters
    ----------
    revision_range : str
        The range of revisions to search for authors, specified in a format recognized by `git log`.
    heading : str, optional
        The heading for the contributors section, default is "Contributors".

    Returns
    -------
    dict
        A dictionary containing the heading, author message, and list of authors who contributed within the specified revision range.

    See Also
    --------
    get_authors : Retrieve a list of authors who contributed within a specified revision range.

    Examples
    --------
    >>> components = build_components("v1.0.0...v2.0.0")
    >>> components["heading"]
    Contributors
    >>> components["author_message"]
    There are 10 contributors.
    >>> components["authors"]
    ['Author1', 'Author2', 'Author3', ...]
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
    Build a formatted string for the contributors section based on a revision range.

    This function creates a formatted string that includes a heading, an author message, and a list of authors who contributed within a specified revision range. The formatting is designed to be suitable for inclusion in documentation.

    Parameters
    ----------
    revision_range : str
        The range of revisions to search for authors, specified in a format recognized by `git log`.
    heading : str, optional
        The heading for the contributors section, default is "Contributors".

    Returns
    -------
    str
        A formatted string containing the contributors section.

    See Also
    --------
    build_components : Build the components for a contributors section.

    Examples
    --------
    >>> contrib_string = build_string("v1.0.0...v2.0.0")
    >>> contrib_string
    Contributors
    ============
    
    There are 10 contributors.
    * Author1
    * Author2
    * Author3
    ...
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
    Main function to document authors based on a revision range.

    This function generates and prints a formatted contributors section for a specified 
    revision range. It serves as the entry point for generating the contributors section 
    output.

    Parameters
    ----------
    revision_range : str
        The range of revisions to search for authors, specified in a format recognized 
        by `git log`.

    Examples
    --------
    >>> main("v1.0.0...v2.0.0")
    Contributors
    ============
    
    There are 10 contributors.
    * Author1
    * Author2
    * Author3
    ...
    """
    # document authors
    text = build_string(revision_range)
    print(text)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Generate author lists for release")
    parser.add_argument("revision_range", help="<revision>..<revision>")
    args = parser.parse_args()
    main(args.revision_range)
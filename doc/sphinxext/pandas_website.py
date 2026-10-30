"""
Sphinx extension for the pages of the pandas website.

The website (https://pandas.pydata.org/) is part of the documentation. Its
pages are the Markdown files in ``doc/source`` (parsed with MyST), which are
rendered with Jinja2 before being parsed, like the rest of the pages (see
``rstjinja`` in ``doc/source/conf.py``). This extension adds to the context
used to render them the data obtained by parsing ``doc/website.yml``, and
enriching it with the context preprocessors listed in its
``main.context_preprocessors`` section (for example, to fetch the latest
releases from GitHub).

Fetching the context requires network access, so it's skipped when the
``pandas_website_context`` config value is false (e.g. when building a single
page); the pages using it are excluded from those builds.
"""

from __future__ import annotations

import collections
import datetime
import importlib
import itertools
import json
import operator
import os
import pathlib
import re
import shutil
import sys
import time
from typing import (
    TYPE_CHECKING,
    Any,
)
import unicodedata

from docutils import nodes
import feedparser
from markdown_it import MarkdownIt
from packaging import version
import requests
import yaml

if TYPE_CHECKING:
    from sphinx.application import Sphinx

CONFIG_PATH = pathlib.Path(__file__).parents[1] / "website.yml"

api_token = os.environ.get("GITHUB_TOKEN")
if api_token is not None:
    GITHUB_API_HEADERS = {"Authorization": f"Bearer {api_token}"}
else:
    GITHUB_API_HEADERS = {}


class Preprocessors:
    """
    Built-in context preprocessors.

    Context preprocessors are functions that receive the context used to
    render the templates, and enriches it with additional information.

    The original context is obtained by parsing ``config.yml``, and
    anything else needed just be added with context preprocessors.
    """

    @staticmethod
    def blog_add_posts(context):
        """
        Given the blog feed defined in the configuration yaml, this context
        preprocessor fetches the posts in the feeds, and returns the relevant
        information for them (sorted from newest to oldest).
        """
        tag_expr = re.compile("<.*?>")
        posts = []
        # posts from the file system
        if context["blog"]["posts_path"]:
            posts_path = context["source_path"] / context["blog"]["posts_path"]
            for fname in posts_path.iterdir():
                if fname.name.startswith("index."):
                    continue
                link = f"{context['blog']['posts_path']}/{fname.stem}.html"
                with fname.open(encoding="utf-8") as f:
                    _, front_matter, body = f.read().split("---\n", 2)
                meta = yaml.safe_load(front_matter)
                html = MarkdownIt("commonmark").render(body)
                title = meta["title"]
                summary = re.sub(tag_expr, "", html)
                try:
                    body_position = summary.index(title) + len(title)
                except ValueError as err:
                    raise ValueError(
                        f'Blog post "{fname}" should have a markdown header '
                        f'corresponding to its "title" front matter "{title}"'
                    ) from err
                summary = " ".join(summary[body_position:].split(" ")[:30])
                posts.append(
                    {
                        "title": title,
                        "author": context["blog"]["author"],
                        "published": datetime.datetime.combine(
                            meta["date"], datetime.time()
                        ),
                        "feed": context["blog"]["feed_name"],
                        "link": link,
                        "description": summary,
                        "summary": summary,
                    }
                )
        # posts from rss feeds
        for feed_url in context["blog"]["feed"]:
            feed_data = feedparser.parse(feed_url)
            for entry in feed_data.entries:
                published = datetime.datetime.fromtimestamp(
                    time.mktime(entry.published_parsed)
                )
                summary = re.sub(tag_expr, "", entry.summary)
                posts.append(
                    {
                        "title": entry.title,
                        "author": entry.author,
                        "published": published,
                        "feed": feed_data["feed"]["title"],
                        "link": entry.link,
                        "description": entry.description,
                        "summary": summary,
                    }
                )
        posts.sort(key=operator.itemgetter("published"), reverse=True)
        context["blog"]["posts"] = posts[: context["blog"]["num_posts"]]
        return context

    @staticmethod
    def maintainers_add_info(context):
        """
        Given the active maintainers defined in the yaml file, it fetches
        the GitHub user information for them.
        """
        repeated = set(context["maintainers"]["active"]) & set(
            context["maintainers"]["inactive"]
        )
        if repeated:
            raise ValueError(f"Maintainers {repeated} are both active and inactive")

        maintainers_info = {}
        for user in (
            context["maintainers"]["active"]
            + context["maintainers"]["inactive"]
            + context["maintainers"]["pandasstubs"]
        ):
            resp = requests.get(
                f"https://api.github.com/users/{user}",
                headers=GITHUB_API_HEADERS,
                timeout=5,
            )
            if resp.status_code == 403:
                sys.stderr.write(
                    "WARN: GitHub API quota exceeded when fetching maintainers\n"
                )
                # if we exceed github api quota, we use the github info
                # of maintainers saved with the website
                resp_bkp = requests.get(
                    context["main"]["production_url"] + "maintainers.json", timeout=5
                )
                resp_bkp.raise_for_status()
                maintainers_info = resp_bkp.json()
                break

            resp.raise_for_status()
            # only keep the fields used in the templates; the full response also
            # includes other fields that change independently (e.g. follower
            # count), causing a diff on every rebuild
            user_json = resp.json()
            maintainers_info[user] = {
                key: user_json[key]
                for key in ("name", "login", "avatar_url", "html_url", "blog")
            }

        context["maintainers"]["github_info"] = maintainers_info

        # save the data fetched from github to use it in case we exceed
        # git github api quota in the future
        with open(
            pathlib.Path(context["target_path"]) / "maintainers.json",
            "w",
            encoding="utf-8",
        ) as f:
            json.dump(maintainers_info, f)

        return context

    @staticmethod
    def home_add_releases(context):
        context["releases"] = []

        github_repo_url = context["main"]["github_repo_url"]
        resp = requests.get(
            f"https://api.github.com/repos/{github_repo_url}/releases",
            headers=GITHUB_API_HEADERS,
            timeout=5,
        )
        if resp.status_code == 403:
            sys.stderr.write("WARN: GitHub API quota exceeded when fetching releases\n")
            resp_bkp = requests.get(
                context["main"]["production_url"] + "releases.json", timeout=5
            )
            resp_bkp.raise_for_status()
            releases = resp_bkp.json()
        else:
            resp.raise_for_status()
            releases = resp.json()

            # only keep the fields actually used when rendering the site; the full
            # response also includes fields like download_count, causing the file
            # to change on every rebuild
            releases = [
                {
                    "tag_name": release["tag_name"],
                    "published_at": release["published_at"],
                    "prerelease": release["prerelease"],
                    "browser_download_url": (
                        release["assets"][0]["browser_download_url"]
                        if release["assets"]
                        else ""
                    ),
                }
                for release in releases
            ]

        with open(
            pathlib.Path(context["target_path"]) / "releases.json",
            "w",
            encoding="utf-8",
        ) as f:
            json.dump(releases, f, default=datetime.datetime.isoformat)

        for release in releases:
            if release["prerelease"]:
                continue
            name = release["tag_name"].lstrip("v")
            parsed_version = version.parse(name)
            published = datetime.datetime.strptime(
                release["published_at"], "%Y-%m-%dT%H:%M:%SZ"
            )
            context["releases"].append(
                {
                    "name": name,
                    "parsed_version": parsed_version,
                    "short_name": f"{parsed_version.major}.{parsed_version.minor}",
                    "tag": release["tag_name"],
                    "published": published,
                    "url": release["browser_download_url"],
                }
            )
        # sorting out obsolete versions
        grouped_releases = itertools.groupby(
            context["releases"],
            key=lambda r: (r["parsed_version"].major, r["parsed_version"].minor),
        )
        context["releases"] = [
            max(release_group, key=lambda r: r["parsed_version"].minor)
            for _, release_group in grouped_releases
        ]
        # sorting releases by version number
        context["releases"].sort(key=lambda r: r["parsed_version"], reverse=True)
        return context

    @staticmethod
    def roadmap_pdeps(context):
        """
        PDEP's (pandas enhancement proposals) are not part of the bar
        navigation. They are included as lists in the "Roadmap" page
        and linked from there. This preprocessor obtains the list of
        PDEP's in different status from the directory tree and GitHub.
        """
        KNOWN_STATUS = {
            "Draft",
            "Under discussion",
            "Accepted",
            "Implemented",
            "Rejected",
            "Withdrawn",
        }
        context["pdeps"] = collections.defaultdict(list)

        # accepted, rejected and implemented
        pdeps_path = (
            pathlib.Path(context["source_path"]) / context["roadmap"]["pdeps_path"]
        )
        for pdep in sorted(pdeps_path.iterdir()):
            if pdep.suffix != ".md":
                continue
            with pdep.open() as f:
                title = f.readline()[2:]  # removing markdown title "# "
                status = None
                for line in f:
                    if line.startswith("- Status: "):
                        status = line.strip().split(": ", 1)[1]
                        break
                if status not in KNOWN_STATUS:
                    raise RuntimeError(
                        f'PDEP "{pdep}" status "{status}" is unknown. '
                        f"Should be one of: {KNOWN_STATUS}"
                    )
            html_file = pdep.with_suffix(".html").name
            context["pdeps"][status].append(
                {
                    "title": title,
                    "url": f"pdeps/{html_file}",
                }
            )

        # under discussion
        github_repo_url = context["main"]["github_repo_url"]
        resp = requests.get(
            "https://api.github.com/search/issues?"
            f"q=is:pr is:open label:PDEP draft:false repo:{github_repo_url}",
            headers=GITHUB_API_HEADERS,
            timeout=5,
        )
        if resp.status_code == 403:
            sys.stderr.write("WARN: GitHub API quota exceeded when fetching pdeps\n")
            resp_bkp = requests.get(
                context["main"]["production_url"] + "pdeps.json", timeout=5
            )
            resp_bkp.raise_for_status()
            pdeps = resp_bkp.json()
        else:
            resp.raise_for_status()
            pdeps = resp.json()

        with open(
            pathlib.Path(context["target_path"]) / "pdeps.json", "w", encoding="utf-8"
        ) as f:
            json.dump(pdeps, f)

        compiled_pattern = re.compile(r"^PDEP-(\d+)")

        def sort_pdep(pdep: dict) -> int:
            title = pdep["title"]
            match = compiled_pattern.match(title)
            if not match:
                msg = f"""Could not find PDEP number in '{title}'. Please make sure to
                write the title as: 'PDEP-num: {title}'."""
                raise ValueError(msg)

            return int(match[1])

        context["pdeps"]["Under discussion"].extend(
            {"title": pdep["title"], "url": pdep["html_url"]}
            for pdep in sorted(pdeps["items"], key=sort_pdep)
        )

        return context


def get_callable(obj_as_str: str) -> object:
    """
    Get a Python object from its string representation.

    For example, for ``sys.stdout.write`` would import the module ``sys``
    and return the ``write`` function.
    """
    components = obj_as_str.split(".")
    attrs = []
    while components:
        try:
            obj = importlib.import_module(".".join(components))
        except ImportError:
            attrs.insert(0, components.pop())
        else:
            break

    if not obj:
        raise ImportError(f'Could not import "{obj_as_str}"')

    for attr in attrs:
        obj = getattr(obj, attr)

    return obj


def get_context(config_fname: pathlib.Path, **kwargs):
    """
    Load the config yaml as the base context, and enrich it with the
    information added by the context preprocessors defined in the file.
    """
    with config_fname.open(encoding="utf-8") as f:
        context = yaml.safe_load(f)

    context["source_path"] = config_fname.parent
    context.update(kwargs)

    preprocessors = (
        get_callable(context_prep)
        for context_prep in context["main"]["context_preprocessors"]
    )
    for preprocessor in preprocessors:
        context = preprocessor(context)
        msg = f"{preprocessor.__name__} is missing the return statement"
        assert context is not None, msg

    return context


def is_website_page(app: Sphinx, docname: str) -> bool:
    """
    Whether ``docname`` is a page of the website (i.e. a Markdown file).
    """
    return app.env.doc2path(docname).suffix == ".md"


def load_context(app: Sphinx) -> None:
    """
    Build the context used to render the website pages, and make it available
    to all pages and templates (as ``html_context``). The preprocessors save
    some of the data as JSON files in the output, to be used if the GitHub API
    quota is exceeded in a later build.
    """
    if not app.config.pandas_website_context:
        return
    outdir = pathlib.Path(app.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    context = get_context(
        CONFIG_PATH, source_path=pathlib.Path(app.srcdir), target_path=outdir
    )
    app.config.html_context.update(context)


def _legacy_heading_id(title: str, used_ids: set[str]) -> str:
    """
    Return the id that the website generator previously used for a heading
    (Python-Markdown's ``toc`` extension), given the ids already used.
    """
    value = unicodedata.normalize("NFKD", title).encode("ascii", "ignore").decode()
    value = re.sub(r"[^\w\s-]", "", value).strip().lower()
    heading_id = re.sub(r"[-\s]+", "-", value)
    while heading_id in used_ids or not heading_id:
        match = re.match(r"^(.*)_([0-9]+)$", heading_id)
        if match:
            heading_id = f"{match.group(1)}_{int(match.group(2)) + 1}"
        else:
            heading_id = f"{heading_id}_1"
    used_ids.add(heading_id)
    return heading_id


def add_legacy_heading_ids(app: Sphinx, doctree: nodes.document) -> None:
    """
    Add to the headings of the website pages the ids they had before the
    website was part of the documentation, so links to their sections work.
    """
    if not is_website_page(app, app.env.docname):
        return
    used_ids: set[str] = set()
    for section in doctree.findall(nodes.section):
        heading_id = _legacy_heading_id(section[0].astext(), used_ids)
        if heading_id not in doctree.ids:
            section["ids"].append(heading_id)
            doctree.ids[heading_id] = section


def copy_static(app: Sphinx, exception: Exception | None) -> None:
    """
    Copy ``static/`` to the output keeping its path, since the website pages
    link to its files relative to the root of the site.
    """
    if exception is None and app.builder.format == "html":
        shutil.copytree(
            pathlib.Path(app.srcdir) / "static",
            pathlib.Path(app.outdir) / "static",
            dirs_exist_ok=True,
        )


def setup(app: Sphinx) -> dict[str, Any]:
    app.add_config_value("pandas_website_context", True, "env", types=[bool])
    app.connect("builder-inited", load_context)
    app.connect("doctree-read", add_legacy_heading_ids)
    app.connect("build-finished", copy_static)
    return {"parallel_read_safe": True, "parallel_write_safe": True}

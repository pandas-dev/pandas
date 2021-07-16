#!/usr/bin/env python3
"""
Simple static site generator for the pandas web.

pandas_web.py takes a directory as parameter, and copies all the files into the
target directory after converting markdown files into html and rendering both
markdown and html files with a context. The context is obtained by parsing
the file ``config.yml`` in the root of the source directory.

The file should contain:
```
main:
  template_path: <path_to_the_jinja2_templates_directory>
  base_template: <template_file_all_other_files_will_extend>
  ignore:
  - <list_of_files_in_the_source_that_will_not_be_copied>
  github_repo_url: <organization/repo-name>
  context_preprocessors:
  - <list_of_functions_that_will_enrich_the_context_parsed_in_this_file>
  markdown_extensions:
  - <list_of_markdown_extensions_that_will_be_loaded>
```

The rest of the items in the file will be added directly to the context.
"""
import argparse
import datetime
import importlib
import operator
import os
import re
import shutil
import sys
import time
import typing

import feedparser
import jinja2
import markdown
import requests
import yaml


class Preprocessors:
    """
    Built-in context preprocessors.

    Context preprocessors are functions that receive the context used to
    render the templates, and enriches it with additional information.

    The original context is obtained by parsing ``config.yml``, and
    anything else needed just be added with context preprocessors.
    """

    @staticmethod
    def navbar_add_info(context):
        """
        Items in the main navigation bar can be direct links, or dropdowns with
        subitems. This context preprocessor adds a boolean field
        ``has_subitems`` that tells which one of them every element is. It
        also adds a ``slug`` field to be used as a CSS id.
        """
        for i, item in enumerate(context["navbar"]):
            context["navbar"][i] = dict(
                item,
                has_subitems=isinstance(item["target"], list),
                slug=(item["name"].replace(" ", "-").lower()),
            )
        return context

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
            posts_path = os.path.join(
                context["source_path"], *context["blog"]["posts_path"].split("/")
            )
            for fname in os.listdir(posts_path):
                if fname.startswith("index."):
                    continue
                link = (
                    f"/{context['blog']['posts_path']}"
                    f"/{os.path.splitext(fname)[0]}.html"
                )
                md = markdown.Markdown(
                    extensions=context["main"]["markdown_extensions"]
                )
                with open(os.path.join(posts_path, fname)) as f:
                    html = md.convert(f.read())
                title = md.Meta["title"][0]
                summary = re.sub(tag_expr, "", html)
                try:
                    body_position = summary.index(title) + len(title)
                except ValueError:
                    raise ValueError(
                        f'Blog post "{fname}" should have a markdown header '
                        f'corresponding to its "Title" element "{title}"'
                    )
                summary = " ".join(summary[body_position:].split(" ")[:30])
                posts.append(
                    {
                        "title": title,
                        "author": context["blog"]["author"],
                        "published": datetime.datetime.strptime(
                            md.Meta["date"][0], "%Y-%m-%d"
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
        context["maintainers"]["people"] = []
        for user in context["maintainers"]["active"]:
            resp = requests.get(f"https://api.github.com/users/{user}")
            if context["ignore_io_errors"] and resp.status_code == 403:
                return context
            resp.raise_for_status()
            context["maintainers"]["people"].append(resp.json())
        return context

    @staticmethod
    def home_add_releases(context):
        context["releases"] = []

        github_repo_url = context["main"]["github_repo_url"]
        resp = requests.get(f"https://api.github.com/repos/{github_repo_url}/releases")
        if context["ignore_io_errors"] and resp.status_code == 403:
            return context
        resp.raise_for_status()

        for release in resp.json():
            if release["prerelease"]:
                continue
            published = datetime.datetime.strptime(
                release["published_at"], "%Y-%m-%dT%H:%M:%SZ"
            )
            context["releases"].append(
                {
                    "name": release["tag_name"].lstrip("v"),
                    "tag": release["tag_name"],
                    "published": published,
                    "url": (
                        release["assets"][0]["browser_download_url"]
                        if release["assets"]
                        else ""
                    ),
                }
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


def get_context(config_fname: str, ignore_io_errors: bool, **kwargs):
    """
    Load the config yaml as the base context, and enrich it with the
    information added by the context preprocessors defined in the file.
    """
    with open(config_fname) as f:
        context = yaml.safe_load(f)

    context["source_path"] = os.path.dirname(config_fname)
    context["ignore_io_errors"] = ignore_io_errors
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


def get_source_files(source_path: str) -> typing.Generator[str, None, None]:
    """
    Generate the list of files present in the source directory.
    """
    for root, dirs, fnames in os.walk(source_path):
        root = os.path.relpath(root, source_path)
        for fname in fnames:
            yield os.path.join(root, fname)


def extend_base_template(content: str, base_template: str) -> str:
    """
    Wrap document to extend the base template, before it is rendered with
    Jinja2.
    """
    result = '{% extends "' + base_template + '" %}'
    result += "{% block body %}"
    result += content
    result += "{% endblock %}"
    return result


def main(
    source_path: str, target_path: str, base_url: str, ignore_io_errors: bool
) -> int:
    """
    Copy every file in the source directory to the target directory.

    For ``.md`` and ``.html`` files, render them with the context
    before copyings them. ``.md`` files are transformed to HTML.
    """
    config_fname = os.path.join(source_path, "config.yml")

    shutil.rmtree(target_path, ignore_errors=True)
    os.makedirs(target_path, exist_ok=True)

    sys.stderr.write("Generating context...\n")
    context = get_context(config_fname, ignore_io_errors, base_url=base_url)
    sys.stderr.write("Context generated\n")

    templates_path = os.path.join(source_path, context["main"]["templates_path"])
    jinja_env = jinja2.Environment(loader=jinja2.FileSystemLoader(templates_path))

    for fname in get_source_files(source_path):
        if os.path.normpath(fname) in context["main"]["ignore"]:
            continue

        sys.stderr.write(f"Processing {fname}\n")
        dirname = os.path.dirname(fname)
        os.makedirs(os.path.join(target_path, dirname), exist_ok=True)

        extension = os.path.splitext(fname)[-1]
        if extension in (".html", ".md"):
            with open(os.path.join(source_path, fname)) as f:
                content = f.read()
            if extension == ".md":
                body = markdown.markdown(
                    content, extensions=context["main"]["markdown_extensions"]
                )
                content = extend_base_template(body, context["main"]["base_template"])
            content = jinja_env.from_string(content).render(**context)
            fname = os.path.splitext(fname)[0] + ".html"
            with open(os.path.join(target_path, fname), "w") as f:
                f.write(content)
        else:
            shutil.copy(
                os.path.join(source_path, fname), os.path.join(target_path, dirname)
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Documentation builder.")
    parser.add_argument(
        "source_path", help="path to the source directory (must contain config.yml)"
    )
    parser.add_argument(
        "--target-path", default="build", help="directory where to write the output"
    )
    parser.add_argument(
        "--base-url", default="", help="base url where the website is served from"
    )
    parser.add_argument(
        "--ignore-io-errors",
        action="store_true",
        help="do not fail if errors happen when fetching "
        "data from http sources, and those fail "
        "(mostly useful to allow github quota errors "
        "when running the script locally)",
    )
    args = parser.parse_args()
    sys.exit(
        main(args.source_path, args.target_path, args.base_url, args.ignore_io_errors)
    )

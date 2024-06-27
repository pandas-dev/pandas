# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os

from ..commands import Command
from ..commands.publish import Publish
from ..console import log
from .. import util


class GithubPages(Command):
    @classmethod
    def setup_arguments(cls, subparsers):
        parser = subparsers.add_parser(
            "gh-pages", help="Publish the results to Github pages",
            description="""
            Publish the results to github pages

            Updates the 'gh-pages' branch in the current repository,
            and pushes it to 'origin'.
            """)
        parser.add_argument(
            "--no-push", action="store_true",
            help="Update local gh-pages branch but don't push")
        parser.add_argument(
            "--rewrite", action="store_true",
            help=("Rewrite gh-pages branch to contain only a single commit, "
                  "instead of adding a new commit"))

        parser.set_defaults(func=cls.run_from_args)

        return parser

    @classmethod
    def run_from_conf_args(cls, conf, args):
        return cls.run(conf=conf, no_push=args.no_push, rewrite=args.rewrite)

    @classmethod
    def run(cls, conf, no_push, rewrite):
        git = util.which('git')

        # Publish
        Publish.run(conf)

        cwd = os.path.abspath(os.getcwd())

        log.info("Updating gh-pages branch")
        try:
            # Create new repo for the html data
            util.check_call([git, 'init'], cwd=conf.html_dir)
            util.check_call([git, 'checkout', '-b', 'gh-pages'], cwd=conf.html_dir)

            if not rewrite:
                # Fetch the old branch
                _, _, retcode = util.check_output([git, 'fetch', cwd, 'gh-pages'],
                                                  cwd=conf.html_dir,
                                                  return_stderr=True,
                                                  valid_return_codes=None)
                if retcode == 0:
                    util.check_call([git, 'reset', '--soft', 'FETCH_HEAD'], cwd=conf.html_dir)

            # We need to tell github this is not a Jekyll document
            with open(os.path.join(conf.html_dir, '.nojekyll'), 'wb') as fd:
                fd.write(b'\n')

            # Add all files
            util.check_call([git, 'add', '--all', '.'], cwd=conf.html_dir)

            # Commit (if needed)
            _, _, retcode = util.check_output([git, "diff-index", "--quiet", "HEAD", "--"],
                                              cwd=conf.html_dir,
                                              return_stderr=True,
                                              valid_return_codes=None)
            if retcode != 0:
                util.check_call([git, 'commit', '-m', 'Generated from sources'], cwd=conf.html_dir)

            # Fetch branch here
            if rewrite:
                util.check_call([git, 'fetch', os.path.abspath(conf.html_dir)])
                util.check_call([git, 'branch', '-f', 'gh-pages', 'FETCH_HEAD'])
            else:
                util.check_call([git, 'fetch', os.path.abspath(conf.html_dir),
                                 'gh-pages:gh-pages'])
        finally:
            # Cleanup the child repo under html
            util.long_path_rmtree(os.path.join(conf.html_dir, '.git'))

        # Push branch
        if not no_push:
            if rewrite:
                log.info("Force-pushing gh-pages branch")
                util.check_call([git, 'push', '-f', 'origin', 'gh-pages'])
            else:
                log.info("Pushing gh-pages branch")
                util.check_call([git, 'push', 'origin', 'gh-pages'])

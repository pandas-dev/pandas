# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Supports git repositories for the benchmarked project.
"""

import os
import re
import shlex

from ..console import log
from ..repo import Repo, NoSuchNameError
from .. import util


class Git(Repo):
    dvcs = "git"

    def __init__(self, url, mirror_path):
        self._git = util.which("git")
        self._path = os.path.abspath(mirror_path)
        self._pulled = False
        # default branch
        self._default_branch = util.git_default_branch()

        if self.is_local_repo(url):
            # Local repository, no need for mirror
            self._path = os.path.abspath(url)
            self._pulled = True
        elif not self.is_local_repo(self._path):
            if os.path.exists(self._path):
                self._raise_bad_mirror_error(self._path)

            # Clone is missing
            log.info("Cloning project")
            self._run_git(['clone', '--mirror', url, self._path], cwd=None)

    @classmethod
    def is_local_repo(cls, path):
        return os.path.isdir(path) and (
            os.path.exists(os.path.join(path, '.git')) or
            os.path.isdir(os.path.join(path, 'objects')))

    @classmethod
    def url_match(cls, url):
        regexes = [
            r'^https?://.*?\.git$',
            r'^git@.*?\.git$']

        for regex in regexes:
            if re.match(regex, url):
                return True

        # Check for a local path
        if cls.is_local_repo(url):
            return True

        return False

    def _run_git(self, args, cwd=True, **kwargs):
        if cwd is True:
            cwd = self._path
        kwargs['cwd'] = cwd
        env = dict(kwargs.pop('env', os.environ))
        if cwd is not None:
            prev = env.get('GIT_CEILING_DIRECTORIES')
            env['GIT_CEILING_DIRECTORIES'] = os.pathsep.join(
                [os.path.join(os.path.abspath(cwd), os.pardir)] +
                ([prev] if prev is not None else []))
        return util.check_output([self._git] + args, env=env, **kwargs)

    def get_new_range_spec(self, latest_result, branch=None):
        return f'{latest_result}..{self.get_branch_name(branch)}'

    def get_range_spec(self, commit_a, commit_b):
        return f'{commit_a}..{commit_b}'

    def pull(self):
        # We assume the remote isn't updated during the run of asv
        # itself.
        if self._pulled:
            return

        log.info("Fetching recent changes")
        self._run_git(['fetch', 'origin'])
        self._pulled = True

    def checkout(self, path, commit_hash):
        def checkout_existing(display_error):
            # Deinit fails if no submodules, so ignore its failure
            self._run_git(['-c', 'protocol.file.allow=always',
                           'submodule', 'deinit', '-f', '.'],
                          cwd=path, display_error=False, valid_return_codes=None)
            self._run_git(['checkout', '-f', commit_hash],
                          cwd=path, display_error=display_error)
            self._run_git(['clean', '-fdx'],
                          cwd=path, display_error=display_error)
            self._run_git(['-c', 'protocol.file.allow=always',
                           'submodule', 'update', '--init', '--recursive'],
                          cwd=path, display_error=display_error)

        if os.path.isdir(path):
            try:
                checkout_existing(display_error=False)
            except util.ProcessError:
                # Remove and try to re-clone
                util.long_path_rmtree(path)

        if not os.path.isdir(path):
            self._run_git(['clone', '--shared', '--recursive', self._path, path],
                          cwd=None)
            checkout_existing(display_error=True)

    def get_date(self, hash):
        return int(self._run_git(
            ['rev-list', '-n', '1', '--format=%at', hash],
            valid_return_codes=(0, 1), dots=False).strip().split()[-1]) * 1000

    def get_hashes_from_range(self, range_spec):
        args = ['rev-list', '--first-parent']
        if range_spec != "":
            args += shlex.split(range_spec) + ["--"]
        output = self._run_git(args, valid_return_codes=(0, 1), dots=False)
        return output.strip().split()

    def get_hash_from_name(self, name):
        if name is None:
            name = self.get_branch_name()

        # In case of annotated tags, return the hash for the commit
        lookup_name = name + '^{commit}'

        try:
            return self._run_git(['rev-parse', lookup_name],
                                 display_error=False,
                                 dots=False).strip().split()[0]
        except util.ProcessError as err:
            if err.stdout.strip() == lookup_name:
                # Name does not exist
                raise NoSuchNameError(name)
            raise

    def get_hash_from_parent(self, name):
        return self.get_hash_from_name(name + '^')

    def get_name_from_hash(self, commit):
        try:
            name = self._run_git(["name-rev", "--name-only",
                                  "--exclude=remotes/*",
                                  "--no-undefined", commit],
                                 display_error=False).strip()
            if not name:
                return None
        except util.ProcessError:
            # Failed to obtain.
            return None

        # Return tags without prefix
        for prefix in ['tags/']:
            if name.startswith(prefix):
                return name[len(prefix):]

        return name

    def get_tags(self):
        tags = {}
        for tag in self._run_git(["tag", "-l", "--sort=taggerdate"]).splitlines():
            tags[tag] = self._run_git(["rev-list", "-n", "1", tag]).strip()
        return tags

    def get_date_from_name(self, name):
        return self.get_date(name + "^{commit}")

    def get_branch_commits(self, branch):
        return self.get_hashes_from_range(self.get_branch_name(branch))

    def get_revisions(self, commits):
        revisions = {}
        for i, commit in enumerate(self._run_git([
            "rev-list", "--all", "--date-order", "--reverse",
        ]).splitlines()):
            if commit in commits:
                revisions[commit] = i
        return revisions

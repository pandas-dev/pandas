# Licensed under a 3-clause BSD style license - see LICENSE.rst

import datetime

from . import util


class NoSuchNameError(RuntimeError):
    """
    Exception raised if requested branch or commit does not exist.
    """
    pass


class Repo:
    """
    Base class for repository handlers.

    There are no commands that modify the source repository.
    """

    # The default branch name when no branch is configured in asv.conf.json
    _default_branch = None

    def __init__(self, url, mirror_path):
        """
        Create a mirror of the repository at `url`, without a working tree.

        Parameters
        ----------
        url : str
            The URL to the repository to clone

        mirror_path : str
            The local path into which to put a mirror of the repository.
            Creating the local mirror is optional, and may be skipped if
            the performance is adequate also without it.

        """
        raise NotImplementedError()

    def _raise_bad_mirror_error(self, path):
        """
        Internal routine for raising an error message if mirror directory already exists,
        but appears incorrect.
        """
        raise util.UserError(f"Directory '{path}' already exists, but is not a mirror "
                             "of the project repository. Please remove it and try again. "
                             "If the benchmark suite is in the project repository, you can "
                             "also adjust the configuration to use the current "
                             "repository (e.g. \"repo\": \".\") instead of a remote URL "
                             "as the source.")

    def checkout(self, path, commit_hash):
        """
        Check out a clean working tree from the current repository
        to the given path

        Parameters
        ----------
        path : str
            The local path to check out into
        commit_hash : str
            The commit hash to check out

        """
        raise NotImplementedError()

    @classmethod
    def url_match(cls, url):
        """
        Returns `True` if the url is of the right type for this kind
        of repository.
        """
        raise NotImplementedError()

    def get_branch_name(self, branch=None):
        """
        Returns the given branch name or the default branch name if branch is
        None or not specified.
        """
        if branch is None:
            return self._default_branch
        else:
            return branch

    def get_range_spec(self, commit_a, commit_b):
        """
        Returns a formatted string giving the results between
        commit_a (exclusive) and commit_b (inclusive).
        """
        raise NotImplementedError()

    def get_new_range_spec(self, latest_result, branch=None):
        """
        Returns a formatted string giving the results between the
        latest result and the newest hash in a given branch.
        If no branch given, use the 'main' branch.
        """
        raise NotImplementedError()

    def get_date(self, hash):
        """
        Get a JavaScript timestamp for a particular commit.
        """
        raise NotImplementedError()

    def get_hashes_from_range(self, range):
        """
        Get a list of commit hashes given a range specifier.  The
        syntax of the range specifier will depend on the DVCS used.
        """
        raise NotImplementedError()

    def get_hash_from_name(self, name):
        """
        Get a hash from a given tag, branch or hash.  The acceptable
        syntax will depend on the DVCS used.
        """
        raise NotImplementedError()

    def get_hash_from_parent(self, name):
        """
        Checkout the parent of the currently checked out commit.
        """
        raise NotImplementedError()

    def get_name_from_hash(self, commit):
        """
        Get a symbolic name for a commit, if it exists.
        The name is resolvable back to the has via get_hash_form_name.
        If there is no symbolic name, returns None.
        """
        raise NotImplementedError()

    def get_decorated_hash(self, commit, hash_length=8):
        """
        Return commit hash with possible branch/tag information added,
        for displaying to the user.
        """
        name = self.get_name_from_hash(commit)
        if name is not None:
            return f"{commit[:hash_length]} <{name}>"
        else:
            return commit[:hash_length]

    def get_tags(self):
        """
        Get a dict of all of the tags defined in the repo and their
        corresponding commit hash
        """
        raise NotImplementedError()

    def get_date_from_name(self, name):
        """
        Get a JavaScript timestamp for a particular name.
        """
        raise NotImplementedError()

    def get_branch_commits(self, branch):
        """
        Returns the ordered list (last commit first) of all commits in
        `branch` following first parent in case of merge
        """
        raise NotImplementedError()

    def get_new_branch_commits(self, branches, existing):
        """
        Return a set of new commits on `branches` that are successors of all
        `existing` commits.

        Returns
        -------
        commit_hashes
            New commit hashes
        all_commit_hashes
            All branch commit hashes
        """
        new = set()
        new_commits = []
        for branch in branches:
            for commit in self.get_branch_commits(branch):
                if commit in existing:
                    break
                if commit not in new:
                    new.add(commit)
                    new_commits.append(commit)
        return new_commits

    def filter_date_period(self, commits, period, old_commits=None):
        """
        Pick a subset of `commits` such that the dates are spaced at least
        by `period` seconds.

        If `old_commits` are given, consider also them in the selection,
        but omit them from the returned list.

        The initial part of the selection remains stable if new commits
        with newer dates are added.
        """
        if old_commits is None:
            old_commits = []

        old_commits = set(old_commits)

        items = []
        for commit in set(commits).union(old_commits):
            items.append((self.get_date(commit), commit))
        items.sort()

        # JS date
        period = period * 1000

        # Pick a subset of commits
        prev_date = None
        selected = set()
        for date, commit in items:
            if prev_date is None or date >= prev_date + period:
                prev_date = date
                if commit not in old_commits:
                    selected.add(commit)

        # Preserve original ordering
        return [commit for commit in commits if commit in selected]


class NoRepository(Repo):
    """
    Project installed in the current environment
    """

    dvcs = "none"
    _default_branch = "master"

    def __init__(self, url=None, path=None):
        self.url = None
        self.path = None

    def _raise_error(self):
        raise ValueError("Using the currently installed project version: "
                         "operations requiring repository are not possible")

    def _check_branch(self, branch):
        if branch is not None:
            self._raise_error()

    @classmethod
    def url_match(cls, url):
        return False

    def checkout(self, path, commit_hash):
        self._check_branch(commit_hash)

    def get_date(self, hash):
        if hash is None:
            return datetime.datetime.now(datetime.timezone.utc)
        self._raise_error()

    def get_hashes_from_range(self, range):
        return [None]

    def get_hash_from_name(self, name):
        return None

    def get_hash_from_parent(self, name):
        self._raise_error()

    def get_name_from_hash(self, commit):
        return None

    def get_tags(self):
        return {}

    def get_date_from_name(self, name):
        self._raise_error()

    def pull(self):
        # Nothing to do
        pass

    def get_branch_commits(self, branch):
        self._raise_error()

    def get_new_branch_commits(self, branches, existing):
        self._raise_error()


def get_repo(conf):
    """
    Get a Repo subclass for the given configuration.

    If the configuration does not explicitly specify a repository
    type, it will attempt to automatically determine one from the
    ``conf.repo`` URL.
    """
    if conf.repo is None:
        return None

    if conf.dvcs is not None:
        for cls in util.iter_subclasses(Repo):
            if getattr(cls, 'dvcs') == conf.dvcs:
                return cls(conf.repo, conf.project)
    else:
        for cls in util.iter_subclasses(Repo):
            if cls.url_match(conf.repo):
                return cls(conf.repo, conf.project)

    raise util.UserError(
        f"Can not determine what kind of DVCS to use for URL '{conf.repo}'")

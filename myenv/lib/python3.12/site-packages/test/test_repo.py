# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import tempfile
import shutil
import datetime
from os.path import join

import pytest

from asv import config, repo, util

try:
    import hglib
except ImportError:
    hglib = None

from . import tools

pytestmark = pytest.mark.skipif(tools.HAS_PYPY, reason="Randomly times out on pypy")


def _test_generic_repo(conf, tmpdir, hash_range, main, branch, is_remote=False):
    workcopy_dir = tempfile.mkdtemp(dir=tmpdir, prefix="workcopy")
    os.rmdir(workcopy_dir)

    # check mirroring fails early if *mirror_dir* exists but is not
    # a mirror
    if is_remote:
        if os.path.isdir(conf.project):
            shutil.rmtree(conf.project)
        os.makedirs(join(conf.project, 'hello'))
        with pytest.raises(util.UserError):
            r = repo.get_repo(conf)
        shutil.rmtree(conf.project)

    # basic checkouts
    r = repo.get_repo(conf)

    # Subrepo creation
    r.checkout(workcopy_dir, main)
    assert os.path.exists(join(workcopy_dir, "setup.py"))

    for filename in ("README", "untracked"):
        with open(join(workcopy_dir, filename), "wb") as fd:
            fd.write(b"foo")

    # After checkout the subrepo has been cleaned
    r.checkout(workcopy_dir, branch)
    assert not os.path.exists(join(workcopy_dir, "untracked"))
    with open(join(workcopy_dir, "README"), "rb") as fd:
        data = fd.read(33)
        assert data == b"This is the asv_test_repo project"

    r.checkout(workcopy_dir, main)

    # check recovering from corruption
    for pth in ['.hg', '.git']:
        pth = os.path.join(workcopy_dir, pth)
        if os.path.isdir(pth):
            shutil.rmtree(pth)
    r.checkout(workcopy_dir, main)

    hashes = r.get_hashes_from_range(hash_range)
    assert len(hashes) == 4

    dates = [r.get_date(hash) for hash in hashes]
    assert dates == sorted(dates)[::-1]

    tags = r.get_tags()
    for tag in tags:
        r.get_date_from_name(tag)


def _test_branches(conf, branch_commits, require_describe=False):
    r = repo.get_repo(conf)

    assert len(conf.branches) == 2

    for branch in conf.branches:
        commits = r.get_branch_commits(branch)

        for commit in branch_commits[branch]:
            assert commit in commits

            name = r.get_name_from_hash(commit)
            if require_describe:
                assert name is not None
            if name is not None:
                assert r.get_hash_from_name(name) == commit
                assert name in r.get_decorated_hash(commit)


def test_repo_git(tmpdir):
    tmpdir = str(tmpdir)

    dvcs = tools.generate_test_repo(tmpdir, list(range(10)),
                                    dvcs_type='git',
                                    extra_branches=[(f'{util.git_default_branch()}~4',
                                                     'some-branch', [11, 12, 13])])

    mirror_dir = join(tmpdir, "repo")

    def test_it(is_remote=False):
        conf = config.Config()

        conf.project = mirror_dir
        conf.repo = dvcs.path
        _test_generic_repo(conf, tmpdir,
                           f'{util.git_default_branch()}~4..{util.git_default_branch()}',
                           f'{util.git_default_branch()}', 'tag5',
                           is_remote=is_remote)

        conf.branches = [f'{util.git_default_branch()}', 'some-branch']
        branch_commits = {
            f'{util.git_default_branch()}':
            [dvcs.get_hash(f'{util.git_default_branch()}'),
             dvcs.get_hash(f'{util.git_default_branch()}~6')],
            'some-branch': [dvcs.get_hash('some-branch'),
                            dvcs.get_hash('some-branch~6')]
        }
        _test_branches(conf, branch_commits, require_describe=True)

    test_it()

    # local repo, so it should not not have cloned it
    assert not os.path.isdir(mirror_dir)

    # try again, pretending the repo is not local
    from asv.plugins.git import Git
    old_local_method = Git.is_local_repo
    old_url_match = Git.url_match
    try:
        Git.is_local_repo = classmethod(lambda cls, path:
                                        path != dvcs.path and
                                        old_local_method(path))
        Git.url_match = classmethod(lambda cls, url: os.path.isdir(url))
        test_it(is_remote=True)
        assert os.path.isdir(mirror_dir)
    finally:
        Git.is_local_repo = old_local_method
        Git.url_match = old_url_match


def test_repo_git_annotated_tag_date(tmpdir):
    tmpdir = str(tmpdir)

    dvcs = tools.generate_test_repo(tmpdir, list(range(5)), dvcs_type='git')

    conf = config.Config()
    conf.project = 'sometest'
    conf.repo = dvcs.path

    r = repo.get_repo(conf)
    d1 = r.get_date('tag1')
    d2 = r.get_date(r.get_hash_from_name('tag1'))
    assert d1 == d2


@pytest.mark.skipif(hglib is None,
                    reason="needs hglib")
def test_repo_hg(tmpdir):
    tmpdir = str(tmpdir)

    conf = config.Config()

    dvcs = tools.generate_test_repo(tmpdir, list(range(10)), dvcs_type='hg',
                                    extra_branches=[('default~4', 'somebranch', [11, 12, 13])])

    mirror_dir = join(tmpdir, "repo")

    def test_it(is_remote=False):
        conf.project = mirror_dir
        conf.repo = dvcs.path
        _test_generic_repo(conf, tmpdir, hash_range="reverse(default~3::default)",
                           main="default", branch="tag5",
                           is_remote=is_remote)

        conf.branches = ['default', 'somebranch']
        branch_commits = {
            'default': [dvcs.get_hash('default'), dvcs.get_hash('default~6')],
            'somebranch': [dvcs.get_hash('somebranch'), dvcs.get_hash('somebranch~6')]
        }
        _test_branches(conf, branch_commits)

    test_it()

    # local repo, so it should not not have cloned it
    assert not os.path.isdir(mirror_dir)

    # try again, pretending the repo is not local
    from asv.plugins.mercurial import Hg
    old_local_method = Hg.is_local_repo
    old_url_match = Hg.url_match
    try:
        Hg.is_local_repo = classmethod(lambda cls, path:
                                       path != dvcs.path and
                                       old_local_method(path))
        Hg.url_match = classmethod(lambda cls, url: os.path.isdir(url))
        test_it(is_remote=True)
        assert os.path.isdir(mirror_dir)
    finally:
        Hg.is_local_repo = old_local_method
        Hg.url_match = old_url_match


def test_get_branch_commits(two_branch_repo_case):
    # Test that get_branch_commits() return an ordered list of commits (last
    # first) and follow first parent in case of merge
    dvcs, main, r, conf = two_branch_repo_case
    expected = {
        main: [
            "Revision 6",
            "Revision 4",
            "Merge stable",
            "Revision 3",
            "Revision 1",
        ],
        "stable": [
            "Revision 5",
            "Merge master",
            "Revision 2",
            "Revision 1",
        ],
    }
    for branch in conf.branches:
        commits = [
            dvcs.get_commit_message(commit_hash)
            for commit_hash in r.get_branch_commits(branch)
        ]
        assert commits == expected[branch]


@pytest.mark.parametrize("existing, expected", [
    # No existing commit, we expect all commits in commit order,
    # master branch first
    ([], ["Revision 6", "Revision 4", "Merge stable", "Revision 3",
          "Revision 1", "Revision 5", "Merge master", "Revision 2"]),

    # New commits on each branch
    (["Revision 4", "Merge master"], ["Revision 6", "Revision 5"]),

    # No new commits
    (["Revision 6", "Revision 5"], []),

    # Missing all commits on one branch (case of new branch added in config)
    (["Revision 6"], ["Revision 5", "Merge master", "Revision 2", "Revision 1"]),
], ids=["all", "new", "no-new", "new-branch-added-in-config"])
def test_get_new_branch_commits(two_branch_repo_case, existing, expected):
    dvcs, main, r, conf = two_branch_repo_case

    existing_commits = set()
    for branch in conf.branches:
        for commit in r.get_branch_commits(branch):
            message = dvcs.get_commit_message(commit)
            if message in existing:
                existing_commits.add(commit)

    assert len(existing_commits) == len(existing)

    new_commits = r.get_new_branch_commits(conf.branches, existing_commits)
    commits = [dvcs.get_commit_message(commit) for commit in new_commits]
    assert commits == expected


def test_git_submodule(tmpdir):
    tmpdir = str(tmpdir)

    # State 0 (no submodule)
    dvcs = tools.generate_test_repo(tmpdir, values=[0], dvcs_type='git')
    sub_dvcs = tools.generate_test_repo(tmpdir, values=[0], dvcs_type='git')
    ssub_dvcs = tools.generate_test_repo(tmpdir, values=[0], dvcs_type='git')
    commit_hash_0 = dvcs.get_hash(f"{util.git_default_branch()}")

    # State 1 (one submodule)
    dvcs.run_git(['-c', 'protocol.file.allow=always',
                  'submodule', 'add', sub_dvcs.path, 'sub1'])
    dvcs.commit('Add sub1')
    commit_hash_1 = dvcs.get_hash(f"{util.git_default_branch()}")

    # State 2 (one submodule with sub-submodule)
    dvcs.run_git(['-c', 'protocol.file.allow=always',
                  'submodule', 'update', '--init'])
    sub1_dvcs = tools.Git(join(dvcs.path, 'sub1'))
    sub_dvcs.run_git(['-c', 'protocol.file.allow=always',
                      'submodule', 'add', ssub_dvcs.path, 'ssub1'])
    sub_dvcs.commit('Add sub1')
    sub1_dvcs.run_git(['pull'])
    dvcs.run_git(['add', 'sub1'])
    dvcs.commit('Update sub1')
    sub1_hash_2 = sub1_dvcs.get_hash(f"{util.git_default_branch()}")
    commit_hash_2 = dvcs.get_hash(f"{util.git_default_branch()}")

    # State 3 (one submodule; sub-submodule removed)
    sub_dvcs.run_git(['rm', '-f', 'ssub1'])
    sub_dvcs.commit('Remove ssub1')
    sub1_dvcs.run_git(['pull'])
    dvcs.run_git(['add', 'sub1'])
    dvcs.commit('Update sub1 again')
    commit_hash_3 = dvcs.get_hash(f"{util.git_default_branch()}")

    # State 4 (back to one submodule with sub-submodule)
    sub1_dvcs.run_git(['checkout', sub1_hash_2])
    dvcs.run_git(['add', 'sub1'])
    dvcs.commit('Update sub1 3rd time')
    commit_hash_4 = dvcs.get_hash(f"{util.git_default_branch()}")

    # State 5 (remove final submodule)
    dvcs.run_git(['rm', '-f', 'sub1'])
    dvcs.commit('Remove sub1')
    commit_hash_5 = dvcs.get_hash(f"{util.git_default_branch()}")
    # Verify clean operation
    conf = config.Config()
    conf.branches = [None]
    conf.repo = dvcs.path
    conf.project = join(tmpdir, "repo")
    r = repo.get_repo(conf)

    checkout_dir = join(tmpdir, "checkout")

    # State 0
    r.checkout(checkout_dir, commit_hash_0)
    assert os.path.isfile(join(checkout_dir, 'README'))
    assert not os.path.exists(join(checkout_dir, 'sub1'))

    # State 1
    r.checkout(checkout_dir, commit_hash_1)
    assert os.path.isfile(join(checkout_dir, 'sub1', 'README'))
    assert not os.path.exists(join(checkout_dir, 'sub1', 'ssub1'))

    # State 2
    r.checkout(checkout_dir, commit_hash_2)
    assert os.path.isfile(join(checkout_dir, 'sub1', 'ssub1', 'README'))

    # State 3
    r.checkout(checkout_dir, commit_hash_3)
    assert os.path.isfile(join(checkout_dir, 'sub1', 'README'))
    assert not os.path.exists(join(checkout_dir, 'sub1', 'ssub1'))

    # State 4
    r.checkout(checkout_dir, commit_hash_4)
    assert os.path.isfile(join(checkout_dir, 'sub1', 'ssub1', 'README'))

    # State 4 (check clean -fdx runs in sub-sub modules)
    garbage_filename = join(checkout_dir, 'sub1', 'ssub1', '.garbage')
    util.write_json(garbage_filename, {})
    assert os.path.isfile(garbage_filename)
    r.checkout(checkout_dir, commit_hash_4)
    assert not os.path.isfile(garbage_filename)

    # State 5
    r.checkout(checkout_dir, commit_hash_5)
    assert os.path.isfile(join(checkout_dir, 'README'))
    assert not os.path.isdir(join(checkout_dir, 'sub1'))


@pytest.mark.parametrize('dvcs_type', [
    "git",
    pytest.param("hg", marks=pytest.mark.skipif(hglib is None, reason="needs hglib"))
])
def test_root_ceiling(dvcs_type, tmpdir):
    # Check that git/hg does not try to look for repository in parent
    # directories.
    tmpdir = str(tmpdir)
    dvcs1 = tools.generate_repo_from_ops(tmpdir, dvcs_type, [("commit", 1)])
    dvcs2 = tools.generate_repo_from_ops(tmpdir, dvcs_type, [("commit", 2)])
    commit1 = dvcs1.get_branch_hashes()[0]
    commit2 = dvcs2.get_branch_hashes()[0]

    conf = config.Config()
    conf.branches = []
    conf.dvcs = dvcs_type
    conf.project = join(tmpdir, "repo")
    conf.repo = dvcs1.path

    r = repo.get_repo(conf)

    # Checkout into a subdir inside another repository
    workcopy_dir = join(dvcs2.path, "workcopy")
    r.checkout(workcopy_dir, commit1)

    # Corrupt the checkout
    for pth in ['.hg', '.git']:
        pth = os.path.join(workcopy_dir, pth)
        if os.path.isdir(pth):
            shutil.rmtree(pth)

    # Operation must fail (commit2 is not in dvcs1), not use the
    # parent repository
    with pytest.raises(Exception):
        r.checkout(workcopy_dir, commit2)


@pytest.mark.parametrize('dvcs_type', [
    "git",
    pytest.param("hg", marks=pytest.mark.skipif(hglib is None, reason="needs hglib"))
])
def test_no_such_name_error(dvcs_type, tmpdir):
    tmpdir = str(tmpdir)
    dvcs = tools.generate_test_repo(tmpdir, values=[0], dvcs_type=dvcs_type)

    conf = config.Config()
    conf.branches = []
    conf.dvcs = dvcs_type
    conf.project = "project"
    conf.repo = dvcs.path

    r = repo.get_repo(conf)

    # Check that NoSuchNameError error gets raised correctly
    assert r.get_hash_from_name(None) == dvcs.get_hash(r._default_branch)
    with pytest.raises(repo.NoSuchNameError):
        r.get_hash_from_name("badbranch")

    if dvcs_type == "git":
        # Corrupted repository/etc should not give NoSuchNameError
        util.long_path_rmtree(join(dvcs.path, ".git"))
        with pytest.raises(Exception) as excinfo:
            r.get_hash_from_name(None)
        assert excinfo.type not in (AssertionError, repo.NoSuchNameError)
    elif dvcs_type == "hg":
        # hglib seems to do some caching, so this doesn't work
        pass


@pytest.mark.parametrize('dvcs_type', [
    "git",
    pytest.param("hg", marks=pytest.mark.skipif(hglib is None, reason="needs hglib"))
])
def test_filter_date_period(tmpdir, dvcs_type):
    tmpdir = str(tmpdir)

    dates = [
        datetime.datetime(2001, 1, 1, tzinfo=datetime.timezone.utc),
        datetime.datetime(2001, 1, 2, tzinfo=datetime.timezone.utc),
        datetime.datetime(2001, 1, 8, tzinfo=datetime.timezone.utc)
    ]

    dvcs = tools.generate_repo_from_ops(
        tmpdir, dvcs_type,
        [("commit", j, dates[j]) for j in range(len(dates))])
    commits = dvcs.get_branch_hashes()[::-1]
    assert len(commits) == len(dates)

    conf = config.Config()
    conf.dvcs = dvcs_type
    conf.repo = dvcs.path
    r = repo.get_repo(conf)

    # Basic filtering
    weekly_commits = r.filter_date_period(commits, 60 * 60 * 24 * 7)
    assert weekly_commits == [commits[0], commits[2]]

    daily_commits = r.filter_date_period(commits, 60 * 60 * 24)
    assert daily_commits == commits

    # Test with old_commits specified
    monthly_commits = r.filter_date_period(commits[1:], 60 * 60 * 24 * 30, commits[:1])
    assert monthly_commits == []

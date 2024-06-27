# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

import pytest

from asv import util

from . import tools


@pytest.mark.skipif(tools.HAS_PYPY, reason="Flaky on pypy")
@pytest.mark.parametrize("rewrite", [False, True], ids=["no-rewrite", "rewrite"])
def test_gh_pages(rewrite, tmpdir, generate_result_dir, monkeypatch):
    tmpdir = os.path.abspath(str(tmpdir))

    monkeypatch.setenv(str('EMAIL'), str('test@asv'))
    monkeypatch.setenv(str('GIT_COMMITTER_NAME'), str('asv test'))
    monkeypatch.setenv(str('GIT_AUTHOR_NAME'), str('asv test'))

    conf, repo, commits = generate_result_dir([1, 2, 3, 4])

    dvcs_dir = os.path.join(tmpdir, 'repo1')
    dvcs_dir2 = os.path.join(tmpdir, 'repo2')

    os.makedirs(dvcs_dir)

    os.chdir(dvcs_dir)

    dvcs = tools.Git(dvcs_dir)
    dvcs.init()

    open(os.path.join(dvcs_dir, 'dummy'), 'wb').close()

    dvcs.add('dummy')
    dvcs.commit('Initial commit')

    if rewrite:
        rewrite_args = ("--rewrite",)
    else:
        rewrite_args = ()

    # Check with no existing gh-pages branch, no push
    tools.run_asv_with_conf(conf, "gh-pages", "--no-push", *rewrite_args)
    dvcs.checkout('gh-pages')
    assert os.path.isfile(os.path.join(dvcs_dir, 'index.html'))
    assert len(dvcs.run_git(['rev-list', 'gh-pages']).splitlines()) == 1
    dvcs.checkout(f"{util.git_default_branch()}")
    assert not os.path.isfile(os.path.join(dvcs_dir, 'index.html'))

    # Check with existing (and checked out) gh-pages branch, with no changes
    tools.run_asv_with_conf(conf, "gh-pages", "--no-push", *rewrite_args)
    dvcs.checkout('gh-pages')
    assert os.path.isfile(os.path.join(dvcs_dir, 'index.html'))
    if rewrite:
        assert len(dvcs.run_git(['rev-list', 'gh-pages']).splitlines()) == 1
    else:
        # Timestamp may have changed
        assert len(dvcs.run_git(['rev-list', 'gh-pages']).splitlines()) <= 2
    dvcs.checkout(f"{util.git_default_branch()}")

    # Check with existing (not checked out) gh-pages branch, with some changes
    benchmarks_json = os.path.join(conf.results_dir, 'benchmarks.json')
    data = util.load_json(benchmarks_json)
    data['time_func']['pretty_name'] = 'something changed'
    util.write_json(benchmarks_json, data)

    prev_len = len(dvcs.run_git(['rev-list', 'gh-pages']).splitlines())
    tools.run_asv_with_conf(conf, "gh-pages", "--no-push", *rewrite_args)
    if not rewrite:
        assert len(dvcs.run_git(['rev-list', 'gh-pages']).splitlines()) == prev_len + 1
    else:
        assert len(dvcs.run_git(['rev-list', 'gh-pages']).splitlines()) == prev_len

    # Check that the push option works
    dvcs.run_git(['branch', '-D', 'gh-pages'])
    dvcs.run_git(['clone', dvcs_dir, dvcs_dir2])

    os.chdir(dvcs_dir2)
    tools.run_asv_with_conf(conf, "gh-pages", *rewrite_args)

    os.chdir(dvcs_dir)
    dvcs.checkout('gh-pages')
    assert os.path.isfile(os.path.join(dvcs_dir, 'index.html'))

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import shutil
import hashlib

from asv import util

from . import tools


def test_update_simple(monkeypatch, generate_result_dir):
    conf, repo, commits = generate_result_dir(5 * [1] + 5 * [10])

    basedir = os.path.abspath(os.path.dirname(conf.results_dir))
    local = os.path.abspath(os.path.dirname(__file__))

    shutil.copyfile(os.path.join(local, 'asv-machine.json'),
                    os.path.join(basedir, 'asv-machine.json'))
    machine_file = 'asv-machine.json'

    conf_values = {}
    for key in ['results_dir', 'html_dir', 'repo', 'project', 'branches']:
        conf_values[key] = getattr(conf, key)

    util.write_json(os.path.join(basedir, 'asv.conf.json'), conf_values,
                    api_version=1)

    # Check renaming of long result files
    machine_dir = os.path.join(basedir, 'results', 'tarzan')

    result_fns = [fn for fn in sorted(os.listdir(machine_dir))
                  if fn != 'machine.json']
    long_result_fn = 'abbacaca-' + 'a' * 128 + '.json'
    hash_result_fn = ('abbacaca-env-' + hashlib.md5(b'a' * 128).hexdigest() + '.json')

    shutil.copyfile(os.path.join(machine_dir, result_fns[0]),
                    os.path.join(machine_dir, long_result_fn))

    old_env_name = util.load_json(os.path.join(machine_dir, result_fns[0]))['env_name']

    # Should succeed
    monkeypatch.chdir(basedir)
    tools.run_asv_with_conf(conf, "update", _machine_file=machine_file)

    # Check file rename
    items = [fn for fn in sorted(os.listdir(machine_dir))
             if fn != 'machine.json']
    assert long_result_fn.lower() not in [x.lower() for x in items]
    assert hash_result_fn.lower() in [x.lower() for x in items]

    # Check env name is preserved
    new_env_name = util.load_json(os.path.join(machine_dir, items[0]))['env_name']
    assert old_env_name == new_env_name

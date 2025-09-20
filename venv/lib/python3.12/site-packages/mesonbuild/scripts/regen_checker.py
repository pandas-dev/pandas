# Copyright 2015-2016 The Meson development team

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from __future__ import annotations

import sys, os
import pickle, subprocess
import typing as T
from ..coredata import CoreData
from ..backend.backends import RegenInfo
from ..mesonlib import OptionKey

# This could also be used for XCode.

def need_regen(regeninfo: RegenInfo, regen_timestamp: float) -> bool:
    for i in regeninfo.depfiles:
        curfile = os.path.join(regeninfo.build_dir, i)
        curtime = os.stat(curfile).st_mtime
        if curtime > regen_timestamp:
            return True
    # The timestamp file gets automatically deleted by MSBuild during a 'Clean' build.
    # We must make sure to recreate it, even if we do not regenerate the solution.
    # Otherwise, Visual Studio will always consider the REGEN project out of date.
    print("Everything is up-to-date, regeneration of build files is not needed.")
    from ..backend.vs2010backend import Vs2010Backend
    Vs2010Backend.touch_regen_timestamp(regeninfo.build_dir)
    return False

def regen(regeninfo: RegenInfo, meson_command: T.List[str], backend: str) -> None:
    cmd = meson_command + ['--internal',
                           'regenerate',
                           regeninfo.build_dir,
                           regeninfo.source_dir,
                           '--backend=' + backend]
    subprocess.check_call(cmd)

def run(args: T.List[str]) -> int:
    private_dir = args[0]
    dumpfile = os.path.join(private_dir, 'regeninfo.dump')
    coredata_file = os.path.join(private_dir, 'coredata.dat')
    with open(dumpfile, 'rb') as f:
        regeninfo = pickle.load(f)
        assert isinstance(regeninfo, RegenInfo)
    with open(coredata_file, 'rb') as f:
        coredata = pickle.load(f)
        assert isinstance(coredata, CoreData)
    backend = coredata.get_option(OptionKey('backend'))
    assert isinstance(backend, str)
    regen_timestamp = os.stat(dumpfile).st_mtime
    if need_regen(regeninfo, regen_timestamp):
        regen(regeninfo, coredata.meson_command, backend)
    return 0

if __name__ == '__main__':
    sys.exit(run(sys.argv[1:]))

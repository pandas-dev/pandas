# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

import pyarrow.fs


def application(env, start_response):
    path = env['PATH_INFO']
    members = path.split('/')
    assert members[0] == ''
    assert len(members) >= 2
    root = members[1]
    if root == 's3':
        # See test_fs::test_uwsgi_integration
        start_response('200 OK', [('Content-Type', 'text/html')])
        # flake8: noqa
        fs = pyarrow.fs.S3FileSystem()
        return [b"Hello World\n"]
    else:
        start_response('404 Not Found', [('Content-Type', 'text/html')])
        return [f"Path {path!r} not found\n".encode()]

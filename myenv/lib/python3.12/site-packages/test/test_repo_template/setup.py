import os

from setuptools import setup

# This file is a template, and will be rendered before executed.
# So the double curly brackets will become single after rendering, and
# when executed, this will work as expected
content = 'env = {{}}\n'.format(repr(dict(os.environ)))  # noqa W291 see above comment
with open('asv_test_repo/build_time_env.py', 'w') as f:
    f.write(content)

setup(name='asv_test_repo',
      version="{version}",
      packages=['asv_test_repo'],

      # The following forces setuptools to generate .egg-info directory,
      # which causes problems in test_environment.py:test_install_success
      include_package_data=True,
      )

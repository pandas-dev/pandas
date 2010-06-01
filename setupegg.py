#!/usr/bin/env python
"""
A setup.py script to use setuptools.
"""

from setup import *
# execfile('setup.py')
setup(name=DISTNAME,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      description=DESCRIPTION,
      license=LICENSE,
      url=URL,
      download_url=DOWNLOAD_URL,
      long_description=LONG_DESCRIPTION,
      classifiers=CLASSIFIERS,
      zip_safe=False,
      platforms='any',
      configuration=configuration)

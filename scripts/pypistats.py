#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calculates the total number of downloads that a particular PyPI package has
received across all versions tracked by PyPI
"""

from datetime import datetime
import locale
import sys
import xmlrpclib
import pandas as pd

locale.setlocale(locale.LC_ALL, '')


class PyPIDownloadAggregator(object):

    def __init__(self, package_name, include_hidden=True):
        self.package_name = package_name
        self.include_hidden = include_hidden
        self.proxy = xmlrpclib.Server('http://pypi.python.org/pypi')
        self._downloads = {}

    @property
    def releases(self):
        """Retrieves the release number for each uploaded release"""

        result = self.proxy.package_releases(self.package_name,
                                             self.include_hidden)

        if len(result) == 0:
            # no matching package--search for possibles, and limit to 15
            # results
            results = self.proxy.search({
                'name': self.package_name,
                'description': self.package_name
            }, 'or')[:15]

            # make sure we only get unique package names
            matches = []
            for match in results:
                name = match['name']
                if name not in matches:
                    matches.append(name)

            # if only one package was found, return it
            if len(matches) == 1:
                self.package_name = matches[0]
                return self.releases

            error = """No such package found: %s

Possible matches include:
%s
""" % (self.package_name, '\n'.join('\t- %s' % n for n in matches))

            sys.exit(error)

        return result

    def get_downloads(self):
        """Calculate the total number of downloads for the package"""
        downloads = {}
        for release in self.releases:
            urls = self.proxy.release_urls(self.package_name, release)
            urls = pd.DataFrame(urls)
            urls['version'] = release
            downloads[release] = urls

        return pd.concat(downloads, ignore_index=True)

if __name__ == '__main__':
    agg = PyPIDownloadAggregator('pandas')

    data = agg.get_downloads()

    to_omit = ['0.2b1', '0.2beta']

    isostrings = data['upload_time'].map(lambda x: x.value)
    data['upload_time'] = pd.to_datetime(isostrings)

    totals = data.groupby('version').downloads.sum()
    rollup = {'0.8.0rc1': '0.8.0',
              '0.8.0rc2': '0.8.0',
              '0.3.0.beta': '0.3.0',
              '0.3.0.beta2': '0.3.0'}
    downloads = totals.groupby(lambda x: rollup.get(x, x)).sum()

    first_upload = data.groupby('version').upload_time.min()

    result = pd.DataFrame({'downloads': totals,
                           'release_date': first_upload})
    result = result.sort('release_date')
    result = result.drop(to_omit + list(rollup.keys()))
    result.index.name = 'release'

    by_date = result.reset_index().set_index('release_date').downloads
    dummy = pd.Series(index=pd.DatetimeIndex([datetime(2012, 12, 27)]))
    by_date = by_date.append(dummy).shift(1).fillna(0)

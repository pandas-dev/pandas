from io import open

import numpy as np
import pytest

from pandas import DataFrame, Index, MultiIndex


@pytest.fixture
def expected_html(datapath):
    def _expected_html(name):
        """
        Read HTML file from formats data directory.

        Parameters
        ----------
        name : str
            The name of the HTML file without the suffix.

        Returns
        -------
        str : contents of HTML file.
        """
        filename = '.'.join([name, 'html'])
        filepath = datapath('io', 'formats', 'data', 'html', filename)
        with open(filepath, encoding='utf-8') as f:
            html = f.read()
        return html.rstrip()
    return _expected_html


@pytest.fixture(params=[True, False])
def index_names_fixture(request):
    return request.param


@pytest.fixture(params=[True, False])
def header_fixture(request):
    return request.param


@pytest.fixture(params=[True, False])
def index_fixture(request):
    return request.param


@pytest.fixture(params=['standard', 'multi'])
def index_type_fixture(request):
    return request.param


@pytest.fixture(params=['standard', 'multi'])
def columns_index_type_fixture(request):
    return request.param


@pytest.fixture(params=['unnamed', 'named'])
def index_naming_fixture(request):
    return request.param


@pytest.fixture(params=['unnamed', 'named'])
def columns_index_naming_fixture(request):
    return request.param


@pytest.fixture(params=[[0, 1]])
def index_labels_fixture(request):
    return request.param


def _df_index(index_labels, index_naming, index_name):
    if index_naming == 'unnamed':
        return Index(index_labels)
    return Index(index_labels, name=index_name)


@pytest.fixture()
def df_index_fixture(index_labels_fixture, index_naming_fixture):
    return _df_index(index_labels_fixture, index_naming_fixture, 'index.name')


@pytest.fixture()
def df_columns_index_fixture(
        index_labels_fixture, columns_index_naming_fixture):
    return _df_index(index_labels_fixture, columns_index_naming_fixture,
                     'columns.name')


@pytest.fixture(params=[[['a'], ['b', 'c']]])
def multi_index_labels_fixture(request):
    return request.param


def _df_multi_index(multi_index_labels, index_naming, index_names):
    if index_naming == 'unnamed':
        return MultiIndex.from_product(multi_index_labels)
    return MultiIndex.from_product(multi_index_labels, names=index_names)


@pytest.fixture(params=[['index.name.0', 'index.name.1']])
def df_multi_index_fixture(
        request, multi_index_labels_fixture, index_naming_fixture):
    names = request.param
    return _df_multi_index(
        multi_index_labels_fixture, index_naming_fixture, names)


@pytest.fixture(params=[['columns.name.0', 'columns.name.1']])
def df_columns_multi_index_fixture(
        request, multi_index_labels_fixture, columns_index_naming_fixture):
    names = request.param
    return _df_multi_index(
        multi_index_labels_fixture, columns_index_naming_fixture, names)


@pytest.fixture()
def df_indexes_fixture(
        index_type_fixture, df_index_fixture, df_multi_index_fixture):
    if index_type_fixture == 'multi':
        return df_multi_index_fixture
    return df_index_fixture


@pytest.fixture()
def df_columns_indexes_fixture(
        columns_index_type_fixture, df_columns_index_fixture,
        df_columns_multi_index_fixture):
    if columns_index_type_fixture == 'multi':
        return df_columns_multi_index_fixture
    return df_columns_index_fixture


@pytest.fixture(params=[np.zeros((2, 2), dtype=int)])
def df_fixture(request, df_indexes_fixture, df_columns_indexes_fixture):
    data = request.param
    return DataFrame(data, index=df_indexes_fixture,
                     columns=df_columns_indexes_fixture)


@pytest.fixture(params=[None])
def expected_fixture(
        request, expected_html, index_type_fixture, index_naming_fixture,
        columns_index_type_fixture, columns_index_naming_fixture,
        index_fixture, header_fixture, index_names_fixture):
    filename_prefix = request.param
    if not index_fixture:
        index_naming_fixture = 'none'
    else:
        if not index_names_fixture:
            index_naming_fixture = 'unnamed'
        index_naming_fixture = index_naming_fixture + '_' + index_type_fixture

    if not header_fixture:
        columns_index_naming_fixture = 'none'
    else:
        if not index_names_fixture:
            columns_index_naming_fixture = 'unnamed'
        columns_index_naming_fixture = (
            columns_index_naming_fixture + '_' + columns_index_type_fixture)

    filename = '_'.join(['index', index_naming_fixture,
                         'columns', columns_index_naming_fixture])
    if filename_prefix:
        filename = filename_prefix + filename
    return expected_html(filename)

<div align="center">
  <img src="https://github.com/pandas-dev/pandas/blob/master/doc/logo/pandas_logo.png"><br>
</div>

-----------------

# pandas: powerful Python data analysis toolkit

<table>
<tr>
  <td>Latest Release</td>
  <td><img src="https://img.shields.io/pypi/v/pandas.svg" alt="latest release" /></td>
</tr>
  <td></td>
  <td><img src="https://anaconda.org/conda-forge/pandas/badges/version.svg" alt="latest release" /></td>
</tr>
<tr>
  <td>Package Status</td>
  <td><img src="https://img.shields.io/pypi/status/pandas.svg" alt="status" /></td>
</tr>
<tr>
  <td>License</td>
  <td><img src="https://img.shields.io/pypi/l/pandas.svg" alt="license" /></td>
</tr>
<tr>
  <td>Build Status</td>
  <td>
    <a href="https://travis-ci.org/pandas-dev/pandas">
    <img src="https://travis-ci.org/pandas-dev/pandas.svg?branch=master" alt="travis build status" />
    </a>
  </td>
</tr>
<tr>
  <td></td>
  <td>
    <a href="https://circleci.com/gh/pandas-dev/pandas">
    <img src="https://circleci.com/gh/circleci/mongofinil/tree/master.svg?style=shield&circle-token=223d8cafa7b02902c3e150242520af8944e34671" alt="circleci build status" />
    </a>
  </td>
</tr>
<tr>
  <td></td>
  <td>
    <a href="https://ci.appveyor.com/project/pandas-dev/pandas">
    <img src="https://ci.appveyor.com/api/projects/status/86vn83mxgnl4xf1s/branch/master?svg=true" alt="appveyor build status" />
    </a>
  </td>
</tr>
<tr>
  <td>Coverage</td>
  <td><img src="https://codecov.io/github/pandas-dev/pandas/coverage.svg?branch=master" alt="coverage" /></td>
</tr>
<tr>
  <td>Conda</td>
  <td>
    <a href="http://pandas.pydata.org">
    <img src="http://pubbadges.s3-website-us-east-1.amazonaws.com/pkgs-downloads-pandas.png" alt="conda default downloads" />
    </a>
  </td>
</tr>
<tr>
  <td>Conda-forge</td>
  <td>
    <a href="http://pandas.pydata.org">
    <img src="https://anaconda.org/conda-forge/pandas/badges/downloads.svg" alt="conda-forge downloads" />
    </a>
  </td>
</tr>
<tr>
  <td>PyPI</td>
  <td>
    <a href="https://pypi.python.org/pypi/pandas/">
    <img src="https://img.shields.io/pypi/dm/pandas.svg" alt="pypi downloads" />
    </a>
  </td>
</tr>
</table>

[![https://gitter.im/pydata/pandas](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/pydata/pandas?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

## Pandas란?

**pandas**는 "관계형"또는 "레이블이 있는" 데이터 작업을 쉽고 직관적으로 처리하기 위해 설계된, 빠르고 유연하며 표현이 풍부한 데이터 구조를 제공하는 Python 패키지입니다. Pandas는 Python에서 **실제** 데이터 분석을 실용적으로 수행하기 위한 기초적인 high-level 구조 블록을 목표로 합니다. 또한 **어떤 언어로도 사용할 수 있고, 가장 강력하고 유연한 오픈 소스 데이터 분석 및 조작 도구** 라는 더 넓은 목표를 가지고 있으며, 이미 이 목표를 향해 나아가고 있습니다.

## 주요 특징
다음은 Pandas의 장점들 입니다:

  - 부동 소수점 데이터뿐만 아니라 부동 소수점 데이터에서 [**누락된 데이터**][missing-data]('NaN')에 대한 처리가 쉽습니다.
  - 크기 가변성 : DataFrame 및 상위 차원 개체에서 열을 [**삽입 및 삭제**][insertion-deletion]할 수 있습니다.
  - 자동 및 명시적인 [**데이터 정렬**][alignment] : 개체를 레이블 집합에 명시적으로 정렬되거나, 사용자가 레이블을 무시하고 계산을 위해'Series', 'DataFrame' 등의 형태로 데이터를 자동 정렬할 수 있습니다.
  - 강력하고 유연한 [**그룹화**][groupby] 기능으로, 데이터의 집계 및 변환을 위해 데이터 세트에 분할-적용-결합 작업을 수행 할 수 있습니다.
  - 다른 Python 및 NumPy의 데이터 구조로 되어있는 비정형 및 다르게 색인된(differently-indexed) 데이터를 DataFrame 개체로 [**쉽게 변환**][conversion]할 수 있습니다.
  - 지능적인 레이블 기반 [**슬라이싱**][slicing], [**고급 인덱싱**][fancy-indexing] 및 대용량 데이터 세트의 [**부분집합화**][subsetting](subsetting).
  - 직관적인 데이터 [**병합**][merging] 및 [**결합**][joining].
  - 데이터 세트에 대해 유연한 [**재형성**][reshape](reshaping) 및 [**고정**][pivot-table](pivoting).
  - 축의 [**계층적**][mi] 레이블링. (틱당 다중 레이블을 가질 수 있습니다.)
  - [**플랫 파일**][flat-files](CSV 및 delimited), [**엑셀 파일**][excel], [**데이터베이스**][db]에서 데이터 로딩 혹은 초고속 [**HDF5 형식**][hdfstore]에서 데이터 로딩 및 저장을 위한 견고한 IO 도구.
  - [**시계열**][timeseries] 관련 기능 : 날짜 범위 생성 및 빈도 변환, moving window 통계, moving window 선형 회귀, 날짜 이동 및 지연 등.


   [missing-data]: http://pandas.pydata.org/pandas-docs/stable/missing_data.html#working-with-missing-data
   [insertion-deletion]: http://pandas.pydata.org/pandas-docs/stable/dsintro.html#column-selection-addition-deletion
   [alignment]: http://pandas.pydata.org/pandas-docs/stable/dsintro.html?highlight=alignment#intro-to-data-structures
   [groupby]: http://pandas.pydata.org/pandas-docs/stable/groupby.html#group-by-split-apply-combine
   [conversion]: http://pandas.pydata.org/pandas-docs/stable/dsintro.html#dataframe
   [slicing]: http://pandas.pydata.org/pandas-docs/stable/indexing.html#slicing-ranges
   [fancy-indexing]: http://pandas.pydata.org/pandas-docs/stable/indexing.html#advanced-indexing-with-ix
   [subsetting]: http://pandas.pydata.org/pandas-docs/stable/indexing.html#boolean-indexing
   [merging]: http://pandas.pydata.org/pandas-docs/stable/merging.html#database-style-dataframe-joining-merging
   [joining]: http://pandas.pydata.org/pandas-docs/stable/merging.html#joining-on-index
   [reshape]: http://pandas.pydata.org/pandas-docs/stable/reshaping.html#reshaping-and-pivot-tables
   [pivot-table]: http://pandas.pydata.org/pandas-docs/stable/reshaping.html#pivot-tables-and-cross-tabulations
   [mi]: http://pandas.pydata.org/pandas-docs/stable/indexing.html#hierarchical-indexing-multiindex
   [flat-files]: http://pandas.pydata.org/pandas-docs/stable/io.html#csv-text-files
   [excel]: http://pandas.pydata.org/pandas-docs/stable/io.html#excel-files
   [db]: http://pandas.pydata.org/pandas-docs/stable/io.html#sql-queries
   [hdfstore]: http://pandas.pydata.org/pandas-docs/stable/io.html#hdf5-pytables
   [timeseries]: http://pandas.pydata.org/pandas-docs/stable/timeseries.html#time-series-date-functionality

## Pandas 얻기
소스 코드는 현재 다음 GitHub 사이트에서 호스팅됩니다 :
http://github.com/pandas-dev/pandas

최신 버전의 이진 설치 프로그램은 [Python 패키지 Index](http://pypi.python.org/pypi/pandas/) 및 conda에서 구할 수 있습니다.

```sh
# conda
conda install pandas
```

```sh
# or PyPI
pip install pandas
```

## Dependencies
- [NumPy](http://www.numpy.org): 1.7.0 or higher
- [python-dateutil](http://labix.org/python-dateutil): 1.5 or higher
- [pytz](http://pytz.sourceforge.net)
    - Needed for time zone support with ``pandas.date_range``

See the [full installation instructions](http://pandas.pydata.org/pandas-docs/stable/install.html#dependencies)
for recommended and optional dependencies.

## Installation from sources
To install pandas from source you need Cython in addition to the normal
dependencies above. Cython can be installed from pypi:

```sh
pip install cython
```

In the `pandas` directory (same one where you found this file after
cloning the git repo), execute:

```sh
python setup.py install
```

or for installing in [development mode](https://pip.pypa.io/en/latest/reference/pip_install.html#editable-installs):

```sh
python setup.py develop
```

Alternatively, you can use `pip` if you want all the dependencies pulled
in automatically (the `-e` option is for installing it in [development
mode](https://pip.pypa.io/en/latest/reference/pip_install.html#editable-installs)):

```sh
pip install -e .
```

See the full instructions for [installing from source](http://pandas.pydata.org/pandas-docs/stable/install.html#installing-from-source).

## License
BSD

## Documentation
The official documentation is hosted on PyData.org: http://pandas.pydata.org/pandas-docs/stable/

The Sphinx documentation should provide a good starting point for learning how
to use the library. Expect the docs to continue to expand as time goes on.

## Background
Work on ``pandas`` started at AQR (a quantitative hedge fund) in 2008 and
has been under active development since then.

## Getting Help

For usage questions, the best place to go to is [StackOverflow](https://stackoverflow.com/questions/tagged/pandas).
Further, general questions and discussions can also take place on the [pydata mailing list](https://groups.google.com/forum/?fromgroups#!forum/pydata).

## Discussion and Development
Most development discussion is taking place on github in this repo. Further, the [pandas-dev mailing list](https://mail.python.org/mailman/listinfo/pandas-dev) can also be used for specialized discussions or design issues, and a [Gitter channel](https://gitter.im/pydata/pandas) is available for quick development related questions.

## Contributing to pandas
All contributions, bug reports, bug fixes, documentation improvements, enhancements and ideas are welcome.

A detailed overview on how to contribute can be found in the **[contributing guide.](http://pandas.pydata.org/pandas-docs/stable/contributing.html)**

If you are simply looking to start working with the pandas codebase, navigate to the [GitHub “issues” tab](https://github.com/pandas-dev/pandas/issues) and start looking through interesting issues. There are a number of issues listed under [Docs](https://github.com/pandas-dev/pandas/issues?labels=Docs&sort=updated&state=open) and [Difficulty Novice](https://github.com/pandas-dev/pandas/issues?q=is%3Aopen+is%3Aissue+label%3A%22Difficulty+Novice%22) where you could start out.

Or maybe through using pandas you have an idea of your own or are looking for something in the documentation and thinking ‘this can be improved’...you can do something about it!

Feel free to ask questions on the [mailing list](https://groups.google.com/forum/?fromgroups#!forum/pydata) or on [Gitter](https://gitter.im/pydata/pandas).

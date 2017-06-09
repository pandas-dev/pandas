<div align="center">
  <img src="https://github.com/pandas-dev/pandas/blob/master/doc/logo/pandas_logo.png"><br>
</div>

-----------------

# pandas : 강력한 파이썬 데이터 분석 툴킷

<table>
<tr>
  <td>최신 버전</td>
  <td><img src="https://img.shields.io/pypi/v/pandas.svg" alt="latest release" /></td>
</tr>
  <td></td>
  <td><img src="https://anaconda.org/conda-forge/pandas/badges/version.svg" alt="latest release" /></td>
</tr>
<tr>
  <td>패키지 상태</td>
  <td><img src="https://img.shields.io/pypi/status/pandas.svg" alt="status" /></td>
</tr>
<tr>
  <td>라이센스</td>
  <td><img src="https://img.shields.io/pypi/l/pandas.svg" alt="license" /></td>
</tr>
<tr>
  <td>빌드 상태</td>
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
  <td>적용 범위</td>
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
  - 자동 및 명시적인 [**데이터 정렬**][alignment] : 개체를 레이블 집합에 명시적으로 정렬되거나, 사용자가 레이블을 무시하고 계산을 위해 'Series', 'DataFrame' 등의 형태로 데이터를 자동 정렬할 수 있습니다.
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
소스 코드는 현재 다음 GitHub 사이트에 호스팅됩니다 :
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

## 종속 프로그램
- [NumPy](http://www.numpy.org): 1.7.0 또는 상위 버전
- [python-dateutil](http://labix.org/python-dateutil): 1.5 또는 상위 버전
- [pytz](http://pytz.sourceforge.net)
    - ``pandas.date_range``를 지원하는 시간대(time zone)가 필요합니다.

권장 혹은 선택적 종속 프로그램은 [전체설치지침](http://pandas.pydata.org/pandas-docs/stable/install.html#dependencies)을 참조하십시오.


## 소스에서 설치
소스로 Pandas를 설치하려면 위의 일반적인 종속 프로그램 외에도 Cython이 필요합니다. Cython은 pypi에서 다음과 같이 설치할 수 있습니다 :

```sh
pip install cython
```

`pandas` 디렉토리(git repo를 복제 한 후, 이 파일을 찾은 디렉토리)에서 다음을 실행하십시오 :

```sh
python setup.py install
```

또는 [개발 모드](https://pip.pypa.io/en/latest/reference/pip_install.html#editable-installs)로 설치할 수 있습니다 :

```sh
python setup.py develop
```

또는, 모든 종속 프로그램을 자동적으로 설치하고 싶다면 `pip`을 사용할 수 있습니다 (`-e` 옵션은 [개발 모드](https://pip.pypa.io/en/latest/reference/pip_install.html#editable-installs)로 설치하기 위함입니다):

```sh
pip install -e .
```

Windows에서는 MinGW를 설치하고, 이를 실행해야합니다 :

```sh
python setup.py build --compiler=mingw32
python setup.py install
```

http://pandas.pydata.org/ 에서 더 많은 정보를 확인하십시오.

## 라이센스
BSD

## 문서
공식 문서는 PyData.org에 호스팅됩니다 : http://pandas.pydata.org/

Sphinx 문서는 라이브러리 사용법을 배우기 위한 좋은 출발점을 제공해야 합니다. 시간이 지남에 따라 문서가 계속 확장 될 것으로 기대합니다.

## 배경
2008년 AQR(A Quantitative Hedge Fund)에서 시작된 ''Pandas''에 대한 작업은 그 후에도 적극적으로 발전해 왔습니다.

## 토론 및 개발
Pandas 개발은 여러 과학관련 파이썬 프로젝트와 관련이 있기 때문에 질의는 SciPy-user 메일링 리스트에서 환영합니다. 특별한 토론 또는 디자인 문제는 PyData 메일링 리스트 / Google 그룹에서 발생되길 권장합니다 :

https://groups.google.com/forum/#!forum/pydata

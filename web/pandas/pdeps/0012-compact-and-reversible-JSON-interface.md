# PDEP-12: Compact and reversible JSON interface

- Created: 16 June 2023
- Status: Rejected
- Discussion:
    [#53252](https://github.com/pandas-dev/pandas/issues/53252)
    [#55038](https://github.com/pandas-dev/pandas/issues/55038)
- Author: [Philippe THOMY](https://github.com/loco-philippe)
- Revision: 3

[TOC]

## Abstract

### Problem description

The `dtype` and "Python type" are not explicitly taken into account in the current JSON interface.

So, the JSON interface is not always reversible and has inconsistencies related to the consideration of the `dtype`.

Another consequence is the partial application of the Table Schema specification in the `orient="table"` option (6 Table Schema data types are taken into account out of the 24 defined).

Some JSON-interface problems are detailed in the [linked NoteBook](https://nbviewer.org/github/loco-philippe/ntv-pandas/blob/main/example/example_json_pandas.ipynb#Current-Json-interface)

### Feature Description

To have a simple, compact and reversible solution, I propose to use the [JSON-NTV format (Named and Typed Value)](https://github.com/loco-philippe/NTV#readme) - which integrates the notion of type - and its JSON-TAB variation for tabular data (the JSON-NTV format is defined in an [IETF Internet-Draft](https://datatracker.ietf.org/doc/draft-thomy-json-ntv/) (not yet an RFC !!) ).

This solution allows to include a large number of types (not necessarily pandas `dtype`) which allows to have:

- a Table Schema JSON interface (`orient="table"`) which respects the Table Schema specification (going from 6 types to 20 types),
- a global JSON interface for all pandas data formats.

#### Global JSON interface example

In the example below, a DataFrame with several data types is converted to JSON.

The DataFrame resulting from this JSON is identical to the initial DataFrame (reversibility).

With the existing JSON interface, this conversion is not possible.

This example uses `ntv_pandas` module defined in the [ntv-pandas repository](https://github.com/loco-philippe/ntv-pandas#readme).

Data example:

```python
In [1]: from shapely.geometry import Point
        from datetime import date
        import pandas as pd
        import ntv_pandas as npd

In [2]: data = {'index':           [100, 200, 300, 400, 500, 600],
                'dates::date':     [date(1964,1,1), date(1985,2,5), date(2022,1,21), date(1964,1,1), date(1985,2,5), date(2022,1,21)],
                'value':           [10, 10, 20, 20, 30, 30],
                'value32':         pd.Series([12, 12, 22, 22, 32, 32], dtype='int32'),
                'res':             [10, 20, 30, 10, 20, 30],
                'coord::point':    [Point(1,2), Point(3,4), Point(5,6), Point(7,8), Point(3,4), Point(5,6)],
                'names':           pd.Series(['john', 'eric', 'judith', 'mila', 'hector', 'maria'], dtype='string'),
                'unique':          True }

In [3]: df = pd.DataFrame(data).set_index('index')

In [4]: df
Out[4]:       dates::date  value  value32  res coord::point   names  unique
        index
        100    1964-01-01     10       12   10  POINT (1 2)    john    True
        200    1985-02-05     10       12   20  POINT (3 4)    eric    True
        300    2022-01-21     20       22   30  POINT (5 6)  judith    True
        400    1964-01-01     20       22   10  POINT (7 8)    mila    True
        500    1985-02-05     30       32   20  POINT (3 4)  hector    True
        600    2022-01-21     30       32   30  POINT (5 6)   maria    True
```

JSON representation

```python
In [5]: df_to_json = npd.to_json(df)
        pprint(df_to_json, width=120)
Out[5]: {':tab': {'coord::point': [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0], [3.0, 4.0], [5.0, 6.0]],
                  'dates::date': ['1964-01-01', '1985-02-05', '2022-01-21', '1964-01-01', '1985-02-05', '2022-01-21'],
                  'index': [100, 200, 300, 400, 500, 600],
                  'names::string': ['john', 'eric', 'judith', 'mila', 'hector', 'maria'],
                  'res': [10, 20, 30, 10, 20, 30],
                  'unique': [True, True, True, True, True, True],
                  'value': [10, 10, 20, 20, 30, 30],
                  'value32::int32': [12, 12, 22, 22, 32, 32]}}
```

Reversibility

```python
In [5]: df_from_json = npd.read_json(df_to_json)
        print('df created from JSON is equal to initial df ? ', df_from_json.equals(df))
Out[5]: df created from JSON is equal to initial df ?  True
```

Several other examples are provided in the [linked NoteBook](https://nbviewer.org/github/loco-philippe/ntv-pandas/blob/main/example/example_ntv_pandas.ipynb)

#### Table Schema JSON interface example

In the example below, a DataFrame with several Table Schema data types is converted to JSON.

The DataFrame resulting from this JSON is identical to the initial DataFrame (reversibility).

With the existing Table Schema JSON interface, this conversion is not possible.

```python
In [1]: from shapely.geometry import Point
        from datetime import date

In [2]: df = pd.DataFrame({
            'end february::date': ['date(2023,2,28)', 'date(2024,2,29)', 'date(2025,2,28)'],
            'coordinates::point': ['Point([2.3, 48.9])', 'Point([5.4, 43.3])', 'Point([4.9, 45.8])'],
            'contact::email':     ['john.doe@table.com', 'lisa.minelli@schema.com', 'walter.white@breaking.com']
            })

In [3]: df
Out[3]: end february::date coordinates::point             contact::email
        0         2023-02-28   POINT (2.3 48.9)         john.doe@table.com
        1         2024-02-29   POINT (5.4 43.3)    lisa.minelli@schema.com
        2         2025-02-28   POINT (4.9 45.8)  walter.white@breaking.com
```

JSON representation

```python
In [4]: df_to_table = npd.to_json(df, table=True)
        pprint(df_to_table, width=140, sort_dicts=False)
Out[4]: {'schema': {'fields': [{'name': 'index', 'type': 'integer'},
                               {'name': 'end february', 'type': 'date'},
                               {'name': 'coordinates', 'type': 'geopoint', 'format': 'array'},
                               {'name': 'contact', 'type': 'string', 'format': 'email'}],
                    'primaryKey': ['index'],
                    'pandas_version': '1.4.0'},
         'data': [{'index': 0, 'end february': '2023-02-28', 'coordinates': [2.3, 48.9], 'contact': 'john.doe@table.com'},
                  {'index': 1, 'end february': '2024-02-29', 'coordinates': [5.4, 43.3], 'contact': 'lisa.minelli@schema.com'},
                  {'index': 2, 'end february': '2025-02-28', 'coordinates': [4.9, 45.8], 'contact': 'walter.white@breaking.com'}]}
```

Reversibility

```python
In [5]: df_from_table = npd.read_json(df_to_table)
        print('df created from JSON is equal to initial df ? ', df_from_table.equals(df))
Out[5]: df created from JSON is equal to initial df ?  True
```

Several other examples are provided in the [linked NoteBook](https://nbviewer.org/github/loco-philippe/ntv-pandas/blob/main/example/example_table_pandas.ipynb)

## Scope

The objective is to make available the proposed JSON interface for any type of data and for `orient="table"` option or a new option `orient="ntv"`.

The proposed interface is compatible with existing data.

## Motivation

### Why extend the `orient=table` option to other data types?

- The Table Schema specification defines 24 data types, 6 are taken into account in the pandas interface

### Why is it important to have a compact and reversible JSON interface ?

- a reversible interface provides an exchange format.
- a textual exchange format facilitates exchanges between platforms (e.g. OpenData)
- a JSON exchange format can be used at API level

### Is it relevant to take an extended type into account ?

- it avoids the addition of an additional data schema
- it increases the semantic scope of the data processed by pandas
- it is an answer to several issues (e.g.  #12997, #14358, #16492, #35420, #35464, #36211, #39537, #49585, #50782, #51375, #52595, #53252)
- the use of a complementary type avoids having to modify the pandas data model

### Is this only useful for pandas ?

- the JSON-TAB format is applicable to tabular data and multi-dimensional data.
- this JSON interface can therefore be used for any application using tabular or multi-dimensional data. This would allow for example reversible data exchanges between pandas - DataFrame and Xarray - DataArray (Xarray issue under construction) [see example DataFrame / DataArray](https://nbviewer.org/github/loco-philippe/NTV/blob/main/example/example_pandas.ipynb#Multidimensional-data).

## Description

The proposed solution is based on several key points:

- data typing
- correspondence between TableSchema and pandas
- JSON format for tabular data
- conversion to and from JSON format

### Data typing

Data types are defined and managed in the NTV project (name, JSON encoder and decoder).

Pandas `dtype` are compatible with NTV types :

| **pandas  dtype**  | **NTV type**   |
|--------------------|------------|
| intxx              | intxx      |
| uintxx             | uintxx     |
| floatxx            | floatxx    |
| datetime[ns]       | datetime   |
| datetime[ns, <tz>] | datetimetz |
| timedelta[ns]      | durationiso|
| string             | string     |
| boolean            | boolean    |

Note:

- datetime with timezone is a single NTV type (string ISO8601)
- `CategoricalDtype` and `SparseDtype` are included in the tabular JSON format
- `object` `dtype` is depending on the context (see below)
- `PeriodDtype` and `IntervalDtype` are to be defined

JSON types (implicit or explicit) are converted in `dtype` following pandas JSON interface:

| **JSON type**  | **pandas  dtype** |
|----------------|-------------------|
| number         | int64 / float64   |
| string         | string / object   |
| array          | object            |
| object         | object            |
| true, false    | boolean           |
| null           | NaT / NaN / None  |

Note:

- if an NTV type is defined, the `dtype` is adjusted accordingly
- the consideration of null type data needs to be clarified

The other NTV types are associated with `object` `dtype`.

### Correspondence between TableSchema and pandas

The TableSchema typing is carried by two attributes `format` and `type`.

The table below shows the correspondence between TableSchema format / type and pandas NTVtype / dtype:

| **format / type**  | **NTV type / dtype** |
|--------------------|----------------------|
| default / datetime |  / datetime64[ns]    |
| default / number   |  / float64           |
| default / integer  |  / int64             |
| default / boolean  |  / bool              |
| default / string   |  / object            |
| default / duration |  / timedelta64[ns]   |
| email   / string   | email   / string     |
| uri     / string   | uri     / string     |
| default / object   | object  / object     |
| default / array    | array   / object     |
| default / date     | date    / object     |
| default / time     | time    / object     |
| default / year     | year    / int64      |
| default / yearmonth| month   / int64      |
| array   / geopoint | point   / object     |
| default / geojson  | geojson / object     |

Note:

- other TableSchema format are defined and are to be studied (uuid, binary, topojson, specific format for geopoint and datation)
- the first six lines correspond to the existing

### JSON format

The JSON format for the TableSchema interface is the existing.

The JSON format for the Global interface is defined in [JSON-TAB](https://github.com/loco-philippe/NTV/blob/v1.1.0/documentation/JSON-TAB-standard.pdf) specification.
It includes the naming rules originally defined in the [JSON-ND project](https://github.com/glenkleidon/JSON-ND) and support for categorical data.
The specification have to be updated to include sparse data.

### Conversion

When data is associated with a non-`object` `dtype`, pandas conversion methods are used.
Otherwise, NTV conversion is used.

#### pandas -> JSON

- `NTV type` is not defined : use `to_json()`
- `NTV type` is defined and `dtype` is not `object` : use `to_json()`
- `NTV type` is defined and `dtype` is `object` : use NTV conversion (if pandas conversion does not exist)

#### JSON -> pandas

- `NTV type` is compatible with a `dtype` : use `read_json()`
- `NTV type` is not compatible with a `dtype` : use NTV conversion (if pandas conversion does not exist)

## Usage and Impact

### Usage

It seems to me that this proposal responds to important issues:

- having an efficient text format for data exchange

    The alternative CSV format is not reversible and obsolete (last revision in 2005). Current CSV tools do not comply with the standard.

- taking into account "semantic" data in pandas objects

- having a complete Table Schema interface

### Compatibility

Interface can be used without NTV type (compatibility with existing data - [see examples](https://nbviewer.org/github/loco-philippe/ntv-pandas/blob/main/example/example_ntv_pandas.ipynb#Appendix-:-Series-tests))

If the interface is available, throw a new `orient` option in the JSON interface, the use of the feature is decoupled from the other features.

### Impacts on the pandas framework

Initially, the impacts are very limited:

- modification of the `name` of `Series` or `DataFrame columns` (no functional impact),
- added an option in the Json interface (e.g. `orient='ntv'`) and added associated methods (no functional interference with the other methods)

In later stages, several developments could be considered:

- validation of the `name` of `Series` or `DataFrame columns` ,
- management of the NTV type as a "complementary-object-dtype"
- functional extensions depending on the NTV type

### Risk to do / risk not to do

The JSON-NTV format and the JSON-TAB format are not (yet) recognized and used formats. The risk for pandas is that this function is not used (no functional impacts).

On the other hand, the early use by pandas will allow a better consideration of the expectations and needs of pandas as well as a reflection on the evolution of the types supported by pandas.

## Implementation

### Modules

Two modules are defined for NTV:

- json-ntv

    this module manages NTV data without dependency to another module

- ntvconnector

    those modules manage the conversion between objects and JSON data. They have dependency with objects modules (e.g. connectors with shapely location have dependency with shapely).

The pandas integration of the JSON interface requires importing only the json-ntv module.

### Implementation options

The interface can be implemented as NTV connector (`SeriesConnector` and `DataFrameConnector`) and as a new pandas JSON interface `orient` option.

Several pandas implementations are possible:

1. External:

    In this implementation, the interface is available only in the NTV side.
    This option means that this evolution of the JSON interface is not useful or strategic for pandas.

2. NTV side:

    In this implementation, the interface is available in the both sides and the conversion is located inside NTV.
    This option is the one that minimizes the impacts on the pandas side

3. pandas side:

    In this implementation, the interface is available in the both sides and the conversion is located inside pandas.
    This option allows pandas to keep control of this evolution

4. pandas restricted:

    In this implementation, the pandas interface and the conversion are located inside pandas and only for non-object `dtype`.
    This option makes it possible to offer a compact and reversible interface while prohibiting the introduction of types incompatible with the existing `dtype`

## F.A.Q.

**Q: Does `orient="table"` not do what you are proposing already?**

**A**: In principle, yes, this option takes into account the notion of type.

But this is very limited (see examples added in the [Notebook](https://nbviewer.org/github/loco-philippe/NTV/blob/main/example/example_pandas.ipynb)) :

- **Types and Json interface**
  - the only way to keep the types in the json interface is to use the `orient='table'` option
  - few dtypes are not allowed in json-table interface : period, timedelta64, interval
  - allowed types are not always kept in json-table interface
  - data with 'object' dtype is kept only id data is string
  - with categorical dtype, the underlying dtype is not included in json interface
- **Data compactness**
  - json-table interface is not compact (in the example in the [Notebook](https://nbviewer.org/github/loco-philippe/NTV/blob/main/example/example_pandas.ipynb#data-compactness))the size is triple or quadruple the size of the compact format
- **Reversibility**
  - Interface is reversible only with few dtypes : int64, float64, bool, string, datetime64 and partially categorical
- **External types**
  - the interface does not accept external types
  - Table-schema defines 20 data types but the `orient="table"` interface takes into account 5 data types (see [table](https://nbviewer.org/github/loco-philippe/NTV/blob/main/example/example_pandas.ipynb#Converting-table-schema-type-to-pandas-dtype))
  - to integrate external types, it is necessary to first create ExtensionArray and ExtensionDtype objects

The current interface is not compatible with the data structure defined by table-schema. For this to be possible, it is necessary to integrate a "type extension" like the one proposed (this has moreover been partially achieved with the notion of `extDtype` found in the interface for several formats).

**Q: In general, we should only have 1 `"table"` format for pandas in read_json/to_json. There is also the issue of backwards compatibility if we do change the format. The fact that the table interface is buggy is not a reason to add a new interface (I'd rather fix those bugs). Can the existing format be adapted in a way that fixes the type issues/issues with roundtripping?**

**A**: I will add two additional remarks:

- the types defined in Tableschema are partially taken into account (examples of types not taken into account in the interface: string-uri, array, date, time, year, geopoint, string-email):
- the `read_json()` interface works too with the following data: `{'simple': [1,2,3] }` (contrary to what is indicated in the documentation) but it is impossible with `to_json()` to recreate this simple json.

I think that the problem cannot be limited to bug fixes and that a clear strategy must be defined for the Json interface in particular with the gradual abandonment in open-data solutions of the obsolete CSV format in favor of a Json format.

As stated, the proposed solution addresses several shortcomings of the current interface and could simply fit into the pandas environment (the other option would be to consider that the Json interface is a peripheral function of pandas and can remain external to pandas) regardless of the `orient='table'` option.

It is nevertheless possible to merge the proposed format and the `orient='table'` format in order to have an explicit management of the notion of `extDtype`

**Q: As far as I can tell, JSON NTV is not in any form a standardised JSON format. I believe that pandas (and geopandas, which is where I came from to this issue) should try to follow either de facto or de jure standards and do not opt in for a file format that does not have any community support at this moment. This can obviously change in the future and that is where this PR should be revised. Why would pandas use this standard?**

**A**: As indicated in the issue (and detailed in [the attached Notebook](https://nbviewer.org/github/loco-philippe/NTV/blob/main/example/example_pandas.ipynb)), the json interface is not reversible (`to_json` then `read_json` does not always return the initial object) and several shortcomings and bugs are present. The main cause of this problem is that the data type is not taken into account in the JSON format (or very partially with the `orient='table'` option).

The proposal made answers this problem ([the example at the beginning of Notebook](https://nbviewer.org/github/loco-philippe/NTV/blob/main/example/example_pandas.ipynb#0---Simple-example) simply and clearly illustrates the interest of the proposal).

Regarding the underlying JSON-NTV format, its impact is quite low for tabular data (it is limited to adding the type in the field name).
Nevertheless, the question is relevant: The JSON-NTV format ([IETF Internet-Draft](https://datatracker.ietf.org/doc/draft-thomy-json-ntv/)) is a shared, documented, supported and implemented format, but indeed the community support is for the moment reduced but it only asks to expand !!

## Synthesis

To conclude,

- if it is important (or strategic) to have a reversible JSON interface for any type of data, the proposal can be allowed,
- if not, a third-party package listed in the [ecosystem](https://pandas.pydata.org/community/ecosystem.html) that reads/writes this format to/from pandas DataFrames should be considered

## Core team decision

Vote was open from september-11 to setpember-26:

- Final tally is 0 approvals, 5 abstentions, 7 disapprove. The quorum has been met. The PDEP fails.

**Disapprove comments** :

- 1 Given the newness of the proposed JSON NTV format, I would support (as described in the PDEP): "if not, a third-party package listed in the ecosystem that reads/writes this format to/from pandas DataFrames should be considered"
- 2 Same reason as -1-, this should be a third party package for now
- 3 Not mature enough, and not clear what the market size would be.
- 4 for the same reason I left in the PDEP: "I think this (JSON-NTV format) does not meet the bar of being a commonly used format for implementation within pandas"
- 5 agree with -4-
- 6 agree with the other core-dev responders. I think work in the existing json interface is extremely valuable. A number of the original issues raised are just bug fixes / extensions of already existing functionality. Trying to start anew is likely not worth the migration effort. That said if a format is well supported in the community we can reconsider in the future (obviously json is well supported but the actual specification detailed here is too new / not accepted as a standard)
- 7 while I do think having a more comprehensive JSON format would be worthwhile, making a new format part of pandas means an implicit endorsement of a standard that is still being reviewed by the broader community.

**Decision**:

- add the `ntv-pandas` package in the [ecosystem](https://pandas.pydata.org/community/ecosystem.html)
- revisit again this PDEP at a later stage, for example in 1/2 to 1 year (based on the evolution of the Internet draft [JSON semantic format (JSON-NTV)](https://www.ietf.org/archive/id/draft-thomy-json-ntv-01.html) and the usage of the [ntv-pandas](https://github.com/loco-philippe/ntv-pandas#readme))

## Timeline

Not applicable

## PDEP History

- 16 June 2023: Initial draft
- 22 July 2023: Add F.A.Q.
- 06 September 2023: Add Table Schema extension
- 01 Octobre: Add Core team decision

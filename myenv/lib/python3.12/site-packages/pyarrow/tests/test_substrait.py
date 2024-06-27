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

import os
import pathlib

import pytest

import pyarrow as pa
import pyarrow.compute as pc
from pyarrow.lib import tobytes
from pyarrow.lib import ArrowInvalid, ArrowNotImplementedError

try:
    import pyarrow.substrait as substrait
except ImportError:
    substrait = None

# Marks all of the tests in this module
# Ignore these with pytest ... -m 'not substrait'
pytestmark = pytest.mark.substrait


def mock_udf_context(batch_length=10):
    from pyarrow._compute import _get_udf_context
    return _get_udf_context(pa.default_memory_pool(), batch_length)


def _write_dummy_data_to_disk(tmpdir, file_name, table):
    path = os.path.join(str(tmpdir), file_name)
    with pa.ipc.RecordBatchFileWriter(path, schema=table.schema) as writer:
        writer.write_table(table)
    return path


@pytest.mark.parametrize("use_threads", [True, False])
def test_run_serialized_query(tmpdir, use_threads):
    substrait_query = """
    {
        "version": { "major": 9999 },
        "relations": [
        {"rel": {
            "read": {
            "base_schema": {
                "struct": {
                "types": [
                            {"i64": {}}
                        ]
                },
                "names": [
                        "foo"
                        ]
            },
            "local_files": {
                "items": [
                {
                    "uri_file": "FILENAME_PLACEHOLDER",
                    "arrow": {}
                }
                ]
            }
            }
        }}
        ]
    }
    """

    file_name = "read_data.arrow"
    table = pa.table([[1, 2, 3, 4, 5]], names=['foo'])
    path = _write_dummy_data_to_disk(tmpdir, file_name, table)
    query = tobytes(substrait_query.replace(
        "FILENAME_PLACEHOLDER", pathlib.Path(path).as_uri()))

    buf = pa._substrait._parse_json_plan(query)

    reader = substrait.run_query(buf, use_threads=use_threads)
    res_tb = reader.read_all()

    assert table.select(["foo"]) == res_tb.select(["foo"])


@pytest.mark.parametrize("query", (pa.py_buffer(b'buffer'), b"bytes", 1))
def test_run_query_input_types(tmpdir, query):

    # Passing unsupported type, like int, will not segfault.
    if not isinstance(query, (pa.Buffer, bytes)):
        msg = f"Expected 'pyarrow.Buffer' or bytes, got '{type(query)}'"
        with pytest.raises(TypeError, match=msg):
            substrait.run_query(query)
        return

    # Otherwise error for invalid query
    msg = "ParseFromZeroCopyStream failed for substrait.Plan"
    with pytest.raises(OSError, match=msg):
        substrait.run_query(query)


def test_invalid_plan():
    query = """
    {
        "relations": [
        ]
    }
    """
    buf = pa._substrait._parse_json_plan(tobytes(query))
    exec_message = "Plan has no relations"
    with pytest.raises(ArrowInvalid, match=exec_message):
        substrait.run_query(buf)


@pytest.mark.parametrize("use_threads", [True, False])
def test_binary_conversion_with_json_options(tmpdir, use_threads):
    substrait_query = """
    {
        "version": { "major": 9999 },
        "relations": [
        {"rel": {
            "read": {
            "base_schema": {
                "struct": {
                "types": [
                            {"i64": {}}
                        ]
                },
                "names": [
                        "bar"
                        ]
            },
            "local_files": {
                "items": [
                {
                    "uri_file": "FILENAME_PLACEHOLDER",
                    "arrow": {},
                    "metadata" : {
                      "created_by" : {},
                    }
                }
                ]
            }
            }
        }}
        ]
    }
    """

    file_name = "binary_json_data.arrow"
    table = pa.table([[1, 2, 3, 4, 5]], names=['bar'])
    path = _write_dummy_data_to_disk(tmpdir, file_name, table)
    query = tobytes(substrait_query.replace(
        "FILENAME_PLACEHOLDER", pathlib.Path(path).as_uri()))
    buf = pa._substrait._parse_json_plan(tobytes(query))

    reader = substrait.run_query(buf, use_threads=use_threads)
    res_tb = reader.read_all()

    assert table.select(["bar"]) == res_tb.select(["bar"])


# Substrait has not finalized what the URI should be for standard functions
# In the meantime, lets just check the suffix
def has_function(fns, ext_file, fn_name):
    suffix = f'{ext_file}#{fn_name}'
    for fn in fns:
        if fn.endswith(suffix):
            return True
    return False


def test_get_supported_functions():
    supported_functions = pa._substrait.get_supported_functions()
    # It probably doesn't make sense to exhaustively verify this list but
    # we can check a sample aggregate and a sample non-aggregate entry
    assert has_function(supported_functions,
                        'functions_arithmetic.yaml', 'add')
    assert has_function(supported_functions,
                        'functions_arithmetic.yaml', 'sum')


@pytest.mark.parametrize("use_threads", [True, False])
def test_named_table(use_threads):
    test_table_1 = pa.Table.from_pydict({"x": [1, 2, 3]})
    test_table_2 = pa.Table.from_pydict({"x": [4, 5, 6]})
    schema_1 = pa.schema([pa.field("x", pa.int64())])

    def table_provider(names, schema):
        if not names:
            raise Exception("No names provided")
        elif names[0] == "t1":
            assert schema == schema_1
            return test_table_1
        elif names[1] == "t2":
            return test_table_2
        else:
            raise Exception("Unrecognized table name")

    substrait_query = """
    {
        "version": { "major": 9999 },
        "relations": [
        {"rel": {
            "read": {
            "base_schema": {
                "struct": {
                "types": [
                            {"i64": {}}
                        ]
                },
                "names": [
                        "x"
                        ]
            },
            "namedTable": {
                    "names": ["t1"]
            }
            }
        }}
        ]
    }
    """

    buf = pa._substrait._parse_json_plan(tobytes(substrait_query))
    reader = pa.substrait.run_query(
        buf, table_provider=table_provider, use_threads=use_threads)
    res_tb = reader.read_all()
    assert res_tb == test_table_1


def test_named_table_invalid_table_name():
    test_table_1 = pa.Table.from_pydict({"x": [1, 2, 3]})

    def table_provider(names, _):
        if not names:
            raise Exception("No names provided")
        elif names[0] == "t1":
            return test_table_1
        else:
            raise Exception("Unrecognized table name")

    substrait_query = """
    {
        "version": { "major": 9999 },
        "relations": [
        {"rel": {
            "read": {
            "base_schema": {
                "struct": {
                "types": [
                            {"i64": {}}
                        ]
                },
                "names": [
                        "x"
                        ]
            },
            "namedTable": {
                    "names": ["t3"]
            }
            }
        }}
        ]
    }
    """

    buf = pa._substrait._parse_json_plan(tobytes(substrait_query))
    exec_message = "Invalid NamedTable Source"
    with pytest.raises(ArrowInvalid, match=exec_message):
        substrait.run_query(buf, table_provider=table_provider)


def test_named_table_empty_names():
    test_table_1 = pa.Table.from_pydict({"x": [1, 2, 3]})

    def table_provider(names, _):
        if not names:
            raise Exception("No names provided")
        elif names[0] == "t1":
            return test_table_1
        else:
            raise Exception("Unrecognized table name")

    substrait_query = """
    {
        "version": { "major": 9999 },
        "relations": [
        {"rel": {
            "read": {
            "base_schema": {
                "struct": {
                "types": [
                            {"i64": {}}
                        ]
                },
                "names": [
                        "x"
                        ]
            },
            "namedTable": {
                    "names": []
            }
            }
        }}
        ]
    }
    """
    query = tobytes(substrait_query)
    buf = pa._substrait._parse_json_plan(tobytes(query))
    exec_message = "names for NamedTable not provided"
    with pytest.raises(ArrowInvalid, match=exec_message):
        substrait.run_query(buf, table_provider=table_provider)


@pytest.mark.parametrize("use_threads", [True, False])
def test_udf_via_substrait(unary_func_fixture, use_threads):
    test_table = pa.Table.from_pydict({"x": [1, 2, 3]})

    def table_provider(names, _):
        if not names:
            raise Exception("No names provided")
        elif names[0] == "t1":
            return test_table
        else:
            raise Exception("Unrecognized table name")

    substrait_query = b"""
    {
  "extensionUris": [
    {
      "extensionUriAnchor": 1
    },
    {
      "extensionUriAnchor": 2,
      "uri": "urn:arrow:substrait_simple_extension_function"
    }
  ],
  "extensions": [
    {
      "extensionFunction": {
        "extensionUriReference": 2,
        "functionAnchor": 1,
        "name": "y=x+1"
      }
    }
  ],
  "relations": [
    {
      "root": {
        "input": {
          "project": {
            "common": {
              "emit": {
                "outputMapping": [
                  1,
                  2,
                ]
              }
            },
            "input": {
              "read": {
                "baseSchema": {
                  "names": [
                    "t",
                  ],
                  "struct": {
                    "types": [
                      {
                        "i64": {
                          "nullability": "NULLABILITY_REQUIRED"
                        }
                      },
                    ],
                    "nullability": "NULLABILITY_REQUIRED"
                  }
                },
                "namedTable": {
                  "names": [
                    "t1"
                  ]
                }
              }
            },
            "expressions": [
              {
                "selection": {
                  "directReference": {
                    "structField": {}
                  },
                  "rootReference": {}
                }
              },
              {
                "scalarFunction": {
                  "functionReference": 1,
                  "outputType": {
                    "i64": {
                      "nullability": "NULLABILITY_NULLABLE"
                    }
                  },
                  "arguments": [
                    {
                      "value": {
                        "selection": {
                          "directReference": {
                            "structField": {}
                          },
                          "rootReference": {}
                        }
                      }
                    }
                  ]
                }
              }
            ]
          }
        },
        "names": [
          "x",
          "y",
        ]
      }
    }
  ]
}
    """

    buf = pa._substrait._parse_json_plan(substrait_query)
    reader = pa.substrait.run_query(
        buf, table_provider=table_provider, use_threads=use_threads)
    res_tb = reader.read_all()

    function, name = unary_func_fixture
    expected_tb = test_table.add_column(1, 'y', function(
        mock_udf_context(10), test_table['x']))
    assert res_tb == expected_tb


def test_udf_via_substrait_wrong_udf_name():
    test_table = pa.Table.from_pydict({"x": [1, 2, 3]})

    def table_provider(names, _):
        if not names:
            raise Exception("No names provided")
        elif names[0] == "t1":
            return test_table
        else:
            raise Exception("Unrecognized table name")

    substrait_query = b"""
    {
  "extensionUris": [
    {
      "extensionUriAnchor": 1
    },
    {
      "extensionUriAnchor": 2,
      "uri": "urn:arrow:substrait_simple_extension_function"
    }
  ],
  "extensions": [
    {
      "extensionFunction": {
        "extensionUriReference": 2,
        "functionAnchor": 1,
        "name": "wrong_udf_name"
      }
    }
  ],
  "relations": [
    {
      "root": {
        "input": {
          "project": {
            "common": {
              "emit": {
                "outputMapping": [
                  1,
                  2,
                ]
              }
            },
            "input": {
              "read": {
                "baseSchema": {
                  "names": [
                    "t",
                  ],
                  "struct": {
                    "types": [
                      {
                        "i64": {
                          "nullability": "NULLABILITY_REQUIRED"
                        }
                      },
                    ],
                    "nullability": "NULLABILITY_REQUIRED"
                  }
                },
                "namedTable": {
                  "names": [
                    "t1"
                  ]
                }
              }
            },
            "expressions": [
              {
                "selection": {
                  "directReference": {
                    "structField": {}
                  },
                  "rootReference": {}
                }
              },
              {
                "scalarFunction": {
                  "functionReference": 1,
                  "outputType": {
                    "i64": {
                      "nullability": "NULLABILITY_NULLABLE"
                    }
                  },
                  "arguments": [
                    {
                      "value": {
                        "selection": {
                          "directReference": {
                            "structField": {}
                          },
                          "rootReference": {}
                        }
                      }
                    }
                  ]
                }
              }
            ]
          }
        },
        "names": [
          "x",
          "y",
        ]
      }
    }
  ]
}
    """

    buf = pa._substrait._parse_json_plan(substrait_query)
    with pytest.raises(pa.ArrowKeyError) as excinfo:
        pa.substrait.run_query(buf, table_provider=table_provider)
    assert "No function registered" in str(excinfo.value)


@pytest.mark.parametrize("use_threads", [True, False])
def test_output_field_names(use_threads):
    in_table = pa.Table.from_pydict({"x": [1, 2, 3]})

    def table_provider(names, schema):
        return in_table

    substrait_query = """
    {
      "version": { "major": 9999 },
      "relations": [
        {
          "root": {
            "input": {
              "read": {
                "base_schema": {
                  "struct": {
                    "types": [{"i64": {}}]
                  },
                  "names": ["x"]
                },
                "namedTable": {
                  "names": ["t1"]
                }
              }
            },
            "names": ["out"]
          }
        }
      ]
    }
    """

    buf = pa._substrait._parse_json_plan(tobytes(substrait_query))
    reader = pa.substrait.run_query(
        buf, table_provider=table_provider, use_threads=use_threads)
    res_tb = reader.read_all()

    expected = pa.Table.from_pydict({"out": [1, 2, 3]})

    assert res_tb == expected


def test_scalar_aggregate_udf_basic(varargs_agg_func_fixture):

    test_table = pa.Table.from_pydict(
        {"k": [1, 1, 2, 2], "v1": [1, 2, 3, 4],
         "v2": [1.0, 1.0, 1.0, 1.0]}
    )

    def table_provider(names, _):
        return test_table

    substrait_query = b"""
{
  "extensionUris": [
    {
      "extensionUriAnchor": 1,
      "uri": "urn:arrow:substrait_simple_extension_function"
    },
  ],
  "extensions": [
    {
      "extensionFunction": {
        "extensionUriReference": 1,
        "functionAnchor": 1,
        "name": "sum_mean"
      }
    }
  ],
  "relations": [
    {
      "root": {
        "input": {
          "extensionSingle": {
            "common": {
              "emit": {
                "outputMapping": [
                  0,
                  1
                ]
              }
            },
            "input": {
              "read": {
                "baseSchema": {
                  "names": [
                    "k",
                    "v1",
                    "v2",
                  ],
                  "struct": {
                    "types": [
                      {
                        "i64": {
                          "nullability": "NULLABILITY_REQUIRED"
                        }
                      },
                      {
                        "i64": {
                          "nullability": "NULLABILITY_NULLABLE"
                        }
                      },
                      {
                        "fp64": {
                          "nullability": "NULLABILITY_NULLABLE"
                        }
                      }
                    ],
                    "nullability": "NULLABILITY_REQUIRED"
                  }
                },
                "namedTable": {
                  "names": ["t1"]
                }
              }
            },
            "detail": {
              "@type": "/arrow.substrait_ext.SegmentedAggregateRel",
              "segmentKeys": [
                {
                  "directReference": {
                    "structField": {}
                  },
                  "rootReference": {}
                }
              ],
              "measures": [
                {
                  "measure": {
                    "functionReference": 1,
                    "phase": "AGGREGATION_PHASE_INITIAL_TO_RESULT",
                    "outputType": {
                      "fp64": {
                        "nullability": "NULLABILITY_NULLABLE"
                      }
                    },
                    "arguments": [
                      {
                        "value": {
                          "selection": {
                            "directReference": {
                              "structField": {
                                "field": 1
                              }
                            },
                            "rootReference": {}
                          }
                        }
                      },
                      {
                        "value": {
                          "selection": {
                            "directReference": {
                              "structField": {
                                "field": 2
                              }
                            },
                            "rootReference": {}
                          }
                        }
                      }
                    ]
                  }
                }
              ]
            }
          }
        },
        "names": [
          "k",
          "v_avg"
        ]
      }
    }
  ],
}
"""
    buf = pa._substrait._parse_json_plan(substrait_query)
    reader = pa.substrait.run_query(
        buf, table_provider=table_provider, use_threads=False)
    res_tb = reader.read_all()

    expected_tb = pa.Table.from_pydict({
        'k': [1, 2],
        'v_avg': [2.5, 4.5]
    })

    assert res_tb == expected_tb


def test_hash_aggregate_udf_basic(varargs_agg_func_fixture):

    test_table = pa.Table.from_pydict(
        {"t": [1, 1, 1, 1, 2, 2, 2, 2],
         "k": [1, 0, 0, 1, 0, 1, 0, 1],
         "v1": [1, 2, 3, 4, 5, 6, 7, 8],
         "v2": [1.0, 1.0, 1.0, 1.0, 2.0, 3.0, 4.0, 5.0]}
    )

    def table_provider(names, _):
        return test_table

    substrait_query = b"""
{
  "extensionUris": [
    {
      "extensionUriAnchor": 1,
      "uri": "urn:arrow:substrait_simple_extension_function"
    },
  ],
  "extensions": [
    {
      "extensionFunction": {
        "extensionUriReference": 1,
        "functionAnchor": 1,
        "name": "sum_mean"
      }
    }
  ],
  "relations": [
    {
      "root": {
        "input": {
          "extensionSingle": {
            "common": {
              "emit": {
                "outputMapping": [
                  0,
                  1,
                  2
                ]
              }
            },
            "input": {
              "read": {
                "baseSchema": {
                  "names": [
                    "t",
                    "k",
                    "v1",
                    "v2",
                  ],
                  "struct": {
                    "types": [
                      {
                        "i64": {
                          "nullability": "NULLABILITY_REQUIRED"
                        }
                      },
                      {
                        "i64": {
                          "nullability": "NULLABILITY_REQUIRED"
                        }
                      },
                      {
                        "i64": {
                          "nullability": "NULLABILITY_NULLABLE"
                        }
                      },
                      {
                        "fp64": {
                          "nullability": "NULLABILITY_NULLABLE"
                        }
                      }
                    ],
                    "nullability": "NULLABILITY_REQUIRED"
                  }
                },
                "namedTable": {
                  "names": ["t1"]
                }
              }
            },
            "detail": {
              "@type": "/arrow.substrait_ext.SegmentedAggregateRel",
              "groupingKeys": [
                {
                  "directReference": {
                    "structField": {
                      "field": 1
                    }
                  },
                  "rootReference": {}
                }
              ],
              "segmentKeys": [
                {
                  "directReference": {
                    "structField": {}
                  },
                  "rootReference": {}
                }
              ],
              "measures": [
                {
                  "measure": {
                    "functionReference": 1,
                    "phase": "AGGREGATION_PHASE_INITIAL_TO_RESULT",
                    "outputType": {
                      "fp64": {
                        "nullability": "NULLABILITY_NULLABLE"
                      }
                    },
                    "arguments": [
                      {
                        "value": {
                          "selection": {
                            "directReference": {
                              "structField": {
                                "field": 2
                              }
                            },
                            "rootReference": {}
                          }
                        }
                      },
                      {
                        "value": {
                          "selection": {
                            "directReference": {
                              "structField": {
                                "field": 3
                              }
                            },
                            "rootReference": {}
                          }
                        }
                      }
                    ]
                  }
                }
              ]
            }
          }
        },
        "names": [
          "t",
          "k",
          "v_avg"
        ]
      }
    }
  ],
}
"""
    buf = pa._substrait._parse_json_plan(substrait_query)
    reader = pa.substrait.run_query(
        buf, table_provider=table_provider, use_threads=False)
    res_tb = reader.read_all()

    expected_tb = pa.Table.from_pydict({
        't': [1, 1, 2, 2],
        'k': [1, 0, 0, 1],
        'v_avg': [3.5, 3.5, 9.0, 11.0]
    })

    # Ordering of k is deterministic because this is running with serial execution
    assert res_tb == expected_tb


@pytest.mark.parametrize("expr", [
    pc.equal(pc.field("x"), 7),
    pc.equal(pc.field("x"), pc.field("y")),
    pc.field("x") > 50
])
def test_serializing_expressions(expr):
    schema = pa.schema([
        pa.field("x", pa.int32()),
        pa.field("y", pa.int32())
    ])

    buf = pa.substrait.serialize_expressions([expr], ["test_expr"], schema)
    returned = pa.substrait.deserialize_expressions(buf)
    assert schema == returned.schema
    assert len(returned.expressions) == 1
    assert "test_expr" in returned.expressions


def test_arrow_specific_types():
    fields = {
        "time_seconds": (pa.time32("s"), 0),
        "time_millis": (pa.time32("ms"), 0),
        "time_nanos": (pa.time64("ns"), 0),
        "date_millis": (pa.date64(), 0),
        "large_string": (pa.large_string(), "test_string"),
        "large_binary": (pa.large_binary(), b"test_string"),
    }
    schema = pa.schema([pa.field(name, typ) for name, (typ, _) in fields.items()])

    def check_round_trip(expr):
        buf = pa.substrait.serialize_expressions([expr], ["test_expr"], schema)
        returned = pa.substrait.deserialize_expressions(buf)
        assert schema == returned.schema

    for name, (typ, val) in fields.items():
        check_round_trip(pc.field(name) == pa.scalar(val, type=typ))


def test_arrow_one_way_types():
    schema = pa.schema(
        [
            pa.field("binary_view", pa.binary_view()),
            pa.field("string_view", pa.string_view()),
            pa.field("dictionary", pa.dictionary(pa.int32(), pa.string())),
            pa.field("ree", pa.run_end_encoded(pa.int32(), pa.string())),
        ]
    )
    alt_schema = pa.schema(
        [
            pa.field("binary_view", pa.binary()),
            pa.field("string_view", pa.string()),
            pa.field("dictionary", pa.string()),
            pa.field("ree", pa.string())
        ]
    )

    def check_one_way(field):
        expr = pc.is_null(pc.field(field.name))
        buf = pa.substrait.serialize_expressions([expr], ["test_expr"], schema)
        returned = pa.substrait.deserialize_expressions(buf)
        assert alt_schema == returned.schema

    for field in schema:
        check_one_way(field)


def test_invalid_expression_ser_des():
    schema = pa.schema([
        pa.field("x", pa.int32()),
        pa.field("y", pa.int32())
    ])
    expr = pc.equal(pc.field("x"), 7)
    bad_expr = pc.equal(pc.field("z"), 7)
    # Invalid number of names
    with pytest.raises(ValueError) as excinfo:
        pa.substrait.serialize_expressions([expr], [], schema)
    assert 'need to have the same length' in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        pa.substrait.serialize_expressions([expr], ["foo", "bar"], schema)
    assert 'need to have the same length' in str(excinfo.value)
    # Expression doesn't match schema
    with pytest.raises(ValueError) as excinfo:
        pa.substrait.serialize_expressions([bad_expr], ["expr"], schema)
    assert 'No match for FieldRef' in str(excinfo.value)


def test_serializing_multiple_expressions():
    schema = pa.schema([
        pa.field("x", pa.int32()),
        pa.field("y", pa.int32())
    ])
    exprs = [pc.equal(pc.field("x"), 7), pc.equal(pc.field("x"), pc.field("y"))]
    buf = pa.substrait.serialize_expressions(exprs, ["first", "second"], schema)
    returned = pa.substrait.deserialize_expressions(buf)
    assert schema == returned.schema
    assert len(returned.expressions) == 2

    norm_exprs = [pc.equal(pc.field(0), 7), pc.equal(pc.field(0), pc.field(1))]
    assert str(returned.expressions["first"]) == str(norm_exprs[0])
    assert str(returned.expressions["second"]) == str(norm_exprs[1])


def test_serializing_with_compute():
    schema = pa.schema([
        pa.field("x", pa.int32()),
        pa.field("y", pa.int32())
    ])
    expr = pc.equal(pc.field("x"), 7)
    expr_norm = pc.equal(pc.field(0), 7)
    buf = expr.to_substrait(schema)
    returned = pa.substrait.deserialize_expressions(buf)

    assert schema == returned.schema
    assert len(returned.expressions) == 1

    assert str(returned.expressions["expression"]) == str(expr_norm)

    # Compute can't deserialize messages with multiple expressions
    buf = pa.substrait.serialize_expressions([expr, expr], ["first", "second"], schema)
    with pytest.raises(ValueError) as excinfo:
        pc.Expression.from_substrait(buf)
    assert 'contained multiple expressions' in str(excinfo.value)

    # Deserialization should be possible regardless of the expression name
    buf = pa.substrait.serialize_expressions([expr], ["weirdname"], schema)
    expr2 = pc.Expression.from_substrait(buf)
    assert str(expr2) == str(expr_norm)


def test_serializing_udfs():
    # Note, UDF in this context means a function that is not
    # recognized by Substrait.  It might still be a builtin pyarrow
    # function.
    schema = pa.schema([
        pa.field("x", pa.uint32())
    ])
    a = pc.scalar(10)
    b = pc.scalar(4)
    exprs = [pc.shift_left(a, b)]

    with pytest.raises(ArrowNotImplementedError):
        pa.substrait.serialize_expressions(exprs, ["expr"], schema)

    buf = pa.substrait.serialize_expressions(
        exprs, ["expr"], schema, allow_arrow_extensions=True)
    returned = pa.substrait.deserialize_expressions(buf)
    assert schema == returned.schema
    assert len(returned.expressions) == 1
    assert str(returned.expressions["expr"]) == str(exprs[0])

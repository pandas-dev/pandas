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

# distutils: language = c++

from pyarrow.includes.common cimport *
from pyarrow.includes.libarrow cimport *


cdef extern from "arrow/acero/options.h" namespace "arrow::acero" nogil:
    cdef enum CJoinType "arrow::acero::JoinType":
        CJoinType_LEFT_SEMI "arrow::acero::JoinType::LEFT_SEMI"
        CJoinType_RIGHT_SEMI "arrow::acero::JoinType::RIGHT_SEMI"
        CJoinType_LEFT_ANTI "arrow::acero::JoinType::LEFT_ANTI"
        CJoinType_RIGHT_ANTI "arrow::acero::JoinType::RIGHT_ANTI"
        CJoinType_INNER "arrow::acero::JoinType::INNER"
        CJoinType_LEFT_OUTER "arrow::acero::JoinType::LEFT_OUTER"
        CJoinType_RIGHT_OUTER "arrow::acero::JoinType::RIGHT_OUTER"
        CJoinType_FULL_OUTER "arrow::acero::JoinType::FULL_OUTER"

    cdef cppclass CExecNodeOptions "arrow::acero::ExecNodeOptions":
        pass

    cdef cppclass CSourceNodeOptions "arrow::acero::SourceNodeOptions"(CExecNodeOptions):
        pass

    cdef cppclass CTableSourceNodeOptions "arrow::acero::TableSourceNodeOptions"(CExecNodeOptions):
        CTableSourceNodeOptions(shared_ptr[CTable] table)
        CTableSourceNodeOptions(shared_ptr[CTable] table, int64_t max_batch_size)

    cdef cppclass CSinkNodeOptions "arrow::acero::SinkNodeOptions"(CExecNodeOptions):
        pass

    cdef cppclass CFilterNodeOptions "arrow::acero::FilterNodeOptions"(CExecNodeOptions):
        CFilterNodeOptions(CExpression)

    cdef cppclass CProjectNodeOptions "arrow::acero::ProjectNodeOptions"(CExecNodeOptions):
        CProjectNodeOptions(vector[CExpression] expressions)
        CProjectNodeOptions(vector[CExpression] expressions,
                            vector[c_string] names)

    cdef cppclass CAggregateNodeOptions "arrow::acero::AggregateNodeOptions"(CExecNodeOptions):
        CAggregateNodeOptions(vector[CAggregate] aggregates, vector[CFieldRef] names)

    cdef cppclass COrderByNodeOptions "arrow::acero::OrderByNodeOptions"(CExecNodeOptions):
        COrderByNodeOptions(COrdering ordering)

    cdef cppclass CHashJoinNodeOptions "arrow::acero::HashJoinNodeOptions"(CExecNodeOptions):
        CHashJoinNodeOptions(CJoinType, vector[CFieldRef] in_left_keys,
                             vector[CFieldRef] in_right_keys)
        CHashJoinNodeOptions(CJoinType, vector[CFieldRef] in_left_keys,
                             vector[CFieldRef] in_right_keys,
                             CExpression filter,
                             c_string output_suffix_for_left,
                             c_string output_suffix_for_right)
        CHashJoinNodeOptions(CJoinType join_type,
                             vector[CFieldRef] left_keys,
                             vector[CFieldRef] right_keys,
                             vector[CFieldRef] left_output,
                             vector[CFieldRef] right_output,
                             CExpression filter,
                             c_string output_suffix_for_left,
                             c_string output_suffix_for_right)

    cdef struct CAsofJoinKeys "arrow::acero::AsofJoinNodeOptions::Keys":
        CFieldRef on_key
        vector[CFieldRef] by_key

    cdef cppclass CAsofJoinNodeOptions "arrow::acero::AsofJoinNodeOptions"(CExecNodeOptions):
        CAsofJoinNodeOptions(vector[CAsofJoinKeys] keys, int64_t tolerance)


cdef extern from "arrow/acero/exec_plan.h" namespace "arrow::acero" nogil:
    cdef cppclass CDeclaration "arrow::acero::Declaration":
        cppclass Input:
            Input(CExecNode*)
            Input(CDeclaration)

        c_string label
        vector[Input] inputs

        CDeclaration()
        CDeclaration(c_string factory_name, CExecNodeOptions options)
        CDeclaration(c_string factory_name, vector[Input] inputs, shared_ptr[CExecNodeOptions] options)

        @staticmethod
        CDeclaration Sequence(vector[CDeclaration] decls)

    cdef cppclass CExecNode "arrow::acero::ExecNode":
        const vector[CExecNode*]& inputs() const
        const shared_ptr[CSchema]& output_schema() const

    CResult[shared_ptr[CTable]] DeclarationToTable(
        CDeclaration declaration, c_bool use_threads
    )
    CResult[shared_ptr[CTable]] DeclarationToTable(
        CDeclaration declaration, c_bool use_threads,
        CMemoryPool* memory_pool, CFunctionRegistry* function_registry
    )
    CResult[unique_ptr[CRecordBatchReader]] DeclarationToReader(
        CDeclaration declaration, c_bool use_threads
    )

    CResult[c_string] DeclarationToString(const CDeclaration& declaration)

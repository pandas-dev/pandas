// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#pragma once

#include <memory>

#include "arrow/acero/exec_plan.h"
#include "arrow/compute/api_aggregate.h"
#include "arrow/engine/substrait/visibility.h"
#include "arrow/type_fwd.h"

namespace arrow {
namespace engine {

/// Execution information resulting from converting a Substrait relation.
struct ARROW_ENGINE_EXPORT DeclarationInfo {
  /// The compute declaration produced thus far.
  acero::Declaration declaration;

  std::shared_ptr<Schema> output_schema;
};

/// Information resulting from converting a Substrait plan
struct ARROW_ENGINE_EXPORT PlanInfo {
  /// The root declaration.
  ///
  /// Only plans containing a single top-level relation are supported and so this will
  /// represent that relation.
  ///
  /// This should technically be a RelRoot but some producers use a simple Rel here and so
  /// Acero currently supports that case.
  DeclarationInfo root;
  /// The names of the output fields
  ///
  /// If `root` was created from a simple Rel then this will be empty
  std::vector<std::string> names;
};

/// An expression whose output has a name
struct ARROW_ENGINE_EXPORT NamedExpression {
  /// An expression
  compute::Expression expression;
  // An optional name to assign to the output, may be the empty string
  std::string name;
};

/// A collection of expressions bound to a common schema
struct ARROW_ENGINE_EXPORT BoundExpressions {
  /// The expressions
  std::vector<NamedExpression> named_expressions;
  /// The schema that all the expressions are bound to
  std::shared_ptr<Schema> schema;
};

}  // namespace engine
}  // namespace arrow

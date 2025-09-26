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

// This API is EXPERIMENTAL.

#pragma once

#include <functional>
#include <string>
#include <vector>

#include "arrow/acero/exec_plan.h"
#include "arrow/acero/options.h"
#include "arrow/compute/type_fwd.h"
#include "arrow/engine/substrait/type_fwd.h"
#include "arrow/engine/substrait/visibility.h"
#include "arrow/type_fwd.h"

namespace arrow {
namespace engine {

/// How strictly to adhere to the input structure when converting between Substrait and
/// Acero representations of a plan. This allows the user to trade conversion accuracy
/// for performance and lenience.
enum class ARROW_ENGINE_EXPORT ConversionStrictness {
  /// When a primitive is used at the input that doesn't have an exact match at the
  /// output, reject the conversion. This effectively asserts that there is no (known)
  /// information loss in the conversion, and that plans should either round-trip back and
  /// forth exactly or not at all. This option is primarily intended for testing and
  /// debugging.
  EXACT_ROUNDTRIP,

  /// When a primitive is used at the input that doesn't have an exact match at the
  /// output, attempt to model it with some collection of primitives at the output. This
  /// means that even if the incoming plan is completely optimal by some metric, the
  /// returned plan is fairly likely to not be optimal anymore, and round-trips back and
  /// forth may make the plan increasingly suboptimal. However, every primitive at the
  /// output can be (manually) traced back to exactly one primitive at the input, which
  /// may be useful when debugging.
  PRESERVE_STRUCTURE,

  /// Behaves like PRESERVE_STRUCTURE, but prefers performance over structural accuracy.
  /// Basic optimizations *may* be applied, in order to attempt to not regress in terms of
  /// plan performance: if the incoming plan was already aggressively optimized, the goal
  /// is for the output plan to not be less performant. In practical use cases, this is
  /// probably the option you want.
  ///
  /// Note that no guarantees are made on top of PRESERVE_STRUCTURE. Past and future
  /// versions of Arrow may even ignore this option entirely and treat it exactly like
  /// PRESERVE_STRUCTURE.
  BEST_EFFORT,
};

using NamedTableProvider = std::function<Result<acero::Declaration>(
    const std::vector<std::string>&, const Schema&)>;
static NamedTableProvider kDefaultNamedTableProvider;

using NamedTapProvider = std::function<Result<acero::Declaration>(
    const std::string&, std::vector<acero::Declaration::Input>, const std::string&,
    std::shared_ptr<Schema>)>;

class ARROW_ENGINE_EXPORT ExtensionDetails {
 public:
  virtual ~ExtensionDetails() = default;
};

class ARROW_ENGINE_EXPORT ExtensionProvider {
 public:
  virtual ~ExtensionProvider() = default;
  virtual Result<DeclarationInfo> MakeRel(const ConversionOptions& conv_opts,
                                          const std::vector<DeclarationInfo>& inputs,
                                          const ExtensionDetails& ext_details,
                                          const ExtensionSet& ext_set) = 0;
};

/// \brief Get the default extension provider
ARROW_ENGINE_EXPORT std::shared_ptr<ExtensionProvider> default_extension_provider();
/// \brief Set the default extension provider
///
/// \param[in] provider the new provider to be set as default
ARROW_ENGINE_EXPORT void set_default_extension_provider(
    const std::shared_ptr<ExtensionProvider>& provider);

ARROW_ENGINE_EXPORT NamedTapProvider default_named_tap_provider();

ARROW_ENGINE_EXPORT void set_default_named_tap_provider(NamedTapProvider provider);

/// Options that control the conversion between Substrait and Acero representations of a
/// plan.
struct ARROW_ENGINE_EXPORT ConversionOptions {
  ConversionOptions()
      : strictness(ConversionStrictness::BEST_EFFORT),
        named_table_provider(kDefaultNamedTableProvider),
        named_tap_provider(default_named_tap_provider()),
        extension_provider(default_extension_provider()),
        allow_arrow_extensions(false) {}

  /// \brief How strictly the converter should adhere to the structure of the input.
  ConversionStrictness strictness;
  /// \brief A custom strategy to be used for providing named tables
  ///
  /// The default behavior will return an invalid status if the plan has any
  /// named table relations.
  NamedTableProvider named_table_provider;
  /// \brief A custom strategy to be used for obtaining a tap declaration
  ///
  /// The default provider returns an error
  NamedTapProvider named_tap_provider;
  /// \brief A custom strategy to be used for providing relation infos.
  ///
  /// The default behavior will provide for relations known to Arrow.
  std::shared_ptr<ExtensionProvider> extension_provider;
  /// \brief If true then Arrow-specific types and functions will be allowed
  ///
  /// Set to false to create plans that are more likely to be compatible with non-Arrow
  /// engines
  bool allow_arrow_extensions;
};

}  // namespace engine
}  // namespace arrow

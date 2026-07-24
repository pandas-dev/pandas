# Copyright 2017 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This import for backward compatibility only.
from functools import wraps  # noqa: F401 pragma: NO COVER

_CREDENTIALS_FILE_WARNING = """\
The `credentials_file` argument is deprecated because of a potential security risk.

The `google.auth.load_credentials_from_file` method does not validate the credential
configuration. The security risk occurs when a credential configuration is accepted
from a source that is not under your control and used without validation on your side.

If you know that you will be loading credential configurations of a
specific type, it is recommended to use a credential-type-specific
load method.

This will ensure that an unexpected credential type with potential for
malicious intent is not loaded unintentionally. You might still have to do
validation for certain credential types. Please follow the recommendations
for that method. For example, if you want to load only service accounts,
you can create the service account credentials explicitly:

```
from google.cloud.vision_v1 import ImageAnnotatorClient
from google.oauth2 import service_account

credentials = service_account.Credentials.from_service_account_file(filename)
client = ImageAnnotatorClient(credentials=credentials)
```

If you are loading your credential configuration from an untrusted source and have
not mitigated the risks (e.g. by validating the configuration yourself), make
these changes as soon as possible to prevent security risks to your environment.

Regardless of the method used, it is always your responsibility to validate
configurations received from external sources.

Refer to https://cloud.google.com/docs/authentication/external/externally-sourced-credentials
for more details.
"""

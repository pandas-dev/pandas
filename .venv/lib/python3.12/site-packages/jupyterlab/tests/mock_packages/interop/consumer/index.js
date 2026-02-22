// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.

const IMockToken = require('@jupyterlab/mock-token').IMockToken;

module.exports = [
  {
    id: 'ext_leaf',
    requires: [IMockToken],
    autoStart: true,
    activate: function (application, mock_token) {
      console.log('hi from consumer', mock_token);
    }
  }
];

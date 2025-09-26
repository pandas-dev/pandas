// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.
const IMockToken = require('@jupyterlab/mock-token').IMockToken;

module.exports = [
  {
    id: 'mock-provider',
    provides: IMockToken,
    autoStart: true,
    activate: function (application) {
      console.log('hi from provider');
      return { value: 'hello' };
    }
  }
];

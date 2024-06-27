// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.

var Widget = require('@lumino/widgets').Widget;

var factory = {
  safe: true,
  mimeTypes: ['text/plain'],
  defaultRank: 1000,
  createRenderer: function () {
    return new Widget();
  }
};

module.exports = {
  id: '@jupyterlab/mock-mime-extension:plugin',
  mimeType: 'text/plain',
  rendererFactory: factory,
  widgetFactoryOptions: {
    name: 'Test',
    fileExtensions: ['.txt'],
    readOnly: true
  }
};

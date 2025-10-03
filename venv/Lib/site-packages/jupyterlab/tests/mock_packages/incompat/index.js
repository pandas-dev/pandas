// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.

module.exports = [
  {
    id: 'mockextension',
    autoStart: true,
    activate: function (application) {
      // eslint-disable-next-line no-console
      console.log('mock extension activated');
      window.commands = application.commands;
    }
  }
];

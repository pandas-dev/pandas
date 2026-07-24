// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.

import { ILauncher } from '@jupyterlab/launcher';

export default [
  {
    id: 'mockextension',
    requires: [ILauncher],
    autoStart: true,
    activate: function (application, launcher) {
      // eslint-disable-next-line no-console
      console.log('mock extension activated', launcher);
      window.commands = application.commands;
    }
  }
];

// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.

import { ILauncher } from '@jupyterlab/launcher';

const plugin = {
  id: 'test_no_hyphens',
  requires: [ILauncher],
  autoStart: true,
  activate: function (application, launcher) {
    // eslint-disable-next-line no-console
    console.log('test_no_hyphens extension activated', launcher);
    window.commands = application.commands;
  }
};

export default plugin;

// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.

import { IDefaultDrive, Drive } from '@jupyterlab/services';

class CustomDrive extends Drive {
  constructor() {
    super();
    console.log('CustomDrive created');
  }
}

export default [
  {
    id: 'mock-service-manager-extension',
    provides: IDefaultDrive,
    autoStart: true,
    activate: function (application) {
      // eslint-disable-next-line no-console
      console.log('service manager mock extension activated');
      const drive = new CustomDrive();
      return drive;
    }
  }
];

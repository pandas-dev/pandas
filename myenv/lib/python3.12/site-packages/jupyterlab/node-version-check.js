/*
 * Copyright (c) Jupyter Development Team.
 * Distributed under the terms of the Modified BSD License.
 */

const pkg = require('./staging/package.json');

function parser(part) {
  return parseInt(part, 10);
}

const engine = pkg.engines.node.replace('>=', '');
const eparts = engine.split('.').map(parser);

const version = process.version.replace('v', '');
const vparts = version.split('.').map(parser);

// eslint-disable-next-line
console.log('Node', process.version);

if (vparts[0] > eparts[0]) {
  process.exit(0);
}

if (vparts[0] < eparts[0]) {
  process.exit(1);
}

if (vparts[1] > eparts[1]) {
  process.exit(0);
}

if (vparts[1] < eparts[1]) {
  process.exit(1);
}

if (vparts[2] < eparts[1]) {
  process.exit(1);
}

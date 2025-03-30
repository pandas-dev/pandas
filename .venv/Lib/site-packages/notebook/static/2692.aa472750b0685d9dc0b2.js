"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[2692],{

/***/ 22692:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {


// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  diagram: () => (/* binding */ diagram)
});

// EXTERNAL MODULE: ../node_modules/mermaid/dist/mermaid-04fb0060.js + 8 modules
var mermaid_04fb0060 = __webpack_require__(24028);
// EXTERNAL MODULE: ../node_modules/dagre-d3-es/src/graphlib/index.js
var graphlib = __webpack_require__(67406);
// EXTERNAL MODULE: ../node_modules/d3/src/index.js + 102 modules
var src = __webpack_require__(23617);
// EXTERNAL MODULE: ../node_modules/dagre-d3-es/src/dagre/index.js + 64 modules
var dagre = __webpack_require__(7259);
;// CONCATENATED MODULE: ../node_modules/mermaid/node_modules/uuid/dist/esm-browser/stringify.js

/**
 * Convert array of 16 byte values to UUID string format of the form:
 * XXXXXXXX-XXXX-XXXX-XXXX-XXXXXXXXXXXX
 */

const byteToHex = [];

for (let i = 0; i < 256; ++i) {
  byteToHex.push((i + 0x100).toString(16).slice(1));
}

function unsafeStringify(arr, offset = 0) {
  // Note: Be careful editing this code!  It's been tuned for performance
  // and works in ways you may not expect. See https://github.com/uuidjs/uuid/pull/434
  return byteToHex[arr[offset + 0]] + byteToHex[arr[offset + 1]] + byteToHex[arr[offset + 2]] + byteToHex[arr[offset + 3]] + '-' + byteToHex[arr[offset + 4]] + byteToHex[arr[offset + 5]] + '-' + byteToHex[arr[offset + 6]] + byteToHex[arr[offset + 7]] + '-' + byteToHex[arr[offset + 8]] + byteToHex[arr[offset + 9]] + '-' + byteToHex[arr[offset + 10]] + byteToHex[arr[offset + 11]] + byteToHex[arr[offset + 12]] + byteToHex[arr[offset + 13]] + byteToHex[arr[offset + 14]] + byteToHex[arr[offset + 15]];
}

function stringify(arr, offset = 0) {
  const uuid = unsafeStringify(arr, offset); // Consistency check for valid UUID.  If this throws, it's likely due to one
  // of the following:
  // - One or more input array values don't map to a hex octet (leading to
  // "undefined" in the uuid)
  // - Invalid input values for the RFC `version` or `variant` fields

  if (!validate(uuid)) {
    throw TypeError('Stringified UUID is invalid');
  }

  return uuid;
}

/* harmony default export */ const esm_browser_stringify = ((/* unused pure expression or super */ null && (stringify)));
;// CONCATENATED MODULE: ../node_modules/mermaid/node_modules/uuid/dist/esm-browser/regex.js
/* harmony default export */ const regex = (/^(?:[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}|00000000-0000-0000-0000-000000000000)$/i);
;// CONCATENATED MODULE: ../node_modules/mermaid/node_modules/uuid/dist/esm-browser/validate.js


function validate_validate(uuid) {
  return typeof uuid === 'string' && regex.test(uuid);
}

/* harmony default export */ const esm_browser_validate = (validate_validate);
;// CONCATENATED MODULE: ../node_modules/mermaid/node_modules/uuid/dist/esm-browser/parse.js


function parse(uuid) {
  if (!esm_browser_validate(uuid)) {
    throw TypeError('Invalid UUID');
  }

  let v;
  const arr = new Uint8Array(16); // Parse ########-....-....-....-............

  arr[0] = (v = parseInt(uuid.slice(0, 8), 16)) >>> 24;
  arr[1] = v >>> 16 & 0xff;
  arr[2] = v >>> 8 & 0xff;
  arr[3] = v & 0xff; // Parse ........-####-....-....-............

  arr[4] = (v = parseInt(uuid.slice(9, 13), 16)) >>> 8;
  arr[5] = v & 0xff; // Parse ........-....-####-....-............

  arr[6] = (v = parseInt(uuid.slice(14, 18), 16)) >>> 8;
  arr[7] = v & 0xff; // Parse ........-....-....-####-............

  arr[8] = (v = parseInt(uuid.slice(19, 23), 16)) >>> 8;
  arr[9] = v & 0xff; // Parse ........-....-....-....-############
  // (Use "/" to avoid 32-bit truncation when bit-shifting high-order bytes)

  arr[10] = (v = parseInt(uuid.slice(24, 36), 16)) / 0x10000000000 & 0xff;
  arr[11] = v / 0x100000000 & 0xff;
  arr[12] = v >>> 24 & 0xff;
  arr[13] = v >>> 16 & 0xff;
  arr[14] = v >>> 8 & 0xff;
  arr[15] = v & 0xff;
  return arr;
}

/* harmony default export */ const esm_browser_parse = (parse);
;// CONCATENATED MODULE: ../node_modules/mermaid/node_modules/uuid/dist/esm-browser/v35.js



function stringToBytes(str) {
  str = unescape(encodeURIComponent(str)); // UTF8 escape

  const bytes = [];

  for (let i = 0; i < str.length; ++i) {
    bytes.push(str.charCodeAt(i));
  }

  return bytes;
}

const DNS = '6ba7b810-9dad-11d1-80b4-00c04fd430c8';
const URL = '6ba7b811-9dad-11d1-80b4-00c04fd430c8';
function v35(name, version, hashfunc) {
  function generateUUID(value, namespace, buf, offset) {
    var _namespace;

    if (typeof value === 'string') {
      value = stringToBytes(value);
    }

    if (typeof namespace === 'string') {
      namespace = esm_browser_parse(namespace);
    }

    if (((_namespace = namespace) === null || _namespace === void 0 ? void 0 : _namespace.length) !== 16) {
      throw TypeError('Namespace must be array-like (16 iterable integer values, 0-255)');
    } // Compute hash of namespace and value, Per 4.3
    // Future: Use spread syntax when supported on all platforms, e.g. `bytes =
    // hashfunc([...namespace, ... value])`


    let bytes = new Uint8Array(16 + value.length);
    bytes.set(namespace);
    bytes.set(value, namespace.length);
    bytes = hashfunc(bytes);
    bytes[6] = bytes[6] & 0x0f | version;
    bytes[8] = bytes[8] & 0x3f | 0x80;

    if (buf) {
      offset = offset || 0;

      for (let i = 0; i < 16; ++i) {
        buf[offset + i] = bytes[i];
      }

      return buf;
    }

    return unsafeStringify(bytes);
  } // Function#name is not settable on some platforms (#270)


  try {
    generateUUID.name = name; // eslint-disable-next-line no-empty
  } catch (err) {} // For CommonJS default export support


  generateUUID.DNS = DNS;
  generateUUID.URL = URL;
  return generateUUID;
}
;// CONCATENATED MODULE: ../node_modules/mermaid/node_modules/uuid/dist/esm-browser/sha1.js
// Adapted from Chris Veness' SHA1 code at
// http://www.movable-type.co.uk/scripts/sha1.html
function f(s, x, y, z) {
  switch (s) {
    case 0:
      return x & y ^ ~x & z;

    case 1:
      return x ^ y ^ z;

    case 2:
      return x & y ^ x & z ^ y & z;

    case 3:
      return x ^ y ^ z;
  }
}

function ROTL(x, n) {
  return x << n | x >>> 32 - n;
}

function sha1(bytes) {
  const K = [0x5a827999, 0x6ed9eba1, 0x8f1bbcdc, 0xca62c1d6];
  const H = [0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476, 0xc3d2e1f0];

  if (typeof bytes === 'string') {
    const msg = unescape(encodeURIComponent(bytes)); // UTF8 escape

    bytes = [];

    for (let i = 0; i < msg.length; ++i) {
      bytes.push(msg.charCodeAt(i));
    }
  } else if (!Array.isArray(bytes)) {
    // Convert Array-like to Array
    bytes = Array.prototype.slice.call(bytes);
  }

  bytes.push(0x80);
  const l = bytes.length / 4 + 2;
  const N = Math.ceil(l / 16);
  const M = new Array(N);

  for (let i = 0; i < N; ++i) {
    const arr = new Uint32Array(16);

    for (let j = 0; j < 16; ++j) {
      arr[j] = bytes[i * 64 + j * 4] << 24 | bytes[i * 64 + j * 4 + 1] << 16 | bytes[i * 64 + j * 4 + 2] << 8 | bytes[i * 64 + j * 4 + 3];
    }

    M[i] = arr;
  }

  M[N - 1][14] = (bytes.length - 1) * 8 / Math.pow(2, 32);
  M[N - 1][14] = Math.floor(M[N - 1][14]);
  M[N - 1][15] = (bytes.length - 1) * 8 & 0xffffffff;

  for (let i = 0; i < N; ++i) {
    const W = new Uint32Array(80);

    for (let t = 0; t < 16; ++t) {
      W[t] = M[i][t];
    }

    for (let t = 16; t < 80; ++t) {
      W[t] = ROTL(W[t - 3] ^ W[t - 8] ^ W[t - 14] ^ W[t - 16], 1);
    }

    let a = H[0];
    let b = H[1];
    let c = H[2];
    let d = H[3];
    let e = H[4];

    for (let t = 0; t < 80; ++t) {
      const s = Math.floor(t / 20);
      const T = ROTL(a, 5) + f(s, b, c, d) + e + K[s] + W[t] >>> 0;
      e = d;
      d = c;
      c = ROTL(b, 30) >>> 0;
      b = a;
      a = T;
    }

    H[0] = H[0] + a >>> 0;
    H[1] = H[1] + b >>> 0;
    H[2] = H[2] + c >>> 0;
    H[3] = H[3] + d >>> 0;
    H[4] = H[4] + e >>> 0;
  }

  return [H[0] >> 24 & 0xff, H[0] >> 16 & 0xff, H[0] >> 8 & 0xff, H[0] & 0xff, H[1] >> 24 & 0xff, H[1] >> 16 & 0xff, H[1] >> 8 & 0xff, H[1] & 0xff, H[2] >> 24 & 0xff, H[2] >> 16 & 0xff, H[2] >> 8 & 0xff, H[2] & 0xff, H[3] >> 24 & 0xff, H[3] >> 16 & 0xff, H[3] >> 8 & 0xff, H[3] & 0xff, H[4] >> 24 & 0xff, H[4] >> 16 & 0xff, H[4] >> 8 & 0xff, H[4] & 0xff];
}

/* harmony default export */ const esm_browser_sha1 = (sha1);
;// CONCATENATED MODULE: ../node_modules/mermaid/node_modules/uuid/dist/esm-browser/v5.js


const v5 = v35('v5', 0x50, esm_browser_sha1);
/* harmony default export */ const esm_browser_v5 = (v5);
// EXTERNAL MODULE: ../node_modules/dayjs/dayjs.min.js
var dayjs_min = __webpack_require__(27693);
// EXTERNAL MODULE: ../node_modules/@braintree/sanitize-url/dist/index.js
var dist = __webpack_require__(7608);
// EXTERNAL MODULE: ../node_modules/dompurify/dist/purify.js
var purify = __webpack_require__(31699);
;// CONCATENATED MODULE: ../node_modules/mermaid/dist/erDiagram-4c71174c.js














var parser = function() {
  var o = function(k, v, o2, l) {
    for (o2 = o2 || {}, l = k.length; l--; o2[k[l]] = v)
      ;
    return o2;
  }, $V0 = [6, 8, 10, 20, 22, 24, 26, 27, 28], $V1 = [1, 10], $V2 = [1, 11], $V3 = [1, 12], $V4 = [1, 13], $V5 = [1, 14], $V6 = [1, 15], $V7 = [1, 21], $V8 = [1, 22], $V9 = [1, 23], $Va = [1, 24], $Vb = [1, 25], $Vc = [6, 8, 10, 13, 15, 18, 19, 20, 22, 24, 26, 27, 28, 41, 42, 43, 44, 45], $Vd = [1, 34], $Ve = [27, 28, 46, 47], $Vf = [41, 42, 43, 44, 45], $Vg = [17, 34], $Vh = [1, 54], $Vi = [1, 53], $Vj = [17, 34, 36, 38];
  var parser2 = {
    trace: function trace() {
    },
    yy: {},
    symbols_: { "error": 2, "start": 3, "ER_DIAGRAM": 4, "document": 5, "EOF": 6, "line": 7, "SPACE": 8, "statement": 9, "NEWLINE": 10, "entityName": 11, "relSpec": 12, ":": 13, "role": 14, "BLOCK_START": 15, "attributes": 16, "BLOCK_STOP": 17, "SQS": 18, "SQE": 19, "title": 20, "title_value": 21, "acc_title": 22, "acc_title_value": 23, "acc_descr": 24, "acc_descr_value": 25, "acc_descr_multiline_value": 26, "ALPHANUM": 27, "ENTITY_NAME": 28, "attribute": 29, "attributeType": 30, "attributeName": 31, "attributeKeyTypeList": 32, "attributeComment": 33, "ATTRIBUTE_WORD": 34, "attributeKeyType": 35, "COMMA": 36, "ATTRIBUTE_KEY": 37, "COMMENT": 38, "cardinality": 39, "relType": 40, "ZERO_OR_ONE": 41, "ZERO_OR_MORE": 42, "ONE_OR_MORE": 43, "ONLY_ONE": 44, "MD_PARENT": 45, "NON_IDENTIFYING": 46, "IDENTIFYING": 47, "WORD": 48, "$accept": 0, "$end": 1 },
    terminals_: { 2: "error", 4: "ER_DIAGRAM", 6: "EOF", 8: "SPACE", 10: "NEWLINE", 13: ":", 15: "BLOCK_START", 17: "BLOCK_STOP", 18: "SQS", 19: "SQE", 20: "title", 21: "title_value", 22: "acc_title", 23: "acc_title_value", 24: "acc_descr", 25: "acc_descr_value", 26: "acc_descr_multiline_value", 27: "ALPHANUM", 28: "ENTITY_NAME", 34: "ATTRIBUTE_WORD", 36: "COMMA", 37: "ATTRIBUTE_KEY", 38: "COMMENT", 41: "ZERO_OR_ONE", 42: "ZERO_OR_MORE", 43: "ONE_OR_MORE", 44: "ONLY_ONE", 45: "MD_PARENT", 46: "NON_IDENTIFYING", 47: "IDENTIFYING", 48: "WORD" },
    productions_: [0, [3, 3], [5, 0], [5, 2], [7, 2], [7, 1], [7, 1], [7, 1], [9, 5], [9, 4], [9, 3], [9, 1], [9, 7], [9, 6], [9, 4], [9, 2], [9, 2], [9, 2], [9, 1], [11, 1], [11, 1], [16, 1], [16, 2], [29, 2], [29, 3], [29, 3], [29, 4], [30, 1], [31, 1], [32, 1], [32, 3], [35, 1], [33, 1], [12, 3], [39, 1], [39, 1], [39, 1], [39, 1], [39, 1], [40, 1], [40, 1], [14, 1], [14, 1], [14, 1]],
    performAction: function anonymous(yytext, yyleng, yylineno, yy, yystate, $$, _$) {
      var $0 = $$.length - 1;
      switch (yystate) {
        case 1:
          break;
        case 2:
          this.$ = [];
          break;
        case 3:
          $$[$0 - 1].push($$[$0]);
          this.$ = $$[$0 - 1];
          break;
        case 4:
        case 5:
          this.$ = $$[$0];
          break;
        case 6:
        case 7:
          this.$ = [];
          break;
        case 8:
          yy.addEntity($$[$0 - 4]);
          yy.addEntity($$[$0 - 2]);
          yy.addRelationship($$[$0 - 4], $$[$0], $$[$0 - 2], $$[$0 - 3]);
          break;
        case 9:
          yy.addEntity($$[$0 - 3]);
          yy.addAttributes($$[$0 - 3], $$[$0 - 1]);
          break;
        case 10:
          yy.addEntity($$[$0 - 2]);
          break;
        case 11:
          yy.addEntity($$[$0]);
          break;
        case 12:
          yy.addEntity($$[$0 - 6], $$[$0 - 4]);
          yy.addAttributes($$[$0 - 6], $$[$0 - 1]);
          break;
        case 13:
          yy.addEntity($$[$0 - 5], $$[$0 - 3]);
          break;
        case 14:
          yy.addEntity($$[$0 - 3], $$[$0 - 1]);
          break;
        case 15:
        case 16:
          this.$ = $$[$0].trim();
          yy.setAccTitle(this.$);
          break;
        case 17:
        case 18:
          this.$ = $$[$0].trim();
          yy.setAccDescription(this.$);
          break;
        case 19:
        case 43:
          this.$ = $$[$0];
          break;
        case 20:
        case 41:
        case 42:
          this.$ = $$[$0].replace(/"/g, "");
          break;
        case 21:
        case 29:
          this.$ = [$$[$0]];
          break;
        case 22:
          $$[$0].push($$[$0 - 1]);
          this.$ = $$[$0];
          break;
        case 23:
          this.$ = { attributeType: $$[$0 - 1], attributeName: $$[$0] };
          break;
        case 24:
          this.$ = { attributeType: $$[$0 - 2], attributeName: $$[$0 - 1], attributeKeyTypeList: $$[$0] };
          break;
        case 25:
          this.$ = { attributeType: $$[$0 - 2], attributeName: $$[$0 - 1], attributeComment: $$[$0] };
          break;
        case 26:
          this.$ = { attributeType: $$[$0 - 3], attributeName: $$[$0 - 2], attributeKeyTypeList: $$[$0 - 1], attributeComment: $$[$0] };
          break;
        case 27:
        case 28:
        case 31:
          this.$ = $$[$0];
          break;
        case 30:
          $$[$0 - 2].push($$[$0]);
          this.$ = $$[$0 - 2];
          break;
        case 32:
          this.$ = $$[$0].replace(/"/g, "");
          break;
        case 33:
          this.$ = { cardA: $$[$0], relType: $$[$0 - 1], cardB: $$[$0 - 2] };
          break;
        case 34:
          this.$ = yy.Cardinality.ZERO_OR_ONE;
          break;
        case 35:
          this.$ = yy.Cardinality.ZERO_OR_MORE;
          break;
        case 36:
          this.$ = yy.Cardinality.ONE_OR_MORE;
          break;
        case 37:
          this.$ = yy.Cardinality.ONLY_ONE;
          break;
        case 38:
          this.$ = yy.Cardinality.MD_PARENT;
          break;
        case 39:
          this.$ = yy.Identification.NON_IDENTIFYING;
          break;
        case 40:
          this.$ = yy.Identification.IDENTIFYING;
          break;
      }
    },
    table: [{ 3: 1, 4: [1, 2] }, { 1: [3] }, o($V0, [2, 2], { 5: 3 }), { 6: [1, 4], 7: 5, 8: [1, 6], 9: 7, 10: [1, 8], 11: 9, 20: $V1, 22: $V2, 24: $V3, 26: $V4, 27: $V5, 28: $V6 }, o($V0, [2, 7], { 1: [2, 1] }), o($V0, [2, 3]), { 9: 16, 11: 9, 20: $V1, 22: $V2, 24: $V3, 26: $V4, 27: $V5, 28: $V6 }, o($V0, [2, 5]), o($V0, [2, 6]), o($V0, [2, 11], { 12: 17, 39: 20, 15: [1, 18], 18: [1, 19], 41: $V7, 42: $V8, 43: $V9, 44: $Va, 45: $Vb }), { 21: [1, 26] }, { 23: [1, 27] }, { 25: [1, 28] }, o($V0, [2, 18]), o($Vc, [2, 19]), o($Vc, [2, 20]), o($V0, [2, 4]), { 11: 29, 27: $V5, 28: $V6 }, { 16: 30, 17: [1, 31], 29: 32, 30: 33, 34: $Vd }, { 11: 35, 27: $V5, 28: $V6 }, { 40: 36, 46: [1, 37], 47: [1, 38] }, o($Ve, [2, 34]), o($Ve, [2, 35]), o($Ve, [2, 36]), o($Ve, [2, 37]), o($Ve, [2, 38]), o($V0, [2, 15]), o($V0, [2, 16]), o($V0, [2, 17]), { 13: [1, 39] }, { 17: [1, 40] }, o($V0, [2, 10]), { 16: 41, 17: [2, 21], 29: 32, 30: 33, 34: $Vd }, { 31: 42, 34: [1, 43] }, { 34: [2, 27] }, { 19: [1, 44] }, { 39: 45, 41: $V7, 42: $V8, 43: $V9, 44: $Va, 45: $Vb }, o($Vf, [2, 39]), o($Vf, [2, 40]), { 14: 46, 27: [1, 49], 28: [1, 48], 48: [1, 47] }, o($V0, [2, 9]), { 17: [2, 22] }, o($Vg, [2, 23], { 32: 50, 33: 51, 35: 52, 37: $Vh, 38: $Vi }), o([17, 34, 37, 38], [2, 28]), o($V0, [2, 14], { 15: [1, 55] }), o([27, 28], [2, 33]), o($V0, [2, 8]), o($V0, [2, 41]), o($V0, [2, 42]), o($V0, [2, 43]), o($Vg, [2, 24], { 33: 56, 36: [1, 57], 38: $Vi }), o($Vg, [2, 25]), o($Vj, [2, 29]), o($Vg, [2, 32]), o($Vj, [2, 31]), { 16: 58, 17: [1, 59], 29: 32, 30: 33, 34: $Vd }, o($Vg, [2, 26]), { 35: 60, 37: $Vh }, { 17: [1, 61] }, o($V0, [2, 13]), o($Vj, [2, 30]), o($V0, [2, 12])],
    defaultActions: { 34: [2, 27], 41: [2, 22] },
    parseError: function parseError(str, hash) {
      if (hash.recoverable) {
        this.trace(str);
      } else {
        var error = new Error(str);
        error.hash = hash;
        throw error;
      }
    },
    parse: function parse(input) {
      var self = this, stack = [0], tstack = [], vstack = [null], lstack = [], table = this.table, yytext = "", yylineno = 0, yyleng = 0, TERROR = 2, EOF = 1;
      var args = lstack.slice.call(arguments, 1);
      var lexer2 = Object.create(this.lexer);
      var sharedState = { yy: {} };
      for (var k in this.yy) {
        if (Object.prototype.hasOwnProperty.call(this.yy, k)) {
          sharedState.yy[k] = this.yy[k];
        }
      }
      lexer2.setInput(input, sharedState.yy);
      sharedState.yy.lexer = lexer2;
      sharedState.yy.parser = this;
      if (typeof lexer2.yylloc == "undefined") {
        lexer2.yylloc = {};
      }
      var yyloc = lexer2.yylloc;
      lstack.push(yyloc);
      var ranges = lexer2.options && lexer2.options.ranges;
      if (typeof sharedState.yy.parseError === "function") {
        this.parseError = sharedState.yy.parseError;
      } else {
        this.parseError = Object.getPrototypeOf(this).parseError;
      }
      function lex() {
        var token;
        token = tstack.pop() || lexer2.lex() || EOF;
        if (typeof token !== "number") {
          if (token instanceof Array) {
            tstack = token;
            token = tstack.pop();
          }
          token = self.symbols_[token] || token;
        }
        return token;
      }
      var symbol, state, action, r, yyval = {}, p, len, newState, expected;
      while (true) {
        state = stack[stack.length - 1];
        if (this.defaultActions[state]) {
          action = this.defaultActions[state];
        } else {
          if (symbol === null || typeof symbol == "undefined") {
            symbol = lex();
          }
          action = table[state] && table[state][symbol];
        }
        if (typeof action === "undefined" || !action.length || !action[0]) {
          var errStr = "";
          expected = [];
          for (p in table[state]) {
            if (this.terminals_[p] && p > TERROR) {
              expected.push("'" + this.terminals_[p] + "'");
            }
          }
          if (lexer2.showPosition) {
            errStr = "Parse error on line " + (yylineno + 1) + ":\n" + lexer2.showPosition() + "\nExpecting " + expected.join(", ") + ", got '" + (this.terminals_[symbol] || symbol) + "'";
          } else {
            errStr = "Parse error on line " + (yylineno + 1) + ": Unexpected " + (symbol == EOF ? "end of input" : "'" + (this.terminals_[symbol] || symbol) + "'");
          }
          this.parseError(errStr, {
            text: lexer2.match,
            token: this.terminals_[symbol] || symbol,
            line: lexer2.yylineno,
            loc: yyloc,
            expected
          });
        }
        if (action[0] instanceof Array && action.length > 1) {
          throw new Error("Parse Error: multiple actions possible at state: " + state + ", token: " + symbol);
        }
        switch (action[0]) {
          case 1:
            stack.push(symbol);
            vstack.push(lexer2.yytext);
            lstack.push(lexer2.yylloc);
            stack.push(action[1]);
            symbol = null;
            {
              yyleng = lexer2.yyleng;
              yytext = lexer2.yytext;
              yylineno = lexer2.yylineno;
              yyloc = lexer2.yylloc;
            }
            break;
          case 2:
            len = this.productions_[action[1]][1];
            yyval.$ = vstack[vstack.length - len];
            yyval._$ = {
              first_line: lstack[lstack.length - (len || 1)].first_line,
              last_line: lstack[lstack.length - 1].last_line,
              first_column: lstack[lstack.length - (len || 1)].first_column,
              last_column: lstack[lstack.length - 1].last_column
            };
            if (ranges) {
              yyval._$.range = [
                lstack[lstack.length - (len || 1)].range[0],
                lstack[lstack.length - 1].range[1]
              ];
            }
            r = this.performAction.apply(yyval, [
              yytext,
              yyleng,
              yylineno,
              sharedState.yy,
              action[1],
              vstack,
              lstack
            ].concat(args));
            if (typeof r !== "undefined") {
              return r;
            }
            if (len) {
              stack = stack.slice(0, -1 * len * 2);
              vstack = vstack.slice(0, -1 * len);
              lstack = lstack.slice(0, -1 * len);
            }
            stack.push(this.productions_[action[1]][0]);
            vstack.push(yyval.$);
            lstack.push(yyval._$);
            newState = table[stack[stack.length - 2]][stack[stack.length - 1]];
            stack.push(newState);
            break;
          case 3:
            return true;
        }
      }
      return true;
    }
  };
  var lexer = function() {
    var lexer2 = {
      EOF: 1,
      parseError: function parseError(str, hash) {
        if (this.yy.parser) {
          this.yy.parser.parseError(str, hash);
        } else {
          throw new Error(str);
        }
      },
      // resets the lexer, sets new input
      setInput: function(input, yy) {
        this.yy = yy || this.yy || {};
        this._input = input;
        this._more = this._backtrack = this.done = false;
        this.yylineno = this.yyleng = 0;
        this.yytext = this.matched = this.match = "";
        this.conditionStack = ["INITIAL"];
        this.yylloc = {
          first_line: 1,
          first_column: 0,
          last_line: 1,
          last_column: 0
        };
        if (this.options.ranges) {
          this.yylloc.range = [0, 0];
        }
        this.offset = 0;
        return this;
      },
      // consumes and returns one char from the input
      input: function() {
        var ch = this._input[0];
        this.yytext += ch;
        this.yyleng++;
        this.offset++;
        this.match += ch;
        this.matched += ch;
        var lines = ch.match(/(?:\r\n?|\n).*/g);
        if (lines) {
          this.yylineno++;
          this.yylloc.last_line++;
        } else {
          this.yylloc.last_column++;
        }
        if (this.options.ranges) {
          this.yylloc.range[1]++;
        }
        this._input = this._input.slice(1);
        return ch;
      },
      // unshifts one char (or a string) into the input
      unput: function(ch) {
        var len = ch.length;
        var lines = ch.split(/(?:\r\n?|\n)/g);
        this._input = ch + this._input;
        this.yytext = this.yytext.substr(0, this.yytext.length - len);
        this.offset -= len;
        var oldLines = this.match.split(/(?:\r\n?|\n)/g);
        this.match = this.match.substr(0, this.match.length - 1);
        this.matched = this.matched.substr(0, this.matched.length - 1);
        if (lines.length - 1) {
          this.yylineno -= lines.length - 1;
        }
        var r = this.yylloc.range;
        this.yylloc = {
          first_line: this.yylloc.first_line,
          last_line: this.yylineno + 1,
          first_column: this.yylloc.first_column,
          last_column: lines ? (lines.length === oldLines.length ? this.yylloc.first_column : 0) + oldLines[oldLines.length - lines.length].length - lines[0].length : this.yylloc.first_column - len
        };
        if (this.options.ranges) {
          this.yylloc.range = [r[0], r[0] + this.yyleng - len];
        }
        this.yyleng = this.yytext.length;
        return this;
      },
      // When called from action, caches matched text and appends it on next action
      more: function() {
        this._more = true;
        return this;
      },
      // When called from action, signals the lexer that this rule fails to match the input, so the next matching rule (regex) should be tested instead.
      reject: function() {
        if (this.options.backtrack_lexer) {
          this._backtrack = true;
        } else {
          return this.parseError("Lexical error on line " + (this.yylineno + 1) + ". You can only invoke reject() in the lexer when the lexer is of the backtracking persuasion (options.backtrack_lexer = true).\n" + this.showPosition(), {
            text: "",
            token: null,
            line: this.yylineno
          });
        }
        return this;
      },
      // retain first n characters of the match
      less: function(n) {
        this.unput(this.match.slice(n));
      },
      // displays already matched input, i.e. for error messages
      pastInput: function() {
        var past = this.matched.substr(0, this.matched.length - this.match.length);
        return (past.length > 20 ? "..." : "") + past.substr(-20).replace(/\n/g, "");
      },
      // displays upcoming input, i.e. for error messages
      upcomingInput: function() {
        var next = this.match;
        if (next.length < 20) {
          next += this._input.substr(0, 20 - next.length);
        }
        return (next.substr(0, 20) + (next.length > 20 ? "..." : "")).replace(/\n/g, "");
      },
      // displays the character position where the lexing error occurred, i.e. for error messages
      showPosition: function() {
        var pre = this.pastInput();
        var c = new Array(pre.length + 1).join("-");
        return pre + this.upcomingInput() + "\n" + c + "^";
      },
      // test the lexed token: return FALSE when not a match, otherwise return token
      test_match: function(match, indexed_rule) {
        var token, lines, backup;
        if (this.options.backtrack_lexer) {
          backup = {
            yylineno: this.yylineno,
            yylloc: {
              first_line: this.yylloc.first_line,
              last_line: this.last_line,
              first_column: this.yylloc.first_column,
              last_column: this.yylloc.last_column
            },
            yytext: this.yytext,
            match: this.match,
            matches: this.matches,
            matched: this.matched,
            yyleng: this.yyleng,
            offset: this.offset,
            _more: this._more,
            _input: this._input,
            yy: this.yy,
            conditionStack: this.conditionStack.slice(0),
            done: this.done
          };
          if (this.options.ranges) {
            backup.yylloc.range = this.yylloc.range.slice(0);
          }
        }
        lines = match[0].match(/(?:\r\n?|\n).*/g);
        if (lines) {
          this.yylineno += lines.length;
        }
        this.yylloc = {
          first_line: this.yylloc.last_line,
          last_line: this.yylineno + 1,
          first_column: this.yylloc.last_column,
          last_column: lines ? lines[lines.length - 1].length - lines[lines.length - 1].match(/\r?\n?/)[0].length : this.yylloc.last_column + match[0].length
        };
        this.yytext += match[0];
        this.match += match[0];
        this.matches = match;
        this.yyleng = this.yytext.length;
        if (this.options.ranges) {
          this.yylloc.range = [this.offset, this.offset += this.yyleng];
        }
        this._more = false;
        this._backtrack = false;
        this._input = this._input.slice(match[0].length);
        this.matched += match[0];
        token = this.performAction.call(this, this.yy, this, indexed_rule, this.conditionStack[this.conditionStack.length - 1]);
        if (this.done && this._input) {
          this.done = false;
        }
        if (token) {
          return token;
        } else if (this._backtrack) {
          for (var k in backup) {
            this[k] = backup[k];
          }
          return false;
        }
        return false;
      },
      // return next match in input
      next: function() {
        if (this.done) {
          return this.EOF;
        }
        if (!this._input) {
          this.done = true;
        }
        var token, match, tempMatch, index;
        if (!this._more) {
          this.yytext = "";
          this.match = "";
        }
        var rules = this._currentRules();
        for (var i = 0; i < rules.length; i++) {
          tempMatch = this._input.match(this.rules[rules[i]]);
          if (tempMatch && (!match || tempMatch[0].length > match[0].length)) {
            match = tempMatch;
            index = i;
            if (this.options.backtrack_lexer) {
              token = this.test_match(tempMatch, rules[i]);
              if (token !== false) {
                return token;
              } else if (this._backtrack) {
                match = false;
                continue;
              } else {
                return false;
              }
            } else if (!this.options.flex) {
              break;
            }
          }
        }
        if (match) {
          token = this.test_match(match, rules[index]);
          if (token !== false) {
            return token;
          }
          return false;
        }
        if (this._input === "") {
          return this.EOF;
        } else {
          return this.parseError("Lexical error on line " + (this.yylineno + 1) + ". Unrecognized text.\n" + this.showPosition(), {
            text: "",
            token: null,
            line: this.yylineno
          });
        }
      },
      // return next match that has a token
      lex: function lex() {
        var r = this.next();
        if (r) {
          return r;
        } else {
          return this.lex();
        }
      },
      // activates a new lexer condition state (pushes the new lexer condition state onto the condition stack)
      begin: function begin(condition) {
        this.conditionStack.push(condition);
      },
      // pop the previously active lexer condition state off the condition stack
      popState: function popState() {
        var n = this.conditionStack.length - 1;
        if (n > 0) {
          return this.conditionStack.pop();
        } else {
          return this.conditionStack[0];
        }
      },
      // produce the lexer rule set which is active for the currently active lexer condition state
      _currentRules: function _currentRules() {
        if (this.conditionStack.length && this.conditionStack[this.conditionStack.length - 1]) {
          return this.conditions[this.conditionStack[this.conditionStack.length - 1]].rules;
        } else {
          return this.conditions["INITIAL"].rules;
        }
      },
      // return the currently active lexer condition state; when an index argument is provided it produces the N-th previous condition state, if available
      topState: function topState(n) {
        n = this.conditionStack.length - 1 - Math.abs(n || 0);
        if (n >= 0) {
          return this.conditionStack[n];
        } else {
          return "INITIAL";
        }
      },
      // alias for begin(condition)
      pushState: function pushState(condition) {
        this.begin(condition);
      },
      // return the number of states currently on the stack
      stateStackSize: function stateStackSize() {
        return this.conditionStack.length;
      },
      options: { "case-insensitive": true },
      performAction: function anonymous(yy, yy_, $avoiding_name_collisions, YY_START) {
        switch ($avoiding_name_collisions) {
          case 0:
            this.begin("acc_title");
            return 22;
          case 1:
            this.popState();
            return "acc_title_value";
          case 2:
            this.begin("acc_descr");
            return 24;
          case 3:
            this.popState();
            return "acc_descr_value";
          case 4:
            this.begin("acc_descr_multiline");
            break;
          case 5:
            this.popState();
            break;
          case 6:
            return "acc_descr_multiline_value";
          case 7:
            return 10;
          case 8:
            break;
          case 9:
            return 8;
          case 10:
            return 28;
          case 11:
            return 48;
          case 12:
            return 4;
          case 13:
            this.begin("block");
            return 15;
          case 14:
            return 36;
          case 15:
            break;
          case 16:
            return 37;
          case 17:
            return 34;
          case 18:
            return 34;
          case 19:
            return 38;
          case 20:
            break;
          case 21:
            this.popState();
            return 17;
          case 22:
            return yy_.yytext[0];
          case 23:
            return 18;
          case 24:
            return 19;
          case 25:
            return 41;
          case 26:
            return 43;
          case 27:
            return 43;
          case 28:
            return 43;
          case 29:
            return 41;
          case 30:
            return 41;
          case 31:
            return 42;
          case 32:
            return 42;
          case 33:
            return 42;
          case 34:
            return 42;
          case 35:
            return 42;
          case 36:
            return 43;
          case 37:
            return 42;
          case 38:
            return 43;
          case 39:
            return 44;
          case 40:
            return 44;
          case 41:
            return 44;
          case 42:
            return 44;
          case 43:
            return 41;
          case 44:
            return 42;
          case 45:
            return 43;
          case 46:
            return 45;
          case 47:
            return 46;
          case 48:
            return 47;
          case 49:
            return 47;
          case 50:
            return 46;
          case 51:
            return 46;
          case 52:
            return 46;
          case 53:
            return 27;
          case 54:
            return yy_.yytext[0];
          case 55:
            return 6;
        }
      },
      rules: [/^(?:accTitle\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*\{\s*)/i, /^(?:[\}])/i, /^(?:[^\}]*)/i, /^(?:[\n]+)/i, /^(?:\s+)/i, /^(?:[\s]+)/i, /^(?:"[^"%\r\n\v\b\\]+")/i, /^(?:"[^"]*")/i, /^(?:erDiagram\b)/i, /^(?:\{)/i, /^(?:,)/i, /^(?:\s+)/i, /^(?:\b((?:PK)|(?:FK)|(?:UK))\b)/i, /^(?:(.*?)[~](.*?)*[~])/i, /^(?:[\*A-Za-z_][A-Za-z0-9\-_\[\]\(\)]*)/i, /^(?:"[^"]*")/i, /^(?:[\n]+)/i, /^(?:\})/i, /^(?:.)/i, /^(?:\[)/i, /^(?:\])/i, /^(?:one or zero\b)/i, /^(?:one or more\b)/i, /^(?:one or many\b)/i, /^(?:1\+)/i, /^(?:\|o\b)/i, /^(?:zero or one\b)/i, /^(?:zero or more\b)/i, /^(?:zero or many\b)/i, /^(?:0\+)/i, /^(?:\}o\b)/i, /^(?:many\(0\))/i, /^(?:many\(1\))/i, /^(?:many\b)/i, /^(?:\}\|)/i, /^(?:one\b)/i, /^(?:only one\b)/i, /^(?:1\b)/i, /^(?:\|\|)/i, /^(?:o\|)/i, /^(?:o\{)/i, /^(?:\|\{)/i, /^(?:\s*u\b)/i, /^(?:\.\.)/i, /^(?:--)/i, /^(?:to\b)/i, /^(?:optionally to\b)/i, /^(?:\.-)/i, /^(?:-\.)/i, /^(?:[A-Za-z_][A-Za-z0-9\-_]*)/i, /^(?:.)/i, /^(?:$)/i],
      conditions: { "acc_descr_multiline": { "rules": [5, 6], "inclusive": false }, "acc_descr": { "rules": [3], "inclusive": false }, "acc_title": { "rules": [1], "inclusive": false }, "block": { "rules": [14, 15, 16, 17, 18, 19, 20, 21, 22], "inclusive": false }, "INITIAL": { "rules": [0, 2, 4, 7, 8, 9, 10, 11, 12, 13, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55], "inclusive": true } }
    };
    return lexer2;
  }();
  parser2.lexer = lexer;
  function Parser() {
    this.yy = {};
  }
  Parser.prototype = parser2;
  parser2.Parser = Parser;
  return new Parser();
}();
parser.parser = parser;
const erParser = parser;
let entities = {};
let relationships = [];
const Cardinality = {
  ZERO_OR_ONE: "ZERO_OR_ONE",
  ZERO_OR_MORE: "ZERO_OR_MORE",
  ONE_OR_MORE: "ONE_OR_MORE",
  ONLY_ONE: "ONLY_ONE",
  MD_PARENT: "MD_PARENT"
};
const Identification = {
  NON_IDENTIFYING: "NON_IDENTIFYING",
  IDENTIFYING: "IDENTIFYING"
};
const addEntity = function(name, alias = void 0) {
  if (entities[name] === void 0) {
    entities[name] = { attributes: [], alias };
    mermaid_04fb0060.l.info("Added new entity :", name);
  } else if (entities[name] && !entities[name].alias && alias) {
    entities[name].alias = alias;
    mermaid_04fb0060.l.info(`Add alias '${alias}' to entity '${name}'`);
  }
  return entities[name];
};
const getEntities = () => entities;
const addAttributes = function(entityName, attribs) {
  let entity = addEntity(entityName);
  let i;
  for (i = attribs.length - 1; i >= 0; i--) {
    entity.attributes.push(attribs[i]);
    mermaid_04fb0060.l.debug("Added attribute ", attribs[i].attributeName);
  }
};
const addRelationship = function(entA, rolA, entB, rSpec) {
  let rel = {
    entityA: entA,
    roleA: rolA,
    entityB: entB,
    relSpec: rSpec
  };
  relationships.push(rel);
  mermaid_04fb0060.l.debug("Added new relationship :", rel);
};
const getRelationships = () => relationships;
const clear = function() {
  entities = {};
  relationships = [];
  (0,mermaid_04fb0060.t)();
};
const erDb = {
  Cardinality,
  Identification,
  getConfig: () => (0,mermaid_04fb0060.c)().er,
  addEntity,
  addAttributes,
  getEntities,
  addRelationship,
  getRelationships,
  clear,
  setAccTitle: mermaid_04fb0060.s,
  getAccTitle: mermaid_04fb0060.g,
  setAccDescription: mermaid_04fb0060.b,
  getAccDescription: mermaid_04fb0060.a,
  setDiagramTitle: mermaid_04fb0060.q,
  getDiagramTitle: mermaid_04fb0060.r
};
const ERMarkers = {
  ONLY_ONE_START: "ONLY_ONE_START",
  ONLY_ONE_END: "ONLY_ONE_END",
  ZERO_OR_ONE_START: "ZERO_OR_ONE_START",
  ZERO_OR_ONE_END: "ZERO_OR_ONE_END",
  ONE_OR_MORE_START: "ONE_OR_MORE_START",
  ONE_OR_MORE_END: "ONE_OR_MORE_END",
  ZERO_OR_MORE_START: "ZERO_OR_MORE_START",
  ZERO_OR_MORE_END: "ZERO_OR_MORE_END",
  MD_PARENT_END: "MD_PARENT_END",
  MD_PARENT_START: "MD_PARENT_START"
};
const insertMarkers = function(elem, conf2) {
  let marker;
  elem.append("defs").append("marker").attr("id", ERMarkers.MD_PARENT_START).attr("refX", 0).attr("refY", 7).attr("markerWidth", 190).attr("markerHeight", 240).attr("orient", "auto").append("path").attr("d", "M 18,7 L9,13 L1,7 L9,1 Z");
  elem.append("defs").append("marker").attr("id", ERMarkers.MD_PARENT_END).attr("refX", 19).attr("refY", 7).attr("markerWidth", 20).attr("markerHeight", 28).attr("orient", "auto").append("path").attr("d", "M 18,7 L9,13 L1,7 L9,1 Z");
  elem.append("defs").append("marker").attr("id", ERMarkers.ONLY_ONE_START).attr("refX", 0).attr("refY", 9).attr("markerWidth", 18).attr("markerHeight", 18).attr("orient", "auto").append("path").attr("stroke", conf2.stroke).attr("fill", "none").attr("d", "M9,0 L9,18 M15,0 L15,18");
  elem.append("defs").append("marker").attr("id", ERMarkers.ONLY_ONE_END).attr("refX", 18).attr("refY", 9).attr("markerWidth", 18).attr("markerHeight", 18).attr("orient", "auto").append("path").attr("stroke", conf2.stroke).attr("fill", "none").attr("d", "M3,0 L3,18 M9,0 L9,18");
  marker = elem.append("defs").append("marker").attr("id", ERMarkers.ZERO_OR_ONE_START).attr("refX", 0).attr("refY", 9).attr("markerWidth", 30).attr("markerHeight", 18).attr("orient", "auto");
  marker.append("circle").attr("stroke", conf2.stroke).attr("fill", "white").attr("cx", 21).attr("cy", 9).attr("r", 6);
  marker.append("path").attr("stroke", conf2.stroke).attr("fill", "none").attr("d", "M9,0 L9,18");
  marker = elem.append("defs").append("marker").attr("id", ERMarkers.ZERO_OR_ONE_END).attr("refX", 30).attr("refY", 9).attr("markerWidth", 30).attr("markerHeight", 18).attr("orient", "auto");
  marker.append("circle").attr("stroke", conf2.stroke).attr("fill", "white").attr("cx", 9).attr("cy", 9).attr("r", 6);
  marker.append("path").attr("stroke", conf2.stroke).attr("fill", "none").attr("d", "M21,0 L21,18");
  elem.append("defs").append("marker").attr("id", ERMarkers.ONE_OR_MORE_START).attr("refX", 18).attr("refY", 18).attr("markerWidth", 45).attr("markerHeight", 36).attr("orient", "auto").append("path").attr("stroke", conf2.stroke).attr("fill", "none").attr("d", "M0,18 Q 18,0 36,18 Q 18,36 0,18 M42,9 L42,27");
  elem.append("defs").append("marker").attr("id", ERMarkers.ONE_OR_MORE_END).attr("refX", 27).attr("refY", 18).attr("markerWidth", 45).attr("markerHeight", 36).attr("orient", "auto").append("path").attr("stroke", conf2.stroke).attr("fill", "none").attr("d", "M3,9 L3,27 M9,18 Q27,0 45,18 Q27,36 9,18");
  marker = elem.append("defs").append("marker").attr("id", ERMarkers.ZERO_OR_MORE_START).attr("refX", 18).attr("refY", 18).attr("markerWidth", 57).attr("markerHeight", 36).attr("orient", "auto");
  marker.append("circle").attr("stroke", conf2.stroke).attr("fill", "white").attr("cx", 48).attr("cy", 18).attr("r", 6);
  marker.append("path").attr("stroke", conf2.stroke).attr("fill", "none").attr("d", "M0,18 Q18,0 36,18 Q18,36 0,18");
  marker = elem.append("defs").append("marker").attr("id", ERMarkers.ZERO_OR_MORE_END).attr("refX", 39).attr("refY", 18).attr("markerWidth", 57).attr("markerHeight", 36).attr("orient", "auto");
  marker.append("circle").attr("stroke", conf2.stroke).attr("fill", "white").attr("cx", 9).attr("cy", 18).attr("r", 6);
  marker.append("path").attr("stroke", conf2.stroke).attr("fill", "none").attr("d", "M21,18 Q39,0 57,18 Q39,36 21,18");
  return;
};
const erMarkers = {
  ERMarkers,
  insertMarkers
};
const BAD_ID_CHARS_REGEXP = /[^\dA-Za-z](\W)*/g;
let conf = {};
let entityNameIds = /* @__PURE__ */ new Map();
const setConf = function(cnf) {
  const keys = Object.keys(cnf);
  for (const key of keys) {
    conf[key] = cnf[key];
  }
};
const drawAttributes = (groupNode, entityTextNode, attributes) => {
  const heightPadding = conf.entityPadding / 3;
  const widthPadding = conf.entityPadding / 3;
  const attrFontSize = conf.fontSize * 0.85;
  const labelBBox = entityTextNode.node().getBBox();
  const attributeNodes = [];
  let hasKeyType = false;
  let hasComment = false;
  let maxTypeWidth = 0;
  let maxNameWidth = 0;
  let maxKeyWidth = 0;
  let maxCommentWidth = 0;
  let cumulativeHeight = labelBBox.height + heightPadding * 2;
  let attrNum = 1;
  attributes.forEach((item) => {
    if (item.attributeKeyTypeList !== void 0 && item.attributeKeyTypeList.length > 0) {
      hasKeyType = true;
    }
    if (item.attributeComment !== void 0) {
      hasComment = true;
    }
  });
  attributes.forEach((item) => {
    const attrPrefix = `${entityTextNode.node().id}-attr-${attrNum}`;
    let nodeHeight = 0;
    const attributeType = (0,mermaid_04fb0060.v)(item.attributeType);
    const typeNode = groupNode.append("text").classed("er entityLabel", true).attr("id", `${attrPrefix}-type`).attr("x", 0).attr("y", 0).style("dominant-baseline", "middle").style("text-anchor", "left").style("font-family", (0,mermaid_04fb0060.c)().fontFamily).style("font-size", attrFontSize + "px").text(attributeType);
    const nameNode = groupNode.append("text").classed("er entityLabel", true).attr("id", `${attrPrefix}-name`).attr("x", 0).attr("y", 0).style("dominant-baseline", "middle").style("text-anchor", "left").style("font-family", (0,mermaid_04fb0060.c)().fontFamily).style("font-size", attrFontSize + "px").text(item.attributeName);
    const attributeNode = {};
    attributeNode.tn = typeNode;
    attributeNode.nn = nameNode;
    const typeBBox = typeNode.node().getBBox();
    const nameBBox = nameNode.node().getBBox();
    maxTypeWidth = Math.max(maxTypeWidth, typeBBox.width);
    maxNameWidth = Math.max(maxNameWidth, nameBBox.width);
    nodeHeight = Math.max(typeBBox.height, nameBBox.height);
    if (hasKeyType) {
      const keyTypeNodeText = item.attributeKeyTypeList !== void 0 ? item.attributeKeyTypeList.join(",") : "";
      const keyTypeNode = groupNode.append("text").classed("er entityLabel", true).attr("id", `${attrPrefix}-key`).attr("x", 0).attr("y", 0).style("dominant-baseline", "middle").style("text-anchor", "left").style("font-family", (0,mermaid_04fb0060.c)().fontFamily).style("font-size", attrFontSize + "px").text(keyTypeNodeText);
      attributeNode.kn = keyTypeNode;
      const keyTypeBBox = keyTypeNode.node().getBBox();
      maxKeyWidth = Math.max(maxKeyWidth, keyTypeBBox.width);
      nodeHeight = Math.max(nodeHeight, keyTypeBBox.height);
    }
    if (hasComment) {
      const commentNode = groupNode.append("text").classed("er entityLabel", true).attr("id", `${attrPrefix}-comment`).attr("x", 0).attr("y", 0).style("dominant-baseline", "middle").style("text-anchor", "left").style("font-family", (0,mermaid_04fb0060.c)().fontFamily).style("font-size", attrFontSize + "px").text(item.attributeComment || "");
      attributeNode.cn = commentNode;
      const commentNodeBBox = commentNode.node().getBBox();
      maxCommentWidth = Math.max(maxCommentWidth, commentNodeBBox.width);
      nodeHeight = Math.max(nodeHeight, commentNodeBBox.height);
    }
    attributeNode.height = nodeHeight;
    attributeNodes.push(attributeNode);
    cumulativeHeight += nodeHeight + heightPadding * 2;
    attrNum += 1;
  });
  let widthPaddingFactor = 4;
  if (hasKeyType) {
    widthPaddingFactor += 2;
  }
  if (hasComment) {
    widthPaddingFactor += 2;
  }
  const maxWidth = maxTypeWidth + maxNameWidth + maxKeyWidth + maxCommentWidth;
  const bBox = {
    width: Math.max(
      conf.minEntityWidth,
      Math.max(
        labelBBox.width + conf.entityPadding * 2,
        maxWidth + widthPadding * widthPaddingFactor
      )
    ),
    height: attributes.length > 0 ? cumulativeHeight : Math.max(conf.minEntityHeight, labelBBox.height + conf.entityPadding * 2)
  };
  if (attributes.length > 0) {
    const spareColumnWidth = Math.max(
      0,
      (bBox.width - maxWidth - widthPadding * widthPaddingFactor) / (widthPaddingFactor / 2)
    );
    entityTextNode.attr(
      "transform",
      "translate(" + bBox.width / 2 + "," + (heightPadding + labelBBox.height / 2) + ")"
    );
    let heightOffset = labelBBox.height + heightPadding * 2;
    let attribStyle = "attributeBoxOdd";
    attributeNodes.forEach((attributeNode) => {
      const alignY = heightOffset + heightPadding + attributeNode.height / 2;
      attributeNode.tn.attr("transform", "translate(" + widthPadding + "," + alignY + ")");
      const typeRect = groupNode.insert("rect", "#" + attributeNode.tn.node().id).classed(`er ${attribStyle}`, true).attr("x", 0).attr("y", heightOffset).attr("width", maxTypeWidth + widthPadding * 2 + spareColumnWidth).attr("height", attributeNode.height + heightPadding * 2);
      const nameXOffset = parseFloat(typeRect.attr("x")) + parseFloat(typeRect.attr("width"));
      attributeNode.nn.attr(
        "transform",
        "translate(" + (nameXOffset + widthPadding) + "," + alignY + ")"
      );
      const nameRect = groupNode.insert("rect", "#" + attributeNode.nn.node().id).classed(`er ${attribStyle}`, true).attr("x", nameXOffset).attr("y", heightOffset).attr("width", maxNameWidth + widthPadding * 2 + spareColumnWidth).attr("height", attributeNode.height + heightPadding * 2);
      let keyTypeAndCommentXOffset = parseFloat(nameRect.attr("x")) + parseFloat(nameRect.attr("width"));
      if (hasKeyType) {
        attributeNode.kn.attr(
          "transform",
          "translate(" + (keyTypeAndCommentXOffset + widthPadding) + "," + alignY + ")"
        );
        const keyTypeRect = groupNode.insert("rect", "#" + attributeNode.kn.node().id).classed(`er ${attribStyle}`, true).attr("x", keyTypeAndCommentXOffset).attr("y", heightOffset).attr("width", maxKeyWidth + widthPadding * 2 + spareColumnWidth).attr("height", attributeNode.height + heightPadding * 2);
        keyTypeAndCommentXOffset = parseFloat(keyTypeRect.attr("x")) + parseFloat(keyTypeRect.attr("width"));
      }
      if (hasComment) {
        attributeNode.cn.attr(
          "transform",
          "translate(" + (keyTypeAndCommentXOffset + widthPadding) + "," + alignY + ")"
        );
        groupNode.insert("rect", "#" + attributeNode.cn.node().id).classed(`er ${attribStyle}`, "true").attr("x", keyTypeAndCommentXOffset).attr("y", heightOffset).attr("width", maxCommentWidth + widthPadding * 2 + spareColumnWidth).attr("height", attributeNode.height + heightPadding * 2);
      }
      heightOffset += attributeNode.height + heightPadding * 2;
      attribStyle = attribStyle === "attributeBoxOdd" ? "attributeBoxEven" : "attributeBoxOdd";
    });
  } else {
    bBox.height = Math.max(conf.minEntityHeight, cumulativeHeight);
    entityTextNode.attr("transform", "translate(" + bBox.width / 2 + "," + bBox.height / 2 + ")");
  }
  return bBox;
};
const drawEntities = function(svgNode, entities2, graph) {
  const keys = Object.keys(entities2);
  let firstOne;
  keys.forEach(function(entityName) {
    const entityId = generateId(entityName, "entity");
    entityNameIds.set(entityName, entityId);
    const groupNode = svgNode.append("g").attr("id", entityId);
    firstOne = firstOne === void 0 ? entityId : firstOne;
    const textId = "text-" + entityId;
    const textNode = groupNode.append("text").classed("er entityLabel", true).attr("id", textId).attr("x", 0).attr("y", 0).style("dominant-baseline", "middle").style("text-anchor", "middle").style("font-family", (0,mermaid_04fb0060.c)().fontFamily).style("font-size", conf.fontSize + "px").text(entities2[entityName].alias ?? entityName);
    const { width: entityWidth, height: entityHeight } = drawAttributes(
      groupNode,
      textNode,
      entities2[entityName].attributes
    );
    const rectNode = groupNode.insert("rect", "#" + textId).classed("er entityBox", true).attr("x", 0).attr("y", 0).attr("width", entityWidth).attr("height", entityHeight);
    const rectBBox = rectNode.node().getBBox();
    graph.setNode(entityId, {
      width: rectBBox.width,
      height: rectBBox.height,
      shape: "rect",
      id: entityId
    });
  });
  return firstOne;
};
const adjustEntities = function(svgNode, graph) {
  graph.nodes().forEach(function(v) {
    if (v !== void 0 && graph.node(v) !== void 0) {
      svgNode.select("#" + v).attr(
        "transform",
        "translate(" + (graph.node(v).x - graph.node(v).width / 2) + "," + (graph.node(v).y - graph.node(v).height / 2) + " )"
      );
    }
  });
};
const getEdgeName = function(rel) {
  return (rel.entityA + rel.roleA + rel.entityB).replace(/\s/g, "");
};
const addRelationships = function(relationships2, g) {
  relationships2.forEach(function(r) {
    g.setEdge(
      entityNameIds.get(r.entityA),
      entityNameIds.get(r.entityB),
      { relationship: r },
      getEdgeName(r)
    );
  });
  return relationships2;
};
let relCnt = 0;
const drawRelationshipFromLayout = function(svg, rel, g, insert, diagObj) {
  relCnt++;
  const edge = g.edge(
    entityNameIds.get(rel.entityA),
    entityNameIds.get(rel.entityB),
    getEdgeName(rel)
  );
  const lineFunction = (0,src/* line */.jvg)().x(function(d) {
    return d.x;
  }).y(function(d) {
    return d.y;
  }).curve(src/* curveBasis */.$0Z);
  const svgPath = svg.insert("path", "#" + insert).classed("er relationshipLine", true).attr("d", lineFunction(edge.points)).style("stroke", conf.stroke).style("fill", "none");
  if (rel.relSpec.relType === diagObj.db.Identification.NON_IDENTIFYING) {
    svgPath.attr("stroke-dasharray", "8,8");
  }
  let url = "";
  if (conf.arrowMarkerAbsolute) {
    url = window.location.protocol + "//" + window.location.host + window.location.pathname + window.location.search;
    url = url.replace(/\(/g, "\\(");
    url = url.replace(/\)/g, "\\)");
  }
  switch (rel.relSpec.cardA) {
    case diagObj.db.Cardinality.ZERO_OR_ONE:
      svgPath.attr("marker-end", "url(" + url + "#" + erMarkers.ERMarkers.ZERO_OR_ONE_END + ")");
      break;
    case diagObj.db.Cardinality.ZERO_OR_MORE:
      svgPath.attr("marker-end", "url(" + url + "#" + erMarkers.ERMarkers.ZERO_OR_MORE_END + ")");
      break;
    case diagObj.db.Cardinality.ONE_OR_MORE:
      svgPath.attr("marker-end", "url(" + url + "#" + erMarkers.ERMarkers.ONE_OR_MORE_END + ")");
      break;
    case diagObj.db.Cardinality.ONLY_ONE:
      svgPath.attr("marker-end", "url(" + url + "#" + erMarkers.ERMarkers.ONLY_ONE_END + ")");
      break;
    case diagObj.db.Cardinality.MD_PARENT:
      svgPath.attr("marker-end", "url(" + url + "#" + erMarkers.ERMarkers.MD_PARENT_END + ")");
      break;
  }
  switch (rel.relSpec.cardB) {
    case diagObj.db.Cardinality.ZERO_OR_ONE:
      svgPath.attr(
        "marker-start",
        "url(" + url + "#" + erMarkers.ERMarkers.ZERO_OR_ONE_START + ")"
      );
      break;
    case diagObj.db.Cardinality.ZERO_OR_MORE:
      svgPath.attr(
        "marker-start",
        "url(" + url + "#" + erMarkers.ERMarkers.ZERO_OR_MORE_START + ")"
      );
      break;
    case diagObj.db.Cardinality.ONE_OR_MORE:
      svgPath.attr(
        "marker-start",
        "url(" + url + "#" + erMarkers.ERMarkers.ONE_OR_MORE_START + ")"
      );
      break;
    case diagObj.db.Cardinality.ONLY_ONE:
      svgPath.attr("marker-start", "url(" + url + "#" + erMarkers.ERMarkers.ONLY_ONE_START + ")");
      break;
    case diagObj.db.Cardinality.MD_PARENT:
      svgPath.attr("marker-start", "url(" + url + "#" + erMarkers.ERMarkers.MD_PARENT_START + ")");
      break;
  }
  const len = svgPath.node().getTotalLength();
  const labelPoint = svgPath.node().getPointAtLength(len * 0.5);
  const labelId = "rel" + relCnt;
  const labelNode = svg.append("text").classed("er relationshipLabel", true).attr("id", labelId).attr("x", labelPoint.x).attr("y", labelPoint.y).style("text-anchor", "middle").style("dominant-baseline", "middle").style("font-family", (0,mermaid_04fb0060.c)().fontFamily).style("font-size", conf.fontSize + "px").text(rel.roleA);
  const labelBBox = labelNode.node().getBBox();
  svg.insert("rect", "#" + labelId).classed("er relationshipLabelBox", true).attr("x", labelPoint.x - labelBBox.width / 2).attr("y", labelPoint.y - labelBBox.height / 2).attr("width", labelBBox.width).attr("height", labelBBox.height);
};
const draw = function(text, id, _version, diagObj) {
  conf = (0,mermaid_04fb0060.c)().er;
  mermaid_04fb0060.l.info("Drawing ER diagram");
  const securityLevel = (0,mermaid_04fb0060.c)().securityLevel;
  let sandboxElement;
  if (securityLevel === "sandbox") {
    sandboxElement = (0,src/* select */.Ys)("#i" + id);
  }
  const root = securityLevel === "sandbox" ? (0,src/* select */.Ys)(sandboxElement.nodes()[0].contentDocument.body) : (0,src/* select */.Ys)("body");
  const svg = root.select(`[id='${id}']`);
  erMarkers.insertMarkers(svg, conf);
  let g;
  g = new graphlib/* Graph */.k({
    multigraph: true,
    directed: true,
    compound: false
  }).setGraph({
    rankdir: conf.layoutDirection,
    marginx: 20,
    marginy: 20,
    nodesep: 100,
    edgesep: 100,
    ranksep: 100
  }).setDefaultEdgeLabel(function() {
    return {};
  });
  const firstEntity = drawEntities(svg, diagObj.db.getEntities(), g);
  const relationships2 = addRelationships(diagObj.db.getRelationships(), g);
  (0,dagre/* layout */.bK)(g);
  adjustEntities(svg, g);
  relationships2.forEach(function(rel) {
    drawRelationshipFromLayout(svg, rel, g, firstEntity, diagObj);
  });
  const padding = conf.diagramPadding;
  mermaid_04fb0060.u.insertTitle(svg, "entityTitleText", conf.titleTopMargin, diagObj.db.getDiagramTitle());
  const svgBounds = svg.node().getBBox();
  const width = svgBounds.width + padding * 2;
  const height = svgBounds.height + padding * 2;
  (0,mermaid_04fb0060.i)(svg, height, width, conf.useMaxWidth);
  svg.attr("viewBox", `${svgBounds.x - padding} ${svgBounds.y - padding} ${width} ${height}`);
};
const MERMAID_ERDIAGRAM_UUID = "28e9f9db-3c8d-5aa5-9faf-44286ae5937c";
function generateId(str = "", prefix = "") {
  const simplifiedStr = str.replace(BAD_ID_CHARS_REGEXP, "");
  return `${strWithHyphen(prefix)}${strWithHyphen(simplifiedStr)}${esm_browser_v5(
    str,
    MERMAID_ERDIAGRAM_UUID
  )}`;
}
function strWithHyphen(str = "") {
  return str.length > 0 ? `${str}-` : "";
}
const erRenderer = {
  setConf,
  draw
};
const getStyles = (options) => `
  .entityBox {
    fill: ${options.mainBkg};
    stroke: ${options.nodeBorder};
  }

  .attributeBoxOdd {
    fill: ${options.attributeBackgroundColorOdd};
    stroke: ${options.nodeBorder};
  }

  .attributeBoxEven {
    fill:  ${options.attributeBackgroundColorEven};
    stroke: ${options.nodeBorder};
  }

  .relationshipLabelBox {
    fill: ${options.tertiaryColor};
    opacity: 0.7;
    background-color: ${options.tertiaryColor};
      rect {
        opacity: 0.5;
      }
  }

    .relationshipLine {
      stroke: ${options.lineColor};
    }

  .entityTitleText {
    text-anchor: middle;
    font-size: 18px;
    fill: ${options.textColor};
  }    
  #MD_PARENT_START {
    fill: #f5f5f5 !important;
    stroke: ${options.lineColor} !important;
    stroke-width: 1;
  }
  #MD_PARENT_END {
    fill: #f5f5f5 !important;
    stroke: ${options.lineColor} !important;
    stroke-width: 1;
  }
  
`;
const erStyles = getStyles;
const diagram = {
  parser: erParser,
  db: erDb,
  renderer: erRenderer,
  styles: erStyles
};



/***/ })

}]);
//# sourceMappingURL=2692.aa472750b0685d9dc0b2.js.map?v=aa472750b0685d9dc0b2
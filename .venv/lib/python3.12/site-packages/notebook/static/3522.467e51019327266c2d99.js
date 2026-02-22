"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[3522],{

/***/ 40448:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   q: () => (/* binding */ getDiagramElement)
/* harmony export */ });
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(74999);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(35321);


// src/rendering-util/insertElementsForSize.js

var getDiagramElement = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((id, securityLevel) => {
  let sandboxElement;
  if (securityLevel === "sandbox") {
    sandboxElement = (0,d3__WEBPACK_IMPORTED_MODULE_1__/* .select */ .Ys)("#i" + id);
  }
  const root = securityLevel === "sandbox" ? (0,d3__WEBPACK_IMPORTED_MODULE_1__/* .select */ .Ys)(sandboxElement.nodes()[0].contentDocument.body) : (0,d3__WEBPACK_IMPORTED_MODULE_1__/* .select */ .Ys)("body");
  const svg = root.select(`[id="${id}"]`);
  return svg;
}, "getDiagramElement");




/***/ }),

/***/ 51755:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   j: () => (/* binding */ setupViewPortForSVG)
/* harmony export */ });
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(74999);



// src/rendering-util/setupViewPortForSVG.ts
var setupViewPortForSVG = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((svg, padding, cssDiagram, useMaxWidth) => {
  svg.attr("class", cssDiagram);
  const { width, height, x, y } = calculateDimensionsWithPadding(svg, padding);
  (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .configureSvgSize */ .v2)(svg, height, width, useMaxWidth);
  const viewBox = createViewBox(x, y, width, height, padding);
  svg.attr("viewBox", viewBox);
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .log */ .cM.debug(`viewBox configured: ${viewBox} with padding: ${padding}`);
}, "setupViewPortForSVG");
var calculateDimensionsWithPadding = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((svg, padding) => {
  const bounds = svg.node()?.getBBox() || { width: 0, height: 0, x: 0, y: 0 };
  return {
    width: bounds.width + padding * 2,
    height: bounds.height + padding * 2,
    x: bounds.x,
    y: bounds.y
  };
}, "calculateDimensionsWithPadding");
var createViewBox = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((x, y, width, height, padding) => {
  return `${x - padding} ${y - padding} ${width} ${height}`;
}, "createViewBox");




/***/ }),

/***/ 3522:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {


// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  diagram: () => (/* binding */ diagram)
});

// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-55IACEB6.mjs
var chunk_55IACEB6 = __webpack_require__(40448);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-QN33PNHL.mjs
var chunk_QN33PNHL = __webpack_require__(51755);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-N4CR4FBY.mjs
var chunk_N4CR4FBY = __webpack_require__(54160);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-QXUST7PY.mjs
var chunk_QXUST7PY = __webpack_require__(55195);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-HN2XXSSU.mjs
var chunk_HN2XXSSU = __webpack_require__(67894);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-JZLCHNYA.mjs
var chunk_JZLCHNYA = __webpack_require__(26212);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-CVBHYZKI.mjs
var chunk_CVBHYZKI = __webpack_require__(22462);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-ATLVNIR6.mjs
var chunk_ATLVNIR6 = __webpack_require__(37756);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-JA3XYJ7Z.mjs + 12 modules
var chunk_JA3XYJ7Z = __webpack_require__(2111);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-S3R3BYOJ.mjs
var chunk_S3R3BYOJ = __webpack_require__(17175);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-ABZYJK2D.mjs + 3 modules
var chunk_ABZYJK2D = __webpack_require__(61805);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-AGHRB4JF.mjs
var chunk_AGHRB4JF = __webpack_require__(74999);
;// CONCATENATED MODULE: ../node_modules/mermaid/node_modules/uuid/dist/esm-browser/native.js
const randomUUID = typeof crypto !== 'undefined' && crypto.randomUUID && crypto.randomUUID.bind(crypto);
/* harmony default export */ const esm_browser_native = ({ randomUUID });

;// CONCATENATED MODULE: ../node_modules/mermaid/node_modules/uuid/dist/esm-browser/rng.js
let getRandomValues;
const rnds8 = new Uint8Array(16);
function rng() {
    if (!getRandomValues) {
        if (typeof crypto === 'undefined' || !crypto.getRandomValues) {
            throw new Error('crypto.getRandomValues() not supported. See https://github.com/uuidjs/uuid#getrandomvalues-not-supported');
        }
        getRandomValues = crypto.getRandomValues.bind(crypto);
    }
    return getRandomValues(rnds8);
}

;// CONCATENATED MODULE: ../node_modules/mermaid/node_modules/uuid/dist/esm-browser/stringify.js

const byteToHex = [];
for (let i = 0; i < 256; ++i) {
    byteToHex.push((i + 0x100).toString(16).slice(1));
}
function unsafeStringify(arr, offset = 0) {
    return (byteToHex[arr[offset + 0]] +
        byteToHex[arr[offset + 1]] +
        byteToHex[arr[offset + 2]] +
        byteToHex[arr[offset + 3]] +
        '-' +
        byteToHex[arr[offset + 4]] +
        byteToHex[arr[offset + 5]] +
        '-' +
        byteToHex[arr[offset + 6]] +
        byteToHex[arr[offset + 7]] +
        '-' +
        byteToHex[arr[offset + 8]] +
        byteToHex[arr[offset + 9]] +
        '-' +
        byteToHex[arr[offset + 10]] +
        byteToHex[arr[offset + 11]] +
        byteToHex[arr[offset + 12]] +
        byteToHex[arr[offset + 13]] +
        byteToHex[arr[offset + 14]] +
        byteToHex[arr[offset + 15]]).toLowerCase();
}
function stringify(arr, offset = 0) {
    const uuid = unsafeStringify(arr, offset);
    if (!validate(uuid)) {
        throw TypeError('Stringified UUID is invalid');
    }
    return uuid;
}
/* harmony default export */ const esm_browser_stringify = ((/* unused pure expression or super */ null && (stringify)));

;// CONCATENATED MODULE: ../node_modules/mermaid/node_modules/uuid/dist/esm-browser/v4.js



function v4(options, buf, offset) {
    if (esm_browser_native.randomUUID && !buf && !options) {
        return esm_browser_native.randomUUID();
    }
    options = options || {};
    const rnds = options.random ?? options.rng?.() ?? rng();
    if (rnds.length < 16) {
        throw new Error('Random bytes length must be >= 16');
    }
    rnds[6] = (rnds[6] & 0x0f) | 0x40;
    rnds[8] = (rnds[8] & 0x3f) | 0x80;
    if (buf) {
        offset = offset || 0;
        if (offset < 0 || offset + 16 > buf.length) {
            throw new RangeError(`UUID byte range ${offset}:${offset + 15} is out of buffer bounds`);
        }
        for (let i = 0; i < 16; ++i) {
            buf[offset + i] = rnds[i];
        }
        return buf;
    }
    return unsafeStringify(rnds);
}
/* harmony default export */ const esm_browser_v4 = (v4);

// EXTERNAL MODULE: ../node_modules/khroma/dist/methods/is_dark.js + 2 modules
var is_dark = __webpack_require__(28186);
// EXTERNAL MODULE: ../node_modules/khroma/dist/methods/lighten.js
var lighten = __webpack_require__(28482);
// EXTERNAL MODULE: ../node_modules/khroma/dist/methods/darken.js
var darken = __webpack_require__(57838);
;// CONCATENATED MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/mindmap-definition-VGOIOE7T.mjs













// src/diagrams/mindmap/parser/mindmap.jison
var parser = (function() {
  var o = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(k, v, o2, l) {
    for (o2 = o2 || {}, l = k.length; l--; o2[k[l]] = v) ;
    return o2;
  }, "o"), $V0 = [1, 4], $V1 = [1, 13], $V2 = [1, 12], $V3 = [1, 15], $V4 = [1, 16], $V5 = [1, 20], $V6 = [1, 19], $V7 = [6, 7, 8], $V8 = [1, 26], $V9 = [1, 24], $Va = [1, 25], $Vb = [6, 7, 11], $Vc = [1, 6, 13, 15, 16, 19, 22], $Vd = [1, 33], $Ve = [1, 34], $Vf = [1, 6, 7, 11, 13, 15, 16, 19, 22];
  var parser2 = {
    trace: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function trace() {
    }, "trace"),
    yy: {},
    symbols_: { "error": 2, "start": 3, "mindMap": 4, "spaceLines": 5, "SPACELINE": 6, "NL": 7, "MINDMAP": 8, "document": 9, "stop": 10, "EOF": 11, "statement": 12, "SPACELIST": 13, "node": 14, "ICON": 15, "CLASS": 16, "nodeWithId": 17, "nodeWithoutId": 18, "NODE_DSTART": 19, "NODE_DESCR": 20, "NODE_DEND": 21, "NODE_ID": 22, "$accept": 0, "$end": 1 },
    terminals_: { 2: "error", 6: "SPACELINE", 7: "NL", 8: "MINDMAP", 11: "EOF", 13: "SPACELIST", 15: "ICON", 16: "CLASS", 19: "NODE_DSTART", 20: "NODE_DESCR", 21: "NODE_DEND", 22: "NODE_ID" },
    productions_: [0, [3, 1], [3, 2], [5, 1], [5, 2], [5, 2], [4, 2], [4, 3], [10, 1], [10, 1], [10, 1], [10, 2], [10, 2], [9, 3], [9, 2], [12, 2], [12, 2], [12, 2], [12, 1], [12, 1], [12, 1], [12, 1], [12, 1], [14, 1], [14, 1], [18, 3], [17, 1], [17, 4]],
    performAction: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function anonymous(yytext, yyleng, yylineno, yy, yystate, $$, _$) {
      var $0 = $$.length - 1;
      switch (yystate) {
        case 6:
        case 7:
          return yy;
          break;
        case 8:
          yy.getLogger().trace("Stop NL ");
          break;
        case 9:
          yy.getLogger().trace("Stop EOF ");
          break;
        case 11:
          yy.getLogger().trace("Stop NL2 ");
          break;
        case 12:
          yy.getLogger().trace("Stop EOF2 ");
          break;
        case 15:
          yy.getLogger().info("Node: ", $$[$0].id);
          yy.addNode($$[$0 - 1].length, $$[$0].id, $$[$0].descr, $$[$0].type);
          break;
        case 16:
          yy.getLogger().trace("Icon: ", $$[$0]);
          yy.decorateNode({ icon: $$[$0] });
          break;
        case 17:
        case 21:
          yy.decorateNode({ class: $$[$0] });
          break;
        case 18:
          yy.getLogger().trace("SPACELIST");
          break;
        case 19:
          yy.getLogger().trace("Node: ", $$[$0].id);
          yy.addNode(0, $$[$0].id, $$[$0].descr, $$[$0].type);
          break;
        case 20:
          yy.decorateNode({ icon: $$[$0] });
          break;
        case 25:
          yy.getLogger().trace("node found ..", $$[$0 - 2]);
          this.$ = { id: $$[$0 - 1], descr: $$[$0 - 1], type: yy.getType($$[$0 - 2], $$[$0]) };
          break;
        case 26:
          this.$ = { id: $$[$0], descr: $$[$0], type: yy.nodeType.DEFAULT };
          break;
        case 27:
          yy.getLogger().trace("node found ..", $$[$0 - 3]);
          this.$ = { id: $$[$0 - 3], descr: $$[$0 - 1], type: yy.getType($$[$0 - 2], $$[$0]) };
          break;
      }
    }, "anonymous"),
    table: [{ 3: 1, 4: 2, 5: 3, 6: [1, 5], 8: $V0 }, { 1: [3] }, { 1: [2, 1] }, { 4: 6, 6: [1, 7], 7: [1, 8], 8: $V0 }, { 6: $V1, 7: [1, 10], 9: 9, 12: 11, 13: $V2, 14: 14, 15: $V3, 16: $V4, 17: 17, 18: 18, 19: $V5, 22: $V6 }, o($V7, [2, 3]), { 1: [2, 2] }, o($V7, [2, 4]), o($V7, [2, 5]), { 1: [2, 6], 6: $V1, 12: 21, 13: $V2, 14: 14, 15: $V3, 16: $V4, 17: 17, 18: 18, 19: $V5, 22: $V6 }, { 6: $V1, 9: 22, 12: 11, 13: $V2, 14: 14, 15: $V3, 16: $V4, 17: 17, 18: 18, 19: $V5, 22: $V6 }, { 6: $V8, 7: $V9, 10: 23, 11: $Va }, o($Vb, [2, 22], { 17: 17, 18: 18, 14: 27, 15: [1, 28], 16: [1, 29], 19: $V5, 22: $V6 }), o($Vb, [2, 18]), o($Vb, [2, 19]), o($Vb, [2, 20]), o($Vb, [2, 21]), o($Vb, [2, 23]), o($Vb, [2, 24]), o($Vb, [2, 26], { 19: [1, 30] }), { 20: [1, 31] }, { 6: $V8, 7: $V9, 10: 32, 11: $Va }, { 1: [2, 7], 6: $V1, 12: 21, 13: $V2, 14: 14, 15: $V3, 16: $V4, 17: 17, 18: 18, 19: $V5, 22: $V6 }, o($Vc, [2, 14], { 7: $Vd, 11: $Ve }), o($Vf, [2, 8]), o($Vf, [2, 9]), o($Vf, [2, 10]), o($Vb, [2, 15]), o($Vb, [2, 16]), o($Vb, [2, 17]), { 20: [1, 35] }, { 21: [1, 36] }, o($Vc, [2, 13], { 7: $Vd, 11: $Ve }), o($Vf, [2, 11]), o($Vf, [2, 12]), { 21: [1, 37] }, o($Vb, [2, 25]), o($Vb, [2, 27])],
    defaultActions: { 2: [2, 1], 6: [2, 2] },
    parseError: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function parseError(str, hash) {
      if (hash.recoverable) {
        this.trace(str);
      } else {
        var error = new Error(str);
        error.hash = hash;
        throw error;
      }
    }, "parseError"),
    parse: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function parse(input) {
      var self = this, stack = [0], tstack = [], vstack = [null], lstack = [], table = this.table, yytext = "", yylineno = 0, yyleng = 0, recovering = 0, TERROR = 2, EOF = 1;
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
      function popStack(n) {
        stack.length = stack.length - 2 * n;
        vstack.length = vstack.length - n;
        lstack.length = lstack.length - n;
      }
      (0,chunk_AGHRB4JF/* __name */.eW)(popStack, "popStack");
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
      (0,chunk_AGHRB4JF/* __name */.eW)(lex, "lex");
      var symbol, preErrorSymbol, state, action, a, r, yyval = {}, p, len, newState, expected;
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
            if (!preErrorSymbol) {
              yyleng = lexer2.yyleng;
              yytext = lexer2.yytext;
              yylineno = lexer2.yylineno;
              yyloc = lexer2.yylloc;
              if (recovering > 0) {
                recovering--;
              }
            } else {
              symbol = preErrorSymbol;
              preErrorSymbol = null;
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
    }, "parse")
  };
  var lexer = /* @__PURE__ */ (function() {
    var lexer2 = {
      EOF: 1,
      parseError: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function parseError(str, hash) {
        if (this.yy.parser) {
          this.yy.parser.parseError(str, hash);
        } else {
          throw new Error(str);
        }
      }, "parseError"),
      // resets the lexer, sets new input
      setInput: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(input, yy) {
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
      }, "setInput"),
      // consumes and returns one char from the input
      input: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
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
      }, "input"),
      // unshifts one char (or a string) into the input
      unput: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(ch) {
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
      }, "unput"),
      // When called from action, caches matched text and appends it on next action
      more: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
        this._more = true;
        return this;
      }, "more"),
      // When called from action, signals the lexer that this rule fails to match the input, so the next matching rule (regex) should be tested instead.
      reject: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
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
      }, "reject"),
      // retain first n characters of the match
      less: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(n) {
        this.unput(this.match.slice(n));
      }, "less"),
      // displays already matched input, i.e. for error messages
      pastInput: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
        var past = this.matched.substr(0, this.matched.length - this.match.length);
        return (past.length > 20 ? "..." : "") + past.substr(-20).replace(/\n/g, "");
      }, "pastInput"),
      // displays upcoming input, i.e. for error messages
      upcomingInput: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
        var next = this.match;
        if (next.length < 20) {
          next += this._input.substr(0, 20 - next.length);
        }
        return (next.substr(0, 20) + (next.length > 20 ? "..." : "")).replace(/\n/g, "");
      }, "upcomingInput"),
      // displays the character position where the lexing error occurred, i.e. for error messages
      showPosition: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
        var pre = this.pastInput();
        var c = new Array(pre.length + 1).join("-");
        return pre + this.upcomingInput() + "\n" + c + "^";
      }, "showPosition"),
      // test the lexed token: return FALSE when not a match, otherwise return token
      test_match: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(match, indexed_rule) {
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
      }, "test_match"),
      // return next match in input
      next: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
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
      }, "next"),
      // return next match that has a token
      lex: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function lex() {
        var r = this.next();
        if (r) {
          return r;
        } else {
          return this.lex();
        }
      }, "lex"),
      // activates a new lexer condition state (pushes the new lexer condition state onto the condition stack)
      begin: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function begin(condition) {
        this.conditionStack.push(condition);
      }, "begin"),
      // pop the previously active lexer condition state off the condition stack
      popState: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function popState() {
        var n = this.conditionStack.length - 1;
        if (n > 0) {
          return this.conditionStack.pop();
        } else {
          return this.conditionStack[0];
        }
      }, "popState"),
      // produce the lexer rule set which is active for the currently active lexer condition state
      _currentRules: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function _currentRules() {
        if (this.conditionStack.length && this.conditionStack[this.conditionStack.length - 1]) {
          return this.conditions[this.conditionStack[this.conditionStack.length - 1]].rules;
        } else {
          return this.conditions["INITIAL"].rules;
        }
      }, "_currentRules"),
      // return the currently active lexer condition state; when an index argument is provided it produces the N-th previous condition state, if available
      topState: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function topState(n) {
        n = this.conditionStack.length - 1 - Math.abs(n || 0);
        if (n >= 0) {
          return this.conditionStack[n];
        } else {
          return "INITIAL";
        }
      }, "topState"),
      // alias for begin(condition)
      pushState: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function pushState(condition) {
        this.begin(condition);
      }, "pushState"),
      // return the number of states currently on the stack
      stateStackSize: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function stateStackSize() {
        return this.conditionStack.length;
      }, "stateStackSize"),
      options: { "case-insensitive": true },
      performAction: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function anonymous(yy, yy_, $avoiding_name_collisions, YY_START) {
        var YYSTATE = YY_START;
        switch ($avoiding_name_collisions) {
          case 0:
            yy.getLogger().trace("Found comment", yy_.yytext);
            return 6;
            break;
          case 1:
            return 8;
            break;
          case 2:
            this.begin("CLASS");
            break;
          case 3:
            this.popState();
            return 16;
            break;
          case 4:
            this.popState();
            break;
          case 5:
            yy.getLogger().trace("Begin icon");
            this.begin("ICON");
            break;
          case 6:
            yy.getLogger().trace("SPACELINE");
            return 6;
            break;
          case 7:
            return 7;
            break;
          case 8:
            return 15;
            break;
          case 9:
            yy.getLogger().trace("end icon");
            this.popState();
            break;
          case 10:
            yy.getLogger().trace("Exploding node");
            this.begin("NODE");
            return 19;
            break;
          case 11:
            yy.getLogger().trace("Cloud");
            this.begin("NODE");
            return 19;
            break;
          case 12:
            yy.getLogger().trace("Explosion Bang");
            this.begin("NODE");
            return 19;
            break;
          case 13:
            yy.getLogger().trace("Cloud Bang");
            this.begin("NODE");
            return 19;
            break;
          case 14:
            this.begin("NODE");
            return 19;
            break;
          case 15:
            this.begin("NODE");
            return 19;
            break;
          case 16:
            this.begin("NODE");
            return 19;
            break;
          case 17:
            this.begin("NODE");
            return 19;
            break;
          case 18:
            return 13;
            break;
          case 19:
            return 22;
            break;
          case 20:
            return 11;
            break;
          case 21:
            this.begin("NSTR2");
            break;
          case 22:
            return "NODE_DESCR";
            break;
          case 23:
            this.popState();
            break;
          case 24:
            yy.getLogger().trace("Starting NSTR");
            this.begin("NSTR");
            break;
          case 25:
            yy.getLogger().trace("description:", yy_.yytext);
            return "NODE_DESCR";
            break;
          case 26:
            this.popState();
            break;
          case 27:
            this.popState();
            yy.getLogger().trace("node end ))");
            return "NODE_DEND";
            break;
          case 28:
            this.popState();
            yy.getLogger().trace("node end )");
            return "NODE_DEND";
            break;
          case 29:
            this.popState();
            yy.getLogger().trace("node end ...", yy_.yytext);
            return "NODE_DEND";
            break;
          case 30:
            this.popState();
            yy.getLogger().trace("node end ((");
            return "NODE_DEND";
            break;
          case 31:
            this.popState();
            yy.getLogger().trace("node end (-");
            return "NODE_DEND";
            break;
          case 32:
            this.popState();
            yy.getLogger().trace("node end (-");
            return "NODE_DEND";
            break;
          case 33:
            this.popState();
            yy.getLogger().trace("node end ((");
            return "NODE_DEND";
            break;
          case 34:
            this.popState();
            yy.getLogger().trace("node end ((");
            return "NODE_DEND";
            break;
          case 35:
            yy.getLogger().trace("Long description:", yy_.yytext);
            return 20;
            break;
          case 36:
            yy.getLogger().trace("Long description:", yy_.yytext);
            return 20;
            break;
        }
      }, "anonymous"),
      rules: [/^(?:\s*%%.*)/i, /^(?:mindmap\b)/i, /^(?::::)/i, /^(?:.+)/i, /^(?:\n)/i, /^(?:::icon\()/i, /^(?:[\s]+[\n])/i, /^(?:[\n]+)/i, /^(?:[^\)]+)/i, /^(?:\))/i, /^(?:-\))/i, /^(?:\(-)/i, /^(?:\)\))/i, /^(?:\))/i, /^(?:\(\()/i, /^(?:\{\{)/i, /^(?:\()/i, /^(?:\[)/i, /^(?:[\s]+)/i, /^(?:[^\(\[\n\)\{\}]+)/i, /^(?:$)/i, /^(?:["][`])/i, /^(?:[^`"]+)/i, /^(?:[`]["])/i, /^(?:["])/i, /^(?:[^"]+)/i, /^(?:["])/i, /^(?:[\)]\))/i, /^(?:[\)])/i, /^(?:[\]])/i, /^(?:\}\})/i, /^(?:\(-)/i, /^(?:-\))/i, /^(?:\(\()/i, /^(?:\()/i, /^(?:[^\)\]\(\}]+)/i, /^(?:.+(?!\(\())/i],
      conditions: { "CLASS": { "rules": [3, 4], "inclusive": false }, "ICON": { "rules": [8, 9], "inclusive": false }, "NSTR2": { "rules": [22, 23], "inclusive": false }, "NSTR": { "rules": [25, 26], "inclusive": false }, "NODE": { "rules": [21, 24, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36], "inclusive": false }, "INITIAL": { "rules": [0, 1, 2, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], "inclusive": true } }
    };
    return lexer2;
  })();
  parser2.lexer = lexer;
  function Parser() {
    this.yy = {};
  }
  (0,chunk_AGHRB4JF/* __name */.eW)(Parser, "Parser");
  Parser.prototype = parser2;
  parser2.Parser = Parser;
  return new Parser();
})();
parser.parser = parser;
var mindmap_default = parser;

// src/diagrams/mindmap/mindmapDb.ts

var nodeType = {
  DEFAULT: 0,
  NO_BORDER: 0,
  ROUNDED_RECT: 1,
  RECT: 2,
  CIRCLE: 3,
  CLOUD: 4,
  BANG: 5,
  HEXAGON: 6
};
var MindmapDB = class {
  constructor() {
    this.nodes = [];
    this.count = 0;
    this.elements = {};
    this.getLogger = this.getLogger.bind(this);
    this.nodeType = nodeType;
    this.clear();
    this.getType = this.getType.bind(this);
    this.getElementById = this.getElementById.bind(this);
    this.getParent = this.getParent.bind(this);
    this.getMindmap = this.getMindmap.bind(this);
    this.addNode = this.addNode.bind(this);
    this.decorateNode = this.decorateNode.bind(this);
  }
  static {
    (0,chunk_AGHRB4JF/* __name */.eW)(this, "MindmapDB");
  }
  clear() {
    this.nodes = [];
    this.count = 0;
    this.elements = {};
    this.baseLevel = void 0;
  }
  getParent(level) {
    for (let i = this.nodes.length - 1; i >= 0; i--) {
      if (this.nodes[i].level < level) {
        return this.nodes[i];
      }
    }
    return null;
  }
  getMindmap() {
    return this.nodes.length > 0 ? this.nodes[0] : null;
  }
  addNode(level, id, descr, type) {
    chunk_AGHRB4JF/* log */.cM.info("addNode", level, id, descr, type);
    let isRoot = false;
    if (this.nodes.length === 0) {
      this.baseLevel = level;
      level = 0;
      isRoot = true;
    } else if (this.baseLevel !== void 0) {
      level = level - this.baseLevel;
      isRoot = false;
    }
    const conf = (0,chunk_ABZYJK2D/* getConfig2 */.nV)();
    let padding = conf.mindmap?.padding ?? chunk_ABZYJK2D/* defaultConfig_default */.vZ.mindmap.padding;
    switch (type) {
      case this.nodeType.ROUNDED_RECT:
      case this.nodeType.RECT:
      case this.nodeType.HEXAGON:
        padding *= 2;
        break;
    }
    const node = {
      id: this.count++,
      nodeId: (0,chunk_ABZYJK2D/* sanitizeText */.oO)(id, conf),
      level,
      descr: (0,chunk_ABZYJK2D/* sanitizeText */.oO)(descr, conf),
      type,
      children: [],
      width: conf.mindmap?.maxNodeWidth ?? chunk_ABZYJK2D/* defaultConfig_default */.vZ.mindmap.maxNodeWidth,
      padding,
      isRoot
    };
    const parent = this.getParent(level);
    if (parent) {
      parent.children.push(node);
      this.nodes.push(node);
    } else {
      if (isRoot) {
        this.nodes.push(node);
      } else {
        throw new Error(
          `There can be only one root. No parent could be found for ("${node.descr}")`
        );
      }
    }
  }
  getType(startStr, endStr) {
    chunk_AGHRB4JF/* log */.cM.debug("In get type", startStr, endStr);
    switch (startStr) {
      case "[":
        return this.nodeType.RECT;
      case "(":
        return endStr === ")" ? this.nodeType.ROUNDED_RECT : this.nodeType.CLOUD;
      case "((":
        return this.nodeType.CIRCLE;
      case ")":
        return this.nodeType.CLOUD;
      case "))":
        return this.nodeType.BANG;
      case "{{":
        return this.nodeType.HEXAGON;
      default:
        return this.nodeType.DEFAULT;
    }
  }
  setElementForId(id, element) {
    this.elements[id] = element;
  }
  getElementById(id) {
    return this.elements[id];
  }
  decorateNode(decoration) {
    if (!decoration) {
      return;
    }
    const config = (0,chunk_ABZYJK2D/* getConfig2 */.nV)();
    const node = this.nodes[this.nodes.length - 1];
    if (decoration.icon) {
      node.icon = (0,chunk_ABZYJK2D/* sanitizeText */.oO)(decoration.icon, config);
    }
    if (decoration.class) {
      node.class = (0,chunk_ABZYJK2D/* sanitizeText */.oO)(decoration.class, config);
    }
  }
  type2Str(type) {
    switch (type) {
      case this.nodeType.DEFAULT:
        return "no-border";
      case this.nodeType.RECT:
        return "rect";
      case this.nodeType.ROUNDED_RECT:
        return "rounded-rect";
      case this.nodeType.CIRCLE:
        return "circle";
      case this.nodeType.CLOUD:
        return "cloud";
      case this.nodeType.BANG:
        return "bang";
      case this.nodeType.HEXAGON:
        return "hexgon";
      // cspell: disable-line
      default:
        return "no-border";
    }
  }
  /**
   * Assign section numbers to nodes based on their position relative to root
   * @param node - The mindmap node to process
   * @param sectionNumber - The section number to assign (undefined for root)
   */
  assignSections(node, sectionNumber) {
    if (node.level === 0) {
      node.section = void 0;
    } else {
      node.section = sectionNumber;
    }
    if (node.children) {
      for (const [index, child] of node.children.entries()) {
        const childSectionNumber = node.level === 0 ? index : sectionNumber;
        this.assignSections(child, childSectionNumber);
      }
    }
  }
  /**
   * Convert mindmap tree structure to flat array of nodes
   * @param node - The mindmap node to process
   * @param processedNodes - Array to collect processed nodes
   */
  flattenNodes(node, processedNodes) {
    const cssClasses = ["mindmap-node"];
    if (node.isRoot === true) {
      cssClasses.push("section-root", "section--1");
    } else if (node.section !== void 0) {
      cssClasses.push(`section-${node.section}`);
    }
    if (node.class) {
      cssClasses.push(node.class);
    }
    const classes = cssClasses.join(" ");
    const getShapeFromType = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((type) => {
      switch (type) {
        case nodeType.CIRCLE:
          return "mindmapCircle";
        case nodeType.RECT:
          return "rect";
        case nodeType.ROUNDED_RECT:
          return "rounded";
        case nodeType.CLOUD:
          return "cloud";
        case nodeType.BANG:
          return "bang";
        case nodeType.HEXAGON:
          return "hexagon";
        case nodeType.DEFAULT:
          return "defaultMindmapNode";
        case nodeType.NO_BORDER:
        default:
          return "rect";
      }
    }, "getShapeFromType");
    const processedNode = {
      id: node.id.toString(),
      domId: "node_" + node.id.toString(),
      label: node.descr,
      isGroup: false,
      shape: getShapeFromType(node.type),
      width: node.width,
      height: node.height ?? 0,
      padding: node.padding,
      cssClasses: classes,
      cssStyles: [],
      look: "default",
      icon: node.icon,
      x: node.x,
      y: node.y,
      // Mindmap-specific properties
      level: node.level,
      nodeId: node.nodeId,
      type: node.type,
      section: node.section
    };
    processedNodes.push(processedNode);
    if (node.children) {
      for (const child of node.children) {
        this.flattenNodes(child, processedNodes);
      }
    }
  }
  /**
   * Generate edges from parent-child relationships in mindmap tree
   * @param node - The mindmap node to process
   * @param edges - Array to collect edges
   */
  generateEdges(node, edges) {
    if (!node.children) {
      return;
    }
    for (const child of node.children) {
      let edgeClasses = "edge";
      if (child.section !== void 0) {
        edgeClasses += ` section-edge-${child.section}`;
      }
      const edgeDepth = node.level + 1;
      edgeClasses += ` edge-depth-${edgeDepth}`;
      const edge = {
        id: `edge_${node.id}_${child.id}`,
        start: node.id.toString(),
        end: child.id.toString(),
        type: "normal",
        curve: "basis",
        thickness: "normal",
        look: "default",
        classes: edgeClasses,
        // Store mindmap-specific data
        depth: node.level,
        section: child.section
      };
      edges.push(edge);
      this.generateEdges(child, edges);
    }
  }
  /**
   * Get structured data for layout algorithms
   * Following the pattern established by ER diagrams
   * @returns Structured data containing nodes, edges, and config
   */
  getData() {
    const mindmapRoot = this.getMindmap();
    const config = (0,chunk_ABZYJK2D/* getConfig2 */.nV)();
    const userDefinedConfig = (0,chunk_ABZYJK2D/* getUserDefinedConfig */.qO)();
    const hasUserDefinedLayout = userDefinedConfig.layout !== void 0;
    const finalConfig = config;
    if (!hasUserDefinedLayout) {
      finalConfig.layout = "cose-bilkent";
    }
    if (!mindmapRoot) {
      return {
        nodes: [],
        edges: [],
        config: finalConfig
      };
    }
    chunk_AGHRB4JF/* log */.cM.debug("getData: mindmapRoot", mindmapRoot, config);
    this.assignSections(mindmapRoot);
    const processedNodes = [];
    const processedEdges = [];
    this.flattenNodes(mindmapRoot, processedNodes);
    this.generateEdges(mindmapRoot, processedEdges);
    chunk_AGHRB4JF/* log */.cM.debug(
      `getData: processed ${processedNodes.length} nodes and ${processedEdges.length} edges`
    );
    const shapes = /* @__PURE__ */ new Map();
    for (const node of processedNodes) {
      shapes.set(node.id, {
        shape: node.shape,
        width: node.width,
        height: node.height,
        padding: node.padding
      });
    }
    return {
      nodes: processedNodes,
      edges: processedEdges,
      config: finalConfig,
      // Store the root node for mindmap-specific layout algorithms
      rootNode: mindmapRoot,
      // Properties required by dagre layout algorithm
      markers: ["point"],
      // Mindmaps don't use markers
      direction: "TB",
      // Top-to-bottom direction for mindmaps
      nodeSpacing: 50,
      // Default spacing between nodes
      rankSpacing: 50,
      // Default spacing between ranks
      // Add shapes for ELK compatibility
      shapes: Object.fromEntries(shapes),
      // Additional properties that layout algorithms might expect
      type: "mindmap",
      diagramId: "mindmap-" + esm_browser_v4()
    };
  }
  // Expose logger to grammar
  getLogger() {
    return chunk_AGHRB4JF/* log */.cM;
  }
};

// src/diagrams/mindmap/mindmapRenderer.ts
var draw = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async (text, id, _version, diagObj) => {
  chunk_AGHRB4JF/* log */.cM.debug("Rendering mindmap diagram\n" + text);
  const db = diagObj.db;
  const data4Layout = db.getData();
  const svg = (0,chunk_55IACEB6/* getDiagramElement */.q)(id, data4Layout.config.securityLevel);
  data4Layout.type = diagObj.type;
  data4Layout.layoutAlgorithm = (0,chunk_N4CR4FBY/* getRegisteredLayoutAlgorithm */._b)(data4Layout.config.layout, {
    fallback: "cose-bilkent"
  });
  data4Layout.diagramId = id;
  const mm = db.getMindmap();
  if (!mm) {
    return;
  }
  data4Layout.nodes.forEach((node) => {
    if (node.shape === "rounded") {
      node.radius = 15;
      node.taper = 15;
      node.stroke = "none";
      node.width = 0;
      node.padding = 15;
    } else if (node.shape === "circle") {
      node.padding = 10;
    } else if (node.shape === "rect") {
      node.width = 0;
      node.padding = 10;
    }
  });
  await (0,chunk_N4CR4FBY/* render */.sY)(data4Layout, svg);
  (0,chunk_QN33PNHL/* setupViewPortForSVG */.j)(
    svg,
    data4Layout.config.mindmap?.padding ?? chunk_ABZYJK2D/* defaultConfig_default */.vZ.mindmap.padding,
    "mindmapDiagram",
    data4Layout.config.mindmap?.useMaxWidth ?? chunk_ABZYJK2D/* defaultConfig_default */.vZ.mindmap.useMaxWidth
  );
}, "draw");
var mindmapRenderer_default = {
  draw
};

// src/diagrams/mindmap/styles.ts

var genSections = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((options) => {
  let sections = "";
  for (let i = 0; i < options.THEME_COLOR_LIMIT; i++) {
    options["lineColor" + i] = options["lineColor" + i] || options["cScaleInv" + i];
    if ((0,is_dark/* default */.Z)(options["lineColor" + i])) {
      options["lineColor" + i] = (0,lighten/* default */.Z)(options["lineColor" + i], 20);
    } else {
      options["lineColor" + i] = (0,darken/* default */.Z)(options["lineColor" + i], 20);
    }
  }
  for (let i = 0; i < options.THEME_COLOR_LIMIT; i++) {
    const sw = "" + (17 - 3 * i);
    sections += `
    .section-${i - 1} rect, .section-${i - 1} path, .section-${i - 1} circle, .section-${i - 1} polygon, .section-${i - 1} path  {
      fill: ${options["cScale" + i]};
    }
    .section-${i - 1} text {
     fill: ${options["cScaleLabel" + i]};
    }
    .node-icon-${i - 1} {
      font-size: 40px;
      color: ${options["cScaleLabel" + i]};
    }
    .section-edge-${i - 1}{
      stroke: ${options["cScale" + i]};
    }
    .edge-depth-${i - 1}{
      stroke-width: ${sw};
    }
    .section-${i - 1} line {
      stroke: ${options["cScaleInv" + i]} ;
      stroke-width: 3;
    }

    .disabled, .disabled circle, .disabled text {
      fill: lightgray;
    }
    .disabled text {
      fill: #efefef;
    }
    `;
  }
  return sections;
}, "genSections");
var getStyles = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((options) => `
  .edge {
    stroke-width: 3;
  }
  ${genSections(options)}
  .section-root rect, .section-root path, .section-root circle, .section-root polygon  {
    fill: ${options.git0};
  }
  .section-root text {
    fill: ${options.gitBranchLabel0};
  }
  .section-root span {
    color: ${options.gitBranchLabel0};
  }
  .section-2 span {
    color: ${options.gitBranchLabel0};
  }
  .icon-container {
    height:100%;
    display: flex;
    justify-content: center;
    align-items: center;
  }
  .edge {
    fill: none;
  }
  .mindmap-node-label {
    dy: 1em;
    alignment-baseline: middle;
    text-anchor: middle;
    dominant-baseline: middle;
    text-align: center;
  }
`, "getStyles");
var styles_default = getStyles;

// src/diagrams/mindmap/mindmap-definition.ts
var diagram = {
  get db() {
    return new MindmapDB();
  },
  renderer: mindmapRenderer_default,
  parser: mindmap_default,
  styles: styles_default
};



/***/ })

}]);
//# sourceMappingURL=3522.467e51019327266c2d99.js.map?v=467e51019327266c2d99
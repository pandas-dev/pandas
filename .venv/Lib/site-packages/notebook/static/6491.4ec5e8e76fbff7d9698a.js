"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[6491],{

/***/ 86491:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _chunk_C3MQ5ANM_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(39769);
/* harmony import */ var _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(42626);
/* harmony import */ var _chunk_7B677QYD_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(30103);
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(86906);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(83619);





// src/diagrams/xychart/parser/xychart.jison
var parser = function() {
  var o = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(k, v, o2, l) {
    for (o2 = o2 || {}, l = k.length; l--; o2[k[l]] = v) ;
    return o2;
  }, "o"), $V0 = [1, 10, 12, 14, 16, 18, 19, 21, 23], $V1 = [2, 6], $V2 = [1, 3], $V3 = [1, 5], $V4 = [1, 6], $V5 = [1, 7], $V6 = [1, 5, 10, 12, 14, 16, 18, 19, 21, 23, 34, 35, 36], $V7 = [1, 25], $V8 = [1, 26], $V9 = [1, 28], $Va = [1, 29], $Vb = [1, 30], $Vc = [1, 31], $Vd = [1, 32], $Ve = [1, 33], $Vf = [1, 34], $Vg = [1, 35], $Vh = [1, 36], $Vi = [1, 37], $Vj = [1, 43], $Vk = [1, 42], $Vl = [1, 47], $Vm = [1, 50], $Vn = [1, 10, 12, 14, 16, 18, 19, 21, 23, 34, 35, 36], $Vo = [1, 10, 12, 14, 16, 18, 19, 21, 23, 24, 26, 27, 28, 34, 35, 36], $Vp = [1, 10, 12, 14, 16, 18, 19, 21, 23, 24, 26, 27, 28, 34, 35, 36, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50], $Vq = [1, 64];
  var parser2 = {
    trace: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function trace() {
    }, "trace"),
    yy: {},
    symbols_: { "error": 2, "start": 3, "eol": 4, "XYCHART": 5, "chartConfig": 6, "document": 7, "CHART_ORIENTATION": 8, "statement": 9, "title": 10, "text": 11, "X_AXIS": 12, "parseXAxis": 13, "Y_AXIS": 14, "parseYAxis": 15, "LINE": 16, "plotData": 17, "BAR": 18, "acc_title": 19, "acc_title_value": 20, "acc_descr": 21, "acc_descr_value": 22, "acc_descr_multiline_value": 23, "SQUARE_BRACES_START": 24, "commaSeparatedNumbers": 25, "SQUARE_BRACES_END": 26, "NUMBER_WITH_DECIMAL": 27, "COMMA": 28, "xAxisData": 29, "bandData": 30, "ARROW_DELIMITER": 31, "commaSeparatedTexts": 32, "yAxisData": 33, "NEWLINE": 34, "SEMI": 35, "EOF": 36, "alphaNum": 37, "STR": 38, "MD_STR": 39, "alphaNumToken": 40, "AMP": 41, "NUM": 42, "ALPHA": 43, "PLUS": 44, "EQUALS": 45, "MULT": 46, "DOT": 47, "BRKT": 48, "MINUS": 49, "UNDERSCORE": 50, "$accept": 0, "$end": 1 },
    terminals_: { 2: "error", 5: "XYCHART", 8: "CHART_ORIENTATION", 10: "title", 12: "X_AXIS", 14: "Y_AXIS", 16: "LINE", 18: "BAR", 19: "acc_title", 20: "acc_title_value", 21: "acc_descr", 22: "acc_descr_value", 23: "acc_descr_multiline_value", 24: "SQUARE_BRACES_START", 26: "SQUARE_BRACES_END", 27: "NUMBER_WITH_DECIMAL", 28: "COMMA", 31: "ARROW_DELIMITER", 34: "NEWLINE", 35: "SEMI", 36: "EOF", 38: "STR", 39: "MD_STR", 41: "AMP", 42: "NUM", 43: "ALPHA", 44: "PLUS", 45: "EQUALS", 46: "MULT", 47: "DOT", 48: "BRKT", 49: "MINUS", 50: "UNDERSCORE" },
    productions_: [0, [3, 2], [3, 3], [3, 2], [3, 1], [6, 1], [7, 0], [7, 2], [9, 2], [9, 2], [9, 2], [9, 2], [9, 2], [9, 3], [9, 2], [9, 3], [9, 2], [9, 2], [9, 1], [17, 3], [25, 3], [25, 1], [13, 1], [13, 2], [13, 1], [29, 1], [29, 3], [30, 3], [32, 3], [32, 1], [15, 1], [15, 2], [15, 1], [33, 3], [4, 1], [4, 1], [4, 1], [11, 1], [11, 1], [11, 1], [37, 1], [37, 2], [40, 1], [40, 1], [40, 1], [40, 1], [40, 1], [40, 1], [40, 1], [40, 1], [40, 1], [40, 1]],
    performAction: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function anonymous(yytext, yyleng, yylineno, yy, yystate, $$, _$) {
      var $0 = $$.length - 1;
      switch (yystate) {
        case 5:
          yy.setOrientation($$[$0]);
          break;
        case 9:
          yy.setDiagramTitle($$[$0].text.trim());
          break;
        case 12:
          yy.setLineData({ text: "", type: "text" }, $$[$0]);
          break;
        case 13:
          yy.setLineData($$[$0 - 1], $$[$0]);
          break;
        case 14:
          yy.setBarData({ text: "", type: "text" }, $$[$0]);
          break;
        case 15:
          yy.setBarData($$[$0 - 1], $$[$0]);
          break;
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
          this.$ = $$[$0 - 1];
          break;
        case 20:
          this.$ = [Number($$[$0 - 2]), ...$$[$0]];
          break;
        case 21:
          this.$ = [Number($$[$0])];
          break;
        case 22:
          yy.setXAxisTitle($$[$0]);
          break;
        case 23:
          yy.setXAxisTitle($$[$0 - 1]);
          break;
        case 24:
          yy.setXAxisTitle({ type: "text", text: "" });
          break;
        case 25:
          yy.setXAxisBand($$[$0]);
          break;
        case 26:
          yy.setXAxisRangeData(Number($$[$0 - 2]), Number($$[$0]));
          break;
        case 27:
          this.$ = $$[$0 - 1];
          break;
        case 28:
          this.$ = [$$[$0 - 2], ...$$[$0]];
          break;
        case 29:
          this.$ = [$$[$0]];
          break;
        case 30:
          yy.setYAxisTitle($$[$0]);
          break;
        case 31:
          yy.setYAxisTitle($$[$0 - 1]);
          break;
        case 32:
          yy.setYAxisTitle({ type: "text", text: "" });
          break;
        case 33:
          yy.setYAxisRangeData(Number($$[$0 - 2]), Number($$[$0]));
          break;
        case 37:
          this.$ = { text: $$[$0], type: "text" };
          break;
        case 38:
          this.$ = { text: $$[$0], type: "text" };
          break;
        case 39:
          this.$ = { text: $$[$0], type: "markdown" };
          break;
        case 40:
          this.$ = $$[$0];
          break;
        case 41:
          this.$ = $$[$0 - 1] + "" + $$[$0];
          break;
      }
    }, "anonymous"),
    table: [o($V0, $V1, { 3: 1, 4: 2, 7: 4, 5: $V2, 34: $V3, 35: $V4, 36: $V5 }), { 1: [3] }, o($V0, $V1, { 4: 2, 7: 4, 3: 8, 5: $V2, 34: $V3, 35: $V4, 36: $V5 }), o($V0, $V1, { 4: 2, 7: 4, 6: 9, 3: 10, 5: $V2, 8: [1, 11], 34: $V3, 35: $V4, 36: $V5 }), { 1: [2, 4], 9: 12, 10: [1, 13], 12: [1, 14], 14: [1, 15], 16: [1, 16], 18: [1, 17], 19: [1, 18], 21: [1, 19], 23: [1, 20] }, o($V6, [2, 34]), o($V6, [2, 35]), o($V6, [2, 36]), { 1: [2, 1] }, o($V0, $V1, { 4: 2, 7: 4, 3: 21, 5: $V2, 34: $V3, 35: $V4, 36: $V5 }), { 1: [2, 3] }, o($V6, [2, 5]), o($V0, [2, 7], { 4: 22, 34: $V3, 35: $V4, 36: $V5 }), { 11: 23, 37: 24, 38: $V7, 39: $V8, 40: 27, 41: $V9, 42: $Va, 43: $Vb, 44: $Vc, 45: $Vd, 46: $Ve, 47: $Vf, 48: $Vg, 49: $Vh, 50: $Vi }, { 11: 39, 13: 38, 24: $Vj, 27: $Vk, 29: 40, 30: 41, 37: 24, 38: $V7, 39: $V8, 40: 27, 41: $V9, 42: $Va, 43: $Vb, 44: $Vc, 45: $Vd, 46: $Ve, 47: $Vf, 48: $Vg, 49: $Vh, 50: $Vi }, { 11: 45, 15: 44, 27: $Vl, 33: 46, 37: 24, 38: $V7, 39: $V8, 40: 27, 41: $V9, 42: $Va, 43: $Vb, 44: $Vc, 45: $Vd, 46: $Ve, 47: $Vf, 48: $Vg, 49: $Vh, 50: $Vi }, { 11: 49, 17: 48, 24: $Vm, 37: 24, 38: $V7, 39: $V8, 40: 27, 41: $V9, 42: $Va, 43: $Vb, 44: $Vc, 45: $Vd, 46: $Ve, 47: $Vf, 48: $Vg, 49: $Vh, 50: $Vi }, { 11: 52, 17: 51, 24: $Vm, 37: 24, 38: $V7, 39: $V8, 40: 27, 41: $V9, 42: $Va, 43: $Vb, 44: $Vc, 45: $Vd, 46: $Ve, 47: $Vf, 48: $Vg, 49: $Vh, 50: $Vi }, { 20: [1, 53] }, { 22: [1, 54] }, o($Vn, [2, 18]), { 1: [2, 2] }, o($Vn, [2, 8]), o($Vn, [2, 9]), o($Vo, [2, 37], { 40: 55, 41: $V9, 42: $Va, 43: $Vb, 44: $Vc, 45: $Vd, 46: $Ve, 47: $Vf, 48: $Vg, 49: $Vh, 50: $Vi }), o($Vo, [2, 38]), o($Vo, [2, 39]), o($Vp, [2, 40]), o($Vp, [2, 42]), o($Vp, [2, 43]), o($Vp, [2, 44]), o($Vp, [2, 45]), o($Vp, [2, 46]), o($Vp, [2, 47]), o($Vp, [2, 48]), o($Vp, [2, 49]), o($Vp, [2, 50]), o($Vp, [2, 51]), o($Vn, [2, 10]), o($Vn, [2, 22], { 30: 41, 29: 56, 24: $Vj, 27: $Vk }), o($Vn, [2, 24]), o($Vn, [2, 25]), { 31: [1, 57] }, { 11: 59, 32: 58, 37: 24, 38: $V7, 39: $V8, 40: 27, 41: $V9, 42: $Va, 43: $Vb, 44: $Vc, 45: $Vd, 46: $Ve, 47: $Vf, 48: $Vg, 49: $Vh, 50: $Vi }, o($Vn, [2, 11]), o($Vn, [2, 30], { 33: 60, 27: $Vl }), o($Vn, [2, 32]), { 31: [1, 61] }, o($Vn, [2, 12]), { 17: 62, 24: $Vm }, { 25: 63, 27: $Vq }, o($Vn, [2, 14]), { 17: 65, 24: $Vm }, o($Vn, [2, 16]), o($Vn, [2, 17]), o($Vp, [2, 41]), o($Vn, [2, 23]), { 27: [1, 66] }, { 26: [1, 67] }, { 26: [2, 29], 28: [1, 68] }, o($Vn, [2, 31]), { 27: [1, 69] }, o($Vn, [2, 13]), { 26: [1, 70] }, { 26: [2, 21], 28: [1, 71] }, o($Vn, [2, 15]), o($Vn, [2, 26]), o($Vn, [2, 27]), { 11: 59, 32: 72, 37: 24, 38: $V7, 39: $V8, 40: 27, 41: $V9, 42: $Va, 43: $Vb, 44: $Vc, 45: $Vd, 46: $Ve, 47: $Vf, 48: $Vg, 49: $Vh, 50: $Vi }, o($Vn, [2, 33]), o($Vn, [2, 19]), { 25: 73, 27: $Vq }, { 26: [2, 28] }, { 26: [2, 20] }],
    defaultActions: { 8: [2, 1], 10: [2, 3], 21: [2, 2], 72: [2, 28], 73: [2, 20] },
    parseError: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function parseError(str, hash) {
      if (hash.recoverable) {
        this.trace(str);
      } else {
        var error = new Error(str);
        error.hash = hash;
        throw error;
      }
    }, "parseError"),
    parse: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function parse(input) {
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
      (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(popStack, "popStack");
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
      (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(lex, "lex");
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
  var lexer = /* @__PURE__ */ function() {
    var lexer2 = {
      EOF: 1,
      parseError: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function parseError(str, hash) {
        if (this.yy.parser) {
          this.yy.parser.parseError(str, hash);
        } else {
          throw new Error(str);
        }
      }, "parseError"),
      // resets the lexer, sets new input
      setInput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(input, yy) {
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
      input: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
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
      unput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(ch) {
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
      more: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        this._more = true;
        return this;
      }, "more"),
      // When called from action, signals the lexer that this rule fails to match the input, so the next matching rule (regex) should be tested instead.
      reject: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
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
      less: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(n) {
        this.unput(this.match.slice(n));
      }, "less"),
      // displays already matched input, i.e. for error messages
      pastInput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        var past = this.matched.substr(0, this.matched.length - this.match.length);
        return (past.length > 20 ? "..." : "") + past.substr(-20).replace(/\n/g, "");
      }, "pastInput"),
      // displays upcoming input, i.e. for error messages
      upcomingInput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        var next = this.match;
        if (next.length < 20) {
          next += this._input.substr(0, 20 - next.length);
        }
        return (next.substr(0, 20) + (next.length > 20 ? "..." : "")).replace(/\n/g, "");
      }, "upcomingInput"),
      // displays the character position where the lexing error occurred, i.e. for error messages
      showPosition: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        var pre = this.pastInput();
        var c = new Array(pre.length + 1).join("-");
        return pre + this.upcomingInput() + "\n" + c + "^";
      }, "showPosition"),
      // test the lexed token: return FALSE when not a match, otherwise return token
      test_match: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(match, indexed_rule) {
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
      next: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
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
      lex: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function lex() {
        var r = this.next();
        if (r) {
          return r;
        } else {
          return this.lex();
        }
      }, "lex"),
      // activates a new lexer condition state (pushes the new lexer condition state onto the condition stack)
      begin: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function begin(condition) {
        this.conditionStack.push(condition);
      }, "begin"),
      // pop the previously active lexer condition state off the condition stack
      popState: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function popState() {
        var n = this.conditionStack.length - 1;
        if (n > 0) {
          return this.conditionStack.pop();
        } else {
          return this.conditionStack[0];
        }
      }, "popState"),
      // produce the lexer rule set which is active for the currently active lexer condition state
      _currentRules: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function _currentRules() {
        if (this.conditionStack.length && this.conditionStack[this.conditionStack.length - 1]) {
          return this.conditions[this.conditionStack[this.conditionStack.length - 1]].rules;
        } else {
          return this.conditions["INITIAL"].rules;
        }
      }, "_currentRules"),
      // return the currently active lexer condition state; when an index argument is provided it produces the N-th previous condition state, if available
      topState: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function topState(n) {
        n = this.conditionStack.length - 1 - Math.abs(n || 0);
        if (n >= 0) {
          return this.conditionStack[n];
        } else {
          return "INITIAL";
        }
      }, "topState"),
      // alias for begin(condition)
      pushState: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function pushState(condition) {
        this.begin(condition);
      }, "pushState"),
      // return the number of states currently on the stack
      stateStackSize: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function stateStackSize() {
        return this.conditionStack.length;
      }, "stateStackSize"),
      options: { "case-insensitive": true },
      performAction: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function anonymous(yy, yy_, $avoiding_name_collisions, YY_START) {
        var YYSTATE = YY_START;
        switch ($avoiding_name_collisions) {
          case 0:
            break;
          case 1:
            break;
          case 2:
            this.popState();
            return 34;
            break;
          case 3:
            this.popState();
            return 34;
            break;
          case 4:
            return 34;
            break;
          case 5:
            break;
          case 6:
            return 10;
            break;
          case 7:
            this.pushState("acc_title");
            return 19;
            break;
          case 8:
            this.popState();
            return "acc_title_value";
            break;
          case 9:
            this.pushState("acc_descr");
            return 21;
            break;
          case 10:
            this.popState();
            return "acc_descr_value";
            break;
          case 11:
            this.pushState("acc_descr_multiline");
            break;
          case 12:
            this.popState();
            break;
          case 13:
            return "acc_descr_multiline_value";
            break;
          case 14:
            return 5;
            break;
          case 15:
            return 8;
            break;
          case 16:
            this.pushState("axis_data");
            return "X_AXIS";
            break;
          case 17:
            this.pushState("axis_data");
            return "Y_AXIS";
            break;
          case 18:
            this.pushState("axis_band_data");
            return 24;
            break;
          case 19:
            return 31;
            break;
          case 20:
            this.pushState("data");
            return 16;
            break;
          case 21:
            this.pushState("data");
            return 18;
            break;
          case 22:
            this.pushState("data_inner");
            return 24;
            break;
          case 23:
            return 27;
            break;
          case 24:
            this.popState();
            return 26;
            break;
          case 25:
            this.popState();
            break;
          case 26:
            this.pushState("string");
            break;
          case 27:
            this.popState();
            break;
          case 28:
            return "STR";
            break;
          case 29:
            return 24;
            break;
          case 30:
            return 26;
            break;
          case 31:
            return 43;
            break;
          case 32:
            return "COLON";
            break;
          case 33:
            return 44;
            break;
          case 34:
            return 28;
            break;
          case 35:
            return 45;
            break;
          case 36:
            return 46;
            break;
          case 37:
            return 48;
            break;
          case 38:
            return 50;
            break;
          case 39:
            return 47;
            break;
          case 40:
            return 41;
            break;
          case 41:
            return 49;
            break;
          case 42:
            return 42;
            break;
          case 43:
            break;
          case 44:
            return 35;
            break;
          case 45:
            return 36;
            break;
        }
      }, "anonymous"),
      rules: [/^(?:%%(?!\{)[^\n]*)/i, /^(?:[^\}]%%[^\n]*)/i, /^(?:(\r?\n))/i, /^(?:(\r?\n))/i, /^(?:[\n\r]+)/i, /^(?:%%[^\n]*)/i, /^(?:title\b)/i, /^(?:accTitle\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*\{\s*)/i, /^(?:\{)/i, /^(?:[^\}]*)/i, /^(?:xychart-beta\b)/i, /^(?:(?:vertical|horizontal))/i, /^(?:x-axis\b)/i, /^(?:y-axis\b)/i, /^(?:\[)/i, /^(?:-->)/i, /^(?:line\b)/i, /^(?:bar\b)/i, /^(?:\[)/i, /^(?:[+-]?(?:\d+(?:\.\d+)?|\.\d+))/i, /^(?:\])/i, /^(?:(?:`\)                                    \{ this\.pushState\(md_string\); \}\n<md_string>\(\?:\(\?!`"\)\.\)\+                  \{ return MD_STR; \}\n<md_string>\(\?:`))/i, /^(?:["])/i, /^(?:["])/i, /^(?:[^"]*)/i, /^(?:\[)/i, /^(?:\])/i, /^(?:[A-Za-z]+)/i, /^(?::)/i, /^(?:\+)/i, /^(?:,)/i, /^(?:=)/i, /^(?:\*)/i, /^(?:#)/i, /^(?:[\_])/i, /^(?:\.)/i, /^(?:&)/i, /^(?:-)/i, /^(?:[0-9]+)/i, /^(?:\s+)/i, /^(?:;)/i, /^(?:$)/i],
      conditions: { "data_inner": { "rules": [0, 1, 4, 5, 6, 7, 9, 11, 14, 15, 16, 17, 20, 21, 23, 24, 25, 26, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45], "inclusive": true }, "data": { "rules": [0, 1, 3, 4, 5, 6, 7, 9, 11, 14, 15, 16, 17, 20, 21, 22, 25, 26, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45], "inclusive": true }, "axis_band_data": { "rules": [0, 1, 4, 5, 6, 7, 9, 11, 14, 15, 16, 17, 20, 21, 24, 25, 26, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45], "inclusive": true }, "axis_data": { "rules": [0, 1, 2, 4, 5, 6, 7, 9, 11, 14, 15, 16, 17, 18, 19, 20, 21, 23, 25, 26, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45], "inclusive": true }, "acc_descr_multiline": { "rules": [12, 13], "inclusive": false }, "acc_descr": { "rules": [10], "inclusive": false }, "acc_title": { "rules": [8], "inclusive": false }, "title": { "rules": [], "inclusive": false }, "md_string": { "rules": [], "inclusive": false }, "string": { "rules": [27, 28], "inclusive": false }, "INITIAL": { "rules": [0, 1, 4, 5, 6, 7, 9, 11, 14, 15, 16, 17, 20, 21, 25, 26, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45], "inclusive": true } }
    };
    return lexer2;
  }();
  parser2.lexer = lexer;
  function Parser() {
    this.yy = {};
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(Parser, "Parser");
  Parser.prototype = parser2;
  parser2.Parser = Parser;
  return new Parser();
}();
parser.parser = parser;
var xychart_default = parser;

// src/diagrams/xychart/chartBuilder/interfaces.ts
function isBarPlot(data) {
  return data.type === "bar";
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(isBarPlot, "isBarPlot");
function isBandAxisData(data) {
  return data.type === "band";
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(isBandAxisData, "isBandAxisData");
function isLinearAxisData(data) {
  return data.type === "linear";
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(isLinearAxisData, "isLinearAxisData");

// src/diagrams/xychart/chartBuilder/textDimensionCalculator.ts
var TextDimensionCalculatorWithFont = class {
  constructor(parentGroup) {
    this.parentGroup = parentGroup;
  }
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(this, "TextDimensionCalculatorWithFont");
  }
  getMaxDimension(texts, fontSize) {
    if (!this.parentGroup) {
      return {
        width: texts.reduce((acc, cur) => Math.max(cur.length, acc), 0) * fontSize,
        height: fontSize
      };
    }
    const dimension = {
      width: 0,
      height: 0
    };
    const elem = this.parentGroup.append("g").attr("visibility", "hidden").attr("font-size", fontSize);
    for (const t of texts) {
      const bbox = (0,_chunk_C3MQ5ANM_mjs__WEBPACK_IMPORTED_MODULE_0__/* .computeDimensionOfText */ .QA)(elem, 1, t);
      const width = bbox ? bbox.width : t.length * fontSize;
      const height = bbox ? bbox.height : fontSize;
      dimension.width = Math.max(dimension.width, width);
      dimension.height = Math.max(dimension.height, height);
    }
    elem.remove();
    return dimension;
  }
};

// src/diagrams/xychart/chartBuilder/components/axis/bandAxis.ts


// src/diagrams/xychart/chartBuilder/components/axis/baseAxis.ts
var BAR_WIDTH_TO_TICK_WIDTH_RATIO = 0.7;
var MAX_OUTER_PADDING_PERCENT_FOR_WRT_LABEL = 0.2;
var BaseAxis = class {
  constructor(axisConfig, title, textDimensionCalculator, axisThemeConfig) {
    this.axisConfig = axisConfig;
    this.title = title;
    this.textDimensionCalculator = textDimensionCalculator;
    this.axisThemeConfig = axisThemeConfig;
    this.boundingRect = { x: 0, y: 0, width: 0, height: 0 };
    this.axisPosition = "left";
    this.showTitle = false;
    this.showLabel = false;
    this.showTick = false;
    this.showAxisLine = false;
    this.outerPadding = 0;
    this.titleTextHeight = 0;
    this.labelTextHeight = 0;
    this.range = [0, 10];
    this.boundingRect = { x: 0, y: 0, width: 0, height: 0 };
    this.axisPosition = "left";
  }
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(this, "BaseAxis");
  }
  setRange(range) {
    this.range = range;
    if (this.axisPosition === "left" || this.axisPosition === "right") {
      this.boundingRect.height = range[1] - range[0];
    } else {
      this.boundingRect.width = range[1] - range[0];
    }
    this.recalculateScale();
  }
  getRange() {
    return [this.range[0] + this.outerPadding, this.range[1] - this.outerPadding];
  }
  setAxisPosition(axisPosition) {
    this.axisPosition = axisPosition;
    this.setRange(this.range);
  }
  getTickDistance() {
    const range = this.getRange();
    return Math.abs(range[0] - range[1]) / this.getTickValues().length;
  }
  getAxisOuterPadding() {
    return this.outerPadding;
  }
  getLabelDimension() {
    return this.textDimensionCalculator.getMaxDimension(
      this.getTickValues().map((tick) => tick.toString()),
      this.axisConfig.labelFontSize
    );
  }
  recalculateOuterPaddingToDrawBar() {
    if (BAR_WIDTH_TO_TICK_WIDTH_RATIO * this.getTickDistance() > this.outerPadding * 2) {
      this.outerPadding = Math.floor(BAR_WIDTH_TO_TICK_WIDTH_RATIO * this.getTickDistance() / 2);
    }
    this.recalculateScale();
  }
  calculateSpaceIfDrawnHorizontally(availableSpace) {
    let availableHeight = availableSpace.height;
    if (this.axisConfig.showAxisLine && availableHeight > this.axisConfig.axisLineWidth) {
      availableHeight -= this.axisConfig.axisLineWidth;
      this.showAxisLine = true;
    }
    if (this.axisConfig.showLabel) {
      const spaceRequired = this.getLabelDimension();
      const maxPadding = MAX_OUTER_PADDING_PERCENT_FOR_WRT_LABEL * availableSpace.width;
      this.outerPadding = Math.min(spaceRequired.width / 2, maxPadding);
      const heightRequired = spaceRequired.height + this.axisConfig.labelPadding * 2;
      this.labelTextHeight = spaceRequired.height;
      if (heightRequired <= availableHeight) {
        availableHeight -= heightRequired;
        this.showLabel = true;
      }
    }
    if (this.axisConfig.showTick && availableHeight >= this.axisConfig.tickLength) {
      this.showTick = true;
      availableHeight -= this.axisConfig.tickLength;
    }
    if (this.axisConfig.showTitle && this.title) {
      const spaceRequired = this.textDimensionCalculator.getMaxDimension(
        [this.title],
        this.axisConfig.titleFontSize
      );
      const heightRequired = spaceRequired.height + this.axisConfig.titlePadding * 2;
      this.titleTextHeight = spaceRequired.height;
      if (heightRequired <= availableHeight) {
        availableHeight -= heightRequired;
        this.showTitle = true;
      }
    }
    this.boundingRect.width = availableSpace.width;
    this.boundingRect.height = availableSpace.height - availableHeight;
  }
  calculateSpaceIfDrawnVertical(availableSpace) {
    let availableWidth = availableSpace.width;
    if (this.axisConfig.showAxisLine && availableWidth > this.axisConfig.axisLineWidth) {
      availableWidth -= this.axisConfig.axisLineWidth;
      this.showAxisLine = true;
    }
    if (this.axisConfig.showLabel) {
      const spaceRequired = this.getLabelDimension();
      const maxPadding = MAX_OUTER_PADDING_PERCENT_FOR_WRT_LABEL * availableSpace.height;
      this.outerPadding = Math.min(spaceRequired.height / 2, maxPadding);
      const widthRequired = spaceRequired.width + this.axisConfig.labelPadding * 2;
      if (widthRequired <= availableWidth) {
        availableWidth -= widthRequired;
        this.showLabel = true;
      }
    }
    if (this.axisConfig.showTick && availableWidth >= this.axisConfig.tickLength) {
      this.showTick = true;
      availableWidth -= this.axisConfig.tickLength;
    }
    if (this.axisConfig.showTitle && this.title) {
      const spaceRequired = this.textDimensionCalculator.getMaxDimension(
        [this.title],
        this.axisConfig.titleFontSize
      );
      const widthRequired = spaceRequired.height + this.axisConfig.titlePadding * 2;
      this.titleTextHeight = spaceRequired.height;
      if (widthRequired <= availableWidth) {
        availableWidth -= widthRequired;
        this.showTitle = true;
      }
    }
    this.boundingRect.width = availableSpace.width - availableWidth;
    this.boundingRect.height = availableSpace.height;
  }
  calculateSpace(availableSpace) {
    if (this.axisPosition === "left" || this.axisPosition === "right") {
      this.calculateSpaceIfDrawnVertical(availableSpace);
    } else {
      this.calculateSpaceIfDrawnHorizontally(availableSpace);
    }
    this.recalculateScale();
    return {
      width: this.boundingRect.width,
      height: this.boundingRect.height
    };
  }
  setBoundingBoxXY(point) {
    this.boundingRect.x = point.x;
    this.boundingRect.y = point.y;
  }
  getDrawableElementsForLeftAxis() {
    const drawableElement = [];
    if (this.showAxisLine) {
      const x = this.boundingRect.x + this.boundingRect.width - this.axisConfig.axisLineWidth / 2;
      drawableElement.push({
        type: "path",
        groupTexts: ["left-axis", "axisl-line"],
        data: [
          {
            path: `M ${x},${this.boundingRect.y} L ${x},${this.boundingRect.y + this.boundingRect.height} `,
            strokeFill: this.axisThemeConfig.axisLineColor,
            strokeWidth: this.axisConfig.axisLineWidth
          }
        ]
      });
    }
    if (this.showLabel) {
      drawableElement.push({
        type: "text",
        groupTexts: ["left-axis", "label"],
        data: this.getTickValues().map((tick) => ({
          text: tick.toString(),
          x: this.boundingRect.x + this.boundingRect.width - (this.showLabel ? this.axisConfig.labelPadding : 0) - (this.showTick ? this.axisConfig.tickLength : 0) - (this.showAxisLine ? this.axisConfig.axisLineWidth : 0),
          y: this.getScaleValue(tick),
          fill: this.axisThemeConfig.labelColor,
          fontSize: this.axisConfig.labelFontSize,
          rotation: 0,
          verticalPos: "middle",
          horizontalPos: "right"
        }))
      });
    }
    if (this.showTick) {
      const x = this.boundingRect.x + this.boundingRect.width - (this.showAxisLine ? this.axisConfig.axisLineWidth : 0);
      drawableElement.push({
        type: "path",
        groupTexts: ["left-axis", "ticks"],
        data: this.getTickValues().map((tick) => ({
          path: `M ${x},${this.getScaleValue(tick)} L ${x - this.axisConfig.tickLength},${this.getScaleValue(tick)}`,
          strokeFill: this.axisThemeConfig.tickColor,
          strokeWidth: this.axisConfig.tickWidth
        }))
      });
    }
    if (this.showTitle) {
      drawableElement.push({
        type: "text",
        groupTexts: ["left-axis", "title"],
        data: [
          {
            text: this.title,
            x: this.boundingRect.x + this.axisConfig.titlePadding,
            y: this.boundingRect.y + this.boundingRect.height / 2,
            fill: this.axisThemeConfig.titleColor,
            fontSize: this.axisConfig.titleFontSize,
            rotation: 270,
            verticalPos: "top",
            horizontalPos: "center"
          }
        ]
      });
    }
    return drawableElement;
  }
  getDrawableElementsForBottomAxis() {
    const drawableElement = [];
    if (this.showAxisLine) {
      const y = this.boundingRect.y + this.axisConfig.axisLineWidth / 2;
      drawableElement.push({
        type: "path",
        groupTexts: ["bottom-axis", "axis-line"],
        data: [
          {
            path: `M ${this.boundingRect.x},${y} L ${this.boundingRect.x + this.boundingRect.width},${y}`,
            strokeFill: this.axisThemeConfig.axisLineColor,
            strokeWidth: this.axisConfig.axisLineWidth
          }
        ]
      });
    }
    if (this.showLabel) {
      drawableElement.push({
        type: "text",
        groupTexts: ["bottom-axis", "label"],
        data: this.getTickValues().map((tick) => ({
          text: tick.toString(),
          x: this.getScaleValue(tick),
          y: this.boundingRect.y + this.axisConfig.labelPadding + (this.showTick ? this.axisConfig.tickLength : 0) + (this.showAxisLine ? this.axisConfig.axisLineWidth : 0),
          fill: this.axisThemeConfig.labelColor,
          fontSize: this.axisConfig.labelFontSize,
          rotation: 0,
          verticalPos: "top",
          horizontalPos: "center"
        }))
      });
    }
    if (this.showTick) {
      const y = this.boundingRect.y + (this.showAxisLine ? this.axisConfig.axisLineWidth : 0);
      drawableElement.push({
        type: "path",
        groupTexts: ["bottom-axis", "ticks"],
        data: this.getTickValues().map((tick) => ({
          path: `M ${this.getScaleValue(tick)},${y} L ${this.getScaleValue(tick)},${y + this.axisConfig.tickLength}`,
          strokeFill: this.axisThemeConfig.tickColor,
          strokeWidth: this.axisConfig.tickWidth
        }))
      });
    }
    if (this.showTitle) {
      drawableElement.push({
        type: "text",
        groupTexts: ["bottom-axis", "title"],
        data: [
          {
            text: this.title,
            x: this.range[0] + (this.range[1] - this.range[0]) / 2,
            y: this.boundingRect.y + this.boundingRect.height - this.axisConfig.titlePadding - this.titleTextHeight,
            fill: this.axisThemeConfig.titleColor,
            fontSize: this.axisConfig.titleFontSize,
            rotation: 0,
            verticalPos: "top",
            horizontalPos: "center"
          }
        ]
      });
    }
    return drawableElement;
  }
  getDrawableElementsForTopAxis() {
    const drawableElement = [];
    if (this.showAxisLine) {
      const y = this.boundingRect.y + this.boundingRect.height - this.axisConfig.axisLineWidth / 2;
      drawableElement.push({
        type: "path",
        groupTexts: ["top-axis", "axis-line"],
        data: [
          {
            path: `M ${this.boundingRect.x},${y} L ${this.boundingRect.x + this.boundingRect.width},${y}`,
            strokeFill: this.axisThemeConfig.axisLineColor,
            strokeWidth: this.axisConfig.axisLineWidth
          }
        ]
      });
    }
    if (this.showLabel) {
      drawableElement.push({
        type: "text",
        groupTexts: ["top-axis", "label"],
        data: this.getTickValues().map((tick) => ({
          text: tick.toString(),
          x: this.getScaleValue(tick),
          y: this.boundingRect.y + (this.showTitle ? this.titleTextHeight + this.axisConfig.titlePadding * 2 : 0) + this.axisConfig.labelPadding,
          fill: this.axisThemeConfig.labelColor,
          fontSize: this.axisConfig.labelFontSize,
          rotation: 0,
          verticalPos: "top",
          horizontalPos: "center"
        }))
      });
    }
    if (this.showTick) {
      const y = this.boundingRect.y;
      drawableElement.push({
        type: "path",
        groupTexts: ["top-axis", "ticks"],
        data: this.getTickValues().map((tick) => ({
          path: `M ${this.getScaleValue(tick)},${y + this.boundingRect.height - (this.showAxisLine ? this.axisConfig.axisLineWidth : 0)} L ${this.getScaleValue(tick)},${y + this.boundingRect.height - this.axisConfig.tickLength - (this.showAxisLine ? this.axisConfig.axisLineWidth : 0)}`,
          strokeFill: this.axisThemeConfig.tickColor,
          strokeWidth: this.axisConfig.tickWidth
        }))
      });
    }
    if (this.showTitle) {
      drawableElement.push({
        type: "text",
        groupTexts: ["top-axis", "title"],
        data: [
          {
            text: this.title,
            x: this.boundingRect.x + this.boundingRect.width / 2,
            y: this.boundingRect.y + this.axisConfig.titlePadding,
            fill: this.axisThemeConfig.titleColor,
            fontSize: this.axisConfig.titleFontSize,
            rotation: 0,
            verticalPos: "top",
            horizontalPos: "center"
          }
        ]
      });
    }
    return drawableElement;
  }
  getDrawableElements() {
    if (this.axisPosition === "left") {
      return this.getDrawableElementsForLeftAxis();
    }
    if (this.axisPosition === "right") {
      throw Error("Drawing of right axis is not implemented");
    }
    if (this.axisPosition === "bottom") {
      return this.getDrawableElementsForBottomAxis();
    }
    if (this.axisPosition === "top") {
      return this.getDrawableElementsForTopAxis();
    }
    return [];
  }
};

// src/diagrams/xychart/chartBuilder/components/axis/bandAxis.ts
var BandAxis = class extends BaseAxis {
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(this, "BandAxis");
  }
  constructor(axisConfig, axisThemeConfig, categories, title, textDimensionCalculator) {
    super(axisConfig, title, textDimensionCalculator, axisThemeConfig);
    this.categories = categories;
    this.scale = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .scaleBand */ .tiA)().domain(this.categories).range(this.getRange());
  }
  setRange(range) {
    super.setRange(range);
  }
  recalculateScale() {
    this.scale = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .scaleBand */ .tiA)().domain(this.categories).range(this.getRange()).paddingInner(1).paddingOuter(0).align(0.5);
    _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.trace("BandAxis axis final categories, range: ", this.categories, this.getRange());
  }
  getTickValues() {
    return this.categories;
  }
  getScaleValue(value) {
    return this.scale(value) ?? this.getRange()[0];
  }
};

// src/diagrams/xychart/chartBuilder/components/axis/linearAxis.ts

var LinearAxis = class extends BaseAxis {
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(this, "LinearAxis");
  }
  constructor(axisConfig, axisThemeConfig, domain, title, textDimensionCalculator) {
    super(axisConfig, title, textDimensionCalculator, axisThemeConfig);
    this.domain = domain;
    this.scale = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .scaleLinear */ .BYU)().domain(this.domain).range(this.getRange());
  }
  getTickValues() {
    return this.scale.ticks();
  }
  recalculateScale() {
    const domain = [...this.domain];
    if (this.axisPosition === "left") {
      domain.reverse();
    }
    this.scale = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .scaleLinear */ .BYU)().domain(domain).range(this.getRange());
  }
  getScaleValue(value) {
    return this.scale(value);
  }
};

// src/diagrams/xychart/chartBuilder/components/axis/index.ts
function getAxis(data, axisConfig, axisThemeConfig, tmpSVGGroup2) {
  const textDimensionCalculator = new TextDimensionCalculatorWithFont(tmpSVGGroup2);
  if (isBandAxisData(data)) {
    return new BandAxis(
      axisConfig,
      axisThemeConfig,
      data.categories,
      data.title,
      textDimensionCalculator
    );
  }
  return new LinearAxis(
    axisConfig,
    axisThemeConfig,
    [data.min, data.max],
    data.title,
    textDimensionCalculator
  );
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getAxis, "getAxis");

// src/diagrams/xychart/chartBuilder/components/chartTitle.ts
var ChartTitle = class {
  constructor(textDimensionCalculator, chartConfig, chartData, chartThemeConfig) {
    this.textDimensionCalculator = textDimensionCalculator;
    this.chartConfig = chartConfig;
    this.chartData = chartData;
    this.chartThemeConfig = chartThemeConfig;
    this.boundingRect = {
      x: 0,
      y: 0,
      width: 0,
      height: 0
    };
    this.showChartTitle = false;
  }
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(this, "ChartTitle");
  }
  setBoundingBoxXY(point) {
    this.boundingRect.x = point.x;
    this.boundingRect.y = point.y;
  }
  calculateSpace(availableSpace) {
    const titleDimension = this.textDimensionCalculator.getMaxDimension(
      [this.chartData.title],
      this.chartConfig.titleFontSize
    );
    const widthRequired = Math.max(titleDimension.width, availableSpace.width);
    const heightRequired = titleDimension.height + 2 * this.chartConfig.titlePadding;
    if (titleDimension.width <= widthRequired && titleDimension.height <= heightRequired && this.chartConfig.showTitle && this.chartData.title) {
      this.boundingRect.width = widthRequired;
      this.boundingRect.height = heightRequired;
      this.showChartTitle = true;
    }
    return {
      width: this.boundingRect.width,
      height: this.boundingRect.height
    };
  }
  getDrawableElements() {
    const drawableElem = [];
    if (this.showChartTitle) {
      drawableElem.push({
        groupTexts: ["chart-title"],
        type: "text",
        data: [
          {
            fontSize: this.chartConfig.titleFontSize,
            text: this.chartData.title,
            verticalPos: "middle",
            horizontalPos: "center",
            x: this.boundingRect.x + this.boundingRect.width / 2,
            y: this.boundingRect.y + this.boundingRect.height / 2,
            fill: this.chartThemeConfig.titleColor,
            rotation: 0
          }
        ]
      });
    }
    return drawableElem;
  }
};
function getChartTitleComponent(chartConfig, chartData, chartThemeConfig, tmpSVGGroup2) {
  const textDimensionCalculator = new TextDimensionCalculatorWithFont(tmpSVGGroup2);
  return new ChartTitle(textDimensionCalculator, chartConfig, chartData, chartThemeConfig);
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getChartTitleComponent, "getChartTitleComponent");

// src/diagrams/xychart/chartBuilder/components/plot/linePlot.ts

var LinePlot = class {
  constructor(plotData, xAxis, yAxis, orientation, plotIndex2) {
    this.plotData = plotData;
    this.xAxis = xAxis;
    this.yAxis = yAxis;
    this.orientation = orientation;
    this.plotIndex = plotIndex2;
  }
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(this, "LinePlot");
  }
  getDrawableElement() {
    const finalData = this.plotData.data.map((d) => [
      this.xAxis.getScaleValue(d[0]),
      this.yAxis.getScaleValue(d[1])
    ]);
    let path;
    if (this.orientation === "horizontal") {
      path = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .line */ .jvg)().y((d) => d[0]).x((d) => d[1])(finalData);
    } else {
      path = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .line */ .jvg)().x((d) => d[0]).y((d) => d[1])(finalData);
    }
    if (!path) {
      return [];
    }
    return [
      {
        groupTexts: ["plot", `line-plot-${this.plotIndex}`],
        type: "path",
        data: [
          {
            path,
            strokeFill: this.plotData.strokeFill,
            strokeWidth: this.plotData.strokeWidth
          }
        ]
      }
    ];
  }
};

// src/diagrams/xychart/chartBuilder/components/plot/barPlot.ts
var BarPlot = class {
  constructor(barData, boundingRect, xAxis, yAxis, orientation, plotIndex2) {
    this.barData = barData;
    this.boundingRect = boundingRect;
    this.xAxis = xAxis;
    this.yAxis = yAxis;
    this.orientation = orientation;
    this.plotIndex = plotIndex2;
  }
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(this, "BarPlot");
  }
  getDrawableElement() {
    const finalData = this.barData.data.map((d) => [
      this.xAxis.getScaleValue(d[0]),
      this.yAxis.getScaleValue(d[1])
    ]);
    const barPaddingPercent = 0.05;
    const barWidth = Math.min(this.xAxis.getAxisOuterPadding() * 2, this.xAxis.getTickDistance()) * (1 - barPaddingPercent);
    const barWidthHalf = barWidth / 2;
    if (this.orientation === "horizontal") {
      return [
        {
          groupTexts: ["plot", `bar-plot-${this.plotIndex}`],
          type: "rect",
          data: finalData.map((data) => ({
            x: this.boundingRect.x,
            y: data[0] - barWidthHalf,
            height: barWidth,
            width: data[1] - this.boundingRect.x,
            fill: this.barData.fill,
            strokeWidth: 0,
            strokeFill: this.barData.fill
          }))
        }
      ];
    }
    return [
      {
        groupTexts: ["plot", `bar-plot-${this.plotIndex}`],
        type: "rect",
        data: finalData.map((data) => ({
          x: data[0] - barWidthHalf,
          y: data[1],
          width: barWidth,
          height: this.boundingRect.y + this.boundingRect.height - data[1],
          fill: this.barData.fill,
          strokeWidth: 0,
          strokeFill: this.barData.fill
        }))
      }
    ];
  }
};

// src/diagrams/xychart/chartBuilder/components/plot/index.ts
var BasePlot = class {
  constructor(chartConfig, chartData, chartThemeConfig) {
    this.chartConfig = chartConfig;
    this.chartData = chartData;
    this.chartThemeConfig = chartThemeConfig;
    this.boundingRect = {
      x: 0,
      y: 0,
      width: 0,
      height: 0
    };
  }
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(this, "BasePlot");
  }
  setAxes(xAxis, yAxis) {
    this.xAxis = xAxis;
    this.yAxis = yAxis;
  }
  setBoundingBoxXY(point) {
    this.boundingRect.x = point.x;
    this.boundingRect.y = point.y;
  }
  calculateSpace(availableSpace) {
    this.boundingRect.width = availableSpace.width;
    this.boundingRect.height = availableSpace.height;
    return {
      width: this.boundingRect.width,
      height: this.boundingRect.height
    };
  }
  getDrawableElements() {
    if (!(this.xAxis && this.yAxis)) {
      throw Error("Axes must be passed to render Plots");
    }
    const drawableElem = [];
    for (const [i, plot] of this.chartData.plots.entries()) {
      switch (plot.type) {
        case "line":
          {
            const linePlot = new LinePlot(
              plot,
              this.xAxis,
              this.yAxis,
              this.chartConfig.chartOrientation,
              i
            );
            drawableElem.push(...linePlot.getDrawableElement());
          }
          break;
        case "bar":
          {
            const barPlot = new BarPlot(
              plot,
              this.boundingRect,
              this.xAxis,
              this.yAxis,
              this.chartConfig.chartOrientation,
              i
            );
            drawableElem.push(...barPlot.getDrawableElement());
          }
          break;
      }
    }
    return drawableElem;
  }
};
function getPlotComponent(chartConfig, chartData, chartThemeConfig) {
  return new BasePlot(chartConfig, chartData, chartThemeConfig);
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getPlotComponent, "getPlotComponent");

// src/diagrams/xychart/chartBuilder/orchestrator.ts
var Orchestrator = class {
  constructor(chartConfig, chartData, chartThemeConfig, tmpSVGGroup2) {
    this.chartConfig = chartConfig;
    this.chartData = chartData;
    this.componentStore = {
      title: getChartTitleComponent(chartConfig, chartData, chartThemeConfig, tmpSVGGroup2),
      plot: getPlotComponent(chartConfig, chartData, chartThemeConfig),
      xAxis: getAxis(
        chartData.xAxis,
        chartConfig.xAxis,
        {
          titleColor: chartThemeConfig.xAxisTitleColor,
          labelColor: chartThemeConfig.xAxisLabelColor,
          tickColor: chartThemeConfig.xAxisTickColor,
          axisLineColor: chartThemeConfig.xAxisLineColor
        },
        tmpSVGGroup2
      ),
      yAxis: getAxis(
        chartData.yAxis,
        chartConfig.yAxis,
        {
          titleColor: chartThemeConfig.yAxisTitleColor,
          labelColor: chartThemeConfig.yAxisLabelColor,
          tickColor: chartThemeConfig.yAxisTickColor,
          axisLineColor: chartThemeConfig.yAxisLineColor
        },
        tmpSVGGroup2
      )
    };
  }
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(this, "Orchestrator");
  }
  calculateVerticalSpace() {
    let availableWidth = this.chartConfig.width;
    let availableHeight = this.chartConfig.height;
    let plotX = 0;
    let plotY = 0;
    let chartWidth = Math.floor(availableWidth * this.chartConfig.plotReservedSpacePercent / 100);
    let chartHeight = Math.floor(
      availableHeight * this.chartConfig.plotReservedSpacePercent / 100
    );
    let spaceUsed = this.componentStore.plot.calculateSpace({
      width: chartWidth,
      height: chartHeight
    });
    availableWidth -= spaceUsed.width;
    availableHeight -= spaceUsed.height;
    spaceUsed = this.componentStore.title.calculateSpace({
      width: this.chartConfig.width,
      height: availableHeight
    });
    plotY = spaceUsed.height;
    availableHeight -= spaceUsed.height;
    this.componentStore.xAxis.setAxisPosition("bottom");
    spaceUsed = this.componentStore.xAxis.calculateSpace({
      width: availableWidth,
      height: availableHeight
    });
    availableHeight -= spaceUsed.height;
    this.componentStore.yAxis.setAxisPosition("left");
    spaceUsed = this.componentStore.yAxis.calculateSpace({
      width: availableWidth,
      height: availableHeight
    });
    plotX = spaceUsed.width;
    availableWidth -= spaceUsed.width;
    if (availableWidth > 0) {
      chartWidth += availableWidth;
      availableWidth = 0;
    }
    if (availableHeight > 0) {
      chartHeight += availableHeight;
      availableHeight = 0;
    }
    this.componentStore.plot.calculateSpace({
      width: chartWidth,
      height: chartHeight
    });
    this.componentStore.plot.setBoundingBoxXY({ x: plotX, y: plotY });
    this.componentStore.xAxis.setRange([plotX, plotX + chartWidth]);
    this.componentStore.xAxis.setBoundingBoxXY({ x: plotX, y: plotY + chartHeight });
    this.componentStore.yAxis.setRange([plotY, plotY + chartHeight]);
    this.componentStore.yAxis.setBoundingBoxXY({ x: 0, y: plotY });
    if (this.chartData.plots.some((p) => isBarPlot(p))) {
      this.componentStore.xAxis.recalculateOuterPaddingToDrawBar();
    }
  }
  calculateHorizontalSpace() {
    let availableWidth = this.chartConfig.width;
    let availableHeight = this.chartConfig.height;
    let titleYEnd = 0;
    let plotX = 0;
    let plotY = 0;
    let chartWidth = Math.floor(availableWidth * this.chartConfig.plotReservedSpacePercent / 100);
    let chartHeight = Math.floor(
      availableHeight * this.chartConfig.plotReservedSpacePercent / 100
    );
    let spaceUsed = this.componentStore.plot.calculateSpace({
      width: chartWidth,
      height: chartHeight
    });
    availableWidth -= spaceUsed.width;
    availableHeight -= spaceUsed.height;
    spaceUsed = this.componentStore.title.calculateSpace({
      width: this.chartConfig.width,
      height: availableHeight
    });
    titleYEnd = spaceUsed.height;
    availableHeight -= spaceUsed.height;
    this.componentStore.xAxis.setAxisPosition("left");
    spaceUsed = this.componentStore.xAxis.calculateSpace({
      width: availableWidth,
      height: availableHeight
    });
    availableWidth -= spaceUsed.width;
    plotX = spaceUsed.width;
    this.componentStore.yAxis.setAxisPosition("top");
    spaceUsed = this.componentStore.yAxis.calculateSpace({
      width: availableWidth,
      height: availableHeight
    });
    availableHeight -= spaceUsed.height;
    plotY = titleYEnd + spaceUsed.height;
    if (availableWidth > 0) {
      chartWidth += availableWidth;
      availableWidth = 0;
    }
    if (availableHeight > 0) {
      chartHeight += availableHeight;
      availableHeight = 0;
    }
    this.componentStore.plot.calculateSpace({
      width: chartWidth,
      height: chartHeight
    });
    this.componentStore.plot.setBoundingBoxXY({ x: plotX, y: plotY });
    this.componentStore.yAxis.setRange([plotX, plotX + chartWidth]);
    this.componentStore.yAxis.setBoundingBoxXY({ x: plotX, y: titleYEnd });
    this.componentStore.xAxis.setRange([plotY, plotY + chartHeight]);
    this.componentStore.xAxis.setBoundingBoxXY({ x: 0, y: plotY });
    if (this.chartData.plots.some((p) => isBarPlot(p))) {
      this.componentStore.xAxis.recalculateOuterPaddingToDrawBar();
    }
  }
  calculateSpace() {
    if (this.chartConfig.chartOrientation === "horizontal") {
      this.calculateHorizontalSpace();
    } else {
      this.calculateVerticalSpace();
    }
  }
  getDrawableElement() {
    this.calculateSpace();
    const drawableElem = [];
    this.componentStore.plot.setAxes(this.componentStore.xAxis, this.componentStore.yAxis);
    for (const component of Object.values(this.componentStore)) {
      drawableElem.push(...component.getDrawableElements());
    }
    return drawableElem;
  }
};

// src/diagrams/xychart/chartBuilder/index.ts
var XYChartBuilder = class {
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(this, "XYChartBuilder");
  }
  static build(config, chartData, chartThemeConfig, tmpSVGGroup2) {
    const orchestrator = new Orchestrator(config, chartData, chartThemeConfig, tmpSVGGroup2);
    return orchestrator.getDrawableElement();
  }
};

// src/diagrams/xychart/xychartDb.ts
var plotIndex = 0;
var tmpSVGGroup;
var xyChartConfig = getChartDefaultConfig();
var xyChartThemeConfig = getChartDefaultThemeConfig();
var xyChartData = getChartDefaultData();
var plotColorPalette = xyChartThemeConfig.plotColorPalette.split(",").map((color) => color.trim());
var hasSetXAxis = false;
var hasSetYAxis = false;
function getChartDefaultThemeConfig() {
  const defaultThemeVariables = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getThemeVariables */ .xN)();
  const config = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig */ .iE)();
  return (0,_chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_1__/* .cleanAndMerge */ .Rb)(defaultThemeVariables.xyChart, config.themeVariables.xyChart);
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getChartDefaultThemeConfig, "getChartDefaultThemeConfig");
function getChartDefaultConfig() {
  const config = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig */ .iE)();
  return (0,_chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_1__/* .cleanAndMerge */ .Rb)(
    _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .defaultConfig_default */ .vZ.xyChart,
    config.xyChart
  );
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getChartDefaultConfig, "getChartDefaultConfig");
function getChartDefaultData() {
  return {
    yAxis: {
      type: "linear",
      title: "",
      min: Infinity,
      max: -Infinity
    },
    xAxis: {
      type: "band",
      title: "",
      categories: []
    },
    title: "",
    plots: []
  };
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getChartDefaultData, "getChartDefaultData");
function textSanitizer(text) {
  const config = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig */ .iE)();
  return (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .sanitizeText */ .oO)(text.trim(), config);
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(textSanitizer, "textSanitizer");
function setTmpSVGG(SVGG) {
  tmpSVGGroup = SVGG;
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(setTmpSVGG, "setTmpSVGG");
function setOrientation(orientation) {
  if (orientation === "horizontal") {
    xyChartConfig.chartOrientation = "horizontal";
  } else {
    xyChartConfig.chartOrientation = "vertical";
  }
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(setOrientation, "setOrientation");
function setXAxisTitle(title) {
  xyChartData.xAxis.title = textSanitizer(title.text);
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(setXAxisTitle, "setXAxisTitle");
function setXAxisRangeData(min, max) {
  xyChartData.xAxis = { type: "linear", title: xyChartData.xAxis.title, min, max };
  hasSetXAxis = true;
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(setXAxisRangeData, "setXAxisRangeData");
function setXAxisBand(categories) {
  xyChartData.xAxis = {
    type: "band",
    title: xyChartData.xAxis.title,
    categories: categories.map((c) => textSanitizer(c.text))
  };
  hasSetXAxis = true;
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(setXAxisBand, "setXAxisBand");
function setYAxisTitle(title) {
  xyChartData.yAxis.title = textSanitizer(title.text);
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(setYAxisTitle, "setYAxisTitle");
function setYAxisRangeData(min, max) {
  xyChartData.yAxis = { type: "linear", title: xyChartData.yAxis.title, min, max };
  hasSetYAxis = true;
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(setYAxisRangeData, "setYAxisRangeData");
function setYAxisRangeFromPlotData(data) {
  const minValue = Math.min(...data);
  const maxValue = Math.max(...data);
  const prevMinValue = isLinearAxisData(xyChartData.yAxis) ? xyChartData.yAxis.min : Infinity;
  const prevMaxValue = isLinearAxisData(xyChartData.yAxis) ? xyChartData.yAxis.max : -Infinity;
  xyChartData.yAxis = {
    type: "linear",
    title: xyChartData.yAxis.title,
    min: Math.min(prevMinValue, minValue),
    max: Math.max(prevMaxValue, maxValue)
  };
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(setYAxisRangeFromPlotData, "setYAxisRangeFromPlotData");
function transformDataWithoutCategory(data) {
  let retData = [];
  if (data.length === 0) {
    return retData;
  }
  if (!hasSetXAxis) {
    const prevMinValue = isLinearAxisData(xyChartData.xAxis) ? xyChartData.xAxis.min : Infinity;
    const prevMaxValue = isLinearAxisData(xyChartData.xAxis) ? xyChartData.xAxis.max : -Infinity;
    setXAxisRangeData(Math.min(prevMinValue, 1), Math.max(prevMaxValue, data.length));
  }
  if (!hasSetYAxis) {
    setYAxisRangeFromPlotData(data);
  }
  if (isBandAxisData(xyChartData.xAxis)) {
    retData = xyChartData.xAxis.categories.map((c, i) => [c, data[i]]);
  }
  if (isLinearAxisData(xyChartData.xAxis)) {
    const min = xyChartData.xAxis.min;
    const max = xyChartData.xAxis.max;
    const step = (max - min) / (data.length - 1);
    const categories = [];
    for (let i = min; i <= max; i += step) {
      categories.push(`${i}`);
    }
    retData = categories.map((c, i) => [c, data[i]]);
  }
  return retData;
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(transformDataWithoutCategory, "transformDataWithoutCategory");
function getPlotColorFromPalette(plotIndex2) {
  return plotColorPalette[plotIndex2 === 0 ? 0 : plotIndex2 % plotColorPalette.length];
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getPlotColorFromPalette, "getPlotColorFromPalette");
function setLineData(title, data) {
  const plotData = transformDataWithoutCategory(data);
  xyChartData.plots.push({
    type: "line",
    strokeFill: getPlotColorFromPalette(plotIndex),
    strokeWidth: 2,
    data: plotData
  });
  plotIndex++;
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(setLineData, "setLineData");
function setBarData(title, data) {
  const plotData = transformDataWithoutCategory(data);
  xyChartData.plots.push({
    type: "bar",
    fill: getPlotColorFromPalette(plotIndex),
    data: plotData
  });
  plotIndex++;
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(setBarData, "setBarData");
function getDrawableElem() {
  if (xyChartData.plots.length === 0) {
    throw Error("No Plot to render, please provide a plot with some data");
  }
  xyChartData.title = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getDiagramTitle */ .Kr)();
  return XYChartBuilder.build(xyChartConfig, xyChartData, xyChartThemeConfig, tmpSVGGroup);
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getDrawableElem, "getDrawableElem");
function getChartThemeConfig() {
  return xyChartThemeConfig;
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getChartThemeConfig, "getChartThemeConfig");
function getChartConfig() {
  return xyChartConfig;
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getChartConfig, "getChartConfig");
var clear2 = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .clear */ .ZH)();
  plotIndex = 0;
  xyChartConfig = getChartDefaultConfig();
  xyChartData = getChartDefaultData();
  xyChartThemeConfig = getChartDefaultThemeConfig();
  plotColorPalette = xyChartThemeConfig.plotColorPalette.split(",").map((color) => color.trim());
  hasSetXAxis = false;
  hasSetYAxis = false;
}, "clear");
var xychartDb_default = {
  getDrawableElem,
  clear: clear2,
  setAccTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccTitle */ .GN,
  getAccTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccTitle */ .eu,
  setDiagramTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setDiagramTitle */ .g2,
  getDiagramTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getDiagramTitle */ .Kr,
  getAccDescription: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccDescription */ .Mx,
  setAccDescription: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccDescription */ .U$,
  setOrientation,
  setXAxisTitle,
  setXAxisRangeData,
  setXAxisBand,
  setYAxisTitle,
  setYAxisRangeData,
  setLineData,
  setBarData,
  setTmpSVGG,
  getChartThemeConfig,
  getChartConfig
};

// src/diagrams/xychart/xychartRenderer.ts
var draw = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((txt, id, _version, diagObj) => {
  const db = diagObj.db;
  const themeConfig = db.getChartThemeConfig();
  const chartConfig = db.getChartConfig();
  function getDominantBaseLine(horizontalPos) {
    return horizontalPos === "top" ? "text-before-edge" : "middle";
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getDominantBaseLine, "getDominantBaseLine");
  function getTextAnchor(verticalPos) {
    return verticalPos === "left" ? "start" : verticalPos === "right" ? "end" : "middle";
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getTextAnchor, "getTextAnchor");
  function getTextTransformation(data) {
    return `translate(${data.x}, ${data.y}) rotate(${data.rotation || 0})`;
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getTextTransformation, "getTextTransformation");
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug("Rendering xychart chart\n" + txt);
  const svg = (0,_chunk_7B677QYD_mjs__WEBPACK_IMPORTED_MODULE_2__/* .selectSvgElement */ .P)(id);
  const group = svg.append("g").attr("class", "main");
  const background = group.append("rect").attr("width", chartConfig.width).attr("height", chartConfig.height).attr("class", "background");
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .configureSvgSize */ .v2)(svg, chartConfig.height, chartConfig.width, true);
  svg.attr("viewBox", `0 0 ${chartConfig.width} ${chartConfig.height}`);
  background.attr("fill", themeConfig.backgroundColor);
  db.setTmpSVGG(svg.append("g").attr("class", "mermaid-tmp-group"));
  const shapes = db.getDrawableElem();
  const groups = {};
  function getGroup(gList) {
    let elem = group;
    let prefix = "";
    for (const [i] of gList.entries()) {
      let parent = group;
      if (i > 0 && groups[prefix]) {
        parent = groups[prefix];
      }
      prefix += gList[i];
      elem = groups[prefix];
      if (!elem) {
        elem = groups[prefix] = parent.append("g").attr("class", gList[i]);
      }
    }
    return elem;
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(getGroup, "getGroup");
  for (const shape of shapes) {
    if (shape.data.length === 0) {
      continue;
    }
    const shapeGroup = getGroup(shape.groupTexts);
    switch (shape.type) {
      case "rect":
        shapeGroup.selectAll("rect").data(shape.data).enter().append("rect").attr("x", (data) => data.x).attr("y", (data) => data.y).attr("width", (data) => data.width).attr("height", (data) => data.height).attr("fill", (data) => data.fill).attr("stroke", (data) => data.strokeFill).attr("stroke-width", (data) => data.strokeWidth);
        break;
      case "text":
        shapeGroup.selectAll("text").data(shape.data).enter().append("text").attr("x", 0).attr("y", 0).attr("fill", (data) => data.fill).attr("font-size", (data) => data.fontSize).attr("dominant-baseline", (data) => getDominantBaseLine(data.verticalPos)).attr("text-anchor", (data) => getTextAnchor(data.horizontalPos)).attr("transform", (data) => getTextTransformation(data)).text((data) => data.text);
        break;
      case "path":
        shapeGroup.selectAll("path").data(shape.data).enter().append("path").attr("d", (data) => data.path).attr("fill", (data) => data.fill ? data.fill : "none").attr("stroke", (data) => data.strokeFill).attr("stroke-width", (data) => data.strokeWidth);
        break;
    }
  }
}, "draw");
var xychartRenderer_default = {
  draw
};

// src/diagrams/xychart/xychartDiagram.ts
var diagram = {
  parser: xychart_default,
  db: xychartDb_default,
  renderer: xychartRenderer_default
};



/***/ })

}]);
//# sourceMappingURL=6491.4ec5e8e76fbff7d9698a.js.map?v=4ec5e8e76fbff7d9698a
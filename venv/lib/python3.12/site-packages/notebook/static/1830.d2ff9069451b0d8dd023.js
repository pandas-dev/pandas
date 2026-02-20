(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[1830],{

/***/ 11830:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(17175);
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(74999);
/* harmony import */ var _braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(7608);
/* harmony import */ var dayjs__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(59395);
/* harmony import */ var dayjs__WEBPACK_IMPORTED_MODULE_4___default = /*#__PURE__*/__webpack_require__.n(dayjs__WEBPACK_IMPORTED_MODULE_4__);
/* harmony import */ var dayjs_plugin_isoWeek_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(62539);
/* harmony import */ var dayjs_plugin_isoWeek_js__WEBPACK_IMPORTED_MODULE_5___default = /*#__PURE__*/__webpack_require__.n(dayjs_plugin_isoWeek_js__WEBPACK_IMPORTED_MODULE_5__);
/* harmony import */ var dayjs_plugin_customParseFormat_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(14881);
/* harmony import */ var dayjs_plugin_customParseFormat_js__WEBPACK_IMPORTED_MODULE_6___default = /*#__PURE__*/__webpack_require__.n(dayjs_plugin_customParseFormat_js__WEBPACK_IMPORTED_MODULE_6__);
/* harmony import */ var dayjs_plugin_advancedFormat_js__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(85747);
/* harmony import */ var dayjs_plugin_advancedFormat_js__WEBPACK_IMPORTED_MODULE_7___default = /*#__PURE__*/__webpack_require__.n(dayjs_plugin_advancedFormat_js__WEBPACK_IMPORTED_MODULE_7__);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(35321);




// src/diagrams/gantt/parser/gantt.jison
var parser = (function() {
  var o = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(k, v, o2, l) {
    for (o2 = o2 || {}, l = k.length; l--; o2[k[l]] = v) ;
    return o2;
  }, "o"), $V0 = [6, 8, 10, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 35, 36, 38, 40], $V1 = [1, 26], $V2 = [1, 27], $V3 = [1, 28], $V4 = [1, 29], $V5 = [1, 30], $V6 = [1, 31], $V7 = [1, 32], $V8 = [1, 33], $V9 = [1, 34], $Va = [1, 9], $Vb = [1, 10], $Vc = [1, 11], $Vd = [1, 12], $Ve = [1, 13], $Vf = [1, 14], $Vg = [1, 15], $Vh = [1, 16], $Vi = [1, 19], $Vj = [1, 20], $Vk = [1, 21], $Vl = [1, 22], $Vm = [1, 23], $Vn = [1, 25], $Vo = [1, 35];
  var parser2 = {
    trace: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function trace() {
    }, "trace"),
    yy: {},
    symbols_: { "error": 2, "start": 3, "gantt": 4, "document": 5, "EOF": 6, "line": 7, "SPACE": 8, "statement": 9, "NL": 10, "weekday": 11, "weekday_monday": 12, "weekday_tuesday": 13, "weekday_wednesday": 14, "weekday_thursday": 15, "weekday_friday": 16, "weekday_saturday": 17, "weekday_sunday": 18, "weekend": 19, "weekend_friday": 20, "weekend_saturday": 21, "dateFormat": 22, "inclusiveEndDates": 23, "topAxis": 24, "axisFormat": 25, "tickInterval": 26, "excludes": 27, "includes": 28, "todayMarker": 29, "title": 30, "acc_title": 31, "acc_title_value": 32, "acc_descr": 33, "acc_descr_value": 34, "acc_descr_multiline_value": 35, "section": 36, "clickStatement": 37, "taskTxt": 38, "taskData": 39, "click": 40, "callbackname": 41, "callbackargs": 42, "href": 43, "clickStatementDebug": 44, "$accept": 0, "$end": 1 },
    terminals_: { 2: "error", 4: "gantt", 6: "EOF", 8: "SPACE", 10: "NL", 12: "weekday_monday", 13: "weekday_tuesday", 14: "weekday_wednesday", 15: "weekday_thursday", 16: "weekday_friday", 17: "weekday_saturday", 18: "weekday_sunday", 20: "weekend_friday", 21: "weekend_saturday", 22: "dateFormat", 23: "inclusiveEndDates", 24: "topAxis", 25: "axisFormat", 26: "tickInterval", 27: "excludes", 28: "includes", 29: "todayMarker", 30: "title", 31: "acc_title", 32: "acc_title_value", 33: "acc_descr", 34: "acc_descr_value", 35: "acc_descr_multiline_value", 36: "section", 38: "taskTxt", 39: "taskData", 40: "click", 41: "callbackname", 42: "callbackargs", 43: "href" },
    productions_: [0, [3, 3], [5, 0], [5, 2], [7, 2], [7, 1], [7, 1], [7, 1], [11, 1], [11, 1], [11, 1], [11, 1], [11, 1], [11, 1], [11, 1], [19, 1], [19, 1], [9, 1], [9, 1], [9, 1], [9, 1], [9, 1], [9, 1], [9, 1], [9, 1], [9, 1], [9, 1], [9, 1], [9, 2], [9, 2], [9, 1], [9, 1], [9, 1], [9, 2], [37, 2], [37, 3], [37, 3], [37, 4], [37, 3], [37, 4], [37, 2], [44, 2], [44, 3], [44, 3], [44, 4], [44, 3], [44, 4], [44, 2]],
    performAction: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function anonymous(yytext, yyleng, yylineno, yy, yystate, $$, _$) {
      var $0 = $$.length - 1;
      switch (yystate) {
        case 1:
          return $$[$0 - 1];
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
          yy.setWeekday("monday");
          break;
        case 9:
          yy.setWeekday("tuesday");
          break;
        case 10:
          yy.setWeekday("wednesday");
          break;
        case 11:
          yy.setWeekday("thursday");
          break;
        case 12:
          yy.setWeekday("friday");
          break;
        case 13:
          yy.setWeekday("saturday");
          break;
        case 14:
          yy.setWeekday("sunday");
          break;
        case 15:
          yy.setWeekend("friday");
          break;
        case 16:
          yy.setWeekend("saturday");
          break;
        case 17:
          yy.setDateFormat($$[$0].substr(11));
          this.$ = $$[$0].substr(11);
          break;
        case 18:
          yy.enableInclusiveEndDates();
          this.$ = $$[$0].substr(18);
          break;
        case 19:
          yy.TopAxis();
          this.$ = $$[$0].substr(8);
          break;
        case 20:
          yy.setAxisFormat($$[$0].substr(11));
          this.$ = $$[$0].substr(11);
          break;
        case 21:
          yy.setTickInterval($$[$0].substr(13));
          this.$ = $$[$0].substr(13);
          break;
        case 22:
          yy.setExcludes($$[$0].substr(9));
          this.$ = $$[$0].substr(9);
          break;
        case 23:
          yy.setIncludes($$[$0].substr(9));
          this.$ = $$[$0].substr(9);
          break;
        case 24:
          yy.setTodayMarker($$[$0].substr(12));
          this.$ = $$[$0].substr(12);
          break;
        case 27:
          yy.setDiagramTitle($$[$0].substr(6));
          this.$ = $$[$0].substr(6);
          break;
        case 28:
          this.$ = $$[$0].trim();
          yy.setAccTitle(this.$);
          break;
        case 29:
        case 30:
          this.$ = $$[$0].trim();
          yy.setAccDescription(this.$);
          break;
        case 31:
          yy.addSection($$[$0].substr(8));
          this.$ = $$[$0].substr(8);
          break;
        case 33:
          yy.addTask($$[$0 - 1], $$[$0]);
          this.$ = "task";
          break;
        case 34:
          this.$ = $$[$0 - 1];
          yy.setClickEvent($$[$0 - 1], $$[$0], null);
          break;
        case 35:
          this.$ = $$[$0 - 2];
          yy.setClickEvent($$[$0 - 2], $$[$0 - 1], $$[$0]);
          break;
        case 36:
          this.$ = $$[$0 - 2];
          yy.setClickEvent($$[$0 - 2], $$[$0 - 1], null);
          yy.setLink($$[$0 - 2], $$[$0]);
          break;
        case 37:
          this.$ = $$[$0 - 3];
          yy.setClickEvent($$[$0 - 3], $$[$0 - 2], $$[$0 - 1]);
          yy.setLink($$[$0 - 3], $$[$0]);
          break;
        case 38:
          this.$ = $$[$0 - 2];
          yy.setClickEvent($$[$0 - 2], $$[$0], null);
          yy.setLink($$[$0 - 2], $$[$0 - 1]);
          break;
        case 39:
          this.$ = $$[$0 - 3];
          yy.setClickEvent($$[$0 - 3], $$[$0 - 1], $$[$0]);
          yy.setLink($$[$0 - 3], $$[$0 - 2]);
          break;
        case 40:
          this.$ = $$[$0 - 1];
          yy.setLink($$[$0 - 1], $$[$0]);
          break;
        case 41:
        case 47:
          this.$ = $$[$0 - 1] + " " + $$[$0];
          break;
        case 42:
        case 43:
        case 45:
          this.$ = $$[$0 - 2] + " " + $$[$0 - 1] + " " + $$[$0];
          break;
        case 44:
        case 46:
          this.$ = $$[$0 - 3] + " " + $$[$0 - 2] + " " + $$[$0 - 1] + " " + $$[$0];
          break;
      }
    }, "anonymous"),
    table: [{ 3: 1, 4: [1, 2] }, { 1: [3] }, o($V0, [2, 2], { 5: 3 }), { 6: [1, 4], 7: 5, 8: [1, 6], 9: 7, 10: [1, 8], 11: 17, 12: $V1, 13: $V2, 14: $V3, 15: $V4, 16: $V5, 17: $V6, 18: $V7, 19: 18, 20: $V8, 21: $V9, 22: $Va, 23: $Vb, 24: $Vc, 25: $Vd, 26: $Ve, 27: $Vf, 28: $Vg, 29: $Vh, 30: $Vi, 31: $Vj, 33: $Vk, 35: $Vl, 36: $Vm, 37: 24, 38: $Vn, 40: $Vo }, o($V0, [2, 7], { 1: [2, 1] }), o($V0, [2, 3]), { 9: 36, 11: 17, 12: $V1, 13: $V2, 14: $V3, 15: $V4, 16: $V5, 17: $V6, 18: $V7, 19: 18, 20: $V8, 21: $V9, 22: $Va, 23: $Vb, 24: $Vc, 25: $Vd, 26: $Ve, 27: $Vf, 28: $Vg, 29: $Vh, 30: $Vi, 31: $Vj, 33: $Vk, 35: $Vl, 36: $Vm, 37: 24, 38: $Vn, 40: $Vo }, o($V0, [2, 5]), o($V0, [2, 6]), o($V0, [2, 17]), o($V0, [2, 18]), o($V0, [2, 19]), o($V0, [2, 20]), o($V0, [2, 21]), o($V0, [2, 22]), o($V0, [2, 23]), o($V0, [2, 24]), o($V0, [2, 25]), o($V0, [2, 26]), o($V0, [2, 27]), { 32: [1, 37] }, { 34: [1, 38] }, o($V0, [2, 30]), o($V0, [2, 31]), o($V0, [2, 32]), { 39: [1, 39] }, o($V0, [2, 8]), o($V0, [2, 9]), o($V0, [2, 10]), o($V0, [2, 11]), o($V0, [2, 12]), o($V0, [2, 13]), o($V0, [2, 14]), o($V0, [2, 15]), o($V0, [2, 16]), { 41: [1, 40], 43: [1, 41] }, o($V0, [2, 4]), o($V0, [2, 28]), o($V0, [2, 29]), o($V0, [2, 33]), o($V0, [2, 34], { 42: [1, 42], 43: [1, 43] }), o($V0, [2, 40], { 41: [1, 44] }), o($V0, [2, 35], { 43: [1, 45] }), o($V0, [2, 36]), o($V0, [2, 38], { 42: [1, 46] }), o($V0, [2, 37]), o($V0, [2, 39])],
    defaultActions: {},
    parseError: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function parseError(str, hash) {
      if (hash.recoverable) {
        this.trace(str);
      } else {
        var error = new Error(str);
        error.hash = hash;
        throw error;
      }
    }, "parseError"),
    parse: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function parse(input) {
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
      (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(popStack, "popStack");
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
      (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(lex, "lex");
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
      parseError: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function parseError(str, hash) {
        if (this.yy.parser) {
          this.yy.parser.parseError(str, hash);
        } else {
          throw new Error(str);
        }
      }, "parseError"),
      // resets the lexer, sets new input
      setInput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(input, yy) {
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
      input: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
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
      unput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(ch) {
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
      more: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
        this._more = true;
        return this;
      }, "more"),
      // When called from action, signals the lexer that this rule fails to match the input, so the next matching rule (regex) should be tested instead.
      reject: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
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
      less: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(n) {
        this.unput(this.match.slice(n));
      }, "less"),
      // displays already matched input, i.e. for error messages
      pastInput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
        var past = this.matched.substr(0, this.matched.length - this.match.length);
        return (past.length > 20 ? "..." : "") + past.substr(-20).replace(/\n/g, "");
      }, "pastInput"),
      // displays upcoming input, i.e. for error messages
      upcomingInput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
        var next = this.match;
        if (next.length < 20) {
          next += this._input.substr(0, 20 - next.length);
        }
        return (next.substr(0, 20) + (next.length > 20 ? "..." : "")).replace(/\n/g, "");
      }, "upcomingInput"),
      // displays the character position where the lexing error occurred, i.e. for error messages
      showPosition: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
        var pre = this.pastInput();
        var c = new Array(pre.length + 1).join("-");
        return pre + this.upcomingInput() + "\n" + c + "^";
      }, "showPosition"),
      // test the lexed token: return FALSE when not a match, otherwise return token
      test_match: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(match, indexed_rule) {
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
      next: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
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
      lex: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function lex() {
        var r = this.next();
        if (r) {
          return r;
        } else {
          return this.lex();
        }
      }, "lex"),
      // activates a new lexer condition state (pushes the new lexer condition state onto the condition stack)
      begin: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function begin(condition) {
        this.conditionStack.push(condition);
      }, "begin"),
      // pop the previously active lexer condition state off the condition stack
      popState: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function popState() {
        var n = this.conditionStack.length - 1;
        if (n > 0) {
          return this.conditionStack.pop();
        } else {
          return this.conditionStack[0];
        }
      }, "popState"),
      // produce the lexer rule set which is active for the currently active lexer condition state
      _currentRules: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function _currentRules() {
        if (this.conditionStack.length && this.conditionStack[this.conditionStack.length - 1]) {
          return this.conditions[this.conditionStack[this.conditionStack.length - 1]].rules;
        } else {
          return this.conditions["INITIAL"].rules;
        }
      }, "_currentRules"),
      // return the currently active lexer condition state; when an index argument is provided it produces the N-th previous condition state, if available
      topState: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function topState(n) {
        n = this.conditionStack.length - 1 - Math.abs(n || 0);
        if (n >= 0) {
          return this.conditionStack[n];
        } else {
          return "INITIAL";
        }
      }, "topState"),
      // alias for begin(condition)
      pushState: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function pushState(condition) {
        this.begin(condition);
      }, "pushState"),
      // return the number of states currently on the stack
      stateStackSize: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function stateStackSize() {
        return this.conditionStack.length;
      }, "stateStackSize"),
      options: { "case-insensitive": true },
      performAction: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function anonymous(yy, yy_, $avoiding_name_collisions, YY_START) {
        var YYSTATE = YY_START;
        switch ($avoiding_name_collisions) {
          case 0:
            this.begin("open_directive");
            return "open_directive";
            break;
          case 1:
            this.begin("acc_title");
            return 31;
            break;
          case 2:
            this.popState();
            return "acc_title_value";
            break;
          case 3:
            this.begin("acc_descr");
            return 33;
            break;
          case 4:
            this.popState();
            return "acc_descr_value";
            break;
          case 5:
            this.begin("acc_descr_multiline");
            break;
          case 6:
            this.popState();
            break;
          case 7:
            return "acc_descr_multiline_value";
            break;
          case 8:
            break;
          case 9:
            break;
          case 10:
            break;
          case 11:
            return 10;
            break;
          case 12:
            break;
          case 13:
            break;
          case 14:
            this.begin("href");
            break;
          case 15:
            this.popState();
            break;
          case 16:
            return 43;
            break;
          case 17:
            this.begin("callbackname");
            break;
          case 18:
            this.popState();
            break;
          case 19:
            this.popState();
            this.begin("callbackargs");
            break;
          case 20:
            return 41;
            break;
          case 21:
            this.popState();
            break;
          case 22:
            return 42;
            break;
          case 23:
            this.begin("click");
            break;
          case 24:
            this.popState();
            break;
          case 25:
            return 40;
            break;
          case 26:
            return 4;
            break;
          case 27:
            return 22;
            break;
          case 28:
            return 23;
            break;
          case 29:
            return 24;
            break;
          case 30:
            return 25;
            break;
          case 31:
            return 26;
            break;
          case 32:
            return 28;
            break;
          case 33:
            return 27;
            break;
          case 34:
            return 29;
            break;
          case 35:
            return 12;
            break;
          case 36:
            return 13;
            break;
          case 37:
            return 14;
            break;
          case 38:
            return 15;
            break;
          case 39:
            return 16;
            break;
          case 40:
            return 17;
            break;
          case 41:
            return 18;
            break;
          case 42:
            return 20;
            break;
          case 43:
            return 21;
            break;
          case 44:
            return "date";
            break;
          case 45:
            return 30;
            break;
          case 46:
            return "accDescription";
            break;
          case 47:
            return 36;
            break;
          case 48:
            return 38;
            break;
          case 49:
            return 39;
            break;
          case 50:
            return ":";
            break;
          case 51:
            return 6;
            break;
          case 52:
            return "INVALID";
            break;
        }
      }, "anonymous"),
      rules: [/^(?:%%\{)/i, /^(?:accTitle\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*\{\s*)/i, /^(?:[\}])/i, /^(?:[^\}]*)/i, /^(?:%%(?!\{)*[^\n]*)/i, /^(?:[^\}]%%*[^\n]*)/i, /^(?:%%*[^\n]*[\n]*)/i, /^(?:[\n]+)/i, /^(?:\s+)/i, /^(?:%[^\n]*)/i, /^(?:href[\s]+["])/i, /^(?:["])/i, /^(?:[^"]*)/i, /^(?:call[\s]+)/i, /^(?:\([\s]*\))/i, /^(?:\()/i, /^(?:[^(]*)/i, /^(?:\))/i, /^(?:[^)]*)/i, /^(?:click[\s]+)/i, /^(?:[\s\n])/i, /^(?:[^\s\n]*)/i, /^(?:gantt\b)/i, /^(?:dateFormat\s[^#\n;]+)/i, /^(?:inclusiveEndDates\b)/i, /^(?:topAxis\b)/i, /^(?:axisFormat\s[^#\n;]+)/i, /^(?:tickInterval\s[^#\n;]+)/i, /^(?:includes\s[^#\n;]+)/i, /^(?:excludes\s[^#\n;]+)/i, /^(?:todayMarker\s[^\n;]+)/i, /^(?:weekday\s+monday\b)/i, /^(?:weekday\s+tuesday\b)/i, /^(?:weekday\s+wednesday\b)/i, /^(?:weekday\s+thursday\b)/i, /^(?:weekday\s+friday\b)/i, /^(?:weekday\s+saturday\b)/i, /^(?:weekday\s+sunday\b)/i, /^(?:weekend\s+friday\b)/i, /^(?:weekend\s+saturday\b)/i, /^(?:\d\d\d\d-\d\d-\d\d\b)/i, /^(?:title\s[^\n]+)/i, /^(?:accDescription\s[^#\n;]+)/i, /^(?:section\s[^\n]+)/i, /^(?:[^:\n]+)/i, /^(?::[^#\n;]+)/i, /^(?::)/i, /^(?:$)/i, /^(?:.)/i],
      conditions: { "acc_descr_multiline": { "rules": [6, 7], "inclusive": false }, "acc_descr": { "rules": [4], "inclusive": false }, "acc_title": { "rules": [2], "inclusive": false }, "callbackargs": { "rules": [21, 22], "inclusive": false }, "callbackname": { "rules": [18, 19, 20], "inclusive": false }, "href": { "rules": [15, 16], "inclusive": false }, "click": { "rules": [24, 25], "inclusive": false }, "INITIAL": { "rules": [0, 1, 3, 5, 8, 9, 10, 11, 12, 13, 14, 17, 23, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52], "inclusive": true } }
    };
    return lexer2;
  })();
  parser2.lexer = lexer;
  function Parser() {
    this.yy = {};
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(Parser, "Parser");
  Parser.prototype = parser2;
  parser2.Parser = Parser;
  return new Parser();
})();
parser.parser = parser;
var gantt_default = parser;

// src/diagrams/gantt/ganttDb.js





dayjs__WEBPACK_IMPORTED_MODULE_4___default().extend((dayjs_plugin_isoWeek_js__WEBPACK_IMPORTED_MODULE_5___default()));
dayjs__WEBPACK_IMPORTED_MODULE_4___default().extend((dayjs_plugin_customParseFormat_js__WEBPACK_IMPORTED_MODULE_6___default()));
dayjs__WEBPACK_IMPORTED_MODULE_4___default().extend((dayjs_plugin_advancedFormat_js__WEBPACK_IMPORTED_MODULE_7___default()));
var WEEKEND_START_DAY = { friday: 5, saturday: 6 };
var dateFormat = "";
var axisFormat = "";
var tickInterval = void 0;
var todayMarker = "";
var includes = [];
var excludes = [];
var links = /* @__PURE__ */ new Map();
var sections = [];
var tasks = [];
var currentSection = "";
var displayMode = "";
var tags = ["active", "done", "crit", "milestone", "vert"];
var funs = [];
var inclusiveEndDates = false;
var topAxis = false;
var weekday = "sunday";
var weekend = "saturday";
var lastOrder = 0;
var clear2 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  sections = [];
  tasks = [];
  currentSection = "";
  funs = [];
  taskCnt = 0;
  lastTask = void 0;
  lastTaskID = void 0;
  rawTasks = [];
  dateFormat = "";
  axisFormat = "";
  displayMode = "";
  tickInterval = void 0;
  todayMarker = "";
  includes = [];
  excludes = [];
  inclusiveEndDates = false;
  topAxis = false;
  lastOrder = 0;
  links = /* @__PURE__ */ new Map();
  (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .clear */ .ZH)();
  weekday = "sunday";
  weekend = "saturday";
}, "clear");
var setAxisFormat = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(txt) {
  axisFormat = txt;
}, "setAxisFormat");
var getAxisFormat = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  return axisFormat;
}, "getAxisFormat");
var setTickInterval = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(txt) {
  tickInterval = txt;
}, "setTickInterval");
var getTickInterval = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  return tickInterval;
}, "getTickInterval");
var setTodayMarker = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(txt) {
  todayMarker = txt;
}, "setTodayMarker");
var getTodayMarker = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  return todayMarker;
}, "getTodayMarker");
var setDateFormat = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(txt) {
  dateFormat = txt;
}, "setDateFormat");
var enableInclusiveEndDates = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  inclusiveEndDates = true;
}, "enableInclusiveEndDates");
var endDatesAreInclusive = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  return inclusiveEndDates;
}, "endDatesAreInclusive");
var enableTopAxis = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  topAxis = true;
}, "enableTopAxis");
var topAxisEnabled = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  return topAxis;
}, "topAxisEnabled");
var setDisplayMode = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(txt) {
  displayMode = txt;
}, "setDisplayMode");
var getDisplayMode = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  return displayMode;
}, "getDisplayMode");
var getDateFormat = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  return dateFormat;
}, "getDateFormat");
var setIncludes = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(txt) {
  includes = txt.toLowerCase().split(/[\s,]+/);
}, "setIncludes");
var getIncludes = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  return includes;
}, "getIncludes");
var setExcludes = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(txt) {
  excludes = txt.toLowerCase().split(/[\s,]+/);
}, "setExcludes");
var getExcludes = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  return excludes;
}, "getExcludes");
var getLinks = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  return links;
}, "getLinks");
var addSection = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(txt) {
  currentSection = txt;
  sections.push(txt);
}, "addSection");
var getSections = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  return sections;
}, "getSections");
var getTasks = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  let allItemsProcessed = compileTasks();
  const maxDepth = 10;
  let iterationCount = 0;
  while (!allItemsProcessed && iterationCount < maxDepth) {
    allItemsProcessed = compileTasks();
    iterationCount++;
  }
  tasks = rawTasks;
  return tasks;
}, "getTasks");
var isInvalidDate = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(date, dateFormat2, excludes2, includes2) {
  const formattedDate = date.format(dateFormat2.trim());
  const dateOnly = date.format("YYYY-MM-DD");
  if (includes2.includes(formattedDate) || includes2.includes(dateOnly)) {
    return false;
  }
  if (excludes2.includes("weekends") && (date.isoWeekday() === WEEKEND_START_DAY[weekend] || date.isoWeekday() === WEEKEND_START_DAY[weekend] + 1)) {
    return true;
  }
  if (excludes2.includes(date.format("dddd").toLowerCase())) {
    return true;
  }
  return excludes2.includes(formattedDate) || excludes2.includes(dateOnly);
}, "isInvalidDate");
var setWeekday = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(txt) {
  weekday = txt;
}, "setWeekday");
var getWeekday = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  return weekday;
}, "getWeekday");
var setWeekend = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(startDay) {
  weekend = startDay;
}, "setWeekend");
var checkTaskDates = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(task, dateFormat2, excludes2, includes2) {
  if (!excludes2.length || task.manualEndTime) {
    return;
  }
  let startTime;
  if (task.startTime instanceof Date) {
    startTime = dayjs__WEBPACK_IMPORTED_MODULE_4___default()(task.startTime);
  } else {
    startTime = dayjs__WEBPACK_IMPORTED_MODULE_4___default()(task.startTime, dateFormat2, true);
  }
  startTime = startTime.add(1, "d");
  let originalEndTime;
  if (task.endTime instanceof Date) {
    originalEndTime = dayjs__WEBPACK_IMPORTED_MODULE_4___default()(task.endTime);
  } else {
    originalEndTime = dayjs__WEBPACK_IMPORTED_MODULE_4___default()(task.endTime, dateFormat2, true);
  }
  const [fixedEndTime, renderEndTime] = fixTaskDates(
    startTime,
    originalEndTime,
    dateFormat2,
    excludes2,
    includes2
  );
  task.endTime = fixedEndTime.toDate();
  task.renderEndTime = renderEndTime;
}, "checkTaskDates");
var fixTaskDates = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(startTime, endTime, dateFormat2, excludes2, includes2) {
  let invalid = false;
  let renderEndTime = null;
  while (startTime <= endTime) {
    if (!invalid) {
      renderEndTime = endTime.toDate();
    }
    invalid = isInvalidDate(startTime, dateFormat2, excludes2, includes2);
    if (invalid) {
      endTime = endTime.add(1, "d");
    }
    startTime = startTime.add(1, "d");
  }
  return [endTime, renderEndTime];
}, "fixTaskDates");
var getStartDate = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(prevTime, dateFormat2, str) {
  str = str.trim();
  if ((dateFormat2.trim() === "x" || dateFormat2.trim() === "X") && /^\d+$/.test(str)) {
    return new Date(Number(str));
  }
  const afterRePattern = /^after\s+(?<ids>[\d\w- ]+)/;
  const afterStatement = afterRePattern.exec(str);
  if (afterStatement !== null) {
    let latestTask = null;
    for (const id of afterStatement.groups.ids.split(" ")) {
      let task = findTaskById(id);
      if (task !== void 0 && (!latestTask || task.endTime > latestTask.endTime)) {
        latestTask = task;
      }
    }
    if (latestTask) {
      return latestTask.endTime;
    }
    const today = /* @__PURE__ */ new Date();
    today.setHours(0, 0, 0, 0);
    return today;
  }
  let mDate = dayjs__WEBPACK_IMPORTED_MODULE_4___default()(str, dateFormat2.trim(), true);
  if (mDate.isValid()) {
    return mDate.toDate();
  } else {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .log */ .cM.debug("Invalid date:" + str);
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .log */ .cM.debug("With date format:" + dateFormat2.trim());
    const d = new Date(str);
    if (d === void 0 || isNaN(d.getTime()) || // WebKit browsers can mis-parse invalid dates to be ridiculously
    // huge numbers, e.g. new Date('202304') gets parsed as January 1, 202304.
    // This can cause virtually infinite loops while rendering, so for the
    // purposes of Gantt charts we'll just treat any date beyond 10,000 AD/BC as
    // invalid.
    d.getFullYear() < -1e4 || d.getFullYear() > 1e4) {
      throw new Error("Invalid date:" + str);
    }
    return d;
  }
}, "getStartDate");
var parseDuration = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(str) {
  const statement = /^(\d+(?:\.\d+)?)([Mdhmswy]|ms)$/.exec(str.trim());
  if (statement !== null) {
    return [Number.parseFloat(statement[1]), statement[2]];
  }
  return [NaN, "ms"];
}, "parseDuration");
var getEndDate = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(prevTime, dateFormat2, str, inclusive = false) {
  str = str.trim();
  const untilRePattern = /^until\s+(?<ids>[\d\w- ]+)/;
  const untilStatement = untilRePattern.exec(str);
  if (untilStatement !== null) {
    let earliestTask = null;
    for (const id of untilStatement.groups.ids.split(" ")) {
      let task = findTaskById(id);
      if (task !== void 0 && (!earliestTask || task.startTime < earliestTask.startTime)) {
        earliestTask = task;
      }
    }
    if (earliestTask) {
      return earliestTask.startTime;
    }
    const today = /* @__PURE__ */ new Date();
    today.setHours(0, 0, 0, 0);
    return today;
  }
  let parsedDate = dayjs__WEBPACK_IMPORTED_MODULE_4___default()(str, dateFormat2.trim(), true);
  if (parsedDate.isValid()) {
    if (inclusive) {
      parsedDate = parsedDate.add(1, "d");
    }
    return parsedDate.toDate();
  }
  let endTime = dayjs__WEBPACK_IMPORTED_MODULE_4___default()(prevTime);
  const [durationValue, durationUnit] = parseDuration(str);
  if (!Number.isNaN(durationValue)) {
    const newEndTime = endTime.add(durationValue, durationUnit);
    if (newEndTime.isValid()) {
      endTime = newEndTime;
    }
  }
  return endTime.toDate();
}, "getEndDate");
var taskCnt = 0;
var parseId = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(idStr) {
  if (idStr === void 0) {
    taskCnt = taskCnt + 1;
    return "task" + taskCnt;
  }
  return idStr;
}, "parseId");
var compileData = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(prevTask, dataStr) {
  let ds;
  if (dataStr.substr(0, 1) === ":") {
    ds = dataStr.substr(1, dataStr.length);
  } else {
    ds = dataStr;
  }
  const data = ds.split(",");
  const task = {};
  getTaskTags(data, task, tags);
  for (let i = 0; i < data.length; i++) {
    data[i] = data[i].trim();
  }
  let endTimeData = "";
  switch (data.length) {
    case 1:
      task.id = parseId();
      task.startTime = prevTask.endTime;
      endTimeData = data[0];
      break;
    case 2:
      task.id = parseId();
      task.startTime = getStartDate(void 0, dateFormat, data[0]);
      endTimeData = data[1];
      break;
    case 3:
      task.id = parseId(data[0]);
      task.startTime = getStartDate(void 0, dateFormat, data[1]);
      endTimeData = data[2];
      break;
    default:
  }
  if (endTimeData) {
    task.endTime = getEndDate(task.startTime, dateFormat, endTimeData, inclusiveEndDates);
    task.manualEndTime = dayjs__WEBPACK_IMPORTED_MODULE_4___default()(endTimeData, "YYYY-MM-DD", true).isValid();
    checkTaskDates(task, dateFormat, excludes, includes);
  }
  return task;
}, "compileData");
var parseData = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(prevTaskId, dataStr) {
  let ds;
  if (dataStr.substr(0, 1) === ":") {
    ds = dataStr.substr(1, dataStr.length);
  } else {
    ds = dataStr;
  }
  const data = ds.split(",");
  const task = {};
  getTaskTags(data, task, tags);
  for (let i = 0; i < data.length; i++) {
    data[i] = data[i].trim();
  }
  switch (data.length) {
    case 1:
      task.id = parseId();
      task.startTime = {
        type: "prevTaskEnd",
        id: prevTaskId
      };
      task.endTime = {
        data: data[0]
      };
      break;
    case 2:
      task.id = parseId();
      task.startTime = {
        type: "getStartDate",
        startData: data[0]
      };
      task.endTime = {
        data: data[1]
      };
      break;
    case 3:
      task.id = parseId(data[0]);
      task.startTime = {
        type: "getStartDate",
        startData: data[1]
      };
      task.endTime = {
        data: data[2]
      };
      break;
    default:
  }
  return task;
}, "parseData");
var lastTask;
var lastTaskID;
var rawTasks = [];
var taskDb = {};
var addTask = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(descr, data) {
  const rawTask = {
    section: currentSection,
    type: currentSection,
    processed: false,
    manualEndTime: false,
    renderEndTime: null,
    raw: { data },
    task: descr,
    classes: []
  };
  const taskInfo = parseData(lastTaskID, data);
  rawTask.raw.startTime = taskInfo.startTime;
  rawTask.raw.endTime = taskInfo.endTime;
  rawTask.id = taskInfo.id;
  rawTask.prevTaskId = lastTaskID;
  rawTask.active = taskInfo.active;
  rawTask.done = taskInfo.done;
  rawTask.crit = taskInfo.crit;
  rawTask.milestone = taskInfo.milestone;
  rawTask.vert = taskInfo.vert;
  rawTask.order = lastOrder;
  lastOrder++;
  const pos = rawTasks.push(rawTask);
  lastTaskID = rawTask.id;
  taskDb[rawTask.id] = pos - 1;
}, "addTask");
var findTaskById = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(id) {
  const pos = taskDb[id];
  return rawTasks[pos];
}, "findTaskById");
var addTaskOrg = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(descr, data) {
  const newTask = {
    section: currentSection,
    type: currentSection,
    description: descr,
    task: descr,
    classes: []
  };
  const taskInfo = compileData(lastTask, data);
  newTask.startTime = taskInfo.startTime;
  newTask.endTime = taskInfo.endTime;
  newTask.id = taskInfo.id;
  newTask.active = taskInfo.active;
  newTask.done = taskInfo.done;
  newTask.crit = taskInfo.crit;
  newTask.milestone = taskInfo.milestone;
  newTask.vert = taskInfo.vert;
  lastTask = newTask;
  tasks.push(newTask);
}, "addTaskOrg");
var compileTasks = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  const compileTask = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(pos) {
    const task = rawTasks[pos];
    let startTime = "";
    switch (rawTasks[pos].raw.startTime.type) {
      case "prevTaskEnd": {
        const prevTask = findTaskById(task.prevTaskId);
        task.startTime = prevTask.endTime;
        break;
      }
      case "getStartDate":
        startTime = getStartDate(void 0, dateFormat, rawTasks[pos].raw.startTime.startData);
        if (startTime) {
          rawTasks[pos].startTime = startTime;
        }
        break;
    }
    if (rawTasks[pos].startTime) {
      rawTasks[pos].endTime = getEndDate(
        rawTasks[pos].startTime,
        dateFormat,
        rawTasks[pos].raw.endTime.data,
        inclusiveEndDates
      );
      if (rawTasks[pos].endTime) {
        rawTasks[pos].processed = true;
        rawTasks[pos].manualEndTime = dayjs__WEBPACK_IMPORTED_MODULE_4___default()(
          rawTasks[pos].raw.endTime.data,
          "YYYY-MM-DD",
          true
        ).isValid();
        checkTaskDates(rawTasks[pos], dateFormat, excludes, includes);
      }
    }
    return rawTasks[pos].processed;
  }, "compileTask");
  let allProcessed = true;
  for (const [i, rawTask] of rawTasks.entries()) {
    compileTask(i);
    allProcessed = allProcessed && rawTask.processed;
  }
  return allProcessed;
}, "compileTasks");
var setLink = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(ids, _linkStr) {
  let linkStr = _linkStr;
  if ((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .getConfig2 */ .nV)().securityLevel !== "loose") {
    linkStr = (0,_braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_3__/* .sanitizeUrl */ .N)(_linkStr);
  }
  ids.split(",").forEach(function(id) {
    let rawTask = findTaskById(id);
    if (rawTask !== void 0) {
      pushFun(id, () => {
        window.open(linkStr, "_self");
      });
      links.set(id, linkStr);
    }
  });
  setClass(ids, "clickable");
}, "setLink");
var setClass = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(ids, className) {
  ids.split(",").forEach(function(id) {
    let rawTask = findTaskById(id);
    if (rawTask !== void 0) {
      rawTask.classes.push(className);
    }
  });
}, "setClass");
var setClickFun = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(id, functionName, functionArgs) {
  if ((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .getConfig2 */ .nV)().securityLevel !== "loose") {
    return;
  }
  if (functionName === void 0) {
    return;
  }
  let argList = [];
  if (typeof functionArgs === "string") {
    argList = functionArgs.split(/,(?=(?:(?:[^"]*"){2})*[^"]*$)/);
    for (let i = 0; i < argList.length; i++) {
      let item = argList[i].trim();
      if (item.startsWith('"') && item.endsWith('"')) {
        item = item.substr(1, item.length - 2);
      }
      argList[i] = item;
    }
  }
  if (argList.length === 0) {
    argList.push(id);
  }
  let rawTask = findTaskById(id);
  if (rawTask !== void 0) {
    pushFun(id, () => {
      _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_0__/* .utils_default */ .w8.runFunc(functionName, ...argList);
    });
  }
}, "setClickFun");
var pushFun = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(id, callbackFunction) {
  funs.push(
    function() {
      const elem = document.querySelector(`[id="${id}"]`);
      if (elem !== null) {
        elem.addEventListener("click", function() {
          callbackFunction();
        });
      }
    },
    function() {
      const elem = document.querySelector(`[id="${id}-text"]`);
      if (elem !== null) {
        elem.addEventListener("click", function() {
          callbackFunction();
        });
      }
    }
  );
}, "pushFun");
var setClickEvent = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(ids, functionName, functionArgs) {
  ids.split(",").forEach(function(id) {
    setClickFun(id, functionName, functionArgs);
  });
  setClass(ids, "clickable");
}, "setClickEvent");
var bindFunctions = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(element) {
  funs.forEach(function(fun) {
    fun(element);
  });
}, "bindFunctions");
var ganttDb_default = {
  getConfig: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(() => (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .getConfig2 */ .nV)().gantt, "getConfig"),
  clear: clear2,
  setDateFormat,
  getDateFormat,
  enableInclusiveEndDates,
  endDatesAreInclusive,
  enableTopAxis,
  topAxisEnabled,
  setAxisFormat,
  getAxisFormat,
  setTickInterval,
  getTickInterval,
  setTodayMarker,
  getTodayMarker,
  setAccTitle: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .setAccTitle */ .GN,
  getAccTitle: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .getAccTitle */ .eu,
  setDiagramTitle: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .setDiagramTitle */ .g2,
  getDiagramTitle: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .getDiagramTitle */ .Kr,
  setDisplayMode,
  getDisplayMode,
  setAccDescription: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .setAccDescription */ .U$,
  getAccDescription: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .getAccDescription */ .Mx,
  addSection,
  getSections,
  getTasks,
  addTask,
  findTaskById,
  addTaskOrg,
  setIncludes,
  getIncludes,
  setExcludes,
  getExcludes,
  setClickEvent,
  setLink,
  getLinks,
  bindFunctions,
  parseDuration,
  isInvalidDate,
  setWeekday,
  getWeekday,
  setWeekend
};
function getTaskTags(data, task, tags2) {
  let matchFound = true;
  while (matchFound) {
    matchFound = false;
    tags2.forEach(function(t) {
      const pattern = "^\\s*" + t + "\\s*$";
      const regex = new RegExp(pattern);
      if (data[0].match(regex)) {
        task[t] = true;
        data.shift(1);
        matchFound = true;
      }
    });
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(getTaskTags, "getTaskTags");

// src/diagrams/gantt/ganttRenderer.js


var setConf = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function() {
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .log */ .cM.debug("Something is calling, setConf, remove the call");
}, "setConf");
var mapWeekdayToTimeFunction = {
  monday: d3__WEBPACK_IMPORTED_MODULE_8__/* .timeMonday */ .Ox9,
  tuesday: d3__WEBPACK_IMPORTED_MODULE_8__/* .timeTuesday */ .YDX,
  wednesday: d3__WEBPACK_IMPORTED_MODULE_8__/* .timeWednesday */ .EFj,
  thursday: d3__WEBPACK_IMPORTED_MODULE_8__/* .timeThursday */ .Igq,
  friday: d3__WEBPACK_IMPORTED_MODULE_8__/* .timeFriday */ .y2j,
  saturday: d3__WEBPACK_IMPORTED_MODULE_8__/* .timeSaturday */ .LqH,
  sunday: d3__WEBPACK_IMPORTED_MODULE_8__/* .timeSunday */ .Zyz
};
var getMaxIntersections = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)((tasks2, orderOffset) => {
  let timeline = [...tasks2].map(() => -Infinity);
  let sorted = [...tasks2].sort((a, b) => a.startTime - b.startTime || a.order - b.order);
  let maxIntersections = 0;
  for (const element of sorted) {
    for (let j = 0; j < timeline.length; j++) {
      if (element.startTime >= timeline[j]) {
        timeline[j] = element.endTime;
        element.order = j + orderOffset;
        if (j > maxIntersections) {
          maxIntersections = j;
        }
        break;
      }
    }
  }
  return maxIntersections;
}, "getMaxIntersections");
var w;
var draw = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(function(text, id, version, diagObj) {
  const conf = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .getConfig2 */ .nV)().gantt;
  const securityLevel = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .getConfig2 */ .nV)().securityLevel;
  let sandboxElement;
  if (securityLevel === "sandbox") {
    sandboxElement = (0,d3__WEBPACK_IMPORTED_MODULE_8__/* .select */ .Ys)("#i" + id);
  }
  const root = securityLevel === "sandbox" ? (0,d3__WEBPACK_IMPORTED_MODULE_8__/* .select */ .Ys)(sandboxElement.nodes()[0].contentDocument.body) : (0,d3__WEBPACK_IMPORTED_MODULE_8__/* .select */ .Ys)("body");
  const doc = securityLevel === "sandbox" ? sandboxElement.nodes()[0].contentDocument : document;
  const elem = doc.getElementById(id);
  w = elem.parentElement.offsetWidth;
  if (w === void 0) {
    w = 1200;
  }
  if (conf.useWidth !== void 0) {
    w = conf.useWidth;
  }
  const taskArray = diagObj.db.getTasks();
  let categories = [];
  for (const element of taskArray) {
    categories.push(element.type);
  }
  categories = checkUnique(categories);
  const categoryHeights = {};
  let h = 2 * conf.topPadding;
  if (diagObj.db.getDisplayMode() === "compact" || conf.displayMode === "compact") {
    const categoryElements = {};
    for (const element of taskArray) {
      if (categoryElements[element.section] === void 0) {
        categoryElements[element.section] = [element];
      } else {
        categoryElements[element.section].push(element);
      }
    }
    let intersections = 0;
    for (const category of Object.keys(categoryElements)) {
      const categoryHeight = getMaxIntersections(categoryElements[category], intersections) + 1;
      intersections += categoryHeight;
      h += categoryHeight * (conf.barHeight + conf.barGap);
      categoryHeights[category] = categoryHeight;
    }
  } else {
    h += taskArray.length * (conf.barHeight + conf.barGap);
    for (const category of categories) {
      categoryHeights[category] = taskArray.filter((task) => task.type === category).length;
    }
  }
  elem.setAttribute("viewBox", "0 0 " + w + " " + h);
  const svg = root.select(`[id="${id}"]`);
  const timeScale = (0,d3__WEBPACK_IMPORTED_MODULE_8__/* .scaleTime */ .Xf)().domain([
    (0,d3__WEBPACK_IMPORTED_MODULE_8__/* .min */ .VV$)(taskArray, function(d) {
      return d.startTime;
    }),
    (0,d3__WEBPACK_IMPORTED_MODULE_8__/* .max */ .Fp7)(taskArray, function(d) {
      return d.endTime;
    })
  ]).rangeRound([0, w - conf.leftPadding - conf.rightPadding]);
  function taskCompare(a, b) {
    const taskA = a.startTime;
    const taskB = b.startTime;
    let result = 0;
    if (taskA > taskB) {
      result = 1;
    } else if (taskA < taskB) {
      result = -1;
    }
    return result;
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(taskCompare, "taskCompare");
  taskArray.sort(taskCompare);
  makeGantt(taskArray, w, h);
  (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .configureSvgSize */ .v2)(svg, h, w, conf.useMaxWidth);
  svg.append("text").text(diagObj.db.getDiagramTitle()).attr("x", w / 2).attr("y", conf.titleTopMargin).attr("class", "titleText");
  function makeGantt(tasks2, pageWidth, pageHeight) {
    const barHeight = conf.barHeight;
    const gap = barHeight + conf.barGap;
    const topPadding = conf.topPadding;
    const leftPadding = conf.leftPadding;
    const colorScale = (0,d3__WEBPACK_IMPORTED_MODULE_8__/* .scaleLinear */ .BYU)().domain([0, categories.length]).range(["#00B9FA", "#F95002"]).interpolate(d3__WEBPACK_IMPORTED_MODULE_8__/* .interpolateHcl */ .JHv);
    drawExcludeDays(
      gap,
      topPadding,
      leftPadding,
      pageWidth,
      pageHeight,
      tasks2,
      diagObj.db.getExcludes(),
      diagObj.db.getIncludes()
    );
    makeGrid(leftPadding, topPadding, pageWidth, pageHeight);
    drawRects(tasks2, gap, topPadding, leftPadding, barHeight, colorScale, pageWidth, pageHeight);
    vertLabels(gap, topPadding, leftPadding, barHeight, colorScale);
    drawToday(leftPadding, topPadding, pageWidth, pageHeight);
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(makeGantt, "makeGantt");
  function drawRects(theArray, theGap, theTopPad, theSidePad, theBarHeight, theColorScale, w2) {
    theArray.sort((a, b) => a.vert === b.vert ? 0 : a.vert ? 1 : -1);
    const uniqueTaskOrderIds = [...new Set(theArray.map((item) => item.order))];
    const uniqueTasks = uniqueTaskOrderIds.map((id2) => theArray.find((item) => item.order === id2));
    svg.append("g").selectAll("rect").data(uniqueTasks).enter().append("rect").attr("x", 0).attr("y", function(d, i) {
      i = d.order;
      return i * theGap + theTopPad - 2;
    }).attr("width", function() {
      return w2 - conf.rightPadding / 2;
    }).attr("height", theGap).attr("class", function(d) {
      for (const [i, category] of categories.entries()) {
        if (d.type === category) {
          return "section section" + i % conf.numberSectionStyles;
        }
      }
      return "section section0";
    }).enter();
    const rectangles = svg.append("g").selectAll("rect").data(theArray).enter();
    const links2 = diagObj.db.getLinks();
    rectangles.append("rect").attr("id", function(d) {
      return d.id;
    }).attr("rx", 3).attr("ry", 3).attr("x", function(d) {
      if (d.milestone) {
        return timeScale(d.startTime) + theSidePad + 0.5 * (timeScale(d.endTime) - timeScale(d.startTime)) - 0.5 * theBarHeight;
      }
      return timeScale(d.startTime) + theSidePad;
    }).attr("y", function(d, i) {
      i = d.order;
      if (d.vert) {
        return conf.gridLineStartPadding;
      }
      return i * theGap + theTopPad;
    }).attr("width", function(d) {
      if (d.milestone) {
        return theBarHeight;
      }
      if (d.vert) {
        return 0.08 * theBarHeight;
      }
      return timeScale(d.renderEndTime || d.endTime) - timeScale(d.startTime);
    }).attr("height", function(d) {
      if (d.vert) {
        return taskArray.length * (conf.barHeight + conf.barGap) + conf.barHeight * 2;
      }
      return theBarHeight;
    }).attr("transform-origin", function(d, i) {
      i = d.order;
      return (timeScale(d.startTime) + theSidePad + 0.5 * (timeScale(d.endTime) - timeScale(d.startTime))).toString() + "px " + (i * theGap + theTopPad + 0.5 * theBarHeight).toString() + "px";
    }).attr("class", function(d) {
      const res = "task";
      let classStr = "";
      if (d.classes.length > 0) {
        classStr = d.classes.join(" ");
      }
      let secNum = 0;
      for (const [i, category] of categories.entries()) {
        if (d.type === category) {
          secNum = i % conf.numberSectionStyles;
        }
      }
      let taskClass = "";
      if (d.active) {
        if (d.crit) {
          taskClass += " activeCrit";
        } else {
          taskClass = " active";
        }
      } else if (d.done) {
        if (d.crit) {
          taskClass = " doneCrit";
        } else {
          taskClass = " done";
        }
      } else {
        if (d.crit) {
          taskClass += " crit";
        }
      }
      if (taskClass.length === 0) {
        taskClass = " task";
      }
      if (d.milestone) {
        taskClass = " milestone " + taskClass;
      }
      if (d.vert) {
        taskClass = " vert " + taskClass;
      }
      taskClass += secNum;
      taskClass += " " + classStr;
      return res + taskClass;
    });
    rectangles.append("text").attr("id", function(d) {
      return d.id + "-text";
    }).text(function(d) {
      return d.task;
    }).attr("font-size", conf.fontSize).attr("x", function(d) {
      let startX = timeScale(d.startTime);
      let endX = timeScale(d.renderEndTime || d.endTime);
      if (d.milestone) {
        startX += 0.5 * (timeScale(d.endTime) - timeScale(d.startTime)) - 0.5 * theBarHeight;
        endX = startX + theBarHeight;
      }
      if (d.vert) {
        return timeScale(d.startTime) + theSidePad;
      }
      const textWidth = this.getBBox().width;
      if (textWidth > endX - startX) {
        if (endX + textWidth + 1.5 * conf.leftPadding > w2) {
          return startX + theSidePad - 5;
        } else {
          return endX + theSidePad + 5;
        }
      } else {
        return (endX - startX) / 2 + startX + theSidePad;
      }
    }).attr("y", function(d, i) {
      if (d.vert) {
        return conf.gridLineStartPadding + taskArray.length * (conf.barHeight + conf.barGap) + 60;
      }
      i = d.order;
      return i * theGap + conf.barHeight / 2 + (conf.fontSize / 2 - 2) + theTopPad;
    }).attr("text-height", theBarHeight).attr("class", function(d) {
      const startX = timeScale(d.startTime);
      let endX = timeScale(d.endTime);
      if (d.milestone) {
        endX = startX + theBarHeight;
      }
      const textWidth = this.getBBox().width;
      let classStr = "";
      if (d.classes.length > 0) {
        classStr = d.classes.join(" ");
      }
      let secNum = 0;
      for (const [i, category] of categories.entries()) {
        if (d.type === category) {
          secNum = i % conf.numberSectionStyles;
        }
      }
      let taskType = "";
      if (d.active) {
        if (d.crit) {
          taskType = "activeCritText" + secNum;
        } else {
          taskType = "activeText" + secNum;
        }
      }
      if (d.done) {
        if (d.crit) {
          taskType = taskType + " doneCritText" + secNum;
        } else {
          taskType = taskType + " doneText" + secNum;
        }
      } else {
        if (d.crit) {
          taskType = taskType + " critText" + secNum;
        }
      }
      if (d.milestone) {
        taskType += " milestoneText";
      }
      if (d.vert) {
        taskType += " vertText";
      }
      if (textWidth > endX - startX) {
        if (endX + textWidth + 1.5 * conf.leftPadding > w2) {
          return classStr + " taskTextOutsideLeft taskTextOutside" + secNum + " " + taskType;
        } else {
          return classStr + " taskTextOutsideRight taskTextOutside" + secNum + " " + taskType + " width-" + textWidth;
        }
      } else {
        return classStr + " taskText taskText" + secNum + " " + taskType + " width-" + textWidth;
      }
    });
    const securityLevel2 = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .getConfig2 */ .nV)().securityLevel;
    if (securityLevel2 === "sandbox") {
      let sandboxElement2;
      sandboxElement2 = (0,d3__WEBPACK_IMPORTED_MODULE_8__/* .select */ .Ys)("#i" + id);
      const doc2 = sandboxElement2.nodes()[0].contentDocument;
      rectangles.filter(function(d) {
        return links2.has(d.id);
      }).each(function(o) {
        var taskRect = doc2.querySelector("#" + o.id);
        var taskText = doc2.querySelector("#" + o.id + "-text");
        const oldParent = taskRect.parentNode;
        var Link = doc2.createElement("a");
        Link.setAttribute("xlink:href", links2.get(o.id));
        Link.setAttribute("target", "_top");
        oldParent.appendChild(Link);
        Link.appendChild(taskRect);
        Link.appendChild(taskText);
      });
    }
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(drawRects, "drawRects");
  function drawExcludeDays(theGap, theTopPad, theSidePad, w2, h2, tasks2, excludes2, includes2) {
    if (excludes2.length === 0 && includes2.length === 0) {
      return;
    }
    let minTime;
    let maxTime;
    for (const { startTime, endTime } of tasks2) {
      if (minTime === void 0 || startTime < minTime) {
        minTime = startTime;
      }
      if (maxTime === void 0 || endTime > maxTime) {
        maxTime = endTime;
      }
    }
    if (!minTime || !maxTime) {
      return;
    }
    if (dayjs__WEBPACK_IMPORTED_MODULE_4___default()(maxTime).diff(dayjs__WEBPACK_IMPORTED_MODULE_4___default()(minTime), "year") > 5) {
      _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .log */ .cM.warn(
        "The difference between the min and max time is more than 5 years. This will cause performance issues. Skipping drawing exclude days."
      );
      return;
    }
    const dateFormat2 = diagObj.db.getDateFormat();
    const excludeRanges = [];
    let range = null;
    let d = dayjs__WEBPACK_IMPORTED_MODULE_4___default()(minTime);
    while (d.valueOf() <= maxTime) {
      if (diagObj.db.isInvalidDate(d, dateFormat2, excludes2, includes2)) {
        if (!range) {
          range = {
            start: d,
            end: d
          };
        } else {
          range.end = d;
        }
      } else {
        if (range) {
          excludeRanges.push(range);
          range = null;
        }
      }
      d = d.add(1, "d");
    }
    const rectangles = svg.append("g").selectAll("rect").data(excludeRanges).enter();
    rectangles.append("rect").attr("id", (d2) => "exclude-" + d2.start.format("YYYY-MM-DD")).attr("x", (d2) => timeScale(d2.start.startOf("day")) + theSidePad).attr("y", conf.gridLineStartPadding).attr("width", (d2) => timeScale(d2.end.endOf("day")) - timeScale(d2.start.startOf("day"))).attr("height", h2 - theTopPad - conf.gridLineStartPadding).attr("transform-origin", function(d2, i) {
      return (timeScale(d2.start) + theSidePad + 0.5 * (timeScale(d2.end) - timeScale(d2.start))).toString() + "px " + (i * theGap + 0.5 * h2).toString() + "px";
    }).attr("class", "exclude-range");
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(drawExcludeDays, "drawExcludeDays");
  function makeGrid(theSidePad, theTopPad, w2, h2) {
    const dateFormat2 = diagObj.db.getDateFormat();
    const userAxisFormat = diagObj.db.getAxisFormat();
    let axisFormat2;
    if (userAxisFormat) {
      axisFormat2 = userAxisFormat;
    } else if (dateFormat2 === "D") {
      axisFormat2 = "%d";
    } else {
      axisFormat2 = conf.axisFormat ?? "%Y-%m-%d";
    }
    let bottomXAxis = (0,d3__WEBPACK_IMPORTED_MODULE_8__/* .axisBottom */ .LLu)(timeScale).tickSize(-h2 + theTopPad + conf.gridLineStartPadding).tickFormat((0,d3__WEBPACK_IMPORTED_MODULE_8__/* .timeFormat */ .i$Z)(axisFormat2));
    const reTickInterval = /^([1-9]\d*)(millisecond|second|minute|hour|day|week|month)$/;
    const resultTickInterval = reTickInterval.exec(
      diagObj.db.getTickInterval() || conf.tickInterval
    );
    if (resultTickInterval !== null) {
      const every = resultTickInterval[1];
      const interval = resultTickInterval[2];
      const weekday2 = diagObj.db.getWeekday() || conf.weekday;
      switch (interval) {
        case "millisecond":
          bottomXAxis.ticks(d3__WEBPACK_IMPORTED_MODULE_8__/* .timeMillisecond */ .U8T.every(every));
          break;
        case "second":
          bottomXAxis.ticks(d3__WEBPACK_IMPORTED_MODULE_8__/* .timeSecond */ .S1K.every(every));
          break;
        case "minute":
          bottomXAxis.ticks(d3__WEBPACK_IMPORTED_MODULE_8__/* .timeMinute */ .Z_i.every(every));
          break;
        case "hour":
          bottomXAxis.ticks(d3__WEBPACK_IMPORTED_MODULE_8__/* .timeHour */ .WQD.every(every));
          break;
        case "day":
          bottomXAxis.ticks(d3__WEBPACK_IMPORTED_MODULE_8__/* .timeDay */ .rr1.every(every));
          break;
        case "week":
          bottomXAxis.ticks(mapWeekdayToTimeFunction[weekday2].every(every));
          break;
        case "month":
          bottomXAxis.ticks(d3__WEBPACK_IMPORTED_MODULE_8__/* .timeMonth */ .F0B.every(every));
          break;
      }
    }
    svg.append("g").attr("class", "grid").attr("transform", "translate(" + theSidePad + ", " + (h2 - 50) + ")").call(bottomXAxis).selectAll("text").style("text-anchor", "middle").attr("fill", "#000").attr("stroke", "none").attr("font-size", 10).attr("dy", "1em");
    if (diagObj.db.topAxisEnabled() || conf.topAxis) {
      let topXAxis = (0,d3__WEBPACK_IMPORTED_MODULE_8__/* .axisTop */ .F5q)(timeScale).tickSize(-h2 + theTopPad + conf.gridLineStartPadding).tickFormat((0,d3__WEBPACK_IMPORTED_MODULE_8__/* .timeFormat */ .i$Z)(axisFormat2));
      if (resultTickInterval !== null) {
        const every = resultTickInterval[1];
        const interval = resultTickInterval[2];
        const weekday2 = diagObj.db.getWeekday() || conf.weekday;
        switch (interval) {
          case "millisecond":
            topXAxis.ticks(d3__WEBPACK_IMPORTED_MODULE_8__/* .timeMillisecond */ .U8T.every(every));
            break;
          case "second":
            topXAxis.ticks(d3__WEBPACK_IMPORTED_MODULE_8__/* .timeSecond */ .S1K.every(every));
            break;
          case "minute":
            topXAxis.ticks(d3__WEBPACK_IMPORTED_MODULE_8__/* .timeMinute */ .Z_i.every(every));
            break;
          case "hour":
            topXAxis.ticks(d3__WEBPACK_IMPORTED_MODULE_8__/* .timeHour */ .WQD.every(every));
            break;
          case "day":
            topXAxis.ticks(d3__WEBPACK_IMPORTED_MODULE_8__/* .timeDay */ .rr1.every(every));
            break;
          case "week":
            topXAxis.ticks(mapWeekdayToTimeFunction[weekday2].every(every));
            break;
          case "month":
            topXAxis.ticks(d3__WEBPACK_IMPORTED_MODULE_8__/* .timeMonth */ .F0B.every(every));
            break;
        }
      }
      svg.append("g").attr("class", "grid").attr("transform", "translate(" + theSidePad + ", " + theTopPad + ")").call(topXAxis).selectAll("text").style("text-anchor", "middle").attr("fill", "#000").attr("stroke", "none").attr("font-size", 10);
    }
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(makeGrid, "makeGrid");
  function vertLabels(theGap, theTopPad) {
    let prevGap = 0;
    const numOccurrences = Object.keys(categoryHeights).map((d) => [d, categoryHeights[d]]);
    svg.append("g").selectAll("text").data(numOccurrences).enter().append(function(d) {
      const rows = d[0].split(_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_1__/* .common_default */ .SY.lineBreakRegex);
      const dy = -(rows.length - 1) / 2;
      const svgLabel = doc.createElementNS("http://www.w3.org/2000/svg", "text");
      svgLabel.setAttribute("dy", dy + "em");
      for (const [j, row] of rows.entries()) {
        const tspan = doc.createElementNS("http://www.w3.org/2000/svg", "tspan");
        tspan.setAttribute("alignment-baseline", "central");
        tspan.setAttribute("x", "10");
        if (j > 0) {
          tspan.setAttribute("dy", "1em");
        }
        tspan.textContent = row;
        svgLabel.appendChild(tspan);
      }
      return svgLabel;
    }).attr("x", 10).attr("y", function(d, i) {
      if (i > 0) {
        for (let j = 0; j < i; j++) {
          prevGap += numOccurrences[i - 1][1];
          return d[1] * theGap / 2 + prevGap * theGap + theTopPad;
        }
      } else {
        return d[1] * theGap / 2 + theTopPad;
      }
    }).attr("font-size", conf.sectionFontSize).attr("class", function(d) {
      for (const [i, category] of categories.entries()) {
        if (d[0] === category) {
          return "sectionTitle sectionTitle" + i % conf.numberSectionStyles;
        }
      }
      return "sectionTitle";
    });
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(vertLabels, "vertLabels");
  function drawToday(theSidePad, theTopPad, w2, h2) {
    const todayMarker2 = diagObj.db.getTodayMarker();
    if (todayMarker2 === "off") {
      return;
    }
    const todayG = svg.append("g").attr("class", "today");
    const today = /* @__PURE__ */ new Date();
    const todayLine = todayG.append("line");
    todayLine.attr("x1", timeScale(today) + theSidePad).attr("x2", timeScale(today) + theSidePad).attr("y1", conf.titleTopMargin).attr("y2", h2 - conf.titleTopMargin).attr("class", "today");
    if (todayMarker2 !== "") {
      todayLine.attr("style", todayMarker2.replace(/,/g, ";"));
    }
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(drawToday, "drawToday");
  function checkUnique(arr) {
    const hash = {};
    const result = [];
    for (let i = 0, l = arr.length; i < l; ++i) {
      if (!Object.prototype.hasOwnProperty.call(hash, arr[i])) {
        hash[arr[i]] = true;
        result.push(arr[i]);
      }
    }
    return result;
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)(checkUnique, "checkUnique");
}, "draw");
var ganttRenderer_default = {
  setConf,
  draw
};

// src/diagrams/gantt/styles.js
var getStyles = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .__name */ .eW)((options) => `
  .mermaid-main-font {
        font-family: ${options.fontFamily};
  }

  .exclude-range {
    fill: ${options.excludeBkgColor};
  }

  .section {
    stroke: none;
    opacity: 0.2;
  }

  .section0 {
    fill: ${options.sectionBkgColor};
  }

  .section2 {
    fill: ${options.sectionBkgColor2};
  }

  .section1,
  .section3 {
    fill: ${options.altSectionBkgColor};
    opacity: 0.2;
  }

  .sectionTitle0 {
    fill: ${options.titleColor};
  }

  .sectionTitle1 {
    fill: ${options.titleColor};
  }

  .sectionTitle2 {
    fill: ${options.titleColor};
  }

  .sectionTitle3 {
    fill: ${options.titleColor};
  }

  .sectionTitle {
    text-anchor: start;
    font-family: ${options.fontFamily};
  }


  /* Grid and axis */

  .grid .tick {
    stroke: ${options.gridColor};
    opacity: 0.8;
    shape-rendering: crispEdges;
  }

  .grid .tick text {
    font-family: ${options.fontFamily};
    fill: ${options.textColor};
  }

  .grid path {
    stroke-width: 0;
  }


  /* Today line */

  .today {
    fill: none;
    stroke: ${options.todayLineColor};
    stroke-width: 2px;
  }


  /* Task styling */

  /* Default task */

  .task {
    stroke-width: 2;
  }

  .taskText {
    text-anchor: middle;
    font-family: ${options.fontFamily};
  }

  .taskTextOutsideRight {
    fill: ${options.taskTextDarkColor};
    text-anchor: start;
    font-family: ${options.fontFamily};
  }

  .taskTextOutsideLeft {
    fill: ${options.taskTextDarkColor};
    text-anchor: end;
  }


  /* Special case clickable */

  .task.clickable {
    cursor: pointer;
  }

  .taskText.clickable {
    cursor: pointer;
    fill: ${options.taskTextClickableColor} !important;
    font-weight: bold;
  }

  .taskTextOutsideLeft.clickable {
    cursor: pointer;
    fill: ${options.taskTextClickableColor} !important;
    font-weight: bold;
  }

  .taskTextOutsideRight.clickable {
    cursor: pointer;
    fill: ${options.taskTextClickableColor} !important;
    font-weight: bold;
  }


  /* Specific task settings for the sections*/

  .taskText0,
  .taskText1,
  .taskText2,
  .taskText3 {
    fill: ${options.taskTextColor};
  }

  .task0,
  .task1,
  .task2,
  .task3 {
    fill: ${options.taskBkgColor};
    stroke: ${options.taskBorderColor};
  }

  .taskTextOutside0,
  .taskTextOutside2
  {
    fill: ${options.taskTextOutsideColor};
  }

  .taskTextOutside1,
  .taskTextOutside3 {
    fill: ${options.taskTextOutsideColor};
  }


  /* Active task */

  .active0,
  .active1,
  .active2,
  .active3 {
    fill: ${options.activeTaskBkgColor};
    stroke: ${options.activeTaskBorderColor};
  }

  .activeText0,
  .activeText1,
  .activeText2,
  .activeText3 {
    fill: ${options.taskTextDarkColor} !important;
  }


  /* Completed task */

  .done0,
  .done1,
  .done2,
  .done3 {
    stroke: ${options.doneTaskBorderColor};
    fill: ${options.doneTaskBkgColor};
    stroke-width: 2;
  }

  .doneText0,
  .doneText1,
  .doneText2,
  .doneText3 {
    fill: ${options.taskTextDarkColor} !important;
  }


  /* Tasks on the critical line */

  .crit0,
  .crit1,
  .crit2,
  .crit3 {
    stroke: ${options.critBorderColor};
    fill: ${options.critBkgColor};
    stroke-width: 2;
  }

  .activeCrit0,
  .activeCrit1,
  .activeCrit2,
  .activeCrit3 {
    stroke: ${options.critBorderColor};
    fill: ${options.activeTaskBkgColor};
    stroke-width: 2;
  }

  .doneCrit0,
  .doneCrit1,
  .doneCrit2,
  .doneCrit3 {
    stroke: ${options.critBorderColor};
    fill: ${options.doneTaskBkgColor};
    stroke-width: 2;
    cursor: pointer;
    shape-rendering: crispEdges;
  }

  .milestone {
    transform: rotate(45deg) scale(0.8,0.8);
  }

  .milestoneText {
    font-style: italic;
  }
  .doneCritText0,
  .doneCritText1,
  .doneCritText2,
  .doneCritText3 {
    fill: ${options.taskTextDarkColor} !important;
  }

  .vert {
    stroke: ${options.vertLineColor};
  }

  .vertText {
    font-size: 15px;
    text-anchor: middle;
    fill: ${options.vertLineColor} !important;
  }

  .activeCritText0,
  .activeCritText1,
  .activeCritText2,
  .activeCritText3 {
    fill: ${options.taskTextDarkColor} !important;
  }

  .titleText {
    text-anchor: middle;
    font-size: 18px;
    fill: ${options.titleColor || options.textColor};
    font-family: ${options.fontFamily};
  }
`, "getStyles");
var styles_default = getStyles;

// src/diagrams/gantt/ganttDiagram.ts
var diagram = {
  parser: gantt_default,
  db: ganttDb_default,
  renderer: ganttRenderer_default,
  styles: styles_default
};



/***/ }),

/***/ 85747:
/***/ (function(module) {

!function(e,t){ true?module.exports=t():0}(this,(function(){"use strict";return function(e,t){var r=t.prototype,n=r.format;r.format=function(e){var t=this,r=this.$locale();if(!this.isValid())return n.bind(this)(e);var s=this.$utils(),a=(e||"YYYY-MM-DDTHH:mm:ssZ").replace(/\[([^\]]+)]|Q|wo|ww|w|WW|W|zzz|z|gggg|GGGG|Do|X|x|k{1,2}|S/g,(function(e){switch(e){case"Q":return Math.ceil((t.$M+1)/3);case"Do":return r.ordinal(t.$D);case"gggg":return t.weekYear();case"GGGG":return t.isoWeekYear();case"wo":return r.ordinal(t.week(),"W");case"w":case"ww":return s.s(t.week(),"w"===e?1:2,"0");case"W":case"WW":return s.s(t.isoWeek(),"W"===e?1:2,"0");case"k":case"kk":return s.s(String(0===t.$H?24:t.$H),"k"===e?1:2,"0");case"X":return Math.floor(t.$d.getTime()/1e3);case"x":return t.$d.getTime();case"z":return"["+t.offsetName()+"]";case"zzz":return"["+t.offsetName("long")+"]";default:return e}}));return n.bind(this)(a)}}}));

/***/ }),

/***/ 14881:
/***/ (function(module) {

!function(e,t){ true?module.exports=t():0}(this,(function(){"use strict";var e={LTS:"h:mm:ss A",LT:"h:mm A",L:"MM/DD/YYYY",LL:"MMMM D, YYYY",LLL:"MMMM D, YYYY h:mm A",LLLL:"dddd, MMMM D, YYYY h:mm A"},t=/(\[[^[]*\])|([-_:/.,()\s]+)|(A|a|Q|YYYY|YY?|ww?|MM?M?M?|Do|DD?|hh?|HH?|mm?|ss?|S{1,3}|z|ZZ?)/g,n=/\d/,r=/\d\d/,i=/\d\d?/,o=/\d*[^-_:/,()\s\d]+/,s={},a=function(e){return(e=+e)+(e>68?1900:2e3)};var f=function(e){return function(t){this[e]=+t}},h=[/[+-]\d\d:?(\d\d)?|Z/,function(e){(this.zone||(this.zone={})).offset=function(e){if(!e)return 0;if("Z"===e)return 0;var t=e.match(/([+-]|\d\d)/g),n=60*t[1]+(+t[2]||0);return 0===n?0:"+"===t[0]?-n:n}(e)}],u=function(e){var t=s[e];return t&&(t.indexOf?t:t.s.concat(t.f))},d=function(e,t){var n,r=s.meridiem;if(r){for(var i=1;i<=24;i+=1)if(e.indexOf(r(i,0,t))>-1){n=i>12;break}}else n=e===(t?"pm":"PM");return n},c={A:[o,function(e){this.afternoon=d(e,!1)}],a:[o,function(e){this.afternoon=d(e,!0)}],Q:[n,function(e){this.month=3*(e-1)+1}],S:[n,function(e){this.milliseconds=100*+e}],SS:[r,function(e){this.milliseconds=10*+e}],SSS:[/\d{3}/,function(e){this.milliseconds=+e}],s:[i,f("seconds")],ss:[i,f("seconds")],m:[i,f("minutes")],mm:[i,f("minutes")],H:[i,f("hours")],h:[i,f("hours")],HH:[i,f("hours")],hh:[i,f("hours")],D:[i,f("day")],DD:[r,f("day")],Do:[o,function(e){var t=s.ordinal,n=e.match(/\d+/);if(this.day=n[0],t)for(var r=1;r<=31;r+=1)t(r).replace(/\[|\]/g,"")===e&&(this.day=r)}],w:[i,f("week")],ww:[r,f("week")],M:[i,f("month")],MM:[r,f("month")],MMM:[o,function(e){var t=u("months"),n=(u("monthsShort")||t.map((function(e){return e.slice(0,3)}))).indexOf(e)+1;if(n<1)throw new Error;this.month=n%12||n}],MMMM:[o,function(e){var t=u("months").indexOf(e)+1;if(t<1)throw new Error;this.month=t%12||t}],Y:[/[+-]?\d+/,f("year")],YY:[r,function(e){this.year=a(e)}],YYYY:[/\d{4}/,f("year")],Z:h,ZZ:h};function l(n){var r,i;r=n,i=s&&s.formats;for(var o=(n=r.replace(/(\[[^\]]+])|(LTS?|l{1,4}|L{1,4})/g,(function(t,n,r){var o=r&&r.toUpperCase();return n||i[r]||e[r]||i[o].replace(/(\[[^\]]+])|(MMMM|MM|DD|dddd)/g,(function(e,t,n){return t||n.slice(1)}))}))).match(t),a=o.length,f=0;f<a;f+=1){var h=o[f],u=c[h],d=u&&u[0],l=u&&u[1];o[f]=l?{regex:d,parser:l}:h.replace(/^\[|\]$/g,"")}return function(e){for(var t={},n=0,r=0;n<a;n+=1){var i=o[n];if("string"==typeof i)r+=i.length;else{var s=i.regex,f=i.parser,h=e.slice(r),u=s.exec(h)[0];f.call(t,u),e=e.replace(u,"")}}return function(e){var t=e.afternoon;if(void 0!==t){var n=e.hours;t?n<12&&(e.hours+=12):12===n&&(e.hours=0),delete e.afternoon}}(t),t}}return function(e,t,n){n.p.customParseFormat=!0,e&&e.parseTwoDigitYear&&(a=e.parseTwoDigitYear);var r=t.prototype,i=r.parse;r.parse=function(e){var t=e.date,r=e.utc,o=e.args;this.$u=r;var a=o[1];if("string"==typeof a){var f=!0===o[2],h=!0===o[3],u=f||h,d=o[2];h&&(d=o[2]),s=this.$locale(),!f&&d&&(s=n.Ls[d]),this.$d=function(e,t,n,r){try{if(["x","X"].indexOf(t)>-1)return new Date(("X"===t?1e3:1)*e);var i=l(t)(e),o=i.year,s=i.month,a=i.day,f=i.hours,h=i.minutes,u=i.seconds,d=i.milliseconds,c=i.zone,m=i.week,M=new Date,Y=a||(o||s?1:M.getDate()),p=o||M.getFullYear(),v=0;o&&!s||(v=s>0?s-1:M.getMonth());var D,w=f||0,g=h||0,y=u||0,L=d||0;return c?new Date(Date.UTC(p,v,Y,w,g,y,L+60*c.offset*1e3)):n?new Date(Date.UTC(p,v,Y,w,g,y,L)):(D=new Date(p,v,Y,w,g,y,L),m&&(D=r(D).week(m).toDate()),D)}catch(e){return new Date("")}}(t,a,r,n),this.init(),d&&!0!==d&&(this.$L=this.locale(d).$L),u&&t!=this.format(a)&&(this.$d=new Date("")),s={}}else if(a instanceof Array)for(var c=a.length,m=1;m<=c;m+=1){o[1]=a[m-1];var M=n.apply(this,o);if(M.isValid()){this.$d=M.$d,this.$L=M.$L,this.init();break}m===c&&(this.$d=new Date(""))}else i.call(this,e)}}}));

/***/ }),

/***/ 62539:
/***/ (function(module) {

!function(e,t){ true?module.exports=t():0}(this,(function(){"use strict";var e="day";return function(t,i,s){var a=function(t){return t.add(4-t.isoWeekday(),e)},d=i.prototype;d.isoWeekYear=function(){return a(this).year()},d.isoWeek=function(t){if(!this.$utils().u(t))return this.add(7*(t-this.isoWeek()),e);var i,d,n,o,r=a(this),u=(i=this.isoWeekYear(),d=this.$u,n=(d?s.utc:s)().year(i).startOf("year"),o=4-n.isoWeekday(),n.isoWeekday()>4&&(o+=7),n.add(o,e));return r.diff(u,"week")+1},d.isoWeekday=function(e){return this.$utils().u(e)?this.day()||7:this.day(this.day()%7?e:e-7)};var n=d.startOf;d.startOf=function(e,t){var i=this.$utils(),s=!!i.u(t)||t;return"isoweek"===i.p(e)?s?this.date(this.date()-(this.isoWeekday()-1)).startOf("day"):this.date(this.date()-1-(this.isoWeekday()-1)+7).endOf("day"):n.bind(this)(e,t)}}}));

/***/ })

}]);
//# sourceMappingURL=1830.d2ff9069451b0d8dd023.js.map?v=d2ff9069451b0d8dd023
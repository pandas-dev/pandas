"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[369],{

/***/ 28862:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   G: () => (/* binding */ getIconStyles)
/* harmony export */ });
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(74999);


// src/diagrams/globalStyles.ts
var getIconStyles = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(() => `
  /* Font Awesome icon styling - consolidated */
  .label-icon {
    display: inline-block;
    height: 1em;
    overflow: visible;
    vertical-align: -0.125em;
  }
  
  .node .label-icon path {
    fill: currentColor;
    stroke: revert;
    stroke-width: revert;
  }
`, "getIconStyles");




/***/ }),

/***/ 30369:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _chunk_EXTU4WIE_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(98154);
/* harmony import */ var _chunk_FMBD7UC4_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(28862);
/* harmony import */ var _chunk_MI3HLSF2_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(71262);
/* harmony import */ var _chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(26212);
/* harmony import */ var _chunk_CVBHYZKI_mjs__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(22462);
/* harmony import */ var _chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(37756);
/* harmony import */ var _chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(2111);
/* harmony import */ var _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(17175);
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(74999);
/* harmony import */ var khroma__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(28186);
/* harmony import */ var khroma__WEBPACK_IMPORTED_MODULE_11__ = __webpack_require__(28482);
/* harmony import */ var khroma__WEBPACK_IMPORTED_MODULE_12__ = __webpack_require__(57838);











// src/diagrams/kanban/parser/kanban.jison
var parser = (function() {
  var o = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function(k, v, o2, l) {
    for (o2 = o2 || {}, l = k.length; l--; o2[k[l]] = v) ;
    return o2;
  }, "o"), $V0 = [1, 4], $V1 = [1, 13], $V2 = [1, 12], $V3 = [1, 15], $V4 = [1, 16], $V5 = [1, 20], $V6 = [1, 19], $V7 = [6, 7, 8], $V8 = [1, 26], $V9 = [1, 24], $Va = [1, 25], $Vb = [6, 7, 11], $Vc = [1, 31], $Vd = [6, 7, 11, 24], $Ve = [1, 6, 13, 16, 17, 20, 23], $Vf = [1, 35], $Vg = [1, 36], $Vh = [1, 6, 7, 11, 13, 16, 17, 20, 23], $Vi = [1, 38];
  var parser2 = {
    trace: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function trace() {
    }, "trace"),
    yy: {},
    symbols_: { "error": 2, "start": 3, "mindMap": 4, "spaceLines": 5, "SPACELINE": 6, "NL": 7, "KANBAN": 8, "document": 9, "stop": 10, "EOF": 11, "statement": 12, "SPACELIST": 13, "node": 14, "shapeData": 15, "ICON": 16, "CLASS": 17, "nodeWithId": 18, "nodeWithoutId": 19, "NODE_DSTART": 20, "NODE_DESCR": 21, "NODE_DEND": 22, "NODE_ID": 23, "SHAPE_DATA": 24, "$accept": 0, "$end": 1 },
    terminals_: { 2: "error", 6: "SPACELINE", 7: "NL", 8: "KANBAN", 11: "EOF", 13: "SPACELIST", 16: "ICON", 17: "CLASS", 20: "NODE_DSTART", 21: "NODE_DESCR", 22: "NODE_DEND", 23: "NODE_ID", 24: "SHAPE_DATA" },
    productions_: [0, [3, 1], [3, 2], [5, 1], [5, 2], [5, 2], [4, 2], [4, 3], [10, 1], [10, 1], [10, 1], [10, 2], [10, 2], [9, 3], [9, 2], [12, 3], [12, 2], [12, 2], [12, 2], [12, 1], [12, 2], [12, 1], [12, 1], [12, 1], [12, 1], [14, 1], [14, 1], [19, 3], [18, 1], [18, 4], [15, 2], [15, 1]],
    performAction: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function anonymous(yytext, yyleng, yylineno, yy, yystate, $$, _$) {
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
          yy.getLogger().info("Node: ", $$[$0 - 1].id);
          yy.addNode($$[$0 - 2].length, $$[$0 - 1].id, $$[$0 - 1].descr, $$[$0 - 1].type, $$[$0]);
          break;
        case 16:
          yy.getLogger().info("Node: ", $$[$0].id);
          yy.addNode($$[$0 - 1].length, $$[$0].id, $$[$0].descr, $$[$0].type);
          break;
        case 17:
          yy.getLogger().trace("Icon: ", $$[$0]);
          yy.decorateNode({ icon: $$[$0] });
          break;
        case 18:
        case 23:
          yy.decorateNode({ class: $$[$0] });
          break;
        case 19:
          yy.getLogger().trace("SPACELIST");
          break;
        case 20:
          yy.getLogger().trace("Node: ", $$[$0 - 1].id);
          yy.addNode(0, $$[$0 - 1].id, $$[$0 - 1].descr, $$[$0 - 1].type, $$[$0]);
          break;
        case 21:
          yy.getLogger().trace("Node: ", $$[$0].id);
          yy.addNode(0, $$[$0].id, $$[$0].descr, $$[$0].type);
          break;
        case 22:
          yy.decorateNode({ icon: $$[$0] });
          break;
        case 27:
          yy.getLogger().trace("node found ..", $$[$0 - 2]);
          this.$ = { id: $$[$0 - 1], descr: $$[$0 - 1], type: yy.getType($$[$0 - 2], $$[$0]) };
          break;
        case 28:
          this.$ = { id: $$[$0], descr: $$[$0], type: 0 };
          break;
        case 29:
          yy.getLogger().trace("node found ..", $$[$0 - 3]);
          this.$ = { id: $$[$0 - 3], descr: $$[$0 - 1], type: yy.getType($$[$0 - 2], $$[$0]) };
          break;
        case 30:
          this.$ = $$[$0 - 1] + $$[$0];
          break;
        case 31:
          this.$ = $$[$0];
          break;
      }
    }, "anonymous"),
    table: [{ 3: 1, 4: 2, 5: 3, 6: [1, 5], 8: $V0 }, { 1: [3] }, { 1: [2, 1] }, { 4: 6, 6: [1, 7], 7: [1, 8], 8: $V0 }, { 6: $V1, 7: [1, 10], 9: 9, 12: 11, 13: $V2, 14: 14, 16: $V3, 17: $V4, 18: 17, 19: 18, 20: $V5, 23: $V6 }, o($V7, [2, 3]), { 1: [2, 2] }, o($V7, [2, 4]), o($V7, [2, 5]), { 1: [2, 6], 6: $V1, 12: 21, 13: $V2, 14: 14, 16: $V3, 17: $V4, 18: 17, 19: 18, 20: $V5, 23: $V6 }, { 6: $V1, 9: 22, 12: 11, 13: $V2, 14: 14, 16: $V3, 17: $V4, 18: 17, 19: 18, 20: $V5, 23: $V6 }, { 6: $V8, 7: $V9, 10: 23, 11: $Va }, o($Vb, [2, 24], { 18: 17, 19: 18, 14: 27, 16: [1, 28], 17: [1, 29], 20: $V5, 23: $V6 }), o($Vb, [2, 19]), o($Vb, [2, 21], { 15: 30, 24: $Vc }), o($Vb, [2, 22]), o($Vb, [2, 23]), o($Vd, [2, 25]), o($Vd, [2, 26]), o($Vd, [2, 28], { 20: [1, 32] }), { 21: [1, 33] }, { 6: $V8, 7: $V9, 10: 34, 11: $Va }, { 1: [2, 7], 6: $V1, 12: 21, 13: $V2, 14: 14, 16: $V3, 17: $V4, 18: 17, 19: 18, 20: $V5, 23: $V6 }, o($Ve, [2, 14], { 7: $Vf, 11: $Vg }), o($Vh, [2, 8]), o($Vh, [2, 9]), o($Vh, [2, 10]), o($Vb, [2, 16], { 15: 37, 24: $Vc }), o($Vb, [2, 17]), o($Vb, [2, 18]), o($Vb, [2, 20], { 24: $Vi }), o($Vd, [2, 31]), { 21: [1, 39] }, { 22: [1, 40] }, o($Ve, [2, 13], { 7: $Vf, 11: $Vg }), o($Vh, [2, 11]), o($Vh, [2, 12]), o($Vb, [2, 15], { 24: $Vi }), o($Vd, [2, 30]), { 22: [1, 41] }, o($Vd, [2, 27]), o($Vd, [2, 29])],
    defaultActions: { 2: [2, 1], 6: [2, 2] },
    parseError: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function parseError(str, hash) {
      if (hash.recoverable) {
        this.trace(str);
      } else {
        var error = new Error(str);
        error.hash = hash;
        throw error;
      }
    }, "parseError"),
    parse: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function parse(input) {
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
      (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(popStack, "popStack");
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
      (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(lex, "lex");
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
      parseError: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function parseError(str, hash) {
        if (this.yy.parser) {
          this.yy.parser.parseError(str, hash);
        } else {
          throw new Error(str);
        }
      }, "parseError"),
      // resets the lexer, sets new input
      setInput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function(input, yy) {
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
      input: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
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
      unput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function(ch) {
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
      more: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
        this._more = true;
        return this;
      }, "more"),
      // When called from action, signals the lexer that this rule fails to match the input, so the next matching rule (regex) should be tested instead.
      reject: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
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
      less: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function(n) {
        this.unput(this.match.slice(n));
      }, "less"),
      // displays already matched input, i.e. for error messages
      pastInput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
        var past = this.matched.substr(0, this.matched.length - this.match.length);
        return (past.length > 20 ? "..." : "") + past.substr(-20).replace(/\n/g, "");
      }, "pastInput"),
      // displays upcoming input, i.e. for error messages
      upcomingInput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
        var next = this.match;
        if (next.length < 20) {
          next += this._input.substr(0, 20 - next.length);
        }
        return (next.substr(0, 20) + (next.length > 20 ? "..." : "")).replace(/\n/g, "");
      }, "upcomingInput"),
      // displays the character position where the lexing error occurred, i.e. for error messages
      showPosition: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
        var pre = this.pastInput();
        var c = new Array(pre.length + 1).join("-");
        return pre + this.upcomingInput() + "\n" + c + "^";
      }, "showPosition"),
      // test the lexed token: return FALSE when not a match, otherwise return token
      test_match: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function(match, indexed_rule) {
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
      next: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
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
      lex: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function lex() {
        var r = this.next();
        if (r) {
          return r;
        } else {
          return this.lex();
        }
      }, "lex"),
      // activates a new lexer condition state (pushes the new lexer condition state onto the condition stack)
      begin: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function begin(condition) {
        this.conditionStack.push(condition);
      }, "begin"),
      // pop the previously active lexer condition state off the condition stack
      popState: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function popState() {
        var n = this.conditionStack.length - 1;
        if (n > 0) {
          return this.conditionStack.pop();
        } else {
          return this.conditionStack[0];
        }
      }, "popState"),
      // produce the lexer rule set which is active for the currently active lexer condition state
      _currentRules: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function _currentRules() {
        if (this.conditionStack.length && this.conditionStack[this.conditionStack.length - 1]) {
          return this.conditions[this.conditionStack[this.conditionStack.length - 1]].rules;
        } else {
          return this.conditions["INITIAL"].rules;
        }
      }, "_currentRules"),
      // return the currently active lexer condition state; when an index argument is provided it produces the N-th previous condition state, if available
      topState: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function topState(n) {
        n = this.conditionStack.length - 1 - Math.abs(n || 0);
        if (n >= 0) {
          return this.conditionStack[n];
        } else {
          return "INITIAL";
        }
      }, "topState"),
      // alias for begin(condition)
      pushState: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function pushState(condition) {
        this.begin(condition);
      }, "pushState"),
      // return the number of states currently on the stack
      stateStackSize: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function stateStackSize() {
        return this.conditionStack.length;
      }, "stateStackSize"),
      options: { "case-insensitive": true },
      performAction: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function anonymous(yy, yy_, $avoiding_name_collisions, YY_START) {
        var YYSTATE = YY_START;
        switch ($avoiding_name_collisions) {
          case 0:
            this.pushState("shapeData");
            yy_.yytext = "";
            return 24;
            break;
          case 1:
            this.pushState("shapeDataStr");
            return 24;
            break;
          case 2:
            this.popState();
            return 24;
            break;
          case 3:
            const re = /\n\s*/g;
            yy_.yytext = yy_.yytext.replace(re, "<br/>");
            return 24;
            break;
          case 4:
            return 24;
            break;
          case 5:
            this.popState();
            break;
          case 6:
            yy.getLogger().trace("Found comment", yy_.yytext);
            return 6;
            break;
          case 7:
            return 8;
            break;
          case 8:
            this.begin("CLASS");
            break;
          case 9:
            this.popState();
            return 17;
            break;
          case 10:
            this.popState();
            break;
          case 11:
            yy.getLogger().trace("Begin icon");
            this.begin("ICON");
            break;
          case 12:
            yy.getLogger().trace("SPACELINE");
            return 6;
            break;
          case 13:
            return 7;
            break;
          case 14:
            return 16;
            break;
          case 15:
            yy.getLogger().trace("end icon");
            this.popState();
            break;
          case 16:
            yy.getLogger().trace("Exploding node");
            this.begin("NODE");
            return 20;
            break;
          case 17:
            yy.getLogger().trace("Cloud");
            this.begin("NODE");
            return 20;
            break;
          case 18:
            yy.getLogger().trace("Explosion Bang");
            this.begin("NODE");
            return 20;
            break;
          case 19:
            yy.getLogger().trace("Cloud Bang");
            this.begin("NODE");
            return 20;
            break;
          case 20:
            this.begin("NODE");
            return 20;
            break;
          case 21:
            this.begin("NODE");
            return 20;
            break;
          case 22:
            this.begin("NODE");
            return 20;
            break;
          case 23:
            this.begin("NODE");
            return 20;
            break;
          case 24:
            return 13;
            break;
          case 25:
            return 23;
            break;
          case 26:
            return 11;
            break;
          case 27:
            this.begin("NSTR2");
            break;
          case 28:
            return "NODE_DESCR";
            break;
          case 29:
            this.popState();
            break;
          case 30:
            yy.getLogger().trace("Starting NSTR");
            this.begin("NSTR");
            break;
          case 31:
            yy.getLogger().trace("description:", yy_.yytext);
            return "NODE_DESCR";
            break;
          case 32:
            this.popState();
            break;
          case 33:
            this.popState();
            yy.getLogger().trace("node end ))");
            return "NODE_DEND";
            break;
          case 34:
            this.popState();
            yy.getLogger().trace("node end )");
            return "NODE_DEND";
            break;
          case 35:
            this.popState();
            yy.getLogger().trace("node end ...", yy_.yytext);
            return "NODE_DEND";
            break;
          case 36:
            this.popState();
            yy.getLogger().trace("node end ((");
            return "NODE_DEND";
            break;
          case 37:
            this.popState();
            yy.getLogger().trace("node end (-");
            return "NODE_DEND";
            break;
          case 38:
            this.popState();
            yy.getLogger().trace("node end (-");
            return "NODE_DEND";
            break;
          case 39:
            this.popState();
            yy.getLogger().trace("node end ((");
            return "NODE_DEND";
            break;
          case 40:
            this.popState();
            yy.getLogger().trace("node end ((");
            return "NODE_DEND";
            break;
          case 41:
            yy.getLogger().trace("Long description:", yy_.yytext);
            return 21;
            break;
          case 42:
            yy.getLogger().trace("Long description:", yy_.yytext);
            return 21;
            break;
        }
      }, "anonymous"),
      rules: [/^(?:@\{)/i, /^(?:["])/i, /^(?:["])/i, /^(?:[^\"]+)/i, /^(?:[^}^"]+)/i, /^(?:\})/i, /^(?:\s*%%.*)/i, /^(?:kanban\b)/i, /^(?::::)/i, /^(?:.+)/i, /^(?:\n)/i, /^(?:::icon\()/i, /^(?:[\s]+[\n])/i, /^(?:[\n]+)/i, /^(?:[^\)]+)/i, /^(?:\))/i, /^(?:-\))/i, /^(?:\(-)/i, /^(?:\)\))/i, /^(?:\))/i, /^(?:\(\()/i, /^(?:\{\{)/i, /^(?:\()/i, /^(?:\[)/i, /^(?:[\s]+)/i, /^(?:[^\(\[\n\)\{\}@]+)/i, /^(?:$)/i, /^(?:["][`])/i, /^(?:[^`"]+)/i, /^(?:[`]["])/i, /^(?:["])/i, /^(?:[^"]+)/i, /^(?:["])/i, /^(?:[\)]\))/i, /^(?:[\)])/i, /^(?:[\]])/i, /^(?:\}\})/i, /^(?:\(-)/i, /^(?:-\))/i, /^(?:\(\()/i, /^(?:\()/i, /^(?:[^\)\]\(\}]+)/i, /^(?:.+(?!\(\())/i],
      conditions: { "shapeDataEndBracket": { "rules": [], "inclusive": false }, "shapeDataStr": { "rules": [2, 3], "inclusive": false }, "shapeData": { "rules": [1, 4, 5], "inclusive": false }, "CLASS": { "rules": [9, 10], "inclusive": false }, "ICON": { "rules": [14, 15], "inclusive": false }, "NSTR2": { "rules": [28, 29], "inclusive": false }, "NSTR": { "rules": [31, 32], "inclusive": false }, "NODE": { "rules": [27, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42], "inclusive": false }, "INITIAL": { "rules": [0, 6, 7, 8, 11, 12, 13, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26], "inclusive": true } }
    };
    return lexer2;
  })();
  parser2.lexer = lexer;
  function Parser() {
    this.yy = {};
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(Parser, "Parser");
  Parser.prototype = parser2;
  parser2.Parser = Parser;
  return new Parser();
})();
parser.parser = parser;
var kanban_default = parser;

// src/diagrams/kanban/kanbanDb.ts
var nodes = [];
var sections = [];
var cnt = 0;
var elements = {};
var clear = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(() => {
  nodes = [];
  sections = [];
  cnt = 0;
  elements = {};
}, "clear");
var getSection = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)((level) => {
  if (nodes.length === 0) {
    return null;
  }
  const sectionLevel = nodes[0].level;
  let lastSection = null;
  for (let i = nodes.length - 1; i >= 0; i--) {
    if (nodes[i].level === sectionLevel && !lastSection) {
      lastSection = nodes[i];
    }
    if (nodes[i].level < sectionLevel) {
      throw new Error('Items without section detected, found section ("' + nodes[i].label + '")');
    }
  }
  if (level === lastSection?.level) {
    return null;
  }
  return lastSection;
}, "getSection");
var getSections = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
  return sections;
}, "getSections");
var getData = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(function() {
  const edges = [];
  const _nodes = [];
  const sections2 = getSections();
  const conf = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .getConfig2 */ .nV)();
  for (const section of sections2) {
    const node = {
      id: section.id,
      label: (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .sanitizeText */ .oO)(section.label ?? "", conf),
      isGroup: true,
      ticket: section.ticket,
      shape: "kanbanSection",
      level: section.level,
      look: conf.look
    };
    _nodes.push(node);
    const children = nodes.filter((n) => n.parentId === section.id);
    for (const item of children) {
      const childNode = {
        id: item.id,
        parentId: section.id,
        label: (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .sanitizeText */ .oO)(item.label ?? "", conf),
        isGroup: false,
        ticket: item?.ticket,
        priority: item?.priority,
        assigned: item?.assigned,
        icon: item?.icon,
        shape: "kanbanItem",
        level: item.level,
        rx: 5,
        ry: 5,
        cssStyles: ["text-align: left"]
      };
      _nodes.push(childNode);
    }
  }
  return { nodes: _nodes, edges, other: {}, config: (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .getConfig2 */ .nV)() };
}, "getData");
var addNode = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)((level, id, descr, type, shapeData) => {
  const conf = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .getConfig2 */ .nV)();
  let padding = conf.mindmap?.padding ?? _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .defaultConfig_default */ .vZ.mindmap.padding;
  switch (type) {
    case nodeType.ROUNDED_RECT:
    case nodeType.RECT:
    case nodeType.HEXAGON:
      padding *= 2;
  }
  const node = {
    id: (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .sanitizeText */ .oO)(id, conf) || "kbn" + cnt++,
    level,
    label: (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .sanitizeText */ .oO)(descr, conf),
    width: conf.mindmap?.maxNodeWidth ?? _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .defaultConfig_default */ .vZ.mindmap.maxNodeWidth,
    padding,
    isGroup: false
  };
  if (shapeData !== void 0) {
    let yamlData;
    if (!shapeData.includes("\n")) {
      yamlData = "{\n" + shapeData + "\n}";
    } else {
      yamlData = shapeData + "\n";
    }
    const doc = (0,_chunk_MI3HLSF2_mjs__WEBPACK_IMPORTED_MODULE_2__/* .load */ .z)(yamlData, { schema: _chunk_MI3HLSF2_mjs__WEBPACK_IMPORTED_MODULE_2__/* .JSON_SCHEMA */ .A });
    if (doc.shape && (doc.shape !== doc.shape.toLowerCase() || doc.shape.includes("_"))) {
      throw new Error(`No such shape: ${doc.shape}. Shape names should be lowercase.`);
    }
    if (doc?.shape && doc.shape === "kanbanItem") {
      node.shape = doc?.shape;
    }
    if (doc?.label) {
      node.label = doc?.label;
    }
    if (doc?.icon) {
      node.icon = doc?.icon.toString();
    }
    if (doc?.assigned) {
      node.assigned = doc?.assigned.toString();
    }
    if (doc?.ticket) {
      node.ticket = doc?.ticket.toString();
    }
    if (doc?.priority) {
      node.priority = doc?.priority;
    }
  }
  const section = getSection(level);
  if (section) {
    node.parentId = section.id || "kbn" + cnt++;
  } else {
    sections.push(node);
  }
  nodes.push(node);
}, "addNode");
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
var getType = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)((startStr, endStr) => {
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .log */ .cM.debug("In get type", startStr, endStr);
  switch (startStr) {
    case "[":
      return nodeType.RECT;
    case "(":
      return endStr === ")" ? nodeType.ROUNDED_RECT : nodeType.CLOUD;
    case "((":
      return nodeType.CIRCLE;
    case ")":
      return nodeType.CLOUD;
    case "))":
      return nodeType.BANG;
    case "{{":
      return nodeType.HEXAGON;
    default:
      return nodeType.DEFAULT;
  }
}, "getType");
var setElementForId = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)((id, element) => {
  elements[id] = element;
}, "setElementForId");
var decorateNode = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)((decoration) => {
  if (!decoration) {
    return;
  }
  const config = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .getConfig2 */ .nV)();
  const node = nodes[nodes.length - 1];
  if (decoration.icon) {
    node.icon = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .sanitizeText */ .oO)(decoration.icon, config);
  }
  if (decoration.class) {
    node.cssClasses = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .sanitizeText */ .oO)(decoration.class, config);
  }
}, "decorateNode");
var type2Str = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)((type) => {
  switch (type) {
    case nodeType.DEFAULT:
      return "no-border";
    case nodeType.RECT:
      return "rect";
    case nodeType.ROUNDED_RECT:
      return "rounded-rect";
    case nodeType.CIRCLE:
      return "circle";
    case nodeType.CLOUD:
      return "cloud";
    case nodeType.BANG:
      return "bang";
    case nodeType.HEXAGON:
      return "hexgon";
    // cspell: disable-line
    default:
      return "no-border";
  }
}, "type2Str");
var getLogger = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(() => _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .log */ .cM, "getLogger");
var getElementById = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)((id) => elements[id], "getElementById");
var db = {
  clear,
  addNode,
  getSections,
  getData,
  nodeType,
  getType,
  setElementForId,
  decorateNode,
  type2Str,
  getLogger,
  getElementById
};
var kanbanDb_default = db;

// src/diagrams/kanban/kanbanRenderer.ts
var draw = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)(async (text, id, _version, diagObj) => {
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .log */ .cM.debug("Rendering kanban diagram\n" + text);
  const db2 = diagObj.db;
  const data4Layout = db2.getData();
  const conf = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .getConfig2 */ .nV)();
  conf.htmlLabels = false;
  const svg = (0,_chunk_EXTU4WIE_mjs__WEBPACK_IMPORTED_MODULE_0__/* .selectSvgElement */ .P)(id);
  const sectionsElem = svg.append("g");
  sectionsElem.attr("class", "sections");
  const nodesElem = svg.append("g");
  nodesElem.attr("class", "items");
  const sections2 = data4Layout.nodes.filter(
    // TODO: TypeScript 5.5 will infer this predicate automatically
    (node) => node.isGroup
  );
  let cnt2 = 0;
  const padding = 10;
  const sectionObjects = [];
  let maxLabelHeight = 25;
  for (const section of sections2) {
    const WIDTH = conf?.kanban?.sectionWidth || 200;
    cnt2 = cnt2 + 1;
    section.x = WIDTH * cnt2 + (cnt2 - 1) * padding / 2;
    section.width = WIDTH;
    section.y = 0;
    section.height = WIDTH * 3;
    section.rx = 5;
    section.ry = 5;
    section.cssClasses = section.cssClasses + " section-" + cnt2;
    const sectionObj = await (0,_chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_3__/* .insertCluster */ .us)(sectionsElem, section);
    maxLabelHeight = Math.max(maxLabelHeight, sectionObj?.labelBBox?.height);
    sectionObjects.push(sectionObj);
  }
  let i = 0;
  for (const section of sections2) {
    const sectionObj = sectionObjects[i];
    i = i + 1;
    const WIDTH = conf?.kanban?.sectionWidth || 200;
    const top = -WIDTH * 3 / 2 + maxLabelHeight;
    let y = top;
    const sectionItems = data4Layout.nodes.filter((node) => node.parentId === section.id);
    for (const item of sectionItems) {
      if (item.isGroup) {
        throw new Error("Groups within groups are not allowed in Kanban diagrams");
      }
      item.x = section.x;
      item.width = WIDTH - 1.5 * padding;
      const nodeEl = await (0,_chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_3__/* .insertNode */ .Lf)(nodesElem, item, { config: conf });
      const bbox = nodeEl.node().getBBox();
      item.y = y + bbox.height / 2;
      await (0,_chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_3__/* .positionNode */ .aH)(item);
      y = item.y + bbox.height / 2 + padding / 2;
    }
    const rect = sectionObj.cluster.select("rect");
    const height = Math.max(y - top + 3 * padding, 50) + (maxLabelHeight - 25);
    rect.attr("height", height);
  }
  (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .setupGraphViewbox */ .j7)(
    void 0,
    svg,
    conf.mindmap?.padding ?? _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .defaultConfig_default */ .vZ.kanban.padding,
    conf.mindmap?.useMaxWidth ?? _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_8__/* .defaultConfig_default */ .vZ.kanban.useMaxWidth
  );
}, "draw");
var kanbanRenderer_default = {
  draw
};

// src/diagrams/kanban/styles.ts

var genSections = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)((options) => {
  let sections2 = "";
  for (let i = 0; i < options.THEME_COLOR_LIMIT; i++) {
    options["lineColor" + i] = options["lineColor" + i] || options["cScaleInv" + i];
    if ((0,khroma__WEBPACK_IMPORTED_MODULE_10__/* ["default"] */ .Z)(options["lineColor" + i])) {
      options["lineColor" + i] = (0,khroma__WEBPACK_IMPORTED_MODULE_11__/* ["default"] */ .Z)(options["lineColor" + i], 20);
    } else {
      options["lineColor" + i] = (0,khroma__WEBPACK_IMPORTED_MODULE_12__/* ["default"] */ .Z)(options["lineColor" + i], 20);
    }
  }
  const adjuster = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)((color, level) => options.darkMode ? (0,khroma__WEBPACK_IMPORTED_MODULE_12__/* ["default"] */ .Z)(color, level) : (0,khroma__WEBPACK_IMPORTED_MODULE_11__/* ["default"] */ .Z)(color, level), "adjuster");
  for (let i = 0; i < options.THEME_COLOR_LIMIT; i++) {
    const sw = "" + (17 - 3 * i);
    sections2 += `
    .section-${i - 1} rect, .section-${i - 1} path, .section-${i - 1} circle, .section-${i - 1} polygon, .section-${i - 1} path  {
      fill: ${adjuster(options["cScale" + i], 10)};
      stroke: ${adjuster(options["cScale" + i], 10)};

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

  .node rect,
  .node circle,
  .node ellipse,
  .node polygon,
  .node path {
    fill: ${options.background};
    stroke: ${options.nodeBorder};
    stroke-width: 1px;
  }

  .kanban-ticket-link {
    fill: ${options.background};
    stroke: ${options.nodeBorder};
    text-decoration: underline;
  }
    `;
  }
  return sections2;
}, "genSections");
var getStyles = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_9__/* .__name */ .eW)((options) => `
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
  .icon-container {
    height:100%;
    display: flex;
    justify-content: center;
    align-items: center;
  }
  .edge {
    fill: none;
  }
  .cluster-label, .label {
    color: ${options.textColor};
    fill: ${options.textColor};
    }
  .kanban-label {
    dy: 1em;
    alignment-baseline: middle;
    text-anchor: middle;
    dominant-baseline: middle;
    text-align: center;
  }
    ${(0,_chunk_FMBD7UC4_mjs__WEBPACK_IMPORTED_MODULE_1__/* .getIconStyles */ .G)()}
`, "getStyles");
var styles_default = getStyles;

// src/diagrams/kanban/kanban-definition.ts
var diagram = {
  db: kanbanDb_default,
  renderer: kanbanRenderer_default,
  parser: kanban_default,
  styles: styles_default
};



/***/ })

}]);
//# sourceMappingURL=369.5cecdf753e161a6bb7fe.js.map?v=5cecdf753e161a6bb7fe
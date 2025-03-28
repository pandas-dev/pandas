"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[5343],{

/***/ 75343:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   D: () => (/* binding */ DEFAULT_STATE_TYPE),
/* harmony export */   S: () => (/* binding */ STMT_RELATION),
/* harmony export */   a: () => (/* binding */ DIVIDER_TYPE),
/* harmony export */   b: () => (/* binding */ STMT_STATE),
/* harmony export */   c: () => (/* binding */ DEFAULT_NESTED_DOC_DIR),
/* harmony export */   d: () => (/* binding */ db),
/* harmony export */   p: () => (/* binding */ parser$1),
/* harmony export */   s: () => (/* binding */ styles)
/* harmony export */ });
/* harmony import */ var _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(24028);

var parser = function() {
  var o = function(k, v, o2, l) {
    for (o2 = o2 || {}, l = k.length; l--; o2[k[l]] = v)
      ;
    return o2;
  }, $V0 = [1, 2], $V1 = [1, 3], $V2 = [1, 4], $V3 = [2, 4], $V4 = [1, 9], $V5 = [1, 11], $V6 = [1, 15], $V7 = [1, 16], $V8 = [1, 17], $V9 = [1, 18], $Va = [1, 30], $Vb = [1, 19], $Vc = [1, 20], $Vd = [1, 21], $Ve = [1, 22], $Vf = [1, 23], $Vg = [1, 25], $Vh = [1, 26], $Vi = [1, 27], $Vj = [1, 28], $Vk = [1, 29], $Vl = [1, 32], $Vm = [1, 33], $Vn = [1, 34], $Vo = [1, 35], $Vp = [1, 31], $Vq = [1, 4, 5, 15, 16, 18, 20, 21, 23, 24, 25, 26, 27, 28, 32, 34, 36, 37, 41, 44, 45, 46, 47, 50], $Vr = [1, 4, 5, 13, 14, 15, 16, 18, 20, 21, 23, 24, 25, 26, 27, 28, 32, 34, 36, 37, 41, 44, 45, 46, 47, 50], $Vs = [4, 5, 15, 16, 18, 20, 21, 23, 24, 25, 26, 27, 28, 32, 34, 36, 37, 41, 44, 45, 46, 47, 50];
  var parser2 = {
    trace: function trace() {
    },
    yy: {},
    symbols_: { "error": 2, "start": 3, "SPACE": 4, "NL": 5, "SD": 6, "document": 7, "line": 8, "statement": 9, "classDefStatement": 10, "cssClassStatement": 11, "idStatement": 12, "DESCR": 13, "-->": 14, "HIDE_EMPTY": 15, "scale": 16, "WIDTH": 17, "COMPOSIT_STATE": 18, "STRUCT_START": 19, "STRUCT_STOP": 20, "STATE_DESCR": 21, "AS": 22, "ID": 23, "FORK": 24, "JOIN": 25, "CHOICE": 26, "CONCURRENT": 27, "note": 28, "notePosition": 29, "NOTE_TEXT": 30, "direction": 31, "acc_title": 32, "acc_title_value": 33, "acc_descr": 34, "acc_descr_value": 35, "acc_descr_multiline_value": 36, "classDef": 37, "CLASSDEF_ID": 38, "CLASSDEF_STYLEOPTS": 39, "DEFAULT": 40, "class": 41, "CLASSENTITY_IDS": 42, "STYLECLASS": 43, "direction_tb": 44, "direction_bt": 45, "direction_rl": 46, "direction_lr": 47, "eol": 48, ";": 49, "EDGE_STATE": 50, "STYLE_SEPARATOR": 51, "left_of": 52, "right_of": 53, "$accept": 0, "$end": 1 },
    terminals_: { 2: "error", 4: "SPACE", 5: "NL", 6: "SD", 13: "DESCR", 14: "-->", 15: "HIDE_EMPTY", 16: "scale", 17: "WIDTH", 18: "COMPOSIT_STATE", 19: "STRUCT_START", 20: "STRUCT_STOP", 21: "STATE_DESCR", 22: "AS", 23: "ID", 24: "FORK", 25: "JOIN", 26: "CHOICE", 27: "CONCURRENT", 28: "note", 30: "NOTE_TEXT", 32: "acc_title", 33: "acc_title_value", 34: "acc_descr", 35: "acc_descr_value", 36: "acc_descr_multiline_value", 37: "classDef", 38: "CLASSDEF_ID", 39: "CLASSDEF_STYLEOPTS", 40: "DEFAULT", 41: "class", 42: "CLASSENTITY_IDS", 43: "STYLECLASS", 44: "direction_tb", 45: "direction_bt", 46: "direction_rl", 47: "direction_lr", 49: ";", 50: "EDGE_STATE", 51: "STYLE_SEPARATOR", 52: "left_of", 53: "right_of" },
    productions_: [0, [3, 2], [3, 2], [3, 2], [7, 0], [7, 2], [8, 2], [8, 1], [8, 1], [9, 1], [9, 1], [9, 1], [9, 2], [9, 3], [9, 4], [9, 1], [9, 2], [9, 1], [9, 4], [9, 3], [9, 6], [9, 1], [9, 1], [9, 1], [9, 1], [9, 4], [9, 4], [9, 1], [9, 2], [9, 2], [9, 1], [10, 3], [10, 3], [11, 3], [31, 1], [31, 1], [31, 1], [31, 1], [48, 1], [48, 1], [12, 1], [12, 1], [12, 3], [12, 3], [29, 1], [29, 1]],
    performAction: function anonymous(yytext, yyleng, yylineno, yy, yystate, $$, _$) {
      var $0 = $$.length - 1;
      switch (yystate) {
        case 3:
          yy.setRootDoc($$[$0]);
          return $$[$0];
        case 4:
          this.$ = [];
          break;
        case 5:
          if ($$[$0] != "nl") {
            $$[$0 - 1].push($$[$0]);
            this.$ = $$[$0 - 1];
          }
          break;
        case 6:
        case 7:
          this.$ = $$[$0];
          break;
        case 8:
          this.$ = "nl";
          break;
        case 11:
          this.$ = $$[$0];
          break;
        case 12:
          const stateStmt = $$[$0 - 1];
          stateStmt.description = yy.trimColon($$[$0]);
          this.$ = stateStmt;
          break;
        case 13:
          this.$ = { stmt: "relation", state1: $$[$0 - 2], state2: $$[$0] };
          break;
        case 14:
          const relDescription = yy.trimColon($$[$0]);
          this.$ = { stmt: "relation", state1: $$[$0 - 3], state2: $$[$0 - 1], description: relDescription };
          break;
        case 18:
          this.$ = { stmt: "state", id: $$[$0 - 3], type: "default", description: "", doc: $$[$0 - 1] };
          break;
        case 19:
          var id = $$[$0];
          var description = $$[$0 - 2].trim();
          if ($$[$0].match(":")) {
            var parts = $$[$0].split(":");
            id = parts[0];
            description = [description, parts[1]];
          }
          this.$ = { stmt: "state", id, type: "default", description };
          break;
        case 20:
          this.$ = { stmt: "state", id: $$[$0 - 3], type: "default", description: $$[$0 - 5], doc: $$[$0 - 1] };
          break;
        case 21:
          this.$ = { stmt: "state", id: $$[$0], type: "fork" };
          break;
        case 22:
          this.$ = { stmt: "state", id: $$[$0], type: "join" };
          break;
        case 23:
          this.$ = { stmt: "state", id: $$[$0], type: "choice" };
          break;
        case 24:
          this.$ = { stmt: "state", id: yy.getDividerId(), type: "divider" };
          break;
        case 25:
          this.$ = { stmt: "state", id: $$[$0 - 1].trim(), note: { position: $$[$0 - 2].trim(), text: $$[$0].trim() } };
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
        case 32:
          this.$ = { stmt: "classDef", id: $$[$0 - 1].trim(), classes: $$[$0].trim() };
          break;
        case 33:
          this.$ = { stmt: "applyClass", id: $$[$0 - 1].trim(), styleClass: $$[$0].trim() };
          break;
        case 34:
          yy.setDirection("TB");
          this.$ = { stmt: "dir", value: "TB" };
          break;
        case 35:
          yy.setDirection("BT");
          this.$ = { stmt: "dir", value: "BT" };
          break;
        case 36:
          yy.setDirection("RL");
          this.$ = { stmt: "dir", value: "RL" };
          break;
        case 37:
          yy.setDirection("LR");
          this.$ = { stmt: "dir", value: "LR" };
          break;
        case 40:
        case 41:
          this.$ = { stmt: "state", id: $$[$0].trim(), type: "default", description: "" };
          break;
        case 42:
          this.$ = { stmt: "state", id: $$[$0 - 2].trim(), classes: [$$[$0].trim()], type: "default", description: "" };
          break;
        case 43:
          this.$ = { stmt: "state", id: $$[$0 - 2].trim(), classes: [$$[$0].trim()], type: "default", description: "" };
          break;
      }
    },
    table: [{ 3: 1, 4: $V0, 5: $V1, 6: $V2 }, { 1: [3] }, { 3: 5, 4: $V0, 5: $V1, 6: $V2 }, { 3: 6, 4: $V0, 5: $V1, 6: $V2 }, o([1, 4, 5, 15, 16, 18, 21, 23, 24, 25, 26, 27, 28, 32, 34, 36, 37, 41, 44, 45, 46, 47, 50], $V3, { 7: 7 }), { 1: [2, 1] }, { 1: [2, 2] }, { 1: [2, 3], 4: $V4, 5: $V5, 8: 8, 9: 10, 10: 12, 11: 13, 12: 14, 15: $V6, 16: $V7, 18: $V8, 21: $V9, 23: $Va, 24: $Vb, 25: $Vc, 26: $Vd, 27: $Ve, 28: $Vf, 31: 24, 32: $Vg, 34: $Vh, 36: $Vi, 37: $Vj, 41: $Vk, 44: $Vl, 45: $Vm, 46: $Vn, 47: $Vo, 50: $Vp }, o($Vq, [2, 5]), { 9: 36, 10: 12, 11: 13, 12: 14, 15: $V6, 16: $V7, 18: $V8, 21: $V9, 23: $Va, 24: $Vb, 25: $Vc, 26: $Vd, 27: $Ve, 28: $Vf, 31: 24, 32: $Vg, 34: $Vh, 36: $Vi, 37: $Vj, 41: $Vk, 44: $Vl, 45: $Vm, 46: $Vn, 47: $Vo, 50: $Vp }, o($Vq, [2, 7]), o($Vq, [2, 8]), o($Vq, [2, 9]), o($Vq, [2, 10]), o($Vq, [2, 11], { 13: [1, 37], 14: [1, 38] }), o($Vq, [2, 15]), { 17: [1, 39] }, o($Vq, [2, 17], { 19: [1, 40] }), { 22: [1, 41] }, o($Vq, [2, 21]), o($Vq, [2, 22]), o($Vq, [2, 23]), o($Vq, [2, 24]), { 29: 42, 30: [1, 43], 52: [1, 44], 53: [1, 45] }, o($Vq, [2, 27]), { 33: [1, 46] }, { 35: [1, 47] }, o($Vq, [2, 30]), { 38: [1, 48], 40: [1, 49] }, { 42: [1, 50] }, o($Vr, [2, 40], { 51: [1, 51] }), o($Vr, [2, 41], { 51: [1, 52] }), o($Vq, [2, 34]), o($Vq, [2, 35]), o($Vq, [2, 36]), o($Vq, [2, 37]), o($Vq, [2, 6]), o($Vq, [2, 12]), { 12: 53, 23: $Va, 50: $Vp }, o($Vq, [2, 16]), o($Vs, $V3, { 7: 54 }), { 23: [1, 55] }, { 23: [1, 56] }, { 22: [1, 57] }, { 23: [2, 44] }, { 23: [2, 45] }, o($Vq, [2, 28]), o($Vq, [2, 29]), { 39: [1, 58] }, { 39: [1, 59] }, { 43: [1, 60] }, { 23: [1, 61] }, { 23: [1, 62] }, o($Vq, [2, 13], { 13: [1, 63] }), { 4: $V4, 5: $V5, 8: 8, 9: 10, 10: 12, 11: 13, 12: 14, 15: $V6, 16: $V7, 18: $V8, 20: [1, 64], 21: $V9, 23: $Va, 24: $Vb, 25: $Vc, 26: $Vd, 27: $Ve, 28: $Vf, 31: 24, 32: $Vg, 34: $Vh, 36: $Vi, 37: $Vj, 41: $Vk, 44: $Vl, 45: $Vm, 46: $Vn, 47: $Vo, 50: $Vp }, o($Vq, [2, 19], { 19: [1, 65] }), { 30: [1, 66] }, { 23: [1, 67] }, o($Vq, [2, 31]), o($Vq, [2, 32]), o($Vq, [2, 33]), o($Vr, [2, 42]), o($Vr, [2, 43]), o($Vq, [2, 14]), o($Vq, [2, 18]), o($Vs, $V3, { 7: 68 }), o($Vq, [2, 25]), o($Vq, [2, 26]), { 4: $V4, 5: $V5, 8: 8, 9: 10, 10: 12, 11: 13, 12: 14, 15: $V6, 16: $V7, 18: $V8, 20: [1, 69], 21: $V9, 23: $Va, 24: $Vb, 25: $Vc, 26: $Vd, 27: $Ve, 28: $Vf, 31: 24, 32: $Vg, 34: $Vh, 36: $Vi, 37: $Vj, 41: $Vk, 44: $Vl, 45: $Vm, 46: $Vn, 47: $Vo, 50: $Vp }, o($Vq, [2, 20])],
    defaultActions: { 5: [2, 1], 6: [2, 2], 44: [2, 44], 45: [2, 45] },
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
            return 40;
          case 1:
            return 44;
          case 2:
            return 45;
          case 3:
            return 46;
          case 4:
            return 47;
          case 5:
            break;
          case 6:
            break;
          case 7:
            return 5;
          case 8:
            break;
          case 9:
            break;
          case 10:
            break;
          case 11:
            break;
          case 12:
            this.pushState("SCALE");
            return 16;
          case 13:
            return 17;
          case 14:
            this.popState();
            break;
          case 15:
            this.begin("acc_title");
            return 32;
          case 16:
            this.popState();
            return "acc_title_value";
          case 17:
            this.begin("acc_descr");
            return 34;
          case 18:
            this.popState();
            return "acc_descr_value";
          case 19:
            this.begin("acc_descr_multiline");
            break;
          case 20:
            this.popState();
            break;
          case 21:
            return "acc_descr_multiline_value";
          case 22:
            this.pushState("CLASSDEF");
            return 37;
          case 23:
            this.popState();
            this.pushState("CLASSDEFID");
            return "DEFAULT_CLASSDEF_ID";
          case 24:
            this.popState();
            this.pushState("CLASSDEFID");
            return 38;
          case 25:
            this.popState();
            return 39;
          case 26:
            this.pushState("CLASS");
            return 41;
          case 27:
            this.popState();
            this.pushState("CLASS_STYLE");
            return 42;
          case 28:
            this.popState();
            return 43;
          case 29:
            this.pushState("SCALE");
            return 16;
          case 30:
            return 17;
          case 31:
            this.popState();
            break;
          case 32:
            this.pushState("STATE");
            break;
          case 33:
            this.popState();
            yy_.yytext = yy_.yytext.slice(0, -8).trim();
            return 24;
          case 34:
            this.popState();
            yy_.yytext = yy_.yytext.slice(0, -8).trim();
            return 25;
          case 35:
            this.popState();
            yy_.yytext = yy_.yytext.slice(0, -10).trim();
            return 26;
          case 36:
            this.popState();
            yy_.yytext = yy_.yytext.slice(0, -8).trim();
            return 24;
          case 37:
            this.popState();
            yy_.yytext = yy_.yytext.slice(0, -8).trim();
            return 25;
          case 38:
            this.popState();
            yy_.yytext = yy_.yytext.slice(0, -10).trim();
            return 26;
          case 39:
            return 44;
          case 40:
            return 45;
          case 41:
            return 46;
          case 42:
            return 47;
          case 43:
            this.pushState("STATE_STRING");
            break;
          case 44:
            this.pushState("STATE_ID");
            return "AS";
          case 45:
            this.popState();
            return "ID";
          case 46:
            this.popState();
            break;
          case 47:
            return "STATE_DESCR";
          case 48:
            return 18;
          case 49:
            this.popState();
            break;
          case 50:
            this.popState();
            this.pushState("struct");
            return 19;
          case 51:
            break;
          case 52:
            this.popState();
            return 20;
          case 53:
            break;
          case 54:
            this.begin("NOTE");
            return 28;
          case 55:
            this.popState();
            this.pushState("NOTE_ID");
            return 52;
          case 56:
            this.popState();
            this.pushState("NOTE_ID");
            return 53;
          case 57:
            this.popState();
            this.pushState("FLOATING_NOTE");
            break;
          case 58:
            this.popState();
            this.pushState("FLOATING_NOTE_ID");
            return "AS";
          case 59:
            break;
          case 60:
            return "NOTE_TEXT";
          case 61:
            this.popState();
            return "ID";
          case 62:
            this.popState();
            this.pushState("NOTE_TEXT");
            return 23;
          case 63:
            this.popState();
            yy_.yytext = yy_.yytext.substr(2).trim();
            return 30;
          case 64:
            this.popState();
            yy_.yytext = yy_.yytext.slice(0, -8).trim();
            return 30;
          case 65:
            return 6;
          case 66:
            return 6;
          case 67:
            return 15;
          case 68:
            return 50;
          case 69:
            return 23;
          case 70:
            yy_.yytext = yy_.yytext.trim();
            return 13;
          case 71:
            return 14;
          case 72:
            return 27;
          case 73:
            return 51;
          case 74:
            return 5;
          case 75:
            return "INVALID";
        }
      },
      rules: [/^(?:default\b)/i, /^(?:.*direction\s+TB[^\n]*)/i, /^(?:.*direction\s+BT[^\n]*)/i, /^(?:.*direction\s+RL[^\n]*)/i, /^(?:.*direction\s+LR[^\n]*)/i, /^(?:%%(?!\{)[^\n]*)/i, /^(?:[^\}]%%[^\n]*)/i, /^(?:[\n]+)/i, /^(?:[\s]+)/i, /^(?:((?!\n)\s)+)/i, /^(?:#[^\n]*)/i, /^(?:%[^\n]*)/i, /^(?:scale\s+)/i, /^(?:\d+)/i, /^(?:\s+width\b)/i, /^(?:accTitle\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*\{\s*)/i, /^(?:[\}])/i, /^(?:[^\}]*)/i, /^(?:classDef\s+)/i, /^(?:DEFAULT\s+)/i, /^(?:\w+\s+)/i, /^(?:[^\n]*)/i, /^(?:class\s+)/i, /^(?:(\w+)+((,\s*\w+)*))/i, /^(?:[^\n]*)/i, /^(?:scale\s+)/i, /^(?:\d+)/i, /^(?:\s+width\b)/i, /^(?:state\s+)/i, /^(?:.*<<fork>>)/i, /^(?:.*<<join>>)/i, /^(?:.*<<choice>>)/i, /^(?:.*\[\[fork\]\])/i, /^(?:.*\[\[join\]\])/i, /^(?:.*\[\[choice\]\])/i, /^(?:.*direction\s+TB[^\n]*)/i, /^(?:.*direction\s+BT[^\n]*)/i, /^(?:.*direction\s+RL[^\n]*)/i, /^(?:.*direction\s+LR[^\n]*)/i, /^(?:["])/i, /^(?:\s*as\s+)/i, /^(?:[^\n\{]*)/i, /^(?:["])/i, /^(?:[^"]*)/i, /^(?:[^\n\s\{]+)/i, /^(?:\n)/i, /^(?:\{)/i, /^(?:%%(?!\{)[^\n]*)/i, /^(?:\})/i, /^(?:[\n])/i, /^(?:note\s+)/i, /^(?:left of\b)/i, /^(?:right of\b)/i, /^(?:")/i, /^(?:\s*as\s*)/i, /^(?:["])/i, /^(?:[^"]*)/i, /^(?:[^\n]*)/i, /^(?:\s*[^:\n\s\-]+)/i, /^(?:\s*:[^:\n;]+)/i, /^(?:[\s\S]*?end note\b)/i, /^(?:stateDiagram\s+)/i, /^(?:stateDiagram-v2\s+)/i, /^(?:hide empty description\b)/i, /^(?:\[\*\])/i, /^(?:[^:\n\s\-\{]+)/i, /^(?:\s*:[^:\n;]+)/i, /^(?:-->)/i, /^(?:--)/i, /^(?::::)/i, /^(?:$)/i, /^(?:.)/i],
      conditions: { "LINE": { "rules": [9, 10], "inclusive": false }, "struct": { "rules": [9, 10, 22, 26, 32, 39, 40, 41, 42, 51, 52, 53, 54, 68, 69, 70, 71, 72], "inclusive": false }, "FLOATING_NOTE_ID": { "rules": [61], "inclusive": false }, "FLOATING_NOTE": { "rules": [58, 59, 60], "inclusive": false }, "NOTE_TEXT": { "rules": [63, 64], "inclusive": false }, "NOTE_ID": { "rules": [62], "inclusive": false }, "NOTE": { "rules": [55, 56, 57], "inclusive": false }, "CLASS_STYLE": { "rules": [28], "inclusive": false }, "CLASS": { "rules": [27], "inclusive": false }, "CLASSDEFID": { "rules": [25], "inclusive": false }, "CLASSDEF": { "rules": [23, 24], "inclusive": false }, "acc_descr_multiline": { "rules": [20, 21], "inclusive": false }, "acc_descr": { "rules": [18], "inclusive": false }, "acc_title": { "rules": [16], "inclusive": false }, "SCALE": { "rules": [13, 14, 30, 31], "inclusive": false }, "ALIAS": { "rules": [], "inclusive": false }, "STATE_ID": { "rules": [45], "inclusive": false }, "STATE_STRING": { "rules": [46, 47], "inclusive": false }, "FORK_STATE": { "rules": [], "inclusive": false }, "STATE": { "rules": [9, 10, 33, 34, 35, 36, 37, 38, 43, 44, 48, 49, 50], "inclusive": false }, "ID": { "rules": [9, 10], "inclusive": false }, "INITIAL": { "rules": [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 15, 17, 19, 22, 26, 29, 32, 50, 54, 65, 66, 67, 68, 69, 70, 71, 73, 74, 75], "inclusive": true } }
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
const parser$1 = parser;
const DEFAULT_DIAGRAM_DIRECTION = "LR";
const DEFAULT_NESTED_DOC_DIR = "TB";
const STMT_STATE = "state";
const STMT_RELATION = "relation";
const STMT_CLASSDEF = "classDef";
const STMT_APPLYCLASS = "applyClass";
const DEFAULT_STATE_TYPE = "default";
const DIVIDER_TYPE = "divider";
const START_NODE = "[*]";
const START_TYPE = "start";
const END_NODE = START_NODE;
const END_TYPE = "end";
const COLOR_KEYWORD = "color";
const FILL_KEYWORD = "fill";
const BG_FILL = "bgFill";
const STYLECLASS_SEP = ",";
function newClassesList() {
  return {};
}
let direction = DEFAULT_DIAGRAM_DIRECTION;
let rootDoc = [];
let classes = newClassesList();
const newDoc = () => {
  return {
    relations: [],
    states: {},
    documents: {}
  };
};
let documents = {
  root: newDoc()
};
let currentDocument = documents.root;
let startEndCount = 0;
let dividerCnt = 0;
const lineType = {
  LINE: 0,
  DOTTED_LINE: 1
};
const relationType = {
  AGGREGATION: 0,
  EXTENSION: 1,
  COMPOSITION: 2,
  DEPENDENCY: 3
};
const clone = (o) => JSON.parse(JSON.stringify(o));
const setRootDoc = (o) => {
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.l.info("Setting root doc", o);
  rootDoc = o;
};
const getRootDoc = () => rootDoc;
const docTranslator = (parent, node, first) => {
  if (node.stmt === STMT_RELATION) {
    docTranslator(parent, node.state1, true);
    docTranslator(parent, node.state2, false);
  } else {
    if (node.stmt === STMT_STATE) {
      if (node.id === "[*]") {
        node.id = first ? parent.id + "_start" : parent.id + "_end";
        node.start = first;
      } else {
        node.id = node.id.trim();
      }
    }
    if (node.doc) {
      const doc = [];
      let currentDoc = [];
      let i;
      for (i = 0; i < node.doc.length; i++) {
        if (node.doc[i].type === DIVIDER_TYPE) {
          const newNode = clone(node.doc[i]);
          newNode.doc = clone(currentDoc);
          doc.push(newNode);
          currentDoc = [];
        } else {
          currentDoc.push(node.doc[i]);
        }
      }
      if (doc.length > 0 && currentDoc.length > 0) {
        const newNode = {
          stmt: STMT_STATE,
          id: (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.G)(),
          type: "divider",
          doc: clone(currentDoc)
        };
        doc.push(clone(newNode));
        node.doc = doc;
      }
      node.doc.forEach((docNode) => docTranslator(node, docNode, true));
    }
  }
};
const getRootDocV2 = () => {
  docTranslator({ id: "root" }, { id: "root", doc: rootDoc }, true);
  return { id: "root", doc: rootDoc };
};
const extract = (_doc) => {
  let doc;
  if (_doc.doc) {
    doc = _doc.doc;
  } else {
    doc = _doc;
  }
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.l.info(doc);
  clear(true);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.l.info("Extract", doc);
  doc.forEach((item) => {
    switch (item.stmt) {
      case STMT_STATE:
        addState(
          item.id.trim(),
          item.type,
          item.doc,
          item.description,
          item.note,
          item.classes,
          item.styles,
          item.textStyles
        );
        break;
      case STMT_RELATION:
        addRelation(item.state1, item.state2, item.description);
        break;
      case STMT_CLASSDEF:
        addStyleClass(item.id.trim(), item.classes);
        break;
      case STMT_APPLYCLASS:
        setCssClass(item.id.trim(), item.styleClass);
        break;
    }
  });
};
const addState = function(id, type = DEFAULT_STATE_TYPE, doc = null, descr = null, note = null, classes2 = null, styles2 = null, textStyles = null) {
  const trimmedId = id == null ? void 0 : id.trim();
  if (currentDocument.states[trimmedId] === void 0) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.l.info("Adding state ", trimmedId, descr);
    currentDocument.states[trimmedId] = {
      id: trimmedId,
      descriptions: [],
      type,
      doc,
      note,
      classes: [],
      styles: [],
      textStyles: []
    };
  } else {
    if (!currentDocument.states[trimmedId].doc) {
      currentDocument.states[trimmedId].doc = doc;
    }
    if (!currentDocument.states[trimmedId].type) {
      currentDocument.states[trimmedId].type = type;
    }
  }
  if (descr) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.l.info("Setting state description", trimmedId, descr);
    if (typeof descr === "string") {
      addDescription(trimmedId, descr.trim());
    }
    if (typeof descr === "object") {
      descr.forEach((des) => addDescription(trimmedId, des.trim()));
    }
  }
  if (note) {
    currentDocument.states[trimmedId].note = note;
    currentDocument.states[trimmedId].note.text = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.e.sanitizeText(
      currentDocument.states[trimmedId].note.text,
      (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.c)()
    );
  }
  if (classes2) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.l.info("Setting state classes", trimmedId, classes2);
    const classesList = typeof classes2 === "string" ? [classes2] : classes2;
    classesList.forEach((klass) => setCssClass(trimmedId, klass.trim()));
  }
  if (styles2) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.l.info("Setting state styles", trimmedId, styles2);
    const stylesList = typeof styles2 === "string" ? [styles2] : styles2;
    stylesList.forEach((style) => setStyle(trimmedId, style.trim()));
  }
  if (textStyles) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.l.info("Setting state styles", trimmedId, styles2);
    const textStylesList = typeof textStyles === "string" ? [textStyles] : textStyles;
    textStylesList.forEach((textStyle) => setTextStyle(trimmedId, textStyle.trim()));
  }
};
const clear = function(saveCommon) {
  documents = {
    root: newDoc()
  };
  currentDocument = documents.root;
  startEndCount = 0;
  classes = newClassesList();
  if (!saveCommon) {
    (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.t)();
  }
};
const getState = function(id) {
  return currentDocument.states[id];
};
const getStates = function() {
  return currentDocument.states;
};
const logDocuments = function() {
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.l.info("Documents = ", documents);
};
const getRelations = function() {
  return currentDocument.relations;
};
function startIdIfNeeded(id = "") {
  let fixedId = id;
  if (id === START_NODE) {
    startEndCount++;
    fixedId = `${START_TYPE}${startEndCount}`;
  }
  return fixedId;
}
function startTypeIfNeeded(id = "", type = DEFAULT_STATE_TYPE) {
  return id === START_NODE ? START_TYPE : type;
}
function endIdIfNeeded(id = "") {
  let fixedId = id;
  if (id === END_NODE) {
    startEndCount++;
    fixedId = `${END_TYPE}${startEndCount}`;
  }
  return fixedId;
}
function endTypeIfNeeded(id = "", type = DEFAULT_STATE_TYPE) {
  return id === END_NODE ? END_TYPE : type;
}
function addRelationObjs(item1, item2, relationTitle) {
  let id1 = startIdIfNeeded(item1.id.trim());
  let type1 = startTypeIfNeeded(item1.id.trim(), item1.type);
  let id2 = startIdIfNeeded(item2.id.trim());
  let type2 = startTypeIfNeeded(item2.id.trim(), item2.type);
  addState(
    id1,
    type1,
    item1.doc,
    item1.description,
    item1.note,
    item1.classes,
    item1.styles,
    item1.textStyles
  );
  addState(
    id2,
    type2,
    item2.doc,
    item2.description,
    item2.note,
    item2.classes,
    item2.styles,
    item2.textStyles
  );
  currentDocument.relations.push({
    id1,
    id2,
    relationTitle: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.e.sanitizeText(relationTitle, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.c)())
  });
}
const addRelation = function(item1, item2, title) {
  if (typeof item1 === "object") {
    addRelationObjs(item1, item2, title);
  } else {
    const id1 = startIdIfNeeded(item1.trim());
    const type1 = startTypeIfNeeded(item1);
    const id2 = endIdIfNeeded(item2.trim());
    const type2 = endTypeIfNeeded(item2);
    addState(id1, type1);
    addState(id2, type2);
    currentDocument.relations.push({
      id1,
      id2,
      title: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.e.sanitizeText(title, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.c)())
    });
  }
};
const addDescription = function(id, descr) {
  const theState = currentDocument.states[id];
  const _descr = descr.startsWith(":") ? descr.replace(":", "").trim() : descr;
  theState.descriptions.push(_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.e.sanitizeText(_descr, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.c)()));
};
const cleanupLabel = function(label) {
  if (label.substring(0, 1) === ":") {
    return label.substr(2).trim();
  } else {
    return label.trim();
  }
};
const getDividerId = () => {
  dividerCnt++;
  return "divider-id-" + dividerCnt;
};
const addStyleClass = function(id, styleAttributes = "") {
  if (classes[id] === void 0) {
    classes[id] = { id, styles: [], textStyles: [] };
  }
  const foundClass = classes[id];
  if (styleAttributes !== void 0 && styleAttributes !== null) {
    styleAttributes.split(STYLECLASS_SEP).forEach((attrib) => {
      const fixedAttrib = attrib.replace(/([^;]*);/, "$1").trim();
      if (attrib.match(COLOR_KEYWORD)) {
        const newStyle1 = fixedAttrib.replace(FILL_KEYWORD, BG_FILL);
        const newStyle2 = newStyle1.replace(COLOR_KEYWORD, FILL_KEYWORD);
        foundClass.textStyles.push(newStyle2);
      }
      foundClass.styles.push(fixedAttrib);
    });
  }
};
const getClasses = function() {
  return classes;
};
const setCssClass = function(itemIds, cssClassName) {
  itemIds.split(",").forEach(function(id) {
    let foundState = getState(id);
    if (foundState === void 0) {
      const trimmedId = id.trim();
      addState(trimmedId);
      foundState = getState(trimmedId);
    }
    foundState.classes.push(cssClassName);
  });
};
const setStyle = function(itemId, styleText) {
  const item = getState(itemId);
  if (item !== void 0) {
    item.textStyles.push(styleText);
  }
};
const setTextStyle = function(itemId, cssClassName) {
  const item = getState(itemId);
  if (item !== void 0) {
    item.textStyles.push(cssClassName);
  }
};
const getDirection = () => direction;
const setDirection = (dir) => {
  direction = dir;
};
const trimColon = (str) => str && str[0] === ":" ? str.substr(1).trim() : str.trim();
const db = {
  getConfig: () => (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.c)().state,
  addState,
  clear,
  getState,
  getStates,
  getRelations,
  getClasses,
  getDirection,
  addRelation,
  getDividerId,
  setDirection,
  cleanupLabel,
  lineType,
  relationType,
  logDocuments,
  getRootDoc,
  setRootDoc,
  getRootDocV2,
  extract,
  trimColon,
  getAccTitle: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.g,
  setAccTitle: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.s,
  getAccDescription: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.a,
  setAccDescription: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.b,
  addStyleClass,
  setCssClass,
  addDescription,
  setDiagramTitle: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.q,
  getDiagramTitle: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_0__.r
};
const getStyles = (options) => `
defs #statediagram-barbEnd {
    fill: ${options.transitionColor};
    stroke: ${options.transitionColor};
  }
g.stateGroup text {
  fill: ${options.nodeBorder};
  stroke: none;
  font-size: 10px;
}
g.stateGroup text {
  fill: ${options.textColor};
  stroke: none;
  font-size: 10px;

}
g.stateGroup .state-title {
  font-weight: bolder;
  fill: ${options.stateLabelColor};
}

g.stateGroup rect {
  fill: ${options.mainBkg};
  stroke: ${options.nodeBorder};
}

g.stateGroup line {
  stroke: ${options.lineColor};
  stroke-width: 1;
}

.transition {
  stroke: ${options.transitionColor};
  stroke-width: 1;
  fill: none;
}

.stateGroup .composit {
  fill: ${options.background};
  border-bottom: 1px
}

.stateGroup .alt-composit {
  fill: #e0e0e0;
  border-bottom: 1px
}

.state-note {
  stroke: ${options.noteBorderColor};
  fill: ${options.noteBkgColor};

  text {
    fill: ${options.noteTextColor};
    stroke: none;
    font-size: 10px;
  }
}

.stateLabel .box {
  stroke: none;
  stroke-width: 0;
  fill: ${options.mainBkg};
  opacity: 0.5;
}

.edgeLabel .label rect {
  fill: ${options.labelBackgroundColor};
  opacity: 0.5;
}
.edgeLabel .label text {
  fill: ${options.transitionLabelColor || options.tertiaryTextColor};
}
.label div .edgeLabel {
  color: ${options.transitionLabelColor || options.tertiaryTextColor};
}

.stateLabel text {
  fill: ${options.stateLabelColor};
  font-size: 10px;
  font-weight: bold;
}

.node circle.state-start {
  fill: ${options.specialStateColor};
  stroke: ${options.specialStateColor};
}

.node .fork-join {
  fill: ${options.specialStateColor};
  stroke: ${options.specialStateColor};
}

.node circle.state-end {
  fill: ${options.innerEndBackground};
  stroke: ${options.background};
  stroke-width: 1.5
}
.end-state-inner {
  fill: ${options.compositeBackground || options.background};
  // stroke: ${options.background};
  stroke-width: 1.5
}

.node rect {
  fill: ${options.stateBkg || options.mainBkg};
  stroke: ${options.stateBorder || options.nodeBorder};
  stroke-width: 1px;
}
.node polygon {
  fill: ${options.mainBkg};
  stroke: ${options.stateBorder || options.nodeBorder};;
  stroke-width: 1px;
}
#statediagram-barbEnd {
  fill: ${options.lineColor};
}

.statediagram-cluster rect {
  fill: ${options.compositeTitleBackground};
  stroke: ${options.stateBorder || options.nodeBorder};
  stroke-width: 1px;
}

.cluster-label, .nodeLabel {
  color: ${options.stateLabelColor};
}

.statediagram-cluster rect.outer {
  rx: 5px;
  ry: 5px;
}
.statediagram-state .divider {
  stroke: ${options.stateBorder || options.nodeBorder};
}

.statediagram-state .title-state {
  rx: 5px;
  ry: 5px;
}
.statediagram-cluster.statediagram-cluster .inner {
  fill: ${options.compositeBackground || options.background};
}
.statediagram-cluster.statediagram-cluster-alt .inner {
  fill: ${options.altBackground ? options.altBackground : "#efefef"};
}

.statediagram-cluster .inner {
  rx:0;
  ry:0;
}

.statediagram-state rect.basic {
  rx: 5px;
  ry: 5px;
}
.statediagram-state rect.divider {
  stroke-dasharray: 10,10;
  fill: ${options.altBackground ? options.altBackground : "#efefef"};
}

.note-edge {
  stroke-dasharray: 5;
}

.statediagram-note rect {
  fill: ${options.noteBkgColor};
  stroke: ${options.noteBorderColor};
  stroke-width: 1px;
  rx: 0;
  ry: 0;
}
.statediagram-note rect {
  fill: ${options.noteBkgColor};
  stroke: ${options.noteBorderColor};
  stroke-width: 1px;
  rx: 0;
  ry: 0;
}

.statediagram-note text {
  fill: ${options.noteTextColor};
}

.statediagram-note .nodeLabel {
  color: ${options.noteTextColor};
}
.statediagram .edgeLabel {
  color: red; // ${options.noteTextColor};
}

#dependencyStart, #dependencyEnd {
  fill: ${options.lineColor};
  stroke: ${options.lineColor};
  stroke-width: 1;
}

.statediagramTitleText {
  text-anchor: middle;
  font-size: 18px;
  fill: ${options.textColor};
}
`;
const styles = getStyles;



/***/ })

}]);
//# sourceMappingURL=5343.9b1635b1e24801d5d086.js.map?v=9b1635b1e24801d5d086
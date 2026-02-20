"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[1225],{

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

/***/ 52315:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   AD: () => (/* binding */ getTextObj),
/* harmony export */   AE: () => (/* binding */ drawImage),
/* harmony export */   Mu: () => (/* binding */ drawRect),
/* harmony export */   O: () => (/* binding */ drawBackgroundRect),
/* harmony export */   kc: () => (/* binding */ getNoteRect),
/* harmony export */   rB: () => (/* binding */ drawEmbeddedImage),
/* harmony export */   yU: () => (/* binding */ drawText)
/* harmony export */ });
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(74999);
/* harmony import */ var _braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(7608);



// src/diagrams/common/svgDrawCommon.ts

var drawRect = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((element, rectData) => {
  const rectElement = element.append("rect");
  rectElement.attr("x", rectData.x);
  rectElement.attr("y", rectData.y);
  rectElement.attr("fill", rectData.fill);
  rectElement.attr("stroke", rectData.stroke);
  rectElement.attr("width", rectData.width);
  rectElement.attr("height", rectData.height);
  if (rectData.name) {
    rectElement.attr("name", rectData.name);
  }
  if (rectData.rx) {
    rectElement.attr("rx", rectData.rx);
  }
  if (rectData.ry) {
    rectElement.attr("ry", rectData.ry);
  }
  if (rectData.attrs !== void 0) {
    for (const attrKey in rectData.attrs) {
      rectElement.attr(attrKey, rectData.attrs[attrKey]);
    }
  }
  if (rectData.class) {
    rectElement.attr("class", rectData.class);
  }
  return rectElement;
}, "drawRect");
var drawBackgroundRect = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((element, bounds) => {
  const rectData = {
    x: bounds.startx,
    y: bounds.starty,
    width: bounds.stopx - bounds.startx,
    height: bounds.stopy - bounds.starty,
    fill: bounds.fill,
    stroke: bounds.stroke,
    class: "rect"
  };
  const rectElement = drawRect(element, rectData);
  rectElement.lower();
}, "drawBackgroundRect");
var drawText = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((element, textData) => {
  const nText = textData.text.replace(_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .lineBreakRegex */ .Vw, " ");
  const textElem = element.append("text");
  textElem.attr("x", textData.x);
  textElem.attr("y", textData.y);
  textElem.attr("class", "legend");
  textElem.style("text-anchor", textData.anchor);
  if (textData.class) {
    textElem.attr("class", textData.class);
  }
  const tspan = textElem.append("tspan");
  tspan.attr("x", textData.x + textData.textMargin * 2);
  tspan.text(nText);
  return textElem;
}, "drawText");
var drawImage = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((elem, x, y, link) => {
  const imageElement = elem.append("image");
  imageElement.attr("x", x);
  imageElement.attr("y", y);
  const sanitizedLink = (0,_braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_2__/* .sanitizeUrl */ .N)(link);
  imageElement.attr("xlink:href", sanitizedLink);
}, "drawImage");
var drawEmbeddedImage = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((element, x, y, link) => {
  const imageElement = element.append("use");
  imageElement.attr("x", x);
  imageElement.attr("y", y);
  const sanitizedLink = (0,_braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_2__/* .sanitizeUrl */ .N)(link);
  imageElement.attr("xlink:href", `#${sanitizedLink}`);
}, "drawEmbeddedImage");
var getNoteRect = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(() => {
  const noteRectData = {
    x: 0,
    y: 0,
    width: 100,
    height: 100,
    fill: "#EDF2AE",
    stroke: "#666",
    anchor: "start",
    rx: 0,
    ry: 0
  };
  return noteRectData;
}, "getNoteRect");
var getTextObj = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(() => {
  const testObject = {
    x: 0,
    y: 0,
    width: 100,
    height: 100,
    "text-anchor": "start",
    style: "#666",
    textMargin: 0,
    rx: 0,
    ry: 0,
    tspan: true
  };
  return testObject;
}, "getTextObj");




/***/ }),

/***/ 34885:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _chunk_TZMSLE5B_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(52315);
/* harmony import */ var _chunk_FMBD7UC4_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(28862);
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(74999);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(35321);





// src/diagrams/user-journey/parser/journey.jison
var parser = (function() {
  var o = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(k, v, o2, l) {
    for (o2 = o2 || {}, l = k.length; l--; o2[k[l]] = v) ;
    return o2;
  }, "o"), $V0 = [6, 8, 10, 11, 12, 14, 16, 17, 18], $V1 = [1, 9], $V2 = [1, 10], $V3 = [1, 11], $V4 = [1, 12], $V5 = [1, 13], $V6 = [1, 14];
  var parser2 = {
    trace: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function trace() {
    }, "trace"),
    yy: {},
    symbols_: { "error": 2, "start": 3, "journey": 4, "document": 5, "EOF": 6, "line": 7, "SPACE": 8, "statement": 9, "NEWLINE": 10, "title": 11, "acc_title": 12, "acc_title_value": 13, "acc_descr": 14, "acc_descr_value": 15, "acc_descr_multiline_value": 16, "section": 17, "taskName": 18, "taskData": 19, "$accept": 0, "$end": 1 },
    terminals_: { 2: "error", 4: "journey", 6: "EOF", 8: "SPACE", 10: "NEWLINE", 11: "title", 12: "acc_title", 13: "acc_title_value", 14: "acc_descr", 15: "acc_descr_value", 16: "acc_descr_multiline_value", 17: "section", 18: "taskName", 19: "taskData" },
    productions_: [0, [3, 3], [5, 0], [5, 2], [7, 2], [7, 1], [7, 1], [7, 1], [9, 1], [9, 2], [9, 2], [9, 1], [9, 1], [9, 2]],
    performAction: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function anonymous(yytext, yyleng, yylineno, yy, yystate, $$, _$) {
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
          yy.setDiagramTitle($$[$0].substr(6));
          this.$ = $$[$0].substr(6);
          break;
        case 9:
          this.$ = $$[$0].trim();
          yy.setAccTitle(this.$);
          break;
        case 10:
        case 11:
          this.$ = $$[$0].trim();
          yy.setAccDescription(this.$);
          break;
        case 12:
          yy.addSection($$[$0].substr(8));
          this.$ = $$[$0].substr(8);
          break;
        case 13:
          yy.addTask($$[$0 - 1], $$[$0]);
          this.$ = "task";
          break;
      }
    }, "anonymous"),
    table: [{ 3: 1, 4: [1, 2] }, { 1: [3] }, o($V0, [2, 2], { 5: 3 }), { 6: [1, 4], 7: 5, 8: [1, 6], 9: 7, 10: [1, 8], 11: $V1, 12: $V2, 14: $V3, 16: $V4, 17: $V5, 18: $V6 }, o($V0, [2, 7], { 1: [2, 1] }), o($V0, [2, 3]), { 9: 15, 11: $V1, 12: $V2, 14: $V3, 16: $V4, 17: $V5, 18: $V6 }, o($V0, [2, 5]), o($V0, [2, 6]), o($V0, [2, 8]), { 13: [1, 16] }, { 15: [1, 17] }, o($V0, [2, 11]), o($V0, [2, 12]), { 19: [1, 18] }, o($V0, [2, 4]), o($V0, [2, 9]), o($V0, [2, 10]), o($V0, [2, 13])],
    defaultActions: {},
    parseError: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function parseError(str, hash) {
      if (hash.recoverable) {
        this.trace(str);
      } else {
        var error = new Error(str);
        error.hash = hash;
        throw error;
      }
    }, "parseError"),
    parse: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function parse(input) {
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
      (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(popStack, "popStack");
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
      (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(lex, "lex");
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
      parseError: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function parseError(str, hash) {
        if (this.yy.parser) {
          this.yy.parser.parseError(str, hash);
        } else {
          throw new Error(str);
        }
      }, "parseError"),
      // resets the lexer, sets new input
      setInput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(input, yy) {
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
      input: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
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
      unput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(ch) {
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
      more: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        this._more = true;
        return this;
      }, "more"),
      // When called from action, signals the lexer that this rule fails to match the input, so the next matching rule (regex) should be tested instead.
      reject: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
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
      less: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(n) {
        this.unput(this.match.slice(n));
      }, "less"),
      // displays already matched input, i.e. for error messages
      pastInput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        var past = this.matched.substr(0, this.matched.length - this.match.length);
        return (past.length > 20 ? "..." : "") + past.substr(-20).replace(/\n/g, "");
      }, "pastInput"),
      // displays upcoming input, i.e. for error messages
      upcomingInput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        var next = this.match;
        if (next.length < 20) {
          next += this._input.substr(0, 20 - next.length);
        }
        return (next.substr(0, 20) + (next.length > 20 ? "..." : "")).replace(/\n/g, "");
      }, "upcomingInput"),
      // displays the character position where the lexing error occurred, i.e. for error messages
      showPosition: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        var pre = this.pastInput();
        var c = new Array(pre.length + 1).join("-");
        return pre + this.upcomingInput() + "\n" + c + "^";
      }, "showPosition"),
      // test the lexed token: return FALSE when not a match, otherwise return token
      test_match: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(match, indexed_rule) {
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
      next: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
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
      lex: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function lex() {
        var r = this.next();
        if (r) {
          return r;
        } else {
          return this.lex();
        }
      }, "lex"),
      // activates a new lexer condition state (pushes the new lexer condition state onto the condition stack)
      begin: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function begin(condition) {
        this.conditionStack.push(condition);
      }, "begin"),
      // pop the previously active lexer condition state off the condition stack
      popState: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function popState() {
        var n = this.conditionStack.length - 1;
        if (n > 0) {
          return this.conditionStack.pop();
        } else {
          return this.conditionStack[0];
        }
      }, "popState"),
      // produce the lexer rule set which is active for the currently active lexer condition state
      _currentRules: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function _currentRules() {
        if (this.conditionStack.length && this.conditionStack[this.conditionStack.length - 1]) {
          return this.conditions[this.conditionStack[this.conditionStack.length - 1]].rules;
        } else {
          return this.conditions["INITIAL"].rules;
        }
      }, "_currentRules"),
      // return the currently active lexer condition state; when an index argument is provided it produces the N-th previous condition state, if available
      topState: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function topState(n) {
        n = this.conditionStack.length - 1 - Math.abs(n || 0);
        if (n >= 0) {
          return this.conditionStack[n];
        } else {
          return "INITIAL";
        }
      }, "topState"),
      // alias for begin(condition)
      pushState: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function pushState(condition) {
        this.begin(condition);
      }, "pushState"),
      // return the number of states currently on the stack
      stateStackSize: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function stateStackSize() {
        return this.conditionStack.length;
      }, "stateStackSize"),
      options: { "case-insensitive": true },
      performAction: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function anonymous(yy, yy_, $avoiding_name_collisions, YY_START) {
        var YYSTATE = YY_START;
        switch ($avoiding_name_collisions) {
          case 0:
            break;
          case 1:
            break;
          case 2:
            return 10;
            break;
          case 3:
            break;
          case 4:
            break;
          case 5:
            return 4;
            break;
          case 6:
            return 11;
            break;
          case 7:
            this.begin("acc_title");
            return 12;
            break;
          case 8:
            this.popState();
            return "acc_title_value";
            break;
          case 9:
            this.begin("acc_descr");
            return 14;
            break;
          case 10:
            this.popState();
            return "acc_descr_value";
            break;
          case 11:
            this.begin("acc_descr_multiline");
            break;
          case 12:
            this.popState();
            break;
          case 13:
            return "acc_descr_multiline_value";
            break;
          case 14:
            return 17;
            break;
          case 15:
            return 18;
            break;
          case 16:
            return 19;
            break;
          case 17:
            return ":";
            break;
          case 18:
            return 6;
            break;
          case 19:
            return "INVALID";
            break;
        }
      }, "anonymous"),
      rules: [/^(?:%(?!\{)[^\n]*)/i, /^(?:[^\}]%%[^\n]*)/i, /^(?:[\n]+)/i, /^(?:\s+)/i, /^(?:#[^\n]*)/i, /^(?:journey\b)/i, /^(?:title\s[^#\n;]+)/i, /^(?:accTitle\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*\{\s*)/i, /^(?:[\}])/i, /^(?:[^\}]*)/i, /^(?:section\s[^#:\n;]+)/i, /^(?:[^#:\n;]+)/i, /^(?::[^#\n;]+)/i, /^(?::)/i, /^(?:$)/i, /^(?:.)/i],
      conditions: { "acc_descr_multiline": { "rules": [12, 13], "inclusive": false }, "acc_descr": { "rules": [10], "inclusive": false }, "acc_title": { "rules": [8], "inclusive": false }, "INITIAL": { "rules": [0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 14, 15, 16, 17, 18, 19], "inclusive": true } }
    };
    return lexer2;
  })();
  parser2.lexer = lexer;
  function Parser() {
    this.yy = {};
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(Parser, "Parser");
  Parser.prototype = parser2;
  parser2.Parser = Parser;
  return new Parser();
})();
parser.parser = parser;
var journey_default = parser;

// src/diagrams/user-journey/journeyDb.js
var currentSection = "";
var sections = [];
var tasks = [];
var rawTasks = [];
var clear2 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
  sections.length = 0;
  tasks.length = 0;
  currentSection = "";
  rawTasks.length = 0;
  (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .clear */ .ZH)();
}, "clear");
var addSection = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(txt) {
  currentSection = txt;
  sections.push(txt);
}, "addSection");
var getSections = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
  return sections;
}, "getSections");
var getTasks = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
  let allItemsProcessed = compileTasks();
  const maxDepth = 100;
  let iterationCount = 0;
  while (!allItemsProcessed && iterationCount < maxDepth) {
    allItemsProcessed = compileTasks();
    iterationCount++;
  }
  tasks.push(...rawTasks);
  return tasks;
}, "getTasks");
var updateActors = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
  const tempActors = [];
  tasks.forEach((task) => {
    if (task.people) {
      tempActors.push(...task.people);
    }
  });
  const unique = new Set(tempActors);
  return [...unique].sort();
}, "updateActors");
var addTask = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(descr, taskData) {
  const pieces = taskData.substr(1).split(":");
  let score = 0;
  let peeps = [];
  if (pieces.length === 1) {
    score = Number(pieces[0]);
    peeps = [];
  } else {
    score = Number(pieces[0]);
    peeps = pieces[1].split(",");
  }
  const peopleList = peeps.map((s) => s.trim());
  const rawTask = {
    section: currentSection,
    type: currentSection,
    people: peopleList,
    task: descr,
    score
  };
  rawTasks.push(rawTask);
}, "addTask");
var addTaskOrg = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(descr) {
  const newTask = {
    section: currentSection,
    type: currentSection,
    description: descr,
    task: descr,
    classes: []
  };
  tasks.push(newTask);
}, "addTaskOrg");
var compileTasks = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
  const compileTask = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(pos) {
    return rawTasks[pos].processed;
  }, "compileTask");
  let allProcessed = true;
  for (const [i, rawTask] of rawTasks.entries()) {
    compileTask(i);
    allProcessed = allProcessed && rawTask.processed;
  }
  return allProcessed;
}, "compileTasks");
var getActors = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
  return updateActors();
}, "getActors");
var journeyDb_default = {
  getConfig: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getConfig2 */ .nV)().journey, "getConfig"),
  clear: clear2,
  setDiagramTitle: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .setDiagramTitle */ .g2,
  getDiagramTitle: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getDiagramTitle */ .Kr,
  setAccTitle: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .setAccTitle */ .GN,
  getAccTitle: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getAccTitle */ .eu,
  setAccDescription: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .setAccDescription */ .U$,
  getAccDescription: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getAccDescription */ .Mx,
  addSection,
  getSections,
  getTasks,
  addTask,
  addTaskOrg,
  getActors
};

// src/diagrams/user-journey/styles.js
var getStyles = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((options) => `.label {
    font-family: ${options.fontFamily};
    color: ${options.textColor};
  }
  .mouth {
    stroke: #666;
  }

  line {
    stroke: ${options.textColor}
  }

  .legend {
    fill: ${options.textColor};
    font-family: ${options.fontFamily};
  }

  .label text {
    fill: #333;
  }
  .label {
    color: ${options.textColor}
  }

  .face {
    ${options.faceColor ? `fill: ${options.faceColor}` : "fill: #FFF8DC"};
    stroke: #999;
  }

  .node rect,
  .node circle,
  .node ellipse,
  .node polygon,
  .node path {
    fill: ${options.mainBkg};
    stroke: ${options.nodeBorder};
    stroke-width: 1px;
  }

  .node .label {
    text-align: center;
  }
  .node.clickable {
    cursor: pointer;
  }

  .arrowheadPath {
    fill: ${options.arrowheadColor};
  }

  .edgePath .path {
    stroke: ${options.lineColor};
    stroke-width: 1.5px;
  }

  .flowchart-link {
    stroke: ${options.lineColor};
    fill: none;
  }

  .edgeLabel {
    background-color: ${options.edgeLabelBackground};
    rect {
      opacity: 0.5;
    }
    text-align: center;
  }

  .cluster rect {
  }

  .cluster text {
    fill: ${options.titleColor};
  }

  div.mermaidTooltip {
    position: absolute;
    text-align: center;
    max-width: 200px;
    padding: 2px;
    font-family: ${options.fontFamily};
    font-size: 12px;
    background: ${options.tertiaryColor};
    border: 1px solid ${options.border2};
    border-radius: 2px;
    pointer-events: none;
    z-index: 100;
  }

  .task-type-0, .section-type-0  {
    ${options.fillType0 ? `fill: ${options.fillType0}` : ""};
  }
  .task-type-1, .section-type-1  {
    ${options.fillType0 ? `fill: ${options.fillType1}` : ""};
  }
  .task-type-2, .section-type-2  {
    ${options.fillType0 ? `fill: ${options.fillType2}` : ""};
  }
  .task-type-3, .section-type-3  {
    ${options.fillType0 ? `fill: ${options.fillType3}` : ""};
  }
  .task-type-4, .section-type-4  {
    ${options.fillType0 ? `fill: ${options.fillType4}` : ""};
  }
  .task-type-5, .section-type-5  {
    ${options.fillType0 ? `fill: ${options.fillType5}` : ""};
  }
  .task-type-6, .section-type-6  {
    ${options.fillType0 ? `fill: ${options.fillType6}` : ""};
  }
  .task-type-7, .section-type-7  {
    ${options.fillType0 ? `fill: ${options.fillType7}` : ""};
  }

  .actor-0 {
    ${options.actor0 ? `fill: ${options.actor0}` : ""};
  }
  .actor-1 {
    ${options.actor1 ? `fill: ${options.actor1}` : ""};
  }
  .actor-2 {
    ${options.actor2 ? `fill: ${options.actor2}` : ""};
  }
  .actor-3 {
    ${options.actor3 ? `fill: ${options.actor3}` : ""};
  }
  .actor-4 {
    ${options.actor4 ? `fill: ${options.actor4}` : ""};
  }
  .actor-5 {
    ${options.actor5 ? `fill: ${options.actor5}` : ""};
  }
  ${(0,_chunk_FMBD7UC4_mjs__WEBPACK_IMPORTED_MODULE_1__/* .getIconStyles */ .G)()}
`, "getStyles");
var styles_default = getStyles;

// src/diagrams/user-journey/journeyRenderer.ts


// src/diagrams/user-journey/svgDraw.js

var drawRect2 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, rectData) {
  return (0,_chunk_TZMSLE5B_mjs__WEBPACK_IMPORTED_MODULE_0__/* .drawRect */ .Mu)(elem, rectData);
}, "drawRect");
var drawFace = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(element, faceData) {
  const radius = 15;
  const circleElement = element.append("circle").attr("cx", faceData.cx).attr("cy", faceData.cy).attr("class", "face").attr("r", radius).attr("stroke-width", 2).attr("overflow", "visible");
  const face = element.append("g");
  face.append("circle").attr("cx", faceData.cx - radius / 3).attr("cy", faceData.cy - radius / 3).attr("r", 1.5).attr("stroke-width", 2).attr("fill", "#666").attr("stroke", "#666");
  face.append("circle").attr("cx", faceData.cx + radius / 3).attr("cy", faceData.cy - radius / 3).attr("r", 1.5).attr("stroke-width", 2).attr("fill", "#666").attr("stroke", "#666");
  function smile(face2) {
    const arc = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .arc */ .Nb1)().startAngle(Math.PI / 2).endAngle(3 * (Math.PI / 2)).innerRadius(radius / 2).outerRadius(radius / 2.2);
    face2.append("path").attr("class", "mouth").attr("d", arc).attr("transform", "translate(" + faceData.cx + "," + (faceData.cy + 2) + ")");
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(smile, "smile");
  function sad(face2) {
    const arc = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .arc */ .Nb1)().startAngle(3 * Math.PI / 2).endAngle(5 * (Math.PI / 2)).innerRadius(radius / 2).outerRadius(radius / 2.2);
    face2.append("path").attr("class", "mouth").attr("d", arc).attr("transform", "translate(" + faceData.cx + "," + (faceData.cy + 7) + ")");
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(sad, "sad");
  function ambivalent(face2) {
    face2.append("line").attr("class", "mouth").attr("stroke", 2).attr("x1", faceData.cx - 5).attr("y1", faceData.cy + 7).attr("x2", faceData.cx + 5).attr("y2", faceData.cy + 7).attr("class", "mouth").attr("stroke-width", "1px").attr("stroke", "#666");
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(ambivalent, "ambivalent");
  if (faceData.score > 3) {
    smile(face);
  } else if (faceData.score < 3) {
    sad(face);
  } else {
    ambivalent(face);
  }
  return circleElement;
}, "drawFace");
var drawCircle = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(element, circleData) {
  const circleElement = element.append("circle");
  circleElement.attr("cx", circleData.cx);
  circleElement.attr("cy", circleData.cy);
  circleElement.attr("class", "actor-" + circleData.pos);
  circleElement.attr("fill", circleData.fill);
  circleElement.attr("stroke", circleData.stroke);
  circleElement.attr("r", circleData.r);
  if (circleElement.class !== void 0) {
    circleElement.attr("class", circleElement.class);
  }
  if (circleData.title !== void 0) {
    circleElement.append("title").text(circleData.title);
  }
  return circleElement;
}, "drawCircle");
var drawText2 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, textData) {
  return (0,_chunk_TZMSLE5B_mjs__WEBPACK_IMPORTED_MODULE_0__/* .drawText */ .yU)(elem, textData);
}, "drawText");
var drawLabel = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, txtObject) {
  function genPoints(x, y, width, height, cut) {
    return x + "," + y + " " + (x + width) + "," + y + " " + (x + width) + "," + (y + height - cut) + " " + (x + width - cut * 1.2) + "," + (y + height) + " " + x + "," + (y + height);
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(genPoints, "genPoints");
  const polygon = elem.append("polygon");
  polygon.attr("points", genPoints(txtObject.x, txtObject.y, 50, 20, 7));
  polygon.attr("class", "labelBox");
  txtObject.y = txtObject.y + txtObject.labelMargin;
  txtObject.x = txtObject.x + 0.5 * txtObject.labelMargin;
  drawText2(elem, txtObject);
}, "drawLabel");
var drawSection = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, section, conf2) {
  const g = elem.append("g");
  const rect = (0,_chunk_TZMSLE5B_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getNoteRect */ .kc)();
  rect.x = section.x;
  rect.y = section.y;
  rect.fill = section.fill;
  rect.width = conf2.width * section.taskCount + // width of the tasks
  conf2.diagramMarginX * (section.taskCount - 1);
  rect.height = conf2.height;
  rect.class = "journey-section section-type-" + section.num;
  rect.rx = 3;
  rect.ry = 3;
  drawRect2(g, rect);
  _drawTextCandidateFunc(conf2)(
    section.text,
    g,
    rect.x,
    rect.y,
    rect.width,
    rect.height,
    { class: "journey-section section-type-" + section.num },
    conf2,
    section.colour
  );
}, "drawSection");
var taskCount = -1;
var drawTask = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, task, conf2) {
  const center = task.x + conf2.width / 2;
  const g = elem.append("g");
  taskCount++;
  const maxHeight = 300 + 5 * 30;
  g.append("line").attr("id", "task" + taskCount).attr("x1", center).attr("y1", task.y).attr("x2", center).attr("y2", maxHeight).attr("class", "task-line").attr("stroke-width", "1px").attr("stroke-dasharray", "4 2").attr("stroke", "#666");
  drawFace(g, {
    cx: center,
    cy: 300 + (5 - task.score) * 30,
    score: task.score
  });
  const rect = (0,_chunk_TZMSLE5B_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getNoteRect */ .kc)();
  rect.x = task.x;
  rect.y = task.y;
  rect.fill = task.fill;
  rect.width = conf2.width;
  rect.height = conf2.height;
  rect.class = "task task-type-" + task.num;
  rect.rx = 3;
  rect.ry = 3;
  drawRect2(g, rect);
  let xPos = task.x + 14;
  task.people.forEach((person) => {
    const colour = task.actors[person].color;
    const circle = {
      cx: xPos,
      cy: task.y,
      r: 7,
      fill: colour,
      stroke: "#000",
      title: person,
      pos: task.actors[person].position
    };
    drawCircle(g, circle);
    xPos += 10;
  });
  _drawTextCandidateFunc(conf2)(
    task.task,
    g,
    rect.x,
    rect.y,
    rect.width,
    rect.height,
    { class: "task" },
    conf2,
    task.colour
  );
}, "drawTask");
var drawBackgroundRect2 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(elem, bounds2) {
  (0,_chunk_TZMSLE5B_mjs__WEBPACK_IMPORTED_MODULE_0__/* .drawBackgroundRect */ .O)(elem, bounds2);
}, "drawBackgroundRect");
var _drawTextCandidateFunc = /* @__PURE__ */ (function() {
  function byText(content, g, x, y, width, height, textAttrs, colour) {
    const text = g.append("text").attr("x", x + width / 2).attr("y", y + height / 2 + 5).style("font-color", colour).style("text-anchor", "middle").text(content);
    _setTextAttrs(text, textAttrs);
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(byText, "byText");
  function byTspan(content, g, x, y, width, height, textAttrs, conf2, colour) {
    const { taskFontSize, taskFontFamily } = conf2;
    const lines = content.split(/<br\s*\/?>/gi);
    for (let i = 0; i < lines.length; i++) {
      const dy = i * taskFontSize - taskFontSize * (lines.length - 1) / 2;
      const text = g.append("text").attr("x", x + width / 2).attr("y", y).attr("fill", colour).style("text-anchor", "middle").style("font-size", taskFontSize).style("font-family", taskFontFamily);
      text.append("tspan").attr("x", x + width / 2).attr("dy", dy).text(lines[i]);
      text.attr("y", y + height / 2).attr("dominant-baseline", "central").attr("alignment-baseline", "central");
      _setTextAttrs(text, textAttrs);
    }
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(byTspan, "byTspan");
  function byFo(content, g, x, y, width, height, textAttrs, conf2) {
    const body = g.append("switch");
    const f = body.append("foreignObject").attr("x", x).attr("y", y).attr("width", width).attr("height", height).attr("position", "fixed");
    const text = f.append("xhtml:div").style("display", "table").style("height", "100%").style("width", "100%");
    text.append("div").attr("class", "label").style("display", "table-cell").style("text-align", "center").style("vertical-align", "middle").text(content);
    byTspan(content, body, x, y, width, height, textAttrs, conf2);
    _setTextAttrs(text, textAttrs);
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(byFo, "byFo");
  function _setTextAttrs(toText, fromTextAttrsDict) {
    for (const key in fromTextAttrsDict) {
      if (key in fromTextAttrsDict) {
        toText.attr(key, fromTextAttrsDict[key]);
      }
    }
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(_setTextAttrs, "_setTextAttrs");
  return function(conf2) {
    return conf2.textPlacement === "fo" ? byFo : conf2.textPlacement === "old" ? byText : byTspan;
  };
})();
var initGraphics = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(graphics) {
  graphics.append("defs").append("marker").attr("id", "arrowhead").attr("refX", 5).attr("refY", 2).attr("markerWidth", 6).attr("markerHeight", 4).attr("orient", "auto").append("path").attr("d", "M 0,0 V 4 L6,2 Z");
}, "initGraphics");
var svgDraw_default = {
  drawRect: drawRect2,
  drawCircle,
  drawSection,
  drawText: drawText2,
  drawLabel,
  drawTask,
  drawBackgroundRect: drawBackgroundRect2,
  initGraphics
};

// src/diagrams/user-journey/journeyRenderer.ts
var setConf = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(cnf) {
  const keys = Object.keys(cnf);
  keys.forEach(function(key) {
    conf[key] = cnf[key];
  });
}, "setConf");
var actors = {};
var maxWidth = 0;
function drawActorLegend(diagram2) {
  const conf2 = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getConfig2 */ .nV)().journey;
  const maxLabelWidth = conf2.maxLabelWidth;
  maxWidth = 0;
  let yPos = 60;
  Object.keys(actors).forEach((person) => {
    const colour = actors[person].color;
    const circleData = {
      cx: 20,
      cy: yPos,
      r: 7,
      fill: colour,
      stroke: "#000",
      pos: actors[person].position
    };
    svgDraw_default.drawCircle(diagram2, circleData);
    let measureText = diagram2.append("text").attr("visibility", "hidden").text(person);
    const fullTextWidth = measureText.node().getBoundingClientRect().width;
    measureText.remove();
    let lines = [];
    if (fullTextWidth <= maxLabelWidth) {
      lines = [person];
    } else {
      const words = person.split(" ");
      let currentLine = "";
      measureText = diagram2.append("text").attr("visibility", "hidden");
      words.forEach((word) => {
        const testLine = currentLine ? `${currentLine} ${word}` : word;
        measureText.text(testLine);
        const textWidth = measureText.node().getBoundingClientRect().width;
        if (textWidth > maxLabelWidth) {
          if (currentLine) {
            lines.push(currentLine);
          }
          currentLine = word;
          measureText.text(word);
          if (measureText.node().getBoundingClientRect().width > maxLabelWidth) {
            let brokenWord = "";
            for (const char of word) {
              brokenWord += char;
              measureText.text(brokenWord + "-");
              if (measureText.node().getBoundingClientRect().width > maxLabelWidth) {
                lines.push(brokenWord.slice(0, -1) + "-");
                brokenWord = char;
              }
            }
            currentLine = brokenWord;
          }
        } else {
          currentLine = testLine;
        }
      });
      if (currentLine) {
        lines.push(currentLine);
      }
      measureText.remove();
    }
    lines.forEach((line, index) => {
      const labelData = {
        x: 40,
        y: yPos + 7 + index * 20,
        fill: "#666",
        text: line,
        textMargin: conf2.boxTextMargin ?? 5
      };
      const textElement = svgDraw_default.drawText(diagram2, labelData);
      const lineWidth = textElement.node().getBoundingClientRect().width;
      if (lineWidth > maxWidth && lineWidth > conf2.leftMargin - lineWidth) {
        maxWidth = lineWidth;
      }
    });
    yPos += Math.max(20, lines.length * 20);
  });
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(drawActorLegend, "drawActorLegend");
var conf = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getConfig2 */ .nV)().journey;
var leftMargin = 0;
var draw = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(text, id, version, diagObj) {
  const configObject = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getConfig2 */ .nV)();
  const titleColor = configObject.journey.titleColor;
  const titleFontSize = configObject.journey.titleFontSize;
  const titleFontFamily = configObject.journey.titleFontFamily;
  const securityLevel = configObject.securityLevel;
  let sandboxElement;
  if (securityLevel === "sandbox") {
    sandboxElement = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .select */ .Ys)("#i" + id);
  }
  const root = securityLevel === "sandbox" ? (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .select */ .Ys)(sandboxElement.nodes()[0].contentDocument.body) : (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .select */ .Ys)("body");
  bounds.init();
  const diagram2 = root.select("#" + id);
  svgDraw_default.initGraphics(diagram2);
  const tasks2 = diagObj.db.getTasks();
  const title = diagObj.db.getDiagramTitle();
  const actorNames = diagObj.db.getActors();
  for (const member in actors) {
    delete actors[member];
  }
  let actorPos = 0;
  actorNames.forEach((actorName) => {
    actors[actorName] = {
      color: conf.actorColours[actorPos % conf.actorColours.length],
      position: actorPos
    };
    actorPos++;
  });
  drawActorLegend(diagram2);
  leftMargin = conf.leftMargin + maxWidth;
  bounds.insert(0, 0, leftMargin, Object.keys(actors).length * 50);
  drawTasks(diagram2, tasks2, 0);
  const box = bounds.getBounds();
  if (title) {
    diagram2.append("text").text(title).attr("x", leftMargin).attr("font-size", titleFontSize).attr("font-weight", "bold").attr("y", 25).attr("fill", titleColor).attr("font-family", titleFontFamily);
  }
  const height = box.stopy - box.starty + 2 * conf.diagramMarginY;
  const width = leftMargin + box.stopx + 2 * conf.diagramMarginX;
  (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .configureSvgSize */ .v2)(diagram2, height, width, conf.useMaxWidth);
  diagram2.append("line").attr("x1", leftMargin).attr("y1", conf.height * 4).attr("x2", width - leftMargin - 4).attr("y2", conf.height * 4).attr("stroke-width", 4).attr("stroke", "black").attr("marker-end", "url(#arrowhead)");
  const extraVertForTitle = title ? 70 : 0;
  diagram2.attr("viewBox", `${box.startx} -25 ${width} ${height + extraVertForTitle}`);
  diagram2.attr("preserveAspectRatio", "xMinYMin meet");
  diagram2.attr("height", height + extraVertForTitle + 25);
}, "draw");
var bounds = {
  data: {
    startx: void 0,
    stopx: void 0,
    starty: void 0,
    stopy: void 0
  },
  verticalPos: 0,
  sequenceItems: [],
  init: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
    this.sequenceItems = [];
    this.data = {
      startx: void 0,
      stopx: void 0,
      starty: void 0,
      stopy: void 0
    };
    this.verticalPos = 0;
  }, "init"),
  updateVal: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(obj, key, val, fun) {
    if (obj[key] === void 0) {
      obj[key] = val;
    } else {
      obj[key] = fun(val, obj[key]);
    }
  }, "updateVal"),
  updateBounds: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(startx, starty, stopx, stopy) {
    const conf2 = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getConfig2 */ .nV)().journey;
    const _self = this;
    let cnt = 0;
    function updateFn(type) {
      return /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function updateItemBounds(item) {
        cnt++;
        const n = _self.sequenceItems.length - cnt + 1;
        _self.updateVal(item, "starty", starty - n * conf2.boxMargin, Math.min);
        _self.updateVal(item, "stopy", stopy + n * conf2.boxMargin, Math.max);
        _self.updateVal(bounds.data, "startx", startx - n * conf2.boxMargin, Math.min);
        _self.updateVal(bounds.data, "stopx", stopx + n * conf2.boxMargin, Math.max);
        if (!(type === "activation")) {
          _self.updateVal(item, "startx", startx - n * conf2.boxMargin, Math.min);
          _self.updateVal(item, "stopx", stopx + n * conf2.boxMargin, Math.max);
          _self.updateVal(bounds.data, "starty", starty - n * conf2.boxMargin, Math.min);
          _self.updateVal(bounds.data, "stopy", stopy + n * conf2.boxMargin, Math.max);
        }
      }, "updateItemBounds");
    }
    (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(updateFn, "updateFn");
    this.sequenceItems.forEach(updateFn());
  }, "updateBounds"),
  insert: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(startx, starty, stopx, stopy) {
    const _startx = Math.min(startx, stopx);
    const _stopx = Math.max(startx, stopx);
    const _starty = Math.min(starty, stopy);
    const _stopy = Math.max(starty, stopy);
    this.updateVal(bounds.data, "startx", _startx, Math.min);
    this.updateVal(bounds.data, "starty", _starty, Math.min);
    this.updateVal(bounds.data, "stopx", _stopx, Math.max);
    this.updateVal(bounds.data, "stopy", _stopy, Math.max);
    this.updateBounds(_startx, _starty, _stopx, _stopy);
  }, "insert"),
  bumpVerticalPos: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(bump) {
    this.verticalPos = this.verticalPos + bump;
    this.data.stopy = this.verticalPos;
  }, "bumpVerticalPos"),
  getVerticalPos: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
    return this.verticalPos;
  }, "getVerticalPos"),
  getBounds: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
    return this.data;
  }, "getBounds")
};
var fills = conf.sectionFills;
var textColours = conf.sectionColours;
var drawTasks = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(diagram2, tasks2, verticalPos) {
  const conf2 = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getConfig2 */ .nV)().journey;
  let lastSection = "";
  const sectionVHeight = conf2.height * 2 + conf2.diagramMarginY;
  const taskPos = verticalPos + sectionVHeight;
  let sectionNumber = 0;
  let fill = "#CCC";
  let colour = "black";
  let num = 0;
  for (const [i, task] of tasks2.entries()) {
    if (lastSection !== task.section) {
      fill = fills[sectionNumber % fills.length];
      num = sectionNumber % fills.length;
      colour = textColours[sectionNumber % textColours.length];
      let taskInSectionCount = 0;
      const currentSection2 = task.section;
      for (let taskIndex = i; taskIndex < tasks2.length; taskIndex++) {
        if (tasks2[taskIndex].section == currentSection2) {
          taskInSectionCount = taskInSectionCount + 1;
        } else {
          break;
        }
      }
      const section = {
        x: i * conf2.taskMargin + i * conf2.width + leftMargin,
        y: 50,
        text: task.section,
        fill,
        num,
        colour,
        taskCount: taskInSectionCount
      };
      svgDraw_default.drawSection(diagram2, section, conf2);
      lastSection = task.section;
      sectionNumber++;
    }
    const taskActors = task.people.reduce((acc, actorName) => {
      if (actors[actorName]) {
        acc[actorName] = actors[actorName];
      }
      return acc;
    }, {});
    task.x = i * conf2.taskMargin + i * conf2.width + leftMargin;
    task.y = taskPos;
    task.width = conf2.diagramMarginX;
    task.height = conf2.diagramMarginY;
    task.colour = colour;
    task.fill = fill;
    task.num = num;
    task.actors = taskActors;
    svgDraw_default.drawTask(diagram2, task, conf2);
    bounds.insert(task.x, task.y, task.x + task.width + conf2.taskMargin, 300 + 5 * 30);
  }
}, "drawTasks");
var journeyRenderer_default = {
  setConf,
  draw
};

// src/diagrams/user-journey/journeyDiagram.ts
var diagram = {
  parser: journey_default,
  db: journeyDb_default,
  renderer: journeyRenderer_default,
  styles: styles_default,
  init: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((cnf) => {
    journeyRenderer_default.setConf(cnf.journey);
    journeyDb_default.clear();
  }, "init")
};



/***/ })

}]);
//# sourceMappingURL=1225.a84f9ad316be9c1538e1.js.map?v=a84f9ad316be9c1538e1
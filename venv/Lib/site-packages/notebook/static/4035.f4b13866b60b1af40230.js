"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[4035],{

/***/ 74035:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(24028);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(23617);
/* harmony import */ var _svgDrawCommon_92a8ff2b_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(26476);
/* harmony import */ var dayjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(27693);
/* harmony import */ var dayjs__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(dayjs__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(7608);
/* harmony import */ var dompurify__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(31699);
/* harmony import */ var dompurify__WEBPACK_IMPORTED_MODULE_3___default = /*#__PURE__*/__webpack_require__.n(dompurify__WEBPACK_IMPORTED_MODULE_3__);












var parser = function() {
  var o = function(k, v, o2, l) {
    for (o2 = o2 || {}, l = k.length; l--; o2[k[l]] = v)
      ;
    return o2;
  }, $V0 = [6, 8, 10, 11, 12, 14, 16, 17, 18], $V1 = [1, 9], $V2 = [1, 10], $V3 = [1, 11], $V4 = [1, 12], $V5 = [1, 13], $V6 = [1, 14];
  var parser2 = {
    trace: function trace() {
    },
    yy: {},
    symbols_: { "error": 2, "start": 3, "journey": 4, "document": 5, "EOF": 6, "line": 7, "SPACE": 8, "statement": 9, "NEWLINE": 10, "title": 11, "acc_title": 12, "acc_title_value": 13, "acc_descr": 14, "acc_descr_value": 15, "acc_descr_multiline_value": 16, "section": 17, "taskName": 18, "taskData": 19, "$accept": 0, "$end": 1 },
    terminals_: { 2: "error", 4: "journey", 6: "EOF", 8: "SPACE", 10: "NEWLINE", 11: "title", 12: "acc_title", 13: "acc_title_value", 14: "acc_descr", 15: "acc_descr_value", 16: "acc_descr_multiline_value", 17: "section", 18: "taskName", 19: "taskData" },
    productions_: [0, [3, 3], [5, 0], [5, 2], [7, 2], [7, 1], [7, 1], [7, 1], [9, 1], [9, 2], [9, 2], [9, 1], [9, 1], [9, 2]],
    performAction: function anonymous(yytext, yyleng, yylineno, yy, yystate, $$, _$) {
      var $0 = $$.length - 1;
      switch (yystate) {
        case 1:
          return $$[$0 - 1];
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
    },
    table: [{ 3: 1, 4: [1, 2] }, { 1: [3] }, o($V0, [2, 2], { 5: 3 }), { 6: [1, 4], 7: 5, 8: [1, 6], 9: 7, 10: [1, 8], 11: $V1, 12: $V2, 14: $V3, 16: $V4, 17: $V5, 18: $V6 }, o($V0, [2, 7], { 1: [2, 1] }), o($V0, [2, 3]), { 9: 15, 11: $V1, 12: $V2, 14: $V3, 16: $V4, 17: $V5, 18: $V6 }, o($V0, [2, 5]), o($V0, [2, 6]), o($V0, [2, 8]), { 13: [1, 16] }, { 15: [1, 17] }, o($V0, [2, 11]), o($V0, [2, 12]), { 19: [1, 18] }, o($V0, [2, 4]), o($V0, [2, 9]), o($V0, [2, 10]), o($V0, [2, 13])],
    defaultActions: {},
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
            break;
          case 1:
            break;
          case 2:
            return 10;
          case 3:
            break;
          case 4:
            break;
          case 5:
            return 4;
          case 6:
            return 11;
          case 7:
            this.begin("acc_title");
            return 12;
          case 8:
            this.popState();
            return "acc_title_value";
          case 9:
            this.begin("acc_descr");
            return 14;
          case 10:
            this.popState();
            return "acc_descr_value";
          case 11:
            this.begin("acc_descr_multiline");
            break;
          case 12:
            this.popState();
            break;
          case 13:
            return "acc_descr_multiline_value";
          case 14:
            return 17;
          case 15:
            return 18;
          case 16:
            return 19;
          case 17:
            return ":";
          case 18:
            return 6;
          case 19:
            return "INVALID";
        }
      },
      rules: [/^(?:%(?!\{)[^\n]*)/i, /^(?:[^\}]%%[^\n]*)/i, /^(?:[\n]+)/i, /^(?:\s+)/i, /^(?:#[^\n]*)/i, /^(?:journey\b)/i, /^(?:title\s[^#\n;]+)/i, /^(?:accTitle\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*\{\s*)/i, /^(?:[\}])/i, /^(?:[^\}]*)/i, /^(?:section\s[^#:\n;]+)/i, /^(?:[^#:\n;]+)/i, /^(?::[^#\n;]+)/i, /^(?::)/i, /^(?:$)/i, /^(?:.)/i],
      conditions: { "acc_descr_multiline": { "rules": [12, 13], "inclusive": false }, "acc_descr": { "rules": [10], "inclusive": false }, "acc_title": { "rules": [8], "inclusive": false }, "INITIAL": { "rules": [0, 1, 2, 3, 4, 5, 6, 7, 9, 11, 14, 15, 16, 17, 18, 19], "inclusive": true } }
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
let currentSection = "";
const sections = [];
const tasks = [];
const rawTasks = [];
const clear = function() {
  sections.length = 0;
  tasks.length = 0;
  currentSection = "";
  rawTasks.length = 0;
  (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.t)();
};
const addSection = function(txt) {
  currentSection = txt;
  sections.push(txt);
};
const getSections = function() {
  return sections;
};
const getTasks = function() {
  let allItemsProcessed = compileTasks();
  const maxDepth = 100;
  let iterationCount = 0;
  while (!allItemsProcessed && iterationCount < maxDepth) {
    allItemsProcessed = compileTasks();
    iterationCount++;
  }
  tasks.push(...rawTasks);
  return tasks;
};
const updateActors = function() {
  const tempActors = [];
  tasks.forEach((task) => {
    if (task.people) {
      tempActors.push(...task.people);
    }
  });
  const unique = new Set(tempActors);
  return [...unique].sort();
};
const addTask = function(descr, taskData) {
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
};
const addTaskOrg = function(descr) {
  const newTask = {
    section: currentSection,
    type: currentSection,
    description: descr,
    task: descr,
    classes: []
  };
  tasks.push(newTask);
};
const compileTasks = function() {
  const compileTask = function(pos) {
    return rawTasks[pos].processed;
  };
  let allProcessed = true;
  for (const [i, rawTask] of rawTasks.entries()) {
    compileTask(i);
    allProcessed = allProcessed && rawTask.processed;
  }
  return allProcessed;
};
const getActors = function() {
  return updateActors();
};
const db = {
  getConfig: () => (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().journey,
  clear,
  setDiagramTitle: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.q,
  getDiagramTitle: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.r,
  setAccTitle: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.s,
  getAccTitle: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.g,
  setAccDescription: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.b,
  getAccDescription: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.a,
  addSection,
  getSections,
  getTasks,
  addTask,
  addTaskOrg,
  getActors
};
const getStyles = (options) => `.label {
    font-family: 'trebuchet ms', verdana, arial, sans-serif;
    font-family: var(--mermaid-font-family);
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
    font-family: 'trebuchet ms', verdana, arial, sans-serif;
    font-family: var(--mermaid-font-family);
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
`;
const styles = getStyles;
const drawRect = function(elem, rectData) {
  return (0,_svgDrawCommon_92a8ff2b_js__WEBPACK_IMPORTED_MODULE_5__.d)(elem, rectData);
};
const drawFace = function(element, faceData) {
  const radius = 15;
  const circleElement = element.append("circle").attr("cx", faceData.cx).attr("cy", faceData.cy).attr("class", "face").attr("r", radius).attr("stroke-width", 2).attr("overflow", "visible");
  const face = element.append("g");
  face.append("circle").attr("cx", faceData.cx - radius / 3).attr("cy", faceData.cy - radius / 3).attr("r", 1.5).attr("stroke-width", 2).attr("fill", "#666").attr("stroke", "#666");
  face.append("circle").attr("cx", faceData.cx + radius / 3).attr("cy", faceData.cy - radius / 3).attr("r", 1.5).attr("stroke-width", 2).attr("fill", "#666").attr("stroke", "#666");
  function smile(face2) {
    const arc$1 = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .arc */ .Nb1)().startAngle(Math.PI / 2).endAngle(3 * (Math.PI / 2)).innerRadius(radius / 2).outerRadius(radius / 2.2);
    face2.append("path").attr("class", "mouth").attr("d", arc$1).attr("transform", "translate(" + faceData.cx + "," + (faceData.cy + 2) + ")");
  }
  function sad(face2) {
    const arc$1 = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .arc */ .Nb1)().startAngle(3 * Math.PI / 2).endAngle(5 * (Math.PI / 2)).innerRadius(radius / 2).outerRadius(radius / 2.2);
    face2.append("path").attr("class", "mouth").attr("d", arc$1).attr("transform", "translate(" + faceData.cx + "," + (faceData.cy + 7) + ")");
  }
  function ambivalent(face2) {
    face2.append("line").attr("class", "mouth").attr("stroke", 2).attr("x1", faceData.cx - 5).attr("y1", faceData.cy + 7).attr("x2", faceData.cx + 5).attr("y2", faceData.cy + 7).attr("class", "mouth").attr("stroke-width", "1px").attr("stroke", "#666");
  }
  if (faceData.score > 3) {
    smile(face);
  } else if (faceData.score < 3) {
    sad(face);
  } else {
    ambivalent(face);
  }
  return circleElement;
};
const drawCircle = function(element, circleData) {
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
};
const drawText = function(elem, textData) {
  return (0,_svgDrawCommon_92a8ff2b_js__WEBPACK_IMPORTED_MODULE_5__.f)(elem, textData);
};
const drawLabel = function(elem, txtObject) {
  function genPoints(x, y, width, height, cut) {
    return x + "," + y + " " + (x + width) + "," + y + " " + (x + width) + "," + (y + height - cut) + " " + (x + width - cut * 1.2) + "," + (y + height) + " " + x + "," + (y + height);
  }
  const polygon = elem.append("polygon");
  polygon.attr("points", genPoints(txtObject.x, txtObject.y, 50, 20, 7));
  polygon.attr("class", "labelBox");
  txtObject.y = txtObject.y + txtObject.labelMargin;
  txtObject.x = txtObject.x + 0.5 * txtObject.labelMargin;
  drawText(elem, txtObject);
};
const drawSection = function(elem, section, conf2) {
  const g = elem.append("g");
  const rect = (0,_svgDrawCommon_92a8ff2b_js__WEBPACK_IMPORTED_MODULE_5__.g)();
  rect.x = section.x;
  rect.y = section.y;
  rect.fill = section.fill;
  rect.width = conf2.width * section.taskCount + // width of the tasks
  conf2.diagramMarginX * (section.taskCount - 1);
  rect.height = conf2.height;
  rect.class = "journey-section section-type-" + section.num;
  rect.rx = 3;
  rect.ry = 3;
  drawRect(g, rect);
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
};
let taskCount = -1;
const drawTask = function(elem, task, conf2) {
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
  const rect = (0,_svgDrawCommon_92a8ff2b_js__WEBPACK_IMPORTED_MODULE_5__.g)();
  rect.x = task.x;
  rect.y = task.y;
  rect.fill = task.fill;
  rect.width = conf2.width;
  rect.height = conf2.height;
  rect.class = "task task-type-" + task.num;
  rect.rx = 3;
  rect.ry = 3;
  drawRect(g, rect);
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
};
const drawBackgroundRect = function(elem, bounds2) {
  (0,_svgDrawCommon_92a8ff2b_js__WEBPACK_IMPORTED_MODULE_5__.a)(elem, bounds2);
};
const _drawTextCandidateFunc = function() {
  function byText(content, g, x, y, width, height, textAttrs, colour) {
    const text = g.append("text").attr("x", x + width / 2).attr("y", y + height / 2 + 5).style("font-color", colour).style("text-anchor", "middle").text(content);
    _setTextAttrs(text, textAttrs);
  }
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
  function byFo(content, g, x, y, width, height, textAttrs, conf2) {
    const body = g.append("switch");
    const f = body.append("foreignObject").attr("x", x).attr("y", y).attr("width", width).attr("height", height).attr("position", "fixed");
    const text = f.append("xhtml:div").style("display", "table").style("height", "100%").style("width", "100%");
    text.append("div").attr("class", "label").style("display", "table-cell").style("text-align", "center").style("vertical-align", "middle").text(content);
    byTspan(content, body, x, y, width, height, textAttrs, conf2);
    _setTextAttrs(text, textAttrs);
  }
  function _setTextAttrs(toText, fromTextAttrsDict) {
    for (const key in fromTextAttrsDict) {
      if (key in fromTextAttrsDict) {
        toText.attr(key, fromTextAttrsDict[key]);
      }
    }
  }
  return function(conf2) {
    return conf2.textPlacement === "fo" ? byFo : conf2.textPlacement === "old" ? byText : byTspan;
  };
}();
const initGraphics = function(graphics) {
  graphics.append("defs").append("marker").attr("id", "arrowhead").attr("refX", 5).attr("refY", 2).attr("markerWidth", 6).attr("markerHeight", 4).attr("orient", "auto").append("path").attr("d", "M 0,0 V 4 L6,2 Z");
};
const svgDraw = {
  drawRect,
  drawCircle,
  drawSection,
  drawText,
  drawLabel,
  drawTask,
  drawBackgroundRect,
  initGraphics
};
const setConf = function(cnf) {
  const keys = Object.keys(cnf);
  keys.forEach(function(key) {
    conf[key] = cnf[key];
  });
};
const actors = {};
function drawActorLegend(diagram2) {
  const conf2 = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().journey;
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
    svgDraw.drawCircle(diagram2, circleData);
    const labelData = {
      x: 40,
      y: yPos + 7,
      fill: "#666",
      text: person,
      textMargin: conf2.boxTextMargin | 5
    };
    svgDraw.drawText(diagram2, labelData);
    yPos += 20;
  });
}
const conf = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().journey;
const LEFT_MARGIN = conf.leftMargin;
const draw = function(text, id, version, diagObj) {
  const conf2 = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().journey;
  const securityLevel = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().securityLevel;
  let sandboxElement;
  if (securityLevel === "sandbox") {
    sandboxElement = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)("#i" + id);
  }
  const root = securityLevel === "sandbox" ? (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(sandboxElement.nodes()[0].contentDocument.body) : (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)("body");
  bounds.init();
  const diagram2 = root.select("#" + id);
  svgDraw.initGraphics(diagram2);
  const tasks2 = diagObj.db.getTasks();
  const title = diagObj.db.getDiagramTitle();
  const actorNames = diagObj.db.getActors();
  for (const member in actors) {
    delete actors[member];
  }
  let actorPos = 0;
  actorNames.forEach((actorName) => {
    actors[actorName] = {
      color: conf2.actorColours[actorPos % conf2.actorColours.length],
      position: actorPos
    };
    actorPos++;
  });
  drawActorLegend(diagram2);
  bounds.insert(0, 0, LEFT_MARGIN, Object.keys(actors).length * 50);
  drawTasks(diagram2, tasks2, 0);
  const box = bounds.getBounds();
  if (title) {
    diagram2.append("text").text(title).attr("x", LEFT_MARGIN).attr("font-size", "4ex").attr("font-weight", "bold").attr("y", 25);
  }
  const height = box.stopy - box.starty + 2 * conf2.diagramMarginY;
  const width = LEFT_MARGIN + box.stopx + 2 * conf2.diagramMarginX;
  (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.i)(diagram2, height, width, conf2.useMaxWidth);
  diagram2.append("line").attr("x1", LEFT_MARGIN).attr("y1", conf2.height * 4).attr("x2", width - LEFT_MARGIN - 4).attr("y2", conf2.height * 4).attr("stroke-width", 4).attr("stroke", "black").attr("marker-end", "url(#arrowhead)");
  const extraVertForTitle = title ? 70 : 0;
  diagram2.attr("viewBox", `${box.startx} -25 ${width} ${height + extraVertForTitle}`);
  diagram2.attr("preserveAspectRatio", "xMinYMin meet");
  diagram2.attr("height", height + extraVertForTitle + 25);
};
const bounds = {
  data: {
    startx: void 0,
    stopx: void 0,
    starty: void 0,
    stopy: void 0
  },
  verticalPos: 0,
  sequenceItems: [],
  init: function() {
    this.sequenceItems = [];
    this.data = {
      startx: void 0,
      stopx: void 0,
      starty: void 0,
      stopy: void 0
    };
    this.verticalPos = 0;
  },
  updateVal: function(obj, key, val, fun) {
    if (obj[key] === void 0) {
      obj[key] = val;
    } else {
      obj[key] = fun(val, obj[key]);
    }
  },
  updateBounds: function(startx, starty, stopx, stopy) {
    const conf2 = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().journey;
    const _self = this;
    let cnt = 0;
    function updateFn(type) {
      return function updateItemBounds(item) {
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
      };
    }
    this.sequenceItems.forEach(updateFn());
  },
  insert: function(startx, starty, stopx, stopy) {
    const _startx = Math.min(startx, stopx);
    const _stopx = Math.max(startx, stopx);
    const _starty = Math.min(starty, stopy);
    const _stopy = Math.max(starty, stopy);
    this.updateVal(bounds.data, "startx", _startx, Math.min);
    this.updateVal(bounds.data, "starty", _starty, Math.min);
    this.updateVal(bounds.data, "stopx", _stopx, Math.max);
    this.updateVal(bounds.data, "stopy", _stopy, Math.max);
    this.updateBounds(_startx, _starty, _stopx, _stopy);
  },
  bumpVerticalPos: function(bump) {
    this.verticalPos = this.verticalPos + bump;
    this.data.stopy = this.verticalPos;
  },
  getVerticalPos: function() {
    return this.verticalPos;
  },
  getBounds: function() {
    return this.data;
  }
};
const fills = conf.sectionFills;
const textColours = conf.sectionColours;
const drawTasks = function(diagram2, tasks2, verticalPos) {
  const conf2 = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().journey;
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
        x: i * conf2.taskMargin + i * conf2.width + LEFT_MARGIN,
        y: 50,
        text: task.section,
        fill,
        num,
        colour,
        taskCount: taskInSectionCount
      };
      svgDraw.drawSection(diagram2, section, conf2);
      lastSection = task.section;
      sectionNumber++;
    }
    const taskActors = task.people.reduce((acc, actorName) => {
      if (actors[actorName]) {
        acc[actorName] = actors[actorName];
      }
      return acc;
    }, {});
    task.x = i * conf2.taskMargin + i * conf2.width + LEFT_MARGIN;
    task.y = taskPos;
    task.width = conf2.diagramMarginX;
    task.height = conf2.diagramMarginY;
    task.colour = colour;
    task.fill = fill;
    task.num = num;
    task.actors = taskActors;
    svgDraw.drawTask(diagram2, task, conf2);
    bounds.insert(task.x, task.y, task.x + task.width + conf2.taskMargin, 300 + 5 * 30);
  }
};
const renderer = {
  setConf,
  draw
};
const diagram = {
  parser: parser$1,
  db,
  renderer,
  styles,
  init: (cnf) => {
    renderer.setConf(cnf.journey);
    db.clear();
  }
};



/***/ }),

/***/ 26476:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   a: () => (/* binding */ drawBackgroundRect),
/* harmony export */   b: () => (/* binding */ drawEmbeddedImage),
/* harmony export */   c: () => (/* binding */ drawImage),
/* harmony export */   d: () => (/* binding */ drawRect),
/* harmony export */   e: () => (/* binding */ getTextObj),
/* harmony export */   f: () => (/* binding */ drawText),
/* harmony export */   g: () => (/* binding */ getNoteRect)
/* harmony export */ });
/* harmony import */ var _braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(7608);
/* harmony import */ var _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(24028);


const drawRect = (element, rectData) => {
  const rectElement = element.append("rect");
  rectElement.attr("x", rectData.x);
  rectElement.attr("y", rectData.y);
  rectElement.attr("fill", rectData.fill);
  rectElement.attr("stroke", rectData.stroke);
  rectElement.attr("width", rectData.width);
  rectElement.attr("height", rectData.height);
  rectData.rx !== void 0 && rectElement.attr("rx", rectData.rx);
  rectData.ry !== void 0 && rectElement.attr("ry", rectData.ry);
  if (rectData.attrs !== void 0) {
    for (const attrKey in rectData.attrs) {
      rectElement.attr(attrKey, rectData.attrs[attrKey]);
    }
  }
  rectData.class !== void 0 && rectElement.attr("class", rectData.class);
  return rectElement;
};
const drawBackgroundRect = (element, bounds) => {
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
};
const drawText = (element, textData) => {
  const nText = textData.text.replace(_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_1__.H, " ");
  const textElem = element.append("text");
  textElem.attr("x", textData.x);
  textElem.attr("y", textData.y);
  textElem.attr("class", "legend");
  textElem.style("text-anchor", textData.anchor);
  textData.class !== void 0 && textElem.attr("class", textData.class);
  const tspan = textElem.append("tspan");
  tspan.attr("x", textData.x + textData.textMargin * 2);
  tspan.text(nText);
  return textElem;
};
const drawImage = (elem, x, y, link) => {
  const imageElement = elem.append("image");
  imageElement.attr("x", x);
  imageElement.attr("y", y);
  const sanitizedLink = (0,_braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_0__/* .sanitizeUrl */ .Nm)(link);
  imageElement.attr("xlink:href", sanitizedLink);
};
const drawEmbeddedImage = (element, x, y, link) => {
  const imageElement = element.append("use");
  imageElement.attr("x", x);
  imageElement.attr("y", y);
  const sanitizedLink = (0,_braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_0__/* .sanitizeUrl */ .Nm)(link);
  imageElement.attr("xlink:href", `#${sanitizedLink}`);
};
const getNoteRect = () => {
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
};
const getTextObj = () => {
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
};



/***/ })

}]);
//# sourceMappingURL=4035.f4b13866b60b1af40230.js.map?v=f4b13866b60b1af40230
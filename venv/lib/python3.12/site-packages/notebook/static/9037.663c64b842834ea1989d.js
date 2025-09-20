"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[9037],{

/***/ 59037:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   factor: () => (/* binding */ factor)
/* harmony export */ });
/* harmony import */ var _simple_mode_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(4306);


const factor = (0,_simple_mode_js__WEBPACK_IMPORTED_MODULE_0__/* .simpleMode */ .Q)({
    start: [
      // comments
      {regex: /#?!.*/, token: "comment"},
      // strings """, multiline --> state
      {regex: /"""/, token: "string", next: "string3"},
      {regex: /(STRING:)(\s)/, token: ["keyword", null], next: "string2"},
      {regex: /\S*?"/, token: "string", next: "string"},
      // numbers: dec, hex, unicode, bin, fractional, complex
      {regex: /(?:0x[\d,a-f]+)|(?:0o[0-7]+)|(?:0b[0,1]+)|(?:\-?\d+.?\d*)(?=\s)/, token: "number"},
      //{regex: /[+-]?/} //fractional
      // definition: defining word, defined word, etc
      {regex: /((?:GENERIC)|\:?\:)(\s+)(\S+)(\s+)(\()/, token: ["keyword", null, "def", null, "bracket"], next: "stack"},
      // method definition: defining word, type, defined word, etc
      {regex: /(M\:)(\s+)(\S+)(\s+)(\S+)/, token: ["keyword", null, "def", null, "tag"]},
      // vocabulary using --> state
      {regex: /USING\:/, token: "keyword", next: "vocabulary"},
      // vocabulary definition/use
      {regex: /(USE\:|IN\:)(\s+)(\S+)(?=\s|$)/, token: ["keyword", null, "tag"]},
      // definition: a defining word, defined word
      {regex: /(\S+\:)(\s+)(\S+)(?=\s|$)/, token: ["keyword", null, "def"]},
      // "keywords", incl. ; t f . [ ] { } defining words
      {regex: /(?:;|\\|t|f|if|loop|while|until|do|PRIVATE>|<PRIVATE|\.|\S*\[|\]|\S*\{|\})(?=\s|$)/, token: "keyword"},
      // <constructors> and the like
      {regex: /\S+[\)>\.\*\?]+(?=\s|$)/, token: "builtin"},
      {regex: /[\)><]+\S+(?=\s|$)/, token: "builtin"},
      // operators
      {regex: /(?:[\+\-\=\/\*<>])(?=\s|$)/, token: "keyword"},
      // any id (?)
      {regex: /\S+/, token: "variable"},
      {regex: /\s+|./, token: null}
    ],
    vocabulary: [
      {regex: /;/, token: "keyword", next: "start"},
      {regex: /\S+/, token: "tag"},
      {regex: /\s+|./, token: null}
    ],
    string: [
      {regex: /(?:[^\\]|\\.)*?"/, token: "string", next: "start"},
      {regex: /.*/, token: "string"}
    ],
    string2: [
      {regex: /^;/, token: "keyword", next: "start"},
      {regex: /.*/, token: "string"}
    ],
    string3: [
      {regex: /(?:[^\\]|\\.)*?"""/, token: "string", next: "start"},
      {regex: /.*/, token: "string"}
    ],
    stack: [
      {regex: /\)/, token: "bracket", next: "start"},
      {regex: /--/, token: "bracket"},
      {regex: /\S+/, token: "meta"},
      {regex: /\s+|./, token: null}
    ],
    languageData: {
      name: "factor",
      dontIndentStates: ["start", "vocabulary", "string", "string3", "stack"],
      commentTokens: {line: "!"}
    }
  });


/***/ }),

/***/ 4306:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Q: () => (/* binding */ simpleMode)
/* harmony export */ });
function simpleMode(states) {
  ensureState(states, "start");
  var states_ = {}, meta = states.languageData || {}, hasIndentation = false;
  for (var state in states) if (state != meta && states.hasOwnProperty(state)) {
    var list = states_[state] = [], orig = states[state];
    for (var i = 0; i < orig.length; i++) {
      var data = orig[i];
      list.push(new Rule(data, states));
      if (data.indent || data.dedent) hasIndentation = true;
    }
  }
  return {
    name: meta.name,
    startState: function() {
      return {state: "start", pending: null, indent: hasIndentation ? [] : null};
    },
    copyState: function(state) {
      var s = {state: state.state, pending: state.pending, indent: state.indent && state.indent.slice(0)};
      if (state.stack)
        s.stack = state.stack.slice(0);
      return s;
    },
    token: tokenFunction(states_),
    indent: indentFunction(states_, meta),
    mergeTokens: meta.mergeTokens,
    languageData: meta
  }
};

function ensureState(states, name) {
  if (!states.hasOwnProperty(name))
    throw new Error("Undefined state " + name + " in simple mode");
}

function toRegex(val, caret) {
  if (!val) return /(?:)/;
  var flags = "";
  if (val instanceof RegExp) {
    if (val.ignoreCase) flags = "i";
    val = val.source;
  } else {
    val = String(val);
  }
  return new RegExp((caret === false ? "" : "^") + "(?:" + val + ")", flags);
}

function asToken(val) {
  if (!val) return null;
  if (val.apply) return val
  if (typeof val == "string") return val.replace(/\./g, " ");
  var result = [];
  for (var i = 0; i < val.length; i++)
    result.push(val[i] && val[i].replace(/\./g, " "));
  return result;
}

function Rule(data, states) {
  if (data.next || data.push) ensureState(states, data.next || data.push);
  this.regex = toRegex(data.regex);
  this.token = asToken(data.token);
  this.data = data;
}

function tokenFunction(states) {
  return function(stream, state) {
    if (state.pending) {
      var pend = state.pending.shift();
      if (state.pending.length == 0) state.pending = null;
      stream.pos += pend.text.length;
      return pend.token;
    }

    var curState = states[state.state];
    for (var i = 0; i < curState.length; i++) {
      var rule = curState[i];
      var matches = (!rule.data.sol || stream.sol()) && stream.match(rule.regex);
      if (matches) {
        if (rule.data.next) {
          state.state = rule.data.next;
        } else if (rule.data.push) {
          (state.stack || (state.stack = [])).push(state.state);
          state.state = rule.data.push;
        } else if (rule.data.pop && state.stack && state.stack.length) {
          state.state = state.stack.pop();
        }

        if (rule.data.indent)
          state.indent.push(stream.indentation() + stream.indentUnit);
        if (rule.data.dedent)
          state.indent.pop();
        var token = rule.token
        if (token && token.apply) token = token(matches)
        if (matches.length > 2 && rule.token && typeof rule.token != "string") {
          state.pending = [];
          for (var j = 2; j < matches.length; j++)
            if (matches[j])
              state.pending.push({text: matches[j], token: rule.token[j - 1]});
          stream.backUp(matches[0].length - (matches[1] ? matches[1].length : 0));
          return token[0];
        } else if (token && token.join) {
          return token[0];
        } else {
          return token;
        }
      }
    }
    stream.next();
    return null;
  };
}

function indentFunction(states, meta) {
  return function(state, textAfter) {
    if (state.indent == null || meta.dontIndentStates && meta.dontIndentStates.indexOf(state.state) > -1)
      return null

    var pos = state.indent.length - 1, rules = states[state.state];
    scan: for (;;) {
      for (var i = 0; i < rules.length; i++) {
        var rule = rules[i];
        if (rule.data.dedent && rule.data.dedentIfLineStart !== false) {
          var m = rule.regex.exec(textAfter);
          if (m && m[0]) {
            pos--;
            if (rule.next || rule.push) rules = states[rule.next || rule.push];
            textAfter = textAfter.slice(m[0].length);
            continue scan;
          }
        }
      }
      break;
    }
    return pos < 0 ? 0 : state.indent[pos];
  };
}


/***/ })

}]);
//# sourceMappingURL=9037.663c64b842834ea1989d.js.map?v=663c64b842834ea1989d
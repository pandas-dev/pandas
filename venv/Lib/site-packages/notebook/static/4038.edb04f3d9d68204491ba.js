"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[4038],{

/***/ 44038:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   yacas: () => (/* binding */ yacas)
/* harmony export */ });
function words(str) {
  var obj = {}, words = str.split(" ");
  for (var i = 0; i < words.length; ++i) obj[words[i]] = true;
  return obj;
}

var bodiedOps = words("Assert BackQuote D Defun Deriv For ForEach FromFile " +
                      "FromString Function Integrate InverseTaylor Limit " +
                      "LocalSymbols Macro MacroRule MacroRulePattern " +
                      "NIntegrate Rule RulePattern Subst TD TExplicitSum " +
                      "TSum Taylor Taylor1 Taylor2 Taylor3 ToFile " +
                      "ToStdout ToString TraceRule Until While");

// patterns
var pFloatForm  = "(?:(?:\\.\\d+|\\d+\\.\\d*|\\d+)(?:[eE][+-]?\\d+)?)";
var pIdentifier = "(?:[a-zA-Z\\$'][a-zA-Z0-9\\$']*)";

// regular expressions
var reFloatForm    = new RegExp(pFloatForm);
var reIdentifier   = new RegExp(pIdentifier);
var rePattern      = new RegExp(pIdentifier + "?_" + pIdentifier);
var reFunctionLike = new RegExp(pIdentifier + "\\s*\\(");

function tokenBase(stream, state) {
  var ch;

  // get next character
  ch = stream.next();

  // string
  if (ch === '"') {
    state.tokenize = tokenString;
    return state.tokenize(stream, state);
  }

  // comment
  if (ch === '/') {
    if (stream.eat('*')) {
      state.tokenize = tokenComment;
      return state.tokenize(stream, state);
    }
    if (stream.eat("/")) {
      stream.skipToEnd();
      return "comment";
    }
  }

  // go back one character
  stream.backUp(1);

  // update scope info
  var m = stream.match(/^(\w+)\s*\(/, false);
  if (m !== null && bodiedOps.hasOwnProperty(m[1]))
    state.scopes.push('bodied');

  var scope = currentScope(state);

  if (scope === 'bodied' && ch === '[')
    state.scopes.pop();

  if (ch === '[' || ch === '{' || ch === '(')
    state.scopes.push(ch);

  scope = currentScope(state);

  if (scope === '[' && ch === ']' ||
      scope === '{' && ch === '}' ||
      scope === '(' && ch === ')')
    state.scopes.pop();

  if (ch === ';') {
    while (scope === 'bodied') {
      state.scopes.pop();
      scope = currentScope(state);
    }
  }

  // look for ordered rules
  if (stream.match(/\d+ *#/, true, false)) {
    return 'qualifier';
  }

  // look for numbers
  if (stream.match(reFloatForm, true, false)) {
    return 'number';
  }

  // look for placeholders
  if (stream.match(rePattern, true, false)) {
    return 'variableName.special';
  }

  // match all braces separately
  if (stream.match(/(?:\[|\]|{|}|\(|\))/, true, false)) {
    return 'bracket';
  }

  // literals looking like function calls
  if (stream.match(reFunctionLike, true, false)) {
    stream.backUp(1);
    return 'variableName.function';
  }

  // all other identifiers
  if (stream.match(reIdentifier, true, false)) {
    return 'variable';
  }

  // operators; note that operators like @@ or /; are matched separately for each symbol.
  if (stream.match(/(?:\\|\+|\-|\*|\/|,|;|\.|:|@|~|=|>|<|&|\||_|`|'|\^|\?|!|%|#)/, true, false)) {
    return 'operator';
  }

  // everything else is an error
  return 'error';
}

function tokenString(stream, state) {
  var next, end = false, escaped = false;
  while ((next = stream.next()) != null) {
    if (next === '"' && !escaped) {
      end = true;
      break;
    }
    escaped = !escaped && next === '\\';
  }
  if (end && !escaped) {
    state.tokenize = tokenBase;
  }
  return 'string';
};

function tokenComment(stream, state) {
  var prev, next;
  while((next = stream.next()) != null) {
    if (prev === '*' && next === '/') {
      state.tokenize = tokenBase;
      break;
    }
    prev = next;
  }
  return 'comment';
}

function currentScope(state) {
  var scope = null;
  if (state.scopes.length > 0)
    scope = state.scopes[state.scopes.length - 1];
  return scope;
}

const yacas = {
  name: "yacas",
  startState: function() {
    return {
      tokenize: tokenBase,
      scopes: []
    };
  },
  token: function(stream, state) {
    if (stream.eatSpace()) return null;
    return state.tokenize(stream, state);
  },
  indent: function(state, textAfter, cx) {
    if (state.tokenize !== tokenBase && state.tokenize !== null)
      return null;

    var delta = 0;
    if (textAfter === ']' || textAfter === '];' ||
        textAfter === '}' || textAfter === '};' ||
        textAfter === ');')
      delta = -1;

    return (state.scopes.length + delta) * cx.unit;
  },

  languageData: {
    electricInput: /[{}\[\]()\;]/,
    commentTokens: {line: "//", block: {open: "/*", close: "*/"}}
  }
};


/***/ })

}]);
//# sourceMappingURL=4038.edb04f3d9d68204491ba.js.map?v=edb04f3d9d68204491ba
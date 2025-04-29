"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[1837],{

/***/ 81837:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   julia: () => (/* binding */ julia)
/* harmony export */ });
function wordRegexp(words, end, pre) {
  if (typeof pre === "undefined") pre = "";
  if (typeof end === "undefined") { end = "\\b"; }
  return new RegExp("^" + pre + "((" + words.join(")|(") + "))" + end);
}

var octChar = "\\\\[0-7]{1,3}";
var hexChar = "\\\\x[A-Fa-f0-9]{1,2}";
var sChar = "\\\\[abefnrtv0%?'\"\\\\]";
var uChar = "([^\\u0027\\u005C\\uD800-\\uDFFF]|[\\uD800-\\uDFFF][\\uDC00-\\uDFFF])";

var asciiOperatorsList = [
  "[<>]:", "[<>=]=", "<<=?", ">>>?=?", "=>", "--?>", "<--[->]?", "\\/\\/",
  "\\.{2,3}", "[\\.\\\\%*+\\-<>!\\/^|&]=?", "\\?", "\\$", "~", ":"
];
var operators = wordRegexp([
  "[<>]:", "[<>=]=", "[!=]==", "<<=?", ">>>?=?", "=>?", "--?>", "<--[->]?", "\\/\\/",
  "[\\\\%*+\\-<>!\\/^|&\\u00F7\\u22BB]=?", "\\?", "\\$", "~", ":",
  "\\u00D7", "\\u2208", "\\u2209", "\\u220B", "\\u220C", "\\u2218",
  "\\u221A", "\\u221B", "\\u2229", "\\u222A", "\\u2260", "\\u2264",
  "\\u2265", "\\u2286", "\\u2288", "\\u228A", "\\u22C5",
  "\\b(in|isa)\\b(?!\.?\\()"
], "");
var delimiters = /^[;,()[\]{}]/;
var identifiers = /^[_A-Za-z\u00A1-\u2217\u2219-\uFFFF][\w\u00A1-\u2217\u2219-\uFFFF]*!*/;

var chars = wordRegexp([octChar, hexChar, sChar, uChar], "'");

var openersList = ["begin", "function", "type", "struct", "immutable", "let",
                   "macro", "for", "while", "quote", "if", "else", "elseif", "try",
                   "finally", "catch", "do"];

var closersList = ["end", "else", "elseif", "catch", "finally"];

var keywordsList = ["if", "else", "elseif", "while", "for", "begin", "let",
                    "end", "do", "try", "catch", "finally", "return", "break", "continue",
                    "global", "local", "const", "export", "import", "importall", "using",
                    "function", "where", "macro", "module", "baremodule", "struct", "type",
                    "mutable", "immutable", "quote", "typealias", "abstract", "primitive",
                    "bitstype"];

var builtinsList = ["true", "false", "nothing", "NaN", "Inf"];

var openers = wordRegexp(openersList);
var closers = wordRegexp(closersList);
var keywords = wordRegexp(keywordsList);
var builtins = wordRegexp(builtinsList);

var macro = /^@[_A-Za-z\u00A1-\uFFFF][\w\u00A1-\uFFFF]*!*/;
var symbol = /^:[_A-Za-z\u00A1-\uFFFF][\w\u00A1-\uFFFF]*!*/;
var stringPrefixes = /^(`|([_A-Za-z\u00A1-\uFFFF]*"("")?))/;

var macroOperators = wordRegexp(asciiOperatorsList, "", "@");
var symbolOperators = wordRegexp(asciiOperatorsList, "", ":");

function inArray(state) {
  return (state.nestedArrays > 0);
}

function inGenerator(state) {
  return (state.nestedGenerators > 0);
}

function currentScope(state, n) {
  if (typeof(n) === "undefined") { n = 0; }
  if (state.scopes.length <= n) {
    return null;
  }
  return state.scopes[state.scopes.length - (n + 1)];
}

// tokenizers
function tokenBase(stream, state) {
  // Handle multiline comments
  if (stream.match('#=', false)) {
    state.tokenize = tokenComment;
    return state.tokenize(stream, state);
  }

  // Handle scope changes
  var leavingExpr = state.leavingExpr;
  if (stream.sol()) {
    leavingExpr = false;
  }
  state.leavingExpr = false;

  if (leavingExpr) {
    if (stream.match(/^'+/)) {
      return "operator";
    }
  }

  if (stream.match(/\.{4,}/)) {
    return "error";
  } else if (stream.match(/\.{1,3}/)) {
    return "operator";
  }

  if (stream.eatSpace()) {
    return null;
  }

  var ch = stream.peek();

  // Handle single line comments
  if (ch === '#') {
    stream.skipToEnd();
    return "comment";
  }

  if (ch === '[') {
    state.scopes.push('[');
    state.nestedArrays++;
  }

  if (ch === '(') {
    state.scopes.push('(');
    state.nestedGenerators++;
  }

  if (inArray(state) && ch === ']') {
    while (state.scopes.length && currentScope(state) !== "[") { state.scopes.pop(); }
    state.scopes.pop();
    state.nestedArrays--;
    state.leavingExpr = true;
  }

  if (inGenerator(state) && ch === ')') {
    while (state.scopes.length && currentScope(state) !== "(") { state.scopes.pop(); }
    state.scopes.pop();
    state.nestedGenerators--;
    state.leavingExpr = true;
  }

  if (inArray(state)) {
    if (state.lastToken == "end" && stream.match(':')) {
      return "operator";
    }
    if (stream.match('end')) {
      return "number";
    }
  }

  var match;
  if (match = stream.match(openers, false)) {
    state.scopes.push(match[0]);
  }

  if (stream.match(closers, false)) {
    state.scopes.pop();
  }

  // Handle type annotations
  if (stream.match(/^::(?![:\$])/)) {
    state.tokenize = tokenAnnotation;
    return state.tokenize(stream, state);
  }

  // Handle symbols
  if (!leavingExpr && (stream.match(symbol) || stream.match(symbolOperators))) {
    return "builtin";
  }

  // Handle parametric types
  //if (stream.match(/^{[^}]*}(?=\()/)) {
  //  return "builtin";
  //}

  // Handle operators and Delimiters
  if (stream.match(operators)) {
    return "operator";
  }

  // Handle Number Literals
  if (stream.match(/^\.?\d/, false)) {
    var imMatcher = RegExp(/^im\b/);
    var numberLiteral = false;
    if (stream.match(/^0x\.[0-9a-f_]+p[\+\-]?[_\d]+/i)) { numberLiteral = true; }
    // Integers
    if (stream.match(/^0x[0-9a-f_]+/i)) { numberLiteral = true; } // Hex
    if (stream.match(/^0b[01_]+/i)) { numberLiteral = true; } // Binary
    if (stream.match(/^0o[0-7_]+/i)) { numberLiteral = true; } // Octal
    // Floats
    if (stream.match(/^(?:(?:\d[_\d]*)?\.(?!\.)(?:\d[_\d]*)?|\d[_\d]*\.(?!\.)(?:\d[_\d]*))?([Eef][\+\-]?[_\d]+)?/i)) { numberLiteral = true; }
    if (stream.match(/^\d[_\d]*(e[\+\-]?\d+)?/i)) { numberLiteral = true; } // Decimal
    if (numberLiteral) {
      // Integer literals may be "long"
      stream.match(imMatcher);
      state.leavingExpr = true;
      return "number";
    }
  }

  // Handle Chars
  if (stream.match("'")) {
    state.tokenize = tokenChar;
    return state.tokenize(stream, state);
  }

  // Handle Strings
  if (stream.match(stringPrefixes)) {
    state.tokenize = tokenStringFactory(stream.current());
    return state.tokenize(stream, state);
  }

  if (stream.match(macro) || stream.match(macroOperators)) {
    return "meta";
  }

  if (stream.match(delimiters)) {
    return null;
  }

  if (stream.match(keywords)) {
    return "keyword";
  }

  if (stream.match(builtins)) {
    return "builtin";
  }

  var isDefinition = state.isDefinition || state.lastToken == "function" ||
      state.lastToken == "macro" || state.lastToken == "type" ||
      state.lastToken == "struct" || state.lastToken == "immutable";

  if (stream.match(identifiers)) {
    if (isDefinition) {
      if (stream.peek() === '.') {
        state.isDefinition = true;
        return "variable";
      }
      state.isDefinition = false;
      return "def";
    }
    state.leavingExpr = true;
    return "variable";
  }

  // Handle non-detected items
  stream.next();
  return "error";
}

function tokenAnnotation(stream, state) {
  stream.match(/.*?(?=[,;{}()=\s]|$)/);
  if (stream.match('{')) {
    state.nestedParameters++;
  } else if (stream.match('}') && state.nestedParameters > 0) {
    state.nestedParameters--;
  }
  if (state.nestedParameters > 0) {
    stream.match(/.*?(?={|})/) || stream.next();
  } else if (state.nestedParameters == 0) {
    state.tokenize = tokenBase;
  }
  return "builtin";
}

function tokenComment(stream, state) {
  if (stream.match('#=')) {
    state.nestedComments++;
  }
  if (!stream.match(/.*?(?=(#=|=#))/)) {
    stream.skipToEnd();
  }
  if (stream.match('=#')) {
    state.nestedComments--;
    if (state.nestedComments == 0)
      state.tokenize = tokenBase;
  }
  return "comment";
}

function tokenChar(stream, state) {
  var isChar = false, match;
  if (stream.match(chars)) {
    isChar = true;
  } else if (match = stream.match(/\\u([a-f0-9]{1,4})(?=')/i)) {
    var value = parseInt(match[1], 16);
    if (value <= 55295 || value >= 57344) { // (U+0,U+D7FF), (U+E000,U+FFFF)
      isChar = true;
      stream.next();
    }
  } else if (match = stream.match(/\\U([A-Fa-f0-9]{5,8})(?=')/)) {
    var value = parseInt(match[1], 16);
    if (value <= 1114111) { // U+10FFFF
      isChar = true;
      stream.next();
    }
  }
  if (isChar) {
    state.leavingExpr = true;
    state.tokenize = tokenBase;
    return "string";
  }
  if (!stream.match(/^[^']+(?=')/)) { stream.skipToEnd(); }
  if (stream.match("'")) { state.tokenize = tokenBase; }
  return "error";
}

function tokenStringFactory(delimiter) {
  if (delimiter.substr(-3) === '"""') {
    delimiter = '"""';
  } else if (delimiter.substr(-1) === '"') {
    delimiter = '"';
  }
  function tokenString(stream, state) {
    if (stream.eat('\\')) {
      stream.next();
    } else if (stream.match(delimiter)) {
      state.tokenize = tokenBase;
      state.leavingExpr = true;
      return "string";
    } else {
      stream.eat(/[`"]/);
    }
    stream.eatWhile(/[^\\`"]/);
    return "string";
  }
  return tokenString;
}

const julia = {
  name: "julia",
  startState: function() {
    return {
      tokenize: tokenBase,
      scopes: [],
      lastToken: null,
      leavingExpr: false,
      isDefinition: false,
      nestedArrays: 0,
      nestedComments: 0,
      nestedGenerators: 0,
      nestedParameters: 0,
      firstParenPos: -1
    };
  },

  token: function(stream, state) {
    var style = state.tokenize(stream, state);
    var current = stream.current();

    if (current && style) {
      state.lastToken = current;
    }

    return style;
  },

  indent: function(state, textAfter, cx) {
    var delta = 0;
    if ( textAfter === ']' || textAfter === ')' || /^end\b/.test(textAfter) ||
         /^else/.test(textAfter) || /^catch\b/.test(textAfter) || /^elseif\b/.test(textAfter) ||
         /^finally/.test(textAfter) ) {
      delta = -1;
    }
    return (state.scopes.length + delta) * cx.unit;
  },

  languageData: {
    indentOnInput: /^\s*(end|else|catch|finally)\b$/,
    commentTokens: {line: "#", block: {open: "#=", close: "=#"}},
    closeBrackets: {brackets: ["(", "[", "{", '"']},
    autocomplete: keywordsList.concat(builtinsList)
  }
};


/***/ })

}]);
//# sourceMappingURL=1837.6bbfd9967be58e1325f1.js.map?v=6bbfd9967be58e1325f1
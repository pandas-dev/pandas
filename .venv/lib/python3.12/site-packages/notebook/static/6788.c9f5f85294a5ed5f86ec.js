"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[6788],{

/***/ 26788:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   modelica: () => (/* binding */ modelica)
/* harmony export */ });
function words(str) {
  var obj = {}, words = str.split(" ");
  for (var i=0; i<words.length; ++i)
    obj[words[i]] = true;
  return obj;
}

var keywords = words("algorithm and annotation assert block break class connect connector constant constrainedby der discrete each else elseif elsewhen encapsulated end enumeration equation expandable extends external false final flow for function if import impure in initial inner input loop model not operator or outer output package parameter partial protected public pure record redeclare replaceable return stream then true type when while within")
var builtin = words("abs acos actualStream asin atan atan2 cardinality ceil cos cosh delay div edge exp floor getInstanceName homotopy inStream integer log log10 mod pre reinit rem semiLinear sign sin sinh spatialDistribution sqrt tan tanh")
var atoms = words("Real Boolean Integer String")

var completions = [].concat(Object.keys(keywords), Object.keys(builtin), Object.keys(atoms))

var isSingleOperatorChar = /[;=\(:\),{}.*<>+\-\/^\[\]]/;
var isDoubleOperatorChar = /(:=|<=|>=|==|<>|\.\+|\.\-|\.\*|\.\/|\.\^)/;
var isDigit = /[0-9]/;
var isNonDigit = /[_a-zA-Z]/;

function tokenLineComment(stream, state) {
  stream.skipToEnd();
  state.tokenize = null;
  return "comment";
}

function tokenBlockComment(stream, state) {
  var maybeEnd = false, ch;
  while (ch = stream.next()) {
    if (maybeEnd && ch == "/") {
      state.tokenize = null;
      break;
    }
    maybeEnd = (ch == "*");
  }
  return "comment";
}

function tokenString(stream, state) {
  var escaped = false, ch;
  while ((ch = stream.next()) != null) {
    if (ch == '"' && !escaped) {
      state.tokenize = null;
      state.sol = false;
      break;
    }
    escaped = !escaped && ch == "\\";
  }

  return "string";
}

function tokenIdent(stream, state) {
  stream.eatWhile(isDigit);
  while (stream.eat(isDigit) || stream.eat(isNonDigit)) { }


  var cur = stream.current();

  if(state.sol && (cur == "package" || cur == "model" || cur == "when" || cur == "connector")) state.level++;
  else if(state.sol && cur == "end" && state.level > 0) state.level--;

  state.tokenize = null;
  state.sol = false;

  if (keywords.propertyIsEnumerable(cur)) return "keyword";
  else if (builtin.propertyIsEnumerable(cur)) return "builtin";
  else if (atoms.propertyIsEnumerable(cur)) return "atom";
  else return "variable";
}

function tokenQIdent(stream, state) {
  while (stream.eat(/[^']/)) { }

  state.tokenize = null;
  state.sol = false;

  if(stream.eat("'"))
    return "variable";
  else
    return "error";
}

function tokenUnsignedNumber(stream, state) {
  stream.eatWhile(isDigit);
  if (stream.eat('.')) {
    stream.eatWhile(isDigit);
  }
  if (stream.eat('e') || stream.eat('E')) {
    if (!stream.eat('-'))
      stream.eat('+');
    stream.eatWhile(isDigit);
  }

  state.tokenize = null;
  state.sol = false;
  return "number";
}

// Interface
const modelica = {
  name: "modelica",
  startState: function() {
    return {
      tokenize: null,
      level: 0,
      sol: true
    };
  },

  token: function(stream, state) {
    if(state.tokenize != null) {
      return state.tokenize(stream, state);
    }

    if(stream.sol()) {
      state.sol = true;
    }

    // WHITESPACE
    if(stream.eatSpace()) {
      state.tokenize = null;
      return null;
    }

    var ch = stream.next();

    // LINECOMMENT
    if(ch == '/' && stream.eat('/')) {
      state.tokenize = tokenLineComment;
    }
    // BLOCKCOMMENT
    else if(ch == '/' && stream.eat('*')) {
      state.tokenize = tokenBlockComment;
    }
    // TWO SYMBOL TOKENS
    else if(isDoubleOperatorChar.test(ch+stream.peek())) {
      stream.next();
      state.tokenize = null;
      return "operator";
    }
    // SINGLE SYMBOL TOKENS
    else if(isSingleOperatorChar.test(ch)) {
      state.tokenize = null;
      return "operator";
    }
    // IDENT
    else if(isNonDigit.test(ch)) {
      state.tokenize = tokenIdent;
    }
    // Q-IDENT
    else if(ch == "'" && stream.peek() && stream.peek() != "'") {
      state.tokenize = tokenQIdent;
    }
    // STRING
    else if(ch == '"') {
      state.tokenize = tokenString;
    }
    // UNSIGNED_NUMBER
    else if(isDigit.test(ch)) {
      state.tokenize = tokenUnsignedNumber;
    }
    // ERROR
    else {
      state.tokenize = null;
      return "error";
    }

    return state.tokenize(stream, state);
  },

  indent: function(state, textAfter, cx) {
    if (state.tokenize != null) return null;

    var level = state.level;
    if(/(algorithm)/.test(textAfter)) level--;
    if(/(equation)/.test(textAfter)) level--;
    if(/(initial algorithm)/.test(textAfter)) level--;
    if(/(initial equation)/.test(textAfter)) level--;
    if(/(end)/.test(textAfter)) level--;

    if(level > 0)
      return cx.unit*level;
    else
      return 0;
  },

  languageData: {
    commentTokens: {line: "//", block: {open: "/*", close: "*/"}},
    autocomplete: completions
  }
};


/***/ })

}]);
//# sourceMappingURL=6788.c9f5f85294a5ed5f86ec.js.map?v=c9f5f85294a5ed5f86ec
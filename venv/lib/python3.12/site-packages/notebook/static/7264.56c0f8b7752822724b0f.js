"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[7264],{

/***/ 57264:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   eiffel: () => (/* binding */ eiffel)
/* harmony export */ });
function wordObj(words) {
  var o = {};
  for (var i = 0, e = words.length; i < e; ++i) o[words[i]] = true;
  return o;
}
var keywords = wordObj([
  'note',
  'across',
  'when',
  'variant',
  'until',
  'unique',
  'undefine',
  'then',
  'strip',
  'select',
  'retry',
  'rescue',
  'require',
  'rename',
  'reference',
  'redefine',
  'prefix',
  'once',
  'old',
  'obsolete',
  'loop',
  'local',
  'like',
  'is',
  'inspect',
  'infix',
  'include',
  'if',
  'frozen',
  'from',
  'external',
  'export',
  'ensure',
  'end',
  'elseif',
  'else',
  'do',
  'creation',
  'create',
  'check',
  'alias',
  'agent',
  'separate',
  'invariant',
  'inherit',
  'indexing',
  'feature',
  'expanded',
  'deferred',
  'class',
  'Void',
  'True',
  'Result',
  'Precursor',
  'False',
  'Current',
  'create',
  'attached',
  'detachable',
  'as',
  'and',
  'implies',
  'not',
  'or'
]);
var operators = wordObj([":=", "and then","and", "or","<<",">>"]);

function chain(newtok, stream, state) {
  state.tokenize.push(newtok);
  return newtok(stream, state);
}

function tokenBase(stream, state) {
  if (stream.eatSpace()) return null;
  var ch = stream.next();
  if (ch == '"'||ch == "'") {
    return chain(readQuoted(ch, "string"), stream, state);
  } else if (ch == "-"&&stream.eat("-")) {
    stream.skipToEnd();
    return "comment";
  } else if (ch == ":"&&stream.eat("=")) {
    return "operator";
  } else if (/[0-9]/.test(ch)) {
    stream.eatWhile(/[xXbBCc0-9\.]/);
    stream.eat(/[\?\!]/);
    return "variable";
  } else if (/[a-zA-Z_0-9]/.test(ch)) {
    stream.eatWhile(/[a-zA-Z_0-9]/);
    stream.eat(/[\?\!]/);
    return "variable";
  } else if (/[=+\-\/*^%<>~]/.test(ch)) {
    stream.eatWhile(/[=+\-\/*^%<>~]/);
    return "operator";
  } else {
    return null;
  }
}

function readQuoted(quote, style,  unescaped) {
  return function(stream, state) {
    var escaped = false, ch;
    while ((ch = stream.next()) != null) {
      if (ch == quote && (unescaped || !escaped)) {
        state.tokenize.pop();
        break;
      }
      escaped = !escaped && ch == "%";
    }
    return style;
  };
}

const eiffel = {
  name: "eiffel",
  startState: function() {
    return {tokenize: [tokenBase]};
  },

  token: function(stream, state) {
    var style = state.tokenize[state.tokenize.length-1](stream, state);
    if (style == "variable") {
      var word = stream.current();
      style = keywords.propertyIsEnumerable(stream.current()) ? "keyword"
        : operators.propertyIsEnumerable(stream.current()) ? "operator"
        : /^[A-Z][A-Z_0-9]*$/g.test(word) ? "tag"
        : /^0[bB][0-1]+$/g.test(word) ? "number"
        : /^0[cC][0-7]+$/g.test(word) ? "number"
        : /^0[xX][a-fA-F0-9]+$/g.test(word) ? "number"
        : /^([0-9]+\.[0-9]*)|([0-9]*\.[0-9]+)$/g.test(word) ? "number"
        : /^[0-9]+$/g.test(word) ? "number"
        : "variable";
    }
    return style;
  },
  languageData: {
    commentTokens: {line: "--"}
  }
};



/***/ })

}]);
//# sourceMappingURL=7264.56c0f8b7752822724b0f.js.map?v=56c0f8b7752822724b0f
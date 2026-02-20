"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[8433],{

/***/ 28433:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   troff: () => (/* binding */ troff)
/* harmony export */ });
var words = {};

function tokenBase(stream) {
  if (stream.eatSpace()) return null;

  var sol = stream.sol();
  var ch = stream.next();

  if (ch === '\\') {
    if (stream.match('fB') || stream.match('fR') || stream.match('fI') ||
        stream.match('u')  || stream.match('d')  ||
        stream.match('%')  || stream.match('&')) {
      return 'string';
    }
    if (stream.match('m[')) {
      stream.skipTo(']');
      stream.next();
      return 'string';
    }
    if (stream.match('s+') || stream.match('s-')) {
      stream.eatWhile(/[\d-]/);
      return 'string';
    }
    if (stream.match('\(') || stream.match('*\(')) {
      stream.eatWhile(/[\w-]/);
      return 'string';
    }
    return 'string';
  }
  if (sol && (ch === '.' || ch === '\'')) {
    if (stream.eat('\\') && stream.eat('\"')) {
      stream.skipToEnd();
      return 'comment';
    }
  }
  if (sol && ch === '.') {
    if (stream.match('B ') || stream.match('I ') || stream.match('R ')) {
      return 'attribute';
    }
    if (stream.match('TH ') || stream.match('SH ') || stream.match('SS ') || stream.match('HP ')) {
      stream.skipToEnd();
      return 'quote';
    }
    if ((stream.match(/[A-Z]/) && stream.match(/[A-Z]/)) || (stream.match(/[a-z]/) && stream.match(/[a-z]/))) {
      return 'attribute';
    }
  }
  stream.eatWhile(/[\w-]/);
  var cur = stream.current();
  return words.hasOwnProperty(cur) ? words[cur] : null;
}

function tokenize(stream, state) {
  return (state.tokens[0] || tokenBase) (stream, state);
};

const troff = {
  name: "troff",
  startState: function() {return {tokens:[]};},
  token: function(stream, state) {
    return tokenize(stream, state);
  }
};


/***/ })

}]);
//# sourceMappingURL=8433.ed9247b868845dc191b2.js.map?v=ed9247b868845dc191b2
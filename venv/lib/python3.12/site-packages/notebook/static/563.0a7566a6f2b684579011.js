"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[563],{

/***/ 10563:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   apl: () => (/* binding */ apl)
/* harmony export */ });
var builtInFuncs = {
  "+": ["conjugate", "add"],
  "−": ["negate", "subtract"],
  "×": ["signOf", "multiply"],
  "÷": ["reciprocal", "divide"],
  "⌈": ["ceiling", "greaterOf"],
  "⌊": ["floor", "lesserOf"],
  "∣": ["absolute", "residue"],
  "⍳": ["indexGenerate", "indexOf"],
  "?": ["roll", "deal"],
  "⋆": ["exponentiate", "toThePowerOf"],
  "⍟": ["naturalLog", "logToTheBase"],
  "○": ["piTimes", "circularFuncs"],
  "!": ["factorial", "binomial"],
  "⌹": ["matrixInverse", "matrixDivide"],
  "<": [null, "lessThan"],
  "≤": [null, "lessThanOrEqual"],
  "=": [null, "equals"],
  ">": [null, "greaterThan"],
  "≥": [null, "greaterThanOrEqual"],
  "≠": [null, "notEqual"],
  "≡": ["depth", "match"],
  "≢": [null, "notMatch"],
  "∈": ["enlist", "membership"],
  "⍷": [null, "find"],
  "∪": ["unique", "union"],
  "∩": [null, "intersection"],
  "∼": ["not", "without"],
  "∨": [null, "or"],
  "∧": [null, "and"],
  "⍱": [null, "nor"],
  "⍲": [null, "nand"],
  "⍴": ["shapeOf", "reshape"],
  ",": ["ravel", "catenate"],
  "⍪": [null, "firstAxisCatenate"],
  "⌽": ["reverse", "rotate"],
  "⊖": ["axis1Reverse", "axis1Rotate"],
  "⍉": ["transpose", null],
  "↑": ["first", "take"],
  "↓": [null, "drop"],
  "⊂": ["enclose", "partitionWithAxis"],
  "⊃": ["diclose", "pick"],
  "⌷": [null, "index"],
  "⍋": ["gradeUp", null],
  "⍒": ["gradeDown", null],
  "⊤": ["encode", null],
  "⊥": ["decode", null],
  "⍕": ["format", "formatByExample"],
  "⍎": ["execute", null],
  "⊣": ["stop", "left"],
  "⊢": ["pass", "right"]
};

var isOperator = /[\.\/⌿⍀¨⍣]/;
var isNiladic = /⍬/;
var isFunction = /[\+−×÷⌈⌊∣⍳\?⋆⍟○!⌹<≤=>≥≠≡≢∈⍷∪∩∼∨∧⍱⍲⍴,⍪⌽⊖⍉↑↓⊂⊃⌷⍋⍒⊤⊥⍕⍎⊣⊢]/;
var isArrow = /←/;
var isComment = /[⍝#].*$/;

var stringEater = function(type) {
  var prev;
  prev = false;
  return function(c) {
    prev = c;
    if (c === type) {
      return prev === "\\";
    }
    return true;
  };
};

const apl = {
  name: "apl",
  startState: function() {
    return {
      prev: false,
      func: false,
      op: false,
      string: false,
      escape: false
    };
  },
  token: function(stream, state) {
    var ch;
    if (stream.eatSpace()) {
      return null;
    }
    ch = stream.next();
    if (ch === '"' || ch === "'") {
      stream.eatWhile(stringEater(ch));
      stream.next();
      state.prev = true;
      return "string";
    }
    if (/[\[{\(]/.test(ch)) {
      state.prev = false;
      return null;
    }
    if (/[\]}\)]/.test(ch)) {
      state.prev = true;
      return null;
    }
    if (isNiladic.test(ch)) {
      state.prev = false;
      return "atom";
    }
    if (/[¯\d]/.test(ch)) {
      if (state.func) {
        state.func = false;
        state.prev = false;
      } else {
        state.prev = true;
      }
      stream.eatWhile(/[\w\.]/);
      return "number";
    }
    if (isOperator.test(ch)) {
      return "operator"
    }
    if (isArrow.test(ch)) {
      return "operator";
    }
    if (isFunction.test(ch)) {
      state.func = true;
      state.prev = false;
      return builtInFuncs[ch] ? "variableName.function.standard" : "variableName.function"
    }
    if (isComment.test(ch)) {
      stream.skipToEnd();
      return "comment";
    }
    if (ch === "∘" && stream.peek() === ".") {
      stream.next();
      return "variableName.function";
    }
    stream.eatWhile(/[\w\$_]/);
    state.prev = true;
    return "keyword";
  }
}


/***/ })

}]);
//# sourceMappingURL=563.0a7566a6f2b684579011.js.map?v=0a7566a6f2b684579011
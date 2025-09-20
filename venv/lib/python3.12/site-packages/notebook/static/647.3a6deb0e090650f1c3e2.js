"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[647],{

/***/ 90647:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   ebnf: () => (/* binding */ ebnf)
/* harmony export */ });
var commentType = {slash: 0, parenthesis: 1};
var stateType = {comment: 0, _string: 1, characterClass: 2};

const ebnf = {
  name: "ebnf",
  startState: function () {
    return {
      stringType: null,
      commentType: null,
      braced: 0,
      lhs: true,
      localState: null,
      stack: [],
      inDefinition: false
    };
  },
  token: function (stream, state) {
    if (!stream) return;

    //check for state changes
    if (state.stack.length === 0) {
      //strings
      if ((stream.peek() == '"') || (stream.peek() == "'")) {
        state.stringType = stream.peek();
        stream.next(); // Skip quote
        state.stack.unshift(stateType._string);
      } else if (stream.match('/*')) { //comments starting with /*
        state.stack.unshift(stateType.comment);
        state.commentType = commentType.slash;
      } else if (stream.match('(*')) { //comments starting with (*
        state.stack.unshift(stateType.comment);
        state.commentType = commentType.parenthesis;
      }
    }

    //return state
    //stack has
    switch (state.stack[0]) {
    case stateType._string:
      while (state.stack[0] === stateType._string && !stream.eol()) {
        if (stream.peek() === state.stringType) {
          stream.next(); // Skip quote
          state.stack.shift(); // Clear flag
        } else if (stream.peek() === "\\") {
          stream.next();
          stream.next();
        } else {
          stream.match(/^.[^\\\"\']*/);
        }
      }
      return state.lhs ? "property" : "string"; // Token style

    case stateType.comment:
      while (state.stack[0] === stateType.comment && !stream.eol()) {
        if (state.commentType === commentType.slash && stream.match('*/')) {
          state.stack.shift(); // Clear flag
          state.commentType = null;
        } else if (state.commentType === commentType.parenthesis && stream.match('*)')) {
          state.stack.shift(); // Clear flag
          state.commentType = null;
        } else {
          stream.match(/^.[^\*]*/);
        }
      }
      return "comment";

    case stateType.characterClass:
      while (state.stack[0] === stateType.characterClass && !stream.eol()) {
        if (!(stream.match(/^[^\]\\]+/) || stream.match('.'))) {
          state.stack.shift();
        }
      }
      return "operator";
    }

    var peek = stream.peek();

    //no stack
    switch (peek) {
    case "[":
      stream.next();
      state.stack.unshift(stateType.characterClass);
      return "bracket";
    case ":":
    case "|":
    case ";":
      stream.next();
      return "operator";
    case "%":
      if (stream.match("%%")) {
        return "header";
      } else if (stream.match(/[%][A-Za-z]+/)) {
        return "keyword";
      } else if (stream.match(/[%][}]/)) {
        return "bracket";
      }
      break;
    case "/":
      if (stream.match(/[\/][A-Za-z]+/)) {
        return "keyword";
      }
    case "\\":
      if (stream.match(/[\][a-z]+/)) {
        return "string.special";
      }
    case ".":
      if (stream.match(".")) {
        return "atom";
      }
    case "*":
    case "-":
    case "+":
    case "^":
      if (stream.match(peek)) {
        return "atom";
      }
    case "$":
      if (stream.match("$$")) {
        return "builtin";
      } else if (stream.match(/[$][0-9]+/)) {
        return "variableName.special";
      }
    case "<":
      if (stream.match(/<<[a-zA-Z_]+>>/)) {
        return "builtin";
      }
    }

    if (stream.match('//')) {
      stream.skipToEnd();
      return "comment";
    } else if (stream.match('return')) {
      return "operator";
    } else if (stream.match(/^[a-zA-Z_][a-zA-Z0-9_]*/)) {
      if (stream.match(/(?=[\(.])/)) {
        return "variable";
      } else if (stream.match(/(?=[\s\n]*[:=])/)) {
        return "def";
      }
      return "variableName.special";
    } else if (["[", "]", "(", ")"].indexOf(stream.peek()) != -1) {
      stream.next();
      return "bracket";
    } else if (!stream.eatSpace()) {
      stream.next();
    }
    return null;
  }
};


/***/ })

}]);
//# sourceMappingURL=647.3a6deb0e090650f1c3e2.js.map?v=3a6deb0e090650f1c3e2
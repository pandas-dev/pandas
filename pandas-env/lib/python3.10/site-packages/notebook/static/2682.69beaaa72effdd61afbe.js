"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[2682],{

/***/ 92682:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   toml: () => (/* binding */ toml)
/* harmony export */ });
const toml = {
  name: "toml",
  startState: function () {
    return {
      inString: false,
      stringType: "",
      lhs: true,
      inArray: 0
    };
  },
  token: function (stream, state) {
    //check for state changes
    if (!state.inString && ((stream.peek() == '"') || (stream.peek() == "'"))) {
      state.stringType = stream.peek();
      stream.next(); // Skip quote
      state.inString = true; // Update state
    }
    if (stream.sol() && state.inArray === 0) {
      state.lhs = true;
    }
    //return state
    if (state.inString) {
      while (state.inString && !stream.eol()) {
        if (stream.peek() === state.stringType) {
          stream.next(); // Skip quote
          state.inString = false; // Clear flag
        } else if (stream.peek() === '\\') {
          stream.next();
          stream.next();
        } else {
          stream.match(/^.[^\\\"\']*/);
        }
      }
      return state.lhs ? "property" : "string"; // Token style
    } else if (state.inArray && stream.peek() === ']') {
      stream.next();
      state.inArray--;
      return 'bracket';
    } else if (state.lhs && stream.peek() === '[' && stream.skipTo(']')) {
      stream.next();//skip closing ]
      // array of objects has an extra open & close []
      if (stream.peek() === ']') stream.next();
      return "atom";
    } else if (stream.peek() === "#") {
      stream.skipToEnd();
      return "comment";
    } else if (stream.eatSpace()) {
      return null;
    } else if (state.lhs && stream.eatWhile(function (c) { return c != '=' && c != ' '; })) {
      return "property";
    } else if (state.lhs && stream.peek() === "=") {
      stream.next();
      state.lhs = false;
      return null;
    } else if (!state.lhs && stream.match(/^\d\d\d\d[\d\-\:\.T]*Z/)) {
      return 'atom'; //date
    } else if (!state.lhs && (stream.match('true') || stream.match('false'))) {
      return 'atom';
    } else if (!state.lhs && stream.peek() === '[') {
      state.inArray++;
      stream.next();
      return 'bracket';
    } else if (!state.lhs && stream.match(/^\-?\d+(?:\.\d+)?/)) {
      return 'number';
    } else if (!stream.eatSpace()) {
      stream.next();
    }
    return null;
  },
  languageData: {
    commentTokens: { line: '#' },
  },
};


/***/ })

}]);
//# sourceMappingURL=2682.69beaaa72effdd61afbe.js.map?v=69beaaa72effdd61afbe
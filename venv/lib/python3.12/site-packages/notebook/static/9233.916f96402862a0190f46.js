"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[9233],{

/***/ 9233:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   yaml: () => (/* binding */ yaml)
/* harmony export */ });
var cons = ['true', 'false', 'on', 'off', 'yes', 'no'];
var keywordRegex = new RegExp("\\b(("+cons.join(")|(")+"))$", 'i');

const yaml = {
  name: "yaml",
  token: function(stream, state) {
    var ch = stream.peek();
    var esc = state.escaped;
    state.escaped = false;
    /* comments */
    if (ch == "#" && (stream.pos == 0 || /\s/.test(stream.string.charAt(stream.pos - 1)))) {
      stream.skipToEnd();
      return "comment";
    }

    if (stream.match(/^('([^']|\\.)*'?|"([^"]|\\.)*"?)/))
      return "string";

    if (state.literal && stream.indentation() > state.keyCol) {
      stream.skipToEnd(); return "string";
    } else if (state.literal) { state.literal = false; }
    if (stream.sol()) {
      state.keyCol = 0;
      state.pair = false;
      state.pairStart = false;
      /* document start */
      if(stream.match('---')) { return "def"; }
      /* document end */
      if (stream.match('...')) { return "def"; }
      /* array list item */
      if (stream.match(/^\s*-\s+/)) { return 'meta'; }
    }
    /* inline pairs/lists */
    if (stream.match(/^(\{|\}|\[|\])/)) {
      if (ch == '{')
        state.inlinePairs++;
      else if (ch == '}')
        state.inlinePairs--;
      else if (ch == '[')
        state.inlineList++;
      else
        state.inlineList--;
      return 'meta';
    }

    /* list separator */
    if (state.inlineList > 0 && !esc && ch == ',') {
      stream.next();
      return 'meta';
    }
    /* pairs separator */
    if (state.inlinePairs > 0 && !esc && ch == ',') {
      state.keyCol = 0;
      state.pair = false;
      state.pairStart = false;
      stream.next();
      return 'meta';
    }

    /* start of value of a pair */
    if (state.pairStart) {
      /* block literals */
      if (stream.match(/^\s*(\||\>)\s*/)) { state.literal = true; return 'meta'; };
      /* references */
      if (stream.match(/^\s*(\&|\*)[a-z0-9\._-]+\b/i)) { return 'variable'; }
      /* numbers */
      if (state.inlinePairs == 0 && stream.match(/^\s*-?[0-9\.\,]+\s?$/)) { return 'number'; }
      if (state.inlinePairs > 0 && stream.match(/^\s*-?[0-9\.\,]+\s?(?=(,|}))/)) { return 'number'; }
      /* keywords */
      if (stream.match(keywordRegex)) { return 'keyword'; }
    }

    /* pairs (associative arrays) -> key */
    if (!state.pair && stream.match(/^\s*(?:[,\[\]{}&*!|>'"%@`][^\s'":]|[^,\[\]{}#&*!|>'"%@`])[^#]*?(?=\s*:($|\s))/)) {
      state.pair = true;
      state.keyCol = stream.indentation();
      return "atom";
    }
    if (state.pair && stream.match(/^:\s*/)) { state.pairStart = true; return 'meta'; }

    /* nothing found, continue */
    state.pairStart = false;
    state.escaped = (ch == '\\');
    stream.next();
    return null;
  },
  startState: function() {
    return {
      pair: false,
      pairStart: false,
      keyCol: 0,
      inlinePairs: 0,
      inlineList: 0,
      literal: false,
      escaped: false
    };
  },
  languageData: {
    commentTokens: {line: "#"}
  }
};


/***/ })

}]);
//# sourceMappingURL=9233.916f96402862a0190f46.js.map?v=916f96402862a0190f46
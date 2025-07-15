"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[7803],{

/***/ 87803:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   properties: () => (/* binding */ properties)
/* harmony export */ });
const properties = {
  name: "properties",

  token: function(stream, state) {
    var sol = stream.sol() || state.afterSection;
    var eol = stream.eol();

    state.afterSection = false;

    if (sol) {
      if (state.nextMultiline) {
        state.inMultiline = true;
        state.nextMultiline = false;
      } else {
        state.position = "def";
      }
    }

    if (eol && ! state.nextMultiline) {
      state.inMultiline = false;
      state.position = "def";
    }

    if (sol) {
      while(stream.eatSpace()) {}
    }

    var ch = stream.next();

    if (sol && (ch === "#" || ch === "!" || ch === ";")) {
      state.position = "comment";
      stream.skipToEnd();
      return "comment";
    } else if (sol && ch === "[") {
      state.afterSection = true;
      stream.skipTo("]"); stream.eat("]");
      return "header";
    } else if (ch === "=" || ch === ":") {
      state.position = "quote";
      return null;
    } else if (ch === "\\" && state.position === "quote") {
      if (stream.eol()) {  // end of line?
        // Multiline value
        state.nextMultiline = true;
      }
    }

    return state.position;
  },

  startState: function() {
    return {
      position : "def",       // Current position, "def", "quote" or "comment"
      nextMultiline : false,  // Is the next line multiline value
      inMultiline : false,    // Is the current line a multiline value
      afterSection : false    // Did we just open a section
    };
  }

};


/***/ })

}]);
//# sourceMappingURL=7803.0c44e7b8d148353eed87.js.map?v=0c44e7b8d148353eed87
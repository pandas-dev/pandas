"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[6139],{

/***/ 66139:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   mumps: () => (/* binding */ mumps)
/* harmony export */ });
function wordRegexp(words) {
  return new RegExp("^((" + words.join(")|(") + "))\\b", "i");
}

var singleOperators = new RegExp("^[\\+\\-\\*/&#!_?\\\\<>=\\'\\[\\]]");
var doubleOperators = new RegExp("^(('=)|(<=)|(>=)|('>)|('<)|([[)|(]])|(^$))");
var singleDelimiters = new RegExp("^[\\.,:]");
var brackets = new RegExp("[()]");
var identifiers = new RegExp("^[%A-Za-z][A-Za-z0-9]*");
var commandKeywords = ["break","close","do","else","for","goto", "halt", "hang", "if", "job","kill","lock","merge","new","open", "quit", "read", "set", "tcommit", "trollback", "tstart", "use", "view", "write", "xecute", "b","c","d","e","f","g", "h", "i", "j","k","l","m","n","o", "q", "r", "s", "tc", "tro", "ts", "u", "v", "w", "x"];
// The following list includes intrinsic functions _and_ special variables
var intrinsicFuncsWords = ["\\$ascii", "\\$char", "\\$data", "\\$ecode", "\\$estack", "\\$etrap", "\\$extract", "\\$find", "\\$fnumber", "\\$get", "\\$horolog", "\\$io", "\\$increment", "\\$job", "\\$justify", "\\$length", "\\$name", "\\$next", "\\$order", "\\$piece", "\\$qlength", "\\$qsubscript", "\\$query", "\\$quit", "\\$random", "\\$reverse", "\\$select", "\\$stack", "\\$test", "\\$text", "\\$translate", "\\$view", "\\$x", "\\$y", "\\$a", "\\$c", "\\$d", "\\$e", "\\$ec", "\\$es", "\\$et", "\\$f", "\\$fn", "\\$g", "\\$h", "\\$i", "\\$j", "\\$l", "\\$n", "\\$na", "\\$o", "\\$p", "\\$q", "\\$ql", "\\$qs", "\\$r", "\\$re", "\\$s", "\\$st", "\\$t", "\\$tr", "\\$v", "\\$z"];
var intrinsicFuncs = wordRegexp(intrinsicFuncsWords);
var command = wordRegexp(commandKeywords);

function tokenBase(stream, state) {
  if (stream.sol()) {
    state.label = true;
    state.commandMode = 0;
  }

  // The <space> character has meaning in MUMPS. Ignoring consecutive
  // spaces would interfere with interpreting whether the next non-space
  // character belongs to the command or argument context.

  // Examine each character and update a mode variable whose interpretation is:
  //   >0 => command    0 => argument    <0 => command post-conditional
  var ch = stream.peek();

  if (ch == " " || ch == "\t") { // Pre-process <space>
    state.label = false;
    if (state.commandMode == 0)
      state.commandMode = 1;
    else if ((state.commandMode < 0) || (state.commandMode == 2))
      state.commandMode = 0;
  } else if ((ch != ".") && (state.commandMode > 0)) {
    if (ch == ":")
      state.commandMode = -1;   // SIS - Command post-conditional
    else
      state.commandMode = 2;
  }

  // Do not color parameter list as line tag
  if ((ch === "(") || (ch === "\u0009"))
    state.label = false;

  // MUMPS comment starts with ";"
  if (ch === ";") {
    stream.skipToEnd();
    return "comment";
  }

  // Number Literals // SIS/RLM - MUMPS permits canonic number followed by concatenate operator
  if (stream.match(/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?/))
    return "number";

  // Handle Strings
  if (ch == '"') {
    if (stream.skipTo('"')) {
      stream.next();
      return "string";
    } else {
      stream.skipToEnd();
      return "error";
    }
  }

  // Handle operators and Delimiters
  if (stream.match(doubleOperators) || stream.match(singleOperators))
    return "operator";

  // Prevents leading "." in DO block from falling through to error
  if (stream.match(singleDelimiters))
    return null;

  if (brackets.test(ch)) {
    stream.next();
    return "bracket";
  }

  if (state.commandMode > 0 && stream.match(command))
    return "controlKeyword";

  if (stream.match(intrinsicFuncs))
    return "builtin";

  if (stream.match(identifiers))
    return "variable";

  // Detect dollar-sign when not a documented intrinsic function
  // "^" may introduce a GVN or SSVN - Color same as function
  if (ch === "$" || ch === "^") {
    stream.next();
    return "builtin";
  }

  // MUMPS Indirection
  if (ch === "@") {
    stream.next();
    return "string.special";
  }

  if (/[\w%]/.test(ch)) {
    stream.eatWhile(/[\w%]/);
    return "variable";
  }

  // Handle non-detected items
  stream.next();
  return "error";
}

const mumps = {
  name: "mumps",
  startState: function() {
    return {
      label: false,
      commandMode: 0
    };
  },

  token: function(stream, state) {
    var style = tokenBase(stream, state);
    if (state.label) return "tag";
    return style;
  }
};


/***/ })

}]);
//# sourceMappingURL=6139.9b4118bd8223a51fa897.js.map?v=9b4118bd8223a51fa897
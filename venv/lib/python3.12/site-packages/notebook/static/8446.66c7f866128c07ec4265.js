"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[8446],{

/***/ 8446:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   tcl: () => (/* binding */ tcl)
/* harmony export */ });
function parseWords(str) {
  var obj = {}, words = str.split(" ");
  for (var i = 0; i < words.length; ++i) obj[words[i]] = true;
  return obj;
}
var keywords = parseWords("Tcl safe after append array auto_execok auto_import auto_load " +
                          "auto_mkindex auto_mkindex_old auto_qualify auto_reset bgerror " +
                          "binary break catch cd close concat continue dde eof encoding error " +
                          "eval exec exit expr fblocked fconfigure fcopy file fileevent filename " +
                          "filename flush for foreach format gets glob global history http if " +
                          "incr info interp join lappend lindex linsert list llength load lrange " +
                          "lreplace lsearch lset lsort memory msgcat namespace open package parray " +
                          "pid pkg::create pkg_mkIndex proc puts pwd re_syntax read regex regexp " +
                          "registry regsub rename resource return scan seek set socket source split " +
                          "string subst switch tcl_endOfWord tcl_findLibrary tcl_startOfNextWord " +
                          "tcl_wordBreakAfter tcl_startOfPreviousWord tcl_wordBreakBefore tcltest " +
                          "tclvars tell time trace unknown unset update uplevel upvar variable " +
                          "vwait");
var functions = parseWords("if elseif else and not or eq ne in ni for foreach while switch");
var isOperatorChar = /[+\-*&%=<>!?^\/\|]/;
function chain(stream, state, f) {
  state.tokenize = f;
  return f(stream, state);
}
function tokenBase(stream, state) {
  var beforeParams = state.beforeParams;
  state.beforeParams = false;
  var ch = stream.next();
  if ((ch == '"' || ch == "'") && state.inParams) {
    return chain(stream, state, tokenString(ch));
  } else if (/[\[\]{}\(\),;\.]/.test(ch)) {
    if (ch == "(" && beforeParams) state.inParams = true;
    else if (ch == ")") state.inParams = false;
    return null;
  } else if (/\d/.test(ch)) {
    stream.eatWhile(/[\w\.]/);
    return "number";
  } else if (ch == "#") {
    if (stream.eat("*"))
      return chain(stream, state, tokenComment);
    if (ch == "#" && stream.match(/ *\[ *\[/))
      return chain(stream, state, tokenUnparsed);
    stream.skipToEnd();
    return "comment";
  } else if (ch == '"') {
    stream.skipTo(/"/);
    return "comment";
  } else if (ch == "$") {
    stream.eatWhile(/[$_a-z0-9A-Z\.{:]/);
    stream.eatWhile(/}/);
    state.beforeParams = true;
    return "builtin";
  } else if (isOperatorChar.test(ch)) {
    stream.eatWhile(isOperatorChar);
    return "comment";
  } else {
    stream.eatWhile(/[\w\$_{}\xa1-\uffff]/);
    var word = stream.current().toLowerCase();
    if (keywords && keywords.propertyIsEnumerable(word))
      return "keyword";
    if (functions && functions.propertyIsEnumerable(word)) {
      state.beforeParams = true;
      return "keyword";
    }
    return null;
  }
}
function tokenString(quote) {
  return function(stream, state) {
    var escaped = false, next, end = false;
    while ((next = stream.next()) != null) {
      if (next == quote && !escaped) {
        end = true;
        break;
      }
      escaped = !escaped && next == "\\";
    }
    if (end) state.tokenize = tokenBase;
    return "string";
  };
}
function tokenComment(stream, state) {
  var maybeEnd = false, ch;
  while (ch = stream.next()) {
    if (ch == "#" && maybeEnd) {
      state.tokenize = tokenBase;
      break;
    }
    maybeEnd = (ch == "*");
  }
  return "comment";
}
function tokenUnparsed(stream, state) {
  var maybeEnd = 0, ch;
  while (ch = stream.next()) {
    if (ch == "#" && maybeEnd == 2) {
      state.tokenize = tokenBase;
      break;
    }
    if (ch == "]")
      maybeEnd++;
    else if (ch != " ")
      maybeEnd = 0;
  }
  return "meta";
}
const tcl = {
  name: "tcl",
  startState: function() {
    return {
      tokenize: tokenBase,
      beforeParams: false,
      inParams: false
    };
  },
  token: function(stream, state) {
    if (stream.eatSpace()) return null;
    return state.tokenize(stream, state);
  },
  languageData: {
    commentTokens: {line: "#"}
  }
};


/***/ })

}]);
//# sourceMappingURL=8446.66c7f866128c07ec4265.js.map?v=66c7f866128c07ec4265
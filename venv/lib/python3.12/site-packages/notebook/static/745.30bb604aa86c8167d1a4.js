"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[745],{

/***/ 50745:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   sieve: () => (/* binding */ sieve)
/* harmony export */ });
function words(str) {
  var obj = {}, words = str.split(" ");
  for (var i = 0; i < words.length; ++i) obj[words[i]] = true;
  return obj;
}

var keywords = words("if elsif else stop require");
var atoms = words("true false not");

function tokenBase(stream, state) {

  var ch = stream.next();
  if (ch == "/" && stream.eat("*")) {
    state.tokenize = tokenCComment;
    return tokenCComment(stream, state);
  }

  if (ch === '#') {
    stream.skipToEnd();
    return "comment";
  }

  if (ch == "\"") {
    state.tokenize = tokenString(ch);
    return state.tokenize(stream, state);
  }

  if (ch == "(") {
    state._indent.push("(");
    // add virtual angel wings so that editor behaves...
    // ...more sane incase of broken brackets
    state._indent.push("{");
    return null;
  }

  if (ch === "{") {
    state._indent.push("{");
    return null;
  }

  if (ch == ")")  {
    state._indent.pop();
    state._indent.pop();
  }

  if (ch === "}") {
    state._indent.pop();
    return null;
  }

  if (ch == ",")
    return null;

  if (ch == ";")
    return null;


  if (/[{}\(\),;]/.test(ch))
    return null;

  // 1*DIGIT "K" / "M" / "G"
  if (/\d/.test(ch)) {
    stream.eatWhile(/[\d]/);
    stream.eat(/[KkMmGg]/);
    return "number";
  }

  // ":" (ALPHA / "_") *(ALPHA / DIGIT / "_")
  if (ch == ":") {
    stream.eatWhile(/[a-zA-Z_]/);
    stream.eatWhile(/[a-zA-Z0-9_]/);

    return "operator";
  }

  stream.eatWhile(/\w/);
  var cur = stream.current();

  // "text:" *(SP / HTAB) (hash-comment / CRLF)
  // *(multiline-literal / multiline-dotstart)
  // "." CRLF
  if ((cur == "text") && stream.eat(":"))
  {
    state.tokenize = tokenMultiLineString;
    return "string";
  }

  if (keywords.propertyIsEnumerable(cur))
    return "keyword";

  if (atoms.propertyIsEnumerable(cur))
    return "atom";

  return null;
}

function tokenMultiLineString(stream, state)
{
  state._multiLineString = true;
  // the first line is special it may contain a comment
  if (!stream.sol()) {
    stream.eatSpace();

    if (stream.peek() == "#") {
      stream.skipToEnd();
      return "comment";
    }

    stream.skipToEnd();
    return "string";
  }

  if ((stream.next() == ".")  && (stream.eol()))
  {
    state._multiLineString = false;
    state.tokenize = tokenBase;
  }

  return "string";
}

function tokenCComment(stream, state) {
  var maybeEnd = false, ch;
  while ((ch = stream.next()) != null) {
    if (maybeEnd && ch == "/") {
      state.tokenize = tokenBase;
      break;
    }
    maybeEnd = (ch == "*");
  }
  return "comment";
}

function tokenString(quote) {
  return function(stream, state) {
    var escaped = false, ch;
    while ((ch = stream.next()) != null) {
      if (ch == quote && !escaped)
        break;
      escaped = !escaped && ch == "\\";
    }
    if (!escaped) state.tokenize = tokenBase;
    return "string";
  };
}

const sieve = {
  name: "sieve",
  startState: function(base) {
    return {tokenize: tokenBase,
            baseIndent: base || 0,
            _indent: []};
  },

  token: function(stream, state) {
    if (stream.eatSpace())
      return null;

    return (state.tokenize || tokenBase)(stream, state);
  },

  indent: function(state, _textAfter, cx) {
    var length = state._indent.length;
    if (_textAfter && (_textAfter[0] == "}"))
      length--;

    if (length <0)
      length = 0;

    return length * cx.unit;
  },

  languageData: {
    indentOnInput: /^\s*\}$/
  }
};


/***/ })

}]);
//# sourceMappingURL=745.30bb604aa86c8167d1a4.js.map?v=30bb604aa86c8167d1a4
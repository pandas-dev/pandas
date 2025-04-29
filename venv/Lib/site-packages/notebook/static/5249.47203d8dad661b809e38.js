"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[5249],{

/***/ 55249:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   sparql: () => (/* binding */ sparql)
/* harmony export */ });
var curPunc;

function wordRegexp(words) {
  return new RegExp("^(?:" + words.join("|") + ")$", "i");
}
var ops = wordRegexp(["str", "lang", "langmatches", "datatype", "bound", "sameterm", "isiri", "isuri",
                      "iri", "uri", "bnode", "count", "sum", "min", "max", "avg", "sample",
                      "group_concat", "rand", "abs", "ceil", "floor", "round", "concat", "substr", "strlen",
                      "replace", "ucase", "lcase", "encode_for_uri", "contains", "strstarts", "strends",
                      "strbefore", "strafter", "year", "month", "day", "hours", "minutes", "seconds",
                      "timezone", "tz", "now", "uuid", "struuid", "md5", "sha1", "sha256", "sha384",
                      "sha512", "coalesce", "if", "strlang", "strdt", "isnumeric", "regex", "exists",
                      "isblank", "isliteral", "a", "bind"]);
var keywords = wordRegexp(["base", "prefix", "select", "distinct", "reduced", "construct", "describe",
                           "ask", "from", "named", "where", "order", "limit", "offset", "filter", "optional",
                           "graph", "by", "asc", "desc", "as", "having", "undef", "values", "group",
                           "minus", "in", "not", "service", "silent", "using", "insert", "delete", "union",
                           "true", "false", "with",
                           "data", "copy", "to", "move", "add", "create", "drop", "clear", "load", "into"]);
var operatorChars = /[*+\-<>=&|\^\/!\?]/;
var PN_CHARS = "[A-Za-z_\\-0-9]";
var PREFIX_START = new RegExp("[A-Za-z]");
var PREFIX_REMAINDER = new RegExp("((" + PN_CHARS + "|\\.)*(" + PN_CHARS + "))?:");

function tokenBase(stream, state) {
  var ch = stream.next();
  curPunc = null;
  if (ch == "$" || ch == "?") {
    if(ch == "?" && stream.match(/\s/, false)){
      return "operator";
    }
    stream.match(/^[A-Za-z0-9_\u00C0-\u00D6\u00D8-\u00F6\u00F8-\u02FF\u0370-\u037D\u037F-\u1FFF\u200C-\u200D\u2070-\u218F\u2C00-\u2FEF\u3001-\uD7FF\uF900-\uFDCF\uFDF0-\uFFFD][A-Za-z0-9_\u00B7\u00C0-\u00D6\u00D8-\u00F6\u00F8-\u037D\u037F-\u1FFF\u200C-\u200D\u203F-\u2040\u2070-\u218F\u2C00-\u2FEF\u3001-\uD7FF\uF900-\uFDCF\uFDF0-\uFFFD]*/);
    return "variableName.local";
  }
  else if (ch == "<" && !stream.match(/^[\s\u00a0=]/, false)) {
    stream.match(/^[^\s\u00a0>]*>?/);
    return "atom";
  }
  else if (ch == "\"" || ch == "'") {
    state.tokenize = tokenLiteral(ch);
    return state.tokenize(stream, state);
  }
  else if (/[{}\(\),\.;\[\]]/.test(ch)) {
    curPunc = ch;
    return "bracket";
  }
  else if (ch == "#") {
    stream.skipToEnd();
    return "comment";
  }
  else if (operatorChars.test(ch)) {
    return "operator";
  }
  else if (ch == ":") {
    eatPnLocal(stream);
    return "atom";
  }
  else if (ch == "@") {
    stream.eatWhile(/[a-z\d\-]/i);
    return "meta";
  }
  else if (PREFIX_START.test(ch) && stream.match(PREFIX_REMAINDER)) {
    eatPnLocal(stream);
    return "atom";
  }
  stream.eatWhile(/[_\w\d]/);
  var word = stream.current();
  if (ops.test(word))
    return "builtin";
  else if (keywords.test(word))
    return "keyword";
  else
    return "variable";
}

function eatPnLocal(stream) {
  stream.match(/(\.(?=[\w_\-\\%])|[:\w_-]|\\[-\\_~.!$&'()*+,;=/?#@%]|%[a-f\d][a-f\d])+/i);
}

function tokenLiteral(quote) {
  return function(stream, state) {
    var escaped = false, ch;
    while ((ch = stream.next()) != null) {
      if (ch == quote && !escaped) {
        state.tokenize = tokenBase;
        break;
      }
      escaped = !escaped && ch == "\\";
    }
    return "string";
  };
}

function pushContext(state, type, col) {
  state.context = {prev: state.context, indent: state.indent, col: col, type: type};
}
function popContext(state) {
  state.indent = state.context.indent;
  state.context = state.context.prev;
}

const sparql = {
  name: "sparql",

  startState: function() {
    return {tokenize: tokenBase,
            context: null,
            indent: 0,
            col: 0};
  },

  token: function(stream, state) {
    if (stream.sol()) {
      if (state.context && state.context.align == null) state.context.align = false;
      state.indent = stream.indentation();
    }
    if (stream.eatSpace()) return null;
    var style = state.tokenize(stream, state);

    if (style != "comment" && state.context && state.context.align == null && state.context.type != "pattern") {
      state.context.align = true;
    }

    if (curPunc == "(") pushContext(state, ")", stream.column());
    else if (curPunc == "[") pushContext(state, "]", stream.column());
    else if (curPunc == "{") pushContext(state, "}", stream.column());
    else if (/[\]\}\)]/.test(curPunc)) {
      while (state.context && state.context.type == "pattern") popContext(state);
      if (state.context && curPunc == state.context.type) {
        popContext(state);
        if (curPunc == "}" && state.context && state.context.type == "pattern")
          popContext(state);
      }
    }
    else if (curPunc == "." && state.context && state.context.type == "pattern") popContext(state);
    else if (/atom|string|variable/.test(style) && state.context) {
      if (/[\}\]]/.test(state.context.type))
        pushContext(state, "pattern", stream.column());
      else if (state.context.type == "pattern" && !state.context.align) {
        state.context.align = true;
        state.context.col = stream.column();
      }
    }

    return style;
  },

  indent: function(state, textAfter, cx) {
    var firstChar = textAfter && textAfter.charAt(0);
    var context = state.context;
    if (/[\]\}]/.test(firstChar))
      while (context && context.type == "pattern") context = context.prev;

    var closing = context && firstChar == context.type;
    if (!context)
      return 0;
    else if (context.type == "pattern")
      return context.col;
    else if (context.align)
      return context.col + (closing ? 0 : 1);
    else
      return context.indent + (closing ? 0 : cx.unit);
  },

  languageData: {
    commentTokens: {line: "#"}
  }
};



/***/ })

}]);
//# sourceMappingURL=5249.47203d8dad661b809e38.js.map?v=47203d8dad661b809e38
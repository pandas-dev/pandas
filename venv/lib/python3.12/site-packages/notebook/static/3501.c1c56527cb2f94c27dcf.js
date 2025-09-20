"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[3501],{

/***/ 73501:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   d: () => (/* binding */ d)
/* harmony export */ });
function words(str) {
  var obj = {}, words = str.split(" ");
  for (var i = 0; i < words.length; ++i) obj[words[i]] = true;
  return obj;
}

var blockKeywordsStr = "body catch class do else enum for foreach foreach_reverse if in interface mixin " +
    "out scope struct switch try union unittest version while with";

const parserConfig = {
  keywords: words("abstract alias align asm assert auto break case cast cdouble cent cfloat const continue " +
                  "debug default delegate delete deprecated export extern final finally function goto immutable " +
                  "import inout invariant is lazy macro module new nothrow override package pragma private " +
                  "protected public pure ref return shared short static super synchronized template this " +
                  "throw typedef typeid typeof volatile __FILE__ __LINE__ __gshared __traits __vector __parameters " +
                  blockKeywordsStr),
  blockKeywords: words(blockKeywordsStr),
  builtin: words("bool byte char creal dchar double float idouble ifloat int ireal long real short ubyte " +
                 "ucent uint ulong ushort wchar wstring void size_t sizediff_t"),
  atoms: words("exit failure success true false null"),
  hooks: {
    "@": function(stream, _state) {
      stream.eatWhile(/[\w\$_]/);
      return "meta";
    }
  }
}

var statementIndentUnit = parserConfig.statementIndentUnit,
    keywords = parserConfig.keywords,
    builtin = parserConfig.builtin,
    blockKeywords = parserConfig.blockKeywords,
    atoms = parserConfig.atoms,
    hooks = parserConfig.hooks,
    multiLineStrings = parserConfig.multiLineStrings;
var isOperatorChar = /[+\-*&%=<>!?|\/]/;

var curPunc;

function tokenBase(stream, state) {
  var ch = stream.next();
  if (hooks[ch]) {
    var result = hooks[ch](stream, state);
    if (result !== false) return result;
  }
  if (ch == '"' || ch == "'" || ch == "`") {
    state.tokenize = tokenString(ch);
    return state.tokenize(stream, state);
  }
  if (/[\[\]{}\(\),;\:\.]/.test(ch)) {
    curPunc = ch;
    return null;
  }
  if (/\d/.test(ch)) {
    stream.eatWhile(/[\w\.]/);
    return "number";
  }
  if (ch == "/") {
    if (stream.eat("+")) {
      state.tokenize = tokenNestedComment;
      return tokenNestedComment(stream, state);
    }
    if (stream.eat("*")) {
      state.tokenize = tokenComment;
      return tokenComment(stream, state);
    }
    if (stream.eat("/")) {
      stream.skipToEnd();
      return "comment";
    }
  }
  if (isOperatorChar.test(ch)) {
    stream.eatWhile(isOperatorChar);
    return "operator";
  }
  stream.eatWhile(/[\w\$_\xa1-\uffff]/);
  var cur = stream.current();
  if (keywords.propertyIsEnumerable(cur)) {
    if (blockKeywords.propertyIsEnumerable(cur)) curPunc = "newstatement";
    return "keyword";
  }
  if (builtin.propertyIsEnumerable(cur)) {
    if (blockKeywords.propertyIsEnumerable(cur)) curPunc = "newstatement";
    return "builtin";
  }
  if (atoms.propertyIsEnumerable(cur)) return "atom";
  return "variable";
}

function tokenString(quote) {
  return function(stream, state) {
    var escaped = false, next, end = false;
    while ((next = stream.next()) != null) {
      if (next == quote && !escaped) {end = true; break;}
      escaped = !escaped && next == "\\";
    }
    if (end || !(escaped || multiLineStrings))
      state.tokenize = null;
    return "string";
  };
}

function tokenComment(stream, state) {
  var maybeEnd = false, ch;
  while (ch = stream.next()) {
    if (ch == "/" && maybeEnd) {
      state.tokenize = null;
      break;
    }
    maybeEnd = (ch == "*");
  }
  return "comment";
}

function tokenNestedComment(stream, state) {
  var maybeEnd = false, ch;
  while (ch = stream.next()) {
    if (ch == "/" && maybeEnd) {
      state.tokenize = null;
      break;
    }
    maybeEnd = (ch == "+");
  }
  return "comment";
}

function Context(indented, column, type, align, prev) {
  this.indented = indented;
  this.column = column;
  this.type = type;
  this.align = align;
  this.prev = prev;
}
function pushContext(state, col, type) {
  var indent = state.indented;
  if (state.context && state.context.type == "statement")
    indent = state.context.indented;
  return state.context = new Context(indent, col, type, null, state.context);
}
function popContext(state) {
  var t = state.context.type;
  if (t == ")" || t == "]" || t == "}")
    state.indented = state.context.indented;
  return state.context = state.context.prev;
}

// Interface

const d = {
  name: "d",
  startState: function(indentUnit) {
    return {
      tokenize: null,
      context: new Context(-indentUnit, 0, "top", false),
      indented: 0,
      startOfLine: true
    };
  },

  token: function(stream, state) {
    var ctx = state.context;
    if (stream.sol()) {
      if (ctx.align == null) ctx.align = false;
      state.indented = stream.indentation();
      state.startOfLine = true;
    }
    if (stream.eatSpace()) return null;
    curPunc = null;
    var style = (state.tokenize || tokenBase)(stream, state);
    if (style == "comment" || style == "meta") return style;
    if (ctx.align == null) ctx.align = true;

    if ((curPunc == ";" || curPunc == ":" || curPunc == ",") && ctx.type == "statement") popContext(state);
    else if (curPunc == "{") pushContext(state, stream.column(), "}");
    else if (curPunc == "[") pushContext(state, stream.column(), "]");
    else if (curPunc == "(") pushContext(state, stream.column(), ")");
    else if (curPunc == "}") {
      while (ctx.type == "statement") ctx = popContext(state);
      if (ctx.type == "}") ctx = popContext(state);
      while (ctx.type == "statement") ctx = popContext(state);
    }
    else if (curPunc == ctx.type) popContext(state);
    else if (((ctx.type == "}" || ctx.type == "top") && curPunc != ';') || (ctx.type == "statement" && curPunc == "newstatement"))
      pushContext(state, stream.column(), "statement");
    state.startOfLine = false;
    return style;
  },

  indent: function(state, textAfter, cx) {
    if (state.tokenize != tokenBase && state.tokenize != null) return null;
    var ctx = state.context, firstChar = textAfter && textAfter.charAt(0);
    if (ctx.type == "statement" && firstChar == "}") ctx = ctx.prev;
    var closing = firstChar == ctx.type;
    if (ctx.type == "statement") return ctx.indented + (firstChar == "{" ? 0 : statementIndentUnit || cx.unit);
    else if (ctx.align) return ctx.column + (closing ? 0 : 1);
    else return ctx.indented + (closing ? 0 : cx.unit);
  },

  languageData: {
    indentOnInput: /^\s*[{}]$/,
    commentTokens: {line: "//", block: {open: "/*", close: "*/"}}
  }
};


/***/ })

}]);
//# sourceMappingURL=3501.c1c56527cb2f94c27dcf.js.map?v=c1c56527cb2f94c27dcf
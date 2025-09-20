"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[7360],{

/***/ 7360:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   groovy: () => (/* binding */ groovy)
/* harmony export */ });
function words(str) {
  var obj = {}, words = str.split(" ");
  for (var i = 0; i < words.length; ++i) obj[words[i]] = true;
  return obj;
}
var keywords = words(
  "abstract as assert boolean break byte case catch char class const continue def default " +
    "do double else enum extends final finally float for goto if implements import in " +
    "instanceof int interface long native new package private protected public return " +
    "short static strictfp super switch synchronized threadsafe throw throws trait transient " +
    "try void volatile while");
var blockKeywords = words("catch class def do else enum finally for if interface switch trait try while");
var standaloneKeywords = words("return break continue");
var atoms = words("null true false this");

var curPunc;
function tokenBase(stream, state) {
  var ch = stream.next();
  if (ch == '"' || ch == "'") {
    return startString(ch, stream, state);
  }
  if (/[\[\]{}\(\),;\:\.]/.test(ch)) {
    curPunc = ch;
    return null;
  }
  if (/\d/.test(ch)) {
    stream.eatWhile(/[\w\.]/);
    if (stream.eat(/eE/)) { stream.eat(/\+\-/); stream.eatWhile(/\d/); }
    return "number";
  }
  if (ch == "/") {
    if (stream.eat("*")) {
      state.tokenize.push(tokenComment);
      return tokenComment(stream, state);
    }
    if (stream.eat("/")) {
      stream.skipToEnd();
      return "comment";
    }
    if (expectExpression(state.lastToken, false)) {
      return startString(ch, stream, state);
    }
  }
  if (ch == "-" && stream.eat(">")) {
    curPunc = "->";
    return null;
  }
  if (/[+\-*&%=<>!?|\/~]/.test(ch)) {
    stream.eatWhile(/[+\-*&%=<>|~]/);
    return "operator";
  }
  stream.eatWhile(/[\w\$_]/);
  if (ch == "@") { stream.eatWhile(/[\w\$_\.]/); return "meta"; }
  if (state.lastToken == ".") return "property";
  if (stream.eat(":")) { curPunc = "proplabel"; return "property"; }
  var cur = stream.current();
  if (atoms.propertyIsEnumerable(cur)) { return "atom"; }
  if (keywords.propertyIsEnumerable(cur)) {
    if (blockKeywords.propertyIsEnumerable(cur)) curPunc = "newstatement";
    else if (standaloneKeywords.propertyIsEnumerable(cur)) curPunc = "standalone";
    return "keyword";
  }
  return "variable";
}
tokenBase.isBase = true;

function startString(quote, stream, state) {
  var tripleQuoted = false;
  if (quote != "/" && stream.eat(quote)) {
    if (stream.eat(quote)) tripleQuoted = true;
    else return "string";
  }
  function t(stream, state) {
    var escaped = false, next, end = !tripleQuoted;
    while ((next = stream.next()) != null) {
      if (next == quote && !escaped) {
        if (!tripleQuoted) { break; }
        if (stream.match(quote + quote)) { end = true; break; }
      }
      if (quote == '"' && next == "$" && !escaped) {
        if (stream.eat("{")) {
          state.tokenize.push(tokenBaseUntilBrace());
          return "string";
        } else if (stream.match(/^\w/, false)) {
          state.tokenize.push(tokenVariableDeref);
          return "string";
        }
      }
      escaped = !escaped && next == "\\";
    }
    if (end) state.tokenize.pop();
    return "string";
  }
  state.tokenize.push(t);
  return t(stream, state);
}

function tokenBaseUntilBrace() {
  var depth = 1;
  function t(stream, state) {
    if (stream.peek() == "}") {
      depth--;
      if (depth == 0) {
        state.tokenize.pop();
        return state.tokenize[state.tokenize.length-1](stream, state);
      }
    } else if (stream.peek() == "{") {
      depth++;
    }
    return tokenBase(stream, state);
  }
  t.isBase = true;
  return t;
}

function tokenVariableDeref(stream, state) {
  var next = stream.match(/^(\.|[\w\$_]+)/)
  if (!next || !stream.match(next[0] == "." ? /^[\w$_]/ : /^\./)) state.tokenize.pop()
  if (!next) return state.tokenize[state.tokenize.length-1](stream, state)
  return next[0] == "." ? null : "variable"
}

function tokenComment(stream, state) {
  var maybeEnd = false, ch;
  while (ch = stream.next()) {
    if (ch == "/" && maybeEnd) {
      state.tokenize.pop();
      break;
    }
    maybeEnd = (ch == "*");
  }
  return "comment";
}

function expectExpression(last, newline) {
  return !last || last == "operator" || last == "->" || /[\.\[\{\(,;:]/.test(last) ||
    last == "newstatement" || last == "keyword" || last == "proplabel" ||
    (last == "standalone" && !newline);
}

function Context(indented, column, type, align, prev) {
  this.indented = indented;
  this.column = column;
  this.type = type;
  this.align = align;
  this.prev = prev;
}
function pushContext(state, col, type) {
  return state.context = new Context(state.indented, col, type, null, state.context);
}
function popContext(state) {
  var t = state.context.type;
  if (t == ")" || t == "]" || t == "}")
    state.indented = state.context.indented;
  return state.context = state.context.prev;
}

// Interface

const groovy = {
  name: "groovy",
  startState: function(indentUnit) {
    return {
      tokenize: [tokenBase],
      context: new Context(-indentUnit, 0, "top", false),
      indented: 0,
      startOfLine: true,
      lastToken: null
    };
  },

  token: function(stream, state) {
    var ctx = state.context;
    if (stream.sol()) {
      if (ctx.align == null) ctx.align = false;
      state.indented = stream.indentation();
      state.startOfLine = true;
      // Automatic semicolon insertion
      if (ctx.type == "statement" && !expectExpression(state.lastToken, true)) {
        popContext(state); ctx = state.context;
      }
    }
    if (stream.eatSpace()) return null;
    curPunc = null;
    var style = state.tokenize[state.tokenize.length-1](stream, state);
    if (style == "comment") return style;
    if (ctx.align == null) ctx.align = true;

    if ((curPunc == ";" || curPunc == ":") && ctx.type == "statement") popContext(state);
    // Handle indentation for {x -> \n ... }
    else if (curPunc == "->" && ctx.type == "statement" && ctx.prev.type == "}") {
      popContext(state);
      state.context.align = false;
    }
    else if (curPunc == "{") pushContext(state, stream.column(), "}");
    else if (curPunc == "[") pushContext(state, stream.column(), "]");
    else if (curPunc == "(") pushContext(state, stream.column(), ")");
    else if (curPunc == "}") {
      while (ctx.type == "statement") ctx = popContext(state);
      if (ctx.type == "}") ctx = popContext(state);
      while (ctx.type == "statement") ctx = popContext(state);
    }
    else if (curPunc == ctx.type) popContext(state);
    else if (ctx.type == "}" || ctx.type == "top" || (ctx.type == "statement" && curPunc == "newstatement"))
      pushContext(state, stream.column(), "statement");
    state.startOfLine = false;
    state.lastToken = curPunc || style;
    return style;
  },

  indent: function(state, textAfter, cx) {
    if (!state.tokenize[state.tokenize.length-1].isBase) return null;
    var firstChar = textAfter && textAfter.charAt(0), ctx = state.context;
    if (ctx.type == "statement" && !expectExpression(state.lastToken, true)) ctx = ctx.prev;
    var closing = firstChar == ctx.type;
    if (ctx.type == "statement") return ctx.indented + (firstChar == "{" ? 0 : cx.unit);
    else if (ctx.align) return ctx.column + (closing ? 0 : 1);
    else return ctx.indented + (closing ? 0 : cx.unit);
  },

  languageData: {
    indentOnInput: /^\s*[{}]$/,
    commentTokens: {line: "//", block: {open: "/*", close: "*/"}},
    closeBrackets: {brackets: ["(", "[", "{", "'", '"', "'''", '"""']}
  }
};


/***/ })

}]);
//# sourceMappingURL=7360.b3741cc7257cecd9efe9.js.map?v=b3741cc7257cecd9efe9
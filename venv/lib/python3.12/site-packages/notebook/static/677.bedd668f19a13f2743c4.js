"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[677],{

/***/ 30677:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   vhdl: () => (/* binding */ vhdl)
/* harmony export */ });
function words(str) {
  var obj = {}, words = str.split(",");
  for (var i = 0; i < words.length; ++i) {
    var allCaps = words[i].toUpperCase();
    var firstCap = words[i].charAt(0).toUpperCase() + words[i].slice(1);
    obj[words[i]] = true;
    obj[allCaps] = true;
    obj[firstCap] = true;
  }
  return obj;
}

function metaHook(stream) {
  stream.eatWhile(/[\w\$_]/);
  return "meta";
}

var atoms = words("null"),
    hooks = {"`": metaHook, "$": metaHook},
    multiLineStrings = false;

var keywords = words("abs,access,after,alias,all,and,architecture,array,assert,attribute,begin,block," +
                     "body,buffer,bus,case,component,configuration,constant,disconnect,downto,else,elsif,end,end block,end case," +
                     "end component,end for,end generate,end if,end loop,end process,end record,end units,entity,exit,file,for," +
                     "function,generate,generic,generic map,group,guarded,if,impure,in,inertial,inout,is,label,library,linkage," +
                     "literal,loop,map,mod,nand,new,next,nor,null,of,on,open,or,others,out,package,package body,port,port map," +
                     "postponed,procedure,process,pure,range,record,register,reject,rem,report,return,rol,ror,select,severity,signal," +
                     "sla,sll,sra,srl,subtype,then,to,transport,type,unaffected,units,until,use,variable,wait,when,while,with,xnor,xor");

var blockKeywords = words("architecture,entity,begin,case,port,else,elsif,end,for,function,if");

var isOperatorChar = /[&|~><!\)\(*#%@+\/=?\:;}{,\.\^\-\[\]]/;
var curPunc;

function tokenBase(stream, state) {
  var ch = stream.next();
  if (hooks[ch]) {
    var result = hooks[ch](stream, state);
    if (result !== false) return result;
  }
  if (ch == '"') {
    state.tokenize = tokenString2(ch);
    return state.tokenize(stream, state);
  }
  if (ch == "'") {
    state.tokenize = tokenString(ch);
    return state.tokenize(stream, state);
  }
  if (/[\[\]{}\(\),;\:\.]/.test(ch)) {
    curPunc = ch;
    return null;
  }
  if (/[\d']/.test(ch)) {
    stream.eatWhile(/[\w\.']/);
    return "number";
  }
  if (ch == "-") {
    if (stream.eat("-")) {
      stream.skipToEnd();
      return "comment";
    }
  }
  if (isOperatorChar.test(ch)) {
    stream.eatWhile(isOperatorChar);
    return "operator";
  }
  stream.eatWhile(/[\w\$_]/);
  var cur = stream.current();
  if (keywords.propertyIsEnumerable(cur.toLowerCase())) {
    if (blockKeywords.propertyIsEnumerable(cur)) curPunc = "newstatement";
    return "keyword";
  }
  if (atoms.propertyIsEnumerable(cur)) return "atom";
  return "variable";
}

function tokenString(quote) {
  return function(stream, state) {
    var escaped = false, next, end = false;
    while ((next = stream.next()) != null) {
      if (next == quote && !escaped) {end = true; break;}
      escaped = !escaped && next == "--";
    }
    if (end || !(escaped || multiLineStrings))
      state.tokenize = tokenBase;
    return "string";
  };
}
function tokenString2(quote) {
  return function(stream, state) {
    var escaped = false, next, end = false;
    while ((next = stream.next()) != null) {
      if (next == quote && !escaped) {end = true; break;}
      escaped = !escaped && next == "--";
    }
    if (end || !(escaped || multiLineStrings))
      state.tokenize = tokenBase;
    return "string.special";
  };
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
const vhdl = {
  name: "vhdl",
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

    if ((curPunc == ";" || curPunc == ":") && ctx.type == "statement") popContext(state);
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
    return style;
  },

  indent: function(state, textAfter, cx) {
    if (state.tokenize != tokenBase && state.tokenize != null) return 0;
    var firstChar = textAfter && textAfter.charAt(0), ctx = state.context, closing = firstChar == ctx.type;
    if (ctx.type == "statement") return ctx.indented + (firstChar == "{" ? 0 : cx.unit);
    else if (ctx.align) return ctx.column + (closing ? 0 : 1);
    else return ctx.indented + (closing ? 0 : cx.unit);
  },

  languageData: {
    indentOnInput: /^\s*[{}]$/,
    commentTokens: {line: "--"}
  }
}


/***/ })

}]);
//# sourceMappingURL=677.bedd668f19a13f2743c4.js.map?v=bedd668f19a13f2743c4
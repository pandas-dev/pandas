"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[8381],{

/***/ 68381:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   pig: () => (/* binding */ pig)
/* harmony export */ });
function words(str) {
  var obj = {}, words = str.split(" ");
  for (var i = 0; i < words.length; ++i) obj[words[i]] = true;
  return obj;
}

// builtin funcs taken from trunk revision 1303237
var pBuiltins = "ABS ACOS ARITY ASIN ATAN AVG BAGSIZE BINSTORAGE BLOOM BUILDBLOOM CBRT CEIL "
    + "CONCAT COR COS COSH COUNT COUNT_STAR COV CONSTANTSIZE CUBEDIMENSIONS DIFF DISTINCT DOUBLEABS "
    + "DOUBLEAVG DOUBLEBASE DOUBLEMAX DOUBLEMIN DOUBLEROUND DOUBLESUM EXP FLOOR FLOATABS FLOATAVG "
    + "FLOATMAX FLOATMIN FLOATROUND FLOATSUM GENERICINVOKER INDEXOF INTABS INTAVG INTMAX INTMIN "
    + "INTSUM INVOKEFORDOUBLE INVOKEFORFLOAT INVOKEFORINT INVOKEFORLONG INVOKEFORSTRING INVOKER "
    + "ISEMPTY JSONLOADER JSONMETADATA JSONSTORAGE LAST_INDEX_OF LCFIRST LOG LOG10 LOWER LONGABS "
    + "LONGAVG LONGMAX LONGMIN LONGSUM MAX MIN MAPSIZE MONITOREDUDF NONDETERMINISTIC OUTPUTSCHEMA  "
    + "PIGSTORAGE PIGSTREAMING RANDOM REGEX_EXTRACT REGEX_EXTRACT_ALL REPLACE ROUND SIN SINH SIZE "
    + "SQRT STRSPLIT SUBSTRING SUM STRINGCONCAT STRINGMAX STRINGMIN STRINGSIZE TAN TANH TOBAG "
    + "TOKENIZE TOMAP TOP TOTUPLE TRIM TEXTLOADER TUPLESIZE UCFIRST UPPER UTF8STORAGECONVERTER ";

// taken from QueryLexer.g
var pKeywords = "VOID IMPORT RETURNS DEFINE LOAD FILTER FOREACH ORDER CUBE DISTINCT COGROUP "
    + "JOIN CROSS UNION SPLIT INTO IF OTHERWISE ALL AS BY USING INNER OUTER ONSCHEMA PARALLEL "
    + "PARTITION GROUP AND OR NOT GENERATE FLATTEN ASC DESC IS STREAM THROUGH STORE MAPREDUCE "
    + "SHIP CACHE INPUT OUTPUT STDERROR STDIN STDOUT LIMIT SAMPLE LEFT RIGHT FULL EQ GT LT GTE LTE "
    + "NEQ MATCHES TRUE FALSE DUMP";

// data types
var pTypes = "BOOLEAN INT LONG FLOAT DOUBLE CHARARRAY BYTEARRAY BAG TUPLE MAP ";

var builtins = words(pBuiltins), keywords = words(pKeywords), types = words(pTypes)

var isOperatorChar = /[*+\-%<>=&?:\/!|]/;

function chain(stream, state, f) {
  state.tokenize = f;
  return f(stream, state);
}

function tokenComment(stream, state) {
  var isEnd = false;
  var ch;
  while(ch = stream.next()) {
    if(ch == "/" && isEnd) {
      state.tokenize = tokenBase;
      break;
    }
    isEnd = (ch == "*");
  }
  return "comment";
}

function tokenString(quote) {
  return function(stream, state) {
    var escaped = false, next, end = false;
    while((next = stream.next()) != null) {
      if (next == quote && !escaped) {
        end = true; break;
      }
      escaped = !escaped && next == "\\";
    }
    if (end || !escaped)
      state.tokenize = tokenBase;
    return "error";
  };
}


function tokenBase(stream, state) {
  var ch = stream.next();

  // is a start of string?
  if (ch == '"' || ch == "'")
    return chain(stream, state, tokenString(ch));
  // is it one of the special chars
  else if(/[\[\]{}\(\),;\.]/.test(ch))
    return null;
  // is it a number?
  else if(/\d/.test(ch)) {
    stream.eatWhile(/[\w\.]/);
    return "number";
  }
  // multi line comment or operator
  else if (ch == "/") {
    if (stream.eat("*")) {
      return chain(stream, state, tokenComment);
    }
    else {
      stream.eatWhile(isOperatorChar);
      return "operator";
    }
  }
  // single line comment or operator
  else if (ch=="-") {
    if(stream.eat("-")){
      stream.skipToEnd();
      return "comment";
    }
    else {
      stream.eatWhile(isOperatorChar);
      return "operator";
    }
  }
  // is it an operator
  else if (isOperatorChar.test(ch)) {
    stream.eatWhile(isOperatorChar);
    return "operator";
  }
  else {
    // get the while word
    stream.eatWhile(/[\w\$_]/);
    // is it one of the listed keywords?
    if (keywords && keywords.propertyIsEnumerable(stream.current().toUpperCase())) {
      //keywords can be used as variables like flatten(group), group.$0 etc..
      if (!stream.eat(")") && !stream.eat("."))
        return "keyword";
    }
    // is it one of the builtin functions?
    if (builtins && builtins.propertyIsEnumerable(stream.current().toUpperCase()))
      return "builtin";
    // is it one of the listed types?
    if (types && types.propertyIsEnumerable(stream.current().toUpperCase()))
      return "type";
    // default is a 'variable'
    return "variable";
  }
}

// Interface
const pig = {
  name: "pig",

  startState: function() {
    return {
      tokenize: tokenBase,
      startOfLine: true
    };
  },

  token: function(stream, state) {
    if(stream.eatSpace()) return null;
    var style = state.tokenize(stream, state);
    return style;
  },

  languageData: {
    autocomplete: (pBuiltins + pTypes + pKeywords).split(" ")
  }
};


/***/ })

}]);
//# sourceMappingURL=8381.0291906ada65d4e5df4e.js.map?v=0291906ada65d4e5df4e
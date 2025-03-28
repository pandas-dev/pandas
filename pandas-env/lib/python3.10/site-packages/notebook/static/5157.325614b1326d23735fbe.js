"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[5157],{

/***/ 97534:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   asn1: () => (/* binding */ asn1)
/* harmony export */ });
function words(str) {
  var obj = {}, words = str.split(" ");
  for (var i = 0; i < words.length; ++i) obj[words[i]] = true;
  return obj;
}

const defaults = {
  keywords: words("DEFINITIONS OBJECTS IF DERIVED INFORMATION ACTION" +
                  " REPLY ANY NAMED CHARACTERIZED BEHAVIOUR REGISTERED" +
                  " WITH AS IDENTIFIED CONSTRAINED BY PRESENT BEGIN" +
                  " IMPORTS FROM UNITS SYNTAX MIN-ACCESS MAX-ACCESS" +
                  " MINACCESS MAXACCESS REVISION STATUS DESCRIPTION" +
                  " SEQUENCE SET COMPONENTS OF CHOICE DistinguishedName" +
                  " ENUMERATED SIZE MODULE END INDEX AUGMENTS EXTENSIBILITY" +
                  " IMPLIED EXPORTS"),
  cmipVerbs: words("ACTIONS ADD GET NOTIFICATIONS REPLACE REMOVE"),
  compareTypes: words("OPTIONAL DEFAULT MANAGED MODULE-TYPE MODULE_IDENTITY" +
                      " MODULE-COMPLIANCE OBJECT-TYPE OBJECT-IDENTITY" +
                      " OBJECT-COMPLIANCE MODE CONFIRMED CONDITIONAL" +
                      " SUBORDINATE SUPERIOR CLASS TRUE FALSE NULL" +
                      " TEXTUAL-CONVENTION"),
  status: words("current deprecated mandatory obsolete"),
  tags: words("APPLICATION AUTOMATIC EXPLICIT IMPLICIT PRIVATE TAGS" +
              " UNIVERSAL"),
  storage: words("BOOLEAN INTEGER OBJECT IDENTIFIER BIT OCTET STRING" +
                 " UTCTime InterfaceIndex IANAifType CMIP-Attribute" +
                 " REAL PACKAGE PACKAGES IpAddress PhysAddress" +
                 " NetworkAddress BITS BMPString TimeStamp TimeTicks" +
                 " TruthValue RowStatus DisplayString GeneralString" +
                 " GraphicString IA5String NumericString" +
                 " PrintableString SnmpAdminString TeletexString" +
                 " UTF8String VideotexString VisibleString StringStore" +
                 " ISO646String T61String UniversalString Unsigned32" +
                 " Integer32 Gauge Gauge32 Counter Counter32 Counter64"),
  modifier: words("ATTRIBUTE ATTRIBUTES MANDATORY-GROUP MANDATORY-GROUPS" +
                  " GROUP GROUPS ELEMENTS EQUALITY ORDERING SUBSTRINGS" +
                  " DEFINED"),
  accessTypes: words("not-accessible accessible-for-notify read-only" +
                     " read-create read-write"),
  multiLineStrings: true
}

function asn1(parserConfig) {
  var keywords = parserConfig.keywords || defaults.keywords,
      cmipVerbs = parserConfig.cmipVerbs || defaults.cmipVerbs,
      compareTypes = parserConfig.compareTypes || defaults.compareTypes,
      status = parserConfig.status || defaults.status,
      tags = parserConfig.tags || defaults.tags,
      storage = parserConfig.storage || defaults.storage,
      modifier = parserConfig.modifier || defaults.modifier,
      accessTypes = parserConfig.accessTypes|| defaults.accessTypes,
      multiLineStrings = parserConfig.multiLineStrings || defaults.multiLineStrings,
      indentStatements = parserConfig.indentStatements !== false;
  var isOperatorChar = /[\|\^]/;
  var curPunc;

  function tokenBase(stream, state) {
    var ch = stream.next();
    if (ch == '"' || ch == "'") {
      state.tokenize = tokenString(ch);
      return state.tokenize(stream, state);
    }
    if (/[\[\]\(\){}:=,;]/.test(ch)) {
      curPunc = ch;
      return "punctuation";
    }
    if (ch == "-"){
      if (stream.eat("-")) {
        stream.skipToEnd();
        return "comment";
      }
    }
    if (/\d/.test(ch)) {
      stream.eatWhile(/[\w\.]/);
      return "number";
    }
    if (isOperatorChar.test(ch)) {
      stream.eatWhile(isOperatorChar);
      return "operator";
    }

    stream.eatWhile(/[\w\-]/);
    var cur = stream.current();
    if (keywords.propertyIsEnumerable(cur)) return "keyword";
    if (cmipVerbs.propertyIsEnumerable(cur)) return "variableName";
    if (compareTypes.propertyIsEnumerable(cur)) return "atom";
    if (status.propertyIsEnumerable(cur)) return "comment";
    if (tags.propertyIsEnumerable(cur)) return "typeName";
    if (storage.propertyIsEnumerable(cur)) return "modifier";
    if (modifier.propertyIsEnumerable(cur)) return "modifier";
    if (accessTypes.propertyIsEnumerable(cur)) return "modifier";

    return "variableName";
  }

  function tokenString(quote) {
    return function(stream, state) {
      var escaped = false, next, end = false;
      while ((next = stream.next()) != null) {
        if (next == quote && !escaped){
          var afterNext = stream.peek();
          //look if the character if the quote is like the B in '10100010'B
          if (afterNext){
            afterNext = afterNext.toLowerCase();
            if(afterNext == "b" || afterNext == "h" || afterNext == "o")
              stream.next();
          }
          end = true; break;
        }
        escaped = !escaped && next == "\\";
      }
      if (end || !(escaped || multiLineStrings))
        state.tokenize = null;
      return "string";
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

  //Interface
  return {
    name: "asn1",
    startState: function() {
      return {
        tokenize: null,
        context: new Context(-2, 0, "top", false),
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
      if (style == "comment") return style;
      if (ctx.align == null) ctx.align = true;

      if ((curPunc == ";" || curPunc == ":" || curPunc == ",")
          && ctx.type == "statement"){
        popContext(state);
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
      else if (indentStatements && (((ctx.type == "}" || ctx.type == "top")
                                     && curPunc != ';') || (ctx.type == "statement"
                                                            && curPunc == "newstatement")))
        pushContext(state, stream.column(), "statement");

      state.startOfLine = false;
      return style;
    },

    languageData: {
      indentOnInput: /^\s*[{}]$/,
      commentTokens: {line: "--"}
    }
  };
};


/***/ })

}]);
//# sourceMappingURL=5157.325614b1326d23735fbe.js.map?v=325614b1326d23735fbe
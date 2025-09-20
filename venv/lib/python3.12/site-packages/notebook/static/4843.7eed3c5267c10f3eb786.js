"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[4843],{

/***/ 64843:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   mscgen: () => (/* binding */ mscgen),
/* harmony export */   msgenny: () => (/* binding */ msgenny),
/* harmony export */   xu: () => (/* binding */ xu)
/* harmony export */ });
function mkParser(lang) {
  return {
    name: "mscgen",
    startState: startStateFn,
    copyState: copyStateFn,
    token: produceTokenFunction(lang),
    languageData: {
      commentTokens: {line: "#", block: {open: "/*", close: "*/"}}
    }
  }
}

const mscgen = mkParser({
  "keywords" : ["msc"],
  "options" : ["hscale", "width", "arcgradient", "wordwraparcs"],
  "constants" : ["true", "false", "on", "off"],
  "attributes" : ["label", "idurl", "id", "url", "linecolor", "linecolour", "textcolor", "textcolour", "textbgcolor", "textbgcolour", "arclinecolor", "arclinecolour", "arctextcolor", "arctextcolour", "arctextbgcolor", "arctextbgcolour", "arcskip"],
  "brackets" : ["\\{", "\\}"], // [ and  ] are brackets too, but these get handled in with lists
  "arcsWords" : ["note", "abox", "rbox", "box"],
  "arcsOthers" : ["\\|\\|\\|", "\\.\\.\\.", "---", "--", "<->", "==", "<<=>>", "<=>", "\\.\\.", "<<>>", "::", "<:>", "->", "=>>", "=>", ">>", ":>", "<-", "<<=", "<=", "<<", "<:", "x-", "-x"],
  "singlecomment" : ["//", "#"],
  "operators" : ["="]
})

const msgenny = mkParser({
  "keywords" : null,
  "options" : ["hscale", "width", "arcgradient", "wordwraparcs", "wordwrapentities", "watermark"],
  "constants" : ["true", "false", "on", "off", "auto"],
  "attributes" : null,
  "brackets" : ["\\{", "\\}"],
  "arcsWords" : ["note", "abox", "rbox", "box", "alt", "else", "opt", "break", "par", "seq", "strict", "neg", "critical", "ignore", "consider", "assert", "loop", "ref", "exc"],
  "arcsOthers" : ["\\|\\|\\|", "\\.\\.\\.", "---", "--", "<->", "==", "<<=>>", "<=>", "\\.\\.", "<<>>", "::", "<:>", "->", "=>>", "=>", ">>", ":>", "<-", "<<=", "<=", "<<", "<:", "x-", "-x"],
  "singlecomment" : ["//", "#"],
  "operators" : ["="]
})

const xu = mkParser({
  "keywords" : ["msc", "xu"],
  "options" : ["hscale", "width", "arcgradient", "wordwraparcs", "wordwrapentities", "watermark"],
  "constants" : ["true", "false", "on", "off", "auto"],
  "attributes" : ["label", "idurl", "id", "url", "linecolor", "linecolour", "textcolor", "textcolour", "textbgcolor", "textbgcolour", "arclinecolor", "arclinecolour", "arctextcolor", "arctextcolour", "arctextbgcolor", "arctextbgcolour", "arcskip", "title", "deactivate", "activate", "activation"],
  "brackets" : ["\\{", "\\}"],  // [ and  ] are brackets too, but these get handled in with lists
  "arcsWords" : ["note", "abox", "rbox", "box", "alt", "else", "opt", "break", "par", "seq", "strict", "neg", "critical", "ignore", "consider", "assert", "loop", "ref", "exc"],
  "arcsOthers" : ["\\|\\|\\|", "\\.\\.\\.", "---", "--", "<->", "==", "<<=>>", "<=>", "\\.\\.", "<<>>", "::", "<:>", "->", "=>>", "=>", ">>", ":>", "<-", "<<=", "<=", "<<", "<:", "x-", "-x"],
  "singlecomment" : ["//", "#"],
  "operators" : ["="]
})

function wordRegexpBoundary(pWords) {
  return new RegExp("^\\b(" + pWords.join("|") + ")\\b", "i");
}

function wordRegexp(pWords) {
  return new RegExp("^(?:" + pWords.join("|") + ")", "i");
}

function startStateFn() {
  return {
    inComment : false,
    inString : false,
    inAttributeList : false,
    inScript : false
  };
}

function copyStateFn(pState) {
  return {
    inComment : pState.inComment,
    inString : pState.inString,
    inAttributeList : pState.inAttributeList,
    inScript : pState.inScript
  };
}

function produceTokenFunction(pConfig) {
  return function(pStream, pState) {
    if (pStream.match(wordRegexp(pConfig.brackets), true, true)) {
      return "bracket";
    }
    /* comments */
    if (!pState.inComment) {
      if (pStream.match(/\/\*[^\*\/]*/, true, true)) {
        pState.inComment = true;
        return "comment";
      }
      if (pStream.match(wordRegexp(pConfig.singlecomment), true, true)) {
        pStream.skipToEnd();
        return "comment";
      }
    }
    if (pState.inComment) {
      if (pStream.match(/[^\*\/]*\*\//, true, true))
        pState.inComment = false;
      else
        pStream.skipToEnd();
      return "comment";
    }
    /* strings */
    if (!pState.inString && pStream.match(/\"(\\\"|[^\"])*/, true, true)) {
      pState.inString = true;
      return "string";
    }
    if (pState.inString) {
      if (pStream.match(/[^\"]*\"/, true, true))
        pState.inString = false;
      else
        pStream.skipToEnd();
      return "string";
    }
    /* keywords & operators */
    if (!!pConfig.keywords && pStream.match(wordRegexpBoundary(pConfig.keywords), true, true))
      return "keyword";

    if (pStream.match(wordRegexpBoundary(pConfig.options), true, true))
      return "keyword";

    if (pStream.match(wordRegexpBoundary(pConfig.arcsWords), true, true))
      return "keyword";

    if (pStream.match(wordRegexp(pConfig.arcsOthers), true, true))
      return "keyword";

    if (!!pConfig.operators && pStream.match(wordRegexp(pConfig.operators), true, true))
      return "operator";

    if (!!pConfig.constants && pStream.match(wordRegexp(pConfig.constants), true, true))
      return "variable";

    /* attribute lists */
    if (!pConfig.inAttributeList && !!pConfig.attributes && pStream.match('[', true, true)) {
      pConfig.inAttributeList = true;
      return "bracket";
    }
    if (pConfig.inAttributeList) {
      if (pConfig.attributes !== null && pStream.match(wordRegexpBoundary(pConfig.attributes), true, true)) {
        return "attribute";
      }
      if (pStream.match(']', true, true)) {
        pConfig.inAttributeList = false;
        return "bracket";
      }
    }

    pStream.next();
    return null
  };
}


/***/ })

}]);
//# sourceMappingURL=4843.7eed3c5267c10f3eb786.js.map?v=7eed3c5267c10f3eb786
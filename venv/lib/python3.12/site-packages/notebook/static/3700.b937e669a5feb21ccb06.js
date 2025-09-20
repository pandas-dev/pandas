"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[3700],{

/***/ 13700:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   rpmChanges: () => (/* binding */ rpmChanges),
/* harmony export */   rpmSpec: () => (/* binding */ rpmSpec)
/* harmony export */ });
var headerSeparator = /^-+$/;
var headerLine = /^(Mon|Tue|Wed|Thu|Fri|Sat|Sun) (Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)  ?\d{1,2} \d{2}:\d{2}(:\d{2})? [A-Z]{3,4} \d{4} - /;
var simpleEmail = /^[\w+.-]+@[\w.-]+/;

const rpmChanges = {
  name: "rpmchanges",
  token: function(stream) {
    if (stream.sol()) {
      if (stream.match(headerSeparator)) { return 'tag'; }
      if (stream.match(headerLine)) { return 'tag'; }
    }
    if (stream.match(simpleEmail)) { return 'string'; }
    stream.next();
    return null;
  }
}

// Quick and dirty spec file highlighting

var arch = /^(i386|i586|i686|x86_64|ppc64le|ppc64|ppc|ia64|s390x|s390|sparc64|sparcv9|sparc|noarch|alphaev6|alpha|hppa|mipsel)/;

var preamble = /^[a-zA-Z0-9()]+:/;
var section = /^%(debug_package|package|description|prep|build|install|files|clean|changelog|preinstall|preun|postinstall|postun|pretrans|posttrans|pre|post|triggerin|triggerun|verifyscript|check|triggerpostun|triggerprein|trigger)/;
var control_flow_complex = /^%(ifnarch|ifarch|if)/; // rpm control flow macros
var control_flow_simple = /^%(else|endif)/; // rpm control flow macros
var operators = /^(\!|\?|\<\=|\<|\>\=|\>|\=\=|\&\&|\|\|)/; // operators in control flow macros

const rpmSpec = {
  name: "rpmspec",
  startState: function () {
    return {
      controlFlow: false,
      macroParameters: false,
      section: false
    };
  },
  token: function (stream, state) {
    var ch = stream.peek();
    if (ch == "#") { stream.skipToEnd(); return "comment"; }

    if (stream.sol()) {
      if (stream.match(preamble)) { return "header"; }
      if (stream.match(section)) { return "atom"; }
    }

    if (stream.match(/^\$\w+/)) { return "def"; } // Variables like '$RPM_BUILD_ROOT'
    if (stream.match(/^\$\{\w+\}/)) { return "def"; } // Variables like '${RPM_BUILD_ROOT}'

    if (stream.match(control_flow_simple)) { return "keyword"; }
    if (stream.match(control_flow_complex)) {
      state.controlFlow = true;
      return "keyword";
    }
    if (state.controlFlow) {
      if (stream.match(operators)) { return "operator"; }
      if (stream.match(/^(\d+)/)) { return "number"; }
      if (stream.eol()) { state.controlFlow = false; }
    }

    if (stream.match(arch)) {
      if (stream.eol()) { state.controlFlow = false; }
      return "number";
    }

    // Macros like '%make_install' or '%attr(0775,root,root)'
    if (stream.match(/^%[\w]+/)) {
      if (stream.match('(')) { state.macroParameters = true; }
      return "keyword";
    }
    if (state.macroParameters) {
      if (stream.match(/^\d+/)) { return "number";}
      if (stream.match(')')) {
        state.macroParameters = false;
        return "keyword";
      }
    }

    // Macros like '%{defined fedora}'
    if (stream.match(/^%\{\??[\w \-\:\!]+\}/)) {
      if (stream.eol()) { state.controlFlow = false; }
      return "def";
    }

    stream.next();
    return null;
  }
};



/***/ })

}]);
//# sourceMappingURL=3700.b937e669a5feb21ccb06.js.map?v=b937e669a5feb21ccb06
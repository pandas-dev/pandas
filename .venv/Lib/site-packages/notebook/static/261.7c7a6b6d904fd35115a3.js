"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[261],{

/***/ 80261:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(24028);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(23617);
/* harmony import */ var dayjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(27693);
/* harmony import */ var dayjs__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(dayjs__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(7608);
/* harmony import */ var dompurify__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(31699);
/* harmony import */ var dompurify__WEBPACK_IMPORTED_MODULE_3___default = /*#__PURE__*/__webpack_require__.n(dompurify__WEBPACK_IMPORTED_MODULE_3__);











var parser = function() {
  var o = function(k, v, o2, l) {
    for (o2 = o2 || {}, l = k.length; l--; o2[k[l]] = v)
      ;
    return o2;
  }, $V0 = [1, 3], $V1 = [1, 6], $V2 = [1, 4], $V3 = [1, 5], $V4 = [2, 5], $V5 = [1, 12], $V6 = [5, 7, 13, 19, 21, 23, 24, 26, 28, 31, 37, 40, 47], $V7 = [7, 13, 19, 21, 23, 24, 26, 28, 31, 37, 40], $V8 = [7, 12, 13, 19, 21, 23, 24, 26, 28, 31, 37, 40], $V9 = [7, 13, 47], $Va = [1, 42], $Vb = [1, 41], $Vc = [7, 13, 29, 32, 35, 38, 47], $Vd = [1, 55], $Ve = [1, 56], $Vf = [1, 57], $Vg = [7, 13, 32, 35, 42, 47];
  var parser2 = {
    trace: function trace() {
    },
    yy: {},
    symbols_: { "error": 2, "start": 3, "eol": 4, "GG": 5, "document": 6, "EOF": 7, ":": 8, "DIR": 9, "options": 10, "body": 11, "OPT": 12, "NL": 13, "line": 14, "statement": 15, "commitStatement": 16, "mergeStatement": 17, "cherryPickStatement": 18, "acc_title": 19, "acc_title_value": 20, "acc_descr": 21, "acc_descr_value": 22, "acc_descr_multiline_value": 23, "section": 24, "branchStatement": 25, "CHECKOUT": 26, "ref": 27, "BRANCH": 28, "ORDER": 29, "NUM": 30, "CHERRY_PICK": 31, "COMMIT_ID": 32, "STR": 33, "PARENT_COMMIT": 34, "COMMIT_TAG": 35, "EMPTYSTR": 36, "MERGE": 37, "COMMIT_TYPE": 38, "commitType": 39, "COMMIT": 40, "commit_arg": 41, "COMMIT_MSG": 42, "NORMAL": 43, "REVERSE": 44, "HIGHLIGHT": 45, "ID": 46, ";": 47, "$accept": 0, "$end": 1 },
    terminals_: { 2: "error", 5: "GG", 7: "EOF", 8: ":", 9: "DIR", 12: "OPT", 13: "NL", 19: "acc_title", 20: "acc_title_value", 21: "acc_descr", 22: "acc_descr_value", 23: "acc_descr_multiline_value", 24: "section", 26: "CHECKOUT", 28: "BRANCH", 29: "ORDER", 30: "NUM", 31: "CHERRY_PICK", 32: "COMMIT_ID", 33: "STR", 34: "PARENT_COMMIT", 35: "COMMIT_TAG", 36: "EMPTYSTR", 37: "MERGE", 38: "COMMIT_TYPE", 40: "COMMIT", 42: "COMMIT_MSG", 43: "NORMAL", 44: "REVERSE", 45: "HIGHLIGHT", 46: "ID", 47: ";" },
    productions_: [0, [3, 2], [3, 3], [3, 4], [3, 5], [6, 0], [6, 2], [10, 2], [10, 1], [11, 0], [11, 2], [14, 2], [14, 1], [15, 1], [15, 1], [15, 1], [15, 2], [15, 2], [15, 1], [15, 1], [15, 1], [15, 2], [25, 2], [25, 4], [18, 3], [18, 5], [18, 5], [18, 7], [18, 7], [18, 5], [18, 5], [18, 5], [18, 7], [18, 7], [18, 7], [18, 7], [17, 2], [17, 4], [17, 4], [17, 4], [17, 6], [17, 6], [17, 6], [17, 6], [17, 6], [17, 6], [17, 8], [17, 8], [17, 8], [17, 8], [17, 8], [17, 8], [16, 2], [16, 3], [16, 3], [16, 5], [16, 5], [16, 3], [16, 5], [16, 5], [16, 5], [16, 5], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 3], [16, 5], [16, 5], [16, 5], [16, 5], [16, 5], [16, 5], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 7], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [16, 9], [41, 0], [41, 1], [39, 1], [39, 1], [39, 1], [27, 1], [27, 1], [4, 1], [4, 1], [4, 1]],
    performAction: function anonymous(yytext, yyleng, yylineno, yy, yystate, $$, _$) {
      var $0 = $$.length - 1;
      switch (yystate) {
        case 2:
          return $$[$0];
        case 3:
          return $$[$0 - 1];
        case 4:
          yy.setDirection($$[$0 - 3]);
          return $$[$0 - 1];
        case 6:
          yy.setOptions($$[$0 - 1]);
          this.$ = $$[$0];
          break;
        case 7:
          $$[$0 - 1] += $$[$0];
          this.$ = $$[$0 - 1];
          break;
        case 9:
          this.$ = [];
          break;
        case 10:
          $$[$0 - 1].push($$[$0]);
          this.$ = $$[$0 - 1];
          break;
        case 11:
          this.$ = $$[$0 - 1];
          break;
        case 16:
          this.$ = $$[$0].trim();
          yy.setAccTitle(this.$);
          break;
        case 17:
        case 18:
          this.$ = $$[$0].trim();
          yy.setAccDescription(this.$);
          break;
        case 19:
          yy.addSection($$[$0].substr(8));
          this.$ = $$[$0].substr(8);
          break;
        case 21:
          yy.checkout($$[$0]);
          break;
        case 22:
          yy.branch($$[$0]);
          break;
        case 23:
          yy.branch($$[$0 - 2], $$[$0]);
          break;
        case 24:
          yy.cherryPick($$[$0], "", void 0);
          break;
        case 25:
          yy.cherryPick($$[$0 - 2], "", void 0, $$[$0]);
          break;
        case 26:
          yy.cherryPick($$[$0 - 2], "", $$[$0]);
          break;
        case 27:
          yy.cherryPick($$[$0 - 4], "", $$[$0], $$[$0 - 2]);
          break;
        case 28:
          yy.cherryPick($$[$0 - 4], "", $$[$0 - 2], $$[$0]);
          break;
        case 29:
          yy.cherryPick($$[$0], "", $$[$0 - 2]);
          break;
        case 30:
          yy.cherryPick($$[$0], "", "");
          break;
        case 31:
          yy.cherryPick($$[$0 - 2], "", "");
          break;
        case 32:
          yy.cherryPick($$[$0 - 4], "", "", $$[$0 - 2]);
          break;
        case 33:
          yy.cherryPick($$[$0 - 4], "", "", $$[$0]);
          break;
        case 34:
          yy.cherryPick($$[$0 - 2], "", $$[$0 - 4], $$[$0]);
          break;
        case 35:
          yy.cherryPick($$[$0 - 2], "", "", $$[$0]);
          break;
        case 36:
          yy.merge($$[$0], "", "", "");
          break;
        case 37:
          yy.merge($$[$0 - 2], $$[$0], "", "");
          break;
        case 38:
          yy.merge($$[$0 - 2], "", $$[$0], "");
          break;
        case 39:
          yy.merge($$[$0 - 2], "", "", $$[$0]);
          break;
        case 40:
          yy.merge($$[$0 - 4], $$[$0], "", $$[$0 - 2]);
          break;
        case 41:
          yy.merge($$[$0 - 4], "", $$[$0], $$[$0 - 2]);
          break;
        case 42:
          yy.merge($$[$0 - 4], "", $$[$0 - 2], $$[$0]);
          break;
        case 43:
          yy.merge($$[$0 - 4], $$[$0 - 2], $$[$0], "");
          break;
        case 44:
          yy.merge($$[$0 - 4], $$[$0 - 2], "", $$[$0]);
          break;
        case 45:
          yy.merge($$[$0 - 4], $$[$0], $$[$0 - 2], "");
          break;
        case 46:
          yy.merge($$[$0 - 6], $$[$0 - 4], $$[$0 - 2], $$[$0]);
          break;
        case 47:
          yy.merge($$[$0 - 6], $$[$0], $$[$0 - 4], $$[$0 - 2]);
          break;
        case 48:
          yy.merge($$[$0 - 6], $$[$0 - 4], $$[$0], $$[$0 - 2]);
          break;
        case 49:
          yy.merge($$[$0 - 6], $$[$0 - 2], $$[$0 - 4], $$[$0]);
          break;
        case 50:
          yy.merge($$[$0 - 6], $$[$0], $$[$0 - 2], $$[$0 - 4]);
          break;
        case 51:
          yy.merge($$[$0 - 6], $$[$0 - 2], $$[$0], $$[$0 - 4]);
          break;
        case 52:
          yy.commit($$[$0]);
          break;
        case 53:
          yy.commit("", "", yy.commitType.NORMAL, $$[$0]);
          break;
        case 54:
          yy.commit("", "", $$[$0], "");
          break;
        case 55:
          yy.commit("", "", $$[$0], $$[$0 - 2]);
          break;
        case 56:
          yy.commit("", "", $$[$0 - 2], $$[$0]);
          break;
        case 57:
          yy.commit("", $$[$0], yy.commitType.NORMAL, "");
          break;
        case 58:
          yy.commit("", $$[$0 - 2], yy.commitType.NORMAL, $$[$0]);
          break;
        case 59:
          yy.commit("", $$[$0], yy.commitType.NORMAL, $$[$0 - 2]);
          break;
        case 60:
          yy.commit("", $$[$0 - 2], $$[$0], "");
          break;
        case 61:
          yy.commit("", $$[$0], $$[$0 - 2], "");
          break;
        case 62:
          yy.commit("", $$[$0 - 4], $$[$0 - 2], $$[$0]);
          break;
        case 63:
          yy.commit("", $$[$0 - 4], $$[$0], $$[$0 - 2]);
          break;
        case 64:
          yy.commit("", $$[$0 - 2], $$[$0 - 4], $$[$0]);
          break;
        case 65:
          yy.commit("", $$[$0], $$[$0 - 4], $$[$0 - 2]);
          break;
        case 66:
          yy.commit("", $$[$0], $$[$0 - 2], $$[$0 - 4]);
          break;
        case 67:
          yy.commit("", $$[$0 - 2], $$[$0], $$[$0 - 4]);
          break;
        case 68:
          yy.commit($$[$0], "", yy.commitType.NORMAL, "");
          break;
        case 69:
          yy.commit($$[$0], "", yy.commitType.NORMAL, $$[$0 - 2]);
          break;
        case 70:
          yy.commit($$[$0 - 2], "", yy.commitType.NORMAL, $$[$0]);
          break;
        case 71:
          yy.commit($$[$0 - 2], "", $$[$0], "");
          break;
        case 72:
          yy.commit($$[$0], "", $$[$0 - 2], "");
          break;
        case 73:
          yy.commit($$[$0], $$[$0 - 2], yy.commitType.NORMAL, "");
          break;
        case 74:
          yy.commit($$[$0 - 2], $$[$0], yy.commitType.NORMAL, "");
          break;
        case 75:
          yy.commit($$[$0 - 4], "", $$[$0 - 2], $$[$0]);
          break;
        case 76:
          yy.commit($$[$0 - 4], "", $$[$0], $$[$0 - 2]);
          break;
        case 77:
          yy.commit($$[$0 - 2], "", $$[$0 - 4], $$[$0]);
          break;
        case 78:
          yy.commit($$[$0], "", $$[$0 - 4], $$[$0 - 2]);
          break;
        case 79:
          yy.commit($$[$0], "", $$[$0 - 2], $$[$0 - 4]);
          break;
        case 80:
          yy.commit($$[$0 - 2], "", $$[$0], $$[$0 - 4]);
          break;
        case 81:
          yy.commit($$[$0 - 4], $$[$0], $$[$0 - 2], "");
          break;
        case 82:
          yy.commit($$[$0 - 4], $$[$0 - 2], $$[$0], "");
          break;
        case 83:
          yy.commit($$[$0 - 2], $$[$0], $$[$0 - 4], "");
          break;
        case 84:
          yy.commit($$[$0], $$[$0 - 2], $$[$0 - 4], "");
          break;
        case 85:
          yy.commit($$[$0], $$[$0 - 4], $$[$0 - 2], "");
          break;
        case 86:
          yy.commit($$[$0 - 2], $$[$0 - 4], $$[$0], "");
          break;
        case 87:
          yy.commit($$[$0 - 4], $$[$0], yy.commitType.NORMAL, $$[$0 - 2]);
          break;
        case 88:
          yy.commit($$[$0 - 4], $$[$0 - 2], yy.commitType.NORMAL, $$[$0]);
          break;
        case 89:
          yy.commit($$[$0 - 2], $$[$0], yy.commitType.NORMAL, $$[$0 - 4]);
          break;
        case 90:
          yy.commit($$[$0], $$[$0 - 2], yy.commitType.NORMAL, $$[$0 - 4]);
          break;
        case 91:
          yy.commit($$[$0], $$[$0 - 4], yy.commitType.NORMAL, $$[$0 - 2]);
          break;
        case 92:
          yy.commit($$[$0 - 2], $$[$0 - 4], yy.commitType.NORMAL, $$[$0]);
          break;
        case 93:
          yy.commit($$[$0 - 6], $$[$0 - 4], $$[$0 - 2], $$[$0]);
          break;
        case 94:
          yy.commit($$[$0 - 6], $$[$0 - 4], $$[$0], $$[$0 - 2]);
          break;
        case 95:
          yy.commit($$[$0 - 6], $$[$0 - 2], $$[$0 - 4], $$[$0]);
          break;
        case 96:
          yy.commit($$[$0 - 6], $$[$0], $$[$0 - 4], $$[$0 - 2]);
          break;
        case 97:
          yy.commit($$[$0 - 6], $$[$0 - 2], $$[$0], $$[$0 - 4]);
          break;
        case 98:
          yy.commit($$[$0 - 6], $$[$0], $$[$0 - 2], $$[$0 - 4]);
          break;
        case 99:
          yy.commit($$[$0 - 4], $$[$0 - 6], $$[$0 - 2], $$[$0]);
          break;
        case 100:
          yy.commit($$[$0 - 4], $$[$0 - 6], $$[$0], $$[$0 - 2]);
          break;
        case 101:
          yy.commit($$[$0 - 2], $$[$0 - 6], $$[$0 - 4], $$[$0]);
          break;
        case 102:
          yy.commit($$[$0], $$[$0 - 6], $$[$0 - 4], $$[$0 - 2]);
          break;
        case 103:
          yy.commit($$[$0 - 2], $$[$0 - 6], $$[$0], $$[$0 - 4]);
          break;
        case 104:
          yy.commit($$[$0], $$[$0 - 6], $$[$0 - 2], $$[$0 - 4]);
          break;
        case 105:
          yy.commit($$[$0], $$[$0 - 4], $$[$0 - 2], $$[$0 - 6]);
          break;
        case 106:
          yy.commit($$[$0 - 2], $$[$0 - 4], $$[$0], $$[$0 - 6]);
          break;
        case 107:
          yy.commit($$[$0], $$[$0 - 2], $$[$0 - 4], $$[$0 - 6]);
          break;
        case 108:
          yy.commit($$[$0 - 2], $$[$0], $$[$0 - 4], $$[$0 - 6]);
          break;
        case 109:
          yy.commit($$[$0 - 4], $$[$0 - 2], $$[$0], $$[$0 - 6]);
          break;
        case 110:
          yy.commit($$[$0 - 4], $$[$0], $$[$0 - 2], $$[$0 - 6]);
          break;
        case 111:
          yy.commit($$[$0 - 2], $$[$0 - 4], $$[$0 - 6], $$[$0]);
          break;
        case 112:
          yy.commit($$[$0], $$[$0 - 4], $$[$0 - 6], $$[$0 - 2]);
          break;
        case 113:
          yy.commit($$[$0 - 2], $$[$0], $$[$0 - 6], $$[$0 - 4]);
          break;
        case 114:
          yy.commit($$[$0], $$[$0 - 2], $$[$0 - 6], $$[$0 - 4]);
          break;
        case 115:
          yy.commit($$[$0 - 4], $$[$0 - 2], $$[$0 - 6], $$[$0]);
          break;
        case 116:
          yy.commit($$[$0 - 4], $$[$0], $$[$0 - 6], $$[$0 - 2]);
          break;
        case 117:
          this.$ = "";
          break;
        case 118:
          this.$ = $$[$0];
          break;
        case 119:
          this.$ = yy.commitType.NORMAL;
          break;
        case 120:
          this.$ = yy.commitType.REVERSE;
          break;
        case 121:
          this.$ = yy.commitType.HIGHLIGHT;
          break;
      }
    },
    table: [{ 3: 1, 4: 2, 5: $V0, 7: $V1, 13: $V2, 47: $V3 }, { 1: [3] }, { 3: 7, 4: 2, 5: $V0, 7: $V1, 13: $V2, 47: $V3 }, { 6: 8, 7: $V4, 8: [1, 9], 9: [1, 10], 10: 11, 13: $V5 }, o($V6, [2, 124]), o($V6, [2, 125]), o($V6, [2, 126]), { 1: [2, 1] }, { 7: [1, 13] }, { 6: 14, 7: $V4, 10: 11, 13: $V5 }, { 8: [1, 15] }, o($V7, [2, 9], { 11: 16, 12: [1, 17] }), o($V8, [2, 8]), { 1: [2, 2] }, { 7: [1, 18] }, { 6: 19, 7: $V4, 10: 11, 13: $V5 }, { 7: [2, 6], 13: [1, 22], 14: 20, 15: 21, 16: 23, 17: 24, 18: 25, 19: [1, 26], 21: [1, 27], 23: [1, 28], 24: [1, 29], 25: 30, 26: [1, 31], 28: [1, 35], 31: [1, 34], 37: [1, 33], 40: [1, 32] }, o($V8, [2, 7]), { 1: [2, 3] }, { 7: [1, 36] }, o($V7, [2, 10]), { 4: 37, 7: $V1, 13: $V2, 47: $V3 }, o($V7, [2, 12]), o($V9, [2, 13]), o($V9, [2, 14]), o($V9, [2, 15]), { 20: [1, 38] }, { 22: [1, 39] }, o($V9, [2, 18]), o($V9, [2, 19]), o($V9, [2, 20]), { 27: 40, 33: $Va, 46: $Vb }, o($V9, [2, 117], { 41: 43, 32: [1, 46], 33: [1, 48], 35: [1, 44], 38: [1, 45], 42: [1, 47] }), { 27: 49, 33: $Va, 46: $Vb }, { 32: [1, 50], 35: [1, 51] }, { 27: 52, 33: $Va, 46: $Vb }, { 1: [2, 4] }, o($V7, [2, 11]), o($V9, [2, 16]), o($V9, [2, 17]), o($V9, [2, 21]), o($Vc, [2, 122]), o($Vc, [2, 123]), o($V9, [2, 52]), { 33: [1, 53] }, { 39: 54, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 58] }, { 33: [1, 59] }, o($V9, [2, 118]), o($V9, [2, 36], { 32: [1, 60], 35: [1, 62], 38: [1, 61] }), { 33: [1, 63] }, { 33: [1, 64], 36: [1, 65] }, o($V9, [2, 22], { 29: [1, 66] }), o($V9, [2, 53], { 32: [1, 68], 38: [1, 67], 42: [1, 69] }), o($V9, [2, 54], { 32: [1, 71], 35: [1, 70], 42: [1, 72] }), o($Vg, [2, 119]), o($Vg, [2, 120]), o($Vg, [2, 121]), o($V9, [2, 57], { 35: [1, 73], 38: [1, 74], 42: [1, 75] }), o($V9, [2, 68], { 32: [1, 78], 35: [1, 76], 38: [1, 77] }), { 33: [1, 79] }, { 39: 80, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 81] }, o($V9, [2, 24], { 34: [1, 82], 35: [1, 83] }), { 32: [1, 84] }, { 32: [1, 85] }, { 30: [1, 86] }, { 39: 87, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 88] }, { 33: [1, 89] }, { 33: [1, 90] }, { 33: [1, 91] }, { 33: [1, 92] }, { 33: [1, 93] }, { 39: 94, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 95] }, { 33: [1, 96] }, { 39: 97, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 98] }, o($V9, [2, 37], { 35: [1, 100], 38: [1, 99] }), o($V9, [2, 38], { 32: [1, 102], 35: [1, 101] }), o($V9, [2, 39], { 32: [1, 103], 38: [1, 104] }), { 33: [1, 105] }, { 33: [1, 106], 36: [1, 107] }, { 33: [1, 108] }, { 33: [1, 109] }, o($V9, [2, 23]), o($V9, [2, 55], { 32: [1, 110], 42: [1, 111] }), o($V9, [2, 59], { 38: [1, 112], 42: [1, 113] }), o($V9, [2, 69], { 32: [1, 115], 38: [1, 114] }), o($V9, [2, 56], { 32: [1, 116], 42: [1, 117] }), o($V9, [2, 61], { 35: [1, 118], 42: [1, 119] }), o($V9, [2, 72], { 32: [1, 121], 35: [1, 120] }), o($V9, [2, 58], { 38: [1, 122], 42: [1, 123] }), o($V9, [2, 60], { 35: [1, 124], 42: [1, 125] }), o($V9, [2, 73], { 35: [1, 127], 38: [1, 126] }), o($V9, [2, 70], { 32: [1, 129], 38: [1, 128] }), o($V9, [2, 71], { 32: [1, 131], 35: [1, 130] }), o($V9, [2, 74], { 35: [1, 133], 38: [1, 132] }), { 39: 134, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 135] }, { 33: [1, 136] }, { 33: [1, 137] }, { 33: [1, 138] }, { 39: 139, 43: $Vd, 44: $Ve, 45: $Vf }, o($V9, [2, 25], { 35: [1, 140] }), o($V9, [2, 26], { 34: [1, 141] }), o($V9, [2, 31], { 34: [1, 142] }), o($V9, [2, 29], { 34: [1, 143] }), o($V9, [2, 30], { 34: [1, 144] }), { 33: [1, 145] }, { 33: [1, 146] }, { 39: 147, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 148] }, { 39: 149, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 150] }, { 33: [1, 151] }, { 33: [1, 152] }, { 33: [1, 153] }, { 33: [1, 154] }, { 33: [1, 155] }, { 33: [1, 156] }, { 39: 157, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 158] }, { 33: [1, 159] }, { 33: [1, 160] }, { 39: 161, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 162] }, { 39: 163, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 164] }, { 33: [1, 165] }, { 33: [1, 166] }, { 39: 167, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 168] }, o($V9, [2, 43], { 35: [1, 169] }), o($V9, [2, 44], { 38: [1, 170] }), o($V9, [2, 42], { 32: [1, 171] }), o($V9, [2, 45], { 35: [1, 172] }), o($V9, [2, 40], { 38: [1, 173] }), o($V9, [2, 41], { 32: [1, 174] }), { 33: [1, 175], 36: [1, 176] }, { 33: [1, 177] }, { 33: [1, 178] }, { 33: [1, 179] }, { 33: [1, 180] }, o($V9, [2, 66], { 42: [1, 181] }), o($V9, [2, 79], { 32: [1, 182] }), o($V9, [2, 67], { 42: [1, 183] }), o($V9, [2, 90], { 38: [1, 184] }), o($V9, [2, 80], { 32: [1, 185] }), o($V9, [2, 89], { 38: [1, 186] }), o($V9, [2, 65], { 42: [1, 187] }), o($V9, [2, 78], { 32: [1, 188] }), o($V9, [2, 64], { 42: [1, 189] }), o($V9, [2, 84], { 35: [1, 190] }), o($V9, [2, 77], { 32: [1, 191] }), o($V9, [2, 83], { 35: [1, 192] }), o($V9, [2, 63], { 42: [1, 193] }), o($V9, [2, 91], { 38: [1, 194] }), o($V9, [2, 62], { 42: [1, 195] }), o($V9, [2, 85], { 35: [1, 196] }), o($V9, [2, 86], { 35: [1, 197] }), o($V9, [2, 92], { 38: [1, 198] }), o($V9, [2, 76], { 32: [1, 199] }), o($V9, [2, 87], { 38: [1, 200] }), o($V9, [2, 75], { 32: [1, 201] }), o($V9, [2, 81], { 35: [1, 202] }), o($V9, [2, 82], { 35: [1, 203] }), o($V9, [2, 88], { 38: [1, 204] }), { 33: [1, 205] }, { 39: 206, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 207] }, { 33: [1, 208] }, { 39: 209, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 210] }, o($V9, [2, 27]), o($V9, [2, 32]), o($V9, [2, 28]), o($V9, [2, 33]), o($V9, [2, 34]), o($V9, [2, 35]), { 33: [1, 211] }, { 33: [1, 212] }, { 33: [1, 213] }, { 39: 214, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 215] }, { 39: 216, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 217] }, { 33: [1, 218] }, { 33: [1, 219] }, { 33: [1, 220] }, { 33: [1, 221] }, { 33: [1, 222] }, { 33: [1, 223] }, { 39: 224, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 225] }, { 33: [1, 226] }, { 33: [1, 227] }, { 39: 228, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 229] }, { 39: 230, 43: $Vd, 44: $Ve, 45: $Vf }, { 33: [1, 231] }, { 33: [1, 232] }, { 33: [1, 233] }, { 39: 234, 43: $Vd, 44: $Ve, 45: $Vf }, o($V9, [2, 46]), o($V9, [2, 48]), o($V9, [2, 47]), o($V9, [2, 49]), o($V9, [2, 51]), o($V9, [2, 50]), o($V9, [2, 107]), o($V9, [2, 108]), o($V9, [2, 105]), o($V9, [2, 106]), o($V9, [2, 110]), o($V9, [2, 109]), o($V9, [2, 114]), o($V9, [2, 113]), o($V9, [2, 112]), o($V9, [2, 111]), o($V9, [2, 116]), o($V9, [2, 115]), o($V9, [2, 104]), o($V9, [2, 103]), o($V9, [2, 102]), o($V9, [2, 101]), o($V9, [2, 99]), o($V9, [2, 100]), o($V9, [2, 98]), o($V9, [2, 97]), o($V9, [2, 96]), o($V9, [2, 95]), o($V9, [2, 93]), o($V9, [2, 94])],
    defaultActions: { 7: [2, 1], 13: [2, 2], 18: [2, 3], 36: [2, 4] },
    parseError: function parseError(str, hash) {
      if (hash.recoverable) {
        this.trace(str);
      } else {
        var error = new Error(str);
        error.hash = hash;
        throw error;
      }
    },
    parse: function parse(input) {
      var self = this, stack = [0], tstack = [], vstack = [null], lstack = [], table = this.table, yytext = "", yylineno = 0, yyleng = 0, TERROR = 2, EOF = 1;
      var args = lstack.slice.call(arguments, 1);
      var lexer2 = Object.create(this.lexer);
      var sharedState = { yy: {} };
      for (var k in this.yy) {
        if (Object.prototype.hasOwnProperty.call(this.yy, k)) {
          sharedState.yy[k] = this.yy[k];
        }
      }
      lexer2.setInput(input, sharedState.yy);
      sharedState.yy.lexer = lexer2;
      sharedState.yy.parser = this;
      if (typeof lexer2.yylloc == "undefined") {
        lexer2.yylloc = {};
      }
      var yyloc = lexer2.yylloc;
      lstack.push(yyloc);
      var ranges = lexer2.options && lexer2.options.ranges;
      if (typeof sharedState.yy.parseError === "function") {
        this.parseError = sharedState.yy.parseError;
      } else {
        this.parseError = Object.getPrototypeOf(this).parseError;
      }
      function lex() {
        var token;
        token = tstack.pop() || lexer2.lex() || EOF;
        if (typeof token !== "number") {
          if (token instanceof Array) {
            tstack = token;
            token = tstack.pop();
          }
          token = self.symbols_[token] || token;
        }
        return token;
      }
      var symbol, state, action, r, yyval = {}, p, len, newState, expected;
      while (true) {
        state = stack[stack.length - 1];
        if (this.defaultActions[state]) {
          action = this.defaultActions[state];
        } else {
          if (symbol === null || typeof symbol == "undefined") {
            symbol = lex();
          }
          action = table[state] && table[state][symbol];
        }
        if (typeof action === "undefined" || !action.length || !action[0]) {
          var errStr = "";
          expected = [];
          for (p in table[state]) {
            if (this.terminals_[p] && p > TERROR) {
              expected.push("'" + this.terminals_[p] + "'");
            }
          }
          if (lexer2.showPosition) {
            errStr = "Parse error on line " + (yylineno + 1) + ":\n" + lexer2.showPosition() + "\nExpecting " + expected.join(", ") + ", got '" + (this.terminals_[symbol] || symbol) + "'";
          } else {
            errStr = "Parse error on line " + (yylineno + 1) + ": Unexpected " + (symbol == EOF ? "end of input" : "'" + (this.terminals_[symbol] || symbol) + "'");
          }
          this.parseError(errStr, {
            text: lexer2.match,
            token: this.terminals_[symbol] || symbol,
            line: lexer2.yylineno,
            loc: yyloc,
            expected
          });
        }
        if (action[0] instanceof Array && action.length > 1) {
          throw new Error("Parse Error: multiple actions possible at state: " + state + ", token: " + symbol);
        }
        switch (action[0]) {
          case 1:
            stack.push(symbol);
            vstack.push(lexer2.yytext);
            lstack.push(lexer2.yylloc);
            stack.push(action[1]);
            symbol = null;
            {
              yyleng = lexer2.yyleng;
              yytext = lexer2.yytext;
              yylineno = lexer2.yylineno;
              yyloc = lexer2.yylloc;
            }
            break;
          case 2:
            len = this.productions_[action[1]][1];
            yyval.$ = vstack[vstack.length - len];
            yyval._$ = {
              first_line: lstack[lstack.length - (len || 1)].first_line,
              last_line: lstack[lstack.length - 1].last_line,
              first_column: lstack[lstack.length - (len || 1)].first_column,
              last_column: lstack[lstack.length - 1].last_column
            };
            if (ranges) {
              yyval._$.range = [
                lstack[lstack.length - (len || 1)].range[0],
                lstack[lstack.length - 1].range[1]
              ];
            }
            r = this.performAction.apply(yyval, [
              yytext,
              yyleng,
              yylineno,
              sharedState.yy,
              action[1],
              vstack,
              lstack
            ].concat(args));
            if (typeof r !== "undefined") {
              return r;
            }
            if (len) {
              stack = stack.slice(0, -1 * len * 2);
              vstack = vstack.slice(0, -1 * len);
              lstack = lstack.slice(0, -1 * len);
            }
            stack.push(this.productions_[action[1]][0]);
            vstack.push(yyval.$);
            lstack.push(yyval._$);
            newState = table[stack[stack.length - 2]][stack[stack.length - 1]];
            stack.push(newState);
            break;
          case 3:
            return true;
        }
      }
      return true;
    }
  };
  var lexer = function() {
    var lexer2 = {
      EOF: 1,
      parseError: function parseError(str, hash) {
        if (this.yy.parser) {
          this.yy.parser.parseError(str, hash);
        } else {
          throw new Error(str);
        }
      },
      // resets the lexer, sets new input
      setInput: function(input, yy) {
        this.yy = yy || this.yy || {};
        this._input = input;
        this._more = this._backtrack = this.done = false;
        this.yylineno = this.yyleng = 0;
        this.yytext = this.matched = this.match = "";
        this.conditionStack = ["INITIAL"];
        this.yylloc = {
          first_line: 1,
          first_column: 0,
          last_line: 1,
          last_column: 0
        };
        if (this.options.ranges) {
          this.yylloc.range = [0, 0];
        }
        this.offset = 0;
        return this;
      },
      // consumes and returns one char from the input
      input: function() {
        var ch = this._input[0];
        this.yytext += ch;
        this.yyleng++;
        this.offset++;
        this.match += ch;
        this.matched += ch;
        var lines = ch.match(/(?:\r\n?|\n).*/g);
        if (lines) {
          this.yylineno++;
          this.yylloc.last_line++;
        } else {
          this.yylloc.last_column++;
        }
        if (this.options.ranges) {
          this.yylloc.range[1]++;
        }
        this._input = this._input.slice(1);
        return ch;
      },
      // unshifts one char (or a string) into the input
      unput: function(ch) {
        var len = ch.length;
        var lines = ch.split(/(?:\r\n?|\n)/g);
        this._input = ch + this._input;
        this.yytext = this.yytext.substr(0, this.yytext.length - len);
        this.offset -= len;
        var oldLines = this.match.split(/(?:\r\n?|\n)/g);
        this.match = this.match.substr(0, this.match.length - 1);
        this.matched = this.matched.substr(0, this.matched.length - 1);
        if (lines.length - 1) {
          this.yylineno -= lines.length - 1;
        }
        var r = this.yylloc.range;
        this.yylloc = {
          first_line: this.yylloc.first_line,
          last_line: this.yylineno + 1,
          first_column: this.yylloc.first_column,
          last_column: lines ? (lines.length === oldLines.length ? this.yylloc.first_column : 0) + oldLines[oldLines.length - lines.length].length - lines[0].length : this.yylloc.first_column - len
        };
        if (this.options.ranges) {
          this.yylloc.range = [r[0], r[0] + this.yyleng - len];
        }
        this.yyleng = this.yytext.length;
        return this;
      },
      // When called from action, caches matched text and appends it on next action
      more: function() {
        this._more = true;
        return this;
      },
      // When called from action, signals the lexer that this rule fails to match the input, so the next matching rule (regex) should be tested instead.
      reject: function() {
        if (this.options.backtrack_lexer) {
          this._backtrack = true;
        } else {
          return this.parseError("Lexical error on line " + (this.yylineno + 1) + ". You can only invoke reject() in the lexer when the lexer is of the backtracking persuasion (options.backtrack_lexer = true).\n" + this.showPosition(), {
            text: "",
            token: null,
            line: this.yylineno
          });
        }
        return this;
      },
      // retain first n characters of the match
      less: function(n) {
        this.unput(this.match.slice(n));
      },
      // displays already matched input, i.e. for error messages
      pastInput: function() {
        var past = this.matched.substr(0, this.matched.length - this.match.length);
        return (past.length > 20 ? "..." : "") + past.substr(-20).replace(/\n/g, "");
      },
      // displays upcoming input, i.e. for error messages
      upcomingInput: function() {
        var next = this.match;
        if (next.length < 20) {
          next += this._input.substr(0, 20 - next.length);
        }
        return (next.substr(0, 20) + (next.length > 20 ? "..." : "")).replace(/\n/g, "");
      },
      // displays the character position where the lexing error occurred, i.e. for error messages
      showPosition: function() {
        var pre = this.pastInput();
        var c = new Array(pre.length + 1).join("-");
        return pre + this.upcomingInput() + "\n" + c + "^";
      },
      // test the lexed token: return FALSE when not a match, otherwise return token
      test_match: function(match, indexed_rule) {
        var token, lines, backup;
        if (this.options.backtrack_lexer) {
          backup = {
            yylineno: this.yylineno,
            yylloc: {
              first_line: this.yylloc.first_line,
              last_line: this.last_line,
              first_column: this.yylloc.first_column,
              last_column: this.yylloc.last_column
            },
            yytext: this.yytext,
            match: this.match,
            matches: this.matches,
            matched: this.matched,
            yyleng: this.yyleng,
            offset: this.offset,
            _more: this._more,
            _input: this._input,
            yy: this.yy,
            conditionStack: this.conditionStack.slice(0),
            done: this.done
          };
          if (this.options.ranges) {
            backup.yylloc.range = this.yylloc.range.slice(0);
          }
        }
        lines = match[0].match(/(?:\r\n?|\n).*/g);
        if (lines) {
          this.yylineno += lines.length;
        }
        this.yylloc = {
          first_line: this.yylloc.last_line,
          last_line: this.yylineno + 1,
          first_column: this.yylloc.last_column,
          last_column: lines ? lines[lines.length - 1].length - lines[lines.length - 1].match(/\r?\n?/)[0].length : this.yylloc.last_column + match[0].length
        };
        this.yytext += match[0];
        this.match += match[0];
        this.matches = match;
        this.yyleng = this.yytext.length;
        if (this.options.ranges) {
          this.yylloc.range = [this.offset, this.offset += this.yyleng];
        }
        this._more = false;
        this._backtrack = false;
        this._input = this._input.slice(match[0].length);
        this.matched += match[0];
        token = this.performAction.call(this, this.yy, this, indexed_rule, this.conditionStack[this.conditionStack.length - 1]);
        if (this.done && this._input) {
          this.done = false;
        }
        if (token) {
          return token;
        } else if (this._backtrack) {
          for (var k in backup) {
            this[k] = backup[k];
          }
          return false;
        }
        return false;
      },
      // return next match in input
      next: function() {
        if (this.done) {
          return this.EOF;
        }
        if (!this._input) {
          this.done = true;
        }
        var token, match, tempMatch, index;
        if (!this._more) {
          this.yytext = "";
          this.match = "";
        }
        var rules = this._currentRules();
        for (var i = 0; i < rules.length; i++) {
          tempMatch = this._input.match(this.rules[rules[i]]);
          if (tempMatch && (!match || tempMatch[0].length > match[0].length)) {
            match = tempMatch;
            index = i;
            if (this.options.backtrack_lexer) {
              token = this.test_match(tempMatch, rules[i]);
              if (token !== false) {
                return token;
              } else if (this._backtrack) {
                match = false;
                continue;
              } else {
                return false;
              }
            } else if (!this.options.flex) {
              break;
            }
          }
        }
        if (match) {
          token = this.test_match(match, rules[index]);
          if (token !== false) {
            return token;
          }
          return false;
        }
        if (this._input === "") {
          return this.EOF;
        } else {
          return this.parseError("Lexical error on line " + (this.yylineno + 1) + ". Unrecognized text.\n" + this.showPosition(), {
            text: "",
            token: null,
            line: this.yylineno
          });
        }
      },
      // return next match that has a token
      lex: function lex() {
        var r = this.next();
        if (r) {
          return r;
        } else {
          return this.lex();
        }
      },
      // activates a new lexer condition state (pushes the new lexer condition state onto the condition stack)
      begin: function begin(condition) {
        this.conditionStack.push(condition);
      },
      // pop the previously active lexer condition state off the condition stack
      popState: function popState() {
        var n = this.conditionStack.length - 1;
        if (n > 0) {
          return this.conditionStack.pop();
        } else {
          return this.conditionStack[0];
        }
      },
      // produce the lexer rule set which is active for the currently active lexer condition state
      _currentRules: function _currentRules() {
        if (this.conditionStack.length && this.conditionStack[this.conditionStack.length - 1]) {
          return this.conditions[this.conditionStack[this.conditionStack.length - 1]].rules;
        } else {
          return this.conditions["INITIAL"].rules;
        }
      },
      // return the currently active lexer condition state; when an index argument is provided it produces the N-th previous condition state, if available
      topState: function topState(n) {
        n = this.conditionStack.length - 1 - Math.abs(n || 0);
        if (n >= 0) {
          return this.conditionStack[n];
        } else {
          return "INITIAL";
        }
      },
      // alias for begin(condition)
      pushState: function pushState(condition) {
        this.begin(condition);
      },
      // return the number of states currently on the stack
      stateStackSize: function stateStackSize() {
        return this.conditionStack.length;
      },
      options: { "case-insensitive": true },
      performAction: function anonymous(yy, yy_, $avoiding_name_collisions, YY_START) {
        switch ($avoiding_name_collisions) {
          case 0:
            this.begin("acc_title");
            return 19;
          case 1:
            this.popState();
            return "acc_title_value";
          case 2:
            this.begin("acc_descr");
            return 21;
          case 3:
            this.popState();
            return "acc_descr_value";
          case 4:
            this.begin("acc_descr_multiline");
            break;
          case 5:
            this.popState();
            break;
          case 6:
            return "acc_descr_multiline_value";
          case 7:
            return 13;
          case 8:
            break;
          case 9:
            break;
          case 10:
            return 5;
          case 11:
            return 40;
          case 12:
            return 32;
          case 13:
            return 38;
          case 14:
            return 42;
          case 15:
            return 43;
          case 16:
            return 44;
          case 17:
            return 45;
          case 18:
            return 35;
          case 19:
            return 28;
          case 20:
            return 29;
          case 21:
            return 37;
          case 22:
            return 31;
          case 23:
            return 34;
          case 24:
            return 26;
          case 25:
            return 9;
          case 26:
            return 9;
          case 27:
            return 8;
          case 28:
            return "CARET";
          case 29:
            this.begin("options");
            break;
          case 30:
            this.popState();
            break;
          case 31:
            return 12;
          case 32:
            return 36;
          case 33:
            this.begin("string");
            break;
          case 34:
            this.popState();
            break;
          case 35:
            return 33;
          case 36:
            return 30;
          case 37:
            return 46;
          case 38:
            return 7;
        }
      },
      rules: [/^(?:accTitle\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*\{\s*)/i, /^(?:[\}])/i, /^(?:[^\}]*)/i, /^(?:(\r?\n)+)/i, /^(?:#[^\n]*)/i, /^(?:%[^\n]*)/i, /^(?:gitGraph\b)/i, /^(?:commit(?=\s|$))/i, /^(?:id:)/i, /^(?:type:)/i, /^(?:msg:)/i, /^(?:NORMAL\b)/i, /^(?:REVERSE\b)/i, /^(?:HIGHLIGHT\b)/i, /^(?:tag:)/i, /^(?:branch(?=\s|$))/i, /^(?:order:)/i, /^(?:merge(?=\s|$))/i, /^(?:cherry-pick(?=\s|$))/i, /^(?:parent:)/i, /^(?:checkout(?=\s|$))/i, /^(?:LR\b)/i, /^(?:TB\b)/i, /^(?::)/i, /^(?:\^)/i, /^(?:options\r?\n)/i, /^(?:[ \r\n\t]+end\b)/i, /^(?:[\s\S]+(?=[ \r\n\t]+end))/i, /^(?:["]["])/i, /^(?:["])/i, /^(?:["])/i, /^(?:[^"]*)/i, /^(?:[0-9]+(?=\s|$))/i, /^(?:\w([-\./\w]*[-\w])?)/i, /^(?:$)/i, /^(?:\s+)/i],
      conditions: { "acc_descr_multiline": { "rules": [5, 6], "inclusive": false }, "acc_descr": { "rules": [3], "inclusive": false }, "acc_title": { "rules": [1], "inclusive": false }, "options": { "rules": [30, 31], "inclusive": false }, "string": { "rules": [34, 35], "inclusive": false }, "INITIAL": { "rules": [0, 2, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 32, 33, 36, 37, 38, 39], "inclusive": true } }
    };
    return lexer2;
  }();
  parser2.lexer = lexer;
  function Parser() {
    this.yy = {};
  }
  Parser.prototype = parser2;
  parser2.Parser = Parser;
  return new Parser();
}();
parser.parser = parser;
const gitGraphParser = parser;
let mainBranchName = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().gitGraph.mainBranchName;
let mainBranchOrder = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().gitGraph.mainBranchOrder;
let commits = {};
let head = null;
let branchesConfig = {};
branchesConfig[mainBranchName] = { name: mainBranchName, order: mainBranchOrder };
let branches = {};
branches[mainBranchName] = head;
let curBranch = mainBranchName;
let direction = "LR";
let seq = 0;
function getId() {
  return (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.x)({ length: 7 });
}
function uniqBy(list, fn) {
  const recordMap = /* @__PURE__ */ Object.create(null);
  return list.reduce((out, item) => {
    const key = fn(item);
    if (!recordMap[key]) {
      recordMap[key] = true;
      out.push(item);
    }
    return out;
  }, []);
}
const setDirection = function(dir2) {
  direction = dir2;
};
let options = {};
const setOptions = function(rawOptString) {
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.debug("options str", rawOptString);
  rawOptString = rawOptString && rawOptString.trim();
  rawOptString = rawOptString || "{}";
  try {
    options = JSON.parse(rawOptString);
  } catch (e) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.error("error while parsing gitGraph options", e.message);
  }
};
const getOptions = function() {
  return options;
};
const commit = function(msg, id, type, tag) {
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.debug("Entering commit:", msg, id, type, tag);
  id = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.e.sanitizeText(id, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)());
  msg = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.e.sanitizeText(msg, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)());
  tag = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.e.sanitizeText(tag, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)());
  const commit2 = {
    id: id ? id : seq + "-" + getId(),
    message: msg,
    seq: seq++,
    type: type ? type : commitType$1.NORMAL,
    tag: tag ? tag : "",
    parents: head == null ? [] : [head.id],
    branch: curBranch
  };
  head = commit2;
  commits[commit2.id] = commit2;
  branches[curBranch] = commit2.id;
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.debug("in pushCommit " + commit2.id);
};
const branch = function(name, order) {
  name = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.e.sanitizeText(name, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)());
  if (branches[name] === void 0) {
    branches[name] = head != null ? head.id : null;
    branchesConfig[name] = { name, order: order ? parseInt(order, 10) : null };
    checkout(name);
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.debug("in createBranch");
  } else {
    let error = new Error(
      'Trying to create an existing branch. (Help: Either use a new name if you want create a new branch or try using "checkout ' + name + '")'
    );
    error.hash = {
      text: "branch " + name,
      token: "branch " + name,
      line: "1",
      loc: { first_line: 1, last_line: 1, first_column: 1, last_column: 1 },
      expected: ['"checkout ' + name + '"']
    };
    throw error;
  }
};
const merge = function(otherBranch, custom_id, override_type, custom_tag) {
  otherBranch = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.e.sanitizeText(otherBranch, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)());
  custom_id = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.e.sanitizeText(custom_id, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)());
  const currentCommit = commits[branches[curBranch]];
  const otherCommit = commits[branches[otherBranch]];
  if (curBranch === otherBranch) {
    let error = new Error('Incorrect usage of "merge". Cannot merge a branch to itself');
    error.hash = {
      text: "merge " + otherBranch,
      token: "merge " + otherBranch,
      line: "1",
      loc: { first_line: 1, last_line: 1, first_column: 1, last_column: 1 },
      expected: ["branch abc"]
    };
    throw error;
  } else if (currentCommit === void 0 || !currentCommit) {
    let error = new Error(
      'Incorrect usage of "merge". Current branch (' + curBranch + ")has no commits"
    );
    error.hash = {
      text: "merge " + otherBranch,
      token: "merge " + otherBranch,
      line: "1",
      loc: { first_line: 1, last_line: 1, first_column: 1, last_column: 1 },
      expected: ["commit"]
    };
    throw error;
  } else if (branches[otherBranch] === void 0) {
    let error = new Error(
      'Incorrect usage of "merge". Branch to be merged (' + otherBranch + ") does not exist"
    );
    error.hash = {
      text: "merge " + otherBranch,
      token: "merge " + otherBranch,
      line: "1",
      loc: { first_line: 1, last_line: 1, first_column: 1, last_column: 1 },
      expected: ["branch " + otherBranch]
    };
    throw error;
  } else if (otherCommit === void 0 || !otherCommit) {
    let error = new Error(
      'Incorrect usage of "merge". Branch to be merged (' + otherBranch + ") has no commits"
    );
    error.hash = {
      text: "merge " + otherBranch,
      token: "merge " + otherBranch,
      line: "1",
      loc: { first_line: 1, last_line: 1, first_column: 1, last_column: 1 },
      expected: ['"commit"']
    };
    throw error;
  } else if (currentCommit === otherCommit) {
    let error = new Error('Incorrect usage of "merge". Both branches have same head');
    error.hash = {
      text: "merge " + otherBranch,
      token: "merge " + otherBranch,
      line: "1",
      loc: { first_line: 1, last_line: 1, first_column: 1, last_column: 1 },
      expected: ["branch abc"]
    };
    throw error;
  } else if (custom_id && commits[custom_id] !== void 0) {
    let error = new Error(
      'Incorrect usage of "merge". Commit with id:' + custom_id + " already exists, use different custom Id"
    );
    error.hash = {
      text: "merge " + otherBranch + custom_id + override_type + custom_tag,
      token: "merge " + otherBranch + custom_id + override_type + custom_tag,
      line: "1",
      loc: { first_line: 1, last_line: 1, first_column: 1, last_column: 1 },
      expected: [
        "merge " + otherBranch + " " + custom_id + "_UNIQUE " + override_type + " " + custom_tag
      ]
    };
    throw error;
  }
  const commit2 = {
    id: custom_id ? custom_id : seq + "-" + getId(),
    message: "merged branch " + otherBranch + " into " + curBranch,
    seq: seq++,
    parents: [head == null ? null : head.id, branches[otherBranch]],
    branch: curBranch,
    type: commitType$1.MERGE,
    customType: override_type,
    customId: custom_id ? true : false,
    tag: custom_tag ? custom_tag : ""
  };
  head = commit2;
  commits[commit2.id] = commit2;
  branches[curBranch] = commit2.id;
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.debug(branches);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.debug("in mergeBranch");
};
const cherryPick = function(sourceId, targetId, tag, parentCommitId) {
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.debug("Entering cherryPick:", sourceId, targetId, tag);
  sourceId = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.e.sanitizeText(sourceId, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)());
  targetId = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.e.sanitizeText(targetId, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)());
  tag = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.e.sanitizeText(tag, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)());
  parentCommitId = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.e.sanitizeText(parentCommitId, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)());
  if (!sourceId || commits[sourceId] === void 0) {
    let error = new Error(
      'Incorrect usage of "cherryPick". Source commit id should exist and provided'
    );
    error.hash = {
      text: "cherryPick " + sourceId + " " + targetId,
      token: "cherryPick " + sourceId + " " + targetId,
      line: "1",
      loc: { first_line: 1, last_line: 1, first_column: 1, last_column: 1 },
      expected: ["cherry-pick abc"]
    };
    throw error;
  }
  let sourceCommit = commits[sourceId];
  let sourceCommitBranch = sourceCommit.branch;
  if (parentCommitId && !(Array.isArray(sourceCommit.parents) && sourceCommit.parents.includes(parentCommitId))) {
    let error = new Error(
      "Invalid operation: The specified parent commit is not an immediate parent of the cherry-picked commit."
    );
    throw error;
  }
  if (sourceCommit.type === commitType$1.MERGE && !parentCommitId) {
    let error = new Error(
      "Incorrect usage of cherry-pick: If the source commit is a merge commit, an immediate parent commit must be specified."
    );
    throw error;
  }
  if (!targetId || commits[targetId] === void 0) {
    if (sourceCommitBranch === curBranch) {
      let error = new Error(
        'Incorrect usage of "cherryPick". Source commit is already on current branch'
      );
      error.hash = {
        text: "cherryPick " + sourceId + " " + targetId,
        token: "cherryPick " + sourceId + " " + targetId,
        line: "1",
        loc: { first_line: 1, last_line: 1, first_column: 1, last_column: 1 },
        expected: ["cherry-pick abc"]
      };
      throw error;
    }
    const currentCommit = commits[branches[curBranch]];
    if (currentCommit === void 0 || !currentCommit) {
      let error = new Error(
        'Incorrect usage of "cherry-pick". Current branch (' + curBranch + ")has no commits"
      );
      error.hash = {
        text: "cherryPick " + sourceId + " " + targetId,
        token: "cherryPick " + sourceId + " " + targetId,
        line: "1",
        loc: { first_line: 1, last_line: 1, first_column: 1, last_column: 1 },
        expected: ["cherry-pick abc"]
      };
      throw error;
    }
    const commit2 = {
      id: seq + "-" + getId(),
      message: "cherry-picked " + sourceCommit + " into " + curBranch,
      seq: seq++,
      parents: [head == null ? null : head.id, sourceCommit.id],
      branch: curBranch,
      type: commitType$1.CHERRY_PICK,
      tag: tag ?? `cherry-pick:${sourceCommit.id}${sourceCommit.type === commitType$1.MERGE ? `|parent:${parentCommitId}` : ""}`
    };
    head = commit2;
    commits[commit2.id] = commit2;
    branches[curBranch] = commit2.id;
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.debug(branches);
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.debug("in cherryPick");
  }
};
const checkout = function(branch2) {
  branch2 = _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.e.sanitizeText(branch2, (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)());
  if (branches[branch2] === void 0) {
    let error = new Error(
      'Trying to checkout branch which is not yet created. (Help try using "branch ' + branch2 + '")'
    );
    error.hash = {
      text: "checkout " + branch2,
      token: "checkout " + branch2,
      line: "1",
      loc: { first_line: 1, last_line: 1, first_column: 1, last_column: 1 },
      expected: ['"branch ' + branch2 + '"']
    };
    throw error;
  } else {
    curBranch = branch2;
    const id = branches[curBranch];
    head = commits[id];
  }
};
function upsert(arr, key, newVal) {
  const index = arr.indexOf(key);
  if (index === -1) {
    arr.push(newVal);
  } else {
    arr.splice(index, 1, newVal);
  }
}
function prettyPrintCommitHistory(commitArr) {
  const commit2 = commitArr.reduce((out, commit3) => {
    if (out.seq > commit3.seq) {
      return out;
    }
    return commit3;
  }, commitArr[0]);
  let line = "";
  commitArr.forEach(function(c) {
    if (c === commit2) {
      line += "	*";
    } else {
      line += "	|";
    }
  });
  const label = [line, commit2.id, commit2.seq];
  for (let branch2 in branches) {
    if (branches[branch2] === commit2.id) {
      label.push(branch2);
    }
  }
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.debug(label.join(" "));
  if (commit2.parents && commit2.parents.length == 2) {
    const newCommit = commits[commit2.parents[0]];
    upsert(commitArr, commit2, newCommit);
    commitArr.push(commits[commit2.parents[1]]);
  } else if (commit2.parents.length == 0) {
    return;
  } else {
    const nextCommit = commits[commit2.parents];
    upsert(commitArr, commit2, nextCommit);
  }
  commitArr = uniqBy(commitArr, (c) => c.id);
  prettyPrintCommitHistory(commitArr);
}
const prettyPrint = function() {
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.debug(commits);
  const node = getCommitsArray()[0];
  prettyPrintCommitHistory([node]);
};
const clear$1 = function() {
  commits = {};
  head = null;
  let mainBranch = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().gitGraph.mainBranchName;
  let mainBranchOrder2 = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().gitGraph.mainBranchOrder;
  branches = {};
  branches[mainBranch] = null;
  branchesConfig = {};
  branchesConfig[mainBranch] = { name: mainBranch, order: mainBranchOrder2 };
  curBranch = mainBranch;
  seq = 0;
  (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.t)();
};
const getBranchesAsObjArray = function() {
  const branchesArray = Object.values(branchesConfig).map((branchConfig, i) => {
    if (branchConfig.order !== null) {
      return branchConfig;
    }
    return {
      ...branchConfig,
      order: parseFloat(`0.${i}`, 10)
    };
  }).sort((a, b) => a.order - b.order).map(({ name }) => ({ name }));
  return branchesArray;
};
const getBranches = function() {
  return branches;
};
const getCommits = function() {
  return commits;
};
const getCommitsArray = function() {
  const commitArr = Object.keys(commits).map(function(key) {
    return commits[key];
  });
  commitArr.forEach(function(o) {
    _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.debug(o.id);
  });
  commitArr.sort((a, b) => a.seq - b.seq);
  return commitArr;
};
const getCurrentBranch = function() {
  return curBranch;
};
const getDirection = function() {
  return direction;
};
const getHead = function() {
  return head;
};
const commitType$1 = {
  NORMAL: 0,
  REVERSE: 1,
  HIGHLIGHT: 2,
  MERGE: 3,
  CHERRY_PICK: 4
};
const gitGraphDb = {
  getConfig: () => (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().gitGraph,
  setDirection,
  setOptions,
  getOptions,
  commit,
  branch,
  merge,
  cherryPick,
  checkout,
  //reset,
  prettyPrint,
  clear: clear$1,
  getBranchesAsObjArray,
  getBranches,
  getCommits,
  getCommitsArray,
  getCurrentBranch,
  getDirection,
  getHead,
  setAccTitle: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.s,
  getAccTitle: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.g,
  getAccDescription: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.a,
  setAccDescription: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.b,
  setDiagramTitle: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.q,
  getDiagramTitle: _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.r,
  commitType: commitType$1
};
let allCommitsDict = {};
const commitType = {
  NORMAL: 0,
  REVERSE: 1,
  HIGHLIGHT: 2,
  MERGE: 3,
  CHERRY_PICK: 4
};
const THEME_COLOR_LIMIT = 8;
let branchPos = {};
let commitPos = {};
let lanes = [];
let maxPos = 0;
let dir = "LR";
const clear = () => {
  branchPos = {};
  commitPos = {};
  allCommitsDict = {};
  maxPos = 0;
  lanes = [];
  dir = "LR";
};
const drawText = (txt) => {
  const svgLabel = document.createElementNS("http://www.w3.org/2000/svg", "text");
  let rows = [];
  if (typeof txt === "string") {
    rows = txt.split(/\\n|\n|<br\s*\/?>/gi);
  } else if (Array.isArray(txt)) {
    rows = txt;
  } else {
    rows = [];
  }
  for (const row of rows) {
    const tspan = document.createElementNS("http://www.w3.org/2000/svg", "tspan");
    tspan.setAttributeNS("http://www.w3.org/XML/1998/namespace", "xml:space", "preserve");
    tspan.setAttribute("dy", "1em");
    tspan.setAttribute("x", "0");
    tspan.setAttribute("class", "row");
    tspan.textContent = row.trim();
    svgLabel.appendChild(tspan);
  }
  return svgLabel;
};
const drawCommits = (svg, commits2, modifyGraph) => {
  const gitGraphConfig = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().gitGraph;
  const gBullets = svg.append("g").attr("class", "commit-bullets");
  const gLabels = svg.append("g").attr("class", "commit-labels");
  let pos = 0;
  if (dir === "TB") {
    pos = 30;
  }
  const keys = Object.keys(commits2);
  const sortedKeys = keys.sort((a, b) => {
    return commits2[a].seq - commits2[b].seq;
  });
  sortedKeys.forEach((key) => {
    const commit2 = commits2[key];
    const y = dir === "TB" ? pos + 10 : branchPos[commit2.branch].pos;
    const x = dir === "TB" ? branchPos[commit2.branch].pos : pos + 10;
    if (modifyGraph) {
      let typeClass;
      let commitSymbolType = commit2.customType !== void 0 && commit2.customType !== "" ? commit2.customType : commit2.type;
      switch (commitSymbolType) {
        case commitType.NORMAL:
          typeClass = "commit-normal";
          break;
        case commitType.REVERSE:
          typeClass = "commit-reverse";
          break;
        case commitType.HIGHLIGHT:
          typeClass = "commit-highlight";
          break;
        case commitType.MERGE:
          typeClass = "commit-merge";
          break;
        case commitType.CHERRY_PICK:
          typeClass = "commit-cherry-pick";
          break;
        default:
          typeClass = "commit-normal";
      }
      if (commitSymbolType === commitType.HIGHLIGHT) {
        const circle = gBullets.append("rect");
        circle.attr("x", x - 10);
        circle.attr("y", y - 10);
        circle.attr("height", 20);
        circle.attr("width", 20);
        circle.attr(
          "class",
          `commit ${commit2.id} commit-highlight${branchPos[commit2.branch].index % THEME_COLOR_LIMIT} ${typeClass}-outer`
        );
        gBullets.append("rect").attr("x", x - 6).attr("y", y - 6).attr("height", 12).attr("width", 12).attr(
          "class",
          `commit ${commit2.id} commit${branchPos[commit2.branch].index % THEME_COLOR_LIMIT} ${typeClass}-inner`
        );
      } else if (commitSymbolType === commitType.CHERRY_PICK) {
        gBullets.append("circle").attr("cx", x).attr("cy", y).attr("r", 10).attr("class", `commit ${commit2.id} ${typeClass}`);
        gBullets.append("circle").attr("cx", x - 3).attr("cy", y + 2).attr("r", 2.75).attr("fill", "#fff").attr("class", `commit ${commit2.id} ${typeClass}`);
        gBullets.append("circle").attr("cx", x + 3).attr("cy", y + 2).attr("r", 2.75).attr("fill", "#fff").attr("class", `commit ${commit2.id} ${typeClass}`);
        gBullets.append("line").attr("x1", x + 3).attr("y1", y + 1).attr("x2", x).attr("y2", y - 5).attr("stroke", "#fff").attr("class", `commit ${commit2.id} ${typeClass}`);
        gBullets.append("line").attr("x1", x - 3).attr("y1", y + 1).attr("x2", x).attr("y2", y - 5).attr("stroke", "#fff").attr("class", `commit ${commit2.id} ${typeClass}`);
      } else {
        const circle = gBullets.append("circle");
        circle.attr("cx", x);
        circle.attr("cy", y);
        circle.attr("r", commit2.type === commitType.MERGE ? 9 : 10);
        circle.attr(
          "class",
          `commit ${commit2.id} commit${branchPos[commit2.branch].index % THEME_COLOR_LIMIT}`
        );
        if (commitSymbolType === commitType.MERGE) {
          const circle2 = gBullets.append("circle");
          circle2.attr("cx", x);
          circle2.attr("cy", y);
          circle2.attr("r", 6);
          circle2.attr(
            "class",
            `commit ${typeClass} ${commit2.id} commit${branchPos[commit2.branch].index % THEME_COLOR_LIMIT}`
          );
        }
        if (commitSymbolType === commitType.REVERSE) {
          const cross = gBullets.append("path");
          cross.attr("d", `M ${x - 5},${y - 5}L${x + 5},${y + 5}M${x - 5},${y + 5}L${x + 5},${y - 5}`).attr(
            "class",
            `commit ${typeClass} ${commit2.id} commit${branchPos[commit2.branch].index % THEME_COLOR_LIMIT}`
          );
        }
      }
    }
    if (dir === "TB") {
      commitPos[commit2.id] = { x, y: pos + 10 };
    } else {
      commitPos[commit2.id] = { x: pos + 10, y };
    }
    if (modifyGraph) {
      const px = 4;
      const py = 2;
      if (commit2.type !== commitType.CHERRY_PICK && (commit2.customId && commit2.type === commitType.MERGE || commit2.type !== commitType.MERGE) && gitGraphConfig.showCommitLabel) {
        const wrapper = gLabels.append("g");
        const labelBkg = wrapper.insert("rect").attr("class", "commit-label-bkg");
        const text = wrapper.append("text").attr("x", pos).attr("y", y + 25).attr("class", "commit-label").text(commit2.id);
        let bbox = text.node().getBBox();
        labelBkg.attr("x", pos + 10 - bbox.width / 2 - py).attr("y", y + 13.5).attr("width", bbox.width + 2 * py).attr("height", bbox.height + 2 * py);
        if (dir === "TB") {
          labelBkg.attr("x", x - (bbox.width + 4 * px + 5)).attr("y", y - 12);
          text.attr("x", x - (bbox.width + 4 * px)).attr("y", y + bbox.height - 12);
        }
        if (dir !== "TB") {
          text.attr("x", pos + 10 - bbox.width / 2);
        }
        if (gitGraphConfig.rotateCommitLabel) {
          if (dir === "TB") {
            text.attr("transform", "rotate(-45, " + x + ", " + y + ")");
            labelBkg.attr("transform", "rotate(-45, " + x + ", " + y + ")");
          } else {
            let r_x = -7.5 - (bbox.width + 10) / 25 * 9.5;
            let r_y = 10 + bbox.width / 25 * 8.5;
            wrapper.attr(
              "transform",
              "translate(" + r_x + ", " + r_y + ") rotate(-45, " + pos + ", " + y + ")"
            );
          }
        }
      }
      if (commit2.tag) {
        const rect = gLabels.insert("polygon");
        const hole = gLabels.append("circle");
        const tag = gLabels.append("text").attr("y", y - 16).attr("class", "tag-label").text(commit2.tag);
        let tagBbox = tag.node().getBBox();
        tag.attr("x", pos + 10 - tagBbox.width / 2);
        const h2 = tagBbox.height / 2;
        const ly = y - 19.2;
        rect.attr("class", "tag-label-bkg").attr(
          "points",
          `
          ${pos - tagBbox.width / 2 - px / 2},${ly + py}
          ${pos - tagBbox.width / 2 - px / 2},${ly - py}
          ${pos + 10 - tagBbox.width / 2 - px},${ly - h2 - py}
          ${pos + 10 + tagBbox.width / 2 + px},${ly - h2 - py}
          ${pos + 10 + tagBbox.width / 2 + px},${ly + h2 + py}
          ${pos + 10 - tagBbox.width / 2 - px},${ly + h2 + py}`
        );
        hole.attr("cx", pos - tagBbox.width / 2 + px / 2).attr("cy", ly).attr("r", 1.5).attr("class", "tag-hole");
        if (dir === "TB") {
          rect.attr("class", "tag-label-bkg").attr(
            "points",
            `
            ${x},${pos + py}
            ${x},${pos - py}
            ${x + 10},${pos - h2 - py}
            ${x + 10 + tagBbox.width + px},${pos - h2 - py}
            ${x + 10 + tagBbox.width + px},${pos + h2 + py}
            ${x + 10},${pos + h2 + py}`
          ).attr("transform", "translate(12,12) rotate(45, " + x + "," + pos + ")");
          hole.attr("cx", x + px / 2).attr("cy", pos).attr("transform", "translate(12,12) rotate(45, " + x + "," + pos + ")");
          tag.attr("x", x + 5).attr("y", pos + 3).attr("transform", "translate(14,14) rotate(45, " + x + "," + pos + ")");
        }
      }
    }
    pos += 50;
    if (pos > maxPos) {
      maxPos = pos;
    }
  });
};
const shouldRerouteArrow = (commitA, commitB, p1, p2, allCommits) => {
  const commitBIsFurthest = dir === "TB" ? p1.x < p2.x : p1.y < p2.y;
  const branchToGetCurve = commitBIsFurthest ? commitB.branch : commitA.branch;
  const isOnBranchToGetCurve = (x) => x.branch === branchToGetCurve;
  const isBetweenCommits = (x) => x.seq > commitA.seq && x.seq < commitB.seq;
  return Object.values(allCommits).some((commitX) => {
    return isBetweenCommits(commitX) && isOnBranchToGetCurve(commitX);
  });
};
const findLane = (y1, y2, depth = 0) => {
  const candidate = y1 + Math.abs(y1 - y2) / 2;
  if (depth > 5) {
    return candidate;
  }
  let ok = lanes.every((lane) => Math.abs(lane - candidate) >= 10);
  if (ok) {
    lanes.push(candidate);
    return candidate;
  }
  const diff = Math.abs(y1 - y2);
  return findLane(y1, y2 - diff / 5, depth + 1);
};
const drawArrow = (svg, commitA, commitB, allCommits) => {
  const p1 = commitPos[commitA.id];
  const p2 = commitPos[commitB.id];
  const arrowNeedsRerouting = shouldRerouteArrow(commitA, commitB, p1, p2, allCommits);
  let arc = "";
  let arc2 = "";
  let radius = 0;
  let offset = 0;
  let colorClassNum = branchPos[commitB.branch].index;
  let lineDef;
  if (arrowNeedsRerouting) {
    arc = "A 10 10, 0, 0, 0,";
    arc2 = "A 10 10, 0, 0, 1,";
    radius = 10;
    offset = 10;
    const lineY = p1.y < p2.y ? findLane(p1.y, p2.y) : findLane(p2.y, p1.y);
    const lineX = p1.x < p2.x ? findLane(p1.x, p2.x) : findLane(p2.x, p1.x);
    if (dir === "TB") {
      if (p1.x < p2.x) {
        colorClassNum = branchPos[commitB.branch].index;
        lineDef = `M ${p1.x} ${p1.y} L ${lineX - radius} ${p1.y} ${arc2} ${lineX} ${p1.y + offset} L ${lineX} ${p2.y - radius} ${arc} ${lineX + offset} ${p2.y} L ${p2.x} ${p2.y}`;
      } else {
        colorClassNum = branchPos[commitA.branch].index;
        lineDef = `M ${p1.x} ${p1.y} L ${lineX + radius} ${p1.y} ${arc} ${lineX} ${p1.y + offset} L ${lineX} ${p2.y - radius} ${arc2} ${lineX - offset} ${p2.y} L ${p2.x} ${p2.y}`;
      }
    } else {
      if (p1.y < p2.y) {
        colorClassNum = branchPos[commitB.branch].index;
        lineDef = `M ${p1.x} ${p1.y} L ${p1.x} ${lineY - radius} ${arc} ${p1.x + offset} ${lineY} L ${p2.x - radius} ${lineY} ${arc2} ${p2.x} ${lineY + offset} L ${p2.x} ${p2.y}`;
      } else {
        colorClassNum = branchPos[commitA.branch].index;
        lineDef = `M ${p1.x} ${p1.y} L ${p1.x} ${lineY + radius} ${arc2} ${p1.x + offset} ${lineY} L ${p2.x - radius} ${lineY} ${arc} ${p2.x} ${lineY - offset} L ${p2.x} ${p2.y}`;
      }
    }
  } else {
    if (dir === "TB") {
      if (p1.x < p2.x) {
        arc = "A 20 20, 0, 0, 0,";
        arc2 = "A 20 20, 0, 0, 1,";
        radius = 20;
        offset = 20;
        colorClassNum = branchPos[commitB.branch].index;
        lineDef = `M ${p1.x} ${p1.y} L ${p2.x - radius} ${p1.y} ${arc2} ${p2.x} ${p1.y + offset} L ${p2.x} ${p2.y}`;
      }
      if (p1.x > p2.x) {
        arc = "A 20 20, 0, 0, 0,";
        arc2 = "A 20 20, 0, 0, 1,";
        radius = 20;
        offset = 20;
        colorClassNum = branchPos[commitA.branch].index;
        lineDef = `M ${p1.x} ${p1.y} L ${p1.x} ${p2.y - radius} ${arc2} ${p1.x - offset} ${p2.y} L ${p2.x} ${p2.y}`;
      }
      if (p1.x === p2.x) {
        colorClassNum = branchPos[commitA.branch].index;
        lineDef = `M ${p1.x} ${p1.y} L ${p1.x + radius} ${p1.y} ${arc} ${p1.x + offset} ${p2.y + radius} L ${p2.x} ${p2.y}`;
      }
    } else {
      if (p1.y < p2.y) {
        arc = "A 20 20, 0, 0, 0,";
        radius = 20;
        offset = 20;
        colorClassNum = branchPos[commitB.branch].index;
        lineDef = `M ${p1.x} ${p1.y} L ${p1.x} ${p2.y - radius} ${arc} ${p1.x + offset} ${p2.y} L ${p2.x} ${p2.y}`;
      }
      if (p1.y > p2.y) {
        arc = "A 20 20, 0, 0, 0,";
        radius = 20;
        offset = 20;
        colorClassNum = branchPos[commitA.branch].index;
        lineDef = `M ${p1.x} ${p1.y} L ${p2.x - radius} ${p1.y} ${arc} ${p2.x} ${p1.y - offset} L ${p2.x} ${p2.y}`;
      }
      if (p1.y === p2.y) {
        colorClassNum = branchPos[commitA.branch].index;
        lineDef = `M ${p1.x} ${p1.y} L ${p1.x} ${p2.y - radius} ${arc} ${p1.x + offset} ${p2.y} L ${p2.x} ${p2.y}`;
      }
    }
  }
  svg.append("path").attr("d", lineDef).attr("class", "arrow arrow" + colorClassNum % THEME_COLOR_LIMIT);
};
const drawArrows = (svg, commits2) => {
  const gArrows = svg.append("g").attr("class", "commit-arrows");
  Object.keys(commits2).forEach((key) => {
    const commit2 = commits2[key];
    if (commit2.parents && commit2.parents.length > 0) {
      commit2.parents.forEach((parent) => {
        drawArrow(gArrows, commits2[parent], commit2, commits2);
      });
    }
  });
};
const drawBranches = (svg, branches2) => {
  const gitGraphConfig = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)().gitGraph;
  const g = svg.append("g");
  branches2.forEach((branch2, index) => {
    const adjustIndexForTheme = index % THEME_COLOR_LIMIT;
    const pos = branchPos[branch2.name].pos;
    const line = g.append("line");
    line.attr("x1", 0);
    line.attr("y1", pos);
    line.attr("x2", maxPos);
    line.attr("y2", pos);
    line.attr("class", "branch branch" + adjustIndexForTheme);
    if (dir === "TB") {
      line.attr("y1", 30);
      line.attr("x1", pos);
      line.attr("y2", maxPos);
      line.attr("x2", pos);
    }
    lanes.push(pos);
    let name = branch2.name;
    const labelElement = drawText(name);
    const bkg = g.insert("rect");
    const branchLabel = g.insert("g").attr("class", "branchLabel");
    const label = branchLabel.insert("g").attr("class", "label branch-label" + adjustIndexForTheme);
    label.node().appendChild(labelElement);
    let bbox = labelElement.getBBox();
    bkg.attr("class", "branchLabelBkg label" + adjustIndexForTheme).attr("rx", 4).attr("ry", 4).attr("x", -bbox.width - 4 - (gitGraphConfig.rotateCommitLabel === true ? 30 : 0)).attr("y", -bbox.height / 2 + 8).attr("width", bbox.width + 18).attr("height", bbox.height + 4);
    label.attr(
      "transform",
      "translate(" + (-bbox.width - 14 - (gitGraphConfig.rotateCommitLabel === true ? 30 : 0)) + ", " + (pos - bbox.height / 2 - 1) + ")"
    );
    if (dir === "TB") {
      bkg.attr("x", pos - bbox.width / 2 - 10).attr("y", 0);
      label.attr("transform", "translate(" + (pos - bbox.width / 2 - 5) + ", 0)");
    }
    if (dir !== "TB") {
      bkg.attr("transform", "translate(-19, " + (pos - bbox.height / 2) + ")");
    }
  });
};
const draw = function(txt, id, ver, diagObj) {
  clear();
  const conf = (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.c)();
  const gitGraphConfig = conf.gitGraph;
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.l.debug("in gitgraph renderer", txt + "\n", "id:", id, ver);
  allCommitsDict = diagObj.db.getCommits();
  const branches2 = diagObj.db.getBranchesAsObjArray();
  dir = diagObj.db.getDirection();
  const diagram2 = (0,d3__WEBPACK_IMPORTED_MODULE_0__/* .select */ .Ys)(`[id="${id}"]`);
  let pos = 0;
  branches2.forEach((branch2, index) => {
    const labelElement = drawText(branch2.name);
    const g = diagram2.append("g");
    const branchLabel = g.insert("g").attr("class", "branchLabel");
    const label = branchLabel.insert("g").attr("class", "label branch-label");
    label.node().appendChild(labelElement);
    let bbox = labelElement.getBBox();
    branchPos[branch2.name] = { pos, index };
    pos += 50 + (gitGraphConfig.rotateCommitLabel ? 40 : 0) + (dir === "TB" ? bbox.width / 2 : 0);
    label.remove();
    branchLabel.remove();
    g.remove();
  });
  drawCommits(diagram2, allCommitsDict, false);
  if (gitGraphConfig.showBranches) {
    drawBranches(diagram2, branches2);
  }
  drawArrows(diagram2, allCommitsDict);
  drawCommits(diagram2, allCommitsDict, true);
  _mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.u.insertTitle(
    diagram2,
    "gitTitleText",
    gitGraphConfig.titleTopMargin,
    diagObj.db.getDiagramTitle()
  );
  (0,_mermaid_04fb0060_js__WEBPACK_IMPORTED_MODULE_4__.y)(
    void 0,
    diagram2,
    gitGraphConfig.diagramPadding,
    gitGraphConfig.useMaxWidth ?? conf.useMaxWidth
  );
};
const gitGraphRenderer = {
  draw
};
const getStyles = (options2) => `
  .commit-id,
  .commit-msg,
  .branch-label {
    fill: lightgrey;
    color: lightgrey;
    font-family: 'trebuchet ms', verdana, arial, sans-serif;
    font-family: var(--mermaid-font-family);
  }
  ${[0, 1, 2, 3, 4, 5, 6, 7].map(
  (i) => `
        .branch-label${i} { fill: ${options2["gitBranchLabel" + i]}; }
        .commit${i} { stroke: ${options2["git" + i]}; fill: ${options2["git" + i]}; }
        .commit-highlight${i} { stroke: ${options2["gitInv" + i]}; fill: ${options2["gitInv" + i]}; }
        .label${i}  { fill: ${options2["git" + i]}; }
        .arrow${i} { stroke: ${options2["git" + i]}; }
        `
).join("\n")}

  .branch {
    stroke-width: 1;
    stroke: ${options2.lineColor};
    stroke-dasharray: 2;
  }
  .commit-label { font-size: ${options2.commitLabelFontSize}; fill: ${options2.commitLabelColor};}
  .commit-label-bkg { font-size: ${options2.commitLabelFontSize}; fill: ${options2.commitLabelBackground}; opacity: 0.5; }
  .tag-label { font-size: ${options2.tagLabelFontSize}; fill: ${options2.tagLabelColor};}
  .tag-label-bkg { fill: ${options2.tagLabelBackground}; stroke: ${options2.tagLabelBorder}; }
  .tag-hole { fill: ${options2.textColor}; }

  .commit-merge {
    stroke: ${options2.primaryColor};
    fill: ${options2.primaryColor};
  }
  .commit-reverse {
    stroke: ${options2.primaryColor};
    fill: ${options2.primaryColor};
    stroke-width: 3;
  }
  .commit-highlight-outer {
  }
  .commit-highlight-inner {
    stroke: ${options2.primaryColor};
    fill: ${options2.primaryColor};
  }

  .arrow { stroke-width: 8; stroke-linecap: round; fill: none}
  .gitTitleText {
    text-anchor: middle;
    font-size: 18px;
    fill: ${options2.textColor};
  }
`;
const gitGraphStyles = getStyles;
const diagram = {
  parser: gitGraphParser,
  db: gitGraphDb,
  renderer: gitGraphRenderer,
  styles: gitGraphStyles
};



/***/ })

}]);
//# sourceMappingURL=261.7c7a6b6d904fd35115a3.js.map?v=7c7a6b6d904fd35115a3
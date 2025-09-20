"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[4886],{

/***/ 34886:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Ee: () => (/* binding */ styles_default),
/* harmony export */   P0: () => (/* binding */ classDiagram_default),
/* harmony export */   b0: () => (/* binding */ classRenderer_v3_unified_default),
/* harmony export */   dR: () => (/* binding */ ClassDB)
/* harmony export */ });
/* harmony import */ var _chunk_RZ5BOZE2_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(49479);
/* harmony import */ var _chunk_TYCBKAJE_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(83696);
/* harmony import */ var _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(42626);
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(86906);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(83619);





// src/diagrams/class/parser/classDiagram.jison
var parser = function() {
  var o = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(k, v, o2, l) {
    for (o2 = o2 || {}, l = k.length; l--; o2[k[l]] = v) ;
    return o2;
  }, "o"), $V0 = [1, 18], $V1 = [1, 19], $V2 = [1, 20], $V3 = [1, 41], $V4 = [1, 42], $V5 = [1, 26], $V6 = [1, 24], $V7 = [1, 25], $V8 = [1, 32], $V9 = [1, 33], $Va = [1, 34], $Vb = [1, 45], $Vc = [1, 35], $Vd = [1, 36], $Ve = [1, 37], $Vf = [1, 38], $Vg = [1, 27], $Vh = [1, 28], $Vi = [1, 29], $Vj = [1, 30], $Vk = [1, 31], $Vl = [1, 44], $Vm = [1, 46], $Vn = [1, 43], $Vo = [1, 47], $Vp = [1, 9], $Vq = [1, 8, 9], $Vr = [1, 58], $Vs = [1, 59], $Vt = [1, 60], $Vu = [1, 61], $Vv = [1, 62], $Vw = [1, 63], $Vx = [1, 64], $Vy = [1, 8, 9, 41], $Vz = [1, 76], $VA = [1, 8, 9, 12, 13, 22, 39, 41, 44, 66, 67, 68, 69, 70, 71, 72, 77, 79], $VB = [1, 8, 9, 12, 13, 17, 20, 22, 39, 41, 44, 48, 58, 66, 67, 68, 69, 70, 71, 72, 77, 79, 84, 99, 101, 102], $VC = [13, 58, 84, 99, 101, 102], $VD = [13, 58, 71, 72, 84, 99, 101, 102], $VE = [13, 58, 66, 67, 68, 69, 70, 84, 99, 101, 102], $VF = [1, 98], $VG = [1, 115], $VH = [1, 107], $VI = [1, 113], $VJ = [1, 108], $VK = [1, 109], $VL = [1, 110], $VM = [1, 111], $VN = [1, 112], $VO = [1, 114], $VP = [22, 58, 59, 80, 84, 85, 86, 87, 88, 89], $VQ = [1, 8, 9, 39, 41, 44], $VR = [1, 8, 9, 22], $VS = [1, 143], $VT = [1, 8, 9, 59], $VU = [1, 8, 9, 22, 58, 59, 80, 84, 85, 86, 87, 88, 89];
  var parser2 = {
    trace: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function trace() {
    }, "trace"),
    yy: {},
    symbols_: { "error": 2, "start": 3, "mermaidDoc": 4, "statements": 5, "graphConfig": 6, "CLASS_DIAGRAM": 7, "NEWLINE": 8, "EOF": 9, "statement": 10, "classLabel": 11, "SQS": 12, "STR": 13, "SQE": 14, "namespaceName": 15, "alphaNumToken": 16, "DOT": 17, "className": 18, "classLiteralName": 19, "GENERICTYPE": 20, "relationStatement": 21, "LABEL": 22, "namespaceStatement": 23, "classStatement": 24, "memberStatement": 25, "annotationStatement": 26, "clickStatement": 27, "styleStatement": 28, "cssClassStatement": 29, "noteStatement": 30, "classDefStatement": 31, "direction": 32, "acc_title": 33, "acc_title_value": 34, "acc_descr": 35, "acc_descr_value": 36, "acc_descr_multiline_value": 37, "namespaceIdentifier": 38, "STRUCT_START": 39, "classStatements": 40, "STRUCT_STOP": 41, "NAMESPACE": 42, "classIdentifier": 43, "STYLE_SEPARATOR": 44, "members": 45, "CLASS": 46, "ANNOTATION_START": 47, "ANNOTATION_END": 48, "MEMBER": 49, "SEPARATOR": 50, "relation": 51, "NOTE_FOR": 52, "noteText": 53, "NOTE": 54, "CLASSDEF": 55, "classList": 56, "stylesOpt": 57, "ALPHA": 58, "COMMA": 59, "direction_tb": 60, "direction_bt": 61, "direction_rl": 62, "direction_lr": 63, "relationType": 64, "lineType": 65, "AGGREGATION": 66, "EXTENSION": 67, "COMPOSITION": 68, "DEPENDENCY": 69, "LOLLIPOP": 70, "LINE": 71, "DOTTED_LINE": 72, "CALLBACK": 73, "LINK": 74, "LINK_TARGET": 75, "CLICK": 76, "CALLBACK_NAME": 77, "CALLBACK_ARGS": 78, "HREF": 79, "STYLE": 80, "CSSCLASS": 81, "style": 82, "styleComponent": 83, "NUM": 84, "COLON": 85, "UNIT": 86, "SPACE": 87, "BRKT": 88, "PCT": 89, "commentToken": 90, "textToken": 91, "graphCodeTokens": 92, "textNoTagsToken": 93, "TAGSTART": 94, "TAGEND": 95, "==": 96, "--": 97, "DEFAULT": 98, "MINUS": 99, "keywords": 100, "UNICODE_TEXT": 101, "BQUOTE_STR": 102, "$accept": 0, "$end": 1 },
    terminals_: { 2: "error", 7: "CLASS_DIAGRAM", 8: "NEWLINE", 9: "EOF", 12: "SQS", 13: "STR", 14: "SQE", 17: "DOT", 20: "GENERICTYPE", 22: "LABEL", 33: "acc_title", 34: "acc_title_value", 35: "acc_descr", 36: "acc_descr_value", 37: "acc_descr_multiline_value", 39: "STRUCT_START", 41: "STRUCT_STOP", 42: "NAMESPACE", 44: "STYLE_SEPARATOR", 46: "CLASS", 47: "ANNOTATION_START", 48: "ANNOTATION_END", 49: "MEMBER", 50: "SEPARATOR", 52: "NOTE_FOR", 54: "NOTE", 55: "CLASSDEF", 58: "ALPHA", 59: "COMMA", 60: "direction_tb", 61: "direction_bt", 62: "direction_rl", 63: "direction_lr", 66: "AGGREGATION", 67: "EXTENSION", 68: "COMPOSITION", 69: "DEPENDENCY", 70: "LOLLIPOP", 71: "LINE", 72: "DOTTED_LINE", 73: "CALLBACK", 74: "LINK", 75: "LINK_TARGET", 76: "CLICK", 77: "CALLBACK_NAME", 78: "CALLBACK_ARGS", 79: "HREF", 80: "STYLE", 81: "CSSCLASS", 84: "NUM", 85: "COLON", 86: "UNIT", 87: "SPACE", 88: "BRKT", 89: "PCT", 92: "graphCodeTokens", 94: "TAGSTART", 95: "TAGEND", 96: "==", 97: "--", 98: "DEFAULT", 99: "MINUS", 100: "keywords", 101: "UNICODE_TEXT", 102: "BQUOTE_STR" },
    productions_: [0, [3, 1], [3, 1], [4, 1], [6, 4], [5, 1], [5, 2], [5, 3], [11, 3], [15, 1], [15, 3], [15, 2], [18, 1], [18, 3], [18, 1], [18, 2], [18, 2], [18, 2], [10, 1], [10, 2], [10, 1], [10, 1], [10, 1], [10, 1], [10, 1], [10, 1], [10, 1], [10, 1], [10, 1], [10, 1], [10, 2], [10, 2], [10, 1], [23, 4], [23, 5], [38, 2], [40, 1], [40, 2], [40, 3], [24, 1], [24, 3], [24, 4], [24, 6], [43, 2], [43, 3], [26, 4], [45, 1], [45, 2], [25, 1], [25, 2], [25, 1], [25, 1], [21, 3], [21, 4], [21, 4], [21, 5], [30, 3], [30, 2], [31, 3], [56, 1], [56, 3], [32, 1], [32, 1], [32, 1], [32, 1], [51, 3], [51, 2], [51, 2], [51, 1], [64, 1], [64, 1], [64, 1], [64, 1], [64, 1], [65, 1], [65, 1], [27, 3], [27, 4], [27, 3], [27, 4], [27, 4], [27, 5], [27, 3], [27, 4], [27, 4], [27, 5], [27, 4], [27, 5], [27, 5], [27, 6], [28, 3], [29, 3], [57, 1], [57, 3], [82, 1], [82, 2], [83, 1], [83, 1], [83, 1], [83, 1], [83, 1], [83, 1], [83, 1], [83, 1], [83, 1], [90, 1], [90, 1], [91, 1], [91, 1], [91, 1], [91, 1], [91, 1], [91, 1], [91, 1], [93, 1], [93, 1], [93, 1], [93, 1], [16, 1], [16, 1], [16, 1], [16, 1], [19, 1], [53, 1]],
    performAction: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function anonymous(yytext, yyleng, yylineno, yy, yystate, $$, _$) {
      var $0 = $$.length - 1;
      switch (yystate) {
        case 8:
          this.$ = $$[$0 - 1];
          break;
        case 9:
        case 12:
        case 14:
          this.$ = $$[$0];
          break;
        case 10:
        case 13:
          this.$ = $$[$0 - 2] + "." + $$[$0];
          break;
        case 11:
        case 15:
          this.$ = $$[$0 - 1] + $$[$0];
          break;
        case 16:
        case 17:
          this.$ = $$[$0 - 1] + "~" + $$[$0] + "~";
          break;
        case 18:
          yy.addRelation($$[$0]);
          break;
        case 19:
          $$[$0 - 1].title = yy.cleanupLabel($$[$0]);
          yy.addRelation($$[$0 - 1]);
          break;
        case 30:
          this.$ = $$[$0].trim();
          yy.setAccTitle(this.$);
          break;
        case 31:
        case 32:
          this.$ = $$[$0].trim();
          yy.setAccDescription(this.$);
          break;
        case 33:
          yy.addClassesToNamespace($$[$0 - 3], $$[$0 - 1]);
          break;
        case 34:
          yy.addClassesToNamespace($$[$0 - 4], $$[$0 - 1]);
          break;
        case 35:
          this.$ = $$[$0];
          yy.addNamespace($$[$0]);
          break;
        case 36:
          this.$ = [$$[$0]];
          break;
        case 37:
          this.$ = [$$[$0 - 1]];
          break;
        case 38:
          $$[$0].unshift($$[$0 - 2]);
          this.$ = $$[$0];
          break;
        case 40:
          yy.setCssClass($$[$0 - 2], $$[$0]);
          break;
        case 41:
          yy.addMembers($$[$0 - 3], $$[$0 - 1]);
          break;
        case 42:
          yy.setCssClass($$[$0 - 5], $$[$0 - 3]);
          yy.addMembers($$[$0 - 5], $$[$0 - 1]);
          break;
        case 43:
          this.$ = $$[$0];
          yy.addClass($$[$0]);
          break;
        case 44:
          this.$ = $$[$0 - 1];
          yy.addClass($$[$0 - 1]);
          yy.setClassLabel($$[$0 - 1], $$[$0]);
          break;
        case 45:
          yy.addAnnotation($$[$0], $$[$0 - 2]);
          break;
        case 46:
        case 59:
          this.$ = [$$[$0]];
          break;
        case 47:
          $$[$0].push($$[$0 - 1]);
          this.$ = $$[$0];
          break;
        case 48:
          break;
        case 49:
          yy.addMember($$[$0 - 1], yy.cleanupLabel($$[$0]));
          break;
        case 50:
          break;
        case 51:
          break;
        case 52:
          this.$ = { "id1": $$[$0 - 2], "id2": $$[$0], relation: $$[$0 - 1], relationTitle1: "none", relationTitle2: "none" };
          break;
        case 53:
          this.$ = { id1: $$[$0 - 3], id2: $$[$0], relation: $$[$0 - 1], relationTitle1: $$[$0 - 2], relationTitle2: "none" };
          break;
        case 54:
          this.$ = { id1: $$[$0 - 3], id2: $$[$0], relation: $$[$0 - 2], relationTitle1: "none", relationTitle2: $$[$0 - 1] };
          break;
        case 55:
          this.$ = { id1: $$[$0 - 4], id2: $$[$0], relation: $$[$0 - 2], relationTitle1: $$[$0 - 3], relationTitle2: $$[$0 - 1] };
          break;
        case 56:
          yy.addNote($$[$0], $$[$0 - 1]);
          break;
        case 57:
          yy.addNote($$[$0]);
          break;
        case 58:
          this.$ = $$[$0 - 2];
          yy.defineClass($$[$0 - 1], $$[$0]);
          break;
        case 60:
          this.$ = $$[$0 - 2].concat([$$[$0]]);
          break;
        case 61:
          yy.setDirection("TB");
          break;
        case 62:
          yy.setDirection("BT");
          break;
        case 63:
          yy.setDirection("RL");
          break;
        case 64:
          yy.setDirection("LR");
          break;
        case 65:
          this.$ = { type1: $$[$0 - 2], type2: $$[$0], lineType: $$[$0 - 1] };
          break;
        case 66:
          this.$ = { type1: "none", type2: $$[$0], lineType: $$[$0 - 1] };
          break;
        case 67:
          this.$ = { type1: $$[$0 - 1], type2: "none", lineType: $$[$0] };
          break;
        case 68:
          this.$ = { type1: "none", type2: "none", lineType: $$[$0] };
          break;
        case 69:
          this.$ = yy.relationType.AGGREGATION;
          break;
        case 70:
          this.$ = yy.relationType.EXTENSION;
          break;
        case 71:
          this.$ = yy.relationType.COMPOSITION;
          break;
        case 72:
          this.$ = yy.relationType.DEPENDENCY;
          break;
        case 73:
          this.$ = yy.relationType.LOLLIPOP;
          break;
        case 74:
          this.$ = yy.lineType.LINE;
          break;
        case 75:
          this.$ = yy.lineType.DOTTED_LINE;
          break;
        case 76:
        case 82:
          this.$ = $$[$0 - 2];
          yy.setClickEvent($$[$0 - 1], $$[$0]);
          break;
        case 77:
        case 83:
          this.$ = $$[$0 - 3];
          yy.setClickEvent($$[$0 - 2], $$[$0 - 1]);
          yy.setTooltip($$[$0 - 2], $$[$0]);
          break;
        case 78:
          this.$ = $$[$0 - 2];
          yy.setLink($$[$0 - 1], $$[$0]);
          break;
        case 79:
          this.$ = $$[$0 - 3];
          yy.setLink($$[$0 - 2], $$[$0 - 1], $$[$0]);
          break;
        case 80:
          this.$ = $$[$0 - 3];
          yy.setLink($$[$0 - 2], $$[$0 - 1]);
          yy.setTooltip($$[$0 - 2], $$[$0]);
          break;
        case 81:
          this.$ = $$[$0 - 4];
          yy.setLink($$[$0 - 3], $$[$0 - 2], $$[$0]);
          yy.setTooltip($$[$0 - 3], $$[$0 - 1]);
          break;
        case 84:
          this.$ = $$[$0 - 3];
          yy.setClickEvent($$[$0 - 2], $$[$0 - 1], $$[$0]);
          break;
        case 85:
          this.$ = $$[$0 - 4];
          yy.setClickEvent($$[$0 - 3], $$[$0 - 2], $$[$0 - 1]);
          yy.setTooltip($$[$0 - 3], $$[$0]);
          break;
        case 86:
          this.$ = $$[$0 - 3];
          yy.setLink($$[$0 - 2], $$[$0]);
          break;
        case 87:
          this.$ = $$[$0 - 4];
          yy.setLink($$[$0 - 3], $$[$0 - 1], $$[$0]);
          break;
        case 88:
          this.$ = $$[$0 - 4];
          yy.setLink($$[$0 - 3], $$[$0 - 1]);
          yy.setTooltip($$[$0 - 3], $$[$0]);
          break;
        case 89:
          this.$ = $$[$0 - 5];
          yy.setLink($$[$0 - 4], $$[$0 - 2], $$[$0]);
          yy.setTooltip($$[$0 - 4], $$[$0 - 1]);
          break;
        case 90:
          this.$ = $$[$0 - 2];
          yy.setCssStyle($$[$0 - 1], $$[$0]);
          break;
        case 91:
          yy.setCssClass($$[$0 - 1], $$[$0]);
          break;
        case 92:
          this.$ = [$$[$0]];
          break;
        case 93:
          $$[$0 - 2].push($$[$0]);
          this.$ = $$[$0 - 2];
          break;
        case 95:
          this.$ = $$[$0 - 1] + $$[$0];
          break;
      }
    }, "anonymous"),
    table: [{ 3: 1, 4: 2, 5: 3, 6: 4, 7: [1, 6], 10: 5, 16: 39, 18: 21, 19: 40, 21: 7, 23: 8, 24: 9, 25: 10, 26: 11, 27: 12, 28: 13, 29: 14, 30: 15, 31: 16, 32: 17, 33: $V0, 35: $V1, 37: $V2, 38: 22, 42: $V3, 43: 23, 46: $V4, 47: $V5, 49: $V6, 50: $V7, 52: $V8, 54: $V9, 55: $Va, 58: $Vb, 60: $Vc, 61: $Vd, 62: $Ve, 63: $Vf, 73: $Vg, 74: $Vh, 76: $Vi, 80: $Vj, 81: $Vk, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }, { 1: [3] }, { 1: [2, 1] }, { 1: [2, 2] }, { 1: [2, 3] }, o($Vp, [2, 5], { 8: [1, 48] }), { 8: [1, 49] }, o($Vq, [2, 18], { 22: [1, 50] }), o($Vq, [2, 20]), o($Vq, [2, 21]), o($Vq, [2, 22]), o($Vq, [2, 23]), o($Vq, [2, 24]), o($Vq, [2, 25]), o($Vq, [2, 26]), o($Vq, [2, 27]), o($Vq, [2, 28]), o($Vq, [2, 29]), { 34: [1, 51] }, { 36: [1, 52] }, o($Vq, [2, 32]), o($Vq, [2, 48], { 51: 53, 64: 56, 65: 57, 13: [1, 54], 22: [1, 55], 66: $Vr, 67: $Vs, 68: $Vt, 69: $Vu, 70: $Vv, 71: $Vw, 72: $Vx }), { 39: [1, 65] }, o($Vy, [2, 39], { 39: [1, 67], 44: [1, 66] }), o($Vq, [2, 50]), o($Vq, [2, 51]), { 16: 68, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn }, { 16: 39, 18: 69, 19: 40, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }, { 16: 39, 18: 70, 19: 40, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }, { 16: 39, 18: 71, 19: 40, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }, { 58: [1, 72] }, { 13: [1, 73] }, { 16: 39, 18: 74, 19: 40, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }, { 13: $Vz, 53: 75 }, { 56: 77, 58: [1, 78] }, o($Vq, [2, 61]), o($Vq, [2, 62]), o($Vq, [2, 63]), o($Vq, [2, 64]), o($VA, [2, 12], { 16: 39, 19: 40, 18: 80, 17: [1, 79], 20: [1, 81], 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }), o($VA, [2, 14], { 20: [1, 82] }), { 15: 83, 16: 84, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn }, { 16: 39, 18: 85, 19: 40, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }, o($VB, [2, 118]), o($VB, [2, 119]), o($VB, [2, 120]), o($VB, [2, 121]), o([1, 8, 9, 12, 13, 20, 22, 39, 41, 44, 66, 67, 68, 69, 70, 71, 72, 77, 79], [2, 122]), o($Vp, [2, 6], { 10: 5, 21: 7, 23: 8, 24: 9, 25: 10, 26: 11, 27: 12, 28: 13, 29: 14, 30: 15, 31: 16, 32: 17, 18: 21, 38: 22, 43: 23, 16: 39, 19: 40, 5: 86, 33: $V0, 35: $V1, 37: $V2, 42: $V3, 46: $V4, 47: $V5, 49: $V6, 50: $V7, 52: $V8, 54: $V9, 55: $Va, 58: $Vb, 60: $Vc, 61: $Vd, 62: $Ve, 63: $Vf, 73: $Vg, 74: $Vh, 76: $Vi, 80: $Vj, 81: $Vk, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }), { 5: 87, 10: 5, 16: 39, 18: 21, 19: 40, 21: 7, 23: 8, 24: 9, 25: 10, 26: 11, 27: 12, 28: 13, 29: 14, 30: 15, 31: 16, 32: 17, 33: $V0, 35: $V1, 37: $V2, 38: 22, 42: $V3, 43: 23, 46: $V4, 47: $V5, 49: $V6, 50: $V7, 52: $V8, 54: $V9, 55: $Va, 58: $Vb, 60: $Vc, 61: $Vd, 62: $Ve, 63: $Vf, 73: $Vg, 74: $Vh, 76: $Vi, 80: $Vj, 81: $Vk, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }, o($Vq, [2, 19]), o($Vq, [2, 30]), o($Vq, [2, 31]), { 13: [1, 89], 16: 39, 18: 88, 19: 40, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }, { 51: 90, 64: 56, 65: 57, 66: $Vr, 67: $Vs, 68: $Vt, 69: $Vu, 70: $Vv, 71: $Vw, 72: $Vx }, o($Vq, [2, 49]), { 65: 91, 71: $Vw, 72: $Vx }, o($VC, [2, 68], { 64: 92, 66: $Vr, 67: $Vs, 68: $Vt, 69: $Vu, 70: $Vv }), o($VD, [2, 69]), o($VD, [2, 70]), o($VD, [2, 71]), o($VD, [2, 72]), o($VD, [2, 73]), o($VE, [2, 74]), o($VE, [2, 75]), { 8: [1, 94], 24: 95, 40: 93, 43: 23, 46: $V4 }, { 16: 96, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn }, { 45: 97, 49: $VF }, { 48: [1, 99] }, { 13: [1, 100] }, { 13: [1, 101] }, { 77: [1, 102], 79: [1, 103] }, { 22: $VG, 57: 104, 58: $VH, 80: $VI, 82: 105, 83: 106, 84: $VJ, 85: $VK, 86: $VL, 87: $VM, 88: $VN, 89: $VO }, { 58: [1, 116] }, { 13: $Vz, 53: 117 }, o($Vq, [2, 57]), o($Vq, [2, 123]), { 22: $VG, 57: 118, 58: $VH, 59: [1, 119], 80: $VI, 82: 105, 83: 106, 84: $VJ, 85: $VK, 86: $VL, 87: $VM, 88: $VN, 89: $VO }, o($VP, [2, 59]), { 16: 39, 18: 120, 19: 40, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }, o($VA, [2, 15]), o($VA, [2, 16]), o($VA, [2, 17]), { 39: [2, 35] }, { 15: 122, 16: 84, 17: [1, 121], 39: [2, 9], 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn }, o($VQ, [2, 43], { 11: 123, 12: [1, 124] }), o($Vp, [2, 7]), { 9: [1, 125] }, o($VR, [2, 52]), { 16: 39, 18: 126, 19: 40, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }, { 13: [1, 128], 16: 39, 18: 127, 19: 40, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }, o($VC, [2, 67], { 64: 129, 66: $Vr, 67: $Vs, 68: $Vt, 69: $Vu, 70: $Vv }), o($VC, [2, 66]), { 41: [1, 130] }, { 24: 95, 40: 131, 43: 23, 46: $V4 }, { 8: [1, 132], 41: [2, 36] }, o($Vy, [2, 40], { 39: [1, 133] }), { 41: [1, 134] }, { 41: [2, 46], 45: 135, 49: $VF }, { 16: 39, 18: 136, 19: 40, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }, o($Vq, [2, 76], { 13: [1, 137] }), o($Vq, [2, 78], { 13: [1, 139], 75: [1, 138] }), o($Vq, [2, 82], { 13: [1, 140], 78: [1, 141] }), { 13: [1, 142] }, o($Vq, [2, 90], { 59: $VS }), o($VT, [2, 92], { 83: 144, 22: $VG, 58: $VH, 80: $VI, 84: $VJ, 85: $VK, 86: $VL, 87: $VM, 88: $VN, 89: $VO }), o($VU, [2, 94]), o($VU, [2, 96]), o($VU, [2, 97]), o($VU, [2, 98]), o($VU, [2, 99]), o($VU, [2, 100]), o($VU, [2, 101]), o($VU, [2, 102]), o($VU, [2, 103]), o($VU, [2, 104]), o($Vq, [2, 91]), o($Vq, [2, 56]), o($Vq, [2, 58], { 59: $VS }), { 58: [1, 145] }, o($VA, [2, 13]), { 15: 146, 16: 84, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn }, { 39: [2, 11] }, o($VQ, [2, 44]), { 13: [1, 147] }, { 1: [2, 4] }, o($VR, [2, 54]), o($VR, [2, 53]), { 16: 39, 18: 148, 19: 40, 58: $Vb, 84: $Vl, 99: $Vm, 101: $Vn, 102: $Vo }, o($VC, [2, 65]), o($Vq, [2, 33]), { 41: [1, 149] }, { 24: 95, 40: 150, 41: [2, 37], 43: 23, 46: $V4 }, { 45: 151, 49: $VF }, o($Vy, [2, 41]), { 41: [2, 47] }, o($Vq, [2, 45]), o($Vq, [2, 77]), o($Vq, [2, 79]), o($Vq, [2, 80], { 75: [1, 152] }), o($Vq, [2, 83]), o($Vq, [2, 84], { 13: [1, 153] }), o($Vq, [2, 86], { 13: [1, 155], 75: [1, 154] }), { 22: $VG, 58: $VH, 80: $VI, 82: 156, 83: 106, 84: $VJ, 85: $VK, 86: $VL, 87: $VM, 88: $VN, 89: $VO }, o($VU, [2, 95]), o($VP, [2, 60]), { 39: [2, 10] }, { 14: [1, 157] }, o($VR, [2, 55]), o($Vq, [2, 34]), { 41: [2, 38] }, { 41: [1, 158] }, o($Vq, [2, 81]), o($Vq, [2, 85]), o($Vq, [2, 87]), o($Vq, [2, 88], { 75: [1, 159] }), o($VT, [2, 93], { 83: 144, 22: $VG, 58: $VH, 80: $VI, 84: $VJ, 85: $VK, 86: $VL, 87: $VM, 88: $VN, 89: $VO }), o($VQ, [2, 8]), o($Vy, [2, 42]), o($Vq, [2, 89])],
    defaultActions: { 2: [2, 1], 3: [2, 2], 4: [2, 3], 83: [2, 35], 122: [2, 11], 125: [2, 4], 135: [2, 47], 146: [2, 10], 150: [2, 38] },
    parseError: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function parseError(str, hash) {
      if (hash.recoverable) {
        this.trace(str);
      } else {
        var error = new Error(str);
        error.hash = hash;
        throw error;
      }
    }, "parseError"),
    parse: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function parse(input) {
      var self = this, stack = [0], tstack = [], vstack = [null], lstack = [], table = this.table, yytext = "", yylineno = 0, yyleng = 0, recovering = 0, TERROR = 2, EOF = 1;
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
      function popStack(n) {
        stack.length = stack.length - 2 * n;
        vstack.length = vstack.length - n;
        lstack.length = lstack.length - n;
      }
      (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(popStack, "popStack");
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
      (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(lex, "lex");
      var symbol, preErrorSymbol, state, action, a, r, yyval = {}, p, len, newState, expected;
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
            if (!preErrorSymbol) {
              yyleng = lexer2.yyleng;
              yytext = lexer2.yytext;
              yylineno = lexer2.yylineno;
              yyloc = lexer2.yylloc;
              if (recovering > 0) {
                recovering--;
              }
            } else {
              symbol = preErrorSymbol;
              preErrorSymbol = null;
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
    }, "parse")
  };
  var lexer = /* @__PURE__ */ function() {
    var lexer2 = {
      EOF: 1,
      parseError: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function parseError(str, hash) {
        if (this.yy.parser) {
          this.yy.parser.parseError(str, hash);
        } else {
          throw new Error(str);
        }
      }, "parseError"),
      // resets the lexer, sets new input
      setInput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(input, yy) {
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
      }, "setInput"),
      // consumes and returns one char from the input
      input: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
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
      }, "input"),
      // unshifts one char (or a string) into the input
      unput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(ch) {
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
      }, "unput"),
      // When called from action, caches matched text and appends it on next action
      more: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        this._more = true;
        return this;
      }, "more"),
      // When called from action, signals the lexer that this rule fails to match the input, so the next matching rule (regex) should be tested instead.
      reject: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
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
      }, "reject"),
      // retain first n characters of the match
      less: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(n) {
        this.unput(this.match.slice(n));
      }, "less"),
      // displays already matched input, i.e. for error messages
      pastInput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        var past = this.matched.substr(0, this.matched.length - this.match.length);
        return (past.length > 20 ? "..." : "") + past.substr(-20).replace(/\n/g, "");
      }, "pastInput"),
      // displays upcoming input, i.e. for error messages
      upcomingInput: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        var next = this.match;
        if (next.length < 20) {
          next += this._input.substr(0, 20 - next.length);
        }
        return (next.substr(0, 20) + (next.length > 20 ? "..." : "")).replace(/\n/g, "");
      }, "upcomingInput"),
      // displays the character position where the lexing error occurred, i.e. for error messages
      showPosition: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
        var pre = this.pastInput();
        var c = new Array(pre.length + 1).join("-");
        return pre + this.upcomingInput() + "\n" + c + "^";
      }, "showPosition"),
      // test the lexed token: return FALSE when not a match, otherwise return token
      test_match: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(match, indexed_rule) {
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
      }, "test_match"),
      // return next match in input
      next: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function() {
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
      }, "next"),
      // return next match that has a token
      lex: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function lex() {
        var r = this.next();
        if (r) {
          return r;
        } else {
          return this.lex();
        }
      }, "lex"),
      // activates a new lexer condition state (pushes the new lexer condition state onto the condition stack)
      begin: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function begin(condition) {
        this.conditionStack.push(condition);
      }, "begin"),
      // pop the previously active lexer condition state off the condition stack
      popState: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function popState() {
        var n = this.conditionStack.length - 1;
        if (n > 0) {
          return this.conditionStack.pop();
        } else {
          return this.conditionStack[0];
        }
      }, "popState"),
      // produce the lexer rule set which is active for the currently active lexer condition state
      _currentRules: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function _currentRules() {
        if (this.conditionStack.length && this.conditionStack[this.conditionStack.length - 1]) {
          return this.conditions[this.conditionStack[this.conditionStack.length - 1]].rules;
        } else {
          return this.conditions["INITIAL"].rules;
        }
      }, "_currentRules"),
      // return the currently active lexer condition state; when an index argument is provided it produces the N-th previous condition state, if available
      topState: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function topState(n) {
        n = this.conditionStack.length - 1 - Math.abs(n || 0);
        if (n >= 0) {
          return this.conditionStack[n];
        } else {
          return "INITIAL";
        }
      }, "topState"),
      // alias for begin(condition)
      pushState: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function pushState(condition) {
        this.begin(condition);
      }, "pushState"),
      // return the number of states currently on the stack
      stateStackSize: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function stateStackSize() {
        return this.conditionStack.length;
      }, "stateStackSize"),
      options: {},
      performAction: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function anonymous(yy, yy_, $avoiding_name_collisions, YY_START) {
        var YYSTATE = YY_START;
        switch ($avoiding_name_collisions) {
          case 0:
            return 60;
            break;
          case 1:
            return 61;
            break;
          case 2:
            return 62;
            break;
          case 3:
            return 63;
            break;
          case 4:
            break;
          case 5:
            break;
          case 6:
            this.begin("acc_title");
            return 33;
            break;
          case 7:
            this.popState();
            return "acc_title_value";
            break;
          case 8:
            this.begin("acc_descr");
            return 35;
            break;
          case 9:
            this.popState();
            return "acc_descr_value";
            break;
          case 10:
            this.begin("acc_descr_multiline");
            break;
          case 11:
            this.popState();
            break;
          case 12:
            return "acc_descr_multiline_value";
            break;
          case 13:
            return 8;
            break;
          case 14:
            break;
          case 15:
            return 7;
            break;
          case 16:
            return 7;
            break;
          case 17:
            return "EDGE_STATE";
            break;
          case 18:
            this.begin("callback_name");
            break;
          case 19:
            this.popState();
            break;
          case 20:
            this.popState();
            this.begin("callback_args");
            break;
          case 21:
            return 77;
            break;
          case 22:
            this.popState();
            break;
          case 23:
            return 78;
            break;
          case 24:
            this.popState();
            break;
          case 25:
            return "STR";
            break;
          case 26:
            this.begin("string");
            break;
          case 27:
            return 80;
            break;
          case 28:
            return 55;
            break;
          case 29:
            this.begin("namespace");
            return 42;
            break;
          case 30:
            this.popState();
            return 8;
            break;
          case 31:
            break;
          case 32:
            this.begin("namespace-body");
            return 39;
            break;
          case 33:
            this.popState();
            return 41;
            break;
          case 34:
            return "EOF_IN_STRUCT";
            break;
          case 35:
            return 8;
            break;
          case 36:
            break;
          case 37:
            return "EDGE_STATE";
            break;
          case 38:
            this.begin("class");
            return 46;
            break;
          case 39:
            this.popState();
            return 8;
            break;
          case 40:
            break;
          case 41:
            this.popState();
            this.popState();
            return 41;
            break;
          case 42:
            this.begin("class-body");
            return 39;
            break;
          case 43:
            this.popState();
            return 41;
            break;
          case 44:
            return "EOF_IN_STRUCT";
            break;
          case 45:
            return "EDGE_STATE";
            break;
          case 46:
            return "OPEN_IN_STRUCT";
            break;
          case 47:
            break;
          case 48:
            return "MEMBER";
            break;
          case 49:
            return 81;
            break;
          case 50:
            return 73;
            break;
          case 51:
            return 74;
            break;
          case 52:
            return 76;
            break;
          case 53:
            return 52;
            break;
          case 54:
            return 54;
            break;
          case 55:
            return 47;
            break;
          case 56:
            return 48;
            break;
          case 57:
            return 79;
            break;
          case 58:
            this.popState();
            break;
          case 59:
            return "GENERICTYPE";
            break;
          case 60:
            this.begin("generic");
            break;
          case 61:
            this.popState();
            break;
          case 62:
            return "BQUOTE_STR";
            break;
          case 63:
            this.begin("bqstring");
            break;
          case 64:
            return 75;
            break;
          case 65:
            return 75;
            break;
          case 66:
            return 75;
            break;
          case 67:
            return 75;
            break;
          case 68:
            return 67;
            break;
          case 69:
            return 67;
            break;
          case 70:
            return 69;
            break;
          case 71:
            return 69;
            break;
          case 72:
            return 68;
            break;
          case 73:
            return 66;
            break;
          case 74:
            return 70;
            break;
          case 75:
            return 71;
            break;
          case 76:
            return 72;
            break;
          case 77:
            return 22;
            break;
          case 78:
            return 44;
            break;
          case 79:
            return 99;
            break;
          case 80:
            return 17;
            break;
          case 81:
            return "PLUS";
            break;
          case 82:
            return 85;
            break;
          case 83:
            return 59;
            break;
          case 84:
            return 88;
            break;
          case 85:
            return 88;
            break;
          case 86:
            return 89;
            break;
          case 87:
            return "EQUALS";
            break;
          case 88:
            return "EQUALS";
            break;
          case 89:
            return 58;
            break;
          case 90:
            return 12;
            break;
          case 91:
            return 14;
            break;
          case 92:
            return "PUNCTUATION";
            break;
          case 93:
            return 84;
            break;
          case 94:
            return 101;
            break;
          case 95:
            return 87;
            break;
          case 96:
            return 87;
            break;
          case 97:
            return 9;
            break;
        }
      }, "anonymous"),
      rules: [/^(?:.*direction\s+TB[^\n]*)/, /^(?:.*direction\s+BT[^\n]*)/, /^(?:.*direction\s+RL[^\n]*)/, /^(?:.*direction\s+LR[^\n]*)/, /^(?:%%(?!\{)*[^\n]*(\r?\n?)+)/, /^(?:%%[^\n]*(\r?\n)*)/, /^(?:accTitle\s*:\s*)/, /^(?:(?!\n||)*[^\n]*)/, /^(?:accDescr\s*:\s*)/, /^(?:(?!\n||)*[^\n]*)/, /^(?:accDescr\s*\{\s*)/, /^(?:[\}])/, /^(?:[^\}]*)/, /^(?:\s*(\r?\n)+)/, /^(?:\s+)/, /^(?:classDiagram-v2\b)/, /^(?:classDiagram\b)/, /^(?:\[\*\])/, /^(?:call[\s]+)/, /^(?:\([\s]*\))/, /^(?:\()/, /^(?:[^(]*)/, /^(?:\))/, /^(?:[^)]*)/, /^(?:["])/, /^(?:[^"]*)/, /^(?:["])/, /^(?:style\b)/, /^(?:classDef\b)/, /^(?:namespace\b)/, /^(?:\s*(\r?\n)+)/, /^(?:\s+)/, /^(?:[{])/, /^(?:[}])/, /^(?:$)/, /^(?:\s*(\r?\n)+)/, /^(?:\s+)/, /^(?:\[\*\])/, /^(?:class\b)/, /^(?:\s*(\r?\n)+)/, /^(?:\s+)/, /^(?:[}])/, /^(?:[{])/, /^(?:[}])/, /^(?:$)/, /^(?:\[\*\])/, /^(?:[{])/, /^(?:[\n])/, /^(?:[^{}\n]*)/, /^(?:cssClass\b)/, /^(?:callback\b)/, /^(?:link\b)/, /^(?:click\b)/, /^(?:note for\b)/, /^(?:note\b)/, /^(?:<<)/, /^(?:>>)/, /^(?:href\b)/, /^(?:[~])/, /^(?:[^~]*)/, /^(?:~)/, /^(?:[`])/, /^(?:[^`]+)/, /^(?:[`])/, /^(?:_self\b)/, /^(?:_blank\b)/, /^(?:_parent\b)/, /^(?:_top\b)/, /^(?:\s*<\|)/, /^(?:\s*\|>)/, /^(?:\s*>)/, /^(?:\s*<)/, /^(?:\s*\*)/, /^(?:\s*o\b)/, /^(?:\s*\(\))/, /^(?:--)/, /^(?:\.\.)/, /^(?::{1}[^:\n;]+)/, /^(?::{3})/, /^(?:-)/, /^(?:\.)/, /^(?:\+)/, /^(?::)/, /^(?:,)/, /^(?:#)/, /^(?:#)/, /^(?:%)/, /^(?:=)/, /^(?:=)/, /^(?:\w+)/, /^(?:\[)/, /^(?:\])/, /^(?:[!"#$%&'*+,-.`?\\/])/, /^(?:[0-9]+)/, /^(?:[\u00AA\u00B5\u00BA\u00C0-\u00D6\u00D8-\u00F6]|[\u00F8-\u02C1\u02C6-\u02D1\u02E0-\u02E4\u02EC\u02EE\u0370-\u0374\u0376\u0377]|[\u037A-\u037D\u0386\u0388-\u038A\u038C\u038E-\u03A1\u03A3-\u03F5]|[\u03F7-\u0481\u048A-\u0527\u0531-\u0556\u0559\u0561-\u0587\u05D0-\u05EA]|[\u05F0-\u05F2\u0620-\u064A\u066E\u066F\u0671-\u06D3\u06D5\u06E5\u06E6\u06EE]|[\u06EF\u06FA-\u06FC\u06FF\u0710\u0712-\u072F\u074D-\u07A5\u07B1\u07CA-\u07EA]|[\u07F4\u07F5\u07FA\u0800-\u0815\u081A\u0824\u0828\u0840-\u0858\u08A0]|[\u08A2-\u08AC\u0904-\u0939\u093D\u0950\u0958-\u0961\u0971-\u0977]|[\u0979-\u097F\u0985-\u098C\u098F\u0990\u0993-\u09A8\u09AA-\u09B0\u09B2]|[\u09B6-\u09B9\u09BD\u09CE\u09DC\u09DD\u09DF-\u09E1\u09F0\u09F1\u0A05-\u0A0A]|[\u0A0F\u0A10\u0A13-\u0A28\u0A2A-\u0A30\u0A32\u0A33\u0A35\u0A36\u0A38\u0A39]|[\u0A59-\u0A5C\u0A5E\u0A72-\u0A74\u0A85-\u0A8D\u0A8F-\u0A91\u0A93-\u0AA8]|[\u0AAA-\u0AB0\u0AB2\u0AB3\u0AB5-\u0AB9\u0ABD\u0AD0\u0AE0\u0AE1\u0B05-\u0B0C]|[\u0B0F\u0B10\u0B13-\u0B28\u0B2A-\u0B30\u0B32\u0B33\u0B35-\u0B39\u0B3D\u0B5C]|[\u0B5D\u0B5F-\u0B61\u0B71\u0B83\u0B85-\u0B8A\u0B8E-\u0B90\u0B92-\u0B95\u0B99]|[\u0B9A\u0B9C\u0B9E\u0B9F\u0BA3\u0BA4\u0BA8-\u0BAA\u0BAE-\u0BB9\u0BD0]|[\u0C05-\u0C0C\u0C0E-\u0C10\u0C12-\u0C28\u0C2A-\u0C33\u0C35-\u0C39\u0C3D]|[\u0C58\u0C59\u0C60\u0C61\u0C85-\u0C8C\u0C8E-\u0C90\u0C92-\u0CA8\u0CAA-\u0CB3]|[\u0CB5-\u0CB9\u0CBD\u0CDE\u0CE0\u0CE1\u0CF1\u0CF2\u0D05-\u0D0C\u0D0E-\u0D10]|[\u0D12-\u0D3A\u0D3D\u0D4E\u0D60\u0D61\u0D7A-\u0D7F\u0D85-\u0D96\u0D9A-\u0DB1]|[\u0DB3-\u0DBB\u0DBD\u0DC0-\u0DC6\u0E01-\u0E30\u0E32\u0E33\u0E40-\u0E46\u0E81]|[\u0E82\u0E84\u0E87\u0E88\u0E8A\u0E8D\u0E94-\u0E97\u0E99-\u0E9F\u0EA1-\u0EA3]|[\u0EA5\u0EA7\u0EAA\u0EAB\u0EAD-\u0EB0\u0EB2\u0EB3\u0EBD\u0EC0-\u0EC4\u0EC6]|[\u0EDC-\u0EDF\u0F00\u0F40-\u0F47\u0F49-\u0F6C\u0F88-\u0F8C\u1000-\u102A]|[\u103F\u1050-\u1055\u105A-\u105D\u1061\u1065\u1066\u106E-\u1070\u1075-\u1081]|[\u108E\u10A0-\u10C5\u10C7\u10CD\u10D0-\u10FA\u10FC-\u1248\u124A-\u124D]|[\u1250-\u1256\u1258\u125A-\u125D\u1260-\u1288\u128A-\u128D\u1290-\u12B0]|[\u12B2-\u12B5\u12B8-\u12BE\u12C0\u12C2-\u12C5\u12C8-\u12D6\u12D8-\u1310]|[\u1312-\u1315\u1318-\u135A\u1380-\u138F\u13A0-\u13F4\u1401-\u166C]|[\u166F-\u167F\u1681-\u169A\u16A0-\u16EA\u1700-\u170C\u170E-\u1711]|[\u1720-\u1731\u1740-\u1751\u1760-\u176C\u176E-\u1770\u1780-\u17B3\u17D7]|[\u17DC\u1820-\u1877\u1880-\u18A8\u18AA\u18B0-\u18F5\u1900-\u191C]|[\u1950-\u196D\u1970-\u1974\u1980-\u19AB\u19C1-\u19C7\u1A00-\u1A16]|[\u1A20-\u1A54\u1AA7\u1B05-\u1B33\u1B45-\u1B4B\u1B83-\u1BA0\u1BAE\u1BAF]|[\u1BBA-\u1BE5\u1C00-\u1C23\u1C4D-\u1C4F\u1C5A-\u1C7D\u1CE9-\u1CEC]|[\u1CEE-\u1CF1\u1CF5\u1CF6\u1D00-\u1DBF\u1E00-\u1F15\u1F18-\u1F1D]|[\u1F20-\u1F45\u1F48-\u1F4D\u1F50-\u1F57\u1F59\u1F5B\u1F5D\u1F5F-\u1F7D]|[\u1F80-\u1FB4\u1FB6-\u1FBC\u1FBE\u1FC2-\u1FC4\u1FC6-\u1FCC\u1FD0-\u1FD3]|[\u1FD6-\u1FDB\u1FE0-\u1FEC\u1FF2-\u1FF4\u1FF6-\u1FFC\u2071\u207F]|[\u2090-\u209C\u2102\u2107\u210A-\u2113\u2115\u2119-\u211D\u2124\u2126\u2128]|[\u212A-\u212D\u212F-\u2139\u213C-\u213F\u2145-\u2149\u214E\u2183\u2184]|[\u2C00-\u2C2E\u2C30-\u2C5E\u2C60-\u2CE4\u2CEB-\u2CEE\u2CF2\u2CF3]|[\u2D00-\u2D25\u2D27\u2D2D\u2D30-\u2D67\u2D6F\u2D80-\u2D96\u2DA0-\u2DA6]|[\u2DA8-\u2DAE\u2DB0-\u2DB6\u2DB8-\u2DBE\u2DC0-\u2DC6\u2DC8-\u2DCE]|[\u2DD0-\u2DD6\u2DD8-\u2DDE\u2E2F\u3005\u3006\u3031-\u3035\u303B\u303C]|[\u3041-\u3096\u309D-\u309F\u30A1-\u30FA\u30FC-\u30FF\u3105-\u312D]|[\u3131-\u318E\u31A0-\u31BA\u31F0-\u31FF\u3400-\u4DB5\u4E00-\u9FCC]|[\uA000-\uA48C\uA4D0-\uA4FD\uA500-\uA60C\uA610-\uA61F\uA62A\uA62B]|[\uA640-\uA66E\uA67F-\uA697\uA6A0-\uA6E5\uA717-\uA71F\uA722-\uA788]|[\uA78B-\uA78E\uA790-\uA793\uA7A0-\uA7AA\uA7F8-\uA801\uA803-\uA805]|[\uA807-\uA80A\uA80C-\uA822\uA840-\uA873\uA882-\uA8B3\uA8F2-\uA8F7\uA8FB]|[\uA90A-\uA925\uA930-\uA946\uA960-\uA97C\uA984-\uA9B2\uA9CF\uAA00-\uAA28]|[\uAA40-\uAA42\uAA44-\uAA4B\uAA60-\uAA76\uAA7A\uAA80-\uAAAF\uAAB1\uAAB5]|[\uAAB6\uAAB9-\uAABD\uAAC0\uAAC2\uAADB-\uAADD\uAAE0-\uAAEA\uAAF2-\uAAF4]|[\uAB01-\uAB06\uAB09-\uAB0E\uAB11-\uAB16\uAB20-\uAB26\uAB28-\uAB2E]|[\uABC0-\uABE2\uAC00-\uD7A3\uD7B0-\uD7C6\uD7CB-\uD7FB\uF900-\uFA6D]|[\uFA70-\uFAD9\uFB00-\uFB06\uFB13-\uFB17\uFB1D\uFB1F-\uFB28\uFB2A-\uFB36]|[\uFB38-\uFB3C\uFB3E\uFB40\uFB41\uFB43\uFB44\uFB46-\uFBB1\uFBD3-\uFD3D]|[\uFD50-\uFD8F\uFD92-\uFDC7\uFDF0-\uFDFB\uFE70-\uFE74\uFE76-\uFEFC]|[\uFF21-\uFF3A\uFF41-\uFF5A\uFF66-\uFFBE\uFFC2-\uFFC7\uFFCA-\uFFCF]|[\uFFD2-\uFFD7\uFFDA-\uFFDC])/, /^(?:\s)/, /^(?:\s)/, /^(?:$)/],
      conditions: { "namespace-body": { "rules": [26, 33, 34, 35, 36, 37, 38, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "namespace": { "rules": [26, 29, 30, 31, 32, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "class-body": { "rules": [26, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "class": { "rules": [26, 39, 40, 41, 42, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "acc_descr_multiline": { "rules": [11, 12, 26, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "acc_descr": { "rules": [9, 26, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "acc_title": { "rules": [7, 26, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "callback_args": { "rules": [22, 23, 26, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "callback_name": { "rules": [19, 20, 21, 26, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "href": { "rules": [26, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "struct": { "rules": [26, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "generic": { "rules": [26, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "bqstring": { "rules": [26, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "string": { "rules": [24, 25, 26, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 97], "inclusive": false }, "INITIAL": { "rules": [0, 1, 2, 3, 4, 5, 6, 8, 10, 13, 14, 15, 16, 17, 18, 26, 27, 28, 29, 38, 49, 50, 51, 52, 53, 54, 55, 56, 57, 60, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97], "inclusive": true } }
    };
    return lexer2;
  }();
  parser2.lexer = lexer;
  function Parser() {
    this.yy = {};
  }
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(Parser, "Parser");
  Parser.prototype = parser2;
  parser2.Parser = Parser;
  return new Parser();
}();
parser.parser = parser;
var classDiagram_default = parser;

// src/diagrams/class/classDb.ts


// src/diagrams/class/classTypes.ts
var visibilityValues = ["#", "+", "~", "-", ""];
var ClassMember = class {
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(this, "ClassMember");
  }
  constructor(input, memberType) {
    this.memberType = memberType;
    this.visibility = "";
    this.classifier = "";
    this.text = "";
    const sanitizedInput = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .sanitizeText */ .oO)(input, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)());
    this.parseMember(sanitizedInput);
  }
  getDisplayDetails() {
    let displayText = this.visibility + (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .parseGenericTypes */ .UO)(this.id);
    if (this.memberType === "method") {
      displayText += `(${(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .parseGenericTypes */ .UO)(this.parameters.trim())})`;
      if (this.returnType) {
        displayText += " : " + (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .parseGenericTypes */ .UO)(this.returnType);
      }
    }
    displayText = displayText.trim();
    const cssStyle = this.parseClassifier();
    return {
      displayText,
      cssStyle
    };
  }
  parseMember(input) {
    let potentialClassifier = "";
    if (this.memberType === "method") {
      const methodRegEx = /([#+~-])?(.+)\((.*)\)([\s$*])?(.*)([$*])?/;
      const match = methodRegEx.exec(input);
      if (match) {
        const detectedVisibility = match[1] ? match[1].trim() : "";
        if (visibilityValues.includes(detectedVisibility)) {
          this.visibility = detectedVisibility;
        }
        this.id = match[2];
        this.parameters = match[3] ? match[3].trim() : "";
        potentialClassifier = match[4] ? match[4].trim() : "";
        this.returnType = match[5] ? match[5].trim() : "";
        if (potentialClassifier === "") {
          const lastChar = this.returnType.substring(this.returnType.length - 1);
          if (/[$*]/.exec(lastChar)) {
            potentialClassifier = lastChar;
            this.returnType = this.returnType.substring(0, this.returnType.length - 1);
          }
        }
      }
    } else {
      const length = input.length;
      const firstChar = input.substring(0, 1);
      const lastChar = input.substring(length - 1);
      if (visibilityValues.includes(firstChar)) {
        this.visibility = firstChar;
      }
      if (/[$*]/.exec(lastChar)) {
        potentialClassifier = lastChar;
      }
      this.id = input.substring(
        this.visibility === "" ? 0 : 1,
        potentialClassifier === "" ? length : length - 1
      );
    }
    this.classifier = potentialClassifier;
    this.id = this.id.startsWith(" ") ? " " + this.id.trim() : this.id.trim();
    const combinedText = `${this.visibility ? "\\" + this.visibility : ""}${(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .parseGenericTypes */ .UO)(this.id)}${this.memberType === "method" ? `(${(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .parseGenericTypes */ .UO)(this.parameters)})${this.returnType ? " : " + (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .parseGenericTypes */ .UO)(this.returnType) : ""}` : ""}`;
    this.text = combinedText.replaceAll("<", "&lt;").replaceAll(">", "&gt;");
    if (this.text.startsWith("\\&lt;")) {
      this.text = this.text.replace("\\&lt;", "~");
    }
  }
  parseClassifier() {
    switch (this.classifier) {
      case "*":
        return "font-style:italic;";
      case "$":
        return "text-decoration:underline;";
      default:
        return "";
    }
  }
};

// src/diagrams/class/classDb.ts
var MERMAID_DOM_ID_PREFIX = "classId-";
var classCounter = 0;
var sanitizeText2 = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((txt) => _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.sanitizeText(txt, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)()), "sanitizeText");
var ClassDB = class {
  constructor() {
    this.relations = [];
    this.classes = /* @__PURE__ */ new Map();
    this.styleClasses = /* @__PURE__ */ new Map();
    this.notes = [];
    this.interfaces = [];
    // private static classCounter = 0;
    this.namespaces = /* @__PURE__ */ new Map();
    this.namespaceCounter = 0;
    this.functions = [];
    this.lineType = {
      LINE: 0,
      DOTTED_LINE: 1
    };
    this.relationType = {
      AGGREGATION: 0,
      EXTENSION: 1,
      COMPOSITION: 2,
      DEPENDENCY: 3,
      LOLLIPOP: 4
    };
    this.setupToolTips = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((element) => {
      let tooltipElem = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .select */ .Ys)(".mermaidTooltip");
      if ((tooltipElem._groups || tooltipElem)[0][0] === null) {
        tooltipElem = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .select */ .Ys)("body").append("div").attr("class", "mermaidTooltip").style("opacity", 0);
      }
      const svg = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .select */ .Ys)(element).select("svg");
      const nodes = svg.selectAll("g.node");
      nodes.on("mouseover", (event) => {
        const el = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .select */ .Ys)(event.currentTarget);
        const title = el.attr("title");
        if (title === null) {
          return;
        }
        const rect = this.getBoundingClientRect();
        tooltipElem.transition().duration(200).style("opacity", ".9");
        tooltipElem.text(el.attr("title")).style("left", window.scrollX + rect.left + (rect.right - rect.left) / 2 + "px").style("top", window.scrollY + rect.top - 14 + document.body.scrollTop + "px");
        tooltipElem.html(tooltipElem.html().replace(/&lt;br\/&gt;/g, "<br/>"));
        el.classed("hover", true);
      }).on("mouseout", (event) => {
        tooltipElem.transition().duration(500).style("opacity", 0);
        const el = (0,d3__WEBPACK_IMPORTED_MODULE_4__/* .select */ .Ys)(event.currentTarget);
        el.classed("hover", false);
      });
    }, "setupToolTips");
    this.direction = "TB";
    this.setAccTitle = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccTitle */ .GN;
    this.getAccTitle = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccTitle */ .eu;
    this.setAccDescription = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccDescription */ .U$;
    this.getAccDescription = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccDescription */ .Mx;
    this.setDiagramTitle = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setDiagramTitle */ .g2;
    this.getDiagramTitle = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getDiagramTitle */ .Kr;
    this.getConfig = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)().class, "getConfig");
    this.functions.push(this.setupToolTips.bind(this));
    this.clear();
    this.addRelation = this.addRelation.bind(this);
    this.addClassesToNamespace = this.addClassesToNamespace.bind(this);
    this.addNamespace = this.addNamespace.bind(this);
    this.setCssClass = this.setCssClass.bind(this);
    this.addMembers = this.addMembers.bind(this);
    this.addClass = this.addClass.bind(this);
    this.setClassLabel = this.setClassLabel.bind(this);
    this.addAnnotation = this.addAnnotation.bind(this);
    this.addMember = this.addMember.bind(this);
    this.cleanupLabel = this.cleanupLabel.bind(this);
    this.addNote = this.addNote.bind(this);
    this.defineClass = this.defineClass.bind(this);
    this.setDirection = this.setDirection.bind(this);
    this.setLink = this.setLink.bind(this);
    this.bindFunctions = this.bindFunctions.bind(this);
    this.clear = this.clear.bind(this);
    this.setTooltip = this.setTooltip.bind(this);
    this.setClickEvent = this.setClickEvent.bind(this);
    this.setCssStyle = this.setCssStyle.bind(this);
  }
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(this, "ClassDB");
  }
  splitClassNameAndType(_id) {
    const id = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.sanitizeText(_id, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)());
    let genericType = "";
    let className = id;
    if (id.indexOf("~") > 0) {
      const split = id.split("~");
      className = sanitizeText2(split[0]);
      genericType = sanitizeText2(split[1]);
    }
    return { className, type: genericType };
  }
  setClassLabel(_id, label) {
    const id = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.sanitizeText(_id, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)());
    if (label) {
      label = sanitizeText2(label);
    }
    const { className } = this.splitClassNameAndType(id);
    this.classes.get(className).label = label;
    this.classes.get(className).text = `${label}${this.classes.get(className).type ? `<${this.classes.get(className).type}>` : ""}`;
  }
  /**
   * Function called by parser when a node definition has been found.
   *
   * @param id - Id of the class to add
   * @public
   */
  addClass(_id) {
    const id = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.sanitizeText(_id, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)());
    const { className, type } = this.splitClassNameAndType(id);
    if (this.classes.has(className)) {
      return;
    }
    const name = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.sanitizeText(className, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)());
    this.classes.set(name, {
      id: name,
      type,
      label: name,
      text: `${name}${type ? `&lt;${type}&gt;` : ""}`,
      shape: "classBox",
      cssClasses: "default",
      methods: [],
      members: [],
      annotations: [],
      styles: [],
      domId: MERMAID_DOM_ID_PREFIX + name + "-" + classCounter
    });
    classCounter++;
  }
  addInterface(label, classId) {
    const classInterface = {
      id: `interface${this.interfaces.length}`,
      label,
      classId
    };
    this.interfaces.push(classInterface);
  }
  /**
   * Function to lookup domId from id in the graph definition.
   *
   * @param id - class ID to lookup
   * @public
   */
  lookUpDomId(_id) {
    const id = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.sanitizeText(_id, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)());
    if (this.classes.has(id)) {
      return this.classes.get(id).domId;
    }
    throw new Error("Class not found: " + id);
  }
  clear() {
    this.relations = [];
    this.classes = /* @__PURE__ */ new Map();
    this.notes = [];
    this.interfaces = [];
    this.functions = [];
    this.functions.push(this.setupToolTips.bind(this));
    this.namespaces = /* @__PURE__ */ new Map();
    this.namespaceCounter = 0;
    this.direction = "TB";
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .clear */ .ZH)();
  }
  getClass(id) {
    return this.classes.get(id);
  }
  getClasses() {
    return this.classes;
  }
  getRelations() {
    return this.relations;
  }
  getNotes() {
    return this.notes;
  }
  addRelation(classRelation) {
    _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug("Adding relation: " + JSON.stringify(classRelation));
    const invalidTypes = [
      this.relationType.LOLLIPOP,
      this.relationType.AGGREGATION,
      this.relationType.COMPOSITION,
      this.relationType.DEPENDENCY,
      this.relationType.EXTENSION
    ];
    if (classRelation.relation.type1 === this.relationType.LOLLIPOP && !invalidTypes.includes(classRelation.relation.type2)) {
      this.addClass(classRelation.id2);
      this.addInterface(classRelation.id1, classRelation.id2);
      classRelation.id1 = `interface${this.interfaces.length - 1}`;
    } else if (classRelation.relation.type2 === this.relationType.LOLLIPOP && !invalidTypes.includes(classRelation.relation.type1)) {
      this.addClass(classRelation.id1);
      this.addInterface(classRelation.id2, classRelation.id1);
      classRelation.id2 = `interface${this.interfaces.length - 1}`;
    } else {
      this.addClass(classRelation.id1);
      this.addClass(classRelation.id2);
    }
    classRelation.id1 = this.splitClassNameAndType(classRelation.id1).className;
    classRelation.id2 = this.splitClassNameAndType(classRelation.id2).className;
    classRelation.relationTitle1 = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.sanitizeText(
      classRelation.relationTitle1.trim(),
      (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)()
    );
    classRelation.relationTitle2 = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.sanitizeText(
      classRelation.relationTitle2.trim(),
      (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)()
    );
    this.relations.push(classRelation);
  }
  /**
   * Adds an annotation to the specified class Annotations mark special properties of the given type
   * (like 'interface' or 'service')
   *
   * @param className - The class name
   * @param annotation - The name of the annotation without any brackets
   * @public
   */
  addAnnotation(className, annotation) {
    const validatedClassName = this.splitClassNameAndType(className).className;
    this.classes.get(validatedClassName).annotations.push(annotation);
  }
  /**
   * Adds a member to the specified class
   *
   * @param className - The class name
   * @param member - The full name of the member. If the member is enclosed in `<<brackets>>` it is
   *   treated as an annotation If the member is ending with a closing bracket ) it is treated as a
   *   method Otherwise the member will be treated as a normal property
   * @public
   */
  addMember(className, member) {
    this.addClass(className);
    const validatedClassName = this.splitClassNameAndType(className).className;
    const theClass = this.classes.get(validatedClassName);
    if (typeof member === "string") {
      const memberString = member.trim();
      if (memberString.startsWith("<<") && memberString.endsWith(">>")) {
        theClass.annotations.push(sanitizeText2(memberString.substring(2, memberString.length - 2)));
      } else if (memberString.indexOf(")") > 0) {
        theClass.methods.push(new ClassMember(memberString, "method"));
      } else if (memberString) {
        theClass.members.push(new ClassMember(memberString, "attribute"));
      }
    }
  }
  addMembers(className, members) {
    if (Array.isArray(members)) {
      members.reverse();
      members.forEach((member) => this.addMember(className, member));
    }
  }
  addNote(text, className) {
    const note = {
      id: `note${this.notes.length}`,
      class: className,
      text
    };
    this.notes.push(note);
  }
  cleanupLabel(label) {
    if (label.startsWith(":")) {
      label = label.substring(1);
    }
    return sanitizeText2(label.trim());
  }
  /**
   * Called by parser when assigning cssClass to a class
   *
   * @param ids - Comma separated list of ids
   * @param className - Class to add
   */
  setCssClass(ids, className) {
    ids.split(",").forEach((_id) => {
      let id = _id;
      if (/\d/.exec(_id[0])) {
        id = MERMAID_DOM_ID_PREFIX + id;
      }
      const classNode = this.classes.get(id);
      if (classNode) {
        classNode.cssClasses += " " + className;
      }
    });
  }
  defineClass(ids, style) {
    for (const id of ids) {
      let styleClass = this.styleClasses.get(id);
      if (styleClass === void 0) {
        styleClass = { id, styles: [], textStyles: [] };
        this.styleClasses.set(id, styleClass);
      }
      if (style) {
        style.forEach((s) => {
          if (/color/.exec(s)) {
            const newStyle = s.replace("fill", "bgFill");
            styleClass.textStyles.push(newStyle);
          }
          styleClass.styles.push(s);
        });
      }
      this.classes.forEach((value) => {
        if (value.cssClasses.includes(id)) {
          value.styles.push(...style.flatMap((s) => s.split(",")));
        }
      });
    }
  }
  /**
   * Called by parser when a tooltip is found, e.g. a clickable element.
   *
   * @param ids - Comma separated list of ids
   * @param tooltip - Tooltip to add
   */
  setTooltip(ids, tooltip) {
    ids.split(",").forEach((id) => {
      if (tooltip !== void 0) {
        this.classes.get(id).tooltip = sanitizeText2(tooltip);
      }
    });
  }
  getTooltip(id, namespace) {
    if (namespace && this.namespaces.has(namespace)) {
      return this.namespaces.get(namespace).classes.get(id).tooltip;
    }
    return this.classes.get(id).tooltip;
  }
  /**
   * Called by parser when a link is found. Adds the URL to the vertex data.
   *
   * @param ids - Comma separated list of ids
   * @param linkStr - URL to create a link for
   * @param target - Target of the link, _blank by default as originally defined in the svgDraw.js file
   */
  setLink(ids, linkStr, target) {
    const config = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)();
    ids.split(",").forEach((_id) => {
      let id = _id;
      if (/\d/.exec(_id[0])) {
        id = MERMAID_DOM_ID_PREFIX + id;
      }
      const theClass = this.classes.get(id);
      if (theClass) {
        theClass.link = _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.formatUrl(linkStr, config);
        if (config.securityLevel === "sandbox") {
          theClass.linkTarget = "_top";
        } else if (typeof target === "string") {
          theClass.linkTarget = sanitizeText2(target);
        } else {
          theClass.linkTarget = "_blank";
        }
      }
    });
    this.setCssClass(ids, "clickable");
  }
  /**
   * Called by parser when a click definition is found. Registers an event handler.
   *
   * @param ids - Comma separated list of ids
   * @param functionName - Function to be called on click
   * @param functionArgs - Function args the function should be called with
   */
  setClickEvent(ids, functionName, functionArgs) {
    ids.split(",").forEach((id) => {
      this.setClickFunc(id, functionName, functionArgs);
      this.classes.get(id).haveCallback = true;
    });
    this.setCssClass(ids, "clickable");
  }
  setClickFunc(_domId, functionName, functionArgs) {
    const domId = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY.sanitizeText(_domId, (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)());
    const config = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)();
    if (config.securityLevel !== "loose") {
      return;
    }
    if (functionName === void 0) {
      return;
    }
    const id = domId;
    if (this.classes.has(id)) {
      const elemId = this.lookUpDomId(id);
      let argList = [];
      if (typeof functionArgs === "string") {
        argList = functionArgs.split(/,(?=(?:(?:[^"]*"){2})*[^"]*$)/);
        for (let i = 0; i < argList.length; i++) {
          let item = argList[i].trim();
          if (item.startsWith('"') && item.endsWith('"')) {
            item = item.substr(1, item.length - 2);
          }
          argList[i] = item;
        }
      }
      if (argList.length === 0) {
        argList.push(elemId);
      }
      this.functions.push(() => {
        const elem = document.querySelector(`[id="${elemId}"]`);
        if (elem !== null) {
          elem.addEventListener(
            "click",
            () => {
              _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.runFunc(functionName, ...argList);
            },
            false
          );
        }
      });
    }
  }
  bindFunctions(element) {
    this.functions.forEach((fun) => {
      fun(element);
    });
  }
  getDirection() {
    return this.direction;
  }
  setDirection(dir) {
    this.direction = dir;
  }
  /**
   * Function called by parser when a namespace definition has been found.
   *
   * @param id - Id of the namespace to add
   * @public
   */
  addNamespace(id) {
    if (this.namespaces.has(id)) {
      return;
    }
    this.namespaces.set(id, {
      id,
      classes: /* @__PURE__ */ new Map(),
      children: {},
      domId: MERMAID_DOM_ID_PREFIX + id + "-" + this.namespaceCounter
    });
    this.namespaceCounter++;
  }
  getNamespace(name) {
    return this.namespaces.get(name);
  }
  getNamespaces() {
    return this.namespaces;
  }
  /**
   * Function called by parser when a namespace definition has been found.
   *
   * @param id - Id of the namespace to add
   * @param classNames - Ids of the class to add
   * @public
   */
  addClassesToNamespace(id, classNames) {
    if (!this.namespaces.has(id)) {
      return;
    }
    for (const name of classNames) {
      const { className } = this.splitClassNameAndType(name);
      this.classes.get(className).parent = id;
      this.namespaces.get(id).classes.set(className, this.classes.get(className));
    }
  }
  setCssStyle(id, styles) {
    const thisClass = this.classes.get(id);
    if (!styles || !thisClass) {
      return;
    }
    for (const s of styles) {
      if (s.includes(",")) {
        thisClass.styles.push(...s.split(","));
      } else {
        thisClass.styles.push(s);
      }
    }
  }
  /**
   * Gets the arrow marker for a type index
   *
   * @param type - The type to look for
   * @returns The arrow marker
   */
  getArrowMarker(type) {
    let marker;
    switch (type) {
      case 0:
        marker = "aggregation";
        break;
      case 1:
        marker = "extension";
        break;
      case 2:
        marker = "composition";
        break;
      case 3:
        marker = "dependency";
        break;
      case 4:
        marker = "lollipop";
        break;
      default:
        marker = "none";
    }
    return marker;
  }
  getData() {
    const nodes = [];
    const edges = [];
    const config = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)();
    for (const namespaceKey of this.namespaces.keys()) {
      const namespace = this.namespaces.get(namespaceKey);
      if (namespace) {
        const node = {
          id: namespace.id,
          label: namespace.id,
          isGroup: true,
          padding: config.class.padding ?? 16,
          // parent node must be one of [rect, roundedWithTitle, noteGroup, divider]
          shape: "rect",
          cssStyles: ["fill: none", "stroke: black"],
          look: config.look
        };
        nodes.push(node);
      }
    }
    for (const classKey of this.classes.keys()) {
      const classNode = this.classes.get(classKey);
      if (classNode) {
        const node = classNode;
        node.parentId = classNode.parent;
        node.look = config.look;
        nodes.push(node);
      }
    }
    let cnt = 0;
    for (const note of this.notes) {
      cnt++;
      const noteNode = {
        id: note.id,
        label: note.text,
        isGroup: false,
        shape: "note",
        padding: config.class.padding ?? 6,
        cssStyles: [
          "text-align: left",
          "white-space: nowrap",
          `fill: ${config.themeVariables.noteBkgColor}`,
          `stroke: ${config.themeVariables.noteBorderColor}`
        ],
        look: config.look
      };
      nodes.push(noteNode);
      const noteClassId = this.classes.get(note.class)?.id ?? "";
      if (noteClassId) {
        const edge = {
          id: `edgeNote${cnt}`,
          start: note.id,
          end: noteClassId,
          type: "normal",
          thickness: "normal",
          classes: "relation",
          arrowTypeStart: "none",
          arrowTypeEnd: "none",
          arrowheadStyle: "",
          labelStyle: [""],
          style: ["fill: none"],
          pattern: "dotted",
          look: config.look
        };
        edges.push(edge);
      }
    }
    for (const _interface of this.interfaces) {
      const interfaceNode = {
        id: _interface.id,
        label: _interface.label,
        isGroup: false,
        shape: "rect",
        cssStyles: ["opacity: 0;"],
        look: config.look
      };
      nodes.push(interfaceNode);
    }
    cnt = 0;
    for (const classRelation of this.relations) {
      cnt++;
      const edge = {
        id: (0,_chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getEdgeId */ .Ln)(classRelation.id1, classRelation.id2, {
          prefix: "id",
          counter: cnt
        }),
        start: classRelation.id1,
        end: classRelation.id2,
        type: "normal",
        label: classRelation.title,
        labelpos: "c",
        thickness: "normal",
        classes: "relation",
        arrowTypeStart: this.getArrowMarker(classRelation.relation.type1),
        arrowTypeEnd: this.getArrowMarker(classRelation.relation.type2),
        startLabelRight: classRelation.relationTitle1 === "none" ? "" : classRelation.relationTitle1,
        endLabelLeft: classRelation.relationTitle2 === "none" ? "" : classRelation.relationTitle2,
        arrowheadStyle: "",
        labelStyle: ["display: inline-block"],
        style: classRelation.style || "",
        pattern: classRelation.relation.lineType == 1 ? "dashed" : "solid",
        look: config.look
      };
      edges.push(edge);
    }
    return { nodes, edges, other: {}, config, direction: this.getDirection() };
  }
};

// src/diagrams/class/styles.js
var getStyles = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((options) => `g.classGroup text {
  fill: ${options.nodeBorder || options.classText};
  stroke: none;
  font-family: ${options.fontFamily};
  font-size: 10px;

  .title {
    font-weight: bolder;
  }

}

.nodeLabel, .edgeLabel {
  color: ${options.classText};
}
.edgeLabel .label rect {
  fill: ${options.mainBkg};
}
.label text {
  fill: ${options.classText};
}

.labelBkg {
  background: ${options.mainBkg};
}
.edgeLabel .label span {
  background: ${options.mainBkg};
}

.classTitle {
  font-weight: bolder;
}
.node rect,
  .node circle,
  .node ellipse,
  .node polygon,
  .node path {
    fill: ${options.mainBkg};
    stroke: ${options.nodeBorder};
    stroke-width: 1px;
  }


.divider {
  stroke: ${options.nodeBorder};
  stroke-width: 1;
}

g.clickable {
  cursor: pointer;
}

g.classGroup rect {
  fill: ${options.mainBkg};
  stroke: ${options.nodeBorder};
}

g.classGroup line {
  stroke: ${options.nodeBorder};
  stroke-width: 1;
}

.classLabel .box {
  stroke: none;
  stroke-width: 0;
  fill: ${options.mainBkg};
  opacity: 0.5;
}

.classLabel .label {
  fill: ${options.nodeBorder};
  font-size: 10px;
}

.relation {
  stroke: ${options.lineColor};
  stroke-width: 1;
  fill: none;
}

.dashed-line{
  stroke-dasharray: 3;
}

.dotted-line{
  stroke-dasharray: 1 2;
}

#compositionStart, .composition {
  fill: ${options.lineColor} !important;
  stroke: ${options.lineColor} !important;
  stroke-width: 1;
}

#compositionEnd, .composition {
  fill: ${options.lineColor} !important;
  stroke: ${options.lineColor} !important;
  stroke-width: 1;
}

#dependencyStart, .dependency {
  fill: ${options.lineColor} !important;
  stroke: ${options.lineColor} !important;
  stroke-width: 1;
}

#dependencyStart, .dependency {
  fill: ${options.lineColor} !important;
  stroke: ${options.lineColor} !important;
  stroke-width: 1;
}

#extensionStart, .extension {
  fill: transparent !important;
  stroke: ${options.lineColor} !important;
  stroke-width: 1;
}

#extensionEnd, .extension {
  fill: transparent !important;
  stroke: ${options.lineColor} !important;
  stroke-width: 1;
}

#aggregationStart, .aggregation {
  fill: transparent !important;
  stroke: ${options.lineColor} !important;
  stroke-width: 1;
}

#aggregationEnd, .aggregation {
  fill: transparent !important;
  stroke: ${options.lineColor} !important;
  stroke-width: 1;
}

#lollipopStart, .lollipop {
  fill: ${options.mainBkg} !important;
  stroke: ${options.lineColor} !important;
  stroke-width: 1;
}

#lollipopEnd, .lollipop {
  fill: ${options.mainBkg} !important;
  stroke: ${options.lineColor} !important;
  stroke-width: 1;
}

.edgeTerminals {
  font-size: 11px;
  line-height: initial;
}

.classTitleText {
  text-anchor: middle;
  font-size: 18px;
  fill: ${options.textColor};
}
`, "getStyles");
var styles_default = getStyles;

// src/diagrams/class/classRenderer-v3-unified.ts
var getDir = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((parsedItem, defaultDir = "TB") => {
  if (!parsedItem.doc) {
    return defaultDir;
  }
  let dir = defaultDir;
  for (const parsedItemDoc of parsedItem.doc) {
    if (parsedItemDoc.stmt === "dir") {
      dir = parsedItemDoc.value;
    }
  }
  return dir;
}, "getDir");
var getClasses = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(function(text, diagramObj) {
  return diagramObj.db.getClasses();
}, "getClasses");
var draw = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(async function(text, id, _version, diag) {
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.info("REF0:");
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.info("Drawing class diagram (v3)", id);
  const { securityLevel, state: conf, layout } = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig2 */ .nV)();
  const data4Layout = diag.db.getData();
  const svg = (0,_chunk_RZ5BOZE2_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getDiagramElement */ .q)(id, securityLevel);
  data4Layout.type = diag.type;
  data4Layout.layoutAlgorithm = (0,_chunk_TYCBKAJE_mjs__WEBPACK_IMPORTED_MODULE_1__/* .getRegisteredLayoutAlgorithm */ ._b)(layout);
  data4Layout.nodeSpacing = conf?.nodeSpacing || 50;
  data4Layout.rankSpacing = conf?.rankSpacing || 50;
  data4Layout.markers = ["aggregation", "extension", "composition", "dependency", "lollipop"];
  data4Layout.diagramId = id;
  await (0,_chunk_TYCBKAJE_mjs__WEBPACK_IMPORTED_MODULE_1__/* .render */ .sY)(data4Layout, svg);
  const padding = 8;
  _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_2__/* .utils_default */ .w8.insertTitle(
    svg,
    "classDiagramTitleText",
    conf?.titleTopMargin ?? 25,
    diag.db.getDiagramTitle()
  );
  (0,_chunk_RZ5BOZE2_mjs__WEBPACK_IMPORTED_MODULE_0__/* .setupViewPortForSVG */ .j)(svg, padding, "classDiagram", conf?.useMaxWidth ?? true);
}, "draw");
var classRenderer_v3_unified_default = {
  getClasses,
  draw,
  getDir
};




/***/ }),

/***/ 49479:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   j: () => (/* binding */ setupViewPortForSVG),
/* harmony export */   q: () => (/* binding */ getDiagramElement)
/* harmony export */ });
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(86906);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(83619);


// src/rendering-util/insertElementsForSize.js

var getDiagramElement = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((id, securityLevel) => {
  let sandboxElement;
  if (securityLevel === "sandbox") {
    sandboxElement = (0,d3__WEBPACK_IMPORTED_MODULE_1__/* .select */ .Ys)("#i" + id);
  }
  const root = securityLevel === "sandbox" ? (0,d3__WEBPACK_IMPORTED_MODULE_1__/* .select */ .Ys)(sandboxElement.nodes()[0].contentDocument.body) : (0,d3__WEBPACK_IMPORTED_MODULE_1__/* .select */ .Ys)("body");
  const svg = root.select(`[id="${id}"]`);
  return svg;
}, "getDiagramElement");

// src/rendering-util/setupViewPortForSVG.ts
var setupViewPortForSVG = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((svg, padding, cssDiagram, useMaxWidth) => {
  svg.attr("class", cssDiagram);
  const { width, height, x, y } = calculateDimensionsWithPadding(svg, padding);
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .configureSvgSize */ .v2)(svg, height, width, useMaxWidth);
  const viewBox = createViewBox(x, y, width, height, padding);
  svg.attr("viewBox", viewBox);
  _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .log */ .cM.debug(`viewBox configured: ${viewBox} with padding: ${padding}`);
}, "setupViewPortForSVG");
var calculateDimensionsWithPadding = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((svg, padding) => {
  const bounds = svg.node()?.getBBox() || { width: 0, height: 0, x: 0, y: 0 };
  return {
    width: bounds.width + padding * 2,
    height: bounds.height + padding * 2,
    x: bounds.x,
    y: bounds.y
  };
}, "calculateDimensionsWithPadding");
var createViewBox = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((x, y, width, height, padding) => {
  return `${x - padding} ${y - padding} ${width} ${height}`;
}, "createViewBox");




/***/ })

}]);
//# sourceMappingURL=4886.6084c97eb0f7628908ee.js.map?v=6084c97eb0f7628908ee
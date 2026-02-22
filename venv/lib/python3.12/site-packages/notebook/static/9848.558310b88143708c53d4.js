"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[9848],{

/***/ 40448:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   q: () => (/* binding */ getDiagramElement)
/* harmony export */ });
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(74999);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(35321);


// src/rendering-util/insertElementsForSize.js

var getDiagramElement = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((id, securityLevel) => {
  let sandboxElement;
  if (securityLevel === "sandbox") {
    sandboxElement = (0,d3__WEBPACK_IMPORTED_MODULE_1__/* .select */ .Ys)("#i" + id);
  }
  const root = securityLevel === "sandbox" ? (0,d3__WEBPACK_IMPORTED_MODULE_1__/* .select */ .Ys)(sandboxElement.nodes()[0].contentDocument.body) : (0,d3__WEBPACK_IMPORTED_MODULE_1__/* .select */ .Ys)("body");
  const svg = root.select(`[id="${id}"]`);
  return svg;
}, "getDiagramElement");




/***/ }),

/***/ 51755:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   j: () => (/* binding */ setupViewPortForSVG)
/* harmony export */ });
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(74999);



// src/rendering-util/setupViewPortForSVG.ts
var setupViewPortForSVG = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((svg, padding, cssDiagram, useMaxWidth) => {
  svg.attr("class", cssDiagram);
  const { width, height, x, y } = calculateDimensionsWithPadding(svg, padding);
  (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .configureSvgSize */ .v2)(svg, height, width, useMaxWidth);
  const viewBox = createViewBox(x, y, width, height, padding);
  svg.attr("viewBox", viewBox);
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .log */ .cM.debug(`viewBox configured: ${viewBox} with padding: ${padding}`);
}, "setupViewPortForSVG");
var calculateDimensionsWithPadding = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((svg, padding) => {
  const bounds = svg.node()?.getBBox() || { width: 0, height: 0, x: 0, y: 0 };
  return {
    width: bounds.width + padding * 2,
    height: bounds.height + padding * 2,
    x: bounds.x,
    y: bounds.y
  };
}, "calculateDimensionsWithPadding");
var createViewBox = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((x, y, width, height, padding) => {
  return `${x - padding} ${y - padding} ${width} ${height}`;
}, "createViewBox");




/***/ }),

/***/ 89848:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _chunk_55IACEB6_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(40448);
/* harmony import */ var _chunk_QN33PNHL_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(51755);
/* harmony import */ var _chunk_N4CR4FBY_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(54160);
/* harmony import */ var _chunk_QXUST7PY_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(55195);
/* harmony import */ var _chunk_HN2XXSSU_mjs__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(67894);
/* harmony import */ var _chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(26212);
/* harmony import */ var _chunk_CVBHYZKI_mjs__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(22462);
/* harmony import */ var _chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(37756);
/* harmony import */ var _chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(2111);
/* harmony import */ var _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(17175);
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__ = __webpack_require__(74999);













// src/diagrams/requirement/parser/requirementDiagram.jison
var parser = (function() {
  var o = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function(k, v, o2, l) {
    for (o2 = o2 || {}, l = k.length; l--; o2[k[l]] = v) ;
    return o2;
  }, "o"), $V0 = [1, 3], $V1 = [1, 4], $V2 = [1, 5], $V3 = [1, 6], $V4 = [5, 6, 8, 9, 11, 13, 21, 22, 23, 24, 41, 42, 43, 44, 45, 46, 54, 72, 74, 77, 89, 90], $V5 = [1, 22], $V6 = [2, 7], $V7 = [1, 26], $V8 = [1, 27], $V9 = [1, 28], $Va = [1, 29], $Vb = [1, 33], $Vc = [1, 34], $Vd = [1, 35], $Ve = [1, 36], $Vf = [1, 37], $Vg = [1, 38], $Vh = [1, 24], $Vi = [1, 31], $Vj = [1, 32], $Vk = [1, 30], $Vl = [1, 39], $Vm = [1, 40], $Vn = [5, 8, 9, 11, 13, 21, 22, 23, 24, 41, 42, 43, 44, 45, 46, 54, 72, 74, 77, 89, 90], $Vo = [1, 61], $Vp = [89, 90], $Vq = [5, 8, 9, 11, 13, 21, 22, 23, 24, 27, 29, 41, 42, 43, 44, 45, 46, 54, 61, 63, 72, 74, 75, 76, 77, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90], $Vr = [27, 29], $Vs = [1, 70], $Vt = [1, 71], $Vu = [1, 72], $Vv = [1, 73], $Vw = [1, 74], $Vx = [1, 75], $Vy = [1, 76], $Vz = [1, 83], $VA = [1, 80], $VB = [1, 84], $VC = [1, 85], $VD = [1, 86], $VE = [1, 87], $VF = [1, 88], $VG = [1, 89], $VH = [1, 90], $VI = [1, 91], $VJ = [1, 92], $VK = [5, 8, 9, 11, 13, 21, 22, 23, 24, 27, 41, 42, 43, 44, 45, 46, 54, 72, 74, 75, 76, 77, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90], $VL = [63, 64], $VM = [1, 101], $VN = [5, 8, 9, 11, 13, 21, 22, 23, 24, 41, 42, 43, 44, 45, 46, 54, 72, 74, 76, 77, 89, 90], $VO = [5, 8, 9, 11, 13, 21, 22, 23, 24, 41, 42, 43, 44, 45, 46, 54, 72, 74, 75, 76, 77, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90], $VP = [1, 110], $VQ = [1, 106], $VR = [1, 107], $VS = [1, 108], $VT = [1, 109], $VU = [1, 111], $VV = [1, 116], $VW = [1, 117], $VX = [1, 114], $VY = [1, 115];
  var parser2 = {
    trace: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function trace() {
    }, "trace"),
    yy: {},
    symbols_: { "error": 2, "start": 3, "directive": 4, "NEWLINE": 5, "RD": 6, "diagram": 7, "EOF": 8, "acc_title": 9, "acc_title_value": 10, "acc_descr": 11, "acc_descr_value": 12, "acc_descr_multiline_value": 13, "requirementDef": 14, "elementDef": 15, "relationshipDef": 16, "direction": 17, "styleStatement": 18, "classDefStatement": 19, "classStatement": 20, "direction_tb": 21, "direction_bt": 22, "direction_rl": 23, "direction_lr": 24, "requirementType": 25, "requirementName": 26, "STRUCT_START": 27, "requirementBody": 28, "STYLE_SEPARATOR": 29, "idList": 30, "ID": 31, "COLONSEP": 32, "id": 33, "TEXT": 34, "text": 35, "RISK": 36, "riskLevel": 37, "VERIFYMTHD": 38, "verifyType": 39, "STRUCT_STOP": 40, "REQUIREMENT": 41, "FUNCTIONAL_REQUIREMENT": 42, "INTERFACE_REQUIREMENT": 43, "PERFORMANCE_REQUIREMENT": 44, "PHYSICAL_REQUIREMENT": 45, "DESIGN_CONSTRAINT": 46, "LOW_RISK": 47, "MED_RISK": 48, "HIGH_RISK": 49, "VERIFY_ANALYSIS": 50, "VERIFY_DEMONSTRATION": 51, "VERIFY_INSPECTION": 52, "VERIFY_TEST": 53, "ELEMENT": 54, "elementName": 55, "elementBody": 56, "TYPE": 57, "type": 58, "DOCREF": 59, "ref": 60, "END_ARROW_L": 61, "relationship": 62, "LINE": 63, "END_ARROW_R": 64, "CONTAINS": 65, "COPIES": 66, "DERIVES": 67, "SATISFIES": 68, "VERIFIES": 69, "REFINES": 70, "TRACES": 71, "CLASSDEF": 72, "stylesOpt": 73, "CLASS": 74, "ALPHA": 75, "COMMA": 76, "STYLE": 77, "style": 78, "styleComponent": 79, "NUM": 80, "COLON": 81, "UNIT": 82, "SPACE": 83, "BRKT": 84, "PCT": 85, "MINUS": 86, "LABEL": 87, "SEMICOLON": 88, "unqString": 89, "qString": 90, "$accept": 0, "$end": 1 },
    terminals_: { 2: "error", 5: "NEWLINE", 6: "RD", 8: "EOF", 9: "acc_title", 10: "acc_title_value", 11: "acc_descr", 12: "acc_descr_value", 13: "acc_descr_multiline_value", 21: "direction_tb", 22: "direction_bt", 23: "direction_rl", 24: "direction_lr", 27: "STRUCT_START", 29: "STYLE_SEPARATOR", 31: "ID", 32: "COLONSEP", 34: "TEXT", 36: "RISK", 38: "VERIFYMTHD", 40: "STRUCT_STOP", 41: "REQUIREMENT", 42: "FUNCTIONAL_REQUIREMENT", 43: "INTERFACE_REQUIREMENT", 44: "PERFORMANCE_REQUIREMENT", 45: "PHYSICAL_REQUIREMENT", 46: "DESIGN_CONSTRAINT", 47: "LOW_RISK", 48: "MED_RISK", 49: "HIGH_RISK", 50: "VERIFY_ANALYSIS", 51: "VERIFY_DEMONSTRATION", 52: "VERIFY_INSPECTION", 53: "VERIFY_TEST", 54: "ELEMENT", 57: "TYPE", 59: "DOCREF", 61: "END_ARROW_L", 63: "LINE", 64: "END_ARROW_R", 65: "CONTAINS", 66: "COPIES", 67: "DERIVES", 68: "SATISFIES", 69: "VERIFIES", 70: "REFINES", 71: "TRACES", 72: "CLASSDEF", 74: "CLASS", 75: "ALPHA", 76: "COMMA", 77: "STYLE", 80: "NUM", 81: "COLON", 82: "UNIT", 83: "SPACE", 84: "BRKT", 85: "PCT", 86: "MINUS", 87: "LABEL", 88: "SEMICOLON", 89: "unqString", 90: "qString" },
    productions_: [0, [3, 3], [3, 2], [3, 4], [4, 2], [4, 2], [4, 1], [7, 0], [7, 2], [7, 2], [7, 2], [7, 2], [7, 2], [7, 2], [7, 2], [7, 2], [7, 2], [17, 1], [17, 1], [17, 1], [17, 1], [14, 5], [14, 7], [28, 5], [28, 5], [28, 5], [28, 5], [28, 2], [28, 1], [25, 1], [25, 1], [25, 1], [25, 1], [25, 1], [25, 1], [37, 1], [37, 1], [37, 1], [39, 1], [39, 1], [39, 1], [39, 1], [15, 5], [15, 7], [56, 5], [56, 5], [56, 2], [56, 1], [16, 5], [16, 5], [62, 1], [62, 1], [62, 1], [62, 1], [62, 1], [62, 1], [62, 1], [19, 3], [20, 3], [20, 3], [30, 1], [30, 3], [30, 1], [30, 3], [18, 3], [73, 1], [73, 3], [78, 1], [78, 2], [79, 1], [79, 1], [79, 1], [79, 1], [79, 1], [79, 1], [79, 1], [79, 1], [79, 1], [79, 1], [26, 1], [26, 1], [33, 1], [33, 1], [35, 1], [35, 1], [55, 1], [55, 1], [58, 1], [58, 1], [60, 1], [60, 1]],
    performAction: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function anonymous(yytext, yyleng, yylineno, yy, yystate, $$, _$) {
      var $0 = $$.length - 1;
      switch (yystate) {
        case 4:
          this.$ = $$[$0].trim();
          yy.setAccTitle(this.$);
          break;
        case 5:
        case 6:
          this.$ = $$[$0].trim();
          yy.setAccDescription(this.$);
          break;
        case 7:
          this.$ = [];
          break;
        case 17:
          yy.setDirection("TB");
          break;
        case 18:
          yy.setDirection("BT");
          break;
        case 19:
          yy.setDirection("RL");
          break;
        case 20:
          yy.setDirection("LR");
          break;
        case 21:
          yy.addRequirement($$[$0 - 3], $$[$0 - 4]);
          break;
        case 22:
          yy.addRequirement($$[$0 - 5], $$[$0 - 6]);
          yy.setClass([$$[$0 - 5]], $$[$0 - 3]);
          break;
        case 23:
          yy.setNewReqId($$[$0 - 2]);
          break;
        case 24:
          yy.setNewReqText($$[$0 - 2]);
          break;
        case 25:
          yy.setNewReqRisk($$[$0 - 2]);
          break;
        case 26:
          yy.setNewReqVerifyMethod($$[$0 - 2]);
          break;
        case 29:
          this.$ = yy.RequirementType.REQUIREMENT;
          break;
        case 30:
          this.$ = yy.RequirementType.FUNCTIONAL_REQUIREMENT;
          break;
        case 31:
          this.$ = yy.RequirementType.INTERFACE_REQUIREMENT;
          break;
        case 32:
          this.$ = yy.RequirementType.PERFORMANCE_REQUIREMENT;
          break;
        case 33:
          this.$ = yy.RequirementType.PHYSICAL_REQUIREMENT;
          break;
        case 34:
          this.$ = yy.RequirementType.DESIGN_CONSTRAINT;
          break;
        case 35:
          this.$ = yy.RiskLevel.LOW_RISK;
          break;
        case 36:
          this.$ = yy.RiskLevel.MED_RISK;
          break;
        case 37:
          this.$ = yy.RiskLevel.HIGH_RISK;
          break;
        case 38:
          this.$ = yy.VerifyType.VERIFY_ANALYSIS;
          break;
        case 39:
          this.$ = yy.VerifyType.VERIFY_DEMONSTRATION;
          break;
        case 40:
          this.$ = yy.VerifyType.VERIFY_INSPECTION;
          break;
        case 41:
          this.$ = yy.VerifyType.VERIFY_TEST;
          break;
        case 42:
          yy.addElement($$[$0 - 3]);
          break;
        case 43:
          yy.addElement($$[$0 - 5]);
          yy.setClass([$$[$0 - 5]], $$[$0 - 3]);
          break;
        case 44:
          yy.setNewElementType($$[$0 - 2]);
          break;
        case 45:
          yy.setNewElementDocRef($$[$0 - 2]);
          break;
        case 48:
          yy.addRelationship($$[$0 - 2], $$[$0], $$[$0 - 4]);
          break;
        case 49:
          yy.addRelationship($$[$0 - 2], $$[$0 - 4], $$[$0]);
          break;
        case 50:
          this.$ = yy.Relationships.CONTAINS;
          break;
        case 51:
          this.$ = yy.Relationships.COPIES;
          break;
        case 52:
          this.$ = yy.Relationships.DERIVES;
          break;
        case 53:
          this.$ = yy.Relationships.SATISFIES;
          break;
        case 54:
          this.$ = yy.Relationships.VERIFIES;
          break;
        case 55:
          this.$ = yy.Relationships.REFINES;
          break;
        case 56:
          this.$ = yy.Relationships.TRACES;
          break;
        case 57:
          this.$ = $$[$0 - 2];
          yy.defineClass($$[$0 - 1], $$[$0]);
          break;
        case 58:
          yy.setClass($$[$0 - 1], $$[$0]);
          break;
        case 59:
          yy.setClass([$$[$0 - 2]], $$[$0]);
          break;
        case 60:
        case 62:
          this.$ = [$$[$0]];
          break;
        case 61:
        case 63:
          this.$ = $$[$0 - 2].concat([$$[$0]]);
          break;
        case 64:
          this.$ = $$[$0 - 2];
          yy.setCssStyle($$[$0 - 1], $$[$0]);
          break;
        case 65:
          this.$ = [$$[$0]];
          break;
        case 66:
          $$[$0 - 2].push($$[$0]);
          this.$ = $$[$0 - 2];
          break;
        case 68:
          this.$ = $$[$0 - 1] + $$[$0];
          break;
      }
    }, "anonymous"),
    table: [{ 3: 1, 4: 2, 6: $V0, 9: $V1, 11: $V2, 13: $V3 }, { 1: [3] }, { 3: 8, 4: 2, 5: [1, 7], 6: $V0, 9: $V1, 11: $V2, 13: $V3 }, { 5: [1, 9] }, { 10: [1, 10] }, { 12: [1, 11] }, o($V4, [2, 6]), { 3: 12, 4: 2, 6: $V0, 9: $V1, 11: $V2, 13: $V3 }, { 1: [2, 2] }, { 4: 17, 5: $V5, 7: 13, 8: $V6, 9: $V1, 11: $V2, 13: $V3, 14: 14, 15: 15, 16: 16, 17: 18, 18: 19, 19: 20, 20: 21, 21: $V7, 22: $V8, 23: $V9, 24: $Va, 25: 23, 33: 25, 41: $Vb, 42: $Vc, 43: $Vd, 44: $Ve, 45: $Vf, 46: $Vg, 54: $Vh, 72: $Vi, 74: $Vj, 77: $Vk, 89: $Vl, 90: $Vm }, o($V4, [2, 4]), o($V4, [2, 5]), { 1: [2, 1] }, { 8: [1, 41] }, { 4: 17, 5: $V5, 7: 42, 8: $V6, 9: $V1, 11: $V2, 13: $V3, 14: 14, 15: 15, 16: 16, 17: 18, 18: 19, 19: 20, 20: 21, 21: $V7, 22: $V8, 23: $V9, 24: $Va, 25: 23, 33: 25, 41: $Vb, 42: $Vc, 43: $Vd, 44: $Ve, 45: $Vf, 46: $Vg, 54: $Vh, 72: $Vi, 74: $Vj, 77: $Vk, 89: $Vl, 90: $Vm }, { 4: 17, 5: $V5, 7: 43, 8: $V6, 9: $V1, 11: $V2, 13: $V3, 14: 14, 15: 15, 16: 16, 17: 18, 18: 19, 19: 20, 20: 21, 21: $V7, 22: $V8, 23: $V9, 24: $Va, 25: 23, 33: 25, 41: $Vb, 42: $Vc, 43: $Vd, 44: $Ve, 45: $Vf, 46: $Vg, 54: $Vh, 72: $Vi, 74: $Vj, 77: $Vk, 89: $Vl, 90: $Vm }, { 4: 17, 5: $V5, 7: 44, 8: $V6, 9: $V1, 11: $V2, 13: $V3, 14: 14, 15: 15, 16: 16, 17: 18, 18: 19, 19: 20, 20: 21, 21: $V7, 22: $V8, 23: $V9, 24: $Va, 25: 23, 33: 25, 41: $Vb, 42: $Vc, 43: $Vd, 44: $Ve, 45: $Vf, 46: $Vg, 54: $Vh, 72: $Vi, 74: $Vj, 77: $Vk, 89: $Vl, 90: $Vm }, { 4: 17, 5: $V5, 7: 45, 8: $V6, 9: $V1, 11: $V2, 13: $V3, 14: 14, 15: 15, 16: 16, 17: 18, 18: 19, 19: 20, 20: 21, 21: $V7, 22: $V8, 23: $V9, 24: $Va, 25: 23, 33: 25, 41: $Vb, 42: $Vc, 43: $Vd, 44: $Ve, 45: $Vf, 46: $Vg, 54: $Vh, 72: $Vi, 74: $Vj, 77: $Vk, 89: $Vl, 90: $Vm }, { 4: 17, 5: $V5, 7: 46, 8: $V6, 9: $V1, 11: $V2, 13: $V3, 14: 14, 15: 15, 16: 16, 17: 18, 18: 19, 19: 20, 20: 21, 21: $V7, 22: $V8, 23: $V9, 24: $Va, 25: 23, 33: 25, 41: $Vb, 42: $Vc, 43: $Vd, 44: $Ve, 45: $Vf, 46: $Vg, 54: $Vh, 72: $Vi, 74: $Vj, 77: $Vk, 89: $Vl, 90: $Vm }, { 4: 17, 5: $V5, 7: 47, 8: $V6, 9: $V1, 11: $V2, 13: $V3, 14: 14, 15: 15, 16: 16, 17: 18, 18: 19, 19: 20, 20: 21, 21: $V7, 22: $V8, 23: $V9, 24: $Va, 25: 23, 33: 25, 41: $Vb, 42: $Vc, 43: $Vd, 44: $Ve, 45: $Vf, 46: $Vg, 54: $Vh, 72: $Vi, 74: $Vj, 77: $Vk, 89: $Vl, 90: $Vm }, { 4: 17, 5: $V5, 7: 48, 8: $V6, 9: $V1, 11: $V2, 13: $V3, 14: 14, 15: 15, 16: 16, 17: 18, 18: 19, 19: 20, 20: 21, 21: $V7, 22: $V8, 23: $V9, 24: $Va, 25: 23, 33: 25, 41: $Vb, 42: $Vc, 43: $Vd, 44: $Ve, 45: $Vf, 46: $Vg, 54: $Vh, 72: $Vi, 74: $Vj, 77: $Vk, 89: $Vl, 90: $Vm }, { 4: 17, 5: $V5, 7: 49, 8: $V6, 9: $V1, 11: $V2, 13: $V3, 14: 14, 15: 15, 16: 16, 17: 18, 18: 19, 19: 20, 20: 21, 21: $V7, 22: $V8, 23: $V9, 24: $Va, 25: 23, 33: 25, 41: $Vb, 42: $Vc, 43: $Vd, 44: $Ve, 45: $Vf, 46: $Vg, 54: $Vh, 72: $Vi, 74: $Vj, 77: $Vk, 89: $Vl, 90: $Vm }, { 4: 17, 5: $V5, 7: 50, 8: $V6, 9: $V1, 11: $V2, 13: $V3, 14: 14, 15: 15, 16: 16, 17: 18, 18: 19, 19: 20, 20: 21, 21: $V7, 22: $V8, 23: $V9, 24: $Va, 25: 23, 33: 25, 41: $Vb, 42: $Vc, 43: $Vd, 44: $Ve, 45: $Vf, 46: $Vg, 54: $Vh, 72: $Vi, 74: $Vj, 77: $Vk, 89: $Vl, 90: $Vm }, { 26: 51, 89: [1, 52], 90: [1, 53] }, { 55: 54, 89: [1, 55], 90: [1, 56] }, { 29: [1, 59], 61: [1, 57], 63: [1, 58] }, o($Vn, [2, 17]), o($Vn, [2, 18]), o($Vn, [2, 19]), o($Vn, [2, 20]), { 30: 60, 33: 62, 75: $Vo, 89: $Vl, 90: $Vm }, { 30: 63, 33: 62, 75: $Vo, 89: $Vl, 90: $Vm }, { 30: 64, 33: 62, 75: $Vo, 89: $Vl, 90: $Vm }, o($Vp, [2, 29]), o($Vp, [2, 30]), o($Vp, [2, 31]), o($Vp, [2, 32]), o($Vp, [2, 33]), o($Vp, [2, 34]), o($Vq, [2, 81]), o($Vq, [2, 82]), { 1: [2, 3] }, { 8: [2, 8] }, { 8: [2, 9] }, { 8: [2, 10] }, { 8: [2, 11] }, { 8: [2, 12] }, { 8: [2, 13] }, { 8: [2, 14] }, { 8: [2, 15] }, { 8: [2, 16] }, { 27: [1, 65], 29: [1, 66] }, o($Vr, [2, 79]), o($Vr, [2, 80]), { 27: [1, 67], 29: [1, 68] }, o($Vr, [2, 85]), o($Vr, [2, 86]), { 62: 69, 65: $Vs, 66: $Vt, 67: $Vu, 68: $Vv, 69: $Vw, 70: $Vx, 71: $Vy }, { 62: 77, 65: $Vs, 66: $Vt, 67: $Vu, 68: $Vv, 69: $Vw, 70: $Vx, 71: $Vy }, { 30: 78, 33: 62, 75: $Vo, 89: $Vl, 90: $Vm }, { 73: 79, 75: $Vz, 76: $VA, 78: 81, 79: 82, 80: $VB, 81: $VC, 82: $VD, 83: $VE, 84: $VF, 85: $VG, 86: $VH, 87: $VI, 88: $VJ }, o($VK, [2, 60]), o($VK, [2, 62]), { 73: 93, 75: $Vz, 76: $VA, 78: 81, 79: 82, 80: $VB, 81: $VC, 82: $VD, 83: $VE, 84: $VF, 85: $VG, 86: $VH, 87: $VI, 88: $VJ }, { 30: 94, 33: 62, 75: $Vo, 76: $VA, 89: $Vl, 90: $Vm }, { 5: [1, 95] }, { 30: 96, 33: 62, 75: $Vo, 89: $Vl, 90: $Vm }, { 5: [1, 97] }, { 30: 98, 33: 62, 75: $Vo, 89: $Vl, 90: $Vm }, { 63: [1, 99] }, o($VL, [2, 50]), o($VL, [2, 51]), o($VL, [2, 52]), o($VL, [2, 53]), o($VL, [2, 54]), o($VL, [2, 55]), o($VL, [2, 56]), { 64: [1, 100] }, o($Vn, [2, 59], { 76: $VA }), o($Vn, [2, 64], { 76: $VM }), { 33: 103, 75: [1, 102], 89: $Vl, 90: $Vm }, o($VN, [2, 65], { 79: 104, 75: $Vz, 80: $VB, 81: $VC, 82: $VD, 83: $VE, 84: $VF, 85: $VG, 86: $VH, 87: $VI, 88: $VJ }), o($VO, [2, 67]), o($VO, [2, 69]), o($VO, [2, 70]), o($VO, [2, 71]), o($VO, [2, 72]), o($VO, [2, 73]), o($VO, [2, 74]), o($VO, [2, 75]), o($VO, [2, 76]), o($VO, [2, 77]), o($VO, [2, 78]), o($Vn, [2, 57], { 76: $VM }), o($Vn, [2, 58], { 76: $VA }), { 5: $VP, 28: 105, 31: $VQ, 34: $VR, 36: $VS, 38: $VT, 40: $VU }, { 27: [1, 112], 76: $VA }, { 5: $VV, 40: $VW, 56: 113, 57: $VX, 59: $VY }, { 27: [1, 118], 76: $VA }, { 33: 119, 89: $Vl, 90: $Vm }, { 33: 120, 89: $Vl, 90: $Vm }, { 75: $Vz, 78: 121, 79: 82, 80: $VB, 81: $VC, 82: $VD, 83: $VE, 84: $VF, 85: $VG, 86: $VH, 87: $VI, 88: $VJ }, o($VK, [2, 61]), o($VK, [2, 63]), o($VO, [2, 68]), o($Vn, [2, 21]), { 32: [1, 122] }, { 32: [1, 123] }, { 32: [1, 124] }, { 32: [1, 125] }, { 5: $VP, 28: 126, 31: $VQ, 34: $VR, 36: $VS, 38: $VT, 40: $VU }, o($Vn, [2, 28]), { 5: [1, 127] }, o($Vn, [2, 42]), { 32: [1, 128] }, { 32: [1, 129] }, { 5: $VV, 40: $VW, 56: 130, 57: $VX, 59: $VY }, o($Vn, [2, 47]), { 5: [1, 131] }, o($Vn, [2, 48]), o($Vn, [2, 49]), o($VN, [2, 66], { 79: 104, 75: $Vz, 80: $VB, 81: $VC, 82: $VD, 83: $VE, 84: $VF, 85: $VG, 86: $VH, 87: $VI, 88: $VJ }), { 33: 132, 89: $Vl, 90: $Vm }, { 35: 133, 89: [1, 134], 90: [1, 135] }, { 37: 136, 47: [1, 137], 48: [1, 138], 49: [1, 139] }, { 39: 140, 50: [1, 141], 51: [1, 142], 52: [1, 143], 53: [1, 144] }, o($Vn, [2, 27]), { 5: $VP, 28: 145, 31: $VQ, 34: $VR, 36: $VS, 38: $VT, 40: $VU }, { 58: 146, 89: [1, 147], 90: [1, 148] }, { 60: 149, 89: [1, 150], 90: [1, 151] }, o($Vn, [2, 46]), { 5: $VV, 40: $VW, 56: 152, 57: $VX, 59: $VY }, { 5: [1, 153] }, { 5: [1, 154] }, { 5: [2, 83] }, { 5: [2, 84] }, { 5: [1, 155] }, { 5: [2, 35] }, { 5: [2, 36] }, { 5: [2, 37] }, { 5: [1, 156] }, { 5: [2, 38] }, { 5: [2, 39] }, { 5: [2, 40] }, { 5: [2, 41] }, o($Vn, [2, 22]), { 5: [1, 157] }, { 5: [2, 87] }, { 5: [2, 88] }, { 5: [1, 158] }, { 5: [2, 89] }, { 5: [2, 90] }, o($Vn, [2, 43]), { 5: $VP, 28: 159, 31: $VQ, 34: $VR, 36: $VS, 38: $VT, 40: $VU }, { 5: $VP, 28: 160, 31: $VQ, 34: $VR, 36: $VS, 38: $VT, 40: $VU }, { 5: $VP, 28: 161, 31: $VQ, 34: $VR, 36: $VS, 38: $VT, 40: $VU }, { 5: $VP, 28: 162, 31: $VQ, 34: $VR, 36: $VS, 38: $VT, 40: $VU }, { 5: $VV, 40: $VW, 56: 163, 57: $VX, 59: $VY }, { 5: $VV, 40: $VW, 56: 164, 57: $VX, 59: $VY }, o($Vn, [2, 23]), o($Vn, [2, 24]), o($Vn, [2, 25]), o($Vn, [2, 26]), o($Vn, [2, 44]), o($Vn, [2, 45])],
    defaultActions: { 8: [2, 2], 12: [2, 1], 41: [2, 3], 42: [2, 8], 43: [2, 9], 44: [2, 10], 45: [2, 11], 46: [2, 12], 47: [2, 13], 48: [2, 14], 49: [2, 15], 50: [2, 16], 134: [2, 83], 135: [2, 84], 137: [2, 35], 138: [2, 36], 139: [2, 37], 141: [2, 38], 142: [2, 39], 143: [2, 40], 144: [2, 41], 147: [2, 87], 148: [2, 88], 150: [2, 89], 151: [2, 90] },
    parseError: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function parseError(str, hash) {
      if (hash.recoverable) {
        this.trace(str);
      } else {
        var error = new Error(str);
        error.hash = hash;
        throw error;
      }
    }, "parseError"),
    parse: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function parse(input) {
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
      (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(popStack, "popStack");
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
      (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(lex, "lex");
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
  var lexer = /* @__PURE__ */ (function() {
    var lexer2 = {
      EOF: 1,
      parseError: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function parseError(str, hash) {
        if (this.yy.parser) {
          this.yy.parser.parseError(str, hash);
        } else {
          throw new Error(str);
        }
      }, "parseError"),
      // resets the lexer, sets new input
      setInput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function(input, yy) {
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
      input: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function() {
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
      unput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function(ch) {
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
      more: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function() {
        this._more = true;
        return this;
      }, "more"),
      // When called from action, signals the lexer that this rule fails to match the input, so the next matching rule (regex) should be tested instead.
      reject: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function() {
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
      less: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function(n) {
        this.unput(this.match.slice(n));
      }, "less"),
      // displays already matched input, i.e. for error messages
      pastInput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function() {
        var past = this.matched.substr(0, this.matched.length - this.match.length);
        return (past.length > 20 ? "..." : "") + past.substr(-20).replace(/\n/g, "");
      }, "pastInput"),
      // displays upcoming input, i.e. for error messages
      upcomingInput: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function() {
        var next = this.match;
        if (next.length < 20) {
          next += this._input.substr(0, 20 - next.length);
        }
        return (next.substr(0, 20) + (next.length > 20 ? "..." : "")).replace(/\n/g, "");
      }, "upcomingInput"),
      // displays the character position where the lexing error occurred, i.e. for error messages
      showPosition: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function() {
        var pre = this.pastInput();
        var c = new Array(pre.length + 1).join("-");
        return pre + this.upcomingInput() + "\n" + c + "^";
      }, "showPosition"),
      // test the lexed token: return FALSE when not a match, otherwise return token
      test_match: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function(match, indexed_rule) {
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
      next: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function() {
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
      lex: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function lex() {
        var r = this.next();
        if (r) {
          return r;
        } else {
          return this.lex();
        }
      }, "lex"),
      // activates a new lexer condition state (pushes the new lexer condition state onto the condition stack)
      begin: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function begin(condition) {
        this.conditionStack.push(condition);
      }, "begin"),
      // pop the previously active lexer condition state off the condition stack
      popState: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function popState() {
        var n = this.conditionStack.length - 1;
        if (n > 0) {
          return this.conditionStack.pop();
        } else {
          return this.conditionStack[0];
        }
      }, "popState"),
      // produce the lexer rule set which is active for the currently active lexer condition state
      _currentRules: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function _currentRules() {
        if (this.conditionStack.length && this.conditionStack[this.conditionStack.length - 1]) {
          return this.conditions[this.conditionStack[this.conditionStack.length - 1]].rules;
        } else {
          return this.conditions["INITIAL"].rules;
        }
      }, "_currentRules"),
      // return the currently active lexer condition state; when an index argument is provided it produces the N-th previous condition state, if available
      topState: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function topState(n) {
        n = this.conditionStack.length - 1 - Math.abs(n || 0);
        if (n >= 0) {
          return this.conditionStack[n];
        } else {
          return "INITIAL";
        }
      }, "topState"),
      // alias for begin(condition)
      pushState: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function pushState(condition) {
        this.begin(condition);
      }, "pushState"),
      // return the number of states currently on the stack
      stateStackSize: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function stateStackSize() {
        return this.conditionStack.length;
      }, "stateStackSize"),
      options: { "case-insensitive": true },
      performAction: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(function anonymous(yy, yy_, $avoiding_name_collisions, YY_START) {
        var YYSTATE = YY_START;
        switch ($avoiding_name_collisions) {
          case 0:
            return "title";
            break;
          case 1:
            this.begin("acc_title");
            return 9;
            break;
          case 2:
            this.popState();
            return "acc_title_value";
            break;
          case 3:
            this.begin("acc_descr");
            return 11;
            break;
          case 4:
            this.popState();
            return "acc_descr_value";
            break;
          case 5:
            this.begin("acc_descr_multiline");
            break;
          case 6:
            this.popState();
            break;
          case 7:
            return "acc_descr_multiline_value";
            break;
          case 8:
            return 21;
            break;
          case 9:
            return 22;
            break;
          case 10:
            return 23;
            break;
          case 11:
            return 24;
            break;
          case 12:
            return 5;
            break;
          case 13:
            break;
          case 14:
            break;
          case 15:
            break;
          case 16:
            return 8;
            break;
          case 17:
            return 6;
            break;
          case 18:
            return 27;
            break;
          case 19:
            return 40;
            break;
          case 20:
            return 29;
            break;
          case 21:
            return 32;
            break;
          case 22:
            return 31;
            break;
          case 23:
            return 34;
            break;
          case 24:
            return 36;
            break;
          case 25:
            return 38;
            break;
          case 26:
            return 41;
            break;
          case 27:
            return 42;
            break;
          case 28:
            return 43;
            break;
          case 29:
            return 44;
            break;
          case 30:
            return 45;
            break;
          case 31:
            return 46;
            break;
          case 32:
            return 47;
            break;
          case 33:
            return 48;
            break;
          case 34:
            return 49;
            break;
          case 35:
            return 50;
            break;
          case 36:
            return 51;
            break;
          case 37:
            return 52;
            break;
          case 38:
            return 53;
            break;
          case 39:
            return 54;
            break;
          case 40:
            return 65;
            break;
          case 41:
            return 66;
            break;
          case 42:
            return 67;
            break;
          case 43:
            return 68;
            break;
          case 44:
            return 69;
            break;
          case 45:
            return 70;
            break;
          case 46:
            return 71;
            break;
          case 47:
            return 57;
            break;
          case 48:
            return 59;
            break;
          case 49:
            this.begin("style");
            return 77;
            break;
          case 50:
            return 75;
            break;
          case 51:
            return 81;
            break;
          case 52:
            return 88;
            break;
          case 53:
            return "PERCENT";
            break;
          case 54:
            return 86;
            break;
          case 55:
            return 84;
            break;
          case 56:
            break;
          case 57:
            this.begin("string");
            break;
          case 58:
            this.popState();
            break;
          case 59:
            this.begin("style");
            return 72;
            break;
          case 60:
            this.begin("style");
            return 74;
            break;
          case 61:
            return 61;
            break;
          case 62:
            return 64;
            break;
          case 63:
            return 63;
            break;
          case 64:
            this.begin("string");
            break;
          case 65:
            this.popState();
            break;
          case 66:
            return "qString";
            break;
          case 67:
            yy_.yytext = yy_.yytext.trim();
            return 89;
            break;
          case 68:
            return 75;
            break;
          case 69:
            return 80;
            break;
          case 70:
            return 76;
            break;
        }
      }, "anonymous"),
      rules: [/^(?:title\s[^#\n;]+)/i, /^(?:accTitle\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*:\s*)/i, /^(?:(?!\n||)*[^\n]*)/i, /^(?:accDescr\s*\{\s*)/i, /^(?:[\}])/i, /^(?:[^\}]*)/i, /^(?:.*direction\s+TB[^\n]*)/i, /^(?:.*direction\s+BT[^\n]*)/i, /^(?:.*direction\s+RL[^\n]*)/i, /^(?:.*direction\s+LR[^\n]*)/i, /^(?:(\r?\n)+)/i, /^(?:\s+)/i, /^(?:#[^\n]*)/i, /^(?:%[^\n]*)/i, /^(?:$)/i, /^(?:requirementDiagram\b)/i, /^(?:\{)/i, /^(?:\})/i, /^(?::{3})/i, /^(?::)/i, /^(?:id\b)/i, /^(?:text\b)/i, /^(?:risk\b)/i, /^(?:verifyMethod\b)/i, /^(?:requirement\b)/i, /^(?:functionalRequirement\b)/i, /^(?:interfaceRequirement\b)/i, /^(?:performanceRequirement\b)/i, /^(?:physicalRequirement\b)/i, /^(?:designConstraint\b)/i, /^(?:low\b)/i, /^(?:medium\b)/i, /^(?:high\b)/i, /^(?:analysis\b)/i, /^(?:demonstration\b)/i, /^(?:inspection\b)/i, /^(?:test\b)/i, /^(?:element\b)/i, /^(?:contains\b)/i, /^(?:copies\b)/i, /^(?:derives\b)/i, /^(?:satisfies\b)/i, /^(?:verifies\b)/i, /^(?:refines\b)/i, /^(?:traces\b)/i, /^(?:type\b)/i, /^(?:docref\b)/i, /^(?:style\b)/i, /^(?:\w+)/i, /^(?::)/i, /^(?:;)/i, /^(?:%)/i, /^(?:-)/i, /^(?:#)/i, /^(?: )/i, /^(?:["])/i, /^(?:\n)/i, /^(?:classDef\b)/i, /^(?:class\b)/i, /^(?:<-)/i, /^(?:->)/i, /^(?:-)/i, /^(?:["])/i, /^(?:["])/i, /^(?:[^"]*)/i, /^(?:[\w][^:,\r\n\{\<\>\-\=]*)/i, /^(?:\w+)/i, /^(?:[0-9]+)/i, /^(?:,)/i],
      conditions: { "acc_descr_multiline": { "rules": [6, 7, 68, 69, 70], "inclusive": false }, "acc_descr": { "rules": [4, 68, 69, 70], "inclusive": false }, "acc_title": { "rules": [2, 68, 69, 70], "inclusive": false }, "style": { "rules": [50, 51, 52, 53, 54, 55, 56, 57, 58, 68, 69, 70], "inclusive": false }, "unqString": { "rules": [68, 69, 70], "inclusive": false }, "token": { "rules": [68, 69, 70], "inclusive": false }, "string": { "rules": [65, 66, 68, 69, 70], "inclusive": false }, "INITIAL": { "rules": [0, 1, 3, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 59, 60, 61, 62, 63, 64, 67, 68, 69, 70], "inclusive": true } }
    };
    return lexer2;
  })();
  parser2.lexer = lexer;
  function Parser() {
    this.yy = {};
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(Parser, "Parser");
  Parser.prototype = parser2;
  parser2.Parser = Parser;
  return new Parser();
})();
parser.parser = parser;
var requirementDiagram_default = parser;

// src/diagrams/requirement/requirementDb.ts
var RequirementDB = class {
  constructor() {
    this.relations = [];
    this.latestRequirement = this.getInitialRequirement();
    this.requirements = /* @__PURE__ */ new Map();
    this.latestElement = this.getInitialElement();
    this.elements = /* @__PURE__ */ new Map();
    this.classes = /* @__PURE__ */ new Map();
    this.direction = "TB";
    this.RequirementType = {
      REQUIREMENT: "Requirement",
      FUNCTIONAL_REQUIREMENT: "Functional Requirement",
      INTERFACE_REQUIREMENT: "Interface Requirement",
      PERFORMANCE_REQUIREMENT: "Performance Requirement",
      PHYSICAL_REQUIREMENT: "Physical Requirement",
      DESIGN_CONSTRAINT: "Design Constraint"
    };
    this.RiskLevel = {
      LOW_RISK: "Low",
      MED_RISK: "Medium",
      HIGH_RISK: "High"
    };
    this.VerifyType = {
      VERIFY_ANALYSIS: "Analysis",
      VERIFY_DEMONSTRATION: "Demonstration",
      VERIFY_INSPECTION: "Inspection",
      VERIFY_TEST: "Test"
    };
    this.Relationships = {
      CONTAINS: "contains",
      COPIES: "copies",
      DERIVES: "derives",
      SATISFIES: "satisfies",
      VERIFIES: "verifies",
      REFINES: "refines",
      TRACES: "traces"
    };
    this.setAccTitle = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_10__/* .setAccTitle */ .GN;
    this.getAccTitle = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_10__/* .getAccTitle */ .eu;
    this.setAccDescription = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_10__/* .setAccDescription */ .U$;
    this.getAccDescription = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_10__/* .getAccDescription */ .Mx;
    this.setDiagramTitle = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_10__/* .setDiagramTitle */ .g2;
    this.getDiagramTitle = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_10__/* .getDiagramTitle */ .Kr;
    this.getConfig = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(() => (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_10__/* .getConfig2 */ .nV)().requirement, "getConfig");
    this.clear();
    this.setDirection = this.setDirection.bind(this);
    this.addRequirement = this.addRequirement.bind(this);
    this.setNewReqId = this.setNewReqId.bind(this);
    this.setNewReqRisk = this.setNewReqRisk.bind(this);
    this.setNewReqText = this.setNewReqText.bind(this);
    this.setNewReqVerifyMethod = this.setNewReqVerifyMethod.bind(this);
    this.addElement = this.addElement.bind(this);
    this.setNewElementType = this.setNewElementType.bind(this);
    this.setNewElementDocRef = this.setNewElementDocRef.bind(this);
    this.addRelationship = this.addRelationship.bind(this);
    this.setCssStyle = this.setCssStyle.bind(this);
    this.setClass = this.setClass.bind(this);
    this.defineClass = this.defineClass.bind(this);
    this.setAccTitle = this.setAccTitle.bind(this);
    this.setAccDescription = this.setAccDescription.bind(this);
  }
  static {
    (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(this, "RequirementDB");
  }
  getDirection() {
    return this.direction;
  }
  setDirection(dir) {
    this.direction = dir;
  }
  resetLatestRequirement() {
    this.latestRequirement = this.getInitialRequirement();
  }
  resetLatestElement() {
    this.latestElement = this.getInitialElement();
  }
  getInitialRequirement() {
    return {
      requirementId: "",
      text: "",
      risk: "",
      verifyMethod: "",
      name: "",
      type: "",
      cssStyles: [],
      classes: ["default"]
    };
  }
  getInitialElement() {
    return {
      name: "",
      type: "",
      docRef: "",
      cssStyles: [],
      classes: ["default"]
    };
  }
  addRequirement(name, type) {
    if (!this.requirements.has(name)) {
      this.requirements.set(name, {
        name,
        type,
        requirementId: this.latestRequirement.requirementId,
        text: this.latestRequirement.text,
        risk: this.latestRequirement.risk,
        verifyMethod: this.latestRequirement.verifyMethod,
        cssStyles: [],
        classes: ["default"]
      });
    }
    this.resetLatestRequirement();
    return this.requirements.get(name);
  }
  getRequirements() {
    return this.requirements;
  }
  setNewReqId(id) {
    if (this.latestRequirement !== void 0) {
      this.latestRequirement.requirementId = id;
    }
  }
  setNewReqText(text) {
    if (this.latestRequirement !== void 0) {
      this.latestRequirement.text = text;
    }
  }
  setNewReqRisk(risk) {
    if (this.latestRequirement !== void 0) {
      this.latestRequirement.risk = risk;
    }
  }
  setNewReqVerifyMethod(verifyMethod) {
    if (this.latestRequirement !== void 0) {
      this.latestRequirement.verifyMethod = verifyMethod;
    }
  }
  addElement(name) {
    if (!this.elements.has(name)) {
      this.elements.set(name, {
        name,
        type: this.latestElement.type,
        docRef: this.latestElement.docRef,
        cssStyles: [],
        classes: ["default"]
      });
      _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .log */ .cM.info("Added new element: ", name);
    }
    this.resetLatestElement();
    return this.elements.get(name);
  }
  getElements() {
    return this.elements;
  }
  setNewElementType(type) {
    if (this.latestElement !== void 0) {
      this.latestElement.type = type;
    }
  }
  setNewElementDocRef(docRef) {
    if (this.latestElement !== void 0) {
      this.latestElement.docRef = docRef;
    }
  }
  addRelationship(type, src, dst) {
    this.relations.push({
      type,
      src,
      dst
    });
  }
  getRelationships() {
    return this.relations;
  }
  clear() {
    this.relations = [];
    this.resetLatestRequirement();
    this.requirements = /* @__PURE__ */ new Map();
    this.resetLatestElement();
    this.elements = /* @__PURE__ */ new Map();
    this.classes = /* @__PURE__ */ new Map();
    (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_10__/* .clear */ .ZH)();
  }
  setCssStyle(ids, styles) {
    for (const id of ids) {
      const node = this.requirements.get(id) ?? this.elements.get(id);
      if (!styles || !node) {
        return;
      }
      for (const s of styles) {
        if (s.includes(",")) {
          node.cssStyles.push(...s.split(","));
        } else {
          node.cssStyles.push(s);
        }
      }
    }
  }
  setClass(ids, classNames) {
    for (const id of ids) {
      const node = this.requirements.get(id) ?? this.elements.get(id);
      if (node) {
        for (const _class of classNames) {
          node.classes.push(_class);
          const styles = this.classes.get(_class)?.styles;
          if (styles) {
            node.cssStyles.push(...styles);
          }
        }
      }
    }
  }
  defineClass(ids, style) {
    for (const id of ids) {
      let styleClass = this.classes.get(id);
      if (styleClass === void 0) {
        styleClass = { id, styles: [], textStyles: [] };
        this.classes.set(id, styleClass);
      }
      if (style) {
        style.forEach(function(s) {
          if (/color/.exec(s)) {
            const newStyle = s.replace("fill", "bgFill");
            styleClass.textStyles.push(newStyle);
          }
          styleClass.styles.push(s);
        });
      }
      this.requirements.forEach((value) => {
        if (value.classes.includes(id)) {
          value.cssStyles.push(...style.flatMap((s) => s.split(",")));
        }
      });
      this.elements.forEach((value) => {
        if (value.classes.includes(id)) {
          value.cssStyles.push(...style.flatMap((s) => s.split(",")));
        }
      });
    }
  }
  getClasses() {
    return this.classes;
  }
  getData() {
    const config = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_10__/* .getConfig2 */ .nV)();
    const nodes = [];
    const edges = [];
    for (const requirement of this.requirements.values()) {
      const node = requirement;
      node.id = requirement.name;
      node.cssStyles = requirement.cssStyles;
      node.cssClasses = requirement.classes.join(" ");
      node.shape = "requirementBox";
      node.look = config.look;
      nodes.push(node);
    }
    for (const element of this.elements.values()) {
      const node = element;
      node.shape = "requirementBox";
      node.look = config.look;
      node.id = element.name;
      node.cssStyles = element.cssStyles;
      node.cssClasses = element.classes.join(" ");
      nodes.push(node);
    }
    for (const relation of this.relations) {
      let counter = 0;
      const isContains = relation.type === this.Relationships.CONTAINS;
      const edge = {
        id: `${relation.src}-${relation.dst}-${counter}`,
        start: this.requirements.get(relation.src)?.name ?? this.elements.get(relation.src)?.name,
        end: this.requirements.get(relation.dst)?.name ?? this.elements.get(relation.dst)?.name,
        label: `&lt;&lt;${relation.type}&gt;&gt;`,
        classes: "relationshipLine",
        style: ["fill:none", isContains ? "" : "stroke-dasharray: 10,7"],
        labelpos: "c",
        thickness: "normal",
        type: "normal",
        pattern: isContains ? "normal" : "dashed",
        arrowTypeStart: isContains ? "requirement_contains" : "",
        arrowTypeEnd: isContains ? "" : "requirement_arrow",
        look: config.look
      };
      edges.push(edge);
      counter++;
    }
    return { nodes, edges, other: {}, config, direction: this.getDirection() };
  }
};

// src/diagrams/requirement/styles.js
var getStyles = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)((options) => `

  marker {
    fill: ${options.relationColor};
    stroke: ${options.relationColor};
  }

  marker.cross {
    stroke: ${options.lineColor};
  }

  svg {
    font-family: ${options.fontFamily};
    font-size: ${options.fontSize};
  }

  .reqBox {
    fill: ${options.requirementBackground};
    fill-opacity: 1.0;
    stroke: ${options.requirementBorderColor};
    stroke-width: ${options.requirementBorderSize};
  }
  
  .reqTitle, .reqLabel{
    fill:  ${options.requirementTextColor};
  }
  .reqLabelBox {
    fill: ${options.relationLabelBackground};
    fill-opacity: 1.0;
  }

  .req-title-line {
    stroke: ${options.requirementBorderColor};
    stroke-width: ${options.requirementBorderSize};
  }
  .relationshipLine {
    stroke: ${options.relationColor};
    stroke-width: 1;
  }
  .relationshipLabel {
    fill: ${options.relationLabelColor};
  }
  .divider {
    stroke: ${options.nodeBorder};
    stroke-width: 1;
  }
  .label {
    font-family: ${options.fontFamily};
    color: ${options.nodeTextColor || options.textColor};
  }
  .label text,span {
    fill: ${options.nodeTextColor || options.textColor};
    color: ${options.nodeTextColor || options.textColor};
  }
  .labelBkg {
    background-color: ${options.edgeLabelBackground};
  }

`, "getStyles");
var styles_default = getStyles;

// src/diagrams/requirement/requirementRenderer.ts
var requirementRenderer_exports = {};
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__export */ .r2)(requirementRenderer_exports, {
  draw: () => draw
});
var draw = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .__name */ .eW)(async function(text, id, _version, diag) {
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .log */ .cM.info("REF0:");
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_11__/* .log */ .cM.info("Drawing requirement diagram (unified)", id);
  const { securityLevel, state: conf, layout } = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_10__/* .getConfig2 */ .nV)();
  const data4Layout = diag.db.getData();
  const svg = (0,_chunk_55IACEB6_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getDiagramElement */ .q)(id, securityLevel);
  data4Layout.type = diag.type;
  data4Layout.layoutAlgorithm = (0,_chunk_N4CR4FBY_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getRegisteredLayoutAlgorithm */ ._b)(layout);
  data4Layout.nodeSpacing = conf?.nodeSpacing ?? 50;
  data4Layout.rankSpacing = conf?.rankSpacing ?? 50;
  data4Layout.markers = ["requirement_contains", "requirement_arrow"];
  data4Layout.diagramId = id;
  await (0,_chunk_N4CR4FBY_mjs__WEBPACK_IMPORTED_MODULE_2__/* .render */ .sY)(data4Layout, svg);
  const padding = 8;
  _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_9__/* .utils_default */ .w8.insertTitle(
    svg,
    "requirementDiagramTitleText",
    conf?.titleTopMargin ?? 25,
    diag.db.getDiagramTitle()
  );
  (0,_chunk_QN33PNHL_mjs__WEBPACK_IMPORTED_MODULE_1__/* .setupViewPortForSVG */ .j)(svg, padding, "requirementDiagram", conf?.useMaxWidth ?? true);
}, "draw");

// src/diagrams/requirement/requirementDiagram.ts
var diagram = {
  parser: requirementDiagram_default,
  get db() {
    return new RequirementDB();
  },
  renderer: requirementRenderer_exports,
  styles: styles_default
};



/***/ })

}]);
//# sourceMappingURL=9848.558310b88143708c53d4.js.map?v=558310b88143708c53d4
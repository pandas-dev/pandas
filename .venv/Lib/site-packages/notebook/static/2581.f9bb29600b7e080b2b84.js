"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[2581],{

/***/ 80910:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   A: () => (/* binding */ populateCommonDb)
/* harmony export */ });
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(86906);


// src/diagrams/common/populateCommonDb.ts
function populateCommonDb(ast, db) {
  if (ast.accDescr) {
    db.setAccDescription?.(ast.accDescr);
  }
  if (ast.accTitle) {
    db.setAccTitle?.(ast.accTitle);
  }
  if (ast.title) {
    db.setDiagramTitle?.(ast.title);
  }
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(populateCommonDb, "populateCommonDb");




/***/ }),

/***/ 12581:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _chunk_4BMEZGHF_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(80910);
/* harmony import */ var _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(42626);
/* harmony import */ var _chunk_7B677QYD_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(30103);
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(86906);
/* harmony import */ var _mermaid_js_parser__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(13197);





// src/diagrams/packet/db.ts
var defaultPacketData = {
  packet: []
};
var data = structuredClone(defaultPacketData);
var DEFAULT_PACKET_CONFIG = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .defaultConfig_default */ .vZ.packet;
var getConfig2 = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => {
  const config = (0,_chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_1__/* .cleanAndMerge */ .Rb)({
    ...DEFAULT_PACKET_CONFIG,
    ...(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig */ .iE)().packet
  });
  if (config.showBits) {
    config.paddingY += 10;
  }
  return config;
}, "getConfig");
var getPacket = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => data.packet, "getPacket");
var pushWord = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((word) => {
  if (word.length > 0) {
    data.packet.push(word);
  }
}, "pushWord");
var clear2 = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(() => {
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .clear */ .ZH)();
  data = structuredClone(defaultPacketData);
}, "clear");
var db = {
  pushWord,
  getPacket,
  getConfig: getConfig2,
  clear: clear2,
  setAccTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccTitle */ .GN,
  getAccTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccTitle */ .eu,
  setDiagramTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setDiagramTitle */ .g2,
  getDiagramTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getDiagramTitle */ .Kr,
  getAccDescription: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccDescription */ .Mx,
  setAccDescription: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccDescription */ .U$
};

// src/diagrams/packet/parser.ts

var maxPacketSize = 1e4;
var populate = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((ast) => {
  (0,_chunk_4BMEZGHF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .populateCommonDb */ .A)(ast, db);
  let lastByte = -1;
  let word = [];
  let row = 1;
  const { bitsPerRow } = db.getConfig();
  for (let { start, end, label } of ast.blocks) {
    if (end && end < start) {
      throw new Error(`Packet block ${start} - ${end} is invalid. End must be greater than start.`);
    }
    if (start !== lastByte + 1) {
      throw new Error(
        `Packet block ${start} - ${end ?? start} is not contiguous. It should start from ${lastByte + 1}.`
      );
    }
    lastByte = end ?? start;
    _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug(`Packet block ${start} - ${lastByte} with label ${label}`);
    while (word.length <= bitsPerRow + 1 && db.getPacket().length < maxPacketSize) {
      const [block, nextBlock] = getNextFittingBlock({ start, end, label }, row, bitsPerRow);
      word.push(block);
      if (block.end + 1 === row * bitsPerRow) {
        db.pushWord(word);
        word = [];
        row++;
      }
      if (!nextBlock) {
        break;
      }
      ({ start, end, label } = nextBlock);
    }
  }
  db.pushWord(word);
}, "populate");
var getNextFittingBlock = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((block, row, bitsPerRow) => {
  if (block.end === void 0) {
    block.end = block.start;
  }
  if (block.start > block.end) {
    throw new Error(`Block start ${block.start} is greater than block end ${block.end}.`);
  }
  if (block.end + 1 <= row * bitsPerRow) {
    return [block, void 0];
  }
  return [
    {
      start: block.start,
      end: row * bitsPerRow - 1,
      label: block.label
    },
    {
      start: row * bitsPerRow,
      end: block.end,
      label: block.label
    }
  ];
}, "getNextFittingBlock");
var parser = {
  parse: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(async (input) => {
    const ast = await (0,_mermaid_js_parser__WEBPACK_IMPORTED_MODULE_4__/* .parse */ .Qc)("packet", input);
    _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .log */ .cM.debug(ast);
    populate(ast);
  }, "parse")
};

// src/diagrams/packet/renderer.ts
var draw = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((_text, id, _version, diagram2) => {
  const db2 = diagram2.db;
  const config = db2.getConfig();
  const { rowHeight, paddingY, bitWidth, bitsPerRow } = config;
  const words = db2.getPacket();
  const title = db2.getDiagramTitle();
  const totalRowHeight = rowHeight + paddingY;
  const svgHeight = totalRowHeight * (words.length + 1) - (title ? 0 : rowHeight);
  const svgWidth = bitWidth * bitsPerRow + 2;
  const svg = (0,_chunk_7B677QYD_mjs__WEBPACK_IMPORTED_MODULE_2__/* .selectSvgElement */ .P)(id);
  svg.attr("viewbox", `0 0 ${svgWidth} ${svgHeight}`);
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .configureSvgSize */ .v2)(svg, svgHeight, svgWidth, config.useMaxWidth);
  for (const [word, packet] of words.entries()) {
    drawWord(svg, packet, word, config);
  }
  svg.append("text").text(title).attr("x", svgWidth / 2).attr("y", svgHeight - totalRowHeight / 2).attr("dominant-baseline", "middle").attr("text-anchor", "middle").attr("class", "packetTitle");
}, "draw");
var drawWord = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)((svg, word, rowNumber, { rowHeight, paddingX, paddingY, bitWidth, bitsPerRow, showBits }) => {
  const group = svg.append("g");
  const wordY = rowNumber * (rowHeight + paddingY) + paddingY;
  for (const block of word) {
    const blockX = block.start % bitsPerRow * bitWidth + 1;
    const width = (block.end - block.start + 1) * bitWidth - paddingX;
    group.append("rect").attr("x", blockX).attr("y", wordY).attr("width", width).attr("height", rowHeight).attr("class", "packetBlock");
    group.append("text").attr("x", blockX + width / 2).attr("y", wordY + rowHeight / 2).attr("class", "packetLabel").attr("dominant-baseline", "middle").attr("text-anchor", "middle").text(block.label);
    if (!showBits) {
      continue;
    }
    const isSingleBlock = block.end === block.start;
    const bitNumberY = wordY - 2;
    group.append("text").attr("x", blockX + (isSingleBlock ? width / 2 : 0)).attr("y", bitNumberY).attr("class", "packetByte start").attr("dominant-baseline", "auto").attr("text-anchor", isSingleBlock ? "middle" : "start").text(block.start);
    if (!isSingleBlock) {
      group.append("text").attr("x", blockX + width).attr("y", bitNumberY).attr("class", "packetByte end").attr("dominant-baseline", "auto").attr("text-anchor", "end").text(block.end);
    }
  }
}, "drawWord");
var renderer = { draw };

// src/diagrams/packet/styles.ts
var defaultPacketStyleOptions = {
  byteFontSize: "10px",
  startByteColor: "black",
  endByteColor: "black",
  labelColor: "black",
  labelFontSize: "12px",
  titleColor: "black",
  titleFontSize: "14px",
  blockStrokeColor: "black",
  blockStrokeWidth: "1",
  blockFillColor: "#efefef"
};
var styles = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_3__/* .__name */ .eW)(({ packet } = {}) => {
  const options = (0,_chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_1__/* .cleanAndMerge */ .Rb)(defaultPacketStyleOptions, packet);
  return `
	.packetByte {
		font-size: ${options.byteFontSize};
	}
	.packetByte.start {
		fill: ${options.startByteColor};
	}
	.packetByte.end {
		fill: ${options.endByteColor};
	}
	.packetLabel {
		fill: ${options.labelColor};
		font-size: ${options.labelFontSize};
	}
	.packetTitle {
		fill: ${options.titleColor};
		font-size: ${options.titleFontSize};
	}
	.packetBlock {
		stroke: ${options.blockStrokeColor};
		stroke-width: ${options.blockStrokeWidth};
		fill: ${options.blockFillColor};
	}
	`;
}, "styles");

// src/diagrams/packet/diagram.ts
var diagram = {
  parser,
  db,
  renderer,
  styles
};



/***/ })

}]);
//# sourceMappingURL=2581.f9bb29600b7e080b2b84.js.map?v=f9bb29600b7e080b2b84
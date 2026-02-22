"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[5530],{

/***/ 18556:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   A: () => (/* binding */ populateCommonDb)
/* harmony export */ });
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(74999);


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
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(populateCommonDb, "populateCommonDb");




/***/ }),

/***/ 25530:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _chunk_EXTU4WIE_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(98154);
/* harmony import */ var _chunk_4BX2VUAB_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(18556);
/* harmony import */ var _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(17175);
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(74999);
/* harmony import */ var _mermaid_js_parser__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(13197);






// src/diagrams/packet/db.ts
var DEFAULT_PACKET_CONFIG = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .defaultConfig_default */ .vZ.packet;
var PacketDB = class {
  constructor() {
    this.packet = [];
    this.setAccTitle = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccTitle */ .GN;
    this.getAccTitle = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccTitle */ .eu;
    this.setDiagramTitle = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setDiagramTitle */ .g2;
    this.getDiagramTitle = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getDiagramTitle */ .Kr;
    this.getAccDescription = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getAccDescription */ .Mx;
    this.setAccDescription = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .setAccDescription */ .U$;
  }
  static {
    (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)(this, "PacketDB");
  }
  getConfig() {
    const config = (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_2__/* .cleanAndMerge */ .Rb)({
      ...DEFAULT_PACKET_CONFIG,
      ...(0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig */ .iE)().packet
    });
    if (config.showBits) {
      config.paddingY += 10;
    }
    return config;
  }
  getPacket() {
    return this.packet;
  }
  pushWord(word) {
    if (word.length > 0) {
      this.packet.push(word);
    }
  }
  clear() {
    (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .clear */ .ZH)();
    this.packet = [];
  }
};

// src/diagrams/packet/parser.ts

var maxPacketSize = 1e4;
var populate = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)((ast, db) => {
  (0,_chunk_4BX2VUAB_mjs__WEBPACK_IMPORTED_MODULE_1__/* .populateCommonDb */ .A)(ast, db);
  let lastBit = -1;
  let word = [];
  let row = 1;
  const { bitsPerRow } = db.getConfig();
  for (let { start, end, bits, label } of ast.blocks) {
    if (start !== void 0 && end !== void 0 && end < start) {
      throw new Error(`Packet block ${start} - ${end} is invalid. End must be greater than start.`);
    }
    start ??= lastBit + 1;
    if (start !== lastBit + 1) {
      throw new Error(
        `Packet block ${start} - ${end ?? start} is not contiguous. It should start from ${lastBit + 1}.`
      );
    }
    if (bits === 0) {
      throw new Error(`Packet block ${start} is invalid. Cannot have a zero bit field.`);
    }
    end ??= start + (bits ?? 1) - 1;
    bits ??= end - start + 1;
    lastBit = end;
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .log */ .cM.debug(`Packet block ${start} - ${lastBit} with label ${label}`);
    while (word.length <= bitsPerRow + 1 && db.getPacket().length < maxPacketSize) {
      const [block, nextBlock] = getNextFittingBlock({ start, end, bits, label }, row, bitsPerRow);
      word.push(block);
      if (block.end + 1 === row * bitsPerRow) {
        db.pushWord(word);
        word = [];
        row++;
      }
      if (!nextBlock) {
        break;
      }
      ({ start, end, bits, label } = nextBlock);
    }
  }
  db.pushWord(word);
}, "populate");
var getNextFittingBlock = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)((block, row, bitsPerRow) => {
  if (block.start === void 0) {
    throw new Error("start should have been set during first phase");
  }
  if (block.end === void 0) {
    throw new Error("end should have been set during first phase");
  }
  if (block.start > block.end) {
    throw new Error(`Block start ${block.start} is greater than block end ${block.end}.`);
  }
  if (block.end + 1 <= row * bitsPerRow) {
    return [block, void 0];
  }
  const rowEnd = row * bitsPerRow - 1;
  const rowStart = row * bitsPerRow;
  return [
    {
      start: block.start,
      end: rowEnd,
      label: block.label,
      bits: rowEnd - block.start
    },
    {
      start: rowStart,
      end: block.end,
      label: block.label,
      bits: block.end - rowStart
    }
  ];
}, "getNextFittingBlock");
var parser = {
  // @ts-expect-error - PacketDB is not assignable to DiagramDB
  parser: { yy: void 0 },
  parse: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)(async (input) => {
    const ast = await (0,_mermaid_js_parser__WEBPACK_IMPORTED_MODULE_5__/* .parse */ .Qc)("packet", input);
    const db = parser.parser?.yy;
    if (!(db instanceof PacketDB)) {
      throw new Error(
        "parser.parser?.yy was not a PacketDB. This is due to a bug within Mermaid, please report this issue at https://github.com/mermaid-js/mermaid/issues."
      );
    }
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .log */ .cM.debug(ast);
    populate(ast, db);
  }, "parse")
};

// src/diagrams/packet/renderer.ts
var draw = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)((_text, id, _version, diagram2) => {
  const db = diagram2.db;
  const config = db.getConfig();
  const { rowHeight, paddingY, bitWidth, bitsPerRow } = config;
  const words = db.getPacket();
  const title = db.getDiagramTitle();
  const totalRowHeight = rowHeight + paddingY;
  const svgHeight = totalRowHeight * (words.length + 1) - (title ? 0 : rowHeight);
  const svgWidth = bitWidth * bitsPerRow + 2;
  const svg = (0,_chunk_EXTU4WIE_mjs__WEBPACK_IMPORTED_MODULE_0__/* .selectSvgElement */ .P)(id);
  svg.attr("viewbox", `0 0 ${svgWidth} ${svgHeight}`);
  (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .configureSvgSize */ .v2)(svg, svgHeight, svgWidth, config.useMaxWidth);
  for (const [word, packet] of words.entries()) {
    drawWord(svg, packet, word, config);
  }
  svg.append("text").text(title).attr("x", svgWidth / 2).attr("y", svgHeight - totalRowHeight / 2).attr("dominant-baseline", "middle").attr("text-anchor", "middle").attr("class", "packetTitle");
}, "draw");
var drawWord = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)((svg, word, rowNumber, { rowHeight, paddingX, paddingY, bitWidth, bitsPerRow, showBits }) => {
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
var styles = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)(({ packet } = {}) => {
  const options = (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_2__/* .cleanAndMerge */ .Rb)(defaultPacketStyleOptions, packet);
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
  get db() {
    return new PacketDB();
  },
  renderer,
  styles
};



/***/ })

}]);
//# sourceMappingURL=5530.8eb3482278bcfcf70e4a.js.map?v=8eb3482278bcfcf70e4a
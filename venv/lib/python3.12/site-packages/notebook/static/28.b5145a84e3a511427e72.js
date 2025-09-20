"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[28],{

/***/ 21605:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.newState = exports.STATE = exports.AbstractMathItem = exports.protoItem = void 0;
function protoItem(open, math, close, n, start, end, display) {
    if (display === void 0) { display = null; }
    var item = { open: open, math: math, close: close,
        n: n, start: { n: start }, end: { n: end }, display: display };
    return item;
}
exports.protoItem = protoItem;
var AbstractMathItem = (function () {
    function AbstractMathItem(math, jax, display, start, end) {
        if (display === void 0) { display = true; }
        if (start === void 0) { start = { i: 0, n: 0, delim: '' }; }
        if (end === void 0) { end = { i: 0, n: 0, delim: '' }; }
        this.root = null;
        this.typesetRoot = null;
        this.metrics = {};
        this.inputData = {};
        this.outputData = {};
        this._state = exports.STATE.UNPROCESSED;
        this.math = math;
        this.inputJax = jax;
        this.display = display;
        this.start = start;
        this.end = end;
        this.root = null;
        this.typesetRoot = null;
        this.metrics = {};
        this.inputData = {};
        this.outputData = {};
    }
    Object.defineProperty(AbstractMathItem.prototype, "isEscaped", {
        get: function () {
            return this.display === null;
        },
        enumerable: false,
        configurable: true
    });
    AbstractMathItem.prototype.render = function (document) {
        document.renderActions.renderMath(this, document);
    };
    AbstractMathItem.prototype.rerender = function (document, start) {
        if (start === void 0) { start = exports.STATE.RERENDER; }
        if (this.state() >= start) {
            this.state(start - 1);
        }
        document.renderActions.renderMath(this, document, start);
    };
    AbstractMathItem.prototype.convert = function (document, end) {
        if (end === void 0) { end = exports.STATE.LAST; }
        document.renderActions.renderConvert(this, document, end);
    };
    AbstractMathItem.prototype.compile = function (document) {
        if (this.state() < exports.STATE.COMPILED) {
            this.root = this.inputJax.compile(this, document);
            this.state(exports.STATE.COMPILED);
        }
    };
    AbstractMathItem.prototype.typeset = function (document) {
        if (this.state() < exports.STATE.TYPESET) {
            this.typesetRoot = document.outputJax[this.isEscaped ? 'escaped' : 'typeset'](this, document);
            this.state(exports.STATE.TYPESET);
        }
    };
    AbstractMathItem.prototype.updateDocument = function (_document) { };
    AbstractMathItem.prototype.removeFromDocument = function (_restore) {
        if (_restore === void 0) { _restore = false; }
    };
    AbstractMathItem.prototype.setMetrics = function (em, ex, cwidth, lwidth, scale) {
        this.metrics = {
            em: em, ex: ex,
            containerWidth: cwidth,
            lineWidth: lwidth,
            scale: scale
        };
    };
    AbstractMathItem.prototype.state = function (state, restore) {
        if (state === void 0) { state = null; }
        if (restore === void 0) { restore = false; }
        if (state != null) {
            if (state < exports.STATE.INSERTED && this._state >= exports.STATE.INSERTED) {
                this.removeFromDocument(restore);
            }
            if (state < exports.STATE.TYPESET && this._state >= exports.STATE.TYPESET) {
                this.outputData = {};
            }
            if (state < exports.STATE.COMPILED && this._state >= exports.STATE.COMPILED) {
                this.inputData = {};
            }
            this._state = state;
        }
        return this._state;
    };
    AbstractMathItem.prototype.reset = function (restore) {
        if (restore === void 0) { restore = false; }
        this.state(exports.STATE.UNPROCESSED, restore);
    };
    return AbstractMathItem;
}());
exports.AbstractMathItem = AbstractMathItem;
exports.STATE = {
    UNPROCESSED: 0,
    FINDMATH: 10,
    COMPILED: 20,
    CONVERT: 100,
    METRICS: 110,
    RERENDER: 125,
    TYPESET: 150,
    INSERTED: 200,
    LAST: 10000
};
function newState(name, state) {
    if (name in exports.STATE) {
        throw Error('State ' + name + ' already exists');
    }
    exports.STATE[name] = state;
}
exports.newState = newState;
//# sourceMappingURL=MathItem.js.map

/***/ }),

/***/ 47578:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.TeXAtom = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var mo_js_1 = __webpack_require__(19625);
var TeXAtom = (function (_super) {
    __extends(TeXAtom, _super);
    function TeXAtom(factory, attributes, children) {
        var _this = _super.call(this, factory, attributes, children) || this;
        _this.texclass = MmlNode_js_1.TEXCLASS.ORD;
        _this.setProperty('texClass', _this.texClass);
        return _this;
    }
    Object.defineProperty(TeXAtom.prototype, "kind", {
        get: function () {
            return 'TeXAtom';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(TeXAtom.prototype, "arity", {
        get: function () {
            return -1;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(TeXAtom.prototype, "notParent", {
        get: function () {
            return this.childNodes[0] && this.childNodes[0].childNodes.length === 1;
        },
        enumerable: false,
        configurable: true
    });
    TeXAtom.prototype.setTeXclass = function (prev) {
        this.childNodes[0].setTeXclass(null);
        return this.adjustTeXclass(prev);
    };
    TeXAtom.prototype.adjustTeXclass = function (prev) {
        return prev;
    };
    TeXAtom.defaults = __assign({}, MmlNode_js_1.AbstractMmlBaseNode.defaults);
    return TeXAtom;
}(MmlNode_js_1.AbstractMmlBaseNode));
exports.TeXAtom = TeXAtom;
TeXAtom.prototype.adjustTeXclass = mo_js_1.MmlMo.prototype.adjustTeXclass;
//# sourceMappingURL=TeXAtom.js.map

/***/ }),

/***/ 63142:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMaction = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMaction = (function (_super) {
    __extends(MmlMaction, _super);
    function MmlMaction() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMaction.prototype, "kind", {
        get: function () {
            return 'maction';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMaction.prototype, "arity", {
        get: function () {
            return 1;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMaction.prototype, "selected", {
        get: function () {
            var selection = this.attributes.get('selection');
            var i = Math.max(1, Math.min(this.childNodes.length, selection)) - 1;
            return this.childNodes[i] || this.factory.create('mrow');
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMaction.prototype, "isEmbellished", {
        get: function () {
            return this.selected.isEmbellished;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMaction.prototype, "isSpacelike", {
        get: function () {
            return this.selected.isSpacelike;
        },
        enumerable: false,
        configurable: true
    });
    MmlMaction.prototype.core = function () {
        return this.selected.core();
    };
    MmlMaction.prototype.coreMO = function () {
        return this.selected.coreMO();
    };
    MmlMaction.prototype.verifyAttributes = function (options) {
        _super.prototype.verifyAttributes.call(this, options);
        if (this.attributes.get('actiontype') !== 'toggle' &&
            this.attributes.getExplicit('selection') !== undefined) {
            var attributes = this.attributes.getAllAttributes();
            delete attributes.selection;
        }
    };
    MmlMaction.prototype.setTeXclass = function (prev) {
        if (this.attributes.get('actiontype') === 'tooltip' && this.childNodes[1]) {
            this.childNodes[1].setTeXclass(null);
        }
        var selected = this.selected;
        prev = selected.setTeXclass(prev);
        this.updateTeXclass(selected);
        return prev;
    };
    MmlMaction.prototype.nextToggleSelection = function () {
        var selection = Math.max(1, this.attributes.get('selection') + 1);
        if (selection > this.childNodes.length) {
            selection = 1;
        }
        this.attributes.set('selection', selection);
    };
    MmlMaction.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlNode.defaults), { actiontype: 'toggle', selection: 1 });
    return MmlMaction;
}(MmlNode_js_1.AbstractMmlNode));
exports.MmlMaction = MmlMaction;
//# sourceMappingURL=maction.js.map

/***/ }),

/***/ 1334:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMath = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMath = (function (_super) {
    __extends(MmlMath, _super);
    function MmlMath() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMath.prototype, "kind", {
        get: function () {
            return 'math';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMath.prototype, "linebreakContainer", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    MmlMath.prototype.setChildInheritedAttributes = function (attributes, display, level, prime) {
        if (this.attributes.get('mode') === 'display') {
            this.attributes.setInherited('display', 'block');
        }
        attributes = this.addInheritedAttributes(attributes, this.attributes.getAllAttributes());
        display = (!!this.attributes.get('displaystyle') ||
            (!this.attributes.get('displaystyle') && this.attributes.get('display') === 'block'));
        this.attributes.setInherited('displaystyle', display);
        level = (this.attributes.get('scriptlevel') ||
            this.constructor.defaults['scriptlevel']);
        _super.prototype.setChildInheritedAttributes.call(this, attributes, display, level, prime);
    };
    MmlMath.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlLayoutNode.defaults), { mathvariant: 'normal', mathsize: 'normal', mathcolor: '', mathbackground: 'transparent', dir: 'ltr', scriptlevel: 0, displaystyle: false, display: 'inline', maxwidth: '', overflow: 'linebreak', altimg: '', 'altimg-width': '', 'altimg-height': '', 'altimg-valign': '', alttext: '', cdgroup: '', scriptsizemultiplier: 1 / Math.sqrt(2), scriptminsize: '8px', infixlinebreakstyle: 'before', lineleading: '1ex', linebreakmultchar: '\u2062', indentshift: 'auto', indentalign: 'auto', indenttarget: '', indentalignfirst: 'indentalign', indentshiftfirst: 'indentshift', indentalignlast: 'indentalign', indentshiftlast: 'indentshift' });
    return MmlMath;
}(MmlNode_js_1.AbstractMmlLayoutNode));
exports.MmlMath = MmlMath;
//# sourceMappingURL=math.js.map

/***/ }),

/***/ 99031:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMenclose = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMenclose = (function (_super) {
    __extends(MmlMenclose, _super);
    function MmlMenclose() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.texclass = MmlNode_js_1.TEXCLASS.ORD;
        return _this;
    }
    Object.defineProperty(MmlMenclose.prototype, "kind", {
        get: function () {
            return 'menclose';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMenclose.prototype, "arity", {
        get: function () {
            return -1;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMenclose.prototype, "linebreakContininer", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    MmlMenclose.prototype.setTeXclass = function (prev) {
        prev = this.childNodes[0].setTeXclass(prev);
        this.updateTeXclass(this.childNodes[0]);
        return prev;
    };
    MmlMenclose.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlNode.defaults), { notation: 'longdiv' });
    return MmlMenclose;
}(MmlNode_js_1.AbstractMmlNode));
exports.MmlMenclose = MmlMenclose;
//# sourceMappingURL=menclose.js.map

/***/ }),

/***/ 47149:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
var __values = (this && this.__values) || function(o) {
    var s = typeof Symbol === "function" && Symbol.iterator, m = s && o[s], i = 0;
    if (m) return m.call(o);
    if (o && typeof o.length === "number") return {
        next: function () {
            if (o && i >= o.length) o = void 0;
            return { value: o && o[i++], done: !o };
        }
    };
    throw new TypeError(s ? "Object is not iterable." : "Symbol.iterator is not defined.");
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMfenced = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMfenced = (function (_super) {
    __extends(MmlMfenced, _super);
    function MmlMfenced() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.texclass = MmlNode_js_1.TEXCLASS.INNER;
        _this.separators = [];
        _this.open = null;
        _this.close = null;
        return _this;
    }
    Object.defineProperty(MmlMfenced.prototype, "kind", {
        get: function () {
            return 'mfenced';
        },
        enumerable: false,
        configurable: true
    });
    MmlMfenced.prototype.setTeXclass = function (prev) {
        this.getPrevClass(prev);
        if (this.open) {
            prev = this.open.setTeXclass(prev);
        }
        if (this.childNodes[0]) {
            prev = this.childNodes[0].setTeXclass(prev);
        }
        for (var i = 1, m = this.childNodes.length; i < m; i++) {
            if (this.separators[i - 1]) {
                prev = this.separators[i - 1].setTeXclass(prev);
            }
            if (this.childNodes[i]) {
                prev = this.childNodes[i].setTeXclass(prev);
            }
        }
        if (this.close) {
            prev = this.close.setTeXclass(prev);
        }
        this.updateTeXclass(this.open);
        return prev;
    };
    MmlMfenced.prototype.setChildInheritedAttributes = function (attributes, display, level, prime) {
        var e_1, _a;
        this.addFakeNodes();
        try {
            for (var _b = __values([this.open, this.close].concat(this.separators)), _c = _b.next(); !_c.done; _c = _b.next()) {
                var child = _c.value;
                if (child) {
                    child.setInheritedAttributes(attributes, display, level, prime);
                }
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_1) throw e_1.error; }
        }
        _super.prototype.setChildInheritedAttributes.call(this, attributes, display, level, prime);
    };
    MmlMfenced.prototype.addFakeNodes = function () {
        var e_2, _a;
        var _b = this.attributes.getList('open', 'close', 'separators'), open = _b.open, close = _b.close, separators = _b.separators;
        open = open.replace(/[ \t\n\r]/g, '');
        close = close.replace(/[ \t\n\r]/g, '');
        separators = separators.replace(/[ \t\n\r]/g, '');
        if (open) {
            this.open = this.fakeNode(open, { fence: true, form: 'prefix' }, MmlNode_js_1.TEXCLASS.OPEN);
        }
        if (separators) {
            while (separators.length < this.childNodes.length - 1) {
                separators += separators.charAt(separators.length - 1);
            }
            var i = 0;
            try {
                for (var _c = __values(this.childNodes.slice(1)), _d = _c.next(); !_d.done; _d = _c.next()) {
                    var child = _d.value;
                    if (child) {
                        this.separators.push(this.fakeNode(separators.charAt(i++)));
                    }
                }
            }
            catch (e_2_1) { e_2 = { error: e_2_1 }; }
            finally {
                try {
                    if (_d && !_d.done && (_a = _c.return)) _a.call(_c);
                }
                finally { if (e_2) throw e_2.error; }
            }
        }
        if (close) {
            this.close = this.fakeNode(close, { fence: true, form: 'postfix' }, MmlNode_js_1.TEXCLASS.CLOSE);
        }
    };
    MmlMfenced.prototype.fakeNode = function (c, properties, texClass) {
        if (properties === void 0) { properties = {}; }
        if (texClass === void 0) { texClass = null; }
        var text = this.factory.create('text').setText(c);
        var node = this.factory.create('mo', properties, [text]);
        node.texClass = texClass;
        node.parent = this;
        return node;
    };
    MmlMfenced.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlNode.defaults), { open: '(', close: ')', separators: ',' });
    return MmlMfenced;
}(MmlNode_js_1.AbstractMmlNode));
exports.MmlMfenced = MmlMfenced;
//# sourceMappingURL=mfenced.js.map

/***/ }),

/***/ 76198:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
var __values = (this && this.__values) || function(o) {
    var s = typeof Symbol === "function" && Symbol.iterator, m = s && o[s], i = 0;
    if (m) return m.call(o);
    if (o && typeof o.length === "number") return {
        next: function () {
            if (o && i >= o.length) o = void 0;
            return { value: o && o[i++], done: !o };
        }
    };
    throw new TypeError(s ? "Object is not iterable." : "Symbol.iterator is not defined.");
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMfrac = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMfrac = (function (_super) {
    __extends(MmlMfrac, _super);
    function MmlMfrac() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMfrac.prototype, "kind", {
        get: function () {
            return 'mfrac';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMfrac.prototype, "arity", {
        get: function () {
            return 2;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMfrac.prototype, "linebreakContainer", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    MmlMfrac.prototype.setTeXclass = function (prev) {
        var e_1, _a;
        this.getPrevClass(prev);
        try {
            for (var _b = __values(this.childNodes), _c = _b.next(); !_c.done; _c = _b.next()) {
                var child = _c.value;
                child.setTeXclass(null);
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_1) throw e_1.error; }
        }
        return this;
    };
    MmlMfrac.prototype.setChildInheritedAttributes = function (attributes, display, level, prime) {
        if (!display || level > 0) {
            level++;
        }
        this.childNodes[0].setInheritedAttributes(attributes, false, level, prime);
        this.childNodes[1].setInheritedAttributes(attributes, false, level, true);
    };
    MmlMfrac.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlBaseNode.defaults), { linethickness: 'medium', numalign: 'center', denomalign: 'center', bevelled: false });
    return MmlMfrac;
}(MmlNode_js_1.AbstractMmlBaseNode));
exports.MmlMfrac = MmlMfrac;
//# sourceMappingURL=mfrac.js.map

/***/ }),

/***/ 49194:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMglyph = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMglyph = (function (_super) {
    __extends(MmlMglyph, _super);
    function MmlMglyph() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.texclass = MmlNode_js_1.TEXCLASS.ORD;
        return _this;
    }
    Object.defineProperty(MmlMglyph.prototype, "kind", {
        get: function () {
            return 'mglyph';
        },
        enumerable: false,
        configurable: true
    });
    MmlMglyph.prototype.verifyAttributes = function (options) {
        var _a = this.attributes.getList('src', 'fontfamily', 'index'), src = _a.src, fontfamily = _a.fontfamily, index = _a.index;
        if (src === '' && (fontfamily === '' || index === '')) {
            this.mError('mglyph must have either src or fontfamily and index attributes', options, true);
        }
        else {
            _super.prototype.verifyAttributes.call(this, options);
        }
    };
    MmlMglyph.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlTokenNode.defaults), { alt: '', src: '', index: '', width: 'auto', height: 'auto', valign: '0em' });
    return MmlMglyph;
}(MmlNode_js_1.AbstractMmlTokenNode));
exports.MmlMglyph = MmlMglyph;
//# sourceMappingURL=mglyph.js.map

/***/ }),

/***/ 91324:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMi = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMi = (function (_super) {
    __extends(MmlMi, _super);
    function MmlMi() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.texclass = MmlNode_js_1.TEXCLASS.ORD;
        return _this;
    }
    Object.defineProperty(MmlMi.prototype, "kind", {
        get: function () {
            return 'mi';
        },
        enumerable: false,
        configurable: true
    });
    MmlMi.prototype.setInheritedAttributes = function (attributes, display, level, prime) {
        if (attributes === void 0) { attributes = {}; }
        if (display === void 0) { display = false; }
        if (level === void 0) { level = 0; }
        if (prime === void 0) { prime = false; }
        _super.prototype.setInheritedAttributes.call(this, attributes, display, level, prime);
        var text = this.getText();
        if (text.match(MmlMi.singleCharacter) && !attributes.mathvariant) {
            this.attributes.setInherited('mathvariant', 'italic');
        }
    };
    MmlMi.prototype.setTeXclass = function (prev) {
        this.getPrevClass(prev);
        var name = this.getText();
        if (name.length > 1 && name.match(MmlMi.operatorName) &&
            this.attributes.get('mathvariant') === 'normal' &&
            this.getProperty('autoOP') === undefined &&
            this.getProperty('texClass') === undefined) {
            this.texClass = MmlNode_js_1.TEXCLASS.OP;
            this.setProperty('autoOP', true);
        }
        return this;
    };
    MmlMi.defaults = __assign({}, MmlNode_js_1.AbstractMmlTokenNode.defaults);
    MmlMi.operatorName = /^[a-z][a-z0-9]*$/i;
    MmlMi.singleCharacter = /^[\uD800-\uDBFF]?.[\u0300-\u036F\u1AB0-\u1ABE\u1DC0-\u1DFF\u20D0-\u20EF]*$/;
    return MmlMi;
}(MmlNode_js_1.AbstractMmlTokenNode));
exports.MmlMi = MmlMi;
//# sourceMappingURL=mi.js.map

/***/ }),

/***/ 80489:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlNone = exports.MmlMprescripts = exports.MmlMmultiscripts = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var msubsup_js_1 = __webpack_require__(41376);
var MmlMmultiscripts = (function (_super) {
    __extends(MmlMmultiscripts, _super);
    function MmlMmultiscripts() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMmultiscripts.prototype, "kind", {
        get: function () {
            return 'mmultiscripts';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMmultiscripts.prototype, "arity", {
        get: function () {
            return 1;
        },
        enumerable: false,
        configurable: true
    });
    MmlMmultiscripts.prototype.setChildInheritedAttributes = function (attributes, display, level, prime) {
        this.childNodes[0].setInheritedAttributes(attributes, display, level, prime);
        var prescripts = false;
        for (var i = 1, n = 0; i < this.childNodes.length; i++) {
            var child = this.childNodes[i];
            if (child.isKind('mprescripts')) {
                if (!prescripts) {
                    prescripts = true;
                    if (i % 2 === 0) {
                        var mrow = this.factory.create('mrow');
                        this.childNodes.splice(i, 0, mrow);
                        mrow.parent = this;
                        i++;
                    }
                }
            }
            else {
                var primestyle = prime || (n % 2 === 0);
                child.setInheritedAttributes(attributes, false, level + 1, primestyle);
                n++;
            }
        }
        if (this.childNodes.length % 2 === (prescripts ? 1 : 0)) {
            this.appendChild(this.factory.create('mrow'));
            this.childNodes[this.childNodes.length - 1].setInheritedAttributes(attributes, false, level + 1, prime);
        }
    };
    MmlMmultiscripts.prototype.verifyChildren = function (options) {
        var prescripts = false;
        var fix = options['fixMmultiscripts'];
        for (var i = 0; i < this.childNodes.length; i++) {
            var child = this.childNodes[i];
            if (child.isKind('mprescripts')) {
                if (prescripts) {
                    child.mError(child.kind + ' can only appear once in ' + this.kind, options, true);
                }
                else {
                    prescripts = true;
                    if (i % 2 === 0 && !fix) {
                        this.mError('There must be an equal number of prescripts of each type', options);
                    }
                }
            }
        }
        if (this.childNodes.length % 2 === (prescripts ? 1 : 0) && !fix) {
            this.mError('There must be an equal number of scripts of each type', options);
        }
        _super.prototype.verifyChildren.call(this, options);
    };
    MmlMmultiscripts.defaults = __assign({}, msubsup_js_1.MmlMsubsup.defaults);
    return MmlMmultiscripts;
}(msubsup_js_1.MmlMsubsup));
exports.MmlMmultiscripts = MmlMmultiscripts;
var MmlMprescripts = (function (_super) {
    __extends(MmlMprescripts, _super);
    function MmlMprescripts() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMprescripts.prototype, "kind", {
        get: function () {
            return 'mprescripts';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMprescripts.prototype, "arity", {
        get: function () {
            return 0;
        },
        enumerable: false,
        configurable: true
    });
    MmlMprescripts.prototype.verifyTree = function (options) {
        _super.prototype.verifyTree.call(this, options);
        if (this.parent && !this.parent.isKind('mmultiscripts')) {
            this.mError(this.kind + ' must be a child of mmultiscripts', options, true);
        }
    };
    MmlMprescripts.defaults = __assign({}, MmlNode_js_1.AbstractMmlNode.defaults);
    return MmlMprescripts;
}(MmlNode_js_1.AbstractMmlNode));
exports.MmlMprescripts = MmlMprescripts;
var MmlNone = (function (_super) {
    __extends(MmlNone, _super);
    function MmlNone() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlNone.prototype, "kind", {
        get: function () {
            return 'none';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlNone.prototype, "arity", {
        get: function () {
            return 0;
        },
        enumerable: false,
        configurable: true
    });
    MmlNone.prototype.verifyTree = function (options) {
        _super.prototype.verifyTree.call(this, options);
        if (this.parent && !this.parent.isKind('mmultiscripts')) {
            this.mError(this.kind + ' must be a child of mmultiscripts', options, true);
        }
    };
    MmlNone.defaults = __assign({}, MmlNode_js_1.AbstractMmlNode.defaults);
    return MmlNone;
}(MmlNode_js_1.AbstractMmlNode));
exports.MmlNone = MmlNone;
//# sourceMappingURL=mmultiscripts.js.map

/***/ }),

/***/ 14734:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMn = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMn = (function (_super) {
    __extends(MmlMn, _super);
    function MmlMn() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.texclass = MmlNode_js_1.TEXCLASS.ORD;
        return _this;
    }
    Object.defineProperty(MmlMn.prototype, "kind", {
        get: function () {
            return 'mn';
        },
        enumerable: false,
        configurable: true
    });
    MmlMn.defaults = __assign({}, MmlNode_js_1.AbstractMmlTokenNode.defaults);
    return MmlMn;
}(MmlNode_js_1.AbstractMmlTokenNode));
exports.MmlMn = MmlMn;
//# sourceMappingURL=mn.js.map

/***/ }),

/***/ 70596:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMpadded = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMpadded = (function (_super) {
    __extends(MmlMpadded, _super);
    function MmlMpadded() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMpadded.prototype, "kind", {
        get: function () {
            return 'mpadded';
        },
        enumerable: false,
        configurable: true
    });
    MmlMpadded.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlLayoutNode.defaults), { width: '', height: '', depth: '', lspace: 0, voffset: 0 });
    return MmlMpadded;
}(MmlNode_js_1.AbstractMmlLayoutNode));
exports.MmlMpadded = MmlMpadded;
//# sourceMappingURL=mpadded.js.map

/***/ }),

/***/ 79020:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMroot = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMroot = (function (_super) {
    __extends(MmlMroot, _super);
    function MmlMroot() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.texclass = MmlNode_js_1.TEXCLASS.ORD;
        return _this;
    }
    Object.defineProperty(MmlMroot.prototype, "kind", {
        get: function () {
            return 'mroot';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMroot.prototype, "arity", {
        get: function () {
            return 2;
        },
        enumerable: false,
        configurable: true
    });
    MmlMroot.prototype.setTeXclass = function (prev) {
        this.getPrevClass(prev);
        this.childNodes[0].setTeXclass(null);
        this.childNodes[1].setTeXclass(null);
        return this;
    };
    MmlMroot.prototype.setChildInheritedAttributes = function (attributes, display, level, prime) {
        this.childNodes[0].setInheritedAttributes(attributes, display, level, true);
        this.childNodes[1].setInheritedAttributes(attributes, false, level + 2, prime);
    };
    MmlMroot.defaults = __assign({}, MmlNode_js_1.AbstractMmlNode.defaults);
    return MmlMroot;
}(MmlNode_js_1.AbstractMmlNode));
exports.MmlMroot = MmlMroot;
//# sourceMappingURL=mroot.js.map

/***/ }),

/***/ 70938:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
var __values = (this && this.__values) || function(o) {
    var s = typeof Symbol === "function" && Symbol.iterator, m = s && o[s], i = 0;
    if (m) return m.call(o);
    if (o && typeof o.length === "number") return {
        next: function () {
            if (o && i >= o.length) o = void 0;
            return { value: o && o[i++], done: !o };
        }
    };
    throw new TypeError(s ? "Object is not iterable." : "Symbol.iterator is not defined.");
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlInferredMrow = exports.MmlMrow = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMrow = (function (_super) {
    __extends(MmlMrow, _super);
    function MmlMrow() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this._core = null;
        return _this;
    }
    Object.defineProperty(MmlMrow.prototype, "kind", {
        get: function () {
            return 'mrow';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMrow.prototype, "isSpacelike", {
        get: function () {
            var e_1, _a;
            try {
                for (var _b = __values(this.childNodes), _c = _b.next(); !_c.done; _c = _b.next()) {
                    var child = _c.value;
                    if (!child.isSpacelike) {
                        return false;
                    }
                }
            }
            catch (e_1_1) { e_1 = { error: e_1_1 }; }
            finally {
                try {
                    if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
                }
                finally { if (e_1) throw e_1.error; }
            }
            return true;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMrow.prototype, "isEmbellished", {
        get: function () {
            var e_2, _a;
            var embellished = false;
            var i = 0;
            try {
                for (var _b = __values(this.childNodes), _c = _b.next(); !_c.done; _c = _b.next()) {
                    var child = _c.value;
                    if (child) {
                        if (child.isEmbellished) {
                            if (embellished) {
                                return false;
                            }
                            embellished = true;
                            this._core = i;
                        }
                        else if (!child.isSpacelike) {
                            return false;
                        }
                    }
                    i++;
                }
            }
            catch (e_2_1) { e_2 = { error: e_2_1 }; }
            finally {
                try {
                    if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
                }
                finally { if (e_2) throw e_2.error; }
            }
            return embellished;
        },
        enumerable: false,
        configurable: true
    });
    MmlMrow.prototype.core = function () {
        if (!this.isEmbellished || this._core == null) {
            return this;
        }
        return this.childNodes[this._core];
    };
    MmlMrow.prototype.coreMO = function () {
        if (!this.isEmbellished || this._core == null) {
            return this;
        }
        return this.childNodes[this._core].coreMO();
    };
    MmlMrow.prototype.nonSpaceLength = function () {
        var e_3, _a;
        var n = 0;
        try {
            for (var _b = __values(this.childNodes), _c = _b.next(); !_c.done; _c = _b.next()) {
                var child = _c.value;
                if (child && !child.isSpacelike) {
                    n++;
                }
            }
        }
        catch (e_3_1) { e_3 = { error: e_3_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_3) throw e_3.error; }
        }
        return n;
    };
    MmlMrow.prototype.firstNonSpace = function () {
        var e_4, _a;
        try {
            for (var _b = __values(this.childNodes), _c = _b.next(); !_c.done; _c = _b.next()) {
                var child = _c.value;
                if (child && !child.isSpacelike) {
                    return child;
                }
            }
        }
        catch (e_4_1) { e_4 = { error: e_4_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_4) throw e_4.error; }
        }
        return null;
    };
    MmlMrow.prototype.lastNonSpace = function () {
        var i = this.childNodes.length;
        while (--i >= 0) {
            var child = this.childNodes[i];
            if (child && !child.isSpacelike) {
                return child;
            }
        }
        return null;
    };
    MmlMrow.prototype.setTeXclass = function (prev) {
        var e_5, _a, e_6, _b;
        if (this.getProperty('open') != null || this.getProperty('close') != null) {
            this.getPrevClass(prev);
            prev = null;
            try {
                for (var _c = __values(this.childNodes), _d = _c.next(); !_d.done; _d = _c.next()) {
                    var child = _d.value;
                    prev = child.setTeXclass(prev);
                }
            }
            catch (e_5_1) { e_5 = { error: e_5_1 }; }
            finally {
                try {
                    if (_d && !_d.done && (_a = _c.return)) _a.call(_c);
                }
                finally { if (e_5) throw e_5.error; }
            }
            if (this.texClass == null) {
                this.texClass = MmlNode_js_1.TEXCLASS.INNER;
            }
        }
        else {
            try {
                for (var _e = __values(this.childNodes), _f = _e.next(); !_f.done; _f = _e.next()) {
                    var child = _f.value;
                    prev = child.setTeXclass(prev);
                }
            }
            catch (e_6_1) { e_6 = { error: e_6_1 }; }
            finally {
                try {
                    if (_f && !_f.done && (_b = _e.return)) _b.call(_e);
                }
                finally { if (e_6) throw e_6.error; }
            }
            if (this.childNodes[0]) {
                this.updateTeXclass(this.childNodes[0]);
            }
        }
        return prev;
    };
    MmlMrow.defaults = __assign({}, MmlNode_js_1.AbstractMmlNode.defaults);
    return MmlMrow;
}(MmlNode_js_1.AbstractMmlNode));
exports.MmlMrow = MmlMrow;
var MmlInferredMrow = (function (_super) {
    __extends(MmlInferredMrow, _super);
    function MmlInferredMrow() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlInferredMrow.prototype, "kind", {
        get: function () {
            return 'inferredMrow';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlInferredMrow.prototype, "isInferred", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlInferredMrow.prototype, "notParent", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    MmlInferredMrow.prototype.toString = function () {
        return '[' + this.childNodes.join(',') + ']';
    };
    MmlInferredMrow.defaults = MmlMrow.defaults;
    return MmlInferredMrow;
}(MmlMrow));
exports.MmlInferredMrow = MmlInferredMrow;
//# sourceMappingURL=mrow.js.map

/***/ }),

/***/ 75375:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMs = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMs = (function (_super) {
    __extends(MmlMs, _super);
    function MmlMs() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.texclass = MmlNode_js_1.TEXCLASS.ORD;
        return _this;
    }
    Object.defineProperty(MmlMs.prototype, "kind", {
        get: function () {
            return 'ms';
        },
        enumerable: false,
        configurable: true
    });
    MmlMs.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlTokenNode.defaults), { lquote: '"', rquote: '"' });
    return MmlMs;
}(MmlNode_js_1.AbstractMmlTokenNode));
exports.MmlMs = MmlMs;
//# sourceMappingURL=ms.js.map

/***/ }),

/***/ 58321:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMspace = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMspace = (function (_super) {
    __extends(MmlMspace, _super);
    function MmlMspace() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.texclass = MmlNode_js_1.TEXCLASS.NONE;
        return _this;
    }
    MmlMspace.prototype.setTeXclass = function (prev) {
        return prev;
    };
    Object.defineProperty(MmlMspace.prototype, "kind", {
        get: function () {
            return 'mspace';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMspace.prototype, "arity", {
        get: function () {
            return 0;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMspace.prototype, "isSpacelike", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMspace.prototype, "hasNewline", {
        get: function () {
            var attributes = this.attributes;
            return (attributes.getExplicit('width') == null && attributes.getExplicit('height') == null &&
                attributes.getExplicit('depth') == null && attributes.get('linebreak') === 'newline');
        },
        enumerable: false,
        configurable: true
    });
    MmlMspace.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlTokenNode.defaults), { width: '0em', height: '0ex', depth: '0ex', linebreak: 'auto' });
    return MmlMspace;
}(MmlNode_js_1.AbstractMmlTokenNode));
exports.MmlMspace = MmlMspace;
//# sourceMappingURL=mspace.js.map

/***/ }),

/***/ 42061:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMsqrt = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMsqrt = (function (_super) {
    __extends(MmlMsqrt, _super);
    function MmlMsqrt() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.texclass = MmlNode_js_1.TEXCLASS.ORD;
        return _this;
    }
    Object.defineProperty(MmlMsqrt.prototype, "kind", {
        get: function () {
            return 'msqrt';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMsqrt.prototype, "arity", {
        get: function () {
            return -1;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMsqrt.prototype, "linebreakContainer", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    MmlMsqrt.prototype.setTeXclass = function (prev) {
        this.getPrevClass(prev);
        this.childNodes[0].setTeXclass(null);
        return this;
    };
    MmlMsqrt.prototype.setChildInheritedAttributes = function (attributes, display, level, _prime) {
        this.childNodes[0].setInheritedAttributes(attributes, display, level, true);
    };
    MmlMsqrt.defaults = __assign({}, MmlNode_js_1.AbstractMmlNode.defaults);
    return MmlMsqrt;
}(MmlNode_js_1.AbstractMmlNode));
exports.MmlMsqrt = MmlMsqrt;
//# sourceMappingURL=msqrt.js.map

/***/ }),

/***/ 41376:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMsup = exports.MmlMsub = exports.MmlMsubsup = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMsubsup = (function (_super) {
    __extends(MmlMsubsup, _super);
    function MmlMsubsup() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMsubsup.prototype, "kind", {
        get: function () {
            return 'msubsup';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMsubsup.prototype, "arity", {
        get: function () {
            return 3;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMsubsup.prototype, "base", {
        get: function () {
            return 0;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMsubsup.prototype, "sub", {
        get: function () {
            return 1;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMsubsup.prototype, "sup", {
        get: function () {
            return 2;
        },
        enumerable: false,
        configurable: true
    });
    MmlMsubsup.prototype.setChildInheritedAttributes = function (attributes, display, level, prime) {
        var nodes = this.childNodes;
        nodes[0].setInheritedAttributes(attributes, display, level, prime);
        nodes[1].setInheritedAttributes(attributes, false, level + 1, prime || this.sub === 1);
        if (!nodes[2]) {
            return;
        }
        nodes[2].setInheritedAttributes(attributes, false, level + 1, prime || this.sub === 2);
    };
    MmlMsubsup.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlBaseNode.defaults), { subscriptshift: '', superscriptshift: '' });
    return MmlMsubsup;
}(MmlNode_js_1.AbstractMmlBaseNode));
exports.MmlMsubsup = MmlMsubsup;
var MmlMsub = (function (_super) {
    __extends(MmlMsub, _super);
    function MmlMsub() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMsub.prototype, "kind", {
        get: function () {
            return 'msub';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMsub.prototype, "arity", {
        get: function () {
            return 2;
        },
        enumerable: false,
        configurable: true
    });
    MmlMsub.defaults = __assign({}, MmlMsubsup.defaults);
    return MmlMsub;
}(MmlMsubsup));
exports.MmlMsub = MmlMsub;
var MmlMsup = (function (_super) {
    __extends(MmlMsup, _super);
    function MmlMsup() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMsup.prototype, "kind", {
        get: function () {
            return 'msup';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMsup.prototype, "arity", {
        get: function () {
            return 2;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMsup.prototype, "sup", {
        get: function () {
            return 1;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMsup.prototype, "sub", {
        get: function () {
            return 2;
        },
        enumerable: false,
        configurable: true
    });
    MmlMsup.defaults = __assign({}, MmlMsubsup.defaults);
    return MmlMsup;
}(MmlMsubsup));
exports.MmlMsup = MmlMsup;
//# sourceMappingURL=msubsup.js.map

/***/ }),

/***/ 60324:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
var __values = (this && this.__values) || function(o) {
    var s = typeof Symbol === "function" && Symbol.iterator, m = s && o[s], i = 0;
    if (m) return m.call(o);
    if (o && typeof o.length === "number") return {
        next: function () {
            if (o && i >= o.length) o = void 0;
            return { value: o && o[i++], done: !o };
        }
    };
    throw new TypeError(s ? "Object is not iterable." : "Symbol.iterator is not defined.");
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMtable = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var string_js_1 = __webpack_require__(55089);
var MmlMtable = (function (_super) {
    __extends(MmlMtable, _super);
    function MmlMtable() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.properties = {
            useHeight: true
        };
        _this.texclass = MmlNode_js_1.TEXCLASS.ORD;
        return _this;
    }
    Object.defineProperty(MmlMtable.prototype, "kind", {
        get: function () {
            return 'mtable';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMtable.prototype, "linebreakContainer", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    MmlMtable.prototype.setInheritedAttributes = function (attributes, display, level, prime) {
        var e_1, _a;
        try {
            for (var indentAttributes_1 = __values(MmlNode_js_1.indentAttributes), indentAttributes_1_1 = indentAttributes_1.next(); !indentAttributes_1_1.done; indentAttributes_1_1 = indentAttributes_1.next()) {
                var name_1 = indentAttributes_1_1.value;
                if (attributes[name_1]) {
                    this.attributes.setInherited(name_1, attributes[name_1][1]);
                }
                if (this.attributes.getExplicit(name_1) !== undefined) {
                    delete (this.attributes.getAllAttributes())[name_1];
                }
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (indentAttributes_1_1 && !indentAttributes_1_1.done && (_a = indentAttributes_1.return)) _a.call(indentAttributes_1);
            }
            finally { if (e_1) throw e_1.error; }
        }
        _super.prototype.setInheritedAttributes.call(this, attributes, display, level, prime);
    };
    MmlMtable.prototype.setChildInheritedAttributes = function (attributes, display, level, _prime) {
        var e_2, _a, e_3, _b;
        try {
            for (var _c = __values(this.childNodes), _d = _c.next(); !_d.done; _d = _c.next()) {
                var child = _d.value;
                if (!child.isKind('mtr')) {
                    this.replaceChild(this.factory.create('mtr'), child)
                        .appendChild(child);
                }
            }
        }
        catch (e_2_1) { e_2 = { error: e_2_1 }; }
        finally {
            try {
                if (_d && !_d.done && (_a = _c.return)) _a.call(_c);
            }
            finally { if (e_2) throw e_2.error; }
        }
        level = this.getProperty('scriptlevel') || level;
        display = !!(this.attributes.getExplicit('displaystyle') || this.attributes.getDefault('displaystyle'));
        attributes = this.addInheritedAttributes(attributes, {
            columnalign: this.attributes.get('columnalign'),
            rowalign: 'center'
        });
        var cramped = this.attributes.getExplicit('data-cramped');
        var ralign = (0, string_js_1.split)(this.attributes.get('rowalign'));
        try {
            for (var _e = __values(this.childNodes), _f = _e.next(); !_f.done; _f = _e.next()) {
                var child = _f.value;
                attributes.rowalign[1] = ralign.shift() || attributes.rowalign[1];
                child.setInheritedAttributes(attributes, display, level, !!cramped);
            }
        }
        catch (e_3_1) { e_3 = { error: e_3_1 }; }
        finally {
            try {
                if (_f && !_f.done && (_b = _e.return)) _b.call(_e);
            }
            finally { if (e_3) throw e_3.error; }
        }
    };
    MmlMtable.prototype.verifyChildren = function (options) {
        var mtr = null;
        var factory = this.factory;
        for (var i = 0; i < this.childNodes.length; i++) {
            var child = this.childNodes[i];
            if (child.isKind('mtr')) {
                mtr = null;
            }
            else {
                var isMtd = child.isKind('mtd');
                if (mtr) {
                    this.removeChild(child);
                    i--;
                }
                else {
                    mtr = this.replaceChild(factory.create('mtr'), child);
                }
                mtr.appendChild(isMtd ? child : factory.create('mtd', {}, [child]));
                if (!options['fixMtables']) {
                    child.parent.removeChild(child);
                    child.parent = this;
                    isMtd && mtr.appendChild(factory.create('mtd'));
                    var merror = child.mError('Children of ' + this.kind + ' must be mtr or mlabeledtr', options, isMtd);
                    mtr.childNodes[mtr.childNodes.length - 1].appendChild(merror);
                }
            }
        }
        _super.prototype.verifyChildren.call(this, options);
    };
    MmlMtable.prototype.setTeXclass = function (prev) {
        var e_4, _a;
        this.getPrevClass(prev);
        try {
            for (var _b = __values(this.childNodes), _c = _b.next(); !_c.done; _c = _b.next()) {
                var child = _c.value;
                child.setTeXclass(null);
            }
        }
        catch (e_4_1) { e_4 = { error: e_4_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_4) throw e_4.error; }
        }
        return this;
    };
    MmlMtable.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlNode.defaults), { align: 'axis', rowalign: 'baseline', columnalign: 'center', groupalign: '{left}', alignmentscope: true, columnwidth: 'auto', width: 'auto', rowspacing: '1ex', columnspacing: '.8em', rowlines: 'none', columnlines: 'none', frame: 'none', framespacing: '0.4em 0.5ex', equalrows: false, equalcolumns: false, displaystyle: false, side: 'right', minlabelspacing: '0.8em' });
    return MmlMtable;
}(MmlNode_js_1.AbstractMmlNode));
exports.MmlMtable = MmlMtable;
//# sourceMappingURL=mtable.js.map

/***/ }),

/***/ 24955:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMtd = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var Attributes_js_1 = __webpack_require__(95518);
var MmlMtd = (function (_super) {
    __extends(MmlMtd, _super);
    function MmlMtd() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMtd.prototype, "kind", {
        get: function () {
            return 'mtd';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMtd.prototype, "arity", {
        get: function () {
            return -1;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMtd.prototype, "linebreakContainer", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    MmlMtd.prototype.verifyChildren = function (options) {
        if (this.parent && !this.parent.isKind('mtr')) {
            this.mError(this.kind + ' can only be a child of an mtr or mlabeledtr', options, true);
            return;
        }
        _super.prototype.verifyChildren.call(this, options);
    };
    MmlMtd.prototype.setTeXclass = function (prev) {
        this.getPrevClass(prev);
        this.childNodes[0].setTeXclass(null);
        return this;
    };
    MmlMtd.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlBaseNode.defaults), { rowspan: 1, columnspan: 1, rowalign: Attributes_js_1.INHERIT, columnalign: Attributes_js_1.INHERIT, groupalign: Attributes_js_1.INHERIT });
    return MmlMtd;
}(MmlNode_js_1.AbstractMmlBaseNode));
exports.MmlMtd = MmlMtd;
//# sourceMappingURL=mtd.js.map

/***/ }),

/***/ 64957:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMtext = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMtext = (function (_super) {
    __extends(MmlMtext, _super);
    function MmlMtext() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.texclass = MmlNode_js_1.TEXCLASS.ORD;
        return _this;
    }
    Object.defineProperty(MmlMtext.prototype, "kind", {
        get: function () {
            return 'mtext';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMtext.prototype, "isSpacelike", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    MmlMtext.defaults = __assign({}, MmlNode_js_1.AbstractMmlTokenNode.defaults);
    return MmlMtext;
}(MmlNode_js_1.AbstractMmlTokenNode));
exports.MmlMtext = MmlMtext;
//# sourceMappingURL=mtext.js.map

/***/ }),

/***/ 84760:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
var __values = (this && this.__values) || function(o) {
    var s = typeof Symbol === "function" && Symbol.iterator, m = s && o[s], i = 0;
    if (m) return m.call(o);
    if (o && typeof o.length === "number") return {
        next: function () {
            if (o && i >= o.length) o = void 0;
            return { value: o && o[i++], done: !o };
        }
    };
    throw new TypeError(s ? "Object is not iterable." : "Symbol.iterator is not defined.");
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMlabeledtr = exports.MmlMtr = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var Attributes_js_1 = __webpack_require__(95518);
var string_js_1 = __webpack_require__(55089);
var MmlMtr = (function (_super) {
    __extends(MmlMtr, _super);
    function MmlMtr() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMtr.prototype, "kind", {
        get: function () {
            return 'mtr';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMtr.prototype, "linebreakContainer", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    MmlMtr.prototype.setChildInheritedAttributes = function (attributes, display, level, prime) {
        var e_1, _a, e_2, _b;
        try {
            for (var _c = __values(this.childNodes), _d = _c.next(); !_d.done; _d = _c.next()) {
                var child = _d.value;
                if (!child.isKind('mtd')) {
                    this.replaceChild(this.factory.create('mtd'), child)
                        .appendChild(child);
                }
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (_d && !_d.done && (_a = _c.return)) _a.call(_c);
            }
            finally { if (e_1) throw e_1.error; }
        }
        var calign = (0, string_js_1.split)(this.attributes.get('columnalign'));
        if (this.arity === 1) {
            calign.unshift(this.parent.attributes.get('side'));
        }
        attributes = this.addInheritedAttributes(attributes, {
            rowalign: this.attributes.get('rowalign'),
            columnalign: 'center'
        });
        try {
            for (var _e = __values(this.childNodes), _f = _e.next(); !_f.done; _f = _e.next()) {
                var child = _f.value;
                attributes.columnalign[1] = calign.shift() || attributes.columnalign[1];
                child.setInheritedAttributes(attributes, display, level, prime);
            }
        }
        catch (e_2_1) { e_2 = { error: e_2_1 }; }
        finally {
            try {
                if (_f && !_f.done && (_b = _e.return)) _b.call(_e);
            }
            finally { if (e_2) throw e_2.error; }
        }
    };
    MmlMtr.prototype.verifyChildren = function (options) {
        var e_3, _a;
        if (this.parent && !this.parent.isKind('mtable')) {
            this.mError(this.kind + ' can only be a child of an mtable', options, true);
            return;
        }
        try {
            for (var _b = __values(this.childNodes), _c = _b.next(); !_c.done; _c = _b.next()) {
                var child = _c.value;
                if (!child.isKind('mtd')) {
                    var mtd = this.replaceChild(this.factory.create('mtd'), child);
                    mtd.appendChild(child);
                    if (!options['fixMtables']) {
                        child.mError('Children of ' + this.kind + ' must be mtd', options);
                    }
                }
            }
        }
        catch (e_3_1) { e_3 = { error: e_3_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_3) throw e_3.error; }
        }
        _super.prototype.verifyChildren.call(this, options);
    };
    MmlMtr.prototype.setTeXclass = function (prev) {
        var e_4, _a;
        this.getPrevClass(prev);
        try {
            for (var _b = __values(this.childNodes), _c = _b.next(); !_c.done; _c = _b.next()) {
                var child = _c.value;
                child.setTeXclass(null);
            }
        }
        catch (e_4_1) { e_4 = { error: e_4_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_4) throw e_4.error; }
        }
        return this;
    };
    MmlMtr.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlNode.defaults), { rowalign: Attributes_js_1.INHERIT, columnalign: Attributes_js_1.INHERIT, groupalign: Attributes_js_1.INHERIT });
    return MmlMtr;
}(MmlNode_js_1.AbstractMmlNode));
exports.MmlMtr = MmlMtr;
var MmlMlabeledtr = (function (_super) {
    __extends(MmlMlabeledtr, _super);
    function MmlMlabeledtr() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMlabeledtr.prototype, "kind", {
        get: function () {
            return 'mlabeledtr';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMlabeledtr.prototype, "arity", {
        get: function () {
            return 1;
        },
        enumerable: false,
        configurable: true
    });
    return MmlMlabeledtr;
}(MmlMtr));
exports.MmlMlabeledtr = MmlMlabeledtr;
//# sourceMappingURL=mtr.js.map

/***/ }),

/***/ 75579:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlMover = exports.MmlMunder = exports.MmlMunderover = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlMunderover = (function (_super) {
    __extends(MmlMunderover, _super);
    function MmlMunderover() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMunderover.prototype, "kind", {
        get: function () {
            return 'munderover';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMunderover.prototype, "arity", {
        get: function () {
            return 3;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMunderover.prototype, "base", {
        get: function () {
            return 0;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMunderover.prototype, "under", {
        get: function () {
            return 1;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMunderover.prototype, "over", {
        get: function () {
            return 2;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMunderover.prototype, "linebreakContainer", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    MmlMunderover.prototype.setChildInheritedAttributes = function (attributes, display, level, prime) {
        var nodes = this.childNodes;
        nodes[0].setInheritedAttributes(attributes, display, level, prime || !!nodes[this.over]);
        var force = !!(!display && nodes[0].coreMO().attributes.get('movablelimits'));
        var ACCENTS = this.constructor.ACCENTS;
        nodes[1].setInheritedAttributes(attributes, false, this.getScriptlevel(ACCENTS[1], force, level), prime || this.under === 1);
        this.setInheritedAccent(1, ACCENTS[1], display, level, prime, force);
        if (!nodes[2]) {
            return;
        }
        nodes[2].setInheritedAttributes(attributes, false, this.getScriptlevel(ACCENTS[2], force, level), prime || this.under === 2);
        this.setInheritedAccent(2, ACCENTS[2], display, level, prime, force);
    };
    MmlMunderover.prototype.getScriptlevel = function (accent, force, level) {
        if (force || !this.attributes.get(accent)) {
            level++;
        }
        return level;
    };
    MmlMunderover.prototype.setInheritedAccent = function (n, accent, display, level, prime, force) {
        var node = this.childNodes[n];
        if (this.attributes.getExplicit(accent) == null && node.isEmbellished) {
            var value = node.coreMO().attributes.get('accent');
            this.attributes.setInherited(accent, value);
            if (value !== this.attributes.getDefault(accent)) {
                node.setInheritedAttributes({}, display, this.getScriptlevel(accent, force, level), prime);
            }
        }
    };
    MmlMunderover.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlBaseNode.defaults), { accent: false, accentunder: false, align: 'center' });
    MmlMunderover.ACCENTS = ['', 'accentunder', 'accent'];
    return MmlMunderover;
}(MmlNode_js_1.AbstractMmlBaseNode));
exports.MmlMunderover = MmlMunderover;
var MmlMunder = (function (_super) {
    __extends(MmlMunder, _super);
    function MmlMunder() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMunder.prototype, "kind", {
        get: function () {
            return 'munder';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMunder.prototype, "arity", {
        get: function () {
            return 2;
        },
        enumerable: false,
        configurable: true
    });
    MmlMunder.defaults = __assign({}, MmlMunderover.defaults);
    return MmlMunder;
}(MmlMunderover));
exports.MmlMunder = MmlMunder;
var MmlMover = (function (_super) {
    __extends(MmlMover, _super);
    function MmlMover() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlMover.prototype, "kind", {
        get: function () {
            return 'mover';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMover.prototype, "arity", {
        get: function () {
            return 2;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMover.prototype, "over", {
        get: function () {
            return 1;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlMover.prototype, "under", {
        get: function () {
            return 2;
        },
        enumerable: false,
        configurable: true
    });
    MmlMover.defaults = __assign({}, MmlMunderover.defaults);
    MmlMover.ACCENTS = ['', 'accent', 'accentunder'];
    return MmlMover;
}(MmlMunderover));
exports.MmlMover = MmlMover;
//# sourceMappingURL=munderover.js.map

/***/ }),

/***/ 10246:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __extends = (this && this.__extends) || (function () {
    var extendStatics = function (d, b) {
        extendStatics = Object.setPrototypeOf ||
            ({ __proto__: [] } instanceof Array && function (d, b) { d.__proto__ = b; }) ||
            function (d, b) { for (var p in b) if (Object.prototype.hasOwnProperty.call(b, p)) d[p] = b[p]; };
        return extendStatics(d, b);
    };
    return function (d, b) {
        if (typeof b !== "function" && b !== null)
            throw new TypeError("Class extends value " + String(b) + " is not a constructor or null");
        extendStatics(d, b);
        function __() { this.constructor = d; }
        d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
    };
})();
var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlAnnotation = exports.MmlAnnotationXML = exports.MmlSemantics = void 0;
var MmlNode_js_1 = __webpack_require__(83045);
var MmlSemantics = (function (_super) {
    __extends(MmlSemantics, _super);
    function MmlSemantics() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlSemantics.prototype, "kind", {
        get: function () {
            return 'semantics';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlSemantics.prototype, "arity", {
        get: function () {
            return 1;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlSemantics.prototype, "notParent", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    MmlSemantics.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlBaseNode.defaults), { definitionUrl: null, encoding: null });
    return MmlSemantics;
}(MmlNode_js_1.AbstractMmlBaseNode));
exports.MmlSemantics = MmlSemantics;
var MmlAnnotationXML = (function (_super) {
    __extends(MmlAnnotationXML, _super);
    function MmlAnnotationXML() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlAnnotationXML.prototype, "kind", {
        get: function () {
            return 'annotation-xml';
        },
        enumerable: false,
        configurable: true
    });
    MmlAnnotationXML.prototype.setChildInheritedAttributes = function () { };
    MmlAnnotationXML.defaults = __assign(__assign({}, MmlNode_js_1.AbstractMmlNode.defaults), { definitionUrl: null, encoding: null, cd: 'mathmlkeys', name: '', src: null });
    return MmlAnnotationXML;
}(MmlNode_js_1.AbstractMmlNode));
exports.MmlAnnotationXML = MmlAnnotationXML;
var MmlAnnotation = (function (_super) {
    __extends(MmlAnnotation, _super);
    function MmlAnnotation() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.properties = {
            isChars: true
        };
        return _this;
    }
    Object.defineProperty(MmlAnnotation.prototype, "kind", {
        get: function () {
            return 'annotation';
        },
        enumerable: false,
        configurable: true
    });
    MmlAnnotation.defaults = __assign({}, MmlAnnotationXML.defaults);
    return MmlAnnotation;
}(MmlAnnotationXML));
exports.MmlAnnotation = MmlAnnotation;
//# sourceMappingURL=semantics.js.map

/***/ }),

/***/ 72901:
/***/ (function(__unused_webpack_module, exports) {


var __values = (this && this.__values) || function(o) {
    var s = typeof Symbol === "function" && Symbol.iterator, m = s && o[s], i = 0;
    if (m) return m.call(o);
    if (o && typeof o.length === "number") return {
        next: function () {
            if (o && i >= o.length) o = void 0;
            return { value: o && o[i++], done: !o };
        }
    };
    throw new TypeError(s ? "Object is not iterable." : "Symbol.iterator is not defined.");
};
var __read = (this && this.__read) || function (o, n) {
    var m = typeof Symbol === "function" && o[Symbol.iterator];
    if (!m) return o;
    var i = m.call(o), r, ar = [], e;
    try {
        while ((n === void 0 || n-- > 0) && !(r = i.next()).done) ar.push(r.value);
    }
    catch (error) { e = { error: error }; }
    finally {
        try {
            if (r && !r.done && (m = i["return"])) m.call(i);
        }
        finally { if (e) throw e.error; }
    }
    return ar;
};
var __spreadArray = (this && this.__spreadArray) || function (to, from, pack) {
    if (pack || arguments.length === 2) for (var i = 0, l = from.length, ar; i < l; i++) {
        if (ar || !(i in from)) {
            if (!ar) ar = Array.prototype.slice.call(from, 0, i);
            ar[i] = from[i];
        }
    }
    return to.concat(ar || Array.prototype.slice.call(from));
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.AbstractFactory = void 0;
var AbstractFactory = (function () {
    function AbstractFactory(nodes) {
        var e_1, _a;
        if (nodes === void 0) { nodes = null; }
        this.defaultKind = 'unknown';
        this.nodeMap = new Map();
        this.node = {};
        if (nodes === null) {
            nodes = this.constructor.defaultNodes;
        }
        try {
            for (var _b = __values(Object.keys(nodes)), _c = _b.next(); !_c.done; _c = _b.next()) {
                var kind = _c.value;
                this.setNodeClass(kind, nodes[kind]);
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_1) throw e_1.error; }
        }
    }
    AbstractFactory.prototype.create = function (kind) {
        var args = [];
        for (var _i = 1; _i < arguments.length; _i++) {
            args[_i - 1] = arguments[_i];
        }
        return (this.node[kind] || this.node[this.defaultKind]).apply(void 0, __spreadArray([], __read(args), false));
    };
    AbstractFactory.prototype.setNodeClass = function (kind, nodeClass) {
        this.nodeMap.set(kind, nodeClass);
        var THIS = this;
        var KIND = this.nodeMap.get(kind);
        this.node[kind] = function () {
            var args = [];
            for (var _i = 0; _i < arguments.length; _i++) {
                args[_i] = arguments[_i];
            }
            return new (KIND.bind.apply(KIND, __spreadArray([void 0, THIS], __read(args), false)))();
        };
    };
    AbstractFactory.prototype.getNodeClass = function (kind) {
        return this.nodeMap.get(kind);
    };
    AbstractFactory.prototype.deleteNodeClass = function (kind) {
        this.nodeMap.delete(kind);
        delete this.node[kind];
    };
    AbstractFactory.prototype.nodeIsKind = function (node, kind) {
        return (node instanceof this.getNodeClass(kind));
    };
    AbstractFactory.prototype.getKinds = function () {
        return Array.from(this.nodeMap.keys());
    };
    AbstractFactory.defaultNodes = {};
    return AbstractFactory;
}());
exports.AbstractFactory = AbstractFactory;
//# sourceMappingURL=Factory.js.map

/***/ })

}]);
//# sourceMappingURL=28.b5145a84e3a511427e72.js.map?v=b5145a84e3a511427e72
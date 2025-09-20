"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[4105],{

/***/ 74105:
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
exports.AssistiveMmlHandler = exports.AssistiveMmlMathDocumentMixin = exports.AssistiveMmlMathItemMixin = exports.LimitedMmlVisitor = void 0;
var MathItem_js_1 = __webpack_require__(21605);
var SerializedMmlVisitor_js_1 = __webpack_require__(24616);
var Options_js_1 = __webpack_require__(4498);
var LimitedMmlVisitor = (function (_super) {
    __extends(LimitedMmlVisitor, _super);
    function LimitedMmlVisitor() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    LimitedMmlVisitor.prototype.getAttributes = function (node) {
        return _super.prototype.getAttributes.call(this, node).replace(/ ?id=".*?"/, '');
    };
    return LimitedMmlVisitor;
}(SerializedMmlVisitor_js_1.SerializedMmlVisitor));
exports.LimitedMmlVisitor = LimitedMmlVisitor;
(0, MathItem_js_1.newState)('ASSISTIVEMML', 153);
function AssistiveMmlMathItemMixin(BaseMathItem) {
    return (function (_super) {
        __extends(class_1, _super);
        function class_1() {
            return _super !== null && _super.apply(this, arguments) || this;
        }
        class_1.prototype.assistiveMml = function (document, force) {
            if (force === void 0) { force = false; }
            if (this.state() >= MathItem_js_1.STATE.ASSISTIVEMML)
                return;
            if (!this.isEscaped && (document.options.enableAssistiveMml || force)) {
                var adaptor = document.adaptor;
                var mml = document.toMML(this.root).replace(/\n */g, '').replace(/<!--.*?-->/g, '');
                var mmlNodes = adaptor.firstChild(adaptor.body(adaptor.parse(mml, 'text/html')));
                var node = adaptor.node('mjx-assistive-mml', {
                    unselectable: 'on', display: (this.display ? 'block' : 'inline')
                }, [mmlNodes]);
                adaptor.setAttribute(adaptor.firstChild(this.typesetRoot), 'aria-hidden', 'true');
                adaptor.setStyle(this.typesetRoot, 'position', 'relative');
                adaptor.append(this.typesetRoot, node);
            }
            this.state(MathItem_js_1.STATE.ASSISTIVEMML);
        };
        return class_1;
    }(BaseMathItem));
}
exports.AssistiveMmlMathItemMixin = AssistiveMmlMathItemMixin;
function AssistiveMmlMathDocumentMixin(BaseDocument) {
    var _a;
    return _a = (function (_super) {
            __extends(BaseClass, _super);
            function BaseClass() {
                var args = [];
                for (var _i = 0; _i < arguments.length; _i++) {
                    args[_i] = arguments[_i];
                }
                var _this = _super.apply(this, __spreadArray([], __read(args), false)) || this;
                var CLASS = _this.constructor;
                var ProcessBits = CLASS.ProcessBits;
                if (!ProcessBits.has('assistive-mml')) {
                    ProcessBits.allocate('assistive-mml');
                }
                _this.visitor = new LimitedMmlVisitor(_this.mmlFactory);
                _this.options.MathItem =
                    AssistiveMmlMathItemMixin(_this.options.MathItem);
                if ('addStyles' in _this) {
                    _this.addStyles(CLASS.assistiveStyles);
                }
                return _this;
            }
            BaseClass.prototype.toMML = function (node) {
                return this.visitor.visitTree(node);
            };
            BaseClass.prototype.assistiveMml = function () {
                var e_1, _a;
                if (!this.processed.isSet('assistive-mml')) {
                    try {
                        for (var _b = __values(this.math), _c = _b.next(); !_c.done; _c = _b.next()) {
                            var math = _c.value;
                            math.assistiveMml(this);
                        }
                    }
                    catch (e_1_1) { e_1 = { error: e_1_1 }; }
                    finally {
                        try {
                            if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
                        }
                        finally { if (e_1) throw e_1.error; }
                    }
                    this.processed.set('assistive-mml');
                }
                return this;
            };
            BaseClass.prototype.state = function (state, restore) {
                if (restore === void 0) { restore = false; }
                _super.prototype.state.call(this, state, restore);
                if (state < MathItem_js_1.STATE.ASSISTIVEMML) {
                    this.processed.clear('assistive-mml');
                }
                return this;
            };
            return BaseClass;
        }(BaseDocument)),
        _a.OPTIONS = __assign(__assign({}, BaseDocument.OPTIONS), { enableAssistiveMml: true, renderActions: (0, Options_js_1.expandable)(__assign(__assign({}, BaseDocument.OPTIONS.renderActions), { assistiveMml: [MathItem_js_1.STATE.ASSISTIVEMML] })) }),
        _a.assistiveStyles = {
            'mjx-assistive-mml': {
                position: 'absolute !important',
                top: '0px', left: '0px',
                clip: 'rect(1px, 1px, 1px, 1px)',
                padding: '1px 0px 0px 0px !important',
                border: '0px !important',
                display: 'block !important',
                width: 'auto !important',
                overflow: 'hidden !important',
                '-webkit-touch-callout': 'none',
                '-webkit-user-select': 'none',
                '-khtml-user-select': 'none',
                '-moz-user-select': 'none',
                '-ms-user-select': 'none',
                'user-select': 'none'
            },
            'mjx-assistive-mml[display="block"]': {
                width: '100% !important'
            }
        },
        _a;
}
exports.AssistiveMmlMathDocumentMixin = AssistiveMmlMathDocumentMixin;
function AssistiveMmlHandler(handler) {
    handler.documentClass =
        AssistiveMmlMathDocumentMixin(handler.documentClass);
    return handler;
}
exports.AssistiveMmlHandler = AssistiveMmlHandler;
//# sourceMappingURL=assistive-mml.js.map

/***/ }),

/***/ 35659:
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
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.MmlVisitor = void 0;
var MmlFactory_js_1 = __webpack_require__(72666);
var Visitor_js_1 = __webpack_require__(93281);
var MmlVisitor = (function (_super) {
    __extends(MmlVisitor, _super);
    function MmlVisitor(factory) {
        if (factory === void 0) { factory = null; }
        if (!factory) {
            factory = new MmlFactory_js_1.MmlFactory();
        }
        return _super.call(this, factory) || this;
    }
    MmlVisitor.prototype.visitTextNode = function (_node) {
        var _args = [];
        for (var _i = 1; _i < arguments.length; _i++) {
            _args[_i - 1] = arguments[_i];
        }
    };
    MmlVisitor.prototype.visitXMLNode = function (_node) {
        var _args = [];
        for (var _i = 1; _i < arguments.length; _i++) {
            _args[_i - 1] = arguments[_i];
        }
    };
    return MmlVisitor;
}(Visitor_js_1.AbstractVisitor));
exports.MmlVisitor = MmlVisitor;
//# sourceMappingURL=MmlVisitor.js.map

/***/ }),

/***/ 24616:
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
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.SerializedMmlVisitor = exports.toEntity = exports.DATAMJX = void 0;
var MmlVisitor_js_1 = __webpack_require__(35659);
var MmlNode_js_1 = __webpack_require__(83045);
var mi_js_1 = __webpack_require__(91324);
exports.DATAMJX = 'data-mjx-';
var toEntity = function (c) { return '&#x' + c.codePointAt(0).toString(16).toUpperCase() + ';'; };
exports.toEntity = toEntity;
var SerializedMmlVisitor = (function (_super) {
    __extends(SerializedMmlVisitor, _super);
    function SerializedMmlVisitor() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    SerializedMmlVisitor.prototype.visitTree = function (node) {
        return this.visitNode(node, '');
    };
    SerializedMmlVisitor.prototype.visitTextNode = function (node, _space) {
        return this.quoteHTML(node.getText());
    };
    SerializedMmlVisitor.prototype.visitXMLNode = function (node, space) {
        return space + node.getSerializedXML();
    };
    SerializedMmlVisitor.prototype.visitInferredMrowNode = function (node, space) {
        var e_1, _a;
        var mml = [];
        try {
            for (var _b = __values(node.childNodes), _c = _b.next(); !_c.done; _c = _b.next()) {
                var child = _c.value;
                mml.push(this.visitNode(child, space));
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_1) throw e_1.error; }
        }
        return mml.join('\n');
    };
    SerializedMmlVisitor.prototype.visitTeXAtomNode = function (node, space) {
        var children = this.childNodeMml(node, space + '  ', '\n');
        var mml = space + '<mrow' + this.getAttributes(node) + '>' +
            (children.match(/\S/) ? '\n' + children + space : '') + '</mrow>';
        return mml;
    };
    SerializedMmlVisitor.prototype.visitAnnotationNode = function (node, space) {
        return space + '<annotation' + this.getAttributes(node) + '>'
            + this.childNodeMml(node, '', '')
            + '</annotation>';
    };
    SerializedMmlVisitor.prototype.visitDefault = function (node, space) {
        var kind = node.kind;
        var _a = __read((node.isToken || node.childNodes.length === 0 ? ['', ''] : ['\n', space]), 2), nl = _a[0], endspace = _a[1];
        var children = this.childNodeMml(node, space + '  ', nl);
        return space + '<' + kind + this.getAttributes(node) + '>'
            + (children.match(/\S/) ? nl + children + endspace : '')
            + '</' + kind + '>';
    };
    SerializedMmlVisitor.prototype.childNodeMml = function (node, space, nl) {
        var e_2, _a;
        var mml = '';
        try {
            for (var _b = __values(node.childNodes), _c = _b.next(); !_c.done; _c = _b.next()) {
                var child = _c.value;
                mml += this.visitNode(child, space) + nl;
            }
        }
        catch (e_2_1) { e_2 = { error: e_2_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_2) throw e_2.error; }
        }
        return mml;
    };
    SerializedMmlVisitor.prototype.getAttributes = function (node) {
        var e_3, _a;
        var attr = [];
        var defaults = this.constructor.defaultAttributes[node.kind] || {};
        var attributes = Object.assign({}, defaults, this.getDataAttributes(node), node.attributes.getAllAttributes());
        var variants = this.constructor.variants;
        if (attributes.hasOwnProperty('mathvariant') && variants.hasOwnProperty(attributes.mathvariant)) {
            attributes.mathvariant = variants[attributes.mathvariant];
        }
        try {
            for (var _b = __values(Object.keys(attributes)), _c = _b.next(); !_c.done; _c = _b.next()) {
                var name_1 = _c.value;
                var value = String(attributes[name_1]);
                if (value === undefined)
                    continue;
                attr.push(name_1 + '="' + this.quoteHTML(value) + '"');
            }
        }
        catch (e_3_1) { e_3 = { error: e_3_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_3) throw e_3.error; }
        }
        return attr.length ? ' ' + attr.join(' ') : '';
    };
    SerializedMmlVisitor.prototype.getDataAttributes = function (node) {
        var data = {};
        var variant = node.attributes.getExplicit('mathvariant');
        var variants = this.constructor.variants;
        variant && variants.hasOwnProperty(variant) && this.setDataAttribute(data, 'variant', variant);
        node.getProperty('variantForm') && this.setDataAttribute(data, 'alternate', '1');
        node.getProperty('pseudoscript') && this.setDataAttribute(data, 'pseudoscript', 'true');
        node.getProperty('autoOP') === false && this.setDataAttribute(data, 'auto-op', 'false');
        var scriptalign = node.getProperty('scriptalign');
        scriptalign && this.setDataAttribute(data, 'script-align', scriptalign);
        var texclass = node.getProperty('texClass');
        if (texclass !== undefined) {
            var setclass = true;
            if (texclass === MmlNode_js_1.TEXCLASS.OP && node.isKind('mi')) {
                var name_2 = node.getText();
                setclass = !(name_2.length > 1 && name_2.match(mi_js_1.MmlMi.operatorName));
            }
            setclass && this.setDataAttribute(data, 'texclass', texclass < 0 ? 'NONE' : MmlNode_js_1.TEXCLASSNAMES[texclass]);
        }
        node.getProperty('scriptlevel') && node.getProperty('useHeight') === false &&
            this.setDataAttribute(data, 'smallmatrix', 'true');
        return data;
    };
    SerializedMmlVisitor.prototype.setDataAttribute = function (data, name, value) {
        data[exports.DATAMJX + name] = value;
    };
    SerializedMmlVisitor.prototype.quoteHTML = function (value) {
        return value
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;').replace(/>/g, '&gt;')
            .replace(/\"/g, '&quot;')
            .replace(/[\uD800-\uDBFF]./g, exports.toEntity)
            .replace(/[\u0080-\uD7FF\uE000-\uFFFF]/g, exports.toEntity);
    };
    SerializedMmlVisitor.variants = {
        '-tex-calligraphic': 'script',
        '-tex-bold-calligraphic': 'bold-script',
        '-tex-oldstyle': 'normal',
        '-tex-bold-oldstyle': 'bold',
        '-tex-mathit': 'italic'
    };
    SerializedMmlVisitor.defaultAttributes = {
        math: {
            xmlns: 'http://www.w3.org/1998/Math/MathML'
        }
    };
    return SerializedMmlVisitor;
}(MmlVisitor_js_1.MmlVisitor));
exports.SerializedMmlVisitor = SerializedMmlVisitor;
//# sourceMappingURL=SerializedMmlVisitor.js.map

/***/ }),

/***/ 93281:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


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
exports.AbstractVisitor = void 0;
var Node_js_1 = __webpack_require__(85403);
var AbstractVisitor = (function () {
    function AbstractVisitor(factory) {
        var e_1, _a;
        this.nodeHandlers = new Map();
        try {
            for (var _b = __values(factory.getKinds()), _c = _b.next(); !_c.done; _c = _b.next()) {
                var kind = _c.value;
                var method = this[AbstractVisitor.methodName(kind)];
                if (method) {
                    this.nodeHandlers.set(kind, method);
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
    }
    AbstractVisitor.methodName = function (kind) {
        return 'visit' + (kind.charAt(0).toUpperCase() + kind.substr(1)).replace(/[^a-z0-9_]/ig, '_') + 'Node';
    };
    AbstractVisitor.prototype.visitTree = function (tree) {
        var args = [];
        for (var _i = 1; _i < arguments.length; _i++) {
            args[_i - 1] = arguments[_i];
        }
        return this.visitNode.apply(this, __spreadArray([tree], __read(args), false));
    };
    AbstractVisitor.prototype.visitNode = function (node) {
        var args = [];
        for (var _i = 1; _i < arguments.length; _i++) {
            args[_i - 1] = arguments[_i];
        }
        var handler = this.nodeHandlers.get(node.kind) || this.visitDefault;
        return handler.call.apply(handler, __spreadArray([this, node], __read(args), false));
    };
    AbstractVisitor.prototype.visitDefault = function (node) {
        var e_2, _a;
        var args = [];
        for (var _i = 1; _i < arguments.length; _i++) {
            args[_i - 1] = arguments[_i];
        }
        if (node instanceof Node_js_1.AbstractNode) {
            try {
                for (var _b = __values(node.childNodes), _c = _b.next(); !_c.done; _c = _b.next()) {
                    var child = _c.value;
                    this.visitNode.apply(this, __spreadArray([child], __read(args), false));
                }
            }
            catch (e_2_1) { e_2 = { error: e_2_1 }; }
            finally {
                try {
                    if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
                }
                finally { if (e_2) throw e_2.error; }
            }
        }
    };
    AbstractVisitor.prototype.setNodeHandler = function (kind, handler) {
        this.nodeHandlers.set(kind, handler);
    };
    AbstractVisitor.prototype.removeNodeHandler = function (kind) {
        this.nodeHandlers.delete(kind);
    };
    return AbstractVisitor;
}());
exports.AbstractVisitor = AbstractVisitor;
//# sourceMappingURL=Visitor.js.map

/***/ })

}]);
//# sourceMappingURL=4105.5144c29f0bbce103fec4.js.map?v=5144c29f0bbce103fec4
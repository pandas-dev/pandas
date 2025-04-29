"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[7154],{

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

/***/ }),

/***/ 12443:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


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
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.NodeFactory = void 0;
var NodeUtil_js_1 = __importDefault(__webpack_require__(53972));
var NodeFactory = (function () {
    function NodeFactory() {
        this.mmlFactory = null;
        this.factory = { 'node': NodeFactory.createNode,
            'token': NodeFactory.createToken,
            'text': NodeFactory.createText,
            'error': NodeFactory.createError
        };
    }
    NodeFactory.createNode = function (factory, kind, children, def, text) {
        if (children === void 0) { children = []; }
        if (def === void 0) { def = {}; }
        var node = factory.mmlFactory.create(kind);
        node.setChildren(children);
        if (text) {
            node.appendChild(text);
        }
        NodeUtil_js_1.default.setProperties(node, def);
        return node;
    };
    NodeFactory.createToken = function (factory, kind, def, text) {
        if (def === void 0) { def = {}; }
        if (text === void 0) { text = ''; }
        var textNode = factory.create('text', text);
        return factory.create('node', kind, [], def, textNode);
    };
    NodeFactory.createText = function (factory, text) {
        if (text == null) {
            return null;
        }
        return factory.mmlFactory.create('text').setText(text);
    };
    NodeFactory.createError = function (factory, message) {
        var text = factory.create('text', message);
        var mtext = factory.create('node', 'mtext', [], {}, text);
        var error = factory.create('node', 'merror', [mtext], { 'data-mjx-error': message });
        return error;
    };
    NodeFactory.prototype.setMmlFactory = function (mmlFactory) {
        this.mmlFactory = mmlFactory;
    };
    NodeFactory.prototype.set = function (kind, func) {
        this.factory[kind] = func;
    };
    NodeFactory.prototype.setCreators = function (maps) {
        for (var kind in maps) {
            this.set(kind, maps[kind]);
        }
    };
    NodeFactory.prototype.create = function (kind) {
        var rest = [];
        for (var _i = 1; _i < arguments.length; _i++) {
            rest[_i - 1] = arguments[_i];
        }
        var func = this.factory[kind] || this.factory['node'];
        var node = func.apply(void 0, __spreadArray([this, rest[0]], __read(rest.slice(1)), false));
        if (kind === 'node') {
            this.configuration.addNode(rest[0], node);
        }
        return node;
    };
    NodeFactory.prototype.get = function (kind) {
        return this.factory[kind];
    };
    return NodeFactory;
}());
exports.NodeFactory = NodeFactory;
//# sourceMappingURL=NodeFactory.js.map

/***/ }),

/***/ 2362:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


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
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
var NodeUtil_js_1 = __importDefault(__webpack_require__(53972));
var TexConstants_js_1 = __webpack_require__(98927);
var ParseUtil_js_1 = __importDefault(__webpack_require__(55038));
var ParseMethods;
(function (ParseMethods) {
    function variable(parser, c) {
        var def = ParseUtil_js_1.default.getFontDef(parser);
        var env = parser.stack.env;
        if (env.multiLetterIdentifiers && env.font !== '') {
            c = parser.string.substr(parser.i - 1).match(env.multiLetterIdentifiers)[0];
            parser.i += c.length - 1;
            if (def.mathvariant === TexConstants_js_1.TexConstant.Variant.NORMAL && env.noAutoOP && c.length > 1) {
                def.autoOP = false;
            }
        }
        var node = parser.create('token', 'mi', def, c);
        parser.Push(node);
    }
    ParseMethods.variable = variable;
    function digit(parser, c) {
        var mml;
        var pattern = parser.configuration.options['digits'];
        var n = parser.string.slice(parser.i - 1).match(pattern);
        var def = ParseUtil_js_1.default.getFontDef(parser);
        if (n) {
            mml = parser.create('token', 'mn', def, n[0].replace(/[{}]/g, ''));
            parser.i += n[0].length - 1;
        }
        else {
            mml = parser.create('token', 'mo', def, c);
        }
        parser.Push(mml);
    }
    ParseMethods.digit = digit;
    function controlSequence(parser, _c) {
        var name = parser.GetCS();
        parser.parse('macro', [parser, name]);
    }
    ParseMethods.controlSequence = controlSequence;
    function mathchar0mi(parser, mchar) {
        var def = mchar.attributes || { mathvariant: TexConstants_js_1.TexConstant.Variant.ITALIC };
        var node = parser.create('token', 'mi', def, mchar.char);
        parser.Push(node);
    }
    ParseMethods.mathchar0mi = mathchar0mi;
    function mathchar0mo(parser, mchar) {
        var def = mchar.attributes || {};
        def['stretchy'] = false;
        var node = parser.create('token', 'mo', def, mchar.char);
        NodeUtil_js_1.default.setProperty(node, 'fixStretchy', true);
        parser.configuration.addNode('fixStretchy', node);
        parser.Push(node);
    }
    ParseMethods.mathchar0mo = mathchar0mo;
    function mathchar7(parser, mchar) {
        var def = mchar.attributes || { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL };
        if (parser.stack.env['font']) {
            def['mathvariant'] = parser.stack.env['font'];
        }
        var node = parser.create('token', 'mi', def, mchar.char);
        parser.Push(node);
    }
    ParseMethods.mathchar7 = mathchar7;
    function delimiter(parser, delim) {
        var def = delim.attributes || {};
        def = Object.assign({ fence: false, stretchy: false }, def);
        var node = parser.create('token', 'mo', def, delim.char);
        parser.Push(node);
    }
    ParseMethods.delimiter = delimiter;
    function environment(parser, env, func, args) {
        var end = args[0];
        var mml = parser.itemFactory.create('begin').setProperties({ name: env, end: end });
        mml = func.apply(void 0, __spreadArray([parser, mml], __read(args.slice(1)), false));
        parser.Push(mml);
    }
    ParseMethods.environment = environment;
})(ParseMethods || (ParseMethods = {}));
exports["default"] = ParseMethods;
//# sourceMappingURL=ParseMethods.js.map

/***/ }),

/***/ 55661:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


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
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
var StackItemFactory_js_1 = __importDefault(__webpack_require__(42757));
var NodeFactory_js_1 = __webpack_require__(12443);
var NodeUtil_js_1 = __importDefault(__webpack_require__(53972));
var Options_js_1 = __webpack_require__(4498);
var ParseOptions = (function () {
    function ParseOptions(configuration, options) {
        if (options === void 0) { options = []; }
        this.options = {};
        this.packageData = new Map();
        this.parsers = [];
        this.root = null;
        this.nodeLists = {};
        this.error = false;
        this.handlers = configuration.handlers;
        this.nodeFactory = new NodeFactory_js_1.NodeFactory();
        this.nodeFactory.configuration = this;
        this.nodeFactory.setCreators(configuration.nodes);
        this.itemFactory = new StackItemFactory_js_1.default(configuration.items);
        this.itemFactory.configuration = this;
        Options_js_1.defaultOptions.apply(void 0, __spreadArray([this.options], __read(options), false));
        (0, Options_js_1.defaultOptions)(this.options, configuration.options);
    }
    ParseOptions.prototype.pushParser = function (parser) {
        this.parsers.unshift(parser);
    };
    ParseOptions.prototype.popParser = function () {
        this.parsers.shift();
    };
    Object.defineProperty(ParseOptions.prototype, "parser", {
        get: function () {
            return this.parsers[0];
        },
        enumerable: false,
        configurable: true
    });
    ParseOptions.prototype.clear = function () {
        this.parsers = [];
        this.root = null;
        this.nodeLists = {};
        this.error = false;
        this.tags.resetTag();
    };
    ParseOptions.prototype.addNode = function (property, node) {
        var list = this.nodeLists[property];
        if (!list) {
            list = this.nodeLists[property] = [];
        }
        list.push(node);
        if (node.kind !== property) {
            var inlists = (NodeUtil_js_1.default.getProperty(node, 'in-lists') || '');
            var lists = (inlists ? inlists.split(/,/) : []).concat(property).join(',');
            NodeUtil_js_1.default.setProperty(node, 'in-lists', lists);
        }
    };
    ParseOptions.prototype.getList = function (property) {
        var e_1, _a;
        var list = this.nodeLists[property] || [];
        var result = [];
        try {
            for (var list_1 = __values(list), list_1_1 = list_1.next(); !list_1_1.done; list_1_1 = list_1.next()) {
                var node = list_1_1.value;
                if (this.inTree(node)) {
                    result.push(node);
                }
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (list_1_1 && !list_1_1.done && (_a = list_1.return)) _a.call(list_1);
            }
            finally { if (e_1) throw e_1.error; }
        }
        this.nodeLists[property] = result;
        return result;
    };
    ParseOptions.prototype.removeFromList = function (property, nodes) {
        var e_2, _a;
        var list = this.nodeLists[property] || [];
        try {
            for (var nodes_1 = __values(nodes), nodes_1_1 = nodes_1.next(); !nodes_1_1.done; nodes_1_1 = nodes_1.next()) {
                var node = nodes_1_1.value;
                var i = list.indexOf(node);
                if (i >= 0) {
                    list.splice(i, 1);
                }
            }
        }
        catch (e_2_1) { e_2 = { error: e_2_1 }; }
        finally {
            try {
                if (nodes_1_1 && !nodes_1_1.done && (_a = nodes_1.return)) _a.call(nodes_1);
            }
            finally { if (e_2) throw e_2.error; }
        }
    };
    ParseOptions.prototype.inTree = function (node) {
        while (node && node !== this.root) {
            node = node.parent;
        }
        return !!node;
    };
    return ParseOptions;
}());
exports["default"] = ParseOptions;
//# sourceMappingURL=ParseOptions.js.map

/***/ }),

/***/ 26539:
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
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.BaseItem = exports.MmlStack = void 0;
var TexError_js_1 = __importDefault(__webpack_require__(54420));
var MmlStack = (function () {
    function MmlStack(_nodes) {
        this._nodes = _nodes;
    }
    Object.defineProperty(MmlStack.prototype, "nodes", {
        get: function () {
            return this._nodes;
        },
        enumerable: false,
        configurable: true
    });
    MmlStack.prototype.Push = function () {
        var _a;
        var nodes = [];
        for (var _i = 0; _i < arguments.length; _i++) {
            nodes[_i] = arguments[_i];
        }
        (_a = this._nodes).push.apply(_a, __spreadArray([], __read(nodes), false));
    };
    MmlStack.prototype.Pop = function () {
        return this._nodes.pop();
    };
    Object.defineProperty(MmlStack.prototype, "First", {
        get: function () {
            return this._nodes[this.Size() - 1];
        },
        set: function (node) {
            this._nodes[this.Size() - 1] = node;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlStack.prototype, "Last", {
        get: function () {
            return this._nodes[0];
        },
        set: function (node) {
            this._nodes[0] = node;
        },
        enumerable: false,
        configurable: true
    });
    MmlStack.prototype.Peek = function (n) {
        if (n == null) {
            n = 1;
        }
        return this._nodes.slice(this.Size() - n);
    };
    MmlStack.prototype.Size = function () {
        return this._nodes.length;
    };
    MmlStack.prototype.Clear = function () {
        this._nodes = [];
    };
    MmlStack.prototype.toMml = function (inferred, forceRow) {
        if (inferred === void 0) { inferred = true; }
        if (this._nodes.length === 1 && !forceRow) {
            return this.First;
        }
        return this.create('node', inferred ? 'inferredMrow' : 'mrow', this._nodes, {});
    };
    MmlStack.prototype.create = function (kind) {
        var _a;
        var rest = [];
        for (var _i = 1; _i < arguments.length; _i++) {
            rest[_i - 1] = arguments[_i];
        }
        return (_a = this.factory.configuration.nodeFactory).create.apply(_a, __spreadArray([kind], __read(rest), false));
    };
    return MmlStack;
}());
exports.MmlStack = MmlStack;
var BaseItem = (function (_super) {
    __extends(BaseItem, _super);
    function BaseItem(factory) {
        var nodes = [];
        for (var _i = 1; _i < arguments.length; _i++) {
            nodes[_i - 1] = arguments[_i];
        }
        var _this = _super.call(this, nodes) || this;
        _this.factory = factory;
        _this.global = {};
        _this._properties = {};
        if (_this.isOpen) {
            _this._env = {};
        }
        return _this;
    }
    Object.defineProperty(BaseItem.prototype, "kind", {
        get: function () {
            return 'base';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(BaseItem.prototype, "env", {
        get: function () {
            return this._env;
        },
        set: function (value) {
            this._env = value;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(BaseItem.prototype, "copyEnv", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    BaseItem.prototype.getProperty = function (key) {
        return this._properties[key];
    };
    BaseItem.prototype.setProperty = function (key, value) {
        this._properties[key] = value;
        return this;
    };
    Object.defineProperty(BaseItem.prototype, "isOpen", {
        get: function () {
            return false;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(BaseItem.prototype, "isClose", {
        get: function () {
            return false;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(BaseItem.prototype, "isFinal", {
        get: function () {
            return false;
        },
        enumerable: false,
        configurable: true
    });
    BaseItem.prototype.isKind = function (kind) {
        return kind === this.kind;
    };
    BaseItem.prototype.checkItem = function (item) {
        if (item.isKind('over') && this.isOpen) {
            item.setProperty('num', this.toMml(false));
            this.Clear();
        }
        if (item.isKind('cell') && this.isOpen) {
            if (item.getProperty('linebreak')) {
                return BaseItem.fail;
            }
            throw new TexError_js_1.default('Misplaced', 'Misplaced %1', item.getName());
        }
        if (item.isClose && this.getErrors(item.kind)) {
            var _a = __read(this.getErrors(item.kind), 2), id = _a[0], message = _a[1];
            throw new TexError_js_1.default(id, message, item.getName());
        }
        if (!item.isFinal) {
            return BaseItem.success;
        }
        this.Push(item.First);
        return BaseItem.fail;
    };
    BaseItem.prototype.clearEnv = function () {
        var e_1, _a;
        try {
            for (var _b = __values(Object.keys(this.env)), _c = _b.next(); !_c.done; _c = _b.next()) {
                var id = _c.value;
                delete this.env[id];
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_1) throw e_1.error; }
        }
    };
    BaseItem.prototype.setProperties = function (def) {
        Object.assign(this._properties, def);
        return this;
    };
    BaseItem.prototype.getName = function () {
        return this.getProperty('name');
    };
    BaseItem.prototype.toString = function () {
        return this.kind + '[' + this.nodes.join('; ') + ']';
    };
    BaseItem.prototype.getErrors = function (kind) {
        var CLASS = this.constructor;
        return (CLASS.errors || {})[kind] || BaseItem.errors[kind];
    };
    BaseItem.fail = [null, false];
    BaseItem.success = [null, true];
    BaseItem.errors = {
        end: ['MissingBeginExtraEnd', 'Missing \\begin{%1} or extra \\end{%1}'],
        close: ['ExtraCloseMissingOpen', 'Extra close brace or missing open brace'],
        right: ['MissingLeftExtraRight', 'Missing \\left or extra \\right'],
        middle: ['ExtraMiddle', 'Extra \\middle']
    };
    return BaseItem;
}(MmlStack));
exports.BaseItem = BaseItem;
//# sourceMappingURL=StackItem.js.map

/***/ }),

/***/ 42757:
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
var _a;
Object.defineProperty(exports, "__esModule", ({ value: true }));
var StackItem_js_1 = __webpack_require__(26539);
var Factory_js_1 = __webpack_require__(72901);
var DummyItem = (function (_super) {
    __extends(DummyItem, _super);
    function DummyItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    return DummyItem;
}(StackItem_js_1.BaseItem));
var StackItemFactory = (function (_super) {
    __extends(StackItemFactory, _super);
    function StackItemFactory() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.defaultKind = 'dummy';
        _this.configuration = null;
        return _this;
    }
    StackItemFactory.DefaultStackItems = (_a = {},
        _a[DummyItem.prototype.kind] = DummyItem,
        _a);
    return StackItemFactory;
}(Factory_js_1.AbstractFactory));
exports["default"] = StackItemFactory;
//# sourceMappingURL=StackItemFactory.js.map

/***/ }),

/***/ 98927:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.TexConstant = void 0;
var TexConstant;
(function (TexConstant) {
    TexConstant.Variant = {
        NORMAL: 'normal',
        BOLD: 'bold',
        ITALIC: 'italic',
        BOLDITALIC: 'bold-italic',
        DOUBLESTRUCK: 'double-struck',
        FRAKTUR: 'fraktur',
        BOLDFRAKTUR: 'bold-fraktur',
        SCRIPT: 'script',
        BOLDSCRIPT: 'bold-script',
        SANSSERIF: 'sans-serif',
        BOLDSANSSERIF: 'bold-sans-serif',
        SANSSERIFITALIC: 'sans-serif-italic',
        SANSSERIFBOLDITALIC: 'sans-serif-bold-italic',
        MONOSPACE: 'monospace',
        INITIAL: 'inital',
        TAILED: 'tailed',
        LOOPED: 'looped',
        STRETCHED: 'stretched',
        CALLIGRAPHIC: '-tex-calligraphic',
        BOLDCALLIGRAPHIC: '-tex-bold-calligraphic',
        OLDSTYLE: '-tex-oldstyle',
        BOLDOLDSTYLE: '-tex-bold-oldstyle',
        MATHITALIC: '-tex-mathit'
    };
    TexConstant.Form = {
        PREFIX: 'prefix',
        INFIX: 'infix',
        POSTFIX: 'postfix'
    };
    TexConstant.LineBreak = {
        AUTO: 'auto',
        NEWLINE: 'newline',
        NOBREAK: 'nobreak',
        GOODBREAK: 'goodbreak',
        BADBREAK: 'badbreak'
    };
    TexConstant.LineBreakStyle = {
        BEFORE: 'before',
        AFTER: 'after',
        DUPLICATE: 'duplicate',
        INFIXLINBREAKSTYLE: 'infixlinebreakstyle'
    };
    TexConstant.IndentAlign = {
        LEFT: 'left',
        CENTER: 'center',
        RIGHT: 'right',
        AUTO: 'auto',
        ID: 'id',
        INDENTALIGN: 'indentalign'
    };
    TexConstant.IndentShift = {
        INDENTSHIFT: 'indentshift'
    };
    TexConstant.LineThickness = {
        THIN: 'thin',
        MEDIUM: 'medium',
        THICK: 'thick'
    };
    TexConstant.Notation = {
        LONGDIV: 'longdiv',
        ACTUARIAL: 'actuarial',
        PHASORANGLE: 'phasorangle',
        RADICAL: 'radical',
        BOX: 'box',
        ROUNDEDBOX: 'roundedbox',
        CIRCLE: 'circle',
        LEFT: 'left',
        RIGHT: 'right',
        TOP: 'top',
        BOTTOM: 'bottom',
        UPDIAGONALSTRIKE: 'updiagonalstrike',
        DOWNDIAGONALSTRIKE: 'downdiagonalstrike',
        VERTICALSTRIKE: 'verticalstrike',
        HORIZONTALSTRIKE: 'horizontalstrike',
        NORTHEASTARROW: 'northeastarrow',
        MADRUWB: 'madruwb',
        UPDIAGONALARROW: 'updiagonalarrow'
    };
    TexConstant.Align = {
        TOP: 'top',
        BOTTOM: 'bottom',
        CENTER: 'center',
        BASELINE: 'baseline',
        AXIS: 'axis',
        LEFT: 'left',
        RIGHT: 'right'
    };
    TexConstant.Lines = {
        NONE: 'none',
        SOLID: 'solid',
        DASHED: 'dashed'
    };
    TexConstant.Side = {
        LEFT: 'left',
        RIGHT: 'right',
        LEFTOVERLAP: 'leftoverlap',
        RIGHTOVERLAP: 'rightoverlap'
    };
    TexConstant.Width = {
        AUTO: 'auto',
        FIT: 'fit'
    };
    TexConstant.Actiontype = {
        TOGGLE: 'toggle',
        STATUSLINE: 'statusline',
        TOOLTIP: 'tooltip',
        INPUT: 'input'
    };
    TexConstant.Overflow = {
        LINBREAK: 'linebreak',
        SCROLL: 'scroll',
        ELIDE: 'elide',
        TRUNCATE: 'truncate',
        SCALE: 'scale'
    };
    TexConstant.Unit = {
        EM: 'em',
        EX: 'ex',
        PX: 'px',
        IN: 'in',
        CM: 'cm',
        MM: 'mm',
        PT: 'pt',
        PC: 'pc'
    };
})(TexConstant = exports.TexConstant || (exports.TexConstant = {}));
//# sourceMappingURL=TexConstants.js.map

/***/ }),

/***/ 32742:
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
var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    var desc = Object.getOwnPropertyDescriptor(m, k);
    if (!desc || ("get" in desc ? !m.__esModule : desc.writable || desc.configurable)) {
      desc = { enumerable: true, get: function() { return m[k]; } };
    }
    Object.defineProperty(o, k2, desc);
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
    Object.defineProperty(o, "default", { enumerable: true, value: v });
}) : function(o, v) {
    o["default"] = v;
});
var __importStar = (this && this.__importStar) || function (mod) {
    if (mod && mod.__esModule) return mod;
    var result = {};
    if (mod != null) for (var k in mod) if (k !== "default" && Object.prototype.hasOwnProperty.call(mod, k)) __createBinding(result, mod, k);
    __setModuleDefault(result, mod);
    return result;
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
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
var _a;
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.BaseConfiguration = exports.BaseTags = exports.Other = void 0;
var Configuration_js_1 = __webpack_require__(63401);
var MapHandler_js_1 = __webpack_require__(21683);
var TexError_js_1 = __importDefault(__webpack_require__(54420));
var NodeUtil_js_1 = __importDefault(__webpack_require__(53972));
var SymbolMap_js_1 = __webpack_require__(65695);
var bitem = __importStar(__webpack_require__(31201));
var Tags_js_1 = __webpack_require__(75723);
__webpack_require__(23430);
var OperatorDictionary_js_1 = __webpack_require__(99772);
new SymbolMap_js_1.CharacterMap('remap', null, {
    '-': '\u2212',
    '*': '\u2217',
    '`': '\u2018'
});
function Other(parser, char) {
    var font = parser.stack.env['font'];
    var def = font ?
        { mathvariant: parser.stack.env['font'] } : {};
    var remap = MapHandler_js_1.MapHandler.getMap('remap').lookup(char);
    var range = (0, OperatorDictionary_js_1.getRange)(char);
    var type = (range ? range[3] : 'mo');
    var mo = parser.create('token', type, def, (remap ? remap.char : char));
    range[4] && mo.attributes.set('mathvariant', range[4]);
    if (type === 'mo') {
        NodeUtil_js_1.default.setProperty(mo, 'fixStretchy', true);
        parser.configuration.addNode('fixStretchy', mo);
    }
    parser.Push(mo);
}
exports.Other = Other;
function csUndefined(_parser, name) {
    throw new TexError_js_1.default('UndefinedControlSequence', 'Undefined control sequence %1', '\\' + name);
}
function envUndefined(_parser, env) {
    throw new TexError_js_1.default('UnknownEnv', 'Unknown environment \'%1\'', env);
}
function filterNonscript(_a) {
    var e_1, _b;
    var data = _a.data;
    try {
        for (var _c = __values(data.getList('nonscript')), _d = _c.next(); !_d.done; _d = _c.next()) {
            var mml = _d.value;
            if (mml.attributes.get('scriptlevel') > 0) {
                var parent_1 = mml.parent;
                parent_1.childNodes.splice(parent_1.childIndex(mml), 1);
                data.removeFromList(mml.kind, [mml]);
                if (mml.isKind('mrow')) {
                    var mstyle = mml.childNodes[0];
                    data.removeFromList('mstyle', [mstyle]);
                    data.removeFromList('mspace', mstyle.childNodes[0].childNodes);
                }
            }
            else if (mml.isKind('mrow')) {
                mml.parent.replaceChild(mml.childNodes[0], mml);
                data.removeFromList('mrow', [mml]);
            }
        }
    }
    catch (e_1_1) { e_1 = { error: e_1_1 }; }
    finally {
        try {
            if (_d && !_d.done && (_b = _c.return)) _b.call(_c);
        }
        finally { if (e_1) throw e_1.error; }
    }
}
var BaseTags = (function (_super) {
    __extends(BaseTags, _super);
    function BaseTags() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    return BaseTags;
}(Tags_js_1.AbstractTags));
exports.BaseTags = BaseTags;
exports.BaseConfiguration = Configuration_js_1.Configuration.create('base', {
    handler: {
        character: ['command', 'special', 'letter', 'digit'],
        delimiter: ['delimiter'],
        macro: ['delimiter', 'macros', 'mathchar0mi', 'mathchar0mo', 'mathchar7'],
        environment: ['environment']
    },
    fallback: {
        character: Other,
        macro: csUndefined,
        environment: envUndefined
    },
    items: (_a = {},
        _a[bitem.StartItem.prototype.kind] = bitem.StartItem,
        _a[bitem.StopItem.prototype.kind] = bitem.StopItem,
        _a[bitem.OpenItem.prototype.kind] = bitem.OpenItem,
        _a[bitem.CloseItem.prototype.kind] = bitem.CloseItem,
        _a[bitem.PrimeItem.prototype.kind] = bitem.PrimeItem,
        _a[bitem.SubsupItem.prototype.kind] = bitem.SubsupItem,
        _a[bitem.OverItem.prototype.kind] = bitem.OverItem,
        _a[bitem.LeftItem.prototype.kind] = bitem.LeftItem,
        _a[bitem.Middle.prototype.kind] = bitem.Middle,
        _a[bitem.RightItem.prototype.kind] = bitem.RightItem,
        _a[bitem.BeginItem.prototype.kind] = bitem.BeginItem,
        _a[bitem.EndItem.prototype.kind] = bitem.EndItem,
        _a[bitem.StyleItem.prototype.kind] = bitem.StyleItem,
        _a[bitem.PositionItem.prototype.kind] = bitem.PositionItem,
        _a[bitem.CellItem.prototype.kind] = bitem.CellItem,
        _a[bitem.MmlItem.prototype.kind] = bitem.MmlItem,
        _a[bitem.FnItem.prototype.kind] = bitem.FnItem,
        _a[bitem.NotItem.prototype.kind] = bitem.NotItem,
        _a[bitem.NonscriptItem.prototype.kind] = bitem.NonscriptItem,
        _a[bitem.DotsItem.prototype.kind] = bitem.DotsItem,
        _a[bitem.ArrayItem.prototype.kind] = bitem.ArrayItem,
        _a[bitem.EqnArrayItem.prototype.kind] = bitem.EqnArrayItem,
        _a[bitem.EquationItem.prototype.kind] = bitem.EquationItem,
        _a),
    options: {
        maxMacros: 1000,
        baseURL: (typeof (document) === 'undefined' ||
            document.getElementsByTagName('base').length === 0) ?
            '' : String(document.location).replace(/#.*$/, '')
    },
    tags: {
        base: BaseTags
    },
    postprocessors: [[filterNonscript, -4]]
});
//# sourceMappingURL=BaseConfiguration.js.map

/***/ }),

/***/ 31201:
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
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.EquationItem = exports.EqnArrayItem = exports.ArrayItem = exports.DotsItem = exports.NonscriptItem = exports.NotItem = exports.FnItem = exports.MmlItem = exports.CellItem = exports.PositionItem = exports.StyleItem = exports.EndItem = exports.BeginItem = exports.RightItem = exports.Middle = exports.LeftItem = exports.OverItem = exports.SubsupItem = exports.PrimeItem = exports.CloseItem = exports.OpenItem = exports.StopItem = exports.StartItem = void 0;
var MapHandler_js_1 = __webpack_require__(21683);
var Entities_js_1 = __webpack_require__(61051);
var MmlNode_js_1 = __webpack_require__(83045);
var TexError_js_1 = __importDefault(__webpack_require__(54420));
var ParseUtil_js_1 = __importDefault(__webpack_require__(55038));
var NodeUtil_js_1 = __importDefault(__webpack_require__(53972));
var StackItem_js_1 = __webpack_require__(26539);
var StartItem = (function (_super) {
    __extends(StartItem, _super);
    function StartItem(factory, global) {
        var _this = _super.call(this, factory) || this;
        _this.global = global;
        return _this;
    }
    Object.defineProperty(StartItem.prototype, "kind", {
        get: function () {
            return 'start';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(StartItem.prototype, "isOpen", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    StartItem.prototype.checkItem = function (item) {
        if (item.isKind('stop')) {
            var node = this.toMml();
            if (!this.global.isInner) {
                node = this.factory.configuration.tags.finalize(node, this.env);
            }
            return [[this.factory.create('mml', node)], true];
        }
        return _super.prototype.checkItem.call(this, item);
    };
    return StartItem;
}(StackItem_js_1.BaseItem));
exports.StartItem = StartItem;
var StopItem = (function (_super) {
    __extends(StopItem, _super);
    function StopItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(StopItem.prototype, "kind", {
        get: function () {
            return 'stop';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(StopItem.prototype, "isClose", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    return StopItem;
}(StackItem_js_1.BaseItem));
exports.StopItem = StopItem;
var OpenItem = (function (_super) {
    __extends(OpenItem, _super);
    function OpenItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(OpenItem.prototype, "kind", {
        get: function () {
            return 'open';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(OpenItem.prototype, "isOpen", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    OpenItem.prototype.checkItem = function (item) {
        if (item.isKind('close')) {
            var mml = this.toMml();
            var node = this.create('node', 'TeXAtom', [mml]);
            return [[this.factory.create('mml', node)], true];
        }
        return _super.prototype.checkItem.call(this, item);
    };
    OpenItem.errors = Object.assign(Object.create(StackItem_js_1.BaseItem.errors), {
        'stop': ['ExtraOpenMissingClose',
            'Extra open brace or missing close brace']
    });
    return OpenItem;
}(StackItem_js_1.BaseItem));
exports.OpenItem = OpenItem;
var CloseItem = (function (_super) {
    __extends(CloseItem, _super);
    function CloseItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(CloseItem.prototype, "kind", {
        get: function () {
            return 'close';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(CloseItem.prototype, "isClose", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    return CloseItem;
}(StackItem_js_1.BaseItem));
exports.CloseItem = CloseItem;
var PrimeItem = (function (_super) {
    __extends(PrimeItem, _super);
    function PrimeItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(PrimeItem.prototype, "kind", {
        get: function () {
            return 'prime';
        },
        enumerable: false,
        configurable: true
    });
    PrimeItem.prototype.checkItem = function (item) {
        var _a = __read(this.Peek(2), 2), top0 = _a[0], top1 = _a[1];
        if (!NodeUtil_js_1.default.isType(top0, 'msubsup') || NodeUtil_js_1.default.isType(top0, 'msup')) {
            var node = this.create('node', 'msup', [top0, top1]);
            return [[node, item], true];
        }
        NodeUtil_js_1.default.setChild(top0, top0.sup, top1);
        return [[top0, item], true];
    };
    return PrimeItem;
}(StackItem_js_1.BaseItem));
exports.PrimeItem = PrimeItem;
var SubsupItem = (function (_super) {
    __extends(SubsupItem, _super);
    function SubsupItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(SubsupItem.prototype, "kind", {
        get: function () {
            return 'subsup';
        },
        enumerable: false,
        configurable: true
    });
    SubsupItem.prototype.checkItem = function (item) {
        if (item.isKind('open') || item.isKind('left')) {
            return StackItem_js_1.BaseItem.success;
        }
        var top = this.First;
        var position = this.getProperty('position');
        if (item.isKind('mml')) {
            if (this.getProperty('primes')) {
                if (position !== 2) {
                    NodeUtil_js_1.default.setChild(top, 2, this.getProperty('primes'));
                }
                else {
                    NodeUtil_js_1.default.setProperty(this.getProperty('primes'), 'variantForm', true);
                    var node = this.create('node', 'mrow', [this.getProperty('primes'), item.First]);
                    item.First = node;
                }
            }
            NodeUtil_js_1.default.setChild(top, position, item.First);
            if (this.getProperty('movesupsub') != null) {
                NodeUtil_js_1.default.setProperty(top, 'movesupsub', this.getProperty('movesupsub'));
            }
            var result = this.factory.create('mml', top);
            return [[result], true];
        }
        if (_super.prototype.checkItem.call(this, item)[1]) {
            var error = this.getErrors(['', 'sub', 'sup'][position]);
            throw new (TexError_js_1.default.bind.apply(TexError_js_1.default, __spreadArray([void 0, error[0], error[1]], __read(error.splice(2)), false)))();
        }
        return null;
    };
    SubsupItem.errors = Object.assign(Object.create(StackItem_js_1.BaseItem.errors), {
        'stop': ['MissingScript',
            'Missing superscript or subscript argument'],
        'sup': ['MissingOpenForSup',
            'Missing open brace for superscript'],
        'sub': ['MissingOpenForSub',
            'Missing open brace for subscript']
    });
    return SubsupItem;
}(StackItem_js_1.BaseItem));
exports.SubsupItem = SubsupItem;
var OverItem = (function (_super) {
    __extends(OverItem, _super);
    function OverItem(factory) {
        var _this = _super.call(this, factory) || this;
        _this.setProperty('name', '\\over');
        return _this;
    }
    Object.defineProperty(OverItem.prototype, "kind", {
        get: function () {
            return 'over';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(OverItem.prototype, "isClose", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    OverItem.prototype.checkItem = function (item) {
        if (item.isKind('over')) {
            throw new TexError_js_1.default('AmbiguousUseOf', 'Ambiguous use of %1', item.getName());
        }
        if (item.isClose) {
            var mml = this.create('node', 'mfrac', [this.getProperty('num'), this.toMml(false)]);
            if (this.getProperty('thickness') != null) {
                NodeUtil_js_1.default.setAttribute(mml, 'linethickness', this.getProperty('thickness'));
            }
            if (this.getProperty('open') || this.getProperty('close')) {
                NodeUtil_js_1.default.setProperty(mml, 'withDelims', true);
                mml = ParseUtil_js_1.default.fixedFence(this.factory.configuration, this.getProperty('open'), mml, this.getProperty('close'));
            }
            return [[this.factory.create('mml', mml), item], true];
        }
        return _super.prototype.checkItem.call(this, item);
    };
    OverItem.prototype.toString = function () {
        return 'over[' + this.getProperty('num') +
            ' / ' + this.nodes.join('; ') + ']';
    };
    return OverItem;
}(StackItem_js_1.BaseItem));
exports.OverItem = OverItem;
var LeftItem = (function (_super) {
    __extends(LeftItem, _super);
    function LeftItem(factory, delim) {
        var _this = _super.call(this, factory) || this;
        _this.setProperty('delim', delim);
        return _this;
    }
    Object.defineProperty(LeftItem.prototype, "kind", {
        get: function () {
            return 'left';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(LeftItem.prototype, "isOpen", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    LeftItem.prototype.checkItem = function (item) {
        if (item.isKind('right')) {
            return [[this.factory.create('mml', ParseUtil_js_1.default.fenced(this.factory.configuration, this.getProperty('delim'), this.toMml(), item.getProperty('delim'), '', item.getProperty('color')))], true];
        }
        if (item.isKind('middle')) {
            var def = { stretchy: true };
            if (item.getProperty('color')) {
                def.mathcolor = item.getProperty('color');
            }
            this.Push(this.create('node', 'TeXAtom', [], { texClass: MmlNode_js_1.TEXCLASS.CLOSE }), this.create('token', 'mo', def, item.getProperty('delim')), this.create('node', 'TeXAtom', [], { texClass: MmlNode_js_1.TEXCLASS.OPEN }));
            this.env = {};
            return [[this], true];
        }
        return _super.prototype.checkItem.call(this, item);
    };
    LeftItem.errors = Object.assign(Object.create(StackItem_js_1.BaseItem.errors), {
        'stop': ['ExtraLeftMissingRight',
            'Extra \\left or missing \\right']
    });
    return LeftItem;
}(StackItem_js_1.BaseItem));
exports.LeftItem = LeftItem;
var Middle = (function (_super) {
    __extends(Middle, _super);
    function Middle(factory, delim, color) {
        var _this = _super.call(this, factory) || this;
        _this.setProperty('delim', delim);
        color && _this.setProperty('color', color);
        return _this;
    }
    Object.defineProperty(Middle.prototype, "kind", {
        get: function () {
            return 'middle';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(Middle.prototype, "isClose", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    return Middle;
}(StackItem_js_1.BaseItem));
exports.Middle = Middle;
var RightItem = (function (_super) {
    __extends(RightItem, _super);
    function RightItem(factory, delim, color) {
        var _this = _super.call(this, factory) || this;
        _this.setProperty('delim', delim);
        color && _this.setProperty('color', color);
        return _this;
    }
    Object.defineProperty(RightItem.prototype, "kind", {
        get: function () {
            return 'right';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(RightItem.prototype, "isClose", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    return RightItem;
}(StackItem_js_1.BaseItem));
exports.RightItem = RightItem;
var BeginItem = (function (_super) {
    __extends(BeginItem, _super);
    function BeginItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(BeginItem.prototype, "kind", {
        get: function () {
            return 'begin';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(BeginItem.prototype, "isOpen", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    BeginItem.prototype.checkItem = function (item) {
        if (item.isKind('end')) {
            if (item.getName() !== this.getName()) {
                throw new TexError_js_1.default('EnvBadEnd', '\\begin{%1} ended with \\end{%2}', this.getName(), item.getName());
            }
            if (!this.getProperty('end')) {
                return [[this.factory.create('mml', this.toMml())], true];
            }
            return StackItem_js_1.BaseItem.fail;
        }
        if (item.isKind('stop')) {
            throw new TexError_js_1.default('EnvMissingEnd', 'Missing \\end{%1}', this.getName());
        }
        return _super.prototype.checkItem.call(this, item);
    };
    return BeginItem;
}(StackItem_js_1.BaseItem));
exports.BeginItem = BeginItem;
var EndItem = (function (_super) {
    __extends(EndItem, _super);
    function EndItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(EndItem.prototype, "kind", {
        get: function () {
            return 'end';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(EndItem.prototype, "isClose", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    return EndItem;
}(StackItem_js_1.BaseItem));
exports.EndItem = EndItem;
var StyleItem = (function (_super) {
    __extends(StyleItem, _super);
    function StyleItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(StyleItem.prototype, "kind", {
        get: function () {
            return 'style';
        },
        enumerable: false,
        configurable: true
    });
    StyleItem.prototype.checkItem = function (item) {
        if (!item.isClose) {
            return _super.prototype.checkItem.call(this, item);
        }
        var mml = this.create('node', 'mstyle', this.nodes, this.getProperty('styles'));
        return [[this.factory.create('mml', mml), item], true];
    };
    return StyleItem;
}(StackItem_js_1.BaseItem));
exports.StyleItem = StyleItem;
var PositionItem = (function (_super) {
    __extends(PositionItem, _super);
    function PositionItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(PositionItem.prototype, "kind", {
        get: function () {
            return 'position';
        },
        enumerable: false,
        configurable: true
    });
    PositionItem.prototype.checkItem = function (item) {
        if (item.isClose) {
            throw new TexError_js_1.default('MissingBoxFor', 'Missing box for %1', this.getName());
        }
        if (item.isFinal) {
            var mml = item.toMml();
            switch (this.getProperty('move')) {
                case 'vertical':
                    mml = this.create('node', 'mpadded', [mml], { height: this.getProperty('dh'),
                        depth: this.getProperty('dd'),
                        voffset: this.getProperty('dh') });
                    return [[this.factory.create('mml', mml)], true];
                case 'horizontal':
                    return [[this.factory.create('mml', this.getProperty('left')), item,
                            this.factory.create('mml', this.getProperty('right'))], true];
            }
        }
        return _super.prototype.checkItem.call(this, item);
    };
    return PositionItem;
}(StackItem_js_1.BaseItem));
exports.PositionItem = PositionItem;
var CellItem = (function (_super) {
    __extends(CellItem, _super);
    function CellItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(CellItem.prototype, "kind", {
        get: function () {
            return 'cell';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(CellItem.prototype, "isClose", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    return CellItem;
}(StackItem_js_1.BaseItem));
exports.CellItem = CellItem;
var MmlItem = (function (_super) {
    __extends(MmlItem, _super);
    function MmlItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(MmlItem.prototype, "isFinal", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(MmlItem.prototype, "kind", {
        get: function () {
            return 'mml';
        },
        enumerable: false,
        configurable: true
    });
    return MmlItem;
}(StackItem_js_1.BaseItem));
exports.MmlItem = MmlItem;
var FnItem = (function (_super) {
    __extends(FnItem, _super);
    function FnItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(FnItem.prototype, "kind", {
        get: function () {
            return 'fn';
        },
        enumerable: false,
        configurable: true
    });
    FnItem.prototype.checkItem = function (item) {
        var top = this.First;
        if (top) {
            if (item.isOpen) {
                return StackItem_js_1.BaseItem.success;
            }
            if (!item.isKind('fn')) {
                var mml = item.First;
                if (!item.isKind('mml') || !mml) {
                    return [[top, item], true];
                }
                if ((NodeUtil_js_1.default.isType(mml, 'mstyle') && mml.childNodes.length &&
                    NodeUtil_js_1.default.isType(mml.childNodes[0].childNodes[0], 'mspace')) ||
                    NodeUtil_js_1.default.isType(mml, 'mspace')) {
                    return [[top, item], true];
                }
                if (NodeUtil_js_1.default.isEmbellished(mml)) {
                    mml = NodeUtil_js_1.default.getCoreMO(mml);
                }
                var form = NodeUtil_js_1.default.getForm(mml);
                if (form != null && [0, 0, 1, 1, 0, 1, 1, 0, 0, 0][form[2]]) {
                    return [[top, item], true];
                }
            }
            var node = this.create('token', 'mo', { texClass: MmlNode_js_1.TEXCLASS.NONE }, Entities_js_1.entities.ApplyFunction);
            return [[top, node, item], true];
        }
        return _super.prototype.checkItem.apply(this, arguments);
    };
    return FnItem;
}(StackItem_js_1.BaseItem));
exports.FnItem = FnItem;
var NotItem = (function (_super) {
    __extends(NotItem, _super);
    function NotItem() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.remap = MapHandler_js_1.MapHandler.getMap('not_remap');
        return _this;
    }
    Object.defineProperty(NotItem.prototype, "kind", {
        get: function () {
            return 'not';
        },
        enumerable: false,
        configurable: true
    });
    NotItem.prototype.checkItem = function (item) {
        var mml;
        var c;
        var textNode;
        if (item.isKind('open') || item.isKind('left')) {
            return StackItem_js_1.BaseItem.success;
        }
        if (item.isKind('mml') &&
            (NodeUtil_js_1.default.isType(item.First, 'mo') || NodeUtil_js_1.default.isType(item.First, 'mi') ||
                NodeUtil_js_1.default.isType(item.First, 'mtext'))) {
            mml = item.First;
            c = NodeUtil_js_1.default.getText(mml);
            if (c.length === 1 && !NodeUtil_js_1.default.getProperty(mml, 'movesupsub') &&
                NodeUtil_js_1.default.getChildren(mml).length === 1) {
                if (this.remap.contains(c)) {
                    textNode = this.create('text', this.remap.lookup(c).char);
                    NodeUtil_js_1.default.setChild(mml, 0, textNode);
                }
                else {
                    textNode = this.create('text', '\u0338');
                    NodeUtil_js_1.default.appendChildren(mml, [textNode]);
                }
                return [[item], true];
            }
        }
        textNode = this.create('text', '\u29F8');
        var mtextNode = this.create('node', 'mtext', [], {}, textNode);
        var paddedNode = this.create('node', 'mpadded', [mtextNode], { width: 0 });
        mml = this.create('node', 'TeXAtom', [paddedNode], { texClass: MmlNode_js_1.TEXCLASS.REL });
        return [[mml, item], true];
    };
    return NotItem;
}(StackItem_js_1.BaseItem));
exports.NotItem = NotItem;
var NonscriptItem = (function (_super) {
    __extends(NonscriptItem, _super);
    function NonscriptItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(NonscriptItem.prototype, "kind", {
        get: function () {
            return 'nonscript';
        },
        enumerable: false,
        configurable: true
    });
    NonscriptItem.prototype.checkItem = function (item) {
        if (item.isKind('mml') && item.Size() === 1) {
            var mml = item.First;
            if (mml.isKind('mstyle') && mml.notParent) {
                mml = NodeUtil_js_1.default.getChildren(NodeUtil_js_1.default.getChildren(mml)[0])[0];
            }
            if (mml.isKind('mspace')) {
                if (mml !== item.First) {
                    var mrow = this.create('node', 'mrow', [item.Pop()]);
                    item.Push(mrow);
                }
                this.factory.configuration.addNode('nonscript', item.First);
            }
        }
        return [[item], true];
    };
    return NonscriptItem;
}(StackItem_js_1.BaseItem));
exports.NonscriptItem = NonscriptItem;
var DotsItem = (function (_super) {
    __extends(DotsItem, _super);
    function DotsItem() {
        return _super !== null && _super.apply(this, arguments) || this;
    }
    Object.defineProperty(DotsItem.prototype, "kind", {
        get: function () {
            return 'dots';
        },
        enumerable: false,
        configurable: true
    });
    DotsItem.prototype.checkItem = function (item) {
        if (item.isKind('open') || item.isKind('left')) {
            return StackItem_js_1.BaseItem.success;
        }
        var dots = this.getProperty('ldots');
        var top = item.First;
        if (item.isKind('mml') && NodeUtil_js_1.default.isEmbellished(top)) {
            var tclass = NodeUtil_js_1.default.getTexClass(NodeUtil_js_1.default.getCoreMO(top));
            if (tclass === MmlNode_js_1.TEXCLASS.BIN || tclass === MmlNode_js_1.TEXCLASS.REL) {
                dots = this.getProperty('cdots');
            }
        }
        return [[dots, item], true];
    };
    return DotsItem;
}(StackItem_js_1.BaseItem));
exports.DotsItem = DotsItem;
var ArrayItem = (function (_super) {
    __extends(ArrayItem, _super);
    function ArrayItem() {
        var _this = _super !== null && _super.apply(this, arguments) || this;
        _this.table = [];
        _this.row = [];
        _this.frame = [];
        _this.hfill = [];
        _this.arraydef = {};
        _this.dashed = false;
        return _this;
    }
    Object.defineProperty(ArrayItem.prototype, "kind", {
        get: function () {
            return 'array';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(ArrayItem.prototype, "isOpen", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(ArrayItem.prototype, "copyEnv", {
        get: function () {
            return false;
        },
        enumerable: false,
        configurable: true
    });
    ArrayItem.prototype.checkItem = function (item) {
        if (item.isClose && !item.isKind('over')) {
            if (item.getProperty('isEntry')) {
                this.EndEntry();
                this.clearEnv();
                return StackItem_js_1.BaseItem.fail;
            }
            if (item.getProperty('isCR')) {
                this.EndEntry();
                this.EndRow();
                this.clearEnv();
                return StackItem_js_1.BaseItem.fail;
            }
            this.EndTable();
            this.clearEnv();
            var newItem = this.factory.create('mml', this.createMml());
            if (this.getProperty('requireClose')) {
                if (item.isKind('close')) {
                    return [[newItem], true];
                }
                throw new TexError_js_1.default('MissingCloseBrace', 'Missing close brace');
            }
            return [[newItem, item], true];
        }
        return _super.prototype.checkItem.call(this, item);
    };
    ArrayItem.prototype.createMml = function () {
        var scriptlevel = this.arraydef['scriptlevel'];
        delete this.arraydef['scriptlevel'];
        var mml = this.create('node', 'mtable', this.table, this.arraydef);
        if (scriptlevel) {
            mml.setProperty('scriptlevel', scriptlevel);
        }
        if (this.frame.length === 4) {
            NodeUtil_js_1.default.setAttribute(mml, 'frame', this.dashed ? 'dashed' : 'solid');
        }
        else if (this.frame.length) {
            if (this.arraydef['rowlines']) {
                this.arraydef['rowlines'] =
                    this.arraydef['rowlines'].replace(/none( none)+$/, 'none');
            }
            NodeUtil_js_1.default.setAttribute(mml, 'frame', '');
            mml = this.create('node', 'menclose', [mml], { notation: this.frame.join(' ') });
            if ((this.arraydef['columnlines'] || 'none') !== 'none' ||
                (this.arraydef['rowlines'] || 'none') !== 'none') {
                NodeUtil_js_1.default.setAttribute(mml, 'data-padding', 0);
            }
        }
        if (this.getProperty('open') || this.getProperty('close')) {
            mml = ParseUtil_js_1.default.fenced(this.factory.configuration, this.getProperty('open'), mml, this.getProperty('close'));
        }
        return mml;
    };
    ArrayItem.prototype.EndEntry = function () {
        var mtd = this.create('node', 'mtd', this.nodes);
        if (this.hfill.length) {
            if (this.hfill[0] === 0) {
                NodeUtil_js_1.default.setAttribute(mtd, 'columnalign', 'right');
            }
            if (this.hfill[this.hfill.length - 1] === this.Size()) {
                NodeUtil_js_1.default.setAttribute(mtd, 'columnalign', NodeUtil_js_1.default.getAttribute(mtd, 'columnalign') ? 'center' : 'left');
            }
        }
        this.row.push(mtd);
        this.Clear();
        this.hfill = [];
    };
    ArrayItem.prototype.EndRow = function () {
        var node;
        if (this.getProperty('isNumbered') && this.row.length === 3) {
            this.row.unshift(this.row.pop());
            node = this.create('node', 'mlabeledtr', this.row);
        }
        else {
            node = this.create('node', 'mtr', this.row);
        }
        this.table.push(node);
        this.row = [];
    };
    ArrayItem.prototype.EndTable = function () {
        if (this.Size() || this.row.length) {
            this.EndEntry();
            this.EndRow();
        }
        this.checkLines();
    };
    ArrayItem.prototype.checkLines = function () {
        if (this.arraydef['rowlines']) {
            var lines = this.arraydef['rowlines'].split(/ /);
            if (lines.length === this.table.length) {
                this.frame.push('bottom');
                lines.pop();
                this.arraydef['rowlines'] = lines.join(' ');
            }
            else if (lines.length < this.table.length - 1) {
                this.arraydef['rowlines'] += ' none';
            }
        }
        if (this.getProperty('rowspacing')) {
            var rows = this.arraydef['rowspacing'].split(/ /);
            while (rows.length < this.table.length) {
                rows.push(this.getProperty('rowspacing') + 'em');
            }
            this.arraydef['rowspacing'] = rows.join(' ');
        }
    };
    ArrayItem.prototype.addRowSpacing = function (spacing) {
        if (this.arraydef['rowspacing']) {
            var rows = this.arraydef['rowspacing'].split(/ /);
            if (!this.getProperty('rowspacing')) {
                var dimem = ParseUtil_js_1.default.dimen2em(rows[0]);
                this.setProperty('rowspacing', dimem);
            }
            var rowspacing = this.getProperty('rowspacing');
            while (rows.length < this.table.length) {
                rows.push(ParseUtil_js_1.default.Em(rowspacing));
            }
            rows[this.table.length - 1] = ParseUtil_js_1.default.Em(Math.max(0, rowspacing + ParseUtil_js_1.default.dimen2em(spacing)));
            this.arraydef['rowspacing'] = rows.join(' ');
        }
    };
    return ArrayItem;
}(StackItem_js_1.BaseItem));
exports.ArrayItem = ArrayItem;
var EqnArrayItem = (function (_super) {
    __extends(EqnArrayItem, _super);
    function EqnArrayItem(factory) {
        var args = [];
        for (var _i = 1; _i < arguments.length; _i++) {
            args[_i - 1] = arguments[_i];
        }
        var _this = _super.call(this, factory) || this;
        _this.maxrow = 0;
        _this.factory.configuration.tags.start(args[0], args[2], args[1]);
        return _this;
    }
    Object.defineProperty(EqnArrayItem.prototype, "kind", {
        get: function () {
            return 'eqnarray';
        },
        enumerable: false,
        configurable: true
    });
    EqnArrayItem.prototype.EndEntry = function () {
        if (this.row.length) {
            ParseUtil_js_1.default.fixInitialMO(this.factory.configuration, this.nodes);
        }
        var node = this.create('node', 'mtd', this.nodes);
        this.row.push(node);
        this.Clear();
    };
    EqnArrayItem.prototype.EndRow = function () {
        if (this.row.length > this.maxrow) {
            this.maxrow = this.row.length;
        }
        var mtr = 'mtr';
        var tag = this.factory.configuration.tags.getTag();
        if (tag) {
            this.row = [tag].concat(this.row);
            mtr = 'mlabeledtr';
        }
        this.factory.configuration.tags.clearTag();
        var node = this.create('node', mtr, this.row);
        this.table.push(node);
        this.row = [];
    };
    EqnArrayItem.prototype.EndTable = function () {
        _super.prototype.EndTable.call(this);
        this.factory.configuration.tags.end();
        this.extendArray('columnalign', this.maxrow);
        this.extendArray('columnwidth', this.maxrow);
        this.extendArray('columnspacing', this.maxrow - 1);
    };
    EqnArrayItem.prototype.extendArray = function (name, max) {
        if (!this.arraydef[name])
            return;
        var repeat = this.arraydef[name].split(/ /);
        var columns = __spreadArray([], __read(repeat), false);
        if (columns.length > 1) {
            while (columns.length < max) {
                columns.push.apply(columns, __spreadArray([], __read(repeat), false));
            }
            this.arraydef[name] = columns.slice(0, max).join(' ');
        }
    };
    return EqnArrayItem;
}(ArrayItem));
exports.EqnArrayItem = EqnArrayItem;
var EquationItem = (function (_super) {
    __extends(EquationItem, _super);
    function EquationItem(factory) {
        var args = [];
        for (var _i = 1; _i < arguments.length; _i++) {
            args[_i - 1] = arguments[_i];
        }
        var _this = _super.call(this, factory) || this;
        _this.factory.configuration.tags.start('equation', true, args[0]);
        return _this;
    }
    Object.defineProperty(EquationItem.prototype, "kind", {
        get: function () {
            return 'equation';
        },
        enumerable: false,
        configurable: true
    });
    Object.defineProperty(EquationItem.prototype, "isOpen", {
        get: function () {
            return true;
        },
        enumerable: false,
        configurable: true
    });
    EquationItem.prototype.checkItem = function (item) {
        if (item.isKind('end')) {
            var mml = this.toMml();
            var tag = this.factory.configuration.tags.getTag();
            this.factory.configuration.tags.end();
            return [[tag ? this.factory.configuration.tags.enTag(mml, tag) : mml, item], true];
        }
        if (item.isKind('stop')) {
            throw new TexError_js_1.default('EnvMissingEnd', 'Missing \\end{%1}', this.getName());
        }
        return _super.prototype.checkItem.call(this, item);
    };
    return EquationItem;
}(StackItem_js_1.BaseItem));
exports.EquationItem = EquationItem;
//# sourceMappingURL=BaseItems.js.map

/***/ }),

/***/ 23430:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    var desc = Object.getOwnPropertyDescriptor(m, k);
    if (!desc || ("get" in desc ? !m.__esModule : desc.writable || desc.configurable)) {
      desc = { enumerable: true, get: function() { return m[k]; } };
    }
    Object.defineProperty(o, k2, desc);
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
    Object.defineProperty(o, "default", { enumerable: true, value: v });
}) : function(o, v) {
    o["default"] = v;
});
var __importStar = (this && this.__importStar) || function (mod) {
    if (mod && mod.__esModule) return mod;
    var result = {};
    if (mod != null) for (var k in mod) if (k !== "default" && Object.prototype.hasOwnProperty.call(mod, k)) __createBinding(result, mod, k);
    __setModuleDefault(result, mod);
    return result;
};
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
var sm = __importStar(__webpack_require__(65695));
var TexConstants_js_1 = __webpack_require__(98927);
var BaseMethods_js_1 = __importDefault(__webpack_require__(76914));
var ParseMethods_js_1 = __importDefault(__webpack_require__(2362));
var ParseUtil_js_1 = __importDefault(__webpack_require__(55038));
var MmlNode_js_1 = __webpack_require__(83045);
var lengths_js_1 = __webpack_require__(56780);
new sm.RegExpMap('letter', ParseMethods_js_1.default.variable, /[a-z]/i);
new sm.RegExpMap('digit', ParseMethods_js_1.default.digit, /[0-9.,]/);
new sm.RegExpMap('command', ParseMethods_js_1.default.controlSequence, /^\\/);
new sm.MacroMap('special', {
    '{': 'Open',
    '}': 'Close',
    '~': 'Tilde',
    '^': 'Superscript',
    '_': 'Subscript',
    ' ': 'Space',
    '\t': 'Space',
    '\r': 'Space',
    '\n': 'Space',
    '\'': 'Prime',
    '%': 'Comment',
    '&': 'Entry',
    '#': 'Hash',
    '\u00A0': 'Space',
    '\u2019': 'Prime'
}, BaseMethods_js_1.default);
new sm.CharacterMap('mathchar0mi', ParseMethods_js_1.default.mathchar0mi, {
    alpha: '\u03B1',
    beta: '\u03B2',
    gamma: '\u03B3',
    delta: '\u03B4',
    epsilon: '\u03F5',
    zeta: '\u03B6',
    eta: '\u03B7',
    theta: '\u03B8',
    iota: '\u03B9',
    kappa: '\u03BA',
    lambda: '\u03BB',
    mu: '\u03BC',
    nu: '\u03BD',
    xi: '\u03BE',
    omicron: '\u03BF',
    pi: '\u03C0',
    rho: '\u03C1',
    sigma: '\u03C3',
    tau: '\u03C4',
    upsilon: '\u03C5',
    phi: '\u03D5',
    chi: '\u03C7',
    psi: '\u03C8',
    omega: '\u03C9',
    varepsilon: '\u03B5',
    vartheta: '\u03D1',
    varpi: '\u03D6',
    varrho: '\u03F1',
    varsigma: '\u03C2',
    varphi: '\u03C6',
    S: ['\u00A7', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    aleph: ['\u2135', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    hbar: ['\u210F', { variantForm: true }],
    imath: '\u0131',
    jmath: '\u0237',
    ell: '\u2113',
    wp: ['\u2118', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    Re: ['\u211C', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    Im: ['\u2111', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    partial: ['\u2202', { mathvariant: TexConstants_js_1.TexConstant.Variant.ITALIC }],
    infty: ['\u221E', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    prime: ['\u2032', { variantForm: true }],
    emptyset: ['\u2205', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    nabla: ['\u2207', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    top: ['\u22A4', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    bot: ['\u22A5', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    angle: ['\u2220', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    triangle: ['\u25B3', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    backslash: ['\u2216', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    forall: ['\u2200', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    exists: ['\u2203', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    neg: ['\u00AC', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    lnot: ['\u00AC', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    flat: ['\u266D', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    natural: ['\u266E', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    sharp: ['\u266F', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    clubsuit: ['\u2663', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    diamondsuit: ['\u2662', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    heartsuit: ['\u2661', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }],
    spadesuit: ['\u2660', { mathvariant: TexConstants_js_1.TexConstant.Variant.NORMAL }]
});
new sm.CharacterMap('mathchar0mo', ParseMethods_js_1.default.mathchar0mo, {
    surd: '\u221A',
    coprod: ['\u2210', { texClass: MmlNode_js_1.TEXCLASS.OP,
            movesupsub: true }],
    bigvee: ['\u22C1', { texClass: MmlNode_js_1.TEXCLASS.OP,
            movesupsub: true }],
    bigwedge: ['\u22C0', { texClass: MmlNode_js_1.TEXCLASS.OP,
            movesupsub: true }],
    biguplus: ['\u2A04', { texClass: MmlNode_js_1.TEXCLASS.OP,
            movesupsub: true }],
    bigcap: ['\u22C2', { texClass: MmlNode_js_1.TEXCLASS.OP,
            movesupsub: true }],
    bigcup: ['\u22C3', { texClass: MmlNode_js_1.TEXCLASS.OP,
            movesupsub: true }],
    'int': ['\u222B', { texClass: MmlNode_js_1.TEXCLASS.OP }],
    intop: ['\u222B', { texClass: MmlNode_js_1.TEXCLASS.OP,
            movesupsub: true, movablelimits: true }],
    iint: ['\u222C', { texClass: MmlNode_js_1.TEXCLASS.OP }],
    iiint: ['\u222D', { texClass: MmlNode_js_1.TEXCLASS.OP }],
    prod: ['\u220F', { texClass: MmlNode_js_1.TEXCLASS.OP,
            movesupsub: true }],
    sum: ['\u2211', { texClass: MmlNode_js_1.TEXCLASS.OP,
            movesupsub: true }],
    bigotimes: ['\u2A02', { texClass: MmlNode_js_1.TEXCLASS.OP,
            movesupsub: true }],
    bigoplus: ['\u2A01', { texClass: MmlNode_js_1.TEXCLASS.OP,
            movesupsub: true }],
    bigodot: ['\u2A00', { texClass: MmlNode_js_1.TEXCLASS.OP,
            movesupsub: true }],
    oint: ['\u222E', { texClass: MmlNode_js_1.TEXCLASS.OP }],
    bigsqcup: ['\u2A06', { texClass: MmlNode_js_1.TEXCLASS.OP,
            movesupsub: true }],
    smallint: ['\u222B', { largeop: false }],
    triangleleft: '\u25C3',
    triangleright: '\u25B9',
    bigtriangleup: '\u25B3',
    bigtriangledown: '\u25BD',
    wedge: '\u2227',
    land: '\u2227',
    vee: '\u2228',
    lor: '\u2228',
    cap: '\u2229',
    cup: '\u222A',
    ddagger: '\u2021',
    dagger: '\u2020',
    sqcap: '\u2293',
    sqcup: '\u2294',
    uplus: '\u228E',
    amalg: '\u2A3F',
    diamond: '\u22C4',
    bullet: '\u2219',
    wr: '\u2240',
    div: '\u00F7',
    divsymbol: '\u00F7',
    odot: ['\u2299', { largeop: false }],
    oslash: ['\u2298', { largeop: false }],
    otimes: ['\u2297', { largeop: false }],
    ominus: ['\u2296', { largeop: false }],
    oplus: ['\u2295', { largeop: false }],
    mp: '\u2213',
    pm: '\u00B1',
    circ: '\u2218',
    bigcirc: '\u25EF',
    setminus: '\u2216',
    cdot: '\u22C5',
    ast: '\u2217',
    times: '\u00D7',
    star: '\u22C6',
    propto: '\u221D',
    sqsubseteq: '\u2291',
    sqsupseteq: '\u2292',
    parallel: '\u2225',
    mid: '\u2223',
    dashv: '\u22A3',
    vdash: '\u22A2',
    leq: '\u2264',
    le: '\u2264',
    geq: '\u2265',
    ge: '\u2265',
    lt: '\u003C',
    gt: '\u003E',
    succ: '\u227B',
    prec: '\u227A',
    approx: '\u2248',
    succeq: '\u2AB0',
    preceq: '\u2AAF',
    supset: '\u2283',
    subset: '\u2282',
    supseteq: '\u2287',
    subseteq: '\u2286',
    'in': '\u2208',
    ni: '\u220B',
    notin: '\u2209',
    owns: '\u220B',
    gg: '\u226B',
    ll: '\u226A',
    sim: '\u223C',
    simeq: '\u2243',
    perp: '\u22A5',
    equiv: '\u2261',
    asymp: '\u224D',
    smile: '\u2323',
    frown: '\u2322',
    ne: '\u2260',
    neq: '\u2260',
    cong: '\u2245',
    doteq: '\u2250',
    bowtie: '\u22C8',
    models: '\u22A8',
    notChar: '\u29F8',
    Leftrightarrow: '\u21D4',
    Leftarrow: '\u21D0',
    Rightarrow: '\u21D2',
    leftrightarrow: '\u2194',
    leftarrow: '\u2190',
    gets: '\u2190',
    rightarrow: '\u2192',
    to: ['\u2192', { accent: false }],
    mapsto: '\u21A6',
    leftharpoonup: '\u21BC',
    leftharpoondown: '\u21BD',
    rightharpoonup: '\u21C0',
    rightharpoondown: '\u21C1',
    nearrow: '\u2197',
    searrow: '\u2198',
    nwarrow: '\u2196',
    swarrow: '\u2199',
    rightleftharpoons: '\u21CC',
    hookrightarrow: '\u21AA',
    hookleftarrow: '\u21A9',
    longleftarrow: '\u27F5',
    Longleftarrow: '\u27F8',
    longrightarrow: '\u27F6',
    Longrightarrow: '\u27F9',
    Longleftrightarrow: '\u27FA',
    longleftrightarrow: '\u27F7',
    longmapsto: '\u27FC',
    ldots: '\u2026',
    cdots: '\u22EF',
    vdots: '\u22EE',
    ddots: '\u22F1',
    dotsc: '\u2026',
    dotsb: '\u22EF',
    dotsm: '\u22EF',
    dotsi: '\u22EF',
    dotso: '\u2026',
    ldotp: ['\u002E', { texClass: MmlNode_js_1.TEXCLASS.PUNCT }],
    cdotp: ['\u22C5', { texClass: MmlNode_js_1.TEXCLASS.PUNCT }],
    colon: ['\u003A', { texClass: MmlNode_js_1.TEXCLASS.PUNCT }]
});
new sm.CharacterMap('mathchar7', ParseMethods_js_1.default.mathchar7, {
    Gamma: '\u0393',
    Delta: '\u0394',
    Theta: '\u0398',
    Lambda: '\u039B',
    Xi: '\u039E',
    Pi: '\u03A0',
    Sigma: '\u03A3',
    Upsilon: '\u03A5',
    Phi: '\u03A6',
    Psi: '\u03A8',
    Omega: '\u03A9',
    '_': '\u005F',
    '#': '\u0023',
    '$': '\u0024',
    '%': '\u0025',
    '&': '\u0026',
    And: '\u0026'
});
new sm.DelimiterMap('delimiter', ParseMethods_js_1.default.delimiter, {
    '(': '(',
    ')': ')',
    '[': '[',
    ']': ']',
    '<': '\u27E8',
    '>': '\u27E9',
    '\\lt': '\u27E8',
    '\\gt': '\u27E9',
    '/': '/',
    '|': ['|', { texClass: MmlNode_js_1.TEXCLASS.ORD }],
    '.': '',
    '\\\\': '\\',
    '\\lmoustache': '\u23B0',
    '\\rmoustache': '\u23B1',
    '\\lgroup': '\u27EE',
    '\\rgroup': '\u27EF',
    '\\arrowvert': '\u23D0',
    '\\Arrowvert': '\u2016',
    '\\bracevert': '\u23AA',
    '\\Vert': ['\u2016', { texClass: MmlNode_js_1.TEXCLASS.ORD }],
    '\\|': ['\u2016', { texClass: MmlNode_js_1.TEXCLASS.ORD }],
    '\\vert': ['|', { texClass: MmlNode_js_1.TEXCLASS.ORD }],
    '\\uparrow': '\u2191',
    '\\downarrow': '\u2193',
    '\\updownarrow': '\u2195',
    '\\Uparrow': '\u21D1',
    '\\Downarrow': '\u21D3',
    '\\Updownarrow': '\u21D5',
    '\\backslash': '\\',
    '\\rangle': '\u27E9',
    '\\langle': '\u27E8',
    '\\rbrace': '}',
    '\\lbrace': '{',
    '\\}': '}',
    '\\{': '{',
    '\\rceil': '\u2309',
    '\\lceil': '\u2308',
    '\\rfloor': '\u230B',
    '\\lfloor': '\u230A',
    '\\lbrack': '[',
    '\\rbrack': ']'
});
new sm.CommandMap('macros', {
    displaystyle: ['SetStyle', 'D', true, 0],
    textstyle: ['SetStyle', 'T', false, 0],
    scriptstyle: ['SetStyle', 'S', false, 1],
    scriptscriptstyle: ['SetStyle', 'SS', false, 2],
    rm: ['SetFont', TexConstants_js_1.TexConstant.Variant.NORMAL],
    mit: ['SetFont', TexConstants_js_1.TexConstant.Variant.ITALIC],
    oldstyle: ['SetFont', TexConstants_js_1.TexConstant.Variant.OLDSTYLE],
    cal: ['SetFont', TexConstants_js_1.TexConstant.Variant.CALLIGRAPHIC],
    it: ['SetFont', TexConstants_js_1.TexConstant.Variant.MATHITALIC],
    bf: ['SetFont', TexConstants_js_1.TexConstant.Variant.BOLD],
    bbFont: ['SetFont', TexConstants_js_1.TexConstant.Variant.DOUBLESTRUCK],
    scr: ['SetFont', TexConstants_js_1.TexConstant.Variant.SCRIPT],
    frak: ['SetFont', TexConstants_js_1.TexConstant.Variant.FRAKTUR],
    sf: ['SetFont', TexConstants_js_1.TexConstant.Variant.SANSSERIF],
    tt: ['SetFont', TexConstants_js_1.TexConstant.Variant.MONOSPACE],
    mathrm: ['MathFont', TexConstants_js_1.TexConstant.Variant.NORMAL],
    mathup: ['MathFont', TexConstants_js_1.TexConstant.Variant.NORMAL],
    mathnormal: ['MathFont', ''],
    mathbf: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLD],
    mathbfup: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLD],
    mathit: ['MathFont', TexConstants_js_1.TexConstant.Variant.MATHITALIC],
    mathbfit: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLDITALIC],
    mathbb: ['MathFont', TexConstants_js_1.TexConstant.Variant.DOUBLESTRUCK],
    Bbb: ['MathFont', TexConstants_js_1.TexConstant.Variant.DOUBLESTRUCK],
    mathfrak: ['MathFont', TexConstants_js_1.TexConstant.Variant.FRAKTUR],
    mathbffrak: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLDFRAKTUR],
    mathscr: ['MathFont', TexConstants_js_1.TexConstant.Variant.SCRIPT],
    mathbfscr: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLDSCRIPT],
    mathsf: ['MathFont', TexConstants_js_1.TexConstant.Variant.SANSSERIF],
    mathsfup: ['MathFont', TexConstants_js_1.TexConstant.Variant.SANSSERIF],
    mathbfsf: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLDSANSSERIF],
    mathbfsfup: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLDSANSSERIF],
    mathsfit: ['MathFont', TexConstants_js_1.TexConstant.Variant.SANSSERIFITALIC],
    mathbfsfit: ['MathFont', TexConstants_js_1.TexConstant.Variant.SANSSERIFBOLDITALIC],
    mathtt: ['MathFont', TexConstants_js_1.TexConstant.Variant.MONOSPACE],
    mathcal: ['MathFont', TexConstants_js_1.TexConstant.Variant.CALLIGRAPHIC],
    mathbfcal: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLDCALLIGRAPHIC],
    symrm: ['MathFont', TexConstants_js_1.TexConstant.Variant.NORMAL],
    symup: ['MathFont', TexConstants_js_1.TexConstant.Variant.NORMAL],
    symnormal: ['MathFont', ''],
    symbf: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLD],
    symbfup: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLD],
    symit: ['MathFont', TexConstants_js_1.TexConstant.Variant.ITALIC],
    symbfit: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLDITALIC],
    symbb: ['MathFont', TexConstants_js_1.TexConstant.Variant.DOUBLESTRUCK],
    symfrak: ['MathFont', TexConstants_js_1.TexConstant.Variant.FRAKTUR],
    symbffrak: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLDFRAKTUR],
    symscr: ['MathFont', TexConstants_js_1.TexConstant.Variant.SCRIPT],
    symbfscr: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLDSCRIPT],
    symsf: ['MathFont', TexConstants_js_1.TexConstant.Variant.SANSSERIF],
    symsfup: ['MathFont', TexConstants_js_1.TexConstant.Variant.SANSSERIF],
    symbfsf: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLDSANSSERIF],
    symbfsfup: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLDSANSSERIF],
    symsfit: ['MathFont', TexConstants_js_1.TexConstant.Variant.SANSSERIFITALIC],
    symbfsfit: ['MathFont', TexConstants_js_1.TexConstant.Variant.SANSSERIFBOLDITALIC],
    symtt: ['MathFont', TexConstants_js_1.TexConstant.Variant.MONOSPACE],
    symcal: ['MathFont', TexConstants_js_1.TexConstant.Variant.CALLIGRAPHIC],
    symbfcal: ['MathFont', TexConstants_js_1.TexConstant.Variant.BOLDCALLIGRAPHIC],
    textrm: ['HBox', null, TexConstants_js_1.TexConstant.Variant.NORMAL],
    textup: ['HBox', null, TexConstants_js_1.TexConstant.Variant.NORMAL],
    textnormal: ['HBox'],
    textit: ['HBox', null, TexConstants_js_1.TexConstant.Variant.ITALIC],
    textbf: ['HBox', null, TexConstants_js_1.TexConstant.Variant.BOLD],
    textsf: ['HBox', null, TexConstants_js_1.TexConstant.Variant.SANSSERIF],
    texttt: ['HBox', null, TexConstants_js_1.TexConstant.Variant.MONOSPACE],
    tiny: ['SetSize', 0.5],
    Tiny: ['SetSize', 0.6],
    scriptsize: ['SetSize', 0.7],
    small: ['SetSize', 0.85],
    normalsize: ['SetSize', 1.0],
    large: ['SetSize', 1.2],
    Large: ['SetSize', 1.44],
    LARGE: ['SetSize', 1.73],
    huge: ['SetSize', 2.07],
    Huge: ['SetSize', 2.49],
    arcsin: 'NamedFn',
    arccos: 'NamedFn',
    arctan: 'NamedFn',
    arg: 'NamedFn',
    cos: 'NamedFn',
    cosh: 'NamedFn',
    cot: 'NamedFn',
    coth: 'NamedFn',
    csc: 'NamedFn',
    deg: 'NamedFn',
    det: 'NamedOp',
    dim: 'NamedFn',
    exp: 'NamedFn',
    gcd: 'NamedOp',
    hom: 'NamedFn',
    inf: 'NamedOp',
    ker: 'NamedFn',
    lg: 'NamedFn',
    lim: 'NamedOp',
    liminf: ['NamedOp', 'lim&thinsp;inf'],
    limsup: ['NamedOp', 'lim&thinsp;sup'],
    ln: 'NamedFn',
    log: 'NamedFn',
    max: 'NamedOp',
    min: 'NamedOp',
    Pr: 'NamedOp',
    sec: 'NamedFn',
    sin: 'NamedFn',
    sinh: 'NamedFn',
    sup: 'NamedOp',
    tan: 'NamedFn',
    tanh: 'NamedFn',
    limits: ['Limits', 1],
    nolimits: ['Limits', 0],
    overline: ['UnderOver', '2015'],
    underline: ['UnderOver', '2015'],
    overbrace: ['UnderOver', '23DE', 1],
    underbrace: ['UnderOver', '23DF', 1],
    overparen: ['UnderOver', '23DC'],
    underparen: ['UnderOver', '23DD'],
    overrightarrow: ['UnderOver', '2192'],
    underrightarrow: ['UnderOver', '2192'],
    overleftarrow: ['UnderOver', '2190'],
    underleftarrow: ['UnderOver', '2190'],
    overleftrightarrow: ['UnderOver', '2194'],
    underleftrightarrow: ['UnderOver', '2194'],
    overset: 'Overset',
    underset: 'Underset',
    overunderset: 'Overunderset',
    stackrel: ['Macro', '\\mathrel{\\mathop{#2}\\limits^{#1}}', 2],
    stackbin: ['Macro', '\\mathbin{\\mathop{#2}\\limits^{#1}}', 2],
    over: 'Over',
    overwithdelims: 'Over',
    atop: 'Over',
    atopwithdelims: 'Over',
    above: 'Over',
    abovewithdelims: 'Over',
    brace: ['Over', '{', '}'],
    brack: ['Over', '[', ']'],
    choose: ['Over', '(', ')'],
    frac: 'Frac',
    sqrt: 'Sqrt',
    root: 'Root',
    uproot: ['MoveRoot', 'upRoot'],
    leftroot: ['MoveRoot', 'leftRoot'],
    left: 'LeftRight',
    right: 'LeftRight',
    middle: 'LeftRight',
    llap: 'Lap',
    rlap: 'Lap',
    raise: 'RaiseLower',
    lower: 'RaiseLower',
    moveleft: 'MoveLeftRight',
    moveright: 'MoveLeftRight',
    ',': ['Spacer', lengths_js_1.MATHSPACE.thinmathspace],
    ':': ['Spacer', lengths_js_1.MATHSPACE.mediummathspace],
    '>': ['Spacer', lengths_js_1.MATHSPACE.mediummathspace],
    ';': ['Spacer', lengths_js_1.MATHSPACE.thickmathspace],
    '!': ['Spacer', lengths_js_1.MATHSPACE.negativethinmathspace],
    enspace: ['Spacer', .5],
    quad: ['Spacer', 1],
    qquad: ['Spacer', 2],
    thinspace: ['Spacer', lengths_js_1.MATHSPACE.thinmathspace],
    negthinspace: ['Spacer', lengths_js_1.MATHSPACE.negativethinmathspace],
    hskip: 'Hskip',
    hspace: 'Hskip',
    kern: 'Hskip',
    mskip: 'Hskip',
    mspace: 'Hskip',
    mkern: 'Hskip',
    rule: 'rule',
    Rule: ['Rule'],
    Space: ['Rule', 'blank'],
    nonscript: 'Nonscript',
    big: ['MakeBig', MmlNode_js_1.TEXCLASS.ORD, 0.85],
    Big: ['MakeBig', MmlNode_js_1.TEXCLASS.ORD, 1.15],
    bigg: ['MakeBig', MmlNode_js_1.TEXCLASS.ORD, 1.45],
    Bigg: ['MakeBig', MmlNode_js_1.TEXCLASS.ORD, 1.75],
    bigl: ['MakeBig', MmlNode_js_1.TEXCLASS.OPEN, 0.85],
    Bigl: ['MakeBig', MmlNode_js_1.TEXCLASS.OPEN, 1.15],
    biggl: ['MakeBig', MmlNode_js_1.TEXCLASS.OPEN, 1.45],
    Biggl: ['MakeBig', MmlNode_js_1.TEXCLASS.OPEN, 1.75],
    bigr: ['MakeBig', MmlNode_js_1.TEXCLASS.CLOSE, 0.85],
    Bigr: ['MakeBig', MmlNode_js_1.TEXCLASS.CLOSE, 1.15],
    biggr: ['MakeBig', MmlNode_js_1.TEXCLASS.CLOSE, 1.45],
    Biggr: ['MakeBig', MmlNode_js_1.TEXCLASS.CLOSE, 1.75],
    bigm: ['MakeBig', MmlNode_js_1.TEXCLASS.REL, 0.85],
    Bigm: ['MakeBig', MmlNode_js_1.TEXCLASS.REL, 1.15],
    biggm: ['MakeBig', MmlNode_js_1.TEXCLASS.REL, 1.45],
    Biggm: ['MakeBig', MmlNode_js_1.TEXCLASS.REL, 1.75],
    mathord: ['TeXAtom', MmlNode_js_1.TEXCLASS.ORD],
    mathop: ['TeXAtom', MmlNode_js_1.TEXCLASS.OP],
    mathopen: ['TeXAtom', MmlNode_js_1.TEXCLASS.OPEN],
    mathclose: ['TeXAtom', MmlNode_js_1.TEXCLASS.CLOSE],
    mathbin: ['TeXAtom', MmlNode_js_1.TEXCLASS.BIN],
    mathrel: ['TeXAtom', MmlNode_js_1.TEXCLASS.REL],
    mathpunct: ['TeXAtom', MmlNode_js_1.TEXCLASS.PUNCT],
    mathinner: ['TeXAtom', MmlNode_js_1.TEXCLASS.INNER],
    vcenter: ['TeXAtom', MmlNode_js_1.TEXCLASS.VCENTER],
    buildrel: 'BuildRel',
    hbox: ['HBox', 0],
    text: 'HBox',
    mbox: ['HBox', 0],
    fbox: 'FBox',
    boxed: ['Macro', '\\fbox{$\\displaystyle{#1}$}', 1],
    framebox: 'FrameBox',
    strut: 'Strut',
    mathstrut: ['Macro', '\\vphantom{(}'],
    phantom: 'Phantom',
    vphantom: ['Phantom', 1, 0],
    hphantom: ['Phantom', 0, 1],
    smash: 'Smash',
    acute: ['Accent', '00B4'],
    grave: ['Accent', '0060'],
    ddot: ['Accent', '00A8'],
    tilde: ['Accent', '007E'],
    bar: ['Accent', '00AF'],
    breve: ['Accent', '02D8'],
    check: ['Accent', '02C7'],
    hat: ['Accent', '005E'],
    vec: ['Accent', '2192'],
    dot: ['Accent', '02D9'],
    widetilde: ['Accent', '007E', 1],
    widehat: ['Accent', '005E', 1],
    matrix: 'Matrix',
    array: 'Matrix',
    pmatrix: ['Matrix', '(', ')'],
    cases: ['Matrix', '{', '', 'left left', null, '.1em', null,
        true],
    eqalign: ['Matrix', null, null, 'right left',
        (0, lengths_js_1.em)(lengths_js_1.MATHSPACE.thickmathspace), '.5em', 'D'],
    displaylines: ['Matrix', null, null, 'center', null, '.5em', 'D'],
    cr: 'Cr',
    '\\': 'CrLaTeX',
    newline: ['CrLaTeX', true],
    hline: ['HLine', 'solid'],
    hdashline: ['HLine', 'dashed'],
    eqalignno: ['Matrix', null, null, 'right left',
        (0, lengths_js_1.em)(lengths_js_1.MATHSPACE.thickmathspace), '.5em', 'D', null,
        'right'],
    leqalignno: ['Matrix', null, null, 'right left',
        (0, lengths_js_1.em)(lengths_js_1.MATHSPACE.thickmathspace), '.5em', 'D', null,
        'left'],
    hfill: 'HFill',
    hfil: 'HFill',
    hfilll: 'HFill',
    bmod: ['Macro', '\\mmlToken{mo}[lspace="thickmathspace"' +
            ' rspace="thickmathspace"]{mod}'],
    pmod: ['Macro', '\\pod{\\mmlToken{mi}{mod}\\kern 6mu #1}', 1],
    mod: ['Macro', '\\mathchoice{\\kern18mu}{\\kern12mu}' +
            '{\\kern12mu}{\\kern12mu}\\mmlToken{mi}{mod}\\,\\,#1',
        1],
    pod: ['Macro', '\\mathchoice{\\kern18mu}{\\kern8mu}' +
            '{\\kern8mu}{\\kern8mu}(#1)', 1],
    iff: ['Macro', '\\;\\Longleftrightarrow\\;'],
    skew: ['Macro', '{{#2{#3\\mkern#1mu}\\mkern-#1mu}{}}', 3],
    pmb: ['Macro', '\\rlap{#1}\\kern1px{#1}', 1],
    TeX: ['Macro', 'T\\kern-.14em\\lower.5ex{E}\\kern-.115em X'],
    LaTeX: ['Macro', 'L\\kern-.325em\\raise.21em' +
            '{\\scriptstyle{A}}\\kern-.17em\\TeX'],
    ' ': ['Macro', '\\text{ }'],
    not: 'Not',
    dots: 'Dots',
    space: 'Tilde',
    '\u00A0': 'Tilde',
    begin: 'BeginEnd',
    end: 'BeginEnd',
    label: 'HandleLabel',
    ref: 'HandleRef',
    nonumber: 'HandleNoTag',
    mathchoice: 'MathChoice',
    mmlToken: 'MmlToken'
}, BaseMethods_js_1.default);
new sm.EnvironmentMap('environment', ParseMethods_js_1.default.environment, {
    array: ['AlignedArray'],
    equation: ['Equation', null, true],
    eqnarray: ['EqnArray', null, true, true, 'rcl',
        ParseUtil_js_1.default.cols(0, lengths_js_1.MATHSPACE.thickmathspace), '.5em']
}, BaseMethods_js_1.default);
new sm.CharacterMap('not_remap', null, {
    '\u2190': '\u219A',
    '\u2192': '\u219B',
    '\u2194': '\u21AE',
    '\u21D0': '\u21CD',
    '\u21D2': '\u21CF',
    '\u21D4': '\u21CE',
    '\u2208': '\u2209',
    '\u220B': '\u220C',
    '\u2223': '\u2224',
    '\u2225': '\u2226',
    '\u223C': '\u2241',
    '\u007E': '\u2241',
    '\u2243': '\u2244',
    '\u2245': '\u2247',
    '\u2248': '\u2249',
    '\u224D': '\u226D',
    '\u003D': '\u2260',
    '\u2261': '\u2262',
    '\u003C': '\u226E',
    '\u003E': '\u226F',
    '\u2264': '\u2270',
    '\u2265': '\u2271',
    '\u2272': '\u2274',
    '\u2273': '\u2275',
    '\u2276': '\u2278',
    '\u2277': '\u2279',
    '\u227A': '\u2280',
    '\u227B': '\u2281',
    '\u2282': '\u2284',
    '\u2283': '\u2285',
    '\u2286': '\u2288',
    '\u2287': '\u2289',
    '\u22A2': '\u22AC',
    '\u22A8': '\u22AD',
    '\u22A9': '\u22AE',
    '\u22AB': '\u22AF',
    '\u227C': '\u22E0',
    '\u227D': '\u22E1',
    '\u2291': '\u22E2',
    '\u2292': '\u22E3',
    '\u22B2': '\u22EA',
    '\u22B3': '\u22EB',
    '\u22B4': '\u22EC',
    '\u22B5': '\u22ED',
    '\u2203': '\u2204'
});
//# sourceMappingURL=BaseMappings.js.map

/***/ }),

/***/ 76914:
/***/ (function(__unused_webpack_module, exports, __webpack_require__) {


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
var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    var desc = Object.getOwnPropertyDescriptor(m, k);
    if (!desc || ("get" in desc ? !m.__esModule : desc.writable || desc.configurable)) {
      desc = { enumerable: true, get: function() { return m[k]; } };
    }
    Object.defineProperty(o, k2, desc);
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
    Object.defineProperty(o, "default", { enumerable: true, value: v });
}) : function(o, v) {
    o["default"] = v;
});
var __importStar = (this && this.__importStar) || function (mod) {
    if (mod && mod.__esModule) return mod;
    var result = {};
    if (mod != null) for (var k in mod) if (k !== "default" && Object.prototype.hasOwnProperty.call(mod, k)) __createBinding(result, mod, k);
    __setModuleDefault(result, mod);
    return result;
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
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", ({ value: true }));
var sitem = __importStar(__webpack_require__(31201));
var NodeUtil_js_1 = __importDefault(__webpack_require__(53972));
var TexError_js_1 = __importDefault(__webpack_require__(54420));
var TexParser_js_1 = __importDefault(__webpack_require__(94032));
var TexConstants_js_1 = __webpack_require__(98927);
var ParseUtil_js_1 = __importDefault(__webpack_require__(55038));
var MmlNode_js_1 = __webpack_require__(83045);
var Tags_js_1 = __webpack_require__(75723);
var lengths_js_1 = __webpack_require__(56780);
var Entities_js_1 = __webpack_require__(61051);
var Options_js_1 = __webpack_require__(4498);
var BaseMethods = {};
var P_HEIGHT = 1.2 / .85;
var MmlTokenAllow = {
    fontfamily: 1, fontsize: 1, fontweight: 1, fontstyle: 1,
    color: 1, background: 1,
    id: 1, 'class': 1, href: 1, style: 1
};
BaseMethods.Open = function (parser, _c) {
    parser.Push(parser.itemFactory.create('open'));
};
BaseMethods.Close = function (parser, _c) {
    parser.Push(parser.itemFactory.create('close'));
};
BaseMethods.Tilde = function (parser, _c) {
    parser.Push(parser.create('token', 'mtext', {}, Entities_js_1.entities.nbsp));
};
BaseMethods.Space = function (_parser, _c) { };
BaseMethods.Superscript = function (parser, _c) {
    var _a;
    if (parser.GetNext().match(/\d/)) {
        parser.string = parser.string.substr(0, parser.i + 1) +
            ' ' + parser.string.substr(parser.i + 1);
    }
    var primes;
    var base;
    var top = parser.stack.Top();
    if (top.isKind('prime')) {
        _a = __read(top.Peek(2), 2), base = _a[0], primes = _a[1];
        parser.stack.Pop();
    }
    else {
        base = parser.stack.Prev();
        if (!base) {
            base = parser.create('token', 'mi', {}, '');
        }
    }
    var movesupsub = NodeUtil_js_1.default.getProperty(base, 'movesupsub');
    var position = NodeUtil_js_1.default.isType(base, 'msubsup') ? base.sup :
        base.over;
    if ((NodeUtil_js_1.default.isType(base, 'msubsup') && !NodeUtil_js_1.default.isType(base, 'msup') &&
        NodeUtil_js_1.default.getChildAt(base, base.sup)) ||
        (NodeUtil_js_1.default.isType(base, 'munderover') && !NodeUtil_js_1.default.isType(base, 'mover') &&
            NodeUtil_js_1.default.getChildAt(base, base.over) &&
            !NodeUtil_js_1.default.getProperty(base, 'subsupOK'))) {
        throw new TexError_js_1.default('DoubleExponent', 'Double exponent: use braces to clarify');
    }
    if (!NodeUtil_js_1.default.isType(base, 'msubsup') || NodeUtil_js_1.default.isType(base, 'msup')) {
        if (movesupsub) {
            if (!NodeUtil_js_1.default.isType(base, 'munderover') || NodeUtil_js_1.default.isType(base, 'mover') ||
                NodeUtil_js_1.default.getChildAt(base, base.over)) {
                base = parser.create('node', 'munderover', [base], { movesupsub: true });
            }
            position = base.over;
        }
        else {
            base = parser.create('node', 'msubsup', [base]);
            position = base.sup;
        }
    }
    parser.Push(parser.itemFactory.create('subsup', base).setProperties({
        position: position, primes: primes, movesupsub: movesupsub
    }));
};
BaseMethods.Subscript = function (parser, _c) {
    var _a;
    if (parser.GetNext().match(/\d/)) {
        parser.string =
            parser.string.substr(0, parser.i + 1) + ' ' +
                parser.string.substr(parser.i + 1);
    }
    var primes, base;
    var top = parser.stack.Top();
    if (top.isKind('prime')) {
        _a = __read(top.Peek(2), 2), base = _a[0], primes = _a[1];
        parser.stack.Pop();
    }
    else {
        base = parser.stack.Prev();
        if (!base) {
            base = parser.create('token', 'mi', {}, '');
        }
    }
    var movesupsub = NodeUtil_js_1.default.getProperty(base, 'movesupsub');
    var position = NodeUtil_js_1.default.isType(base, 'msubsup') ?
        base.sub : base.under;
    if ((NodeUtil_js_1.default.isType(base, 'msubsup') && !NodeUtil_js_1.default.isType(base, 'msup') &&
        NodeUtil_js_1.default.getChildAt(base, base.sub)) ||
        (NodeUtil_js_1.default.isType(base, 'munderover') && !NodeUtil_js_1.default.isType(base, 'mover') &&
            NodeUtil_js_1.default.getChildAt(base, base.under) &&
            !NodeUtil_js_1.default.getProperty(base, 'subsupOK'))) {
        throw new TexError_js_1.default('DoubleSubscripts', 'Double subscripts: use braces to clarify');
    }
    if (!NodeUtil_js_1.default.isType(base, 'msubsup') || NodeUtil_js_1.default.isType(base, 'msup')) {
        if (movesupsub) {
            if (!NodeUtil_js_1.default.isType(base, 'munderover') || NodeUtil_js_1.default.isType(base, 'mover') ||
                NodeUtil_js_1.default.getChildAt(base, base.under)) {
                base = parser.create('node', 'munderover', [base], { movesupsub: true });
            }
            position = base.under;
        }
        else {
            base = parser.create('node', 'msubsup', [base]);
            position = base.sub;
        }
    }
    parser.Push(parser.itemFactory.create('subsup', base).setProperties({
        position: position, primes: primes, movesupsub: movesupsub
    }));
};
BaseMethods.Prime = function (parser, c) {
    var base = parser.stack.Prev();
    if (!base) {
        base = parser.create('node', 'mi');
    }
    if (NodeUtil_js_1.default.isType(base, 'msubsup') && !NodeUtil_js_1.default.isType(base, 'msup') &&
        NodeUtil_js_1.default.getChildAt(base, base.sup)) {
        throw new TexError_js_1.default('DoubleExponentPrime', 'Prime causes double exponent: use braces to clarify');
    }
    var sup = '';
    parser.i--;
    do {
        sup += Entities_js_1.entities.prime;
        parser.i++, c = parser.GetNext();
    } while (c === '\'' || c === Entities_js_1.entities.rsquo);
    sup = ['', '\u2032', '\u2033', '\u2034', '\u2057'][sup.length] || sup;
    var node = parser.create('token', 'mo', { variantForm: true }, sup);
    parser.Push(parser.itemFactory.create('prime', base, node));
};
BaseMethods.Comment = function (parser, _c) {
    while (parser.i < parser.string.length && parser.string.charAt(parser.i) !== '\n') {
        parser.i++;
    }
};
BaseMethods.Hash = function (_parser, _c) {
    throw new TexError_js_1.default('CantUseHash1', 'You can\'t use \'macro parameter character #\' in math mode');
};
BaseMethods.MathFont = function (parser, name, variant) {
    var text = parser.GetArgument(name);
    var mml = new TexParser_js_1.default(text, __assign(__assign({}, parser.stack.env), { font: variant, multiLetterIdentifiers: /^[a-zA-Z]+/, noAutoOP: true }), parser.configuration).mml();
    parser.Push(parser.create('node', 'TeXAtom', [mml]));
};
BaseMethods.SetFont = function (parser, _name, font) {
    parser.stack.env['font'] = font;
};
BaseMethods.SetStyle = function (parser, _name, texStyle, style, level) {
    parser.stack.env['style'] = texStyle;
    parser.stack.env['level'] = level;
    parser.Push(parser.itemFactory.create('style').setProperty('styles', { displaystyle: style, scriptlevel: level }));
};
BaseMethods.SetSize = function (parser, _name, size) {
    parser.stack.env['size'] = size;
    parser.Push(parser.itemFactory.create('style').setProperty('styles', { mathsize: (0, lengths_js_1.em)(size) }));
};
BaseMethods.Spacer = function (parser, _name, space) {
    var node = parser.create('node', 'mspace', [], { width: (0, lengths_js_1.em)(space) });
    var style = parser.create('node', 'mstyle', [node], { scriptlevel: 0 });
    parser.Push(style);
};
BaseMethods.LeftRight = function (parser, name) {
    var first = name.substr(1);
    parser.Push(parser.itemFactory.create(first, parser.GetDelimiter(name), parser.stack.env.color));
};
BaseMethods.NamedFn = function (parser, name, id) {
    if (!id) {
        id = name.substr(1);
    }
    var mml = parser.create('token', 'mi', { texClass: MmlNode_js_1.TEXCLASS.OP }, id);
    parser.Push(parser.itemFactory.create('fn', mml));
};
BaseMethods.NamedOp = function (parser, name, id) {
    if (!id) {
        id = name.substr(1);
    }
    id = id.replace(/&thinsp;/, '\u2006');
    var mml = parser.create('token', 'mo', {
        movablelimits: true,
        movesupsub: true,
        form: TexConstants_js_1.TexConstant.Form.PREFIX,
        texClass: MmlNode_js_1.TEXCLASS.OP
    }, id);
    parser.Push(mml);
};
BaseMethods.Limits = function (parser, _name, limits) {
    var op = parser.stack.Prev(true);
    if (!op || (NodeUtil_js_1.default.getTexClass(NodeUtil_js_1.default.getCoreMO(op)) !== MmlNode_js_1.TEXCLASS.OP &&
        NodeUtil_js_1.default.getProperty(op, 'movesupsub') == null)) {
        throw new TexError_js_1.default('MisplacedLimits', '%1 is allowed only on operators', parser.currentCS);
    }
    var top = parser.stack.Top();
    var node;
    if (NodeUtil_js_1.default.isType(op, 'munderover') && !limits) {
        node = parser.create('node', 'msubsup');
        NodeUtil_js_1.default.copyChildren(op, node);
        op = top.Last = node;
    }
    else if (NodeUtil_js_1.default.isType(op, 'msubsup') && limits) {
        node = parser.create('node', 'munderover');
        NodeUtil_js_1.default.copyChildren(op, node);
        op = top.Last = node;
    }
    NodeUtil_js_1.default.setProperty(op, 'movesupsub', limits ? true : false);
    NodeUtil_js_1.default.setProperties(NodeUtil_js_1.default.getCoreMO(op), { 'movablelimits': false });
    if (NodeUtil_js_1.default.getAttribute(op, 'movablelimits') ||
        NodeUtil_js_1.default.getProperty(op, 'movablelimits')) {
        NodeUtil_js_1.default.setProperties(op, { 'movablelimits': false });
    }
};
BaseMethods.Over = function (parser, name, open, close) {
    var mml = parser.itemFactory.create('over').setProperty('name', parser.currentCS);
    if (open || close) {
        mml.setProperty('open', open);
        mml.setProperty('close', close);
    }
    else if (name.match(/withdelims$/)) {
        mml.setProperty('open', parser.GetDelimiter(name));
        mml.setProperty('close', parser.GetDelimiter(name));
    }
    if (name.match(/^\\above/)) {
        mml.setProperty('thickness', parser.GetDimen(name));
    }
    else if (name.match(/^\\atop/) || open || close) {
        mml.setProperty('thickness', 0);
    }
    parser.Push(mml);
};
BaseMethods.Frac = function (parser, name) {
    var num = parser.ParseArg(name);
    var den = parser.ParseArg(name);
    var node = parser.create('node', 'mfrac', [num, den]);
    parser.Push(node);
};
BaseMethods.Sqrt = function (parser, name) {
    var n = parser.GetBrackets(name);
    var arg = parser.GetArgument(name);
    if (arg === '\\frac') {
        arg += '{' + parser.GetArgument(arg) + '}{' + parser.GetArgument(arg) + '}';
    }
    var mml = new TexParser_js_1.default(arg, parser.stack.env, parser.configuration).mml();
    if (!n) {
        mml = parser.create('node', 'msqrt', [mml]);
    }
    else {
        mml = parser.create('node', 'mroot', [mml, parseRoot(parser, n)]);
    }
    parser.Push(mml);
};
function parseRoot(parser, n) {
    var env = parser.stack.env;
    var inRoot = env['inRoot'];
    env['inRoot'] = true;
    var newParser = new TexParser_js_1.default(n, env, parser.configuration);
    var node = newParser.mml();
    var global = newParser.stack.global;
    if (global['leftRoot'] || global['upRoot']) {
        var def = {};
        if (global['leftRoot']) {
            def['width'] = global['leftRoot'];
        }
        if (global['upRoot']) {
            def['voffset'] = global['upRoot'];
            def['height'] = global['upRoot'];
        }
        node = parser.create('node', 'mpadded', [node], def);
    }
    env['inRoot'] = inRoot;
    return node;
}
BaseMethods.Root = function (parser, name) {
    var n = parser.GetUpTo(name, '\\of');
    var arg = parser.ParseArg(name);
    var node = parser.create('node', 'mroot', [arg, parseRoot(parser, n)]);
    parser.Push(node);
};
BaseMethods.MoveRoot = function (parser, name, id) {
    if (!parser.stack.env['inRoot']) {
        throw new TexError_js_1.default('MisplacedMoveRoot', '%1 can appear only within a root', parser.currentCS);
    }
    if (parser.stack.global[id]) {
        throw new TexError_js_1.default('MultipleMoveRoot', 'Multiple use of %1', parser.currentCS);
    }
    var n = parser.GetArgument(name);
    if (!n.match(/-?[0-9]+/)) {
        throw new TexError_js_1.default('IntegerArg', 'The argument to %1 must be an integer', parser.currentCS);
    }
    n = (parseInt(n, 10) / 15) + 'em';
    if (n.substr(0, 1) !== '-') {
        n = '+' + n;
    }
    parser.stack.global[id] = n;
};
BaseMethods.Accent = function (parser, name, accent, stretchy) {
    var c = parser.ParseArg(name);
    var def = __assign(__assign({}, ParseUtil_js_1.default.getFontDef(parser)), { accent: true, mathaccent: true });
    var entity = NodeUtil_js_1.default.createEntity(accent);
    var moNode = parser.create('token', 'mo', def, entity);
    var mml = moNode;
    NodeUtil_js_1.default.setAttribute(mml, 'stretchy', stretchy ? true : false);
    var mo = (NodeUtil_js_1.default.isEmbellished(c) ? NodeUtil_js_1.default.getCoreMO(c) : c);
    if (NodeUtil_js_1.default.isType(mo, 'mo') || NodeUtil_js_1.default.getProperty(mo, 'movablelimits')) {
        NodeUtil_js_1.default.setProperties(mo, { 'movablelimits': false });
    }
    var muoNode = parser.create('node', 'munderover');
    NodeUtil_js_1.default.setChild(muoNode, 0, c);
    NodeUtil_js_1.default.setChild(muoNode, 1, null);
    NodeUtil_js_1.default.setChild(muoNode, 2, mml);
    var texAtom = parser.create('node', 'TeXAtom', [muoNode]);
    parser.Push(texAtom);
};
BaseMethods.UnderOver = function (parser, name, c, stack) {
    var entity = NodeUtil_js_1.default.createEntity(c);
    var mo = parser.create('token', 'mo', { stretchy: true, accent: true }, entity);
    var pos = (name.charAt(1) === 'o' ? 'over' : 'under');
    var base = parser.ParseArg(name);
    parser.Push(ParseUtil_js_1.default.underOver(parser, base, mo, pos, stack));
};
BaseMethods.Overset = function (parser, name) {
    var top = parser.ParseArg(name);
    var base = parser.ParseArg(name);
    ParseUtil_js_1.default.checkMovableLimits(base);
    if (top.isKind('mo')) {
        NodeUtil_js_1.default.setAttribute(top, 'accent', false);
    }
    var node = parser.create('node', 'mover', [base, top]);
    parser.Push(node);
};
BaseMethods.Underset = function (parser, name) {
    var bot = parser.ParseArg(name);
    var base = parser.ParseArg(name);
    ParseUtil_js_1.default.checkMovableLimits(base);
    if (bot.isKind('mo')) {
        NodeUtil_js_1.default.setAttribute(bot, 'accent', false);
    }
    var node = parser.create('node', 'munder', [base, bot], { accentunder: false });
    parser.Push(node);
};
BaseMethods.Overunderset = function (parser, name) {
    var top = parser.ParseArg(name);
    var bot = parser.ParseArg(name);
    var base = parser.ParseArg(name);
    ParseUtil_js_1.default.checkMovableLimits(base);
    if (top.isKind('mo')) {
        NodeUtil_js_1.default.setAttribute(top, 'accent', false);
    }
    if (bot.isKind('mo')) {
        NodeUtil_js_1.default.setAttribute(bot, 'accent', false);
    }
    var node = parser.create('node', 'munderover', [base, bot, top], { accent: false, accentunder: false });
    parser.Push(node);
};
BaseMethods.TeXAtom = function (parser, name, mclass) {
    var def = { texClass: mclass };
    var mml;
    var node;
    var parsed;
    if (mclass === MmlNode_js_1.TEXCLASS.OP) {
        def['movesupsub'] = def['movablelimits'] = true;
        var arg = parser.GetArgument(name);
        var match = arg.match(/^\s*\\rm\s+([a-zA-Z0-9 ]+)$/);
        if (match) {
            def['mathvariant'] = TexConstants_js_1.TexConstant.Variant.NORMAL;
            node = parser.create('token', 'mi', def, match[1]);
        }
        else {
            parsed = new TexParser_js_1.default(arg, parser.stack.env, parser.configuration).mml();
            node = parser.create('node', 'TeXAtom', [parsed], def);
        }
        mml = parser.itemFactory.create('fn', node);
    }
    else {
        parsed = parser.ParseArg(name);
        mml = parser.create('node', 'TeXAtom', [parsed], def);
    }
    parser.Push(mml);
};
BaseMethods.MmlToken = function (parser, name) {
    var kind = parser.GetArgument(name);
    var attr = parser.GetBrackets(name, '').replace(/^\s+/, '');
    var text = parser.GetArgument(name);
    var def = {};
    var keep = [];
    var node;
    try {
        node = parser.create('node', kind);
    }
    catch (e) {
        node = null;
    }
    if (!node || !node.isToken) {
        throw new TexError_js_1.default('NotMathMLToken', '%1 is not a token element', kind);
    }
    while (attr !== '') {
        var match = attr.match(/^([a-z]+)\s*=\s*('[^']*'|"[^"]*"|[^ ,]*)\s*,?\s*/i);
        if (!match) {
            throw new TexError_js_1.default('InvalidMathMLAttr', 'Invalid MathML attribute: %1', attr);
        }
        if (!node.attributes.hasDefault(match[1]) && !MmlTokenAllow[match[1]]) {
            throw new TexError_js_1.default('UnknownAttrForElement', '%1 is not a recognized attribute for %2', match[1], kind);
        }
        var value = ParseUtil_js_1.default.MmlFilterAttribute(parser, match[1], match[2].replace(/^(['"])(.*)\1$/, '$2'));
        if (value) {
            if (value.toLowerCase() === 'true') {
                value = true;
            }
            else if (value.toLowerCase() === 'false') {
                value = false;
            }
            def[match[1]] = value;
            keep.push(match[1]);
        }
        attr = attr.substr(match[0].length);
    }
    if (keep.length) {
        def['mjx-keep-attrs'] = keep.join(' ');
    }
    var textNode = parser.create('text', text);
    node.appendChild(textNode);
    NodeUtil_js_1.default.setProperties(node, def);
    parser.Push(node);
};
BaseMethods.Strut = function (parser, _name) {
    var row = parser.create('node', 'mrow');
    var padded = parser.create('node', 'mpadded', [row], { height: '8.6pt', depth: '3pt', width: 0 });
    parser.Push(padded);
};
BaseMethods.Phantom = function (parser, name, v, h) {
    var box = parser.create('node', 'mphantom', [parser.ParseArg(name)]);
    if (v || h) {
        box = parser.create('node', 'mpadded', [box]);
        if (h) {
            NodeUtil_js_1.default.setAttribute(box, 'height', 0);
            NodeUtil_js_1.default.setAttribute(box, 'depth', 0);
        }
        if (v) {
            NodeUtil_js_1.default.setAttribute(box, 'width', 0);
        }
    }
    var atom = parser.create('node', 'TeXAtom', [box]);
    parser.Push(atom);
};
BaseMethods.Smash = function (parser, name) {
    var bt = ParseUtil_js_1.default.trimSpaces(parser.GetBrackets(name, ''));
    var smash = parser.create('node', 'mpadded', [parser.ParseArg(name)]);
    switch (bt) {
        case 'b':
            NodeUtil_js_1.default.setAttribute(smash, 'depth', 0);
            break;
        case 't':
            NodeUtil_js_1.default.setAttribute(smash, 'height', 0);
            break;
        default:
            NodeUtil_js_1.default.setAttribute(smash, 'height', 0);
            NodeUtil_js_1.default.setAttribute(smash, 'depth', 0);
    }
    var atom = parser.create('node', 'TeXAtom', [smash]);
    parser.Push(atom);
};
BaseMethods.Lap = function (parser, name) {
    var mml = parser.create('node', 'mpadded', [parser.ParseArg(name)], { width: 0 });
    if (name === '\\llap') {
        NodeUtil_js_1.default.setAttribute(mml, 'lspace', '-1width');
    }
    var atom = parser.create('node', 'TeXAtom', [mml]);
    parser.Push(atom);
};
BaseMethods.RaiseLower = function (parser, name) {
    var h = parser.GetDimen(name);
    var item = parser.itemFactory.create('position').setProperties({ name: parser.currentCS, move: 'vertical' });
    if (h.charAt(0) === '-') {
        h = h.slice(1);
        name = name.substr(1) === 'raise' ? '\\lower' : '\\raise';
    }
    if (name === '\\lower') {
        item.setProperty('dh', '-' + h);
        item.setProperty('dd', '+' + h);
    }
    else {
        item.setProperty('dh', '+' + h);
        item.setProperty('dd', '-' + h);
    }
    parser.Push(item);
};
BaseMethods.MoveLeftRight = function (parser, name) {
    var h = parser.GetDimen(name);
    var nh = (h.charAt(0) === '-' ? h.slice(1) : '-' + h);
    if (name === '\\moveleft') {
        var tmp = h;
        h = nh;
        nh = tmp;
    }
    parser.Push(parser.itemFactory.create('position').setProperties({
        name: parser.currentCS, move: 'horizontal',
        left: parser.create('node', 'mspace', [], { width: h }),
        right: parser.create('node', 'mspace', [], { width: nh })
    }));
};
BaseMethods.Hskip = function (parser, name) {
    var node = parser.create('node', 'mspace', [], { width: parser.GetDimen(name) });
    parser.Push(node);
};
BaseMethods.Nonscript = function (parser, _name) {
    parser.Push(parser.itemFactory.create('nonscript'));
};
BaseMethods.Rule = function (parser, name, style) {
    var w = parser.GetDimen(name), h = parser.GetDimen(name), d = parser.GetDimen(name);
    var def = { width: w, height: h, depth: d };
    if (style !== 'blank') {
        def['mathbackground'] = (parser.stack.env['color'] || 'black');
    }
    var node = parser.create('node', 'mspace', [], def);
    parser.Push(node);
};
BaseMethods.rule = function (parser, name) {
    var v = parser.GetBrackets(name), w = parser.GetDimen(name), h = parser.GetDimen(name);
    var mml = parser.create('node', 'mspace', [], {
        width: w, height: h,
        mathbackground: (parser.stack.env['color'] || 'black')
    });
    if (v) {
        mml = parser.create('node', 'mpadded', [mml], { voffset: v });
        if (v.match(/^\-/)) {
            NodeUtil_js_1.default.setAttribute(mml, 'height', v);
            NodeUtil_js_1.default.setAttribute(mml, 'depth', '+' + v.substr(1));
        }
        else {
            NodeUtil_js_1.default.setAttribute(mml, 'height', '+' + v);
        }
    }
    parser.Push(mml);
};
BaseMethods.MakeBig = function (parser, name, mclass, size) {
    size *= P_HEIGHT;
    var sizeStr = String(size).replace(/(\.\d\d\d).+/, '$1') + 'em';
    var delim = parser.GetDelimiter(name, true);
    var mo = parser.create('token', 'mo', {
        minsize: sizeStr, maxsize: sizeStr,
        fence: true, stretchy: true, symmetric: true
    }, delim);
    var node = parser.create('node', 'TeXAtom', [mo], { texClass: mclass });
    parser.Push(node);
};
BaseMethods.BuildRel = function (parser, name) {
    var top = parser.ParseUpTo(name, '\\over');
    var bot = parser.ParseArg(name);
    var node = parser.create('node', 'munderover');
    NodeUtil_js_1.default.setChild(node, 0, bot);
    NodeUtil_js_1.default.setChild(node, 1, null);
    NodeUtil_js_1.default.setChild(node, 2, top);
    var atom = parser.create('node', 'TeXAtom', [node], { texClass: MmlNode_js_1.TEXCLASS.REL });
    parser.Push(atom);
};
BaseMethods.HBox = function (parser, name, style, font) {
    parser.PushAll(ParseUtil_js_1.default.internalMath(parser, parser.GetArgument(name), style, font));
};
BaseMethods.FBox = function (parser, name) {
    var internal = ParseUtil_js_1.default.internalMath(parser, parser.GetArgument(name));
    var node = parser.create('node', 'menclose', internal, { notation: 'box' });
    parser.Push(node);
};
BaseMethods.FrameBox = function (parser, name) {
    var width = parser.GetBrackets(name);
    var pos = parser.GetBrackets(name) || 'c';
    var mml = ParseUtil_js_1.default.internalMath(parser, parser.GetArgument(name));
    if (width) {
        mml = [parser.create('node', 'mpadded', mml, {
                width: width,
                'data-align': (0, Options_js_1.lookup)(pos, { l: 'left', r: 'right' }, 'center')
            })];
    }
    var node = parser.create('node', 'TeXAtom', [parser.create('node', 'menclose', mml, { notation: 'box' })], { texClass: MmlNode_js_1.TEXCLASS.ORD });
    parser.Push(node);
};
BaseMethods.Not = function (parser, _name) {
    parser.Push(parser.itemFactory.create('not'));
};
BaseMethods.Dots = function (parser, _name) {
    var ldotsEntity = NodeUtil_js_1.default.createEntity('2026');
    var cdotsEntity = NodeUtil_js_1.default.createEntity('22EF');
    var ldots = parser.create('token', 'mo', { stretchy: false }, ldotsEntity);
    var cdots = parser.create('token', 'mo', { stretchy: false }, cdotsEntity);
    parser.Push(parser.itemFactory.create('dots').setProperties({
        ldots: ldots,
        cdots: cdots
    }));
};
BaseMethods.Matrix = function (parser, _name, open, close, align, spacing, vspacing, style, cases, numbered) {
    var c = parser.GetNext();
    if (c === '') {
        throw new TexError_js_1.default('MissingArgFor', 'Missing argument for %1', parser.currentCS);
    }
    if (c === '{') {
        parser.i++;
    }
    else {
        parser.string = c + '}' + parser.string.slice(parser.i + 1);
        parser.i = 0;
    }
    var array = parser.itemFactory.create('array').setProperty('requireClose', true);
    array.arraydef = {
        rowspacing: (vspacing || '4pt'),
        columnspacing: (spacing || '1em')
    };
    if (cases) {
        array.setProperty('isCases', true);
    }
    if (numbered) {
        array.setProperty('isNumbered', true);
        array.arraydef.side = numbered;
    }
    if (open || close) {
        array.setProperty('open', open);
        array.setProperty('close', close);
    }
    if (style === 'D') {
        array.arraydef.displaystyle = true;
    }
    if (align != null) {
        array.arraydef.columnalign = align;
    }
    parser.Push(array);
};
BaseMethods.Entry = function (parser, name) {
    parser.Push(parser.itemFactory.create('cell').setProperties({ isEntry: true, name: name }));
    var top = parser.stack.Top();
    var env = top.getProperty('casesEnv');
    var cases = top.getProperty('isCases');
    if (!cases && !env)
        return;
    var str = parser.string;
    var braces = 0, close = -1, i = parser.i, m = str.length;
    var end = (env ? new RegExp("^\\\\end\\s*\\{".concat(env.replace(/\*/, '\\*'), "\\}")) : null);
    while (i < m) {
        var c = str.charAt(i);
        if (c === '{') {
            braces++;
            i++;
        }
        else if (c === '}') {
            if (braces === 0) {
                m = 0;
            }
            else {
                braces--;
                if (braces === 0 && close < 0) {
                    close = i - parser.i;
                }
                i++;
            }
        }
        else if (c === '&' && braces === 0) {
            throw new TexError_js_1.default('ExtraAlignTab', 'Extra alignment tab in \\cases text');
        }
        else if (c === '\\') {
            var rest = str.substr(i);
            if (rest.match(/^((\\cr)[^a-zA-Z]|\\\\)/) || (end && rest.match(end))) {
                m = 0;
            }
            else {
                i += 2;
            }
        }
        else {
            i++;
        }
    }
    var text = str.substr(parser.i, i - parser.i);
    if (!text.match(/^\s*\\text[^a-zA-Z]/) || close !== text.replace(/\s+$/, '').length - 1) {
        var internal = ParseUtil_js_1.default.internalMath(parser, ParseUtil_js_1.default.trimSpaces(text), 0);
        parser.PushAll(internal);
        parser.i = i;
    }
};
BaseMethods.Cr = function (parser, name) {
    parser.Push(parser.itemFactory.create('cell').setProperties({ isCR: true, name: name }));
};
BaseMethods.CrLaTeX = function (parser, name, nobrackets) {
    if (nobrackets === void 0) { nobrackets = false; }
    var n;
    if (!nobrackets) {
        if (parser.string.charAt(parser.i) === '*') {
            parser.i++;
        }
        if (parser.string.charAt(parser.i) === '[') {
            var dim = parser.GetBrackets(name, '');
            var _a = __read(ParseUtil_js_1.default.matchDimen(dim), 2), value = _a[0], unit = _a[1];
            if (dim && !value) {
                throw new TexError_js_1.default('BracketMustBeDimension', 'Bracket argument to %1 must be a dimension', parser.currentCS);
            }
            n = value + unit;
        }
    }
    parser.Push(parser.itemFactory.create('cell').setProperties({ isCR: true, name: name, linebreak: true }));
    var top = parser.stack.Top();
    var node;
    if (top instanceof sitem.ArrayItem) {
        if (n) {
            top.addRowSpacing(n);
        }
    }
    else {
        if (n) {
            node = parser.create('node', 'mspace', [], { depth: n });
            parser.Push(node);
        }
        node = parser.create('node', 'mspace', [], { linebreak: TexConstants_js_1.TexConstant.LineBreak.NEWLINE });
        parser.Push(node);
    }
};
BaseMethods.HLine = function (parser, _name, style) {
    if (style == null) {
        style = 'solid';
    }
    var top = parser.stack.Top();
    if (!(top instanceof sitem.ArrayItem) || top.Size()) {
        throw new TexError_js_1.default('Misplaced', 'Misplaced %1', parser.currentCS);
    }
    if (!top.table.length) {
        top.frame.push('top');
    }
    else {
        var lines = (top.arraydef['rowlines'] ? top.arraydef['rowlines'].split(/ /) : []);
        while (lines.length < top.table.length) {
            lines.push('none');
        }
        lines[top.table.length - 1] = style;
        top.arraydef['rowlines'] = lines.join(' ');
    }
};
BaseMethods.HFill = function (parser, _name) {
    var top = parser.stack.Top();
    if (top instanceof sitem.ArrayItem) {
        top.hfill.push(top.Size());
    }
    else {
        throw new TexError_js_1.default('UnsupportedHFill', 'Unsupported use of %1', parser.currentCS);
    }
};
BaseMethods.BeginEnd = function (parser, name) {
    var env = parser.GetArgument(name);
    if (env.match(/\\/i)) {
        throw new TexError_js_1.default('InvalidEnv', 'Invalid environment name \'%1\'', env);
    }
    var macro = parser.configuration.handlers.get('environment').lookup(env);
    if (macro && name === '\\end') {
        if (!macro.args[0]) {
            var mml = parser.itemFactory.create('end').setProperty('name', env);
            parser.Push(mml);
            return;
        }
        parser.stack.env['closing'] = env;
    }
    ParseUtil_js_1.default.checkMaxMacros(parser, false);
    parser.parse('environment', [parser, env]);
};
BaseMethods.Array = function (parser, begin, open, close, align, spacing, vspacing, style, raggedHeight) {
    if (!align) {
        align = parser.GetArgument('\\begin{' + begin.getName() + '}');
    }
    var lines = ('c' + align).replace(/[^clr|:]/g, '').replace(/[^|:]([|:])+/g, '$1');
    align = align.replace(/[^clr]/g, '').split('').join(' ');
    align = align.replace(/l/g, 'left').replace(/r/g, 'right').replace(/c/g, 'center');
    var array = parser.itemFactory.create('array');
    array.arraydef = {
        columnalign: align,
        columnspacing: (spacing || '1em'),
        rowspacing: (vspacing || '4pt')
    };
    if (lines.match(/[|:]/)) {
        if (lines.charAt(0).match(/[|:]/)) {
            array.frame.push('left');
            array.dashed = lines.charAt(0) === ':';
        }
        if (lines.charAt(lines.length - 1).match(/[|:]/)) {
            array.frame.push('right');
        }
        lines = lines.substr(1, lines.length - 2);
        array.arraydef.columnlines =
            lines.split('').join(' ').replace(/[^|: ]/g, 'none').replace(/\|/g, 'solid').replace(/:/g, 'dashed');
    }
    if (open) {
        array.setProperty('open', parser.convertDelimiter(open));
    }
    if (close) {
        array.setProperty('close', parser.convertDelimiter(close));
    }
    if ((style || '').charAt(1) === '\'') {
        array.arraydef['data-cramped'] = true;
        style = style.charAt(0);
    }
    if (style === 'D') {
        array.arraydef['displaystyle'] = true;
    }
    else if (style) {
        array.arraydef['displaystyle'] = false;
    }
    if (style === 'S') {
        array.arraydef['scriptlevel'] = 1;
    }
    if (raggedHeight) {
        array.arraydef['useHeight'] = false;
    }
    parser.Push(begin);
    return array;
};
BaseMethods.AlignedArray = function (parser, begin) {
    var align = parser.GetBrackets('\\begin{' + begin.getName() + '}');
    var item = BaseMethods.Array(parser, begin);
    return ParseUtil_js_1.default.setArrayAlign(item, align);
};
BaseMethods.Equation = function (parser, begin, numbered) {
    parser.Push(begin);
    ParseUtil_js_1.default.checkEqnEnv(parser);
    return parser.itemFactory.create('equation', numbered).
        setProperty('name', begin.getName());
};
BaseMethods.EqnArray = function (parser, begin, numbered, taggable, align, spacing) {
    parser.Push(begin);
    if (taggable) {
        ParseUtil_js_1.default.checkEqnEnv(parser);
    }
    align = align.replace(/[^clr]/g, '').split('').join(' ');
    align = align.replace(/l/g, 'left').replace(/r/g, 'right').replace(/c/g, 'center');
    var newItem = parser.itemFactory.create('eqnarray', begin.getName(), numbered, taggable, parser.stack.global);
    newItem.arraydef = {
        displaystyle: true,
        columnalign: align,
        columnspacing: (spacing || '1em'),
        rowspacing: '3pt',
        side: parser.options['tagSide'],
        minlabelspacing: parser.options['tagIndent']
    };
    return newItem;
};
BaseMethods.HandleNoTag = function (parser, _name) {
    parser.tags.notag();
};
BaseMethods.HandleLabel = function (parser, name) {
    var label = parser.GetArgument(name);
    if (label === '') {
        return;
    }
    if (!parser.tags.refUpdate) {
        if (parser.tags.label) {
            throw new TexError_js_1.default('MultipleCommand', 'Multiple %1', parser.currentCS);
        }
        parser.tags.label = label;
        if ((parser.tags.allLabels[label] || parser.tags.labels[label]) && !parser.options['ignoreDuplicateLabels']) {
            throw new TexError_js_1.default('MultipleLabel', 'Label \'%1\' multiply defined', label);
        }
        parser.tags.labels[label] = new Tags_js_1.Label();
    }
};
BaseMethods.HandleRef = function (parser, name, eqref) {
    var label = parser.GetArgument(name);
    var ref = parser.tags.allLabels[label] || parser.tags.labels[label];
    if (!ref) {
        if (!parser.tags.refUpdate) {
            parser.tags.redo = true;
        }
        ref = new Tags_js_1.Label();
    }
    var tag = ref.tag;
    if (eqref) {
        tag = parser.tags.formatTag(tag);
    }
    var node = parser.create('node', 'mrow', ParseUtil_js_1.default.internalMath(parser, tag), {
        href: parser.tags.formatUrl(ref.id, parser.options.baseURL), 'class': 'MathJax_ref'
    });
    parser.Push(node);
};
BaseMethods.Macro = function (parser, name, macro, argcount, def) {
    if (argcount) {
        var args = [];
        if (def != null) {
            var optional = parser.GetBrackets(name);
            args.push(optional == null ? def : optional);
        }
        for (var i = args.length; i < argcount; i++) {
            args.push(parser.GetArgument(name));
        }
        macro = ParseUtil_js_1.default.substituteArgs(parser, args, macro);
    }
    parser.string = ParseUtil_js_1.default.addArgs(parser, macro, parser.string.slice(parser.i));
    parser.i = 0;
    ParseUtil_js_1.default.checkMaxMacros(parser);
};
BaseMethods.MathChoice = function (parser, name) {
    var D = parser.ParseArg(name);
    var T = parser.ParseArg(name);
    var S = parser.ParseArg(name);
    var SS = parser.ParseArg(name);
    parser.Push(parser.create('node', 'MathChoice', [D, T, S, SS]));
};
exports["default"] = BaseMethods;
//# sourceMappingURL=BaseMethods.js.map

/***/ }),

/***/ 56780:
/***/ ((__unused_webpack_module, exports) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.px = exports.emRounded = exports.em = exports.percent = exports.length2em = exports.MATHSPACE = exports.RELUNITS = exports.UNITS = exports.BIGDIMEN = void 0;
exports.BIGDIMEN = 1000000;
exports.UNITS = {
    px: 1,
    'in': 96,
    cm: 96 / 2.54,
    mm: 96 / 25.4
};
exports.RELUNITS = {
    em: 1,
    ex: .431,
    pt: 1 / 10,
    pc: 12 / 10,
    mu: 1 / 18
};
exports.MATHSPACE = {
    veryverythinmathspace: 1 / 18,
    verythinmathspace: 2 / 18,
    thinmathspace: 3 / 18,
    mediummathspace: 4 / 18,
    thickmathspace: 5 / 18,
    verythickmathspace: 6 / 18,
    veryverythickmathspace: 7 / 18,
    negativeveryverythinmathspace: -1 / 18,
    negativeverythinmathspace: -2 / 18,
    negativethinmathspace: -3 / 18,
    negativemediummathspace: -4 / 18,
    negativethickmathspace: -5 / 18,
    negativeverythickmathspace: -6 / 18,
    negativeveryverythickmathspace: -7 / 18,
    thin: .04,
    medium: .06,
    thick: .1,
    normal: 1,
    big: 2,
    small: 1 / Math.sqrt(2),
    infinity: exports.BIGDIMEN
};
function length2em(length, size, scale, em) {
    if (size === void 0) { size = 0; }
    if (scale === void 0) { scale = 1; }
    if (em === void 0) { em = 16; }
    if (typeof length !== 'string') {
        length = String(length);
    }
    if (length === '' || length == null) {
        return size;
    }
    if (exports.MATHSPACE[length]) {
        return exports.MATHSPACE[length];
    }
    var match = length.match(/^\s*([-+]?(?:\.\d+|\d+(?:\.\d*)?))?(pt|em|ex|mu|px|pc|in|mm|cm|%)?/);
    if (!match) {
        return size;
    }
    var m = parseFloat(match[1] || '1'), unit = match[2];
    if (exports.UNITS.hasOwnProperty(unit)) {
        return m * exports.UNITS[unit] / em / scale;
    }
    if (exports.RELUNITS.hasOwnProperty(unit)) {
        return m * exports.RELUNITS[unit];
    }
    if (unit === '%') {
        return m / 100 * size;
    }
    return m * size;
}
exports.length2em = length2em;
function percent(m) {
    return (100 * m).toFixed(1).replace(/\.?0+$/, '') + '%';
}
exports.percent = percent;
function em(m) {
    if (Math.abs(m) < .001)
        return '0';
    return (m.toFixed(3).replace(/\.?0+$/, '')) + 'em';
}
exports.em = em;
function emRounded(m, em) {
    if (em === void 0) { em = 16; }
    m = (Math.round(m * em) + .05) / em;
    if (Math.abs(m) < .001)
        return '0em';
    return m.toFixed(3).replace(/\.?0+$/, '') + 'em';
}
exports.emRounded = emRounded;
function px(m, M, em) {
    if (M === void 0) { M = -exports.BIGDIMEN; }
    if (em === void 0) { em = 16; }
    m *= em;
    if (M && m < M)
        m = M;
    if (Math.abs(m) < .1)
        return '0';
    return m.toFixed(1).replace(/\.0$/, '') + 'px';
}
exports.px = px;
//# sourceMappingURL=lengths.js.map

/***/ })

}]);
//# sourceMappingURL=7154.cbe00782c73879782ae7.js.map?v=cbe00782c73879782ae7
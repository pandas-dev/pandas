"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[270],{

/***/ 95564:
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
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.HTMLAdaptor = void 0;
var DOMAdaptor_js_1 = __webpack_require__(66791);
var HTMLAdaptor = (function (_super) {
    __extends(HTMLAdaptor, _super);
    function HTMLAdaptor(window) {
        var _this = _super.call(this, window.document) || this;
        _this.window = window;
        _this.parser = new window.DOMParser();
        return _this;
    }
    HTMLAdaptor.prototype.parse = function (text, format) {
        if (format === void 0) { format = 'text/html'; }
        return this.parser.parseFromString(text, format);
    };
    HTMLAdaptor.prototype.create = function (kind, ns) {
        return (ns ?
            this.document.createElementNS(ns, kind) :
            this.document.createElement(kind));
    };
    HTMLAdaptor.prototype.text = function (text) {
        return this.document.createTextNode(text);
    };
    HTMLAdaptor.prototype.head = function (doc) {
        return doc.head || doc;
    };
    HTMLAdaptor.prototype.body = function (doc) {
        return doc.body || doc;
    };
    HTMLAdaptor.prototype.root = function (doc) {
        return doc.documentElement || doc;
    };
    HTMLAdaptor.prototype.doctype = function (doc) {
        return (doc.doctype ? "<!DOCTYPE ".concat(doc.doctype.name, ">") : '');
    };
    HTMLAdaptor.prototype.tags = function (node, name, ns) {
        if (ns === void 0) { ns = null; }
        var nodes = (ns ? node.getElementsByTagNameNS(ns, name) : node.getElementsByTagName(name));
        return Array.from(nodes);
    };
    HTMLAdaptor.prototype.getElements = function (nodes, _document) {
        var e_1, _a;
        var containers = [];
        try {
            for (var nodes_1 = __values(nodes), nodes_1_1 = nodes_1.next(); !nodes_1_1.done; nodes_1_1 = nodes_1.next()) {
                var node = nodes_1_1.value;
                if (typeof (node) === 'string') {
                    containers = containers.concat(Array.from(this.document.querySelectorAll(node)));
                }
                else if (Array.isArray(node)) {
                    containers = containers.concat(Array.from(node));
                }
                else if (node instanceof this.window.NodeList || node instanceof this.window.HTMLCollection) {
                    containers = containers.concat(Array.from(node));
                }
                else {
                    containers.push(node);
                }
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (nodes_1_1 && !nodes_1_1.done && (_a = nodes_1.return)) _a.call(nodes_1);
            }
            finally { if (e_1) throw e_1.error; }
        }
        return containers;
    };
    HTMLAdaptor.prototype.contains = function (container, node) {
        return container.contains(node);
    };
    HTMLAdaptor.prototype.parent = function (node) {
        return node.parentNode;
    };
    HTMLAdaptor.prototype.append = function (node, child) {
        return node.appendChild(child);
    };
    HTMLAdaptor.prototype.insert = function (nchild, ochild) {
        return this.parent(ochild).insertBefore(nchild, ochild);
    };
    HTMLAdaptor.prototype.remove = function (child) {
        return this.parent(child).removeChild(child);
    };
    HTMLAdaptor.prototype.replace = function (nnode, onode) {
        return this.parent(onode).replaceChild(nnode, onode);
    };
    HTMLAdaptor.prototype.clone = function (node) {
        return node.cloneNode(true);
    };
    HTMLAdaptor.prototype.split = function (node, n) {
        return node.splitText(n);
    };
    HTMLAdaptor.prototype.next = function (node) {
        return node.nextSibling;
    };
    HTMLAdaptor.prototype.previous = function (node) {
        return node.previousSibling;
    };
    HTMLAdaptor.prototype.firstChild = function (node) {
        return node.firstChild;
    };
    HTMLAdaptor.prototype.lastChild = function (node) {
        return node.lastChild;
    };
    HTMLAdaptor.prototype.childNodes = function (node) {
        return Array.from(node.childNodes);
    };
    HTMLAdaptor.prototype.childNode = function (node, i) {
        return node.childNodes[i];
    };
    HTMLAdaptor.prototype.kind = function (node) {
        var n = node.nodeType;
        return (n === 1 || n === 3 || n === 8 ? node.nodeName.toLowerCase() : '');
    };
    HTMLAdaptor.prototype.value = function (node) {
        return node.nodeValue || '';
    };
    HTMLAdaptor.prototype.textContent = function (node) {
        return node.textContent;
    };
    HTMLAdaptor.prototype.innerHTML = function (node) {
        return node.innerHTML;
    };
    HTMLAdaptor.prototype.outerHTML = function (node) {
        return node.outerHTML;
    };
    HTMLAdaptor.prototype.serializeXML = function (node) {
        var serializer = new this.window.XMLSerializer();
        return serializer.serializeToString(node);
    };
    HTMLAdaptor.prototype.setAttribute = function (node, name, value, ns) {
        if (ns === void 0) { ns = null; }
        if (!ns) {
            return node.setAttribute(name, value);
        }
        name = ns.replace(/.*\//, '') + ':' + name.replace(/^.*:/, '');
        return node.setAttributeNS(ns, name, value);
    };
    HTMLAdaptor.prototype.getAttribute = function (node, name) {
        return node.getAttribute(name);
    };
    HTMLAdaptor.prototype.removeAttribute = function (node, name) {
        return node.removeAttribute(name);
    };
    HTMLAdaptor.prototype.hasAttribute = function (node, name) {
        return node.hasAttribute(name);
    };
    HTMLAdaptor.prototype.allAttributes = function (node) {
        return Array.from(node.attributes).map(function (x) {
            return { name: x.name, value: x.value };
        });
    };
    HTMLAdaptor.prototype.addClass = function (node, name) {
        if (node.classList) {
            node.classList.add(name);
        }
        else {
            node.className = (node.className + ' ' + name).trim();
        }
    };
    HTMLAdaptor.prototype.removeClass = function (node, name) {
        if (node.classList) {
            node.classList.remove(name);
        }
        else {
            node.className = node.className.split(/ /).filter(function (c) { return c !== name; }).join(' ');
        }
    };
    HTMLAdaptor.prototype.hasClass = function (node, name) {
        if (node.classList) {
            return node.classList.contains(name);
        }
        return node.className.split(/ /).indexOf(name) >= 0;
    };
    HTMLAdaptor.prototype.setStyle = function (node, name, value) {
        node.style[name] = value;
    };
    HTMLAdaptor.prototype.getStyle = function (node, name) {
        return node.style[name];
    };
    HTMLAdaptor.prototype.allStyles = function (node) {
        return node.style.cssText;
    };
    HTMLAdaptor.prototype.insertRules = function (node, rules) {
        var e_2, _a;
        try {
            for (var _b = __values(rules.reverse()), _c = _b.next(); !_c.done; _c = _b.next()) {
                var rule = _c.value;
                try {
                    node.sheet.insertRule(rule, 0);
                }
                catch (e) {
                    console.warn("MathJax: can't insert css rule '".concat(rule, "': ").concat(e.message));
                }
            }
        }
        catch (e_2_1) { e_2 = { error: e_2_1 }; }
        finally {
            try {
                if (_c && !_c.done && (_a = _b.return)) _a.call(_b);
            }
            finally { if (e_2) throw e_2.error; }
        }
    };
    HTMLAdaptor.prototype.fontSize = function (node) {
        var style = this.window.getComputedStyle(node);
        return parseFloat(style.fontSize);
    };
    HTMLAdaptor.prototype.fontFamily = function (node) {
        var style = this.window.getComputedStyle(node);
        return style.fontFamily || '';
    };
    HTMLAdaptor.prototype.nodeSize = function (node, em, local) {
        if (em === void 0) { em = 1; }
        if (local === void 0) { local = false; }
        if (local && node.getBBox) {
            var _a = node.getBBox(), width = _a.width, height = _a.height;
            return [width / em, height / em];
        }
        return [node.offsetWidth / em, node.offsetHeight / em];
    };
    HTMLAdaptor.prototype.nodeBBox = function (node) {
        var _a = node.getBoundingClientRect(), left = _a.left, right = _a.right, top = _a.top, bottom = _a.bottom;
        return { left: left, right: right, top: top, bottom: bottom };
    };
    return HTMLAdaptor;
}(DOMAdaptor_js_1.AbstractDOMAdaptor));
exports.HTMLAdaptor = HTMLAdaptor;
//# sourceMappingURL=HTMLAdaptor.js.map

/***/ }),

/***/ 90270:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {


Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.browserAdaptor = void 0;
var HTMLAdaptor_js_1 = __webpack_require__(95564);
function browserAdaptor() {
    return new HTMLAdaptor_js_1.HTMLAdaptor(window);
}
exports.browserAdaptor = browserAdaptor;
//# sourceMappingURL=browserAdaptor.js.map

/***/ }),

/***/ 66791:
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
Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.AbstractDOMAdaptor = void 0;
var AbstractDOMAdaptor = (function () {
    function AbstractDOMAdaptor(document) {
        if (document === void 0) { document = null; }
        this.document = document;
    }
    AbstractDOMAdaptor.prototype.node = function (kind, def, children, ns) {
        var e_1, _a;
        if (def === void 0) { def = {}; }
        if (children === void 0) { children = []; }
        var node = this.create(kind, ns);
        this.setAttributes(node, def);
        try {
            for (var children_1 = __values(children), children_1_1 = children_1.next(); !children_1_1.done; children_1_1 = children_1.next()) {
                var child = children_1_1.value;
                this.append(node, child);
            }
        }
        catch (e_1_1) { e_1 = { error: e_1_1 }; }
        finally {
            try {
                if (children_1_1 && !children_1_1.done && (_a = children_1.return)) _a.call(children_1);
            }
            finally { if (e_1) throw e_1.error; }
        }
        return node;
    };
    AbstractDOMAdaptor.prototype.setAttributes = function (node, def) {
        var e_2, _a, e_3, _b, e_4, _c;
        if (def.style && typeof (def.style) !== 'string') {
            try {
                for (var _d = __values(Object.keys(def.style)), _e = _d.next(); !_e.done; _e = _d.next()) {
                    var key = _e.value;
                    this.setStyle(node, key.replace(/-([a-z])/g, function (_m, c) { return c.toUpperCase(); }), def.style[key]);
                }
            }
            catch (e_2_1) { e_2 = { error: e_2_1 }; }
            finally {
                try {
                    if (_e && !_e.done && (_a = _d.return)) _a.call(_d);
                }
                finally { if (e_2) throw e_2.error; }
            }
        }
        if (def.properties) {
            try {
                for (var _f = __values(Object.keys(def.properties)), _g = _f.next(); !_g.done; _g = _f.next()) {
                    var key = _g.value;
                    node[key] = def.properties[key];
                }
            }
            catch (e_3_1) { e_3 = { error: e_3_1 }; }
            finally {
                try {
                    if (_g && !_g.done && (_b = _f.return)) _b.call(_f);
                }
                finally { if (e_3) throw e_3.error; }
            }
        }
        try {
            for (var _h = __values(Object.keys(def)), _j = _h.next(); !_j.done; _j = _h.next()) {
                var key = _j.value;
                if ((key !== 'style' || typeof (def.style) === 'string') && key !== 'properties') {
                    this.setAttribute(node, key, def[key]);
                }
            }
        }
        catch (e_4_1) { e_4 = { error: e_4_1 }; }
        finally {
            try {
                if (_j && !_j.done && (_c = _h.return)) _c.call(_h);
            }
            finally { if (e_4) throw e_4.error; }
        }
    };
    AbstractDOMAdaptor.prototype.replace = function (nnode, onode) {
        this.insert(nnode, onode);
        this.remove(onode);
        return onode;
    };
    AbstractDOMAdaptor.prototype.childNode = function (node, i) {
        return this.childNodes(node)[i];
    };
    AbstractDOMAdaptor.prototype.allClasses = function (node) {
        var classes = this.getAttribute(node, 'class');
        return (!classes ? [] :
            classes.replace(/  +/g, ' ').replace(/^ /, '').replace(/ $/, '').split(/ /));
    };
    return AbstractDOMAdaptor;
}());
exports.AbstractDOMAdaptor = AbstractDOMAdaptor;
//# sourceMappingURL=DOMAdaptor.js.map

/***/ })

}]);
//# sourceMappingURL=270.dced80a7f5cbf1705712.js.map?v=dced80a7f5cbf1705712
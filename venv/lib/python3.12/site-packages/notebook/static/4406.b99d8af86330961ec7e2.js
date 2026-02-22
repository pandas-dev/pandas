(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[4406],{

/***/ 2091:
/***/ ((__unused_webpack_module, exports) => {

"use strict";

Object.defineProperty(exports, "__esModule", ({ value: true }));
exports.BLANK_URL = exports.relativeFirstCharacters = exports.whitespaceEscapeCharsRegex = exports.urlSchemeRegex = exports.ctrlCharactersRegex = exports.htmlCtrlEntityRegex = exports.htmlEntitiesRegex = exports.invalidProtocolRegex = void 0;
exports.invalidProtocolRegex = /^([^\w]*)(javascript|data|vbscript)/im;
exports.htmlEntitiesRegex = /&#(\w+)(^\w|;)?/g;
exports.htmlCtrlEntityRegex = /&(newline|tab);/gi;
exports.ctrlCharactersRegex = /[\u0000-\u001F\u007F-\u009F\u2000-\u200D\uFEFF]/gim;
exports.urlSchemeRegex = /^.+(:|&colon;)/gim;
exports.whitespaceEscapeCharsRegex = /(\\|%5[cC])((%(6[eE]|72|74))|[nrt])/g;
exports.relativeFirstCharacters = [".", "/"];
exports.BLANK_URL = "about:blank";


/***/ }),

/***/ 7608:
/***/ ((__unused_webpack_module, exports, __webpack_require__) => {

"use strict";
var __webpack_unused_export__;

__webpack_unused_export__ = ({ value: true });
exports.N = void 0;
var constants_1 = __webpack_require__(2091);
function isRelativeUrlWithoutProtocol(url) {
    return constants_1.relativeFirstCharacters.indexOf(url[0]) > -1;
}
function decodeHtmlCharacters(str) {
    var removedNullByte = str.replace(constants_1.ctrlCharactersRegex, "");
    return removedNullByte.replace(constants_1.htmlEntitiesRegex, function (match, dec) {
        return String.fromCharCode(dec);
    });
}
function isValidUrl(url) {
    return URL.canParse(url);
}
function decodeURI(uri) {
    try {
        return decodeURIComponent(uri);
    }
    catch (e) {
        // Ignoring error
        // It is possible that the URI contains a `%` not associated
        // with URI/URL-encoding.
        return uri;
    }
}
function sanitizeUrl(url) {
    if (!url) {
        return constants_1.BLANK_URL;
    }
    var charsToDecode;
    var decodedUrl = decodeURI(url.trim());
    do {
        decodedUrl = decodeHtmlCharacters(decodedUrl)
            .replace(constants_1.htmlCtrlEntityRegex, "")
            .replace(constants_1.ctrlCharactersRegex, "")
            .replace(constants_1.whitespaceEscapeCharsRegex, "")
            .trim();
        decodedUrl = decodeURI(decodedUrl);
        charsToDecode =
            decodedUrl.match(constants_1.ctrlCharactersRegex) ||
                decodedUrl.match(constants_1.htmlEntitiesRegex) ||
                decodedUrl.match(constants_1.htmlCtrlEntityRegex) ||
                decodedUrl.match(constants_1.whitespaceEscapeCharsRegex);
    } while (charsToDecode && charsToDecode.length > 0);
    var sanitizedUrl = decodedUrl;
    if (!sanitizedUrl) {
        return constants_1.BLANK_URL;
    }
    if (isRelativeUrlWithoutProtocol(sanitizedUrl)) {
        return sanitizedUrl;
    }
    // Remove any leading whitespace before checking the URL scheme
    var trimmedUrl = sanitizedUrl.trimStart();
    var urlSchemeParseResults = trimmedUrl.match(constants_1.urlSchemeRegex);
    if (!urlSchemeParseResults) {
        return sanitizedUrl;
    }
    var urlScheme = urlSchemeParseResults[0].toLowerCase().trim();
    if (constants_1.invalidProtocolRegex.test(urlScheme)) {
        return constants_1.BLANK_URL;
    }
    var backSanitized = trimmedUrl.replace(/\\/g, "/");
    // Handle special cases for mailto: and custom deep-link protocols
    if (urlScheme === "mailto:" || urlScheme.includes("://")) {
        return backSanitized;
    }
    // For http and https URLs, perform additional validation
    if (urlScheme === "http:" || urlScheme === "https:") {
        if (!isValidUrl(backSanitized)) {
            return constants_1.BLANK_URL;
        }
        var url_1 = new URL(backSanitized);
        url_1.protocol = url_1.protocol.toLowerCase();
        url_1.hostname = url_1.hostname.toLowerCase();
        return url_1.toString();
    }
    return backSanitized;
}
exports.N = sanitizeUrl;


/***/ }),

/***/ 4719:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (/* binding */ purify)
/* harmony export */ });
/*! @license DOMPurify 3.2.6 | (c) Cure53 and other contributors | Released under the Apache license 2.0 and Mozilla Public License 2.0 | github.com/cure53/DOMPurify/blob/3.2.6/LICENSE */

const {
  entries,
  setPrototypeOf,
  isFrozen,
  getPrototypeOf,
  getOwnPropertyDescriptor
} = Object;
let {
  freeze,
  seal,
  create
} = Object; // eslint-disable-line import/no-mutable-exports
let {
  apply,
  construct
} = typeof Reflect !== 'undefined' && Reflect;
if (!freeze) {
  freeze = function freeze(x) {
    return x;
  };
}
if (!seal) {
  seal = function seal(x) {
    return x;
  };
}
if (!apply) {
  apply = function apply(fun, thisValue, args) {
    return fun.apply(thisValue, args);
  };
}
if (!construct) {
  construct = function construct(Func, args) {
    return new Func(...args);
  };
}
const arrayForEach = unapply(Array.prototype.forEach);
const arrayLastIndexOf = unapply(Array.prototype.lastIndexOf);
const arrayPop = unapply(Array.prototype.pop);
const arrayPush = unapply(Array.prototype.push);
const arraySplice = unapply(Array.prototype.splice);
const stringToLowerCase = unapply(String.prototype.toLowerCase);
const stringToString = unapply(String.prototype.toString);
const stringMatch = unapply(String.prototype.match);
const stringReplace = unapply(String.prototype.replace);
const stringIndexOf = unapply(String.prototype.indexOf);
const stringTrim = unapply(String.prototype.trim);
const objectHasOwnProperty = unapply(Object.prototype.hasOwnProperty);
const regExpTest = unapply(RegExp.prototype.test);
const typeErrorCreate = unconstruct(TypeError);
/**
 * Creates a new function that calls the given function with a specified thisArg and arguments.
 *
 * @param func - The function to be wrapped and called.
 * @returns A new function that calls the given function with a specified thisArg and arguments.
 */
function unapply(func) {
  return function (thisArg) {
    if (thisArg instanceof RegExp) {
      thisArg.lastIndex = 0;
    }
    for (var _len = arguments.length, args = new Array(_len > 1 ? _len - 1 : 0), _key = 1; _key < _len; _key++) {
      args[_key - 1] = arguments[_key];
    }
    return apply(func, thisArg, args);
  };
}
/**
 * Creates a new function that constructs an instance of the given constructor function with the provided arguments.
 *
 * @param func - The constructor function to be wrapped and called.
 * @returns A new function that constructs an instance of the given constructor function with the provided arguments.
 */
function unconstruct(func) {
  return function () {
    for (var _len2 = arguments.length, args = new Array(_len2), _key2 = 0; _key2 < _len2; _key2++) {
      args[_key2] = arguments[_key2];
    }
    return construct(func, args);
  };
}
/**
 * Add properties to a lookup table
 *
 * @param set - The set to which elements will be added.
 * @param array - The array containing elements to be added to the set.
 * @param transformCaseFunc - An optional function to transform the case of each element before adding to the set.
 * @returns The modified set with added elements.
 */
function addToSet(set, array) {
  let transformCaseFunc = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : stringToLowerCase;
  if (setPrototypeOf) {
    // Make 'in' and truthy checks like Boolean(set.constructor)
    // independent of any properties defined on Object.prototype.
    // Prevent prototype setters from intercepting set as a this value.
    setPrototypeOf(set, null);
  }
  let l = array.length;
  while (l--) {
    let element = array[l];
    if (typeof element === 'string') {
      const lcElement = transformCaseFunc(element);
      if (lcElement !== element) {
        // Config presets (e.g. tags.js, attrs.js) are immutable.
        if (!isFrozen(array)) {
          array[l] = lcElement;
        }
        element = lcElement;
      }
    }
    set[element] = true;
  }
  return set;
}
/**
 * Clean up an array to harden against CSPP
 *
 * @param array - The array to be cleaned.
 * @returns The cleaned version of the array
 */
function cleanArray(array) {
  for (let index = 0; index < array.length; index++) {
    const isPropertyExist = objectHasOwnProperty(array, index);
    if (!isPropertyExist) {
      array[index] = null;
    }
  }
  return array;
}
/**
 * Shallow clone an object
 *
 * @param object - The object to be cloned.
 * @returns A new object that copies the original.
 */
function clone(object) {
  const newObject = create(null);
  for (const [property, value] of entries(object)) {
    const isPropertyExist = objectHasOwnProperty(object, property);
    if (isPropertyExist) {
      if (Array.isArray(value)) {
        newObject[property] = cleanArray(value);
      } else if (value && typeof value === 'object' && value.constructor === Object) {
        newObject[property] = clone(value);
      } else {
        newObject[property] = value;
      }
    }
  }
  return newObject;
}
/**
 * This method automatically checks if the prop is function or getter and behaves accordingly.
 *
 * @param object - The object to look up the getter function in its prototype chain.
 * @param prop - The property name for which to find the getter function.
 * @returns The getter function found in the prototype chain or a fallback function.
 */
function lookupGetter(object, prop) {
  while (object !== null) {
    const desc = getOwnPropertyDescriptor(object, prop);
    if (desc) {
      if (desc.get) {
        return unapply(desc.get);
      }
      if (typeof desc.value === 'function') {
        return unapply(desc.value);
      }
    }
    object = getPrototypeOf(object);
  }
  function fallbackValue() {
    return null;
  }
  return fallbackValue;
}

const html$1 = freeze(['a', 'abbr', 'acronym', 'address', 'area', 'article', 'aside', 'audio', 'b', 'bdi', 'bdo', 'big', 'blink', 'blockquote', 'body', 'br', 'button', 'canvas', 'caption', 'center', 'cite', 'code', 'col', 'colgroup', 'content', 'data', 'datalist', 'dd', 'decorator', 'del', 'details', 'dfn', 'dialog', 'dir', 'div', 'dl', 'dt', 'element', 'em', 'fieldset', 'figcaption', 'figure', 'font', 'footer', 'form', 'h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'head', 'header', 'hgroup', 'hr', 'html', 'i', 'img', 'input', 'ins', 'kbd', 'label', 'legend', 'li', 'main', 'map', 'mark', 'marquee', 'menu', 'menuitem', 'meter', 'nav', 'nobr', 'ol', 'optgroup', 'option', 'output', 'p', 'picture', 'pre', 'progress', 'q', 'rp', 'rt', 'ruby', 's', 'samp', 'section', 'select', 'shadow', 'small', 'source', 'spacer', 'span', 'strike', 'strong', 'style', 'sub', 'summary', 'sup', 'table', 'tbody', 'td', 'template', 'textarea', 'tfoot', 'th', 'thead', 'time', 'tr', 'track', 'tt', 'u', 'ul', 'var', 'video', 'wbr']);
const svg$1 = freeze(['svg', 'a', 'altglyph', 'altglyphdef', 'altglyphitem', 'animatecolor', 'animatemotion', 'animatetransform', 'circle', 'clippath', 'defs', 'desc', 'ellipse', 'filter', 'font', 'g', 'glyph', 'glyphref', 'hkern', 'image', 'line', 'lineargradient', 'marker', 'mask', 'metadata', 'mpath', 'path', 'pattern', 'polygon', 'polyline', 'radialgradient', 'rect', 'stop', 'style', 'switch', 'symbol', 'text', 'textpath', 'title', 'tref', 'tspan', 'view', 'vkern']);
const svgFilters = freeze(['feBlend', 'feColorMatrix', 'feComponentTransfer', 'feComposite', 'feConvolveMatrix', 'feDiffuseLighting', 'feDisplacementMap', 'feDistantLight', 'feDropShadow', 'feFlood', 'feFuncA', 'feFuncB', 'feFuncG', 'feFuncR', 'feGaussianBlur', 'feImage', 'feMerge', 'feMergeNode', 'feMorphology', 'feOffset', 'fePointLight', 'feSpecularLighting', 'feSpotLight', 'feTile', 'feTurbulence']);
// List of SVG elements that are disallowed by default.
// We still need to know them so that we can do namespace
// checks properly in case one wants to add them to
// allow-list.
const svgDisallowed = freeze(['animate', 'color-profile', 'cursor', 'discard', 'font-face', 'font-face-format', 'font-face-name', 'font-face-src', 'font-face-uri', 'foreignobject', 'hatch', 'hatchpath', 'mesh', 'meshgradient', 'meshpatch', 'meshrow', 'missing-glyph', 'script', 'set', 'solidcolor', 'unknown', 'use']);
const mathMl$1 = freeze(['math', 'menclose', 'merror', 'mfenced', 'mfrac', 'mglyph', 'mi', 'mlabeledtr', 'mmultiscripts', 'mn', 'mo', 'mover', 'mpadded', 'mphantom', 'mroot', 'mrow', 'ms', 'mspace', 'msqrt', 'mstyle', 'msub', 'msup', 'msubsup', 'mtable', 'mtd', 'mtext', 'mtr', 'munder', 'munderover', 'mprescripts']);
// Similarly to SVG, we want to know all MathML elements,
// even those that we disallow by default.
const mathMlDisallowed = freeze(['maction', 'maligngroup', 'malignmark', 'mlongdiv', 'mscarries', 'mscarry', 'msgroup', 'mstack', 'msline', 'msrow', 'semantics', 'annotation', 'annotation-xml', 'mprescripts', 'none']);
const text = freeze(['#text']);

const html = freeze(['accept', 'action', 'align', 'alt', 'autocapitalize', 'autocomplete', 'autopictureinpicture', 'autoplay', 'background', 'bgcolor', 'border', 'capture', 'cellpadding', 'cellspacing', 'checked', 'cite', 'class', 'clear', 'color', 'cols', 'colspan', 'controls', 'controlslist', 'coords', 'crossorigin', 'datetime', 'decoding', 'default', 'dir', 'disabled', 'disablepictureinpicture', 'disableremoteplayback', 'download', 'draggable', 'enctype', 'enterkeyhint', 'face', 'for', 'headers', 'height', 'hidden', 'high', 'href', 'hreflang', 'id', 'inputmode', 'integrity', 'ismap', 'kind', 'label', 'lang', 'list', 'loading', 'loop', 'low', 'max', 'maxlength', 'media', 'method', 'min', 'minlength', 'multiple', 'muted', 'name', 'nonce', 'noshade', 'novalidate', 'nowrap', 'open', 'optimum', 'pattern', 'placeholder', 'playsinline', 'popover', 'popovertarget', 'popovertargetaction', 'poster', 'preload', 'pubdate', 'radiogroup', 'readonly', 'rel', 'required', 'rev', 'reversed', 'role', 'rows', 'rowspan', 'spellcheck', 'scope', 'selected', 'shape', 'size', 'sizes', 'span', 'srclang', 'start', 'src', 'srcset', 'step', 'style', 'summary', 'tabindex', 'title', 'translate', 'type', 'usemap', 'valign', 'value', 'width', 'wrap', 'xmlns', 'slot']);
const svg = freeze(['accent-height', 'accumulate', 'additive', 'alignment-baseline', 'amplitude', 'ascent', 'attributename', 'attributetype', 'azimuth', 'basefrequency', 'baseline-shift', 'begin', 'bias', 'by', 'class', 'clip', 'clippathunits', 'clip-path', 'clip-rule', 'color', 'color-interpolation', 'color-interpolation-filters', 'color-profile', 'color-rendering', 'cx', 'cy', 'd', 'dx', 'dy', 'diffuseconstant', 'direction', 'display', 'divisor', 'dur', 'edgemode', 'elevation', 'end', 'exponent', 'fill', 'fill-opacity', 'fill-rule', 'filter', 'filterunits', 'flood-color', 'flood-opacity', 'font-family', 'font-size', 'font-size-adjust', 'font-stretch', 'font-style', 'font-variant', 'font-weight', 'fx', 'fy', 'g1', 'g2', 'glyph-name', 'glyphref', 'gradientunits', 'gradienttransform', 'height', 'href', 'id', 'image-rendering', 'in', 'in2', 'intercept', 'k', 'k1', 'k2', 'k3', 'k4', 'kerning', 'keypoints', 'keysplines', 'keytimes', 'lang', 'lengthadjust', 'letter-spacing', 'kernelmatrix', 'kernelunitlength', 'lighting-color', 'local', 'marker-end', 'marker-mid', 'marker-start', 'markerheight', 'markerunits', 'markerwidth', 'maskcontentunits', 'maskunits', 'max', 'mask', 'media', 'method', 'mode', 'min', 'name', 'numoctaves', 'offset', 'operator', 'opacity', 'order', 'orient', 'orientation', 'origin', 'overflow', 'paint-order', 'path', 'pathlength', 'patterncontentunits', 'patterntransform', 'patternunits', 'points', 'preservealpha', 'preserveaspectratio', 'primitiveunits', 'r', 'rx', 'ry', 'radius', 'refx', 'refy', 'repeatcount', 'repeatdur', 'restart', 'result', 'rotate', 'scale', 'seed', 'shape-rendering', 'slope', 'specularconstant', 'specularexponent', 'spreadmethod', 'startoffset', 'stddeviation', 'stitchtiles', 'stop-color', 'stop-opacity', 'stroke-dasharray', 'stroke-dashoffset', 'stroke-linecap', 'stroke-linejoin', 'stroke-miterlimit', 'stroke-opacity', 'stroke', 'stroke-width', 'style', 'surfacescale', 'systemlanguage', 'tabindex', 'tablevalues', 'targetx', 'targety', 'transform', 'transform-origin', 'text-anchor', 'text-decoration', 'text-rendering', 'textlength', 'type', 'u1', 'u2', 'unicode', 'values', 'viewbox', 'visibility', 'version', 'vert-adv-y', 'vert-origin-x', 'vert-origin-y', 'width', 'word-spacing', 'wrap', 'writing-mode', 'xchannelselector', 'ychannelselector', 'x', 'x1', 'x2', 'xmlns', 'y', 'y1', 'y2', 'z', 'zoomandpan']);
const mathMl = freeze(['accent', 'accentunder', 'align', 'bevelled', 'close', 'columnsalign', 'columnlines', 'columnspan', 'denomalign', 'depth', 'dir', 'display', 'displaystyle', 'encoding', 'fence', 'frame', 'height', 'href', 'id', 'largeop', 'length', 'linethickness', 'lspace', 'lquote', 'mathbackground', 'mathcolor', 'mathsize', 'mathvariant', 'maxsize', 'minsize', 'movablelimits', 'notation', 'numalign', 'open', 'rowalign', 'rowlines', 'rowspacing', 'rowspan', 'rspace', 'rquote', 'scriptlevel', 'scriptminsize', 'scriptsizemultiplier', 'selection', 'separator', 'separators', 'stretchy', 'subscriptshift', 'supscriptshift', 'symmetric', 'voffset', 'width', 'xmlns']);
const xml = freeze(['xlink:href', 'xml:id', 'xlink:title', 'xml:space', 'xmlns:xlink']);

// eslint-disable-next-line unicorn/better-regex
const MUSTACHE_EXPR = seal(/\{\{[\w\W]*|[\w\W]*\}\}/gm); // Specify template detection regex for SAFE_FOR_TEMPLATES mode
const ERB_EXPR = seal(/<%[\w\W]*|[\w\W]*%>/gm);
const TMPLIT_EXPR = seal(/\$\{[\w\W]*/gm); // eslint-disable-line unicorn/better-regex
const DATA_ATTR = seal(/^data-[\-\w.\u00B7-\uFFFF]+$/); // eslint-disable-line no-useless-escape
const ARIA_ATTR = seal(/^aria-[\-\w]+$/); // eslint-disable-line no-useless-escape
const IS_ALLOWED_URI = seal(/^(?:(?:(?:f|ht)tps?|mailto|tel|callto|sms|cid|xmpp|matrix):|[^a-z]|[a-z+.\-]+(?:[^a-z+.\-:]|$))/i // eslint-disable-line no-useless-escape
);
const IS_SCRIPT_OR_DATA = seal(/^(?:\w+script|data):/i);
const ATTR_WHITESPACE = seal(/[\u0000-\u0020\u00A0\u1680\u180E\u2000-\u2029\u205F\u3000]/g // eslint-disable-line no-control-regex
);
const DOCTYPE_NAME = seal(/^html$/i);
const CUSTOM_ELEMENT = seal(/^[a-z][.\w]*(-[.\w]+)+$/i);

var EXPRESSIONS = /*#__PURE__*/Object.freeze({
  __proto__: null,
  ARIA_ATTR: ARIA_ATTR,
  ATTR_WHITESPACE: ATTR_WHITESPACE,
  CUSTOM_ELEMENT: CUSTOM_ELEMENT,
  DATA_ATTR: DATA_ATTR,
  DOCTYPE_NAME: DOCTYPE_NAME,
  ERB_EXPR: ERB_EXPR,
  IS_ALLOWED_URI: IS_ALLOWED_URI,
  IS_SCRIPT_OR_DATA: IS_SCRIPT_OR_DATA,
  MUSTACHE_EXPR: MUSTACHE_EXPR,
  TMPLIT_EXPR: TMPLIT_EXPR
});

/* eslint-disable @typescript-eslint/indent */
// https://developer.mozilla.org/en-US/docs/Web/API/Node/nodeType
const NODE_TYPE = {
  element: 1,
  attribute: 2,
  text: 3,
  cdataSection: 4,
  entityReference: 5,
  // Deprecated
  entityNode: 6,
  // Deprecated
  progressingInstruction: 7,
  comment: 8,
  document: 9,
  documentType: 10,
  documentFragment: 11,
  notation: 12 // Deprecated
};
const getGlobal = function getGlobal() {
  return typeof window === 'undefined' ? null : window;
};
/**
 * Creates a no-op policy for internal use only.
 * Don't export this function outside this module!
 * @param trustedTypes The policy factory.
 * @param purifyHostElement The Script element used to load DOMPurify (to determine policy name suffix).
 * @return The policy created (or null, if Trusted Types
 * are not supported or creating the policy failed).
 */
const _createTrustedTypesPolicy = function _createTrustedTypesPolicy(trustedTypes, purifyHostElement) {
  if (typeof trustedTypes !== 'object' || typeof trustedTypes.createPolicy !== 'function') {
    return null;
  }
  // Allow the callers to control the unique policy name
  // by adding a data-tt-policy-suffix to the script element with the DOMPurify.
  // Policy creation with duplicate names throws in Trusted Types.
  let suffix = null;
  const ATTR_NAME = 'data-tt-policy-suffix';
  if (purifyHostElement && purifyHostElement.hasAttribute(ATTR_NAME)) {
    suffix = purifyHostElement.getAttribute(ATTR_NAME);
  }
  const policyName = 'dompurify' + (suffix ? '#' + suffix : '');
  try {
    return trustedTypes.createPolicy(policyName, {
      createHTML(html) {
        return html;
      },
      createScriptURL(scriptUrl) {
        return scriptUrl;
      }
    });
  } catch (_) {
    // Policy creation failed (most likely another DOMPurify script has
    // already run). Skip creating the policy, as this will only cause errors
    // if TT are enforced.
    console.warn('TrustedTypes policy ' + policyName + ' could not be created.');
    return null;
  }
};
const _createHooksMap = function _createHooksMap() {
  return {
    afterSanitizeAttributes: [],
    afterSanitizeElements: [],
    afterSanitizeShadowDOM: [],
    beforeSanitizeAttributes: [],
    beforeSanitizeElements: [],
    beforeSanitizeShadowDOM: [],
    uponSanitizeAttribute: [],
    uponSanitizeElement: [],
    uponSanitizeShadowNode: []
  };
};
function createDOMPurify() {
  let window = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : getGlobal();
  const DOMPurify = root => createDOMPurify(root);
  DOMPurify.version = '3.2.6';
  DOMPurify.removed = [];
  if (!window || !window.document || window.document.nodeType !== NODE_TYPE.document || !window.Element) {
    // Not running in a browser, provide a factory function
    // so that you can pass your own Window
    DOMPurify.isSupported = false;
    return DOMPurify;
  }
  let {
    document
  } = window;
  const originalDocument = document;
  const currentScript = originalDocument.currentScript;
  const {
    DocumentFragment,
    HTMLTemplateElement,
    Node,
    Element,
    NodeFilter,
    NamedNodeMap = window.NamedNodeMap || window.MozNamedAttrMap,
    HTMLFormElement,
    DOMParser,
    trustedTypes
  } = window;
  const ElementPrototype = Element.prototype;
  const cloneNode = lookupGetter(ElementPrototype, 'cloneNode');
  const remove = lookupGetter(ElementPrototype, 'remove');
  const getNextSibling = lookupGetter(ElementPrototype, 'nextSibling');
  const getChildNodes = lookupGetter(ElementPrototype, 'childNodes');
  const getParentNode = lookupGetter(ElementPrototype, 'parentNode');
  // As per issue #47, the web-components registry is inherited by a
  // new document created via createHTMLDocument. As per the spec
  // (http://w3c.github.io/webcomponents/spec/custom/#creating-and-passing-registries)
  // a new empty registry is used when creating a template contents owner
  // document, so we use that as our parent document to ensure nothing
  // is inherited.
  if (typeof HTMLTemplateElement === 'function') {
    const template = document.createElement('template');
    if (template.content && template.content.ownerDocument) {
      document = template.content.ownerDocument;
    }
  }
  let trustedTypesPolicy;
  let emptyHTML = '';
  const {
    implementation,
    createNodeIterator,
    createDocumentFragment,
    getElementsByTagName
  } = document;
  const {
    importNode
  } = originalDocument;
  let hooks = _createHooksMap();
  /**
   * Expose whether this browser supports running the full DOMPurify.
   */
  DOMPurify.isSupported = typeof entries === 'function' && typeof getParentNode === 'function' && implementation && implementation.createHTMLDocument !== undefined;
  const {
    MUSTACHE_EXPR,
    ERB_EXPR,
    TMPLIT_EXPR,
    DATA_ATTR,
    ARIA_ATTR,
    IS_SCRIPT_OR_DATA,
    ATTR_WHITESPACE,
    CUSTOM_ELEMENT
  } = EXPRESSIONS;
  let {
    IS_ALLOWED_URI: IS_ALLOWED_URI$1
  } = EXPRESSIONS;
  /**
   * We consider the elements and attributes below to be safe. Ideally
   * don't add any new ones but feel free to remove unwanted ones.
   */
  /* allowed element names */
  let ALLOWED_TAGS = null;
  const DEFAULT_ALLOWED_TAGS = addToSet({}, [...html$1, ...svg$1, ...svgFilters, ...mathMl$1, ...text]);
  /* Allowed attribute names */
  let ALLOWED_ATTR = null;
  const DEFAULT_ALLOWED_ATTR = addToSet({}, [...html, ...svg, ...mathMl, ...xml]);
  /*
   * Configure how DOMPurify should handle custom elements and their attributes as well as customized built-in elements.
   * @property {RegExp|Function|null} tagNameCheck one of [null, regexPattern, predicate]. Default: `null` (disallow any custom elements)
   * @property {RegExp|Function|null} attributeNameCheck one of [null, regexPattern, predicate]. Default: `null` (disallow any attributes not on the allow list)
   * @property {boolean} allowCustomizedBuiltInElements allow custom elements derived from built-ins if they pass CUSTOM_ELEMENT_HANDLING.tagNameCheck. Default: `false`.
   */
  let CUSTOM_ELEMENT_HANDLING = Object.seal(create(null, {
    tagNameCheck: {
      writable: true,
      configurable: false,
      enumerable: true,
      value: null
    },
    attributeNameCheck: {
      writable: true,
      configurable: false,
      enumerable: true,
      value: null
    },
    allowCustomizedBuiltInElements: {
      writable: true,
      configurable: false,
      enumerable: true,
      value: false
    }
  }));
  /* Explicitly forbidden tags (overrides ALLOWED_TAGS/ADD_TAGS) */
  let FORBID_TAGS = null;
  /* Explicitly forbidden attributes (overrides ALLOWED_ATTR/ADD_ATTR) */
  let FORBID_ATTR = null;
  /* Decide if ARIA attributes are okay */
  let ALLOW_ARIA_ATTR = true;
  /* Decide if custom data attributes are okay */
  let ALLOW_DATA_ATTR = true;
  /* Decide if unknown protocols are okay */
  let ALLOW_UNKNOWN_PROTOCOLS = false;
  /* Decide if self-closing tags in attributes are allowed.
   * Usually removed due to a mXSS issue in jQuery 3.0 */
  let ALLOW_SELF_CLOSE_IN_ATTR = true;
  /* Output should be safe for common template engines.
   * This means, DOMPurify removes data attributes, mustaches and ERB
   */
  let SAFE_FOR_TEMPLATES = false;
  /* Output should be safe even for XML used within HTML and alike.
   * This means, DOMPurify removes comments when containing risky content.
   */
  let SAFE_FOR_XML = true;
  /* Decide if document with <html>... should be returned */
  let WHOLE_DOCUMENT = false;
  /* Track whether config is already set on this instance of DOMPurify. */
  let SET_CONFIG = false;
  /* Decide if all elements (e.g. style, script) must be children of
   * document.body. By default, browsers might move them to document.head */
  let FORCE_BODY = false;
  /* Decide if a DOM `HTMLBodyElement` should be returned, instead of a html
   * string (or a TrustedHTML object if Trusted Types are supported).
   * If `WHOLE_DOCUMENT` is enabled a `HTMLHtmlElement` will be returned instead
   */
  let RETURN_DOM = false;
  /* Decide if a DOM `DocumentFragment` should be returned, instead of a html
   * string  (or a TrustedHTML object if Trusted Types are supported) */
  let RETURN_DOM_FRAGMENT = false;
  /* Try to return a Trusted Type object instead of a string, return a string in
   * case Trusted Types are not supported  */
  let RETURN_TRUSTED_TYPE = false;
  /* Output should be free from DOM clobbering attacks?
   * This sanitizes markups named with colliding, clobberable built-in DOM APIs.
   */
  let SANITIZE_DOM = true;
  /* Achieve full DOM Clobbering protection by isolating the namespace of named
   * properties and JS variables, mitigating attacks that abuse the HTML/DOM spec rules.
   *
   * HTML/DOM spec rules that enable DOM Clobbering:
   *   - Named Access on Window (§7.3.3)
   *   - DOM Tree Accessors (§3.1.5)
   *   - Form Element Parent-Child Relations (§4.10.3)
   *   - Iframe srcdoc / Nested WindowProxies (§4.8.5)
   *   - HTMLCollection (§4.2.10.2)
   *
   * Namespace isolation is implemented by prefixing `id` and `name` attributes
   * with a constant string, i.e., `user-content-`
   */
  let SANITIZE_NAMED_PROPS = false;
  const SANITIZE_NAMED_PROPS_PREFIX = 'user-content-';
  /* Keep element content when removing element? */
  let KEEP_CONTENT = true;
  /* If a `Node` is passed to sanitize(), then performs sanitization in-place instead
   * of importing it into a new Document and returning a sanitized copy */
  let IN_PLACE = false;
  /* Allow usage of profiles like html, svg and mathMl */
  let USE_PROFILES = {};
  /* Tags to ignore content of when KEEP_CONTENT is true */
  let FORBID_CONTENTS = null;
  const DEFAULT_FORBID_CONTENTS = addToSet({}, ['annotation-xml', 'audio', 'colgroup', 'desc', 'foreignobject', 'head', 'iframe', 'math', 'mi', 'mn', 'mo', 'ms', 'mtext', 'noembed', 'noframes', 'noscript', 'plaintext', 'script', 'style', 'svg', 'template', 'thead', 'title', 'video', 'xmp']);
  /* Tags that are safe for data: URIs */
  let DATA_URI_TAGS = null;
  const DEFAULT_DATA_URI_TAGS = addToSet({}, ['audio', 'video', 'img', 'source', 'image', 'track']);
  /* Attributes safe for values like "javascript:" */
  let URI_SAFE_ATTRIBUTES = null;
  const DEFAULT_URI_SAFE_ATTRIBUTES = addToSet({}, ['alt', 'class', 'for', 'id', 'label', 'name', 'pattern', 'placeholder', 'role', 'summary', 'title', 'value', 'style', 'xmlns']);
  const MATHML_NAMESPACE = 'http://www.w3.org/1998/Math/MathML';
  const SVG_NAMESPACE = 'http://www.w3.org/2000/svg';
  const HTML_NAMESPACE = 'http://www.w3.org/1999/xhtml';
  /* Document namespace */
  let NAMESPACE = HTML_NAMESPACE;
  let IS_EMPTY_INPUT = false;
  /* Allowed XHTML+XML namespaces */
  let ALLOWED_NAMESPACES = null;
  const DEFAULT_ALLOWED_NAMESPACES = addToSet({}, [MATHML_NAMESPACE, SVG_NAMESPACE, HTML_NAMESPACE], stringToString);
  let MATHML_TEXT_INTEGRATION_POINTS = addToSet({}, ['mi', 'mo', 'mn', 'ms', 'mtext']);
  let HTML_INTEGRATION_POINTS = addToSet({}, ['annotation-xml']);
  // Certain elements are allowed in both SVG and HTML
  // namespace. We need to specify them explicitly
  // so that they don't get erroneously deleted from
  // HTML namespace.
  const COMMON_SVG_AND_HTML_ELEMENTS = addToSet({}, ['title', 'style', 'font', 'a', 'script']);
  /* Parsing of strict XHTML documents */
  let PARSER_MEDIA_TYPE = null;
  const SUPPORTED_PARSER_MEDIA_TYPES = ['application/xhtml+xml', 'text/html'];
  const DEFAULT_PARSER_MEDIA_TYPE = 'text/html';
  let transformCaseFunc = null;
  /* Keep a reference to config to pass to hooks */
  let CONFIG = null;
  /* Ideally, do not touch anything below this line */
  /* ______________________________________________ */
  const formElement = document.createElement('form');
  const isRegexOrFunction = function isRegexOrFunction(testValue) {
    return testValue instanceof RegExp || testValue instanceof Function;
  };
  /**
   * _parseConfig
   *
   * @param cfg optional config literal
   */
  // eslint-disable-next-line complexity
  const _parseConfig = function _parseConfig() {
    let cfg = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {};
    if (CONFIG && CONFIG === cfg) {
      return;
    }
    /* Shield configuration object from tampering */
    if (!cfg || typeof cfg !== 'object') {
      cfg = {};
    }
    /* Shield configuration object from prototype pollution */
    cfg = clone(cfg);
    PARSER_MEDIA_TYPE =
    // eslint-disable-next-line unicorn/prefer-includes
    SUPPORTED_PARSER_MEDIA_TYPES.indexOf(cfg.PARSER_MEDIA_TYPE) === -1 ? DEFAULT_PARSER_MEDIA_TYPE : cfg.PARSER_MEDIA_TYPE;
    // HTML tags and attributes are not case-sensitive, converting to lowercase. Keeping XHTML as is.
    transformCaseFunc = PARSER_MEDIA_TYPE === 'application/xhtml+xml' ? stringToString : stringToLowerCase;
    /* Set configuration parameters */
    ALLOWED_TAGS = objectHasOwnProperty(cfg, 'ALLOWED_TAGS') ? addToSet({}, cfg.ALLOWED_TAGS, transformCaseFunc) : DEFAULT_ALLOWED_TAGS;
    ALLOWED_ATTR = objectHasOwnProperty(cfg, 'ALLOWED_ATTR') ? addToSet({}, cfg.ALLOWED_ATTR, transformCaseFunc) : DEFAULT_ALLOWED_ATTR;
    ALLOWED_NAMESPACES = objectHasOwnProperty(cfg, 'ALLOWED_NAMESPACES') ? addToSet({}, cfg.ALLOWED_NAMESPACES, stringToString) : DEFAULT_ALLOWED_NAMESPACES;
    URI_SAFE_ATTRIBUTES = objectHasOwnProperty(cfg, 'ADD_URI_SAFE_ATTR') ? addToSet(clone(DEFAULT_URI_SAFE_ATTRIBUTES), cfg.ADD_URI_SAFE_ATTR, transformCaseFunc) : DEFAULT_URI_SAFE_ATTRIBUTES;
    DATA_URI_TAGS = objectHasOwnProperty(cfg, 'ADD_DATA_URI_TAGS') ? addToSet(clone(DEFAULT_DATA_URI_TAGS), cfg.ADD_DATA_URI_TAGS, transformCaseFunc) : DEFAULT_DATA_URI_TAGS;
    FORBID_CONTENTS = objectHasOwnProperty(cfg, 'FORBID_CONTENTS') ? addToSet({}, cfg.FORBID_CONTENTS, transformCaseFunc) : DEFAULT_FORBID_CONTENTS;
    FORBID_TAGS = objectHasOwnProperty(cfg, 'FORBID_TAGS') ? addToSet({}, cfg.FORBID_TAGS, transformCaseFunc) : clone({});
    FORBID_ATTR = objectHasOwnProperty(cfg, 'FORBID_ATTR') ? addToSet({}, cfg.FORBID_ATTR, transformCaseFunc) : clone({});
    USE_PROFILES = objectHasOwnProperty(cfg, 'USE_PROFILES') ? cfg.USE_PROFILES : false;
    ALLOW_ARIA_ATTR = cfg.ALLOW_ARIA_ATTR !== false; // Default true
    ALLOW_DATA_ATTR = cfg.ALLOW_DATA_ATTR !== false; // Default true
    ALLOW_UNKNOWN_PROTOCOLS = cfg.ALLOW_UNKNOWN_PROTOCOLS || false; // Default false
    ALLOW_SELF_CLOSE_IN_ATTR = cfg.ALLOW_SELF_CLOSE_IN_ATTR !== false; // Default true
    SAFE_FOR_TEMPLATES = cfg.SAFE_FOR_TEMPLATES || false; // Default false
    SAFE_FOR_XML = cfg.SAFE_FOR_XML !== false; // Default true
    WHOLE_DOCUMENT = cfg.WHOLE_DOCUMENT || false; // Default false
    RETURN_DOM = cfg.RETURN_DOM || false; // Default false
    RETURN_DOM_FRAGMENT = cfg.RETURN_DOM_FRAGMENT || false; // Default false
    RETURN_TRUSTED_TYPE = cfg.RETURN_TRUSTED_TYPE || false; // Default false
    FORCE_BODY = cfg.FORCE_BODY || false; // Default false
    SANITIZE_DOM = cfg.SANITIZE_DOM !== false; // Default true
    SANITIZE_NAMED_PROPS = cfg.SANITIZE_NAMED_PROPS || false; // Default false
    KEEP_CONTENT = cfg.KEEP_CONTENT !== false; // Default true
    IN_PLACE = cfg.IN_PLACE || false; // Default false
    IS_ALLOWED_URI$1 = cfg.ALLOWED_URI_REGEXP || IS_ALLOWED_URI;
    NAMESPACE = cfg.NAMESPACE || HTML_NAMESPACE;
    MATHML_TEXT_INTEGRATION_POINTS = cfg.MATHML_TEXT_INTEGRATION_POINTS || MATHML_TEXT_INTEGRATION_POINTS;
    HTML_INTEGRATION_POINTS = cfg.HTML_INTEGRATION_POINTS || HTML_INTEGRATION_POINTS;
    CUSTOM_ELEMENT_HANDLING = cfg.CUSTOM_ELEMENT_HANDLING || {};
    if (cfg.CUSTOM_ELEMENT_HANDLING && isRegexOrFunction(cfg.CUSTOM_ELEMENT_HANDLING.tagNameCheck)) {
      CUSTOM_ELEMENT_HANDLING.tagNameCheck = cfg.CUSTOM_ELEMENT_HANDLING.tagNameCheck;
    }
    if (cfg.CUSTOM_ELEMENT_HANDLING && isRegexOrFunction(cfg.CUSTOM_ELEMENT_HANDLING.attributeNameCheck)) {
      CUSTOM_ELEMENT_HANDLING.attributeNameCheck = cfg.CUSTOM_ELEMENT_HANDLING.attributeNameCheck;
    }
    if (cfg.CUSTOM_ELEMENT_HANDLING && typeof cfg.CUSTOM_ELEMENT_HANDLING.allowCustomizedBuiltInElements === 'boolean') {
      CUSTOM_ELEMENT_HANDLING.allowCustomizedBuiltInElements = cfg.CUSTOM_ELEMENT_HANDLING.allowCustomizedBuiltInElements;
    }
    if (SAFE_FOR_TEMPLATES) {
      ALLOW_DATA_ATTR = false;
    }
    if (RETURN_DOM_FRAGMENT) {
      RETURN_DOM = true;
    }
    /* Parse profile info */
    if (USE_PROFILES) {
      ALLOWED_TAGS = addToSet({}, text);
      ALLOWED_ATTR = [];
      if (USE_PROFILES.html === true) {
        addToSet(ALLOWED_TAGS, html$1);
        addToSet(ALLOWED_ATTR, html);
      }
      if (USE_PROFILES.svg === true) {
        addToSet(ALLOWED_TAGS, svg$1);
        addToSet(ALLOWED_ATTR, svg);
        addToSet(ALLOWED_ATTR, xml);
      }
      if (USE_PROFILES.svgFilters === true) {
        addToSet(ALLOWED_TAGS, svgFilters);
        addToSet(ALLOWED_ATTR, svg);
        addToSet(ALLOWED_ATTR, xml);
      }
      if (USE_PROFILES.mathMl === true) {
        addToSet(ALLOWED_TAGS, mathMl$1);
        addToSet(ALLOWED_ATTR, mathMl);
        addToSet(ALLOWED_ATTR, xml);
      }
    }
    /* Merge configuration parameters */
    if (cfg.ADD_TAGS) {
      if (ALLOWED_TAGS === DEFAULT_ALLOWED_TAGS) {
        ALLOWED_TAGS = clone(ALLOWED_TAGS);
      }
      addToSet(ALLOWED_TAGS, cfg.ADD_TAGS, transformCaseFunc);
    }
    if (cfg.ADD_ATTR) {
      if (ALLOWED_ATTR === DEFAULT_ALLOWED_ATTR) {
        ALLOWED_ATTR = clone(ALLOWED_ATTR);
      }
      addToSet(ALLOWED_ATTR, cfg.ADD_ATTR, transformCaseFunc);
    }
    if (cfg.ADD_URI_SAFE_ATTR) {
      addToSet(URI_SAFE_ATTRIBUTES, cfg.ADD_URI_SAFE_ATTR, transformCaseFunc);
    }
    if (cfg.FORBID_CONTENTS) {
      if (FORBID_CONTENTS === DEFAULT_FORBID_CONTENTS) {
        FORBID_CONTENTS = clone(FORBID_CONTENTS);
      }
      addToSet(FORBID_CONTENTS, cfg.FORBID_CONTENTS, transformCaseFunc);
    }
    /* Add #text in case KEEP_CONTENT is set to true */
    if (KEEP_CONTENT) {
      ALLOWED_TAGS['#text'] = true;
    }
    /* Add html, head and body to ALLOWED_TAGS in case WHOLE_DOCUMENT is true */
    if (WHOLE_DOCUMENT) {
      addToSet(ALLOWED_TAGS, ['html', 'head', 'body']);
    }
    /* Add tbody to ALLOWED_TAGS in case tables are permitted, see #286, #365 */
    if (ALLOWED_TAGS.table) {
      addToSet(ALLOWED_TAGS, ['tbody']);
      delete FORBID_TAGS.tbody;
    }
    if (cfg.TRUSTED_TYPES_POLICY) {
      if (typeof cfg.TRUSTED_TYPES_POLICY.createHTML !== 'function') {
        throw typeErrorCreate('TRUSTED_TYPES_POLICY configuration option must provide a "createHTML" hook.');
      }
      if (typeof cfg.TRUSTED_TYPES_POLICY.createScriptURL !== 'function') {
        throw typeErrorCreate('TRUSTED_TYPES_POLICY configuration option must provide a "createScriptURL" hook.');
      }
      // Overwrite existing TrustedTypes policy.
      trustedTypesPolicy = cfg.TRUSTED_TYPES_POLICY;
      // Sign local variables required by `sanitize`.
      emptyHTML = trustedTypesPolicy.createHTML('');
    } else {
      // Uninitialized policy, attempt to initialize the internal dompurify policy.
      if (trustedTypesPolicy === undefined) {
        trustedTypesPolicy = _createTrustedTypesPolicy(trustedTypes, currentScript);
      }
      // If creating the internal policy succeeded sign internal variables.
      if (trustedTypesPolicy !== null && typeof emptyHTML === 'string') {
        emptyHTML = trustedTypesPolicy.createHTML('');
      }
    }
    // Prevent further manipulation of configuration.
    // Not available in IE8, Safari 5, etc.
    if (freeze) {
      freeze(cfg);
    }
    CONFIG = cfg;
  };
  /* Keep track of all possible SVG and MathML tags
   * so that we can perform the namespace checks
   * correctly. */
  const ALL_SVG_TAGS = addToSet({}, [...svg$1, ...svgFilters, ...svgDisallowed]);
  const ALL_MATHML_TAGS = addToSet({}, [...mathMl$1, ...mathMlDisallowed]);
  /**
   * @param element a DOM element whose namespace is being checked
   * @returns Return false if the element has a
   *  namespace that a spec-compliant parser would never
   *  return. Return true otherwise.
   */
  const _checkValidNamespace = function _checkValidNamespace(element) {
    let parent = getParentNode(element);
    // In JSDOM, if we're inside shadow DOM, then parentNode
    // can be null. We just simulate parent in this case.
    if (!parent || !parent.tagName) {
      parent = {
        namespaceURI: NAMESPACE,
        tagName: 'template'
      };
    }
    const tagName = stringToLowerCase(element.tagName);
    const parentTagName = stringToLowerCase(parent.tagName);
    if (!ALLOWED_NAMESPACES[element.namespaceURI]) {
      return false;
    }
    if (element.namespaceURI === SVG_NAMESPACE) {
      // The only way to switch from HTML namespace to SVG
      // is via <svg>. If it happens via any other tag, then
      // it should be killed.
      if (parent.namespaceURI === HTML_NAMESPACE) {
        return tagName === 'svg';
      }
      // The only way to switch from MathML to SVG is via`
      // svg if parent is either <annotation-xml> or MathML
      // text integration points.
      if (parent.namespaceURI === MATHML_NAMESPACE) {
        return tagName === 'svg' && (parentTagName === 'annotation-xml' || MATHML_TEXT_INTEGRATION_POINTS[parentTagName]);
      }
      // We only allow elements that are defined in SVG
      // spec. All others are disallowed in SVG namespace.
      return Boolean(ALL_SVG_TAGS[tagName]);
    }
    if (element.namespaceURI === MATHML_NAMESPACE) {
      // The only way to switch from HTML namespace to MathML
      // is via <math>. If it happens via any other tag, then
      // it should be killed.
      if (parent.namespaceURI === HTML_NAMESPACE) {
        return tagName === 'math';
      }
      // The only way to switch from SVG to MathML is via
      // <math> and HTML integration points
      if (parent.namespaceURI === SVG_NAMESPACE) {
        return tagName === 'math' && HTML_INTEGRATION_POINTS[parentTagName];
      }
      // We only allow elements that are defined in MathML
      // spec. All others are disallowed in MathML namespace.
      return Boolean(ALL_MATHML_TAGS[tagName]);
    }
    if (element.namespaceURI === HTML_NAMESPACE) {
      // The only way to switch from SVG to HTML is via
      // HTML integration points, and from MathML to HTML
      // is via MathML text integration points
      if (parent.namespaceURI === SVG_NAMESPACE && !HTML_INTEGRATION_POINTS[parentTagName]) {
        return false;
      }
      if (parent.namespaceURI === MATHML_NAMESPACE && !MATHML_TEXT_INTEGRATION_POINTS[parentTagName]) {
        return false;
      }
      // We disallow tags that are specific for MathML
      // or SVG and should never appear in HTML namespace
      return !ALL_MATHML_TAGS[tagName] && (COMMON_SVG_AND_HTML_ELEMENTS[tagName] || !ALL_SVG_TAGS[tagName]);
    }
    // For XHTML and XML documents that support custom namespaces
    if (PARSER_MEDIA_TYPE === 'application/xhtml+xml' && ALLOWED_NAMESPACES[element.namespaceURI]) {
      return true;
    }
    // The code should never reach this place (this means
    // that the element somehow got namespace that is not
    // HTML, SVG, MathML or allowed via ALLOWED_NAMESPACES).
    // Return false just in case.
    return false;
  };
  /**
   * _forceRemove
   *
   * @param node a DOM node
   */
  const _forceRemove = function _forceRemove(node) {
    arrayPush(DOMPurify.removed, {
      element: node
    });
    try {
      // eslint-disable-next-line unicorn/prefer-dom-node-remove
      getParentNode(node).removeChild(node);
    } catch (_) {
      remove(node);
    }
  };
  /**
   * _removeAttribute
   *
   * @param name an Attribute name
   * @param element a DOM node
   */
  const _removeAttribute = function _removeAttribute(name, element) {
    try {
      arrayPush(DOMPurify.removed, {
        attribute: element.getAttributeNode(name),
        from: element
      });
    } catch (_) {
      arrayPush(DOMPurify.removed, {
        attribute: null,
        from: element
      });
    }
    element.removeAttribute(name);
    // We void attribute values for unremovable "is" attributes
    if (name === 'is') {
      if (RETURN_DOM || RETURN_DOM_FRAGMENT) {
        try {
          _forceRemove(element);
        } catch (_) {}
      } else {
        try {
          element.setAttribute(name, '');
        } catch (_) {}
      }
    }
  };
  /**
   * _initDocument
   *
   * @param dirty - a string of dirty markup
   * @return a DOM, filled with the dirty markup
   */
  const _initDocument = function _initDocument(dirty) {
    /* Create a HTML document */
    let doc = null;
    let leadingWhitespace = null;
    if (FORCE_BODY) {
      dirty = '<remove></remove>' + dirty;
    } else {
      /* If FORCE_BODY isn't used, leading whitespace needs to be preserved manually */
      const matches = stringMatch(dirty, /^[\r\n\t ]+/);
      leadingWhitespace = matches && matches[0];
    }
    if (PARSER_MEDIA_TYPE === 'application/xhtml+xml' && NAMESPACE === HTML_NAMESPACE) {
      // Root of XHTML doc must contain xmlns declaration (see https://www.w3.org/TR/xhtml1/normative.html#strict)
      dirty = '<html xmlns="http://www.w3.org/1999/xhtml"><head></head><body>' + dirty + '</body></html>';
    }
    const dirtyPayload = trustedTypesPolicy ? trustedTypesPolicy.createHTML(dirty) : dirty;
    /*
     * Use the DOMParser API by default, fallback later if needs be
     * DOMParser not work for svg when has multiple root element.
     */
    if (NAMESPACE === HTML_NAMESPACE) {
      try {
        doc = new DOMParser().parseFromString(dirtyPayload, PARSER_MEDIA_TYPE);
      } catch (_) {}
    }
    /* Use createHTMLDocument in case DOMParser is not available */
    if (!doc || !doc.documentElement) {
      doc = implementation.createDocument(NAMESPACE, 'template', null);
      try {
        doc.documentElement.innerHTML = IS_EMPTY_INPUT ? emptyHTML : dirtyPayload;
      } catch (_) {
        // Syntax error if dirtyPayload is invalid xml
      }
    }
    const body = doc.body || doc.documentElement;
    if (dirty && leadingWhitespace) {
      body.insertBefore(document.createTextNode(leadingWhitespace), body.childNodes[0] || null);
    }
    /* Work on whole document or just its body */
    if (NAMESPACE === HTML_NAMESPACE) {
      return getElementsByTagName.call(doc, WHOLE_DOCUMENT ? 'html' : 'body')[0];
    }
    return WHOLE_DOCUMENT ? doc.documentElement : body;
  };
  /**
   * Creates a NodeIterator object that you can use to traverse filtered lists of nodes or elements in a document.
   *
   * @param root The root element or node to start traversing on.
   * @return The created NodeIterator
   */
  const _createNodeIterator = function _createNodeIterator(root) {
    return createNodeIterator.call(root.ownerDocument || root, root,
    // eslint-disable-next-line no-bitwise
    NodeFilter.SHOW_ELEMENT | NodeFilter.SHOW_COMMENT | NodeFilter.SHOW_TEXT | NodeFilter.SHOW_PROCESSING_INSTRUCTION | NodeFilter.SHOW_CDATA_SECTION, null);
  };
  /**
   * _isClobbered
   *
   * @param element element to check for clobbering attacks
   * @return true if clobbered, false if safe
   */
  const _isClobbered = function _isClobbered(element) {
    return element instanceof HTMLFormElement && (typeof element.nodeName !== 'string' || typeof element.textContent !== 'string' || typeof element.removeChild !== 'function' || !(element.attributes instanceof NamedNodeMap) || typeof element.removeAttribute !== 'function' || typeof element.setAttribute !== 'function' || typeof element.namespaceURI !== 'string' || typeof element.insertBefore !== 'function' || typeof element.hasChildNodes !== 'function');
  };
  /**
   * Checks whether the given object is a DOM node.
   *
   * @param value object to check whether it's a DOM node
   * @return true is object is a DOM node
   */
  const _isNode = function _isNode(value) {
    return typeof Node === 'function' && value instanceof Node;
  };
  function _executeHooks(hooks, currentNode, data) {
    arrayForEach(hooks, hook => {
      hook.call(DOMPurify, currentNode, data, CONFIG);
    });
  }
  /**
   * _sanitizeElements
   *
   * @protect nodeName
   * @protect textContent
   * @protect removeChild
   * @param currentNode to check for permission to exist
   * @return true if node was killed, false if left alive
   */
  const _sanitizeElements = function _sanitizeElements(currentNode) {
    let content = null;
    /* Execute a hook if present */
    _executeHooks(hooks.beforeSanitizeElements, currentNode, null);
    /* Check if element is clobbered or can clobber */
    if (_isClobbered(currentNode)) {
      _forceRemove(currentNode);
      return true;
    }
    /* Now let's check the element's type and name */
    const tagName = transformCaseFunc(currentNode.nodeName);
    /* Execute a hook if present */
    _executeHooks(hooks.uponSanitizeElement, currentNode, {
      tagName,
      allowedTags: ALLOWED_TAGS
    });
    /* Detect mXSS attempts abusing namespace confusion */
    if (SAFE_FOR_XML && currentNode.hasChildNodes() && !_isNode(currentNode.firstElementChild) && regExpTest(/<[/\w!]/g, currentNode.innerHTML) && regExpTest(/<[/\w!]/g, currentNode.textContent)) {
      _forceRemove(currentNode);
      return true;
    }
    /* Remove any occurrence of processing instructions */
    if (currentNode.nodeType === NODE_TYPE.progressingInstruction) {
      _forceRemove(currentNode);
      return true;
    }
    /* Remove any kind of possibly harmful comments */
    if (SAFE_FOR_XML && currentNode.nodeType === NODE_TYPE.comment && regExpTest(/<[/\w]/g, currentNode.data)) {
      _forceRemove(currentNode);
      return true;
    }
    /* Remove element if anything forbids its presence */
    if (!ALLOWED_TAGS[tagName] || FORBID_TAGS[tagName]) {
      /* Check if we have a custom element to handle */
      if (!FORBID_TAGS[tagName] && _isBasicCustomElement(tagName)) {
        if (CUSTOM_ELEMENT_HANDLING.tagNameCheck instanceof RegExp && regExpTest(CUSTOM_ELEMENT_HANDLING.tagNameCheck, tagName)) {
          return false;
        }
        if (CUSTOM_ELEMENT_HANDLING.tagNameCheck instanceof Function && CUSTOM_ELEMENT_HANDLING.tagNameCheck(tagName)) {
          return false;
        }
      }
      /* Keep content except for bad-listed elements */
      if (KEEP_CONTENT && !FORBID_CONTENTS[tagName]) {
        const parentNode = getParentNode(currentNode) || currentNode.parentNode;
        const childNodes = getChildNodes(currentNode) || currentNode.childNodes;
        if (childNodes && parentNode) {
          const childCount = childNodes.length;
          for (let i = childCount - 1; i >= 0; --i) {
            const childClone = cloneNode(childNodes[i], true);
            childClone.__removalCount = (currentNode.__removalCount || 0) + 1;
            parentNode.insertBefore(childClone, getNextSibling(currentNode));
          }
        }
      }
      _forceRemove(currentNode);
      return true;
    }
    /* Check whether element has a valid namespace */
    if (currentNode instanceof Element && !_checkValidNamespace(currentNode)) {
      _forceRemove(currentNode);
      return true;
    }
    /* Make sure that older browsers don't get fallback-tag mXSS */
    if ((tagName === 'noscript' || tagName === 'noembed' || tagName === 'noframes') && regExpTest(/<\/no(script|embed|frames)/i, currentNode.innerHTML)) {
      _forceRemove(currentNode);
      return true;
    }
    /* Sanitize element content to be template-safe */
    if (SAFE_FOR_TEMPLATES && currentNode.nodeType === NODE_TYPE.text) {
      /* Get the element's text content */
      content = currentNode.textContent;
      arrayForEach([MUSTACHE_EXPR, ERB_EXPR, TMPLIT_EXPR], expr => {
        content = stringReplace(content, expr, ' ');
      });
      if (currentNode.textContent !== content) {
        arrayPush(DOMPurify.removed, {
          element: currentNode.cloneNode()
        });
        currentNode.textContent = content;
      }
    }
    /* Execute a hook if present */
    _executeHooks(hooks.afterSanitizeElements, currentNode, null);
    return false;
  };
  /**
   * _isValidAttribute
   *
   * @param lcTag Lowercase tag name of containing element.
   * @param lcName Lowercase attribute name.
   * @param value Attribute value.
   * @return Returns true if `value` is valid, otherwise false.
   */
  // eslint-disable-next-line complexity
  const _isValidAttribute = function _isValidAttribute(lcTag, lcName, value) {
    /* Make sure attribute cannot clobber */
    if (SANITIZE_DOM && (lcName === 'id' || lcName === 'name') && (value in document || value in formElement)) {
      return false;
    }
    /* Allow valid data-* attributes: At least one character after "-"
        (https://html.spec.whatwg.org/multipage/dom.html#embedding-custom-non-visible-data-with-the-data-*-attributes)
        XML-compatible (https://html.spec.whatwg.org/multipage/infrastructure.html#xml-compatible and http://www.w3.org/TR/xml/#d0e804)
        We don't need to check the value; it's always URI safe. */
    if (ALLOW_DATA_ATTR && !FORBID_ATTR[lcName] && regExpTest(DATA_ATTR, lcName)) ; else if (ALLOW_ARIA_ATTR && regExpTest(ARIA_ATTR, lcName)) ; else if (!ALLOWED_ATTR[lcName] || FORBID_ATTR[lcName]) {
      if (
      // First condition does a very basic check if a) it's basically a valid custom element tagname AND
      // b) if the tagName passes whatever the user has configured for CUSTOM_ELEMENT_HANDLING.tagNameCheck
      // and c) if the attribute name passes whatever the user has configured for CUSTOM_ELEMENT_HANDLING.attributeNameCheck
      _isBasicCustomElement(lcTag) && (CUSTOM_ELEMENT_HANDLING.tagNameCheck instanceof RegExp && regExpTest(CUSTOM_ELEMENT_HANDLING.tagNameCheck, lcTag) || CUSTOM_ELEMENT_HANDLING.tagNameCheck instanceof Function && CUSTOM_ELEMENT_HANDLING.tagNameCheck(lcTag)) && (CUSTOM_ELEMENT_HANDLING.attributeNameCheck instanceof RegExp && regExpTest(CUSTOM_ELEMENT_HANDLING.attributeNameCheck, lcName) || CUSTOM_ELEMENT_HANDLING.attributeNameCheck instanceof Function && CUSTOM_ELEMENT_HANDLING.attributeNameCheck(lcName)) ||
      // Alternative, second condition checks if it's an `is`-attribute, AND
      // the value passes whatever the user has configured for CUSTOM_ELEMENT_HANDLING.tagNameCheck
      lcName === 'is' && CUSTOM_ELEMENT_HANDLING.allowCustomizedBuiltInElements && (CUSTOM_ELEMENT_HANDLING.tagNameCheck instanceof RegExp && regExpTest(CUSTOM_ELEMENT_HANDLING.tagNameCheck, value) || CUSTOM_ELEMENT_HANDLING.tagNameCheck instanceof Function && CUSTOM_ELEMENT_HANDLING.tagNameCheck(value))) ; else {
        return false;
      }
      /* Check value is safe. First, is attr inert? If so, is safe */
    } else if (URI_SAFE_ATTRIBUTES[lcName]) ; else if (regExpTest(IS_ALLOWED_URI$1, stringReplace(value, ATTR_WHITESPACE, ''))) ; else if ((lcName === 'src' || lcName === 'xlink:href' || lcName === 'href') && lcTag !== 'script' && stringIndexOf(value, 'data:') === 0 && DATA_URI_TAGS[lcTag]) ; else if (ALLOW_UNKNOWN_PROTOCOLS && !regExpTest(IS_SCRIPT_OR_DATA, stringReplace(value, ATTR_WHITESPACE, ''))) ; else if (value) {
      return false;
    } else ;
    return true;
  };
  /**
   * _isBasicCustomElement
   * checks if at least one dash is included in tagName, and it's not the first char
   * for more sophisticated checking see https://github.com/sindresorhus/validate-element-name
   *
   * @param tagName name of the tag of the node to sanitize
   * @returns Returns true if the tag name meets the basic criteria for a custom element, otherwise false.
   */
  const _isBasicCustomElement = function _isBasicCustomElement(tagName) {
    return tagName !== 'annotation-xml' && stringMatch(tagName, CUSTOM_ELEMENT);
  };
  /**
   * _sanitizeAttributes
   *
   * @protect attributes
   * @protect nodeName
   * @protect removeAttribute
   * @protect setAttribute
   *
   * @param currentNode to sanitize
   */
  const _sanitizeAttributes = function _sanitizeAttributes(currentNode) {
    /* Execute a hook if present */
    _executeHooks(hooks.beforeSanitizeAttributes, currentNode, null);
    const {
      attributes
    } = currentNode;
    /* Check if we have attributes; if not we might have a text node */
    if (!attributes || _isClobbered(currentNode)) {
      return;
    }
    const hookEvent = {
      attrName: '',
      attrValue: '',
      keepAttr: true,
      allowedAttributes: ALLOWED_ATTR,
      forceKeepAttr: undefined
    };
    let l = attributes.length;
    /* Go backwards over all attributes; safely remove bad ones */
    while (l--) {
      const attr = attributes[l];
      const {
        name,
        namespaceURI,
        value: attrValue
      } = attr;
      const lcName = transformCaseFunc(name);
      const initValue = attrValue;
      let value = name === 'value' ? initValue : stringTrim(initValue);
      /* Execute a hook if present */
      hookEvent.attrName = lcName;
      hookEvent.attrValue = value;
      hookEvent.keepAttr = true;
      hookEvent.forceKeepAttr = undefined; // Allows developers to see this is a property they can set
      _executeHooks(hooks.uponSanitizeAttribute, currentNode, hookEvent);
      value = hookEvent.attrValue;
      /* Full DOM Clobbering protection via namespace isolation,
       * Prefix id and name attributes with `user-content-`
       */
      if (SANITIZE_NAMED_PROPS && (lcName === 'id' || lcName === 'name')) {
        // Remove the attribute with this value
        _removeAttribute(name, currentNode);
        // Prefix the value and later re-create the attribute with the sanitized value
        value = SANITIZE_NAMED_PROPS_PREFIX + value;
      }
      /* Work around a security issue with comments inside attributes */
      if (SAFE_FOR_XML && regExpTest(/((--!?|])>)|<\/(style|title)/i, value)) {
        _removeAttribute(name, currentNode);
        continue;
      }
      /* Did the hooks approve of the attribute? */
      if (hookEvent.forceKeepAttr) {
        continue;
      }
      /* Did the hooks approve of the attribute? */
      if (!hookEvent.keepAttr) {
        _removeAttribute(name, currentNode);
        continue;
      }
      /* Work around a security issue in jQuery 3.0 */
      if (!ALLOW_SELF_CLOSE_IN_ATTR && regExpTest(/\/>/i, value)) {
        _removeAttribute(name, currentNode);
        continue;
      }
      /* Sanitize attribute content to be template-safe */
      if (SAFE_FOR_TEMPLATES) {
        arrayForEach([MUSTACHE_EXPR, ERB_EXPR, TMPLIT_EXPR], expr => {
          value = stringReplace(value, expr, ' ');
        });
      }
      /* Is `value` valid for this attribute? */
      const lcTag = transformCaseFunc(currentNode.nodeName);
      if (!_isValidAttribute(lcTag, lcName, value)) {
        _removeAttribute(name, currentNode);
        continue;
      }
      /* Handle attributes that require Trusted Types */
      if (trustedTypesPolicy && typeof trustedTypes === 'object' && typeof trustedTypes.getAttributeType === 'function') {
        if (namespaceURI) ; else {
          switch (trustedTypes.getAttributeType(lcTag, lcName)) {
            case 'TrustedHTML':
              {
                value = trustedTypesPolicy.createHTML(value);
                break;
              }
            case 'TrustedScriptURL':
              {
                value = trustedTypesPolicy.createScriptURL(value);
                break;
              }
          }
        }
      }
      /* Handle invalid data-* attribute set by try-catching it */
      if (value !== initValue) {
        try {
          if (namespaceURI) {
            currentNode.setAttributeNS(namespaceURI, name, value);
          } else {
            /* Fallback to setAttribute() for browser-unrecognized namespaces e.g. "x-schema". */
            currentNode.setAttribute(name, value);
          }
          if (_isClobbered(currentNode)) {
            _forceRemove(currentNode);
          } else {
            arrayPop(DOMPurify.removed);
          }
        } catch (_) {
          _removeAttribute(name, currentNode);
        }
      }
    }
    /* Execute a hook if present */
    _executeHooks(hooks.afterSanitizeAttributes, currentNode, null);
  };
  /**
   * _sanitizeShadowDOM
   *
   * @param fragment to iterate over recursively
   */
  const _sanitizeShadowDOM = function _sanitizeShadowDOM(fragment) {
    let shadowNode = null;
    const shadowIterator = _createNodeIterator(fragment);
    /* Execute a hook if present */
    _executeHooks(hooks.beforeSanitizeShadowDOM, fragment, null);
    while (shadowNode = shadowIterator.nextNode()) {
      /* Execute a hook if present */
      _executeHooks(hooks.uponSanitizeShadowNode, shadowNode, null);
      /* Sanitize tags and elements */
      _sanitizeElements(shadowNode);
      /* Check attributes next */
      _sanitizeAttributes(shadowNode);
      /* Deep shadow DOM detected */
      if (shadowNode.content instanceof DocumentFragment) {
        _sanitizeShadowDOM(shadowNode.content);
      }
    }
    /* Execute a hook if present */
    _executeHooks(hooks.afterSanitizeShadowDOM, fragment, null);
  };
  // eslint-disable-next-line complexity
  DOMPurify.sanitize = function (dirty) {
    let cfg = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : {};
    let body = null;
    let importedNode = null;
    let currentNode = null;
    let returnNode = null;
    /* Make sure we have a string to sanitize.
      DO NOT return early, as this will return the wrong type if
      the user has requested a DOM object rather than a string */
    IS_EMPTY_INPUT = !dirty;
    if (IS_EMPTY_INPUT) {
      dirty = '<!-->';
    }
    /* Stringify, in case dirty is an object */
    if (typeof dirty !== 'string' && !_isNode(dirty)) {
      if (typeof dirty.toString === 'function') {
        dirty = dirty.toString();
        if (typeof dirty !== 'string') {
          throw typeErrorCreate('dirty is not a string, aborting');
        }
      } else {
        throw typeErrorCreate('toString is not a function');
      }
    }
    /* Return dirty HTML if DOMPurify cannot run */
    if (!DOMPurify.isSupported) {
      return dirty;
    }
    /* Assign config vars */
    if (!SET_CONFIG) {
      _parseConfig(cfg);
    }
    /* Clean up removed elements */
    DOMPurify.removed = [];
    /* Check if dirty is correctly typed for IN_PLACE */
    if (typeof dirty === 'string') {
      IN_PLACE = false;
    }
    if (IN_PLACE) {
      /* Do some early pre-sanitization to avoid unsafe root nodes */
      if (dirty.nodeName) {
        const tagName = transformCaseFunc(dirty.nodeName);
        if (!ALLOWED_TAGS[tagName] || FORBID_TAGS[tagName]) {
          throw typeErrorCreate('root node is forbidden and cannot be sanitized in-place');
        }
      }
    } else if (dirty instanceof Node) {
      /* If dirty is a DOM element, append to an empty document to avoid
         elements being stripped by the parser */
      body = _initDocument('<!---->');
      importedNode = body.ownerDocument.importNode(dirty, true);
      if (importedNode.nodeType === NODE_TYPE.element && importedNode.nodeName === 'BODY') {
        /* Node is already a body, use as is */
        body = importedNode;
      } else if (importedNode.nodeName === 'HTML') {
        body = importedNode;
      } else {
        // eslint-disable-next-line unicorn/prefer-dom-node-append
        body.appendChild(importedNode);
      }
    } else {
      /* Exit directly if we have nothing to do */
      if (!RETURN_DOM && !SAFE_FOR_TEMPLATES && !WHOLE_DOCUMENT &&
      // eslint-disable-next-line unicorn/prefer-includes
      dirty.indexOf('<') === -1) {
        return trustedTypesPolicy && RETURN_TRUSTED_TYPE ? trustedTypesPolicy.createHTML(dirty) : dirty;
      }
      /* Initialize the document to work on */
      body = _initDocument(dirty);
      /* Check we have a DOM node from the data */
      if (!body) {
        return RETURN_DOM ? null : RETURN_TRUSTED_TYPE ? emptyHTML : '';
      }
    }
    /* Remove first element node (ours) if FORCE_BODY is set */
    if (body && FORCE_BODY) {
      _forceRemove(body.firstChild);
    }
    /* Get node iterator */
    const nodeIterator = _createNodeIterator(IN_PLACE ? dirty : body);
    /* Now start iterating over the created document */
    while (currentNode = nodeIterator.nextNode()) {
      /* Sanitize tags and elements */
      _sanitizeElements(currentNode);
      /* Check attributes next */
      _sanitizeAttributes(currentNode);
      /* Shadow DOM detected, sanitize it */
      if (currentNode.content instanceof DocumentFragment) {
        _sanitizeShadowDOM(currentNode.content);
      }
    }
    /* If we sanitized `dirty` in-place, return it. */
    if (IN_PLACE) {
      return dirty;
    }
    /* Return sanitized string or DOM */
    if (RETURN_DOM) {
      if (RETURN_DOM_FRAGMENT) {
        returnNode = createDocumentFragment.call(body.ownerDocument);
        while (body.firstChild) {
          // eslint-disable-next-line unicorn/prefer-dom-node-append
          returnNode.appendChild(body.firstChild);
        }
      } else {
        returnNode = body;
      }
      if (ALLOWED_ATTR.shadowroot || ALLOWED_ATTR.shadowrootmode) {
        /*
          AdoptNode() is not used because internal state is not reset
          (e.g. the past names map of a HTMLFormElement), this is safe
          in theory but we would rather not risk another attack vector.
          The state that is cloned by importNode() is explicitly defined
          by the specs.
        */
        returnNode = importNode.call(originalDocument, returnNode, true);
      }
      return returnNode;
    }
    let serializedHTML = WHOLE_DOCUMENT ? body.outerHTML : body.innerHTML;
    /* Serialize doctype if allowed */
    if (WHOLE_DOCUMENT && ALLOWED_TAGS['!doctype'] && body.ownerDocument && body.ownerDocument.doctype && body.ownerDocument.doctype.name && regExpTest(DOCTYPE_NAME, body.ownerDocument.doctype.name)) {
      serializedHTML = '<!DOCTYPE ' + body.ownerDocument.doctype.name + '>\n' + serializedHTML;
    }
    /* Sanitize final string template-safe */
    if (SAFE_FOR_TEMPLATES) {
      arrayForEach([MUSTACHE_EXPR, ERB_EXPR, TMPLIT_EXPR], expr => {
        serializedHTML = stringReplace(serializedHTML, expr, ' ');
      });
    }
    return trustedTypesPolicy && RETURN_TRUSTED_TYPE ? trustedTypesPolicy.createHTML(serializedHTML) : serializedHTML;
  };
  DOMPurify.setConfig = function () {
    let cfg = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : {};
    _parseConfig(cfg);
    SET_CONFIG = true;
  };
  DOMPurify.clearConfig = function () {
    CONFIG = null;
    SET_CONFIG = false;
  };
  DOMPurify.isValidAttribute = function (tag, attr, value) {
    /* Initialize shared config vars if necessary. */
    if (!CONFIG) {
      _parseConfig({});
    }
    const lcTag = transformCaseFunc(tag);
    const lcName = transformCaseFunc(attr);
    return _isValidAttribute(lcTag, lcName, value);
  };
  DOMPurify.addHook = function (entryPoint, hookFunction) {
    if (typeof hookFunction !== 'function') {
      return;
    }
    arrayPush(hooks[entryPoint], hookFunction);
  };
  DOMPurify.removeHook = function (entryPoint, hookFunction) {
    if (hookFunction !== undefined) {
      const index = arrayLastIndexOf(hooks[entryPoint], hookFunction);
      return index === -1 ? undefined : arraySplice(hooks[entryPoint], index, 1)[0];
    }
    return arrayPop(hooks[entryPoint]);
  };
  DOMPurify.removeHooks = function (entryPoint) {
    hooks[entryPoint] = [];
  };
  DOMPurify.removeAllHooks = function () {
    hooks = _createHooksMap();
  };
  return DOMPurify;
}
var purify = createDOMPurify();


//# sourceMappingURL=purify.es.mjs.map


/***/ }),

/***/ 54758:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ reusable)
});

// EXTERNAL MODULE: ../node_modules/khroma/dist/utils/index.js + 3 modules
var utils = __webpack_require__(90267);
// EXTERNAL MODULE: ../node_modules/khroma/dist/constants.js
var constants = __webpack_require__(62187);
;// CONCATENATED MODULE: ../node_modules/khroma/dist/channels/type.js
/* IMPORT */

/* MAIN */
class Type {
    constructor() {
        /* VARIABLES */
        this.type = constants/* TYPE */.w.ALL;
    }
    /* API */
    get() {
        return this.type;
    }
    set(type) {
        if (this.type && this.type !== type)
            throw new Error('Cannot change both RGB and HSL channels at the same time');
        this.type = type;
    }
    reset() {
        this.type = constants/* TYPE */.w.ALL;
    }
    is(type) {
        return this.type === type;
    }
}
/* EXPORT */
/* harmony default export */ const type = (Type);

;// CONCATENATED MODULE: ../node_modules/khroma/dist/channels/index.js
/* IMPORT */



/* MAIN */
class Channels {
    /* CONSTRUCTOR */
    constructor(data, color) {
        this.color = color;
        this.changed = false;
        this.data = data; //TSC
        this.type = new type();
    }
    /* API */
    set(data, color) {
        this.color = color;
        this.changed = false;
        this.data = data; //TSC
        this.type.type = constants/* TYPE */.w.ALL;
        return this;
    }
    /* HELPERS */
    _ensureHSL() {
        const data = this.data;
        const { h, s, l } = data;
        if (h === undefined)
            data.h = utils/* default */.Z.channel.rgb2hsl(data, 'h');
        if (s === undefined)
            data.s = utils/* default */.Z.channel.rgb2hsl(data, 's');
        if (l === undefined)
            data.l = utils/* default */.Z.channel.rgb2hsl(data, 'l');
    }
    _ensureRGB() {
        const data = this.data;
        const { r, g, b } = data;
        if (r === undefined)
            data.r = utils/* default */.Z.channel.hsl2rgb(data, 'r');
        if (g === undefined)
            data.g = utils/* default */.Z.channel.hsl2rgb(data, 'g');
        if (b === undefined)
            data.b = utils/* default */.Z.channel.hsl2rgb(data, 'b');
    }
    /* GETTERS */
    get r() {
        const data = this.data;
        const r = data.r;
        if (!this.type.is(constants/* TYPE */.w.HSL) && r !== undefined)
            return r;
        this._ensureHSL();
        return utils/* default */.Z.channel.hsl2rgb(data, 'r');
    }
    get g() {
        const data = this.data;
        const g = data.g;
        if (!this.type.is(constants/* TYPE */.w.HSL) && g !== undefined)
            return g;
        this._ensureHSL();
        return utils/* default */.Z.channel.hsl2rgb(data, 'g');
    }
    get b() {
        const data = this.data;
        const b = data.b;
        if (!this.type.is(constants/* TYPE */.w.HSL) && b !== undefined)
            return b;
        this._ensureHSL();
        return utils/* default */.Z.channel.hsl2rgb(data, 'b');
    }
    get h() {
        const data = this.data;
        const h = data.h;
        if (!this.type.is(constants/* TYPE */.w.RGB) && h !== undefined)
            return h;
        this._ensureRGB();
        return utils/* default */.Z.channel.rgb2hsl(data, 'h');
    }
    get s() {
        const data = this.data;
        const s = data.s;
        if (!this.type.is(constants/* TYPE */.w.RGB) && s !== undefined)
            return s;
        this._ensureRGB();
        return utils/* default */.Z.channel.rgb2hsl(data, 's');
    }
    get l() {
        const data = this.data;
        const l = data.l;
        if (!this.type.is(constants/* TYPE */.w.RGB) && l !== undefined)
            return l;
        this._ensureRGB();
        return utils/* default */.Z.channel.rgb2hsl(data, 'l');
    }
    get a() {
        return this.data.a;
    }
    /* SETTERS */
    set r(r) {
        this.type.set(constants/* TYPE */.w.RGB);
        this.changed = true;
        this.data.r = r;
    }
    set g(g) {
        this.type.set(constants/* TYPE */.w.RGB);
        this.changed = true;
        this.data.g = g;
    }
    set b(b) {
        this.type.set(constants/* TYPE */.w.RGB);
        this.changed = true;
        this.data.b = b;
    }
    set h(h) {
        this.type.set(constants/* TYPE */.w.HSL);
        this.changed = true;
        this.data.h = h;
    }
    set s(s) {
        this.type.set(constants/* TYPE */.w.HSL);
        this.changed = true;
        this.data.s = s;
    }
    set l(l) {
        this.type.set(constants/* TYPE */.w.HSL);
        this.changed = true;
        this.data.l = l;
    }
    set a(a) {
        this.changed = true;
        this.data.a = a;
    }
}
/* EXPORT */
/* harmony default export */ const channels = (Channels);

;// CONCATENATED MODULE: ../node_modules/khroma/dist/channels/reusable.js
/* IMPORT */

/* MAIN */
const reusable_channels = new channels({ r: 0, g: 0, b: 0, a: 0 }, 'transparent');
/* EXPORT */
/* harmony default export */ const reusable = (reusable_channels);


/***/ }),

/***/ 42528:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ color)
});

// EXTERNAL MODULE: ../node_modules/khroma/dist/channels/reusable.js + 2 modules
var reusable = __webpack_require__(54758);
// EXTERNAL MODULE: ../node_modules/khroma/dist/constants.js
var constants = __webpack_require__(62187);
;// CONCATENATED MODULE: ../node_modules/khroma/dist/color/hex.js
/* IMPORT */



/* MAIN */
const Hex = {
    /* VARIABLES */
    re: /^#((?:[a-f0-9]{2}){2,4}|[a-f0-9]{3})$/i,
    /* API */
    parse: (color) => {
        if (color.charCodeAt(0) !== 35)
            return; // '#'
        const match = color.match(Hex.re);
        if (!match)
            return;
        const hex = match[1];
        const dec = parseInt(hex, 16);
        const length = hex.length;
        const hasAlpha = length % 4 === 0;
        const isFullLength = length > 4;
        const multiplier = isFullLength ? 1 : 17;
        const bits = isFullLength ? 8 : 4;
        const bitsOffset = hasAlpha ? 0 : -1;
        const mask = isFullLength ? 255 : 15;
        return reusable/* default */.Z.set({
            r: ((dec >> (bits * (bitsOffset + 3))) & mask) * multiplier,
            g: ((dec >> (bits * (bitsOffset + 2))) & mask) * multiplier,
            b: ((dec >> (bits * (bitsOffset + 1))) & mask) * multiplier,
            a: hasAlpha ? (dec & mask) * multiplier / 255 : 1
        }, color);
    },
    stringify: (channels) => {
        const { r, g, b, a } = channels;
        if (a < 1) { // #RRGGBBAA
            return `#${constants/* DEC2HEX */.Q[Math.round(r)]}${constants/* DEC2HEX */.Q[Math.round(g)]}${constants/* DEC2HEX */.Q[Math.round(b)]}${constants/* DEC2HEX */.Q[Math.round(a * 255)]}`;
        }
        else { // #RRGGBB
            return `#${constants/* DEC2HEX */.Q[Math.round(r)]}${constants/* DEC2HEX */.Q[Math.round(g)]}${constants/* DEC2HEX */.Q[Math.round(b)]}`;
        }
    }
};
/* EXPORT */
/* harmony default export */ const color_hex = (Hex);

// EXTERNAL MODULE: ../node_modules/khroma/dist/utils/index.js + 3 modules
var utils = __webpack_require__(90267);
;// CONCATENATED MODULE: ../node_modules/khroma/dist/color/hsl.js
/* IMPORT */


/* MAIN */
const HSL = {
    /* VARIABLES */
    re: /^hsla?\(\s*?(-?(?:\d+(?:\.\d+)?|(?:\.\d+))(?:e-?\d+)?(?:deg|grad|rad|turn)?)\s*?(?:,|\s)\s*?(-?(?:\d+(?:\.\d+)?|(?:\.\d+))(?:e-?\d+)?%)\s*?(?:,|\s)\s*?(-?(?:\d+(?:\.\d+)?|(?:\.\d+))(?:e-?\d+)?%)(?:\s*?(?:,|\/)\s*?\+?(-?(?:\d+(?:\.\d+)?|(?:\.\d+))(?:e-?\d+)?(%)?))?\s*?\)$/i,
    hueRe: /^(.+?)(deg|grad|rad|turn)$/i,
    /* HELPERS */
    _hue2deg: (hue) => {
        const match = hue.match(HSL.hueRe);
        if (match) {
            const [, number, unit] = match;
            switch (unit) {
                case 'grad': return utils/* default */.Z.channel.clamp.h(parseFloat(number) * .9);
                case 'rad': return utils/* default */.Z.channel.clamp.h(parseFloat(number) * 180 / Math.PI);
                case 'turn': return utils/* default */.Z.channel.clamp.h(parseFloat(number) * 360);
            }
        }
        return utils/* default */.Z.channel.clamp.h(parseFloat(hue));
    },
    /* API */
    parse: (color) => {
        const charCode = color.charCodeAt(0);
        if (charCode !== 104 && charCode !== 72)
            return; // 'h'/'H'
        const match = color.match(HSL.re);
        if (!match)
            return;
        const [, h, s, l, a, isAlphaPercentage] = match;
        return reusable/* default */.Z.set({
            h: HSL._hue2deg(h),
            s: utils/* default */.Z.channel.clamp.s(parseFloat(s)),
            l: utils/* default */.Z.channel.clamp.l(parseFloat(l)),
            a: a ? utils/* default */.Z.channel.clamp.a(isAlphaPercentage ? parseFloat(a) / 100 : parseFloat(a)) : 1
        }, color);
    },
    stringify: (channels) => {
        const { h, s, l, a } = channels;
        if (a < 1) { // HSLA
            return `hsla(${utils/* default */.Z.lang.round(h)}, ${utils/* default */.Z.lang.round(s)}%, ${utils/* default */.Z.lang.round(l)}%, ${a})`;
        }
        else { // HSL
            return `hsl(${utils/* default */.Z.lang.round(h)}, ${utils/* default */.Z.lang.round(s)}%, ${utils/* default */.Z.lang.round(l)}%)`;
        }
    }
};
/* EXPORT */
/* harmony default export */ const hsl = (HSL);

;// CONCATENATED MODULE: ../node_modules/khroma/dist/color/keyword.js
/* IMPORT */

/* MAIN */
const Keyword = {
    /* VARIABLES */
    colors: {
        aliceblue: '#f0f8ff',
        antiquewhite: '#faebd7',
        aqua: '#00ffff',
        aquamarine: '#7fffd4',
        azure: '#f0ffff',
        beige: '#f5f5dc',
        bisque: '#ffe4c4',
        black: '#000000',
        blanchedalmond: '#ffebcd',
        blue: '#0000ff',
        blueviolet: '#8a2be2',
        brown: '#a52a2a',
        burlywood: '#deb887',
        cadetblue: '#5f9ea0',
        chartreuse: '#7fff00',
        chocolate: '#d2691e',
        coral: '#ff7f50',
        cornflowerblue: '#6495ed',
        cornsilk: '#fff8dc',
        crimson: '#dc143c',
        cyanaqua: '#00ffff',
        darkblue: '#00008b',
        darkcyan: '#008b8b',
        darkgoldenrod: '#b8860b',
        darkgray: '#a9a9a9',
        darkgreen: '#006400',
        darkgrey: '#a9a9a9',
        darkkhaki: '#bdb76b',
        darkmagenta: '#8b008b',
        darkolivegreen: '#556b2f',
        darkorange: '#ff8c00',
        darkorchid: '#9932cc',
        darkred: '#8b0000',
        darksalmon: '#e9967a',
        darkseagreen: '#8fbc8f',
        darkslateblue: '#483d8b',
        darkslategray: '#2f4f4f',
        darkslategrey: '#2f4f4f',
        darkturquoise: '#00ced1',
        darkviolet: '#9400d3',
        deeppink: '#ff1493',
        deepskyblue: '#00bfff',
        dimgray: '#696969',
        dimgrey: '#696969',
        dodgerblue: '#1e90ff',
        firebrick: '#b22222',
        floralwhite: '#fffaf0',
        forestgreen: '#228b22',
        fuchsia: '#ff00ff',
        gainsboro: '#dcdcdc',
        ghostwhite: '#f8f8ff',
        gold: '#ffd700',
        goldenrod: '#daa520',
        gray: '#808080',
        green: '#008000',
        greenyellow: '#adff2f',
        grey: '#808080',
        honeydew: '#f0fff0',
        hotpink: '#ff69b4',
        indianred: '#cd5c5c',
        indigo: '#4b0082',
        ivory: '#fffff0',
        khaki: '#f0e68c',
        lavender: '#e6e6fa',
        lavenderblush: '#fff0f5',
        lawngreen: '#7cfc00',
        lemonchiffon: '#fffacd',
        lightblue: '#add8e6',
        lightcoral: '#f08080',
        lightcyan: '#e0ffff',
        lightgoldenrodyellow: '#fafad2',
        lightgray: '#d3d3d3',
        lightgreen: '#90ee90',
        lightgrey: '#d3d3d3',
        lightpink: '#ffb6c1',
        lightsalmon: '#ffa07a',
        lightseagreen: '#20b2aa',
        lightskyblue: '#87cefa',
        lightslategray: '#778899',
        lightslategrey: '#778899',
        lightsteelblue: '#b0c4de',
        lightyellow: '#ffffe0',
        lime: '#00ff00',
        limegreen: '#32cd32',
        linen: '#faf0e6',
        magenta: '#ff00ff',
        maroon: '#800000',
        mediumaquamarine: '#66cdaa',
        mediumblue: '#0000cd',
        mediumorchid: '#ba55d3',
        mediumpurple: '#9370db',
        mediumseagreen: '#3cb371',
        mediumslateblue: '#7b68ee',
        mediumspringgreen: '#00fa9a',
        mediumturquoise: '#48d1cc',
        mediumvioletred: '#c71585',
        midnightblue: '#191970',
        mintcream: '#f5fffa',
        mistyrose: '#ffe4e1',
        moccasin: '#ffe4b5',
        navajowhite: '#ffdead',
        navy: '#000080',
        oldlace: '#fdf5e6',
        olive: '#808000',
        olivedrab: '#6b8e23',
        orange: '#ffa500',
        orangered: '#ff4500',
        orchid: '#da70d6',
        palegoldenrod: '#eee8aa',
        palegreen: '#98fb98',
        paleturquoise: '#afeeee',
        palevioletred: '#db7093',
        papayawhip: '#ffefd5',
        peachpuff: '#ffdab9',
        peru: '#cd853f',
        pink: '#ffc0cb',
        plum: '#dda0dd',
        powderblue: '#b0e0e6',
        purple: '#800080',
        rebeccapurple: '#663399',
        red: '#ff0000',
        rosybrown: '#bc8f8f',
        royalblue: '#4169e1',
        saddlebrown: '#8b4513',
        salmon: '#fa8072',
        sandybrown: '#f4a460',
        seagreen: '#2e8b57',
        seashell: '#fff5ee',
        sienna: '#a0522d',
        silver: '#c0c0c0',
        skyblue: '#87ceeb',
        slateblue: '#6a5acd',
        slategray: '#708090',
        slategrey: '#708090',
        snow: '#fffafa',
        springgreen: '#00ff7f',
        tan: '#d2b48c',
        teal: '#008080',
        thistle: '#d8bfd8',
        transparent: '#00000000',
        turquoise: '#40e0d0',
        violet: '#ee82ee',
        wheat: '#f5deb3',
        white: '#ffffff',
        whitesmoke: '#f5f5f5',
        yellow: '#ffff00',
        yellowgreen: '#9acd32'
    },
    /* API */
    parse: (color) => {
        color = color.toLowerCase();
        const hex = Keyword.colors[color];
        if (!hex)
            return;
        return color_hex.parse(hex);
    },
    stringify: (channels) => {
        const hex = color_hex.stringify(channels);
        for (const name in Keyword.colors) {
            if (Keyword.colors[name] === hex)
                return name;
        }
        return;
    }
};
/* EXPORT */
/* harmony default export */ const keyword = (Keyword);

;// CONCATENATED MODULE: ../node_modules/khroma/dist/color/rgb.js
/* IMPORT */


/* MAIN */
const RGB = {
    /* VARIABLES */
    re: /^rgba?\(\s*?(-?(?:\d+(?:\.\d+)?|(?:\.\d+))(?:e\d+)?(%?))\s*?(?:,|\s)\s*?(-?(?:\d+(?:\.\d+)?|(?:\.\d+))(?:e\d+)?(%?))\s*?(?:,|\s)\s*?(-?(?:\d+(?:\.\d+)?|(?:\.\d+))(?:e\d+)?(%?))(?:\s*?(?:,|\/)\s*?\+?(-?(?:\d+(?:\.\d+)?|(?:\.\d+))(?:e\d+)?(%?)))?\s*?\)$/i,
    /* API */
    parse: (color) => {
        const charCode = color.charCodeAt(0);
        if (charCode !== 114 && charCode !== 82)
            return; // 'r'/'R'
        const match = color.match(RGB.re);
        if (!match)
            return;
        const [, r, isRedPercentage, g, isGreenPercentage, b, isBluePercentage, a, isAlphaPercentage] = match;
        return reusable/* default */.Z.set({
            r: utils/* default */.Z.channel.clamp.r(isRedPercentage ? parseFloat(r) * 2.55 : parseFloat(r)),
            g: utils/* default */.Z.channel.clamp.g(isGreenPercentage ? parseFloat(g) * 2.55 : parseFloat(g)),
            b: utils/* default */.Z.channel.clamp.b(isBluePercentage ? parseFloat(b) * 2.55 : parseFloat(b)),
            a: a ? utils/* default */.Z.channel.clamp.a(isAlphaPercentage ? parseFloat(a) / 100 : parseFloat(a)) : 1
        }, color);
    },
    stringify: (channels) => {
        const { r, g, b, a } = channels;
        if (a < 1) { // RGBA
            return `rgba(${utils/* default */.Z.lang.round(r)}, ${utils/* default */.Z.lang.round(g)}, ${utils/* default */.Z.lang.round(b)}, ${utils/* default */.Z.lang.round(a)})`;
        }
        else { // RGB
            return `rgb(${utils/* default */.Z.lang.round(r)}, ${utils/* default */.Z.lang.round(g)}, ${utils/* default */.Z.lang.round(b)})`;
        }
    }
};
/* EXPORT */
/* harmony default export */ const rgb = (RGB);

;// CONCATENATED MODULE: ../node_modules/khroma/dist/color/index.js
/* IMPORT */






/* MAIN */
const Color = {
    /* VARIABLES */
    format: {
        keyword: keyword,
        hex: color_hex,
        rgb: rgb,
        rgba: rgb,
        hsl: hsl,
        hsla: hsl
    },
    /* API */
    parse: (color) => {
        if (typeof color !== 'string')
            return color;
        const channels = color_hex.parse(color) || rgb.parse(color) || hsl.parse(color) || keyword.parse(color); // Color providers ordered with performance in mind
        if (channels)
            return channels;
        throw new Error(`Unsupported color format: "${color}"`);
    },
    stringify: (channels) => {
        // SASS returns a keyword if possible, but we avoid doing that as it's slower and doesn't really add any value
        if (!channels.changed && channels.color)
            return channels.color;
        if (channels.type.is(constants/* TYPE */.w.HSL) || channels.data.r === undefined) {
            return hsl.stringify(channels);
        }
        else if (channels.a < 1 || !Number.isInteger(channels.r) || !Number.isInteger(channels.g) || !Number.isInteger(channels.b)) {
            return rgb.stringify(channels);
        }
        else {
            return color_hex.stringify(channels);
        }
    }
};
/* EXPORT */
/* harmony default export */ const color = (Color);


/***/ }),

/***/ 62187:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Q: () => (/* binding */ DEC2HEX),
/* harmony export */   w: () => (/* binding */ TYPE)
/* harmony export */ });
/* harmony import */ var _utils_index_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(90267);
/* IMPORT */

/* MAIN */
const DEC2HEX = {};
for (let i = 0; i <= 255; i++)
    DEC2HEX[i] = _utils_index_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z.unit.dec2hex(i); // Populating dynamically, striking a balance between code size and performance
const TYPE = {
    ALL: 0,
    RGB: 1,
    HSL: 2
};
/* EXPORT */



/***/ }),

/***/ 60680:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _utils_index_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(90267);
/* harmony import */ var _color_index_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(42528);
/* IMPORT */


/* MAIN */
const adjustChannel = (color, channel, amount) => {
    const channels = _color_index_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z.parse(color);
    const amountCurrent = channels[channel];
    const amountNext = _utils_index_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z.channel.clamp[channel](amountCurrent + amount);
    if (amountCurrent !== amountNext)
        channels[channel] = amountNext;
    return _color_index_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z.stringify(channels);
};
/* EXPORT */
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (adjustChannel);


/***/ }),

/***/ 95323:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _utils_index_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(90267);
/* harmony import */ var _color_index_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(42528);
/* IMPORT */


/* MAIN */
const change = (color, channels) => {
    const ch = _color_index_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z.parse(color);
    for (const c in channels) {
        ch[c] = _utils_index_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z.channel.clamp[c](channels[c]);
    }
    return _color_index_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z.stringify(ch);
};
/* EXPORT */
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (change);


/***/ }),

/***/ 57838:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _adjust_channel_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(60680);
/* IMPORT */

/* MAIN */
const darken = (color, amount) => {
    return (0,_adjust_channel_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(color, 'l', -amount);
};
/* EXPORT */
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (darken);


/***/ }),

/***/ 28186:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ is_dark)
});

// EXTERNAL MODULE: ../node_modules/khroma/dist/utils/index.js + 3 modules
var utils = __webpack_require__(90267);
// EXTERNAL MODULE: ../node_modules/khroma/dist/color/index.js + 4 modules
var dist_color = __webpack_require__(42528);
;// CONCATENATED MODULE: ../node_modules/khroma/dist/methods/luminance.js
/* IMPORT */


/* MAIN */
//SOURCE: https://planetcalc.com/7779
const luminance = (color) => {
    const { r, g, b } = dist_color/* default */.Z.parse(color);
    const luminance = .2126 * utils/* default */.Z.channel.toLinear(r) + .7152 * utils/* default */.Z.channel.toLinear(g) + .0722 * utils/* default */.Z.channel.toLinear(b);
    return utils/* default */.Z.lang.round(luminance);
};
/* EXPORT */
/* harmony default export */ const methods_luminance = (luminance);

;// CONCATENATED MODULE: ../node_modules/khroma/dist/methods/is_light.js
/* IMPORT */

/* MAIN */
const isLight = (color) => {
    return methods_luminance(color) >= .5;
};
/* EXPORT */
/* harmony default export */ const is_light = (isLight);

;// CONCATENATED MODULE: ../node_modules/khroma/dist/methods/is_dark.js
/* IMPORT */

/* MAIN */
const isDark = (color) => {
    return !is_light(color);
};
/* EXPORT */
/* harmony default export */ const is_dark = (isDark);


/***/ }),

/***/ 28482:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _adjust_channel_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(60680);
/* IMPORT */

/* MAIN */
const lighten = (color, amount) => {
    return (0,_adjust_channel_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(color, 'l', amount);
};
/* EXPORT */
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (lighten);


/***/ }),

/***/ 14728:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _utils_index_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(90267);
/* harmony import */ var _channels_reusable_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(54758);
/* harmony import */ var _color_index_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(42528);
/* harmony import */ var _change_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(95323);
/* IMPORT */




/* MAIN */
const rgba = (r, g, b = 0, a = 1) => {
    if (typeof r !== 'number')
        return (0,_change_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(r, { a: g });
    const channels = _channels_reusable_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z.set({
        r: _utils_index_js__WEBPACK_IMPORTED_MODULE_2__/* ["default"] */ .Z.channel.clamp.r(r),
        g: _utils_index_js__WEBPACK_IMPORTED_MODULE_2__/* ["default"] */ .Z.channel.clamp.g(g),
        b: _utils_index_js__WEBPACK_IMPORTED_MODULE_2__/* ["default"] */ .Z.channel.clamp.b(b),
        a: _utils_index_js__WEBPACK_IMPORTED_MODULE_2__/* ["default"] */ .Z.channel.clamp.a(a)
    });
    return _color_index_js__WEBPACK_IMPORTED_MODULE_3__/* ["default"] */ .Z.stringify(channels);
};
/* EXPORT */
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (rgba);


/***/ }),

/***/ 90267:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ utils)
});

;// CONCATENATED MODULE: ../node_modules/khroma/dist/utils/channel.js
/* IMPORT */
/* MAIN */
const Channel = {
    /* CLAMP */
    min: {
        r: 0,
        g: 0,
        b: 0,
        s: 0,
        l: 0,
        a: 0
    },
    max: {
        r: 255,
        g: 255,
        b: 255,
        h: 360,
        s: 100,
        l: 100,
        a: 1
    },
    clamp: {
        r: (r) => r >= 255 ? 255 : (r < 0 ? 0 : r),
        g: (g) => g >= 255 ? 255 : (g < 0 ? 0 : g),
        b: (b) => b >= 255 ? 255 : (b < 0 ? 0 : b),
        h: (h) => h % 360,
        s: (s) => s >= 100 ? 100 : (s < 0 ? 0 : s),
        l: (l) => l >= 100 ? 100 : (l < 0 ? 0 : l),
        a: (a) => a >= 1 ? 1 : (a < 0 ? 0 : a)
    },
    /* CONVERSION */
    //SOURCE: https://planetcalc.com/7779
    toLinear: (c) => {
        const n = c / 255;
        return c > .03928 ? Math.pow(((n + .055) / 1.055), 2.4) : n / 12.92;
    },
    //SOURCE: https://gist.github.com/mjackson/5311256
    hue2rgb: (p, q, t) => {
        if (t < 0)
            t += 1;
        if (t > 1)
            t -= 1;
        if (t < 1 / 6)
            return p + (q - p) * 6 * t;
        if (t < 1 / 2)
            return q;
        if (t < 2 / 3)
            return p + (q - p) * (2 / 3 - t) * 6;
        return p;
    },
    hsl2rgb: ({ h, s, l }, channel) => {
        if (!s)
            return l * 2.55; // Achromatic
        h /= 360;
        s /= 100;
        l /= 100;
        const q = (l < .5) ? l * (1 + s) : (l + s) - (l * s);
        const p = 2 * l - q;
        switch (channel) {
            case 'r': return Channel.hue2rgb(p, q, h + 1 / 3) * 255;
            case 'g': return Channel.hue2rgb(p, q, h) * 255;
            case 'b': return Channel.hue2rgb(p, q, h - 1 / 3) * 255;
        }
    },
    rgb2hsl: ({ r, g, b }, channel) => {
        r /= 255;
        g /= 255;
        b /= 255;
        const max = Math.max(r, g, b);
        const min = Math.min(r, g, b);
        const l = (max + min) / 2;
        if (channel === 'l')
            return l * 100;
        if (max === min)
            return 0; // Achromatic
        const d = max - min;
        const s = (l > .5) ? d / (2 - max - min) : d / (max + min);
        if (channel === 's')
            return s * 100;
        switch (max) {
            case r: return ((g - b) / d + (g < b ? 6 : 0)) * 60;
            case g: return ((b - r) / d + 2) * 60;
            case b: return ((r - g) / d + 4) * 60;
            default: return -1; //TSC: TypeScript is stupid and complains if there isn't this useless default statement
        }
    }
};
/* EXPORT */
/* harmony default export */ const channel = (Channel);

;// CONCATENATED MODULE: ../node_modules/khroma/dist/utils/lang.js
/* MAIN */
const Lang = {
    /* API */
    clamp: (number, lower, upper) => {
        if (lower > upper)
            return Math.min(lower, Math.max(upper, number));
        return Math.min(upper, Math.max(lower, number));
    },
    round: (number) => {
        return Math.round(number * 10000000000) / 10000000000;
    }
};
/* EXPORT */
/* harmony default export */ const lang = (Lang);

;// CONCATENATED MODULE: ../node_modules/khroma/dist/utils/unit.js
/* MAIN */
const Unit = {
    /* API */
    dec2hex: (dec) => {
        const hex = Math.round(dec).toString(16);
        return hex.length > 1 ? hex : `0${hex}`;
    }
};
/* EXPORT */
/* harmony default export */ const unit = (Unit);

;// CONCATENATED MODULE: ../node_modules/khroma/dist/utils/index.js
/* IMPORT */



/* MAIN */
const Utils = {
    channel: channel,
    lang: lang,
    unit: unit
};
/* EXPORT */
/* harmony default export */ const utils = (Utils);


/***/ }),

/***/ 91138:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ _ListCache)
});

;// CONCATENATED MODULE: ../node_modules/lodash-es/_listCacheClear.js
/**
 * Removes all key-value entries from the list cache.
 *
 * @private
 * @name clear
 * @memberOf ListCache
 */
function listCacheClear() {
  this.__data__ = [];
  this.size = 0;
}

/* harmony default export */ const _listCacheClear = (listCacheClear);

// EXTERNAL MODULE: ../node_modules/lodash-es/eq.js
var eq = __webpack_require__(35050);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_assocIndexOf.js


/**
 * Gets the index at which the `key` is found in `array` of key-value pairs.
 *
 * @private
 * @param {Array} array The array to inspect.
 * @param {*} key The key to search for.
 * @returns {number} Returns the index of the matched value, else `-1`.
 */
function assocIndexOf(array, key) {
  var length = array.length;
  while (length--) {
    if ((0,eq/* default */.Z)(array[length][0], key)) {
      return length;
    }
  }
  return -1;
}

/* harmony default export */ const _assocIndexOf = (assocIndexOf);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_listCacheDelete.js


/** Used for built-in method references. */
var arrayProto = Array.prototype;

/** Built-in value references. */
var splice = arrayProto.splice;

/**
 * Removes `key` and its value from the list cache.
 *
 * @private
 * @name delete
 * @memberOf ListCache
 * @param {string} key The key of the value to remove.
 * @returns {boolean} Returns `true` if the entry was removed, else `false`.
 */
function listCacheDelete(key) {
  var data = this.__data__,
      index = _assocIndexOf(data, key);

  if (index < 0) {
    return false;
  }
  var lastIndex = data.length - 1;
  if (index == lastIndex) {
    data.pop();
  } else {
    splice.call(data, index, 1);
  }
  --this.size;
  return true;
}

/* harmony default export */ const _listCacheDelete = (listCacheDelete);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_listCacheGet.js


/**
 * Gets the list cache value for `key`.
 *
 * @private
 * @name get
 * @memberOf ListCache
 * @param {string} key The key of the value to get.
 * @returns {*} Returns the entry value.
 */
function listCacheGet(key) {
  var data = this.__data__,
      index = _assocIndexOf(data, key);

  return index < 0 ? undefined : data[index][1];
}

/* harmony default export */ const _listCacheGet = (listCacheGet);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_listCacheHas.js


/**
 * Checks if a list cache value for `key` exists.
 *
 * @private
 * @name has
 * @memberOf ListCache
 * @param {string} key The key of the entry to check.
 * @returns {boolean} Returns `true` if an entry for `key` exists, else `false`.
 */
function listCacheHas(key) {
  return _assocIndexOf(this.__data__, key) > -1;
}

/* harmony default export */ const _listCacheHas = (listCacheHas);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_listCacheSet.js


/**
 * Sets the list cache `key` to `value`.
 *
 * @private
 * @name set
 * @memberOf ListCache
 * @param {string} key The key of the value to set.
 * @param {*} value The value to set.
 * @returns {Object} Returns the list cache instance.
 */
function listCacheSet(key, value) {
  var data = this.__data__,
      index = _assocIndexOf(data, key);

  if (index < 0) {
    ++this.size;
    data.push([key, value]);
  } else {
    data[index][1] = value;
  }
  return this;
}

/* harmony default export */ const _listCacheSet = (listCacheSet);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_ListCache.js






/**
 * Creates an list cache object.
 *
 * @private
 * @constructor
 * @param {Array} [entries] The key-value pairs to cache.
 */
function ListCache(entries) {
  var index = -1,
      length = entries == null ? 0 : entries.length;

  this.clear();
  while (++index < length) {
    var entry = entries[index];
    this.set(entry[0], entry[1]);
  }
}

// Add methods to `ListCache`.
ListCache.prototype.clear = _listCacheClear;
ListCache.prototype['delete'] = _listCacheDelete;
ListCache.prototype.get = _listCacheGet;
ListCache.prototype.has = _listCacheHas;
ListCache.prototype.set = _listCacheSet;

/* harmony default export */ const _ListCache = (ListCache);


/***/ }),

/***/ 81700:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _getNative_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(26266);
/* harmony import */ var _root_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(94311);



/* Built-in method references that are verified to be native. */
var Map = (0,_getNative_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(_root_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z, 'Map');

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (Map);


/***/ }),

/***/ 24395:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ _MapCache)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/_getNative.js + 4 modules
var _getNative = __webpack_require__(26266);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_nativeCreate.js


/* Built-in method references that are verified to be native. */
var nativeCreate = (0,_getNative/* default */.Z)(Object, 'create');

/* harmony default export */ const _nativeCreate = (nativeCreate);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_hashClear.js


/**
 * Removes all key-value entries from the hash.
 *
 * @private
 * @name clear
 * @memberOf Hash
 */
function hashClear() {
  this.__data__ = _nativeCreate ? _nativeCreate(null) : {};
  this.size = 0;
}

/* harmony default export */ const _hashClear = (hashClear);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_hashDelete.js
/**
 * Removes `key` and its value from the hash.
 *
 * @private
 * @name delete
 * @memberOf Hash
 * @param {Object} hash The hash to modify.
 * @param {string} key The key of the value to remove.
 * @returns {boolean} Returns `true` if the entry was removed, else `false`.
 */
function hashDelete(key) {
  var result = this.has(key) && delete this.__data__[key];
  this.size -= result ? 1 : 0;
  return result;
}

/* harmony default export */ const _hashDelete = (hashDelete);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_hashGet.js


/** Used to stand-in for `undefined` hash values. */
var HASH_UNDEFINED = '__lodash_hash_undefined__';

/** Used for built-in method references. */
var objectProto = Object.prototype;

/** Used to check objects for own properties. */
var _hashGet_hasOwnProperty = objectProto.hasOwnProperty;

/**
 * Gets the hash value for `key`.
 *
 * @private
 * @name get
 * @memberOf Hash
 * @param {string} key The key of the value to get.
 * @returns {*} Returns the entry value.
 */
function hashGet(key) {
  var data = this.__data__;
  if (_nativeCreate) {
    var result = data[key];
    return result === HASH_UNDEFINED ? undefined : result;
  }
  return _hashGet_hasOwnProperty.call(data, key) ? data[key] : undefined;
}

/* harmony default export */ const _hashGet = (hashGet);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_hashHas.js


/** Used for built-in method references. */
var _hashHas_objectProto = Object.prototype;

/** Used to check objects for own properties. */
var _hashHas_hasOwnProperty = _hashHas_objectProto.hasOwnProperty;

/**
 * Checks if a hash value for `key` exists.
 *
 * @private
 * @name has
 * @memberOf Hash
 * @param {string} key The key of the entry to check.
 * @returns {boolean} Returns `true` if an entry for `key` exists, else `false`.
 */
function hashHas(key) {
  var data = this.__data__;
  return _nativeCreate ? (data[key] !== undefined) : _hashHas_hasOwnProperty.call(data, key);
}

/* harmony default export */ const _hashHas = (hashHas);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_hashSet.js


/** Used to stand-in for `undefined` hash values. */
var _hashSet_HASH_UNDEFINED = '__lodash_hash_undefined__';

/**
 * Sets the hash `key` to `value`.
 *
 * @private
 * @name set
 * @memberOf Hash
 * @param {string} key The key of the value to set.
 * @param {*} value The value to set.
 * @returns {Object} Returns the hash instance.
 */
function hashSet(key, value) {
  var data = this.__data__;
  this.size += this.has(key) ? 0 : 1;
  data[key] = (_nativeCreate && value === undefined) ? _hashSet_HASH_UNDEFINED : value;
  return this;
}

/* harmony default export */ const _hashSet = (hashSet);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_Hash.js






/**
 * Creates a hash object.
 *
 * @private
 * @constructor
 * @param {Array} [entries] The key-value pairs to cache.
 */
function Hash(entries) {
  var index = -1,
      length = entries == null ? 0 : entries.length;

  this.clear();
  while (++index < length) {
    var entry = entries[index];
    this.set(entry[0], entry[1]);
  }
}

// Add methods to `Hash`.
Hash.prototype.clear = _hashClear;
Hash.prototype['delete'] = _hashDelete;
Hash.prototype.get = _hashGet;
Hash.prototype.has = _hashHas;
Hash.prototype.set = _hashSet;

/* harmony default export */ const _Hash = (Hash);

// EXTERNAL MODULE: ../node_modules/lodash-es/_ListCache.js + 6 modules
var _ListCache = __webpack_require__(91138);
// EXTERNAL MODULE: ../node_modules/lodash-es/_Map.js
var _Map = __webpack_require__(81700);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_mapCacheClear.js




/**
 * Removes all key-value entries from the map.
 *
 * @private
 * @name clear
 * @memberOf MapCache
 */
function mapCacheClear() {
  this.size = 0;
  this.__data__ = {
    'hash': new _Hash,
    'map': new (_Map/* default */.Z || _ListCache/* default */.Z),
    'string': new _Hash
  };
}

/* harmony default export */ const _mapCacheClear = (mapCacheClear);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_isKeyable.js
/**
 * Checks if `value` is suitable for use as unique object key.
 *
 * @private
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is suitable, else `false`.
 */
function isKeyable(value) {
  var type = typeof value;
  return (type == 'string' || type == 'number' || type == 'symbol' || type == 'boolean')
    ? (value !== '__proto__')
    : (value === null);
}

/* harmony default export */ const _isKeyable = (isKeyable);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_getMapData.js


/**
 * Gets the data for `map`.
 *
 * @private
 * @param {Object} map The map to query.
 * @param {string} key The reference key.
 * @returns {*} Returns the map data.
 */
function getMapData(map, key) {
  var data = map.__data__;
  return _isKeyable(key)
    ? data[typeof key == 'string' ? 'string' : 'hash']
    : data.map;
}

/* harmony default export */ const _getMapData = (getMapData);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_mapCacheDelete.js


/**
 * Removes `key` and its value from the map.
 *
 * @private
 * @name delete
 * @memberOf MapCache
 * @param {string} key The key of the value to remove.
 * @returns {boolean} Returns `true` if the entry was removed, else `false`.
 */
function mapCacheDelete(key) {
  var result = _getMapData(this, key)['delete'](key);
  this.size -= result ? 1 : 0;
  return result;
}

/* harmony default export */ const _mapCacheDelete = (mapCacheDelete);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_mapCacheGet.js


/**
 * Gets the map value for `key`.
 *
 * @private
 * @name get
 * @memberOf MapCache
 * @param {string} key The key of the value to get.
 * @returns {*} Returns the entry value.
 */
function mapCacheGet(key) {
  return _getMapData(this, key).get(key);
}

/* harmony default export */ const _mapCacheGet = (mapCacheGet);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_mapCacheHas.js


/**
 * Checks if a map value for `key` exists.
 *
 * @private
 * @name has
 * @memberOf MapCache
 * @param {string} key The key of the entry to check.
 * @returns {boolean} Returns `true` if an entry for `key` exists, else `false`.
 */
function mapCacheHas(key) {
  return _getMapData(this, key).has(key);
}

/* harmony default export */ const _mapCacheHas = (mapCacheHas);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_mapCacheSet.js


/**
 * Sets the map `key` to `value`.
 *
 * @private
 * @name set
 * @memberOf MapCache
 * @param {string} key The key of the value to set.
 * @param {*} value The value to set.
 * @returns {Object} Returns the map cache instance.
 */
function mapCacheSet(key, value) {
  var data = _getMapData(this, key),
      size = data.size;

  data.set(key, value);
  this.size += data.size == size ? 0 : 1;
  return this;
}

/* harmony default export */ const _mapCacheSet = (mapCacheSet);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_MapCache.js






/**
 * Creates a map cache object to store key-value pairs.
 *
 * @private
 * @constructor
 * @param {Array} [entries] The key-value pairs to cache.
 */
function MapCache(entries) {
  var index = -1,
      length = entries == null ? 0 : entries.length;

  this.clear();
  while (++index < length) {
    var entry = entries[index];
    this.set(entry[0], entry[1]);
  }
}

// Add methods to `MapCache`.
MapCache.prototype.clear = _mapCacheClear;
MapCache.prototype['delete'] = _mapCacheDelete;
MapCache.prototype.get = _mapCacheGet;
MapCache.prototype.has = _mapCacheHas;
MapCache.prototype.set = _mapCacheSet;

/* harmony default export */ const _MapCache = (MapCache);


/***/ }),

/***/ 16889:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _getNative_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(26266);
/* harmony import */ var _root_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(94311);



/* Built-in method references that are verified to be native. */
var Set = (0,_getNative_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(_root_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z, 'Set');

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (Set);


/***/ }),

/***/ 82948:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ _Stack)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/_ListCache.js + 6 modules
var _ListCache = __webpack_require__(91138);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_stackClear.js


/**
 * Removes all key-value entries from the stack.
 *
 * @private
 * @name clear
 * @memberOf Stack
 */
function stackClear() {
  this.__data__ = new _ListCache/* default */.Z;
  this.size = 0;
}

/* harmony default export */ const _stackClear = (stackClear);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_stackDelete.js
/**
 * Removes `key` and its value from the stack.
 *
 * @private
 * @name delete
 * @memberOf Stack
 * @param {string} key The key of the value to remove.
 * @returns {boolean} Returns `true` if the entry was removed, else `false`.
 */
function stackDelete(key) {
  var data = this.__data__,
      result = data['delete'](key);

  this.size = data.size;
  return result;
}

/* harmony default export */ const _stackDelete = (stackDelete);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_stackGet.js
/**
 * Gets the stack value for `key`.
 *
 * @private
 * @name get
 * @memberOf Stack
 * @param {string} key The key of the value to get.
 * @returns {*} Returns the entry value.
 */
function stackGet(key) {
  return this.__data__.get(key);
}

/* harmony default export */ const _stackGet = (stackGet);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_stackHas.js
/**
 * Checks if a stack value for `key` exists.
 *
 * @private
 * @name has
 * @memberOf Stack
 * @param {string} key The key of the entry to check.
 * @returns {boolean} Returns `true` if an entry for `key` exists, else `false`.
 */
function stackHas(key) {
  return this.__data__.has(key);
}

/* harmony default export */ const _stackHas = (stackHas);

// EXTERNAL MODULE: ../node_modules/lodash-es/_Map.js
var _Map = __webpack_require__(81700);
// EXTERNAL MODULE: ../node_modules/lodash-es/_MapCache.js + 14 modules
var _MapCache = __webpack_require__(24395);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_stackSet.js




/** Used as the size to enable large array optimizations. */
var LARGE_ARRAY_SIZE = 200;

/**
 * Sets the stack `key` to `value`.
 *
 * @private
 * @name set
 * @memberOf Stack
 * @param {string} key The key of the value to set.
 * @param {*} value The value to set.
 * @returns {Object} Returns the stack cache instance.
 */
function stackSet(key, value) {
  var data = this.__data__;
  if (data instanceof _ListCache/* default */.Z) {
    var pairs = data.__data__;
    if (!_Map/* default */.Z || (pairs.length < LARGE_ARRAY_SIZE - 1)) {
      pairs.push([key, value]);
      this.size = ++data.size;
      return this;
    }
    data = this.__data__ = new _MapCache/* default */.Z(pairs);
  }
  data.set(key, value);
  this.size = data.size;
  return this;
}

/* harmony default export */ const _stackSet = (stackSet);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_Stack.js







/**
 * Creates a stack cache object to store key-value pairs.
 *
 * @private
 * @constructor
 * @param {Array} [entries] The key-value pairs to cache.
 */
function Stack(entries) {
  var data = this.__data__ = new _ListCache/* default */.Z(entries);
  this.size = data.size;
}

// Add methods to `Stack`.
Stack.prototype.clear = _stackClear;
Stack.prototype['delete'] = _stackDelete;
Stack.prototype.get = _stackGet;
Stack.prototype.has = _stackHas;
Stack.prototype.set = _stackSet;

/* harmony default export */ const _Stack = (Stack);


/***/ }),

/***/ 91642:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _root_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(94311);


/** Built-in value references. */
var Symbol = _root_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z.Symbol;

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (Symbol);


/***/ }),

/***/ 41049:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _root_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(94311);


/** Built-in value references. */
var Uint8Array = _root_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z.Uint8Array;

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (Uint8Array);


/***/ }),

/***/ 40709:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ _arrayLikeKeys)
});

;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseTimes.js
/**
 * The base implementation of `_.times` without support for iteratee shorthands
 * or max array length checks.
 *
 * @private
 * @param {number} n The number of times to invoke `iteratee`.
 * @param {Function} iteratee The function invoked per iteration.
 * @returns {Array} Returns the array of results.
 */
function baseTimes(n, iteratee) {
  var index = -1,
      result = Array(n);

  while (++index < n) {
    result[index] = iteratee(index);
  }
  return result;
}

/* harmony default export */ const _baseTimes = (baseTimes);

// EXTERNAL MODULE: ../node_modules/lodash-es/isArguments.js + 1 modules
var isArguments = __webpack_require__(9028);
// EXTERNAL MODULE: ../node_modules/lodash-es/isArray.js
var isArray = __webpack_require__(64058);
// EXTERNAL MODULE: ../node_modules/lodash-es/isBuffer.js + 1 modules
var isBuffer = __webpack_require__(23230);
// EXTERNAL MODULE: ../node_modules/lodash-es/_isIndex.js
var _isIndex = __webpack_require__(8616);
// EXTERNAL MODULE: ../node_modules/lodash-es/isTypedArray.js + 1 modules
var isTypedArray = __webpack_require__(14923);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_arrayLikeKeys.js







/** Used for built-in method references. */
var objectProto = Object.prototype;

/** Used to check objects for own properties. */
var _arrayLikeKeys_hasOwnProperty = objectProto.hasOwnProperty;

/**
 * Creates an array of the enumerable property names of the array-like `value`.
 *
 * @private
 * @param {*} value The value to query.
 * @param {boolean} inherited Specify returning inherited property names.
 * @returns {Array} Returns the array of property names.
 */
function arrayLikeKeys(value, inherited) {
  var isArr = (0,isArray/* default */.Z)(value),
      isArg = !isArr && (0,isArguments/* default */.Z)(value),
      isBuff = !isArr && !isArg && (0,isBuffer/* default */.Z)(value),
      isType = !isArr && !isArg && !isBuff && (0,isTypedArray/* default */.Z)(value),
      skipIndexes = isArr || isArg || isBuff || isType,
      result = skipIndexes ? _baseTimes(value.length, String) : [],
      length = result.length;

  for (var key in value) {
    if ((inherited || _arrayLikeKeys_hasOwnProperty.call(value, key)) &&
        !(skipIndexes && (
           // Safari 9 has enumerable `arguments.length` in strict mode.
           key == 'length' ||
           // Node.js 0.10 has enumerable non-index properties on buffers.
           (isBuff && (key == 'offset' || key == 'parent')) ||
           // PhantomJS 2 has enumerable non-index properties on typed arrays.
           (isType && (key == 'buffer' || key == 'byteLength' || key == 'byteOffset')) ||
           // Skip index properties.
           (0,_isIndex/* default */.Z)(key, length)
        ))) {
      result.push(key);
    }
  }
  return result;
}

/* harmony default export */ const _arrayLikeKeys = (arrayLikeKeys);


/***/ }),

/***/ 15561:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _baseAssignValue_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(93586);
/* harmony import */ var _eq_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(35050);



/** Used for built-in method references. */
var objectProto = Object.prototype;

/** Used to check objects for own properties. */
var hasOwnProperty = objectProto.hasOwnProperty;

/**
 * Assigns `value` to `key` of `object` if the existing value is not equivalent
 * using [`SameValueZero`](http://ecma-international.org/ecma-262/7.0/#sec-samevaluezero)
 * for equality comparisons.
 *
 * @private
 * @param {Object} object The object to modify.
 * @param {string} key The key of the property to assign.
 * @param {*} value The value to assign.
 */
function assignValue(object, key, value) {
  var objValue = object[key];
  if (!(hasOwnProperty.call(object, key) && (0,_eq_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(objValue, value)) ||
      (value === undefined && !(key in object))) {
    (0,_baseAssignValue_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z)(object, key, value);
  }
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (assignValue);


/***/ }),

/***/ 93586:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _defineProperty_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(30253);


/**
 * The base implementation of `assignValue` and `assignMergeValue` without
 * value checks.
 *
 * @private
 * @param {Object} object The object to modify.
 * @param {string} key The key of the property to assign.
 * @param {*} value The value to assign.
 */
function baseAssignValue(object, key, value) {
  if (key == '__proto__' && _defineProperty_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z) {
    (0,_defineProperty_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(object, key, {
      'configurable': true,
      'enumerable': true,
      'value': value,
      'writable': true
    });
  } else {
    object[key] = value;
  }
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (baseAssignValue);


/***/ }),

/***/ 49399:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ _baseFor)
});

;// CONCATENATED MODULE: ../node_modules/lodash-es/_createBaseFor.js
/**
 * Creates a base function for methods like `_.forIn` and `_.forOwn`.
 *
 * @private
 * @param {boolean} [fromRight] Specify iterating from right to left.
 * @returns {Function} Returns the new base function.
 */
function createBaseFor(fromRight) {
  return function(object, iteratee, keysFunc) {
    var index = -1,
        iterable = Object(object),
        props = keysFunc(object),
        length = props.length;

    while (length--) {
      var key = props[fromRight ? length : ++index];
      if (iteratee(iterable[key], key, iterable) === false) {
        break;
      }
    }
    return object;
  };
}

/* harmony default export */ const _createBaseFor = (createBaseFor);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseFor.js


/**
 * The base implementation of `baseForOwn` which iterates over `object`
 * properties returned by `keysFunc` and invokes `iteratee` for each property.
 * Iteratee functions may exit iteration early by explicitly returning `false`.
 *
 * @private
 * @param {Object} object The object to iterate over.
 * @param {Function} iteratee The function invoked per iteration.
 * @param {Function} keysFunc The function to get the keys of `object`.
 * @returns {Object} Returns `object`.
 */
var baseFor = _createBaseFor();

/* harmony default export */ const _baseFor = (baseFor);


/***/ }),

/***/ 77070:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ _baseGetTag)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/_Symbol.js
var _Symbol = __webpack_require__(91642);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_getRawTag.js


/** Used for built-in method references. */
var objectProto = Object.prototype;

/** Used to check objects for own properties. */
var _getRawTag_hasOwnProperty = objectProto.hasOwnProperty;

/**
 * Used to resolve the
 * [`toStringTag`](http://ecma-international.org/ecma-262/7.0/#sec-object.prototype.tostring)
 * of values.
 */
var nativeObjectToString = objectProto.toString;

/** Built-in value references. */
var symToStringTag = _Symbol/* default */.Z ? _Symbol/* default */.Z.toStringTag : undefined;

/**
 * A specialized version of `baseGetTag` which ignores `Symbol.toStringTag` values.
 *
 * @private
 * @param {*} value The value to query.
 * @returns {string} Returns the raw `toStringTag`.
 */
function getRawTag(value) {
  var isOwn = _getRawTag_hasOwnProperty.call(value, symToStringTag),
      tag = value[symToStringTag];

  try {
    value[symToStringTag] = undefined;
    var unmasked = true;
  } catch (e) {}

  var result = nativeObjectToString.call(value);
  if (unmasked) {
    if (isOwn) {
      value[symToStringTag] = tag;
    } else {
      delete value[symToStringTag];
    }
  }
  return result;
}

/* harmony default export */ const _getRawTag = (getRawTag);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_objectToString.js
/** Used for built-in method references. */
var _objectToString_objectProto = Object.prototype;

/**
 * Used to resolve the
 * [`toStringTag`](http://ecma-international.org/ecma-262/7.0/#sec-object.prototype.tostring)
 * of values.
 */
var _objectToString_nativeObjectToString = _objectToString_objectProto.toString;

/**
 * Converts `value` to a string using `Object.prototype.toString`.
 *
 * @private
 * @param {*} value The value to convert.
 * @returns {string} Returns the converted string.
 */
function objectToString(value) {
  return _objectToString_nativeObjectToString.call(value);
}

/* harmony default export */ const _objectToString = (objectToString);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseGetTag.js




/** `Object#toString` result references. */
var nullTag = '[object Null]',
    undefinedTag = '[object Undefined]';

/** Built-in value references. */
var _baseGetTag_symToStringTag = _Symbol/* default */.Z ? _Symbol/* default */.Z.toStringTag : undefined;

/**
 * The base implementation of `getTag` without fallbacks for buggy environments.
 *
 * @private
 * @param {*} value The value to query.
 * @returns {string} Returns the `toStringTag`.
 */
function baseGetTag(value) {
  if (value == null) {
    return value === undefined ? undefinedTag : nullTag;
  }
  return (_baseGetTag_symToStringTag && _baseGetTag_symToStringTag in Object(value))
    ? _getRawTag(value)
    : _objectToString(value);
}

/* harmony default export */ const _baseGetTag = (baseGetTag);


/***/ }),

/***/ 45934:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ _baseKeys)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/_isPrototype.js
var _isPrototype = __webpack_require__(89418);
// EXTERNAL MODULE: ../node_modules/lodash-es/_overArg.js
var _overArg = __webpack_require__(4883);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_nativeKeys.js


/* Built-in method references for those with the same name as other `lodash` methods. */
var nativeKeys = (0,_overArg/* default */.Z)(Object.keys, Object);

/* harmony default export */ const _nativeKeys = (nativeKeys);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseKeys.js



/** Used for built-in method references. */
var objectProto = Object.prototype;

/** Used to check objects for own properties. */
var _baseKeys_hasOwnProperty = objectProto.hasOwnProperty;

/**
 * The base implementation of `_.keys` which doesn't treat sparse arrays as dense.
 *
 * @private
 * @param {Object} object The object to query.
 * @returns {Array} Returns the array of property names.
 */
function baseKeys(object) {
  if (!(0,_isPrototype/* default */.Z)(object)) {
    return _nativeKeys(object);
  }
  var result = [];
  for (var key in Object(object)) {
    if (_baseKeys_hasOwnProperty.call(object, key) && key != 'constructor') {
      result.push(key);
    }
  }
  return result;
}

/* harmony default export */ const _baseKeys = (baseKeys);


/***/ }),

/***/ 99719:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _identity_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(64056);
/* harmony import */ var _overRest_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(15829);
/* harmony import */ var _setToString_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(71649);




/**
 * The base implementation of `_.rest` which doesn't validate or coerce arguments.
 *
 * @private
 * @param {Function} func The function to apply a rest parameter to.
 * @param {number} [start=func.length-1] The start position of the rest parameter.
 * @returns {Function} Returns the new function.
 */
function baseRest(func, start) {
  return (0,_setToString_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)((0,_overRest_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z)(func, start, _identity_js__WEBPACK_IMPORTED_MODULE_2__/* ["default"] */ .Z), func + '');
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (baseRest);


/***/ }),

/***/ 20274:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/**
 * The base implementation of `_.unary` without support for storing metadata.
 *
 * @private
 * @param {Function} func The function to cap arguments for.
 * @returns {Function} Returns the new capped function.
 */
function baseUnary(func) {
  return function(value) {
    return func(value);
  };
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (baseUnary);


/***/ }),

/***/ 52049:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _Uint8Array_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(41049);


/**
 * Creates a clone of `arrayBuffer`.
 *
 * @private
 * @param {ArrayBuffer} arrayBuffer The array buffer to clone.
 * @returns {ArrayBuffer} Returns the cloned array buffer.
 */
function cloneArrayBuffer(arrayBuffer) {
  var result = new arrayBuffer.constructor(arrayBuffer.byteLength);
  new _Uint8Array_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z(result).set(new _Uint8Array_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z(arrayBuffer));
  return result;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (cloneArrayBuffer);


/***/ }),

/***/ 64405:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _root_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(94311);
/* module decorator */ module = __webpack_require__.hmd(module);


/** Detect free variable `exports`. */
var freeExports = typeof exports == 'object' && exports && !exports.nodeType && exports;

/** Detect free variable `module`. */
var freeModule = freeExports && "object" == 'object' && module && !module.nodeType && module;

/** Detect the popular CommonJS extension `module.exports`. */
var moduleExports = freeModule && freeModule.exports === freeExports;

/** Built-in value references. */
var Buffer = moduleExports ? _root_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z.Buffer : undefined,
    allocUnsafe = Buffer ? Buffer.allocUnsafe : undefined;

/**
 * Creates a clone of  `buffer`.
 *
 * @private
 * @param {Buffer} buffer The buffer to clone.
 * @param {boolean} [isDeep] Specify a deep clone.
 * @returns {Buffer} Returns the cloned buffer.
 */
function cloneBuffer(buffer, isDeep) {
  if (isDeep) {
    return buffer.slice();
  }
  var length = buffer.length,
      result = allocUnsafe ? allocUnsafe(length) : new buffer.constructor(length);

  buffer.copy(result);
  return result;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (cloneBuffer);


/***/ }),

/***/ 61601:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _cloneArrayBuffer_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(52049);


/**
 * Creates a clone of `typedArray`.
 *
 * @private
 * @param {Object} typedArray The typed array to clone.
 * @param {boolean} [isDeep] Specify a deep clone.
 * @returns {Object} Returns the cloned typed array.
 */
function cloneTypedArray(typedArray, isDeep) {
  var buffer = isDeep ? (0,_cloneArrayBuffer_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(typedArray.buffer) : typedArray.buffer;
  return new typedArray.constructor(buffer, typedArray.byteOffset, typedArray.length);
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (cloneTypedArray);


/***/ }),

/***/ 93580:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/**
 * Copies the values of `source` to `array`.
 *
 * @private
 * @param {Array} source The array to copy values from.
 * @param {Array} [array=[]] The array to copy values to.
 * @returns {Array} Returns `array`.
 */
function copyArray(source, array) {
  var index = -1,
      length = source.length;

  array || (array = Array(length));
  while (++index < length) {
    array[index] = source[index];
  }
  return array;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (copyArray);


/***/ }),

/***/ 47313:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _assignValue_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(15561);
/* harmony import */ var _baseAssignValue_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(93586);



/**
 * Copies properties of `source` to `object`.
 *
 * @private
 * @param {Object} source The object to copy properties from.
 * @param {Array} props The property identifiers to copy.
 * @param {Object} [object={}] The object to copy properties to.
 * @param {Function} [customizer] The function to customize copied values.
 * @returns {Object} Returns `object`.
 */
function copyObject(source, props, object, customizer) {
  var isNew = !object;
  object || (object = {});

  var index = -1,
      length = props.length;

  while (++index < length) {
    var key = props[index];

    var newValue = customizer
      ? customizer(object[key], source[key], key, object, source)
      : undefined;

    if (newValue === undefined) {
      newValue = source[key];
    }
    if (isNew) {
      (0,_baseAssignValue_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(object, key, newValue);
    } else {
      (0,_assignValue_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z)(object, key, newValue);
    }
  }
  return object;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (copyObject);


/***/ }),

/***/ 40690:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _baseRest_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(99719);
/* harmony import */ var _isIterateeCall_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(47952);



/**
 * Creates a function like `_.assign`.
 *
 * @private
 * @param {Function} assigner The function to assign values.
 * @returns {Function} Returns the new assigner function.
 */
function createAssigner(assigner) {
  return (0,_baseRest_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(function(object, sources) {
    var index = -1,
        length = sources.length,
        customizer = length > 1 ? sources[length - 1] : undefined,
        guard = length > 2 ? sources[2] : undefined;

    customizer = (assigner.length > 3 && typeof customizer == 'function')
      ? (length--, customizer)
      : undefined;

    if (guard && (0,_isIterateeCall_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z)(sources[0], sources[1], guard)) {
      customizer = length < 3 ? undefined : customizer;
      length = 1;
    }
    object = Object(object);
    while (++index < length) {
      var source = sources[index];
      if (source) {
        assigner(object, source, index, customizer);
      }
    }
    return object;
  });
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (createAssigner);


/***/ }),

/***/ 30253:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _getNative_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(26266);


var defineProperty = (function() {
  try {
    var func = (0,_getNative_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(Object, 'defineProperty');
    func({}, '', {});
    return func;
  } catch (e) {}
}());

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (defineProperty);


/***/ }),

/***/ 89268:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/** Detect free variable `global` from Node.js. */
var freeGlobal = typeof __webpack_require__.g == 'object' && __webpack_require__.g && __webpack_require__.g.Object === Object && __webpack_require__.g;

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (freeGlobal);


/***/ }),

/***/ 26266:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ _getNative)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/isFunction.js
var isFunction = __webpack_require__(48489);
// EXTERNAL MODULE: ../node_modules/lodash-es/_root.js
var _root = __webpack_require__(94311);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_coreJsData.js


/** Used to detect overreaching core-js shims. */
var coreJsData = _root/* default */.Z['__core-js_shared__'];

/* harmony default export */ const _coreJsData = (coreJsData);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_isMasked.js


/** Used to detect methods masquerading as native. */
var maskSrcKey = (function() {
  var uid = /[^.]+$/.exec(_coreJsData && _coreJsData.keys && _coreJsData.keys.IE_PROTO || '');
  return uid ? ('Symbol(src)_1.' + uid) : '';
}());

/**
 * Checks if `func` has its source masked.
 *
 * @private
 * @param {Function} func The function to check.
 * @returns {boolean} Returns `true` if `func` is masked, else `false`.
 */
function isMasked(func) {
  return !!maskSrcKey && (maskSrcKey in func);
}

/* harmony default export */ const _isMasked = (isMasked);

// EXTERNAL MODULE: ../node_modules/lodash-es/isObject.js
var isObject = __webpack_require__(60417);
// EXTERNAL MODULE: ../node_modules/lodash-es/_toSource.js
var _toSource = __webpack_require__(62392);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseIsNative.js





/**
 * Used to match `RegExp`
 * [syntax characters](http://ecma-international.org/ecma-262/7.0/#sec-patterns).
 */
var reRegExpChar = /[\\^$.*+?()[\]{}|]/g;

/** Used to detect host constructors (Safari). */
var reIsHostCtor = /^\[object .+?Constructor\]$/;

/** Used for built-in method references. */
var funcProto = Function.prototype,
    objectProto = Object.prototype;

/** Used to resolve the decompiled source of functions. */
var funcToString = funcProto.toString;

/** Used to check objects for own properties. */
var _baseIsNative_hasOwnProperty = objectProto.hasOwnProperty;

/** Used to detect if a method is native. */
var reIsNative = RegExp('^' +
  funcToString.call(_baseIsNative_hasOwnProperty).replace(reRegExpChar, '\\$&')
  .replace(/hasOwnProperty|(function).*?(?=\\\()| for .+?(?=\\\])/g, '$1.*?') + '$'
);

/**
 * The base implementation of `_.isNative` without bad shim checks.
 *
 * @private
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is a native function,
 *  else `false`.
 */
function baseIsNative(value) {
  if (!(0,isObject/* default */.Z)(value) || _isMasked(value)) {
    return false;
  }
  var pattern = (0,isFunction/* default */.Z)(value) ? reIsNative : reIsHostCtor;
  return pattern.test((0,_toSource/* default */.Z)(value));
}

/* harmony default export */ const _baseIsNative = (baseIsNative);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_getValue.js
/**
 * Gets the value at `key` of `object`.
 *
 * @private
 * @param {Object} [object] The object to query.
 * @param {string} key The key of the property to get.
 * @returns {*} Returns the property value.
 */
function getValue(object, key) {
  return object == null ? undefined : object[key];
}

/* harmony default export */ const _getValue = (getValue);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_getNative.js



/**
 * Gets the native function at `key` of `object`.
 *
 * @private
 * @param {Object} object The object to query.
 * @param {string} key The key of the method to get.
 * @returns {*} Returns the function if it's native, else `undefined`.
 */
function getNative(object, key) {
  var value = _getValue(object, key);
  return _baseIsNative(value) ? value : undefined;
}

/* harmony default export */ const _getNative = (getNative);


/***/ }),

/***/ 60895:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _overArg_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(4883);


/** Built-in value references. */
var getPrototype = (0,_overArg_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(Object.getPrototypeOf, Object);

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (getPrototype);


/***/ }),

/***/ 41182:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ _getTag)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/_getNative.js + 4 modules
var _getNative = __webpack_require__(26266);
// EXTERNAL MODULE: ../node_modules/lodash-es/_root.js
var _root = __webpack_require__(94311);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_DataView.js



/* Built-in method references that are verified to be native. */
var DataView = (0,_getNative/* default */.Z)(_root/* default */.Z, 'DataView');

/* harmony default export */ const _DataView = (DataView);

// EXTERNAL MODULE: ../node_modules/lodash-es/_Map.js
var _Map = __webpack_require__(81700);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_Promise.js



/* Built-in method references that are verified to be native. */
var Promise = (0,_getNative/* default */.Z)(_root/* default */.Z, 'Promise');

/* harmony default export */ const _Promise = (Promise);

// EXTERNAL MODULE: ../node_modules/lodash-es/_Set.js
var _Set = __webpack_require__(16889);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_WeakMap.js



/* Built-in method references that are verified to be native. */
var WeakMap = (0,_getNative/* default */.Z)(_root/* default */.Z, 'WeakMap');

/* harmony default export */ const _WeakMap = (WeakMap);

// EXTERNAL MODULE: ../node_modules/lodash-es/_baseGetTag.js + 2 modules
var _baseGetTag = __webpack_require__(77070);
// EXTERNAL MODULE: ../node_modules/lodash-es/_toSource.js
var _toSource = __webpack_require__(62392);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_getTag.js








/** `Object#toString` result references. */
var mapTag = '[object Map]',
    objectTag = '[object Object]',
    promiseTag = '[object Promise]',
    setTag = '[object Set]',
    weakMapTag = '[object WeakMap]';

var dataViewTag = '[object DataView]';

/** Used to detect maps, sets, and weakmaps. */
var dataViewCtorString = (0,_toSource/* default */.Z)(_DataView),
    mapCtorString = (0,_toSource/* default */.Z)(_Map/* default */.Z),
    promiseCtorString = (0,_toSource/* default */.Z)(_Promise),
    setCtorString = (0,_toSource/* default */.Z)(_Set/* default */.Z),
    weakMapCtorString = (0,_toSource/* default */.Z)(_WeakMap);

/**
 * Gets the `toStringTag` of `value`.
 *
 * @private
 * @param {*} value The value to query.
 * @returns {string} Returns the `toStringTag`.
 */
var getTag = _baseGetTag/* default */.Z;

// Fallback for data views, maps, sets, and weak maps in IE 11 and promises in Node.js < 6.
if ((_DataView && getTag(new _DataView(new ArrayBuffer(1))) != dataViewTag) ||
    (_Map/* default */.Z && getTag(new _Map/* default */.Z) != mapTag) ||
    (_Promise && getTag(_Promise.resolve()) != promiseTag) ||
    (_Set/* default */.Z && getTag(new _Set/* default */.Z) != setTag) ||
    (_WeakMap && getTag(new _WeakMap) != weakMapTag)) {
  getTag = function(value) {
    var result = (0,_baseGetTag/* default */.Z)(value),
        Ctor = result == objectTag ? value.constructor : undefined,
        ctorString = Ctor ? (0,_toSource/* default */.Z)(Ctor) : '';

    if (ctorString) {
      switch (ctorString) {
        case dataViewCtorString: return dataViewTag;
        case mapCtorString: return mapTag;
        case promiseCtorString: return promiseTag;
        case setCtorString: return setTag;
        case weakMapCtorString: return weakMapTag;
      }
    }
    return result;
  };
}

/* harmony default export */ const _getTag = (getTag);


/***/ }),

/***/ 95764:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ _initCloneObject)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/isObject.js
var isObject = __webpack_require__(60417);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseCreate.js


/** Built-in value references. */
var objectCreate = Object.create;

/**
 * The base implementation of `_.create` without support for assigning
 * properties to the created object.
 *
 * @private
 * @param {Object} proto The object to inherit from.
 * @returns {Object} Returns the new object.
 */
var baseCreate = (function() {
  function object() {}
  return function(proto) {
    if (!(0,isObject/* default */.Z)(proto)) {
      return {};
    }
    if (objectCreate) {
      return objectCreate(proto);
    }
    object.prototype = proto;
    var result = new object;
    object.prototype = undefined;
    return result;
  };
}());

/* harmony default export */ const _baseCreate = (baseCreate);

// EXTERNAL MODULE: ../node_modules/lodash-es/_getPrototype.js
var _getPrototype = __webpack_require__(60895);
// EXTERNAL MODULE: ../node_modules/lodash-es/_isPrototype.js
var _isPrototype = __webpack_require__(89418);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_initCloneObject.js




/**
 * Initializes an object clone.
 *
 * @private
 * @param {Object} object The object to clone.
 * @returns {Object} Returns the initialized clone.
 */
function initCloneObject(object) {
  return (typeof object.constructor == 'function' && !(0,_isPrototype/* default */.Z)(object))
    ? _baseCreate((0,_getPrototype/* default */.Z)(object))
    : {};
}

/* harmony default export */ const _initCloneObject = (initCloneObject);


/***/ }),

/***/ 8616:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/** Used as references for various `Number` constants. */
var MAX_SAFE_INTEGER = 9007199254740991;

/** Used to detect unsigned integer values. */
var reIsUint = /^(?:0|[1-9]\d*)$/;

/**
 * Checks if `value` is a valid array-like index.
 *
 * @private
 * @param {*} value The value to check.
 * @param {number} [length=MAX_SAFE_INTEGER] The upper bounds of a valid index.
 * @returns {boolean} Returns `true` if `value` is a valid index, else `false`.
 */
function isIndex(value, length) {
  var type = typeof value;
  length = length == null ? MAX_SAFE_INTEGER : length;

  return !!length &&
    (type == 'number' ||
      (type != 'symbol' && reIsUint.test(value))) &&
        (value > -1 && value % 1 == 0 && value < length);
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (isIndex);


/***/ }),

/***/ 47952:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _eq_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(35050);
/* harmony import */ var _isArrayLike_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(69959);
/* harmony import */ var _isIndex_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(8616);
/* harmony import */ var _isObject_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(60417);





/**
 * Checks if the given arguments are from an iteratee call.
 *
 * @private
 * @param {*} value The potential iteratee value argument.
 * @param {*} index The potential iteratee index or key argument.
 * @param {*} object The potential iteratee object argument.
 * @returns {boolean} Returns `true` if the arguments are from an iteratee call,
 *  else `false`.
 */
function isIterateeCall(value, index, object) {
  if (!(0,_isObject_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(object)) {
    return false;
  }
  var type = typeof index;
  if (type == 'number'
        ? ((0,_isArrayLike_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z)(object) && (0,_isIndex_js__WEBPACK_IMPORTED_MODULE_2__/* ["default"] */ .Z)(index, object.length))
        : (type == 'string' && index in object)
      ) {
    return (0,_eq_js__WEBPACK_IMPORTED_MODULE_3__/* ["default"] */ .Z)(object[index], value);
  }
  return false;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (isIterateeCall);


/***/ }),

/***/ 89418:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/** Used for built-in method references. */
var objectProto = Object.prototype;

/**
 * Checks if `value` is likely a prototype object.
 *
 * @private
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is a prototype, else `false`.
 */
function isPrototype(value) {
  var Ctor = value && value.constructor,
      proto = (typeof Ctor == 'function' && Ctor.prototype) || objectProto;

  return value === proto;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (isPrototype);


/***/ }),

/***/ 53594:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _freeGlobal_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(89268);
/* module decorator */ module = __webpack_require__.hmd(module);


/** Detect free variable `exports`. */
var freeExports = typeof exports == 'object' && exports && !exports.nodeType && exports;

/** Detect free variable `module`. */
var freeModule = freeExports && "object" == 'object' && module && !module.nodeType && module;

/** Detect the popular CommonJS extension `module.exports`. */
var moduleExports = freeModule && freeModule.exports === freeExports;

/** Detect free variable `process` from Node.js. */
var freeProcess = moduleExports && _freeGlobal_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z.process;

/** Used to access faster Node.js helpers. */
var nodeUtil = (function() {
  try {
    // Use `util.types` for Node.js 10+.
    var types = freeModule && freeModule.require && freeModule.require('util').types;

    if (types) {
      return types;
    }

    // Legacy `process.binding('util')` for Node.js < 10.
    return freeProcess && freeProcess.binding && freeProcess.binding('util');
  } catch (e) {}
}());

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (nodeUtil);


/***/ }),

/***/ 4883:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/**
 * Creates a unary function that invokes `func` with its argument transformed.
 *
 * @private
 * @param {Function} func The function to wrap.
 * @param {Function} transform The argument transform.
 * @returns {Function} Returns the new function.
 */
function overArg(func, transform) {
  return function(arg) {
    return func(transform(arg));
  };
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (overArg);


/***/ }),

/***/ 15829:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ _overRest)
});

;// CONCATENATED MODULE: ../node_modules/lodash-es/_apply.js
/**
 * A faster alternative to `Function#apply`, this function invokes `func`
 * with the `this` binding of `thisArg` and the arguments of `args`.
 *
 * @private
 * @param {Function} func The function to invoke.
 * @param {*} thisArg The `this` binding of `func`.
 * @param {Array} args The arguments to invoke `func` with.
 * @returns {*} Returns the result of `func`.
 */
function apply(func, thisArg, args) {
  switch (args.length) {
    case 0: return func.call(thisArg);
    case 1: return func.call(thisArg, args[0]);
    case 2: return func.call(thisArg, args[0], args[1]);
    case 3: return func.call(thisArg, args[0], args[1], args[2]);
  }
  return func.apply(thisArg, args);
}

/* harmony default export */ const _apply = (apply);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_overRest.js


/* Built-in method references for those with the same name as other `lodash` methods. */
var nativeMax = Math.max;

/**
 * A specialized version of `baseRest` which transforms the rest array.
 *
 * @private
 * @param {Function} func The function to apply a rest parameter to.
 * @param {number} [start=func.length-1] The start position of the rest parameter.
 * @param {Function} transform The rest array transform.
 * @returns {Function} Returns the new function.
 */
function overRest(func, start, transform) {
  start = nativeMax(start === undefined ? (func.length - 1) : start, 0);
  return function() {
    var args = arguments,
        index = -1,
        length = nativeMax(args.length - start, 0),
        array = Array(length);

    while (++index < length) {
      array[index] = args[start + index];
    }
    index = -1;
    var otherArgs = Array(start + 1);
    while (++index < start) {
      otherArgs[index] = args[index];
    }
    otherArgs[start] = transform(array);
    return _apply(func, this, otherArgs);
  };
}

/* harmony default export */ const _overRest = (overRest);


/***/ }),

/***/ 94311:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _freeGlobal_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(89268);


/** Detect free variable `self`. */
var freeSelf = typeof self == 'object' && self && self.Object === Object && self;

/** Used as a reference to the global object. */
var root = _freeGlobal_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z || freeSelf || Function('return this')();

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (root);


/***/ }),

/***/ 71649:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ _setToString)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/constant.js
var constant = __webpack_require__(78795);
// EXTERNAL MODULE: ../node_modules/lodash-es/_defineProperty.js
var _defineProperty = __webpack_require__(30253);
// EXTERNAL MODULE: ../node_modules/lodash-es/identity.js
var identity = __webpack_require__(64056);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseSetToString.js




/**
 * The base implementation of `setToString` without support for hot loop shorting.
 *
 * @private
 * @param {Function} func The function to modify.
 * @param {Function} string The `toString` result.
 * @returns {Function} Returns `func`.
 */
var baseSetToString = !_defineProperty/* default */.Z ? identity/* default */.Z : function(func, string) {
  return (0,_defineProperty/* default */.Z)(func, 'toString', {
    'configurable': true,
    'enumerable': false,
    'value': (0,constant/* default */.Z)(string),
    'writable': true
  });
};

/* harmony default export */ const _baseSetToString = (baseSetToString);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_shortOut.js
/** Used to detect hot functions by number of calls within a span of milliseconds. */
var HOT_COUNT = 800,
    HOT_SPAN = 16;

/* Built-in method references for those with the same name as other `lodash` methods. */
var nativeNow = Date.now;

/**
 * Creates a function that'll short out and invoke `identity` instead
 * of `func` when it's called `HOT_COUNT` or more times in `HOT_SPAN`
 * milliseconds.
 *
 * @private
 * @param {Function} func The function to restrict.
 * @returns {Function} Returns the new shortable function.
 */
function shortOut(func) {
  var count = 0,
      lastCalled = 0;

  return function() {
    var stamp = nativeNow(),
        remaining = HOT_SPAN - (stamp - lastCalled);

    lastCalled = stamp;
    if (remaining > 0) {
      if (++count >= HOT_COUNT) {
        return arguments[0];
      }
    } else {
      count = 0;
    }
    return func.apply(undefined, arguments);
  };
}

/* harmony default export */ const _shortOut = (shortOut);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_setToString.js



/**
 * Sets the `toString` method of `func` to return `string`.
 *
 * @private
 * @param {Function} func The function to modify.
 * @param {Function} string The `toString` result.
 * @returns {Function} Returns `func`.
 */
var setToString = _shortOut(_baseSetToString);

/* harmony default export */ const _setToString = (setToString);


/***/ }),

/***/ 62392:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/** Used for built-in method references. */
var funcProto = Function.prototype;

/** Used to resolve the decompiled source of functions. */
var funcToString = funcProto.toString;

/**
 * Converts `func` to its source code.
 *
 * @private
 * @param {Function} func The function to convert.
 * @returns {string} Returns the source code.
 */
function toSource(func) {
  if (func != null) {
    try {
      return funcToString.call(func);
    } catch (e) {}
    try {
      return (func + '');
    } catch (e) {}
  }
  return '';
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (toSource);


/***/ }),

/***/ 78795:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/**
 * Creates a function that returns `value`.
 *
 * @static
 * @memberOf _
 * @since 2.4.0
 * @category Util
 * @param {*} value The value to return from the new function.
 * @returns {Function} Returns the new constant function.
 * @example
 *
 * var objects = _.times(2, _.constant({ 'a': 1 }));
 *
 * console.log(objects);
 * // => [{ 'a': 1 }, { 'a': 1 }]
 *
 * console.log(objects[0] === objects[1]);
 * // => true
 */
function constant(value) {
  return function() {
    return value;
  };
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (constant);


/***/ }),

/***/ 35050:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/**
 * Performs a
 * [`SameValueZero`](http://ecma-international.org/ecma-262/7.0/#sec-samevaluezero)
 * comparison between two values to determine if they are equivalent.
 *
 * @static
 * @memberOf _
 * @since 4.0.0
 * @category Lang
 * @param {*} value The value to compare.
 * @param {*} other The other value to compare.
 * @returns {boolean} Returns `true` if the values are equivalent, else `false`.
 * @example
 *
 * var object = { 'a': 1 };
 * var other = { 'a': 1 };
 *
 * _.eq(object, object);
 * // => true
 *
 * _.eq(object, other);
 * // => false
 *
 * _.eq('a', 'a');
 * // => true
 *
 * _.eq('a', Object('a'));
 * // => false
 *
 * _.eq(NaN, NaN);
 * // => true
 */
function eq(value, other) {
  return value === other || (value !== value && other !== other);
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (eq);


/***/ }),

/***/ 64056:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/**
 * This method returns the first argument it receives.
 *
 * @static
 * @since 0.1.0
 * @memberOf _
 * @category Util
 * @param {*} value Any value.
 * @returns {*} Returns `value`.
 * @example
 *
 * var object = { 'a': 1 };
 *
 * console.log(_.identity(object) === object);
 * // => true
 */
function identity(value) {
  return value;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (identity);


/***/ }),

/***/ 9028:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ lodash_es_isArguments)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/_baseGetTag.js + 2 modules
var _baseGetTag = __webpack_require__(77070);
// EXTERNAL MODULE: ../node_modules/lodash-es/isObjectLike.js
var isObjectLike = __webpack_require__(9615);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseIsArguments.js



/** `Object#toString` result references. */
var argsTag = '[object Arguments]';

/**
 * The base implementation of `_.isArguments`.
 *
 * @private
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is an `arguments` object,
 */
function baseIsArguments(value) {
  return (0,isObjectLike/* default */.Z)(value) && (0,_baseGetTag/* default */.Z)(value) == argsTag;
}

/* harmony default export */ const _baseIsArguments = (baseIsArguments);

;// CONCATENATED MODULE: ../node_modules/lodash-es/isArguments.js



/** Used for built-in method references. */
var objectProto = Object.prototype;

/** Used to check objects for own properties. */
var isArguments_hasOwnProperty = objectProto.hasOwnProperty;

/** Built-in value references. */
var propertyIsEnumerable = objectProto.propertyIsEnumerable;

/**
 * Checks if `value` is likely an `arguments` object.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is an `arguments` object,
 *  else `false`.
 * @example
 *
 * _.isArguments(function() { return arguments; }());
 * // => true
 *
 * _.isArguments([1, 2, 3]);
 * // => false
 */
var isArguments = _baseIsArguments(function() { return arguments; }()) ? _baseIsArguments : function(value) {
  return (0,isObjectLike/* default */.Z)(value) && isArguments_hasOwnProperty.call(value, 'callee') &&
    !propertyIsEnumerable.call(value, 'callee');
};

/* harmony default export */ const lodash_es_isArguments = (isArguments);


/***/ }),

/***/ 64058:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/**
 * Checks if `value` is classified as an `Array` object.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is an array, else `false`.
 * @example
 *
 * _.isArray([1, 2, 3]);
 * // => true
 *
 * _.isArray(document.body.children);
 * // => false
 *
 * _.isArray('abc');
 * // => false
 *
 * _.isArray(_.noop);
 * // => false
 */
var isArray = Array.isArray;

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (isArray);


/***/ }),

/***/ 69959:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _isFunction_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(48489);
/* harmony import */ var _isLength_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(30918);



/**
 * Checks if `value` is array-like. A value is considered array-like if it's
 * not a function and has a `value.length` that's an integer greater than or
 * equal to `0` and less than or equal to `Number.MAX_SAFE_INTEGER`.
 *
 * @static
 * @memberOf _
 * @since 4.0.0
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is array-like, else `false`.
 * @example
 *
 * _.isArrayLike([1, 2, 3]);
 * // => true
 *
 * _.isArrayLike(document.body.children);
 * // => true
 *
 * _.isArrayLike('abc');
 * // => true
 *
 * _.isArrayLike(_.noop);
 * // => false
 */
function isArrayLike(value) {
  return value != null && (0,_isLength_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(value.length) && !(0,_isFunction_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z)(value);
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (isArrayLike);


/***/ }),

/***/ 60492:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _isArrayLike_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(69959);
/* harmony import */ var _isObjectLike_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(9615);



/**
 * This method is like `_.isArrayLike` except that it also checks if `value`
 * is an object.
 *
 * @static
 * @memberOf _
 * @since 4.0.0
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is an array-like object,
 *  else `false`.
 * @example
 *
 * _.isArrayLikeObject([1, 2, 3]);
 * // => true
 *
 * _.isArrayLikeObject(document.body.children);
 * // => true
 *
 * _.isArrayLikeObject('abc');
 * // => false
 *
 * _.isArrayLikeObject(_.noop);
 * // => false
 */
function isArrayLikeObject(value) {
  return (0,_isObjectLike_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(value) && (0,_isArrayLike_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z)(value);
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (isArrayLikeObject);


/***/ }),

/***/ 23230:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ lodash_es_isBuffer)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/_root.js
var _root = __webpack_require__(94311);
;// CONCATENATED MODULE: ../node_modules/lodash-es/stubFalse.js
/**
 * This method returns `false`.
 *
 * @static
 * @memberOf _
 * @since 4.13.0
 * @category Util
 * @returns {boolean} Returns `false`.
 * @example
 *
 * _.times(2, _.stubFalse);
 * // => [false, false]
 */
function stubFalse() {
  return false;
}

/* harmony default export */ const lodash_es_stubFalse = (stubFalse);

;// CONCATENATED MODULE: ../node_modules/lodash-es/isBuffer.js
/* module decorator */ module = __webpack_require__.hmd(module);



/** Detect free variable `exports`. */
var freeExports = typeof exports == 'object' && exports && !exports.nodeType && exports;

/** Detect free variable `module`. */
var freeModule = freeExports && "object" == 'object' && module && !module.nodeType && module;

/** Detect the popular CommonJS extension `module.exports`. */
var moduleExports = freeModule && freeModule.exports === freeExports;

/** Built-in value references. */
var Buffer = moduleExports ? _root/* default */.Z.Buffer : undefined;

/* Built-in method references for those with the same name as other `lodash` methods. */
var nativeIsBuffer = Buffer ? Buffer.isBuffer : undefined;

/**
 * Checks if `value` is a buffer.
 *
 * @static
 * @memberOf _
 * @since 4.3.0
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is a buffer, else `false`.
 * @example
 *
 * _.isBuffer(new Buffer(2));
 * // => true
 *
 * _.isBuffer(new Uint8Array(2));
 * // => false
 */
var isBuffer = nativeIsBuffer || lodash_es_stubFalse;

/* harmony default export */ const lodash_es_isBuffer = (isBuffer);


/***/ }),

/***/ 66400:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _baseKeys_js__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(45934);
/* harmony import */ var _getTag_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(41182);
/* harmony import */ var _isArguments_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(9028);
/* harmony import */ var _isArray_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(64058);
/* harmony import */ var _isArrayLike_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(69959);
/* harmony import */ var _isBuffer_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(23230);
/* harmony import */ var _isPrototype_js__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(89418);
/* harmony import */ var _isTypedArray_js__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(14923);









/** `Object#toString` result references. */
var mapTag = '[object Map]',
    setTag = '[object Set]';

/** Used for built-in method references. */
var objectProto = Object.prototype;

/** Used to check objects for own properties. */
var hasOwnProperty = objectProto.hasOwnProperty;

/**
 * Checks if `value` is an empty object, collection, map, or set.
 *
 * Objects are considered empty if they have no own enumerable string keyed
 * properties.
 *
 * Array-like values such as `arguments` objects, arrays, buffers, strings, or
 * jQuery-like collections are considered empty if they have a `length` of `0`.
 * Similarly, maps and sets are considered empty if they have a `size` of `0`.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is empty, else `false`.
 * @example
 *
 * _.isEmpty(null);
 * // => true
 *
 * _.isEmpty(true);
 * // => true
 *
 * _.isEmpty(1);
 * // => true
 *
 * _.isEmpty([1, 2, 3]);
 * // => false
 *
 * _.isEmpty({ 'a': 1 });
 * // => false
 */
function isEmpty(value) {
  if (value == null) {
    return true;
  }
  if ((0,_isArrayLike_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(value) &&
      ((0,_isArray_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z)(value) || typeof value == 'string' || typeof value.splice == 'function' ||
        (0,_isBuffer_js__WEBPACK_IMPORTED_MODULE_2__/* ["default"] */ .Z)(value) || (0,_isTypedArray_js__WEBPACK_IMPORTED_MODULE_3__/* ["default"] */ .Z)(value) || (0,_isArguments_js__WEBPACK_IMPORTED_MODULE_4__/* ["default"] */ .Z)(value))) {
    return !value.length;
  }
  var tag = (0,_getTag_js__WEBPACK_IMPORTED_MODULE_5__/* ["default"] */ .Z)(value);
  if (tag == mapTag || tag == setTag) {
    return !value.size;
  }
  if ((0,_isPrototype_js__WEBPACK_IMPORTED_MODULE_6__/* ["default"] */ .Z)(value)) {
    return !(0,_baseKeys_js__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z)(value).length;
  }
  for (var key in value) {
    if (hasOwnProperty.call(value, key)) {
      return false;
    }
  }
  return true;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (isEmpty);


/***/ }),

/***/ 48489:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _baseGetTag_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(77070);
/* harmony import */ var _isObject_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(60417);



/** `Object#toString` result references. */
var asyncTag = '[object AsyncFunction]',
    funcTag = '[object Function]',
    genTag = '[object GeneratorFunction]',
    proxyTag = '[object Proxy]';

/**
 * Checks if `value` is classified as a `Function` object.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is a function, else `false`.
 * @example
 *
 * _.isFunction(_);
 * // => true
 *
 * _.isFunction(/abc/);
 * // => false
 */
function isFunction(value) {
  if (!(0,_isObject_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(value)) {
    return false;
  }
  // The use of `Object#toString` avoids issues with the `typeof` operator
  // in Safari 9 which returns 'object' for typed arrays and other constructors.
  var tag = (0,_baseGetTag_js__WEBPACK_IMPORTED_MODULE_1__/* ["default"] */ .Z)(value);
  return tag == funcTag || tag == genTag || tag == asyncTag || tag == proxyTag;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (isFunction);


/***/ }),

/***/ 30918:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/** Used as references for various `Number` constants. */
var MAX_SAFE_INTEGER = 9007199254740991;

/**
 * Checks if `value` is a valid array-like length.
 *
 * **Note:** This method is loosely based on
 * [`ToLength`](http://ecma-international.org/ecma-262/7.0/#sec-tolength).
 *
 * @static
 * @memberOf _
 * @since 4.0.0
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is a valid length, else `false`.
 * @example
 *
 * _.isLength(3);
 * // => true
 *
 * _.isLength(Number.MIN_VALUE);
 * // => false
 *
 * _.isLength(Infinity);
 * // => false
 *
 * _.isLength('3');
 * // => false
 */
function isLength(value) {
  return typeof value == 'number' &&
    value > -1 && value % 1 == 0 && value <= MAX_SAFE_INTEGER;
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (isLength);


/***/ }),

/***/ 60417:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/**
 * Checks if `value` is the
 * [language type](http://www.ecma-international.org/ecma-262/7.0/#sec-ecmascript-language-types)
 * of `Object`. (e.g. arrays, functions, objects, regexes, `new Number(0)`, and `new String('')`)
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is an object, else `false`.
 * @example
 *
 * _.isObject({});
 * // => true
 *
 * _.isObject([1, 2, 3]);
 * // => true
 *
 * _.isObject(_.noop);
 * // => true
 *
 * _.isObject(null);
 * // => false
 */
function isObject(value) {
  var type = typeof value;
  return value != null && (type == 'object' || type == 'function');
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (isObject);


/***/ }),

/***/ 9615:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/**
 * Checks if `value` is object-like. A value is object-like if it's not `null`
 * and has a `typeof` result of "object".
 *
 * @static
 * @memberOf _
 * @since 4.0.0
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is object-like, else `false`.
 * @example
 *
 * _.isObjectLike({});
 * // => true
 *
 * _.isObjectLike([1, 2, 3]);
 * // => true
 *
 * _.isObjectLike(_.noop);
 * // => false
 *
 * _.isObjectLike(null);
 * // => false
 */
function isObjectLike(value) {
  return value != null && typeof value == 'object';
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (isObjectLike);


/***/ }),

/***/ 14923:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ lodash_es_isTypedArray)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/_baseGetTag.js + 2 modules
var _baseGetTag = __webpack_require__(77070);
// EXTERNAL MODULE: ../node_modules/lodash-es/isLength.js
var isLength = __webpack_require__(30918);
// EXTERNAL MODULE: ../node_modules/lodash-es/isObjectLike.js
var isObjectLike = __webpack_require__(9615);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseIsTypedArray.js




/** `Object#toString` result references. */
var argsTag = '[object Arguments]',
    arrayTag = '[object Array]',
    boolTag = '[object Boolean]',
    dateTag = '[object Date]',
    errorTag = '[object Error]',
    funcTag = '[object Function]',
    mapTag = '[object Map]',
    numberTag = '[object Number]',
    objectTag = '[object Object]',
    regexpTag = '[object RegExp]',
    setTag = '[object Set]',
    stringTag = '[object String]',
    weakMapTag = '[object WeakMap]';

var arrayBufferTag = '[object ArrayBuffer]',
    dataViewTag = '[object DataView]',
    float32Tag = '[object Float32Array]',
    float64Tag = '[object Float64Array]',
    int8Tag = '[object Int8Array]',
    int16Tag = '[object Int16Array]',
    int32Tag = '[object Int32Array]',
    uint8Tag = '[object Uint8Array]',
    uint8ClampedTag = '[object Uint8ClampedArray]',
    uint16Tag = '[object Uint16Array]',
    uint32Tag = '[object Uint32Array]';

/** Used to identify `toStringTag` values of typed arrays. */
var typedArrayTags = {};
typedArrayTags[float32Tag] = typedArrayTags[float64Tag] =
typedArrayTags[int8Tag] = typedArrayTags[int16Tag] =
typedArrayTags[int32Tag] = typedArrayTags[uint8Tag] =
typedArrayTags[uint8ClampedTag] = typedArrayTags[uint16Tag] =
typedArrayTags[uint32Tag] = true;
typedArrayTags[argsTag] = typedArrayTags[arrayTag] =
typedArrayTags[arrayBufferTag] = typedArrayTags[boolTag] =
typedArrayTags[dataViewTag] = typedArrayTags[dateTag] =
typedArrayTags[errorTag] = typedArrayTags[funcTag] =
typedArrayTags[mapTag] = typedArrayTags[numberTag] =
typedArrayTags[objectTag] = typedArrayTags[regexpTag] =
typedArrayTags[setTag] = typedArrayTags[stringTag] =
typedArrayTags[weakMapTag] = false;

/**
 * The base implementation of `_.isTypedArray` without Node.js optimizations.
 *
 * @private
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is a typed array, else `false`.
 */
function baseIsTypedArray(value) {
  return (0,isObjectLike/* default */.Z)(value) &&
    (0,isLength/* default */.Z)(value.length) && !!typedArrayTags[(0,_baseGetTag/* default */.Z)(value)];
}

/* harmony default export */ const _baseIsTypedArray = (baseIsTypedArray);

// EXTERNAL MODULE: ../node_modules/lodash-es/_baseUnary.js
var _baseUnary = __webpack_require__(20274);
// EXTERNAL MODULE: ../node_modules/lodash-es/_nodeUtil.js
var _nodeUtil = __webpack_require__(53594);
;// CONCATENATED MODULE: ../node_modules/lodash-es/isTypedArray.js




/* Node.js helper references. */
var nodeIsTypedArray = _nodeUtil/* default */.Z && _nodeUtil/* default */.Z.isTypedArray;

/**
 * Checks if `value` is classified as a typed array.
 *
 * @static
 * @memberOf _
 * @since 3.0.0
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is a typed array, else `false`.
 * @example
 *
 * _.isTypedArray(new Uint8Array);
 * // => true
 *
 * _.isTypedArray([]);
 * // => false
 */
var isTypedArray = nodeIsTypedArray ? (0,_baseUnary/* default */.Z)(nodeIsTypedArray) : _baseIsTypedArray;

/* harmony default export */ const lodash_es_isTypedArray = (isTypedArray);


/***/ }),

/***/ 48441:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ lodash_es_keysIn)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/_arrayLikeKeys.js + 1 modules
var _arrayLikeKeys = __webpack_require__(40709);
// EXTERNAL MODULE: ../node_modules/lodash-es/isObject.js
var isObject = __webpack_require__(60417);
// EXTERNAL MODULE: ../node_modules/lodash-es/_isPrototype.js
var _isPrototype = __webpack_require__(89418);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_nativeKeysIn.js
/**
 * This function is like
 * [`Object.keys`](http://ecma-international.org/ecma-262/7.0/#sec-object.keys)
 * except that it includes inherited enumerable properties.
 *
 * @private
 * @param {Object} object The object to query.
 * @returns {Array} Returns the array of property names.
 */
function nativeKeysIn(object) {
  var result = [];
  if (object != null) {
    for (var key in Object(object)) {
      result.push(key);
    }
  }
  return result;
}

/* harmony default export */ const _nativeKeysIn = (nativeKeysIn);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseKeysIn.js




/** Used for built-in method references. */
var objectProto = Object.prototype;

/** Used to check objects for own properties. */
var _baseKeysIn_hasOwnProperty = objectProto.hasOwnProperty;

/**
 * The base implementation of `_.keysIn` which doesn't treat sparse arrays as dense.
 *
 * @private
 * @param {Object} object The object to query.
 * @returns {Array} Returns the array of property names.
 */
function baseKeysIn(object) {
  if (!(0,isObject/* default */.Z)(object)) {
    return _nativeKeysIn(object);
  }
  var isProto = (0,_isPrototype/* default */.Z)(object),
      result = [];

  for (var key in object) {
    if (!(key == 'constructor' && (isProto || !_baseKeysIn_hasOwnProperty.call(object, key)))) {
      result.push(key);
    }
  }
  return result;
}

/* harmony default export */ const _baseKeysIn = (baseKeysIn);

// EXTERNAL MODULE: ../node_modules/lodash-es/isArrayLike.js
var isArrayLike = __webpack_require__(69959);
;// CONCATENATED MODULE: ../node_modules/lodash-es/keysIn.js




/**
 * Creates an array of the own and inherited enumerable property names of `object`.
 *
 * **Note:** Non-object values are coerced to objects.
 *
 * @static
 * @memberOf _
 * @since 3.0.0
 * @category Object
 * @param {Object} object The object to query.
 * @returns {Array} Returns the array of property names.
 * @example
 *
 * function Foo() {
 *   this.a = 1;
 *   this.b = 2;
 * }
 *
 * Foo.prototype.c = 3;
 *
 * _.keysIn(new Foo);
 * // => ['a', 'b', 'c'] (iteration order is not guaranteed)
 */
function keysIn(object) {
  return (0,isArrayLike/* default */.Z)(object) ? (0,_arrayLikeKeys/* default */.Z)(object, true) : _baseKeysIn(object);
}

/* harmony default export */ const lodash_es_keysIn = (keysIn);


/***/ }),

/***/ 86861:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _MapCache_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(24395);


/** Error message constants. */
var FUNC_ERROR_TEXT = 'Expected a function';

/**
 * Creates a function that memoizes the result of `func`. If `resolver` is
 * provided, it determines the cache key for storing the result based on the
 * arguments provided to the memoized function. By default, the first argument
 * provided to the memoized function is used as the map cache key. The `func`
 * is invoked with the `this` binding of the memoized function.
 *
 * **Note:** The cache is exposed as the `cache` property on the memoized
 * function. Its creation may be customized by replacing the `_.memoize.Cache`
 * constructor with one whose instances implement the
 * [`Map`](http://ecma-international.org/ecma-262/7.0/#sec-properties-of-the-map-prototype-object)
 * method interface of `clear`, `delete`, `get`, `has`, and `set`.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Function
 * @param {Function} func The function to have its output memoized.
 * @param {Function} [resolver] The function to resolve the cache key.
 * @returns {Function} Returns the new memoized function.
 * @example
 *
 * var object = { 'a': 1, 'b': 2 };
 * var other = { 'c': 3, 'd': 4 };
 *
 * var values = _.memoize(_.values);
 * values(object);
 * // => [1, 2]
 *
 * values(other);
 * // => [3, 4]
 *
 * object.a = 2;
 * values(object);
 * // => [1, 2]
 *
 * // Modify the result cache.
 * values.cache.set(object, ['a', 'b']);
 * values(object);
 * // => ['a', 'b']
 *
 * // Replace `_.memoize.Cache`.
 * _.memoize.Cache = WeakMap;
 */
function memoize(func, resolver) {
  if (typeof func != 'function' || (resolver != null && typeof resolver != 'function')) {
    throw new TypeError(FUNC_ERROR_TEXT);
  }
  var memoized = function() {
    var args = arguments,
        key = resolver ? resolver.apply(this, args) : args[0],
        cache = memoized.cache;

    if (cache.has(key)) {
      return cache.get(key);
    }
    var result = func.apply(this, args);
    memoized.cache = cache.set(key, result) || cache;
    return result;
  };
  memoized.cache = new (memoize.Cache || _MapCache_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z);
  return memoized;
}

// Expose `MapCache`.
memoize.Cache = _MapCache_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z;

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (memoize);


/***/ }),

/***/ 20895:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  Z: () => (/* binding */ lodash_es_merge)
});

// EXTERNAL MODULE: ../node_modules/lodash-es/_Stack.js + 5 modules
var _Stack = __webpack_require__(82948);
// EXTERNAL MODULE: ../node_modules/lodash-es/_baseAssignValue.js
var _baseAssignValue = __webpack_require__(93586);
// EXTERNAL MODULE: ../node_modules/lodash-es/eq.js
var eq = __webpack_require__(35050);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_assignMergeValue.js



/**
 * This function is like `assignValue` except that it doesn't assign
 * `undefined` values.
 *
 * @private
 * @param {Object} object The object to modify.
 * @param {string} key The key of the property to assign.
 * @param {*} value The value to assign.
 */
function assignMergeValue(object, key, value) {
  if ((value !== undefined && !(0,eq/* default */.Z)(object[key], value)) ||
      (value === undefined && !(key in object))) {
    (0,_baseAssignValue/* default */.Z)(object, key, value);
  }
}

/* harmony default export */ const _assignMergeValue = (assignMergeValue);

// EXTERNAL MODULE: ../node_modules/lodash-es/_baseFor.js + 1 modules
var _baseFor = __webpack_require__(49399);
// EXTERNAL MODULE: ../node_modules/lodash-es/_cloneBuffer.js
var _cloneBuffer = __webpack_require__(64405);
// EXTERNAL MODULE: ../node_modules/lodash-es/_cloneTypedArray.js
var _cloneTypedArray = __webpack_require__(61601);
// EXTERNAL MODULE: ../node_modules/lodash-es/_copyArray.js
var _copyArray = __webpack_require__(93580);
// EXTERNAL MODULE: ../node_modules/lodash-es/_initCloneObject.js + 1 modules
var _initCloneObject = __webpack_require__(95764);
// EXTERNAL MODULE: ../node_modules/lodash-es/isArguments.js + 1 modules
var isArguments = __webpack_require__(9028);
// EXTERNAL MODULE: ../node_modules/lodash-es/isArray.js
var isArray = __webpack_require__(64058);
// EXTERNAL MODULE: ../node_modules/lodash-es/isArrayLikeObject.js
var isArrayLikeObject = __webpack_require__(60492);
// EXTERNAL MODULE: ../node_modules/lodash-es/isBuffer.js + 1 modules
var isBuffer = __webpack_require__(23230);
// EXTERNAL MODULE: ../node_modules/lodash-es/isFunction.js
var isFunction = __webpack_require__(48489);
// EXTERNAL MODULE: ../node_modules/lodash-es/isObject.js
var isObject = __webpack_require__(60417);
// EXTERNAL MODULE: ../node_modules/lodash-es/_baseGetTag.js + 2 modules
var _baseGetTag = __webpack_require__(77070);
// EXTERNAL MODULE: ../node_modules/lodash-es/_getPrototype.js
var _getPrototype = __webpack_require__(60895);
// EXTERNAL MODULE: ../node_modules/lodash-es/isObjectLike.js
var isObjectLike = __webpack_require__(9615);
;// CONCATENATED MODULE: ../node_modules/lodash-es/isPlainObject.js




/** `Object#toString` result references. */
var objectTag = '[object Object]';

/** Used for built-in method references. */
var funcProto = Function.prototype,
    objectProto = Object.prototype;

/** Used to resolve the decompiled source of functions. */
var funcToString = funcProto.toString;

/** Used to check objects for own properties. */
var isPlainObject_hasOwnProperty = objectProto.hasOwnProperty;

/** Used to infer the `Object` constructor. */
var objectCtorString = funcToString.call(Object);

/**
 * Checks if `value` is a plain object, that is, an object created by the
 * `Object` constructor or one with a `[[Prototype]]` of `null`.
 *
 * @static
 * @memberOf _
 * @since 0.8.0
 * @category Lang
 * @param {*} value The value to check.
 * @returns {boolean} Returns `true` if `value` is a plain object, else `false`.
 * @example
 *
 * function Foo() {
 *   this.a = 1;
 * }
 *
 * _.isPlainObject(new Foo);
 * // => false
 *
 * _.isPlainObject([1, 2, 3]);
 * // => false
 *
 * _.isPlainObject({ 'x': 0, 'y': 0 });
 * // => true
 *
 * _.isPlainObject(Object.create(null));
 * // => true
 */
function isPlainObject(value) {
  if (!(0,isObjectLike/* default */.Z)(value) || (0,_baseGetTag/* default */.Z)(value) != objectTag) {
    return false;
  }
  var proto = (0,_getPrototype/* default */.Z)(value);
  if (proto === null) {
    return true;
  }
  var Ctor = isPlainObject_hasOwnProperty.call(proto, 'constructor') && proto.constructor;
  return typeof Ctor == 'function' && Ctor instanceof Ctor &&
    funcToString.call(Ctor) == objectCtorString;
}

/* harmony default export */ const lodash_es_isPlainObject = (isPlainObject);

// EXTERNAL MODULE: ../node_modules/lodash-es/isTypedArray.js + 1 modules
var isTypedArray = __webpack_require__(14923);
;// CONCATENATED MODULE: ../node_modules/lodash-es/_safeGet.js
/**
 * Gets the value at `key`, unless `key` is "__proto__" or "constructor".
 *
 * @private
 * @param {Object} object The object to query.
 * @param {string} key The key of the property to get.
 * @returns {*} Returns the property value.
 */
function safeGet(object, key) {
  if (key === 'constructor' && typeof object[key] === 'function') {
    return;
  }

  if (key == '__proto__') {
    return;
  }

  return object[key];
}

/* harmony default export */ const _safeGet = (safeGet);

// EXTERNAL MODULE: ../node_modules/lodash-es/_copyObject.js
var _copyObject = __webpack_require__(47313);
// EXTERNAL MODULE: ../node_modules/lodash-es/keysIn.js + 2 modules
var keysIn = __webpack_require__(48441);
;// CONCATENATED MODULE: ../node_modules/lodash-es/toPlainObject.js



/**
 * Converts `value` to a plain object flattening inherited enumerable string
 * keyed properties of `value` to own properties of the plain object.
 *
 * @static
 * @memberOf _
 * @since 3.0.0
 * @category Lang
 * @param {*} value The value to convert.
 * @returns {Object} Returns the converted plain object.
 * @example
 *
 * function Foo() {
 *   this.b = 2;
 * }
 *
 * Foo.prototype.c = 3;
 *
 * _.assign({ 'a': 1 }, new Foo);
 * // => { 'a': 1, 'b': 2 }
 *
 * _.assign({ 'a': 1 }, _.toPlainObject(new Foo));
 * // => { 'a': 1, 'b': 2, 'c': 3 }
 */
function toPlainObject(value) {
  return (0,_copyObject/* default */.Z)(value, (0,keysIn/* default */.Z)(value));
}

/* harmony default export */ const lodash_es_toPlainObject = (toPlainObject);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseMergeDeep.js
















/**
 * A specialized version of `baseMerge` for arrays and objects which performs
 * deep merges and tracks traversed objects enabling objects with circular
 * references to be merged.
 *
 * @private
 * @param {Object} object The destination object.
 * @param {Object} source The source object.
 * @param {string} key The key of the value to merge.
 * @param {number} srcIndex The index of `source`.
 * @param {Function} mergeFunc The function to merge values.
 * @param {Function} [customizer] The function to customize assigned values.
 * @param {Object} [stack] Tracks traversed source values and their merged
 *  counterparts.
 */
function baseMergeDeep(object, source, key, srcIndex, mergeFunc, customizer, stack) {
  var objValue = _safeGet(object, key),
      srcValue = _safeGet(source, key),
      stacked = stack.get(srcValue);

  if (stacked) {
    _assignMergeValue(object, key, stacked);
    return;
  }
  var newValue = customizer
    ? customizer(objValue, srcValue, (key + ''), object, source, stack)
    : undefined;

  var isCommon = newValue === undefined;

  if (isCommon) {
    var isArr = (0,isArray/* default */.Z)(srcValue),
        isBuff = !isArr && (0,isBuffer/* default */.Z)(srcValue),
        isTyped = !isArr && !isBuff && (0,isTypedArray/* default */.Z)(srcValue);

    newValue = srcValue;
    if (isArr || isBuff || isTyped) {
      if ((0,isArray/* default */.Z)(objValue)) {
        newValue = objValue;
      }
      else if ((0,isArrayLikeObject/* default */.Z)(objValue)) {
        newValue = (0,_copyArray/* default */.Z)(objValue);
      }
      else if (isBuff) {
        isCommon = false;
        newValue = (0,_cloneBuffer/* default */.Z)(srcValue, true);
      }
      else if (isTyped) {
        isCommon = false;
        newValue = (0,_cloneTypedArray/* default */.Z)(srcValue, true);
      }
      else {
        newValue = [];
      }
    }
    else if (lodash_es_isPlainObject(srcValue) || (0,isArguments/* default */.Z)(srcValue)) {
      newValue = objValue;
      if ((0,isArguments/* default */.Z)(objValue)) {
        newValue = lodash_es_toPlainObject(objValue);
      }
      else if (!(0,isObject/* default */.Z)(objValue) || (0,isFunction/* default */.Z)(objValue)) {
        newValue = (0,_initCloneObject/* default */.Z)(srcValue);
      }
    }
    else {
      isCommon = false;
    }
  }
  if (isCommon) {
    // Recursively merge objects and arrays (susceptible to call stack limits).
    stack.set(srcValue, newValue);
    mergeFunc(newValue, srcValue, srcIndex, customizer, stack);
    stack['delete'](srcValue);
  }
  _assignMergeValue(object, key, newValue);
}

/* harmony default export */ const _baseMergeDeep = (baseMergeDeep);

;// CONCATENATED MODULE: ../node_modules/lodash-es/_baseMerge.js








/**
 * The base implementation of `_.merge` without support for multiple sources.
 *
 * @private
 * @param {Object} object The destination object.
 * @param {Object} source The source object.
 * @param {number} srcIndex The index of `source`.
 * @param {Function} [customizer] The function to customize merged values.
 * @param {Object} [stack] Tracks traversed source values and their merged
 *  counterparts.
 */
function baseMerge(object, source, srcIndex, customizer, stack) {
  if (object === source) {
    return;
  }
  (0,_baseFor/* default */.Z)(source, function(srcValue, key) {
    stack || (stack = new _Stack/* default */.Z);
    if ((0,isObject/* default */.Z)(srcValue)) {
      _baseMergeDeep(object, source, key, srcIndex, baseMerge, customizer, stack);
    }
    else {
      var newValue = customizer
        ? customizer(_safeGet(object, key), srcValue, (key + ''), object, source, stack)
        : undefined;

      if (newValue === undefined) {
        newValue = srcValue;
      }
      _assignMergeValue(object, key, newValue);
    }
  }, keysIn/* default */.Z);
}

/* harmony default export */ const _baseMerge = (baseMerge);

// EXTERNAL MODULE: ../node_modules/lodash-es/_createAssigner.js
var _createAssigner = __webpack_require__(40690);
;// CONCATENATED MODULE: ../node_modules/lodash-es/merge.js



/**
 * This method is like `_.assign` except that it recursively merges own and
 * inherited enumerable string keyed properties of source objects into the
 * destination object. Source properties that resolve to `undefined` are
 * skipped if a destination value exists. Array and plain object properties
 * are merged recursively. Other objects and value types are overridden by
 * assignment. Source objects are applied from left to right. Subsequent
 * sources overwrite property assignments of previous sources.
 *
 * **Note:** This method mutates `object`.
 *
 * @static
 * @memberOf _
 * @since 0.5.0
 * @category Object
 * @param {Object} object The destination object.
 * @param {...Object} [sources] The source objects.
 * @returns {Object} Returns `object`.
 * @example
 *
 * var object = {
 *   'a': [{ 'b': 2 }, { 'd': 4 }]
 * };
 *
 * var other = {
 *   'a': [{ 'c': 3 }, { 'e': 5 }]
 * };
 *
 * _.merge(object, other);
 * // => { 'a': [{ 'b': 2, 'c': 3 }, { 'd': 4, 'e': 5 }] }
 */
var merge = (0,_createAssigner/* default */.Z)(function(object, source, srcIndex) {
  _baseMerge(object, source, srcIndex);
});

/* harmony default export */ const lodash_es_merge = (merge);


/***/ }),

/***/ 50692:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   X: () => (/* binding */ package_default)
/* harmony export */ });
// package.json
var package_default = {
  name: "mermaid",
  version: "11.12.1",
  description: "Markdown-ish syntax for generating flowcharts, mindmaps, sequence diagrams, class diagrams, gantt charts, git graphs and more.",
  type: "module",
  module: "./dist/mermaid.core.mjs",
  types: "./dist/mermaid.d.ts",
  exports: {
    ".": {
      types: "./dist/mermaid.d.ts",
      import: "./dist/mermaid.core.mjs",
      default: "./dist/mermaid.core.mjs"
    },
    "./*": "./*"
  },
  keywords: [
    "diagram",
    "markdown",
    "flowchart",
    "sequence diagram",
    "gantt",
    "class diagram",
    "git graph",
    "mindmap",
    "packet diagram",
    "c4 diagram",
    "er diagram",
    "pie chart",
    "pie diagram",
    "quadrant chart",
    "requirement diagram",
    "graph"
  ],
  scripts: {
    clean: "rimraf dist",
    dev: "pnpm -w dev",
    "docs:code": "typedoc src/defaultConfig.ts src/config.ts src/mermaid.ts && prettier --write ./src/docs/config/setup",
    "docs:build": "rimraf ../../docs && pnpm docs:code && pnpm docs:spellcheck && tsx scripts/docs.cli.mts",
    "docs:verify": "pnpm docs:code && pnpm docs:spellcheck && tsx scripts/docs.cli.mts --verify",
    "docs:pre:vitepress": "pnpm --filter ./src/docs prefetch && rimraf src/vitepress && pnpm docs:code && tsx scripts/docs.cli.mts --vitepress && pnpm --filter ./src/vitepress install --no-frozen-lockfile --ignore-scripts",
    "docs:build:vitepress": "pnpm docs:pre:vitepress && (cd src/vitepress && pnpm run build) && cpy --flat src/docs/landing/ ./src/vitepress/.vitepress/dist/landing",
    "docs:dev": 'pnpm docs:pre:vitepress && concurrently "pnpm --filter ./src/vitepress dev" "tsx scripts/docs.cli.mts --watch --vitepress"',
    "docs:dev:docker": 'pnpm docs:pre:vitepress && concurrently "pnpm --filter ./src/vitepress dev:docker" "tsx scripts/docs.cli.mts --watch --vitepress"',
    "docs:serve": "pnpm docs:build:vitepress && vitepress serve src/vitepress",
    "docs:spellcheck": 'cspell "src/docs/**/*.md"',
    "docs:release-version": "tsx scripts/update-release-version.mts",
    "docs:verify-version": "tsx scripts/update-release-version.mts --verify",
    "types:build-config": "tsx scripts/create-types-from-json-schema.mts",
    "types:verify-config": "tsx scripts/create-types-from-json-schema.mts --verify",
    checkCircle: "npx madge --circular ./src",
    prepublishOnly: "pnpm docs:verify-version"
  },
  repository: {
    type: "git",
    url: "https://github.com/mermaid-js/mermaid"
  },
  author: "Knut Sveidqvist",
  license: "MIT",
  standard: {
    ignore: [
      "**/parser/*.js",
      "dist/**/*.js",
      "cypress/**/*.js"
    ],
    globals: [
      "page"
    ]
  },
  dependencies: {
    "@braintree/sanitize-url": "^7.1.1",
    "@iconify/utils": "^3.0.1",
    "@mermaid-js/parser": "workspace:^",
    "@types/d3": "^7.4.3",
    cytoscape: "^3.29.3",
    "cytoscape-cose-bilkent": "^4.1.0",
    "cytoscape-fcose": "^2.2.0",
    d3: "^7.9.0",
    "d3-sankey": "^0.12.3",
    "dagre-d3-es": "7.0.13",
    dayjs: "^1.11.18",
    dompurify: "^3.2.5",
    katex: "^0.16.22",
    khroma: "^2.1.0",
    "lodash-es": "^4.17.21",
    marked: "^16.2.1",
    roughjs: "^4.6.6",
    stylis: "^4.3.6",
    "ts-dedent": "^2.2.0",
    uuid: "^11.1.0"
  },
  devDependencies: {
    "@adobe/jsonschema2md": "^8.0.5",
    "@iconify/types": "^2.0.0",
    "@types/cytoscape": "^3.21.9",
    "@types/cytoscape-fcose": "^2.2.4",
    "@types/d3-sankey": "^0.12.4",
    "@types/d3-scale": "^4.0.9",
    "@types/d3-scale-chromatic": "^3.1.0",
    "@types/d3-selection": "^3.0.11",
    "@types/d3-shape": "^3.1.7",
    "@types/jsdom": "^21.1.7",
    "@types/katex": "^0.16.7",
    "@types/lodash-es": "^4.17.12",
    "@types/micromatch": "^4.0.9",
    "@types/stylis": "^4.2.7",
    "@types/uuid": "^10.0.0",
    ajv: "^8.17.1",
    canvas: "^3.1.2",
    chokidar: "3.6.0",
    concurrently: "^9.1.2",
    "csstree-validator": "^4.0.1",
    globby: "^14.1.0",
    jison: "^0.4.18",
    "js-base64": "^3.7.8",
    jsdom: "^26.1.0",
    "json-schema-to-typescript": "^15.0.4",
    micromatch: "^4.0.8",
    "path-browserify": "^1.0.1",
    prettier: "^3.5.3",
    remark: "^15.0.1",
    "remark-frontmatter": "^5.0.0",
    "remark-gfm": "^4.0.1",
    rimraf: "^6.0.1",
    "start-server-and-test": "^2.0.13",
    "type-fest": "^4.35.0",
    typedoc: "^0.28.12",
    "typedoc-plugin-markdown": "^4.8.1",
    typescript: "~5.7.3",
    "unist-util-flatmap": "^1.0.0",
    "unist-util-visit": "^5.0.0",
    vitepress: "^1.6.4",
    "vitepress-plugin-search": "1.0.4-alpha.22"
  },
  files: [
    "dist/",
    "README.md"
  ],
  publishConfig: {
    access: "public"
  }
};




/***/ }),

/***/ 61805:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  cj: () => (/* binding */ UnknownDiagramError),
  XV: () => (/* binding */ addDirective),
  Yc: () => (/* binding */ assignWithDepth_default),
  nH: () => (/* binding */ calculateMathMLDimensions),
  ZH: () => (/* binding */ clear),
  LJ: () => (/* binding */ commonDb_exports),
  SY: () => (/* binding */ common_default),
  v2: () => (/* binding */ configureSvgSize),
  u_: () => (/* binding */ defaultConfig),
  Fy: () => (/* binding */ defaultConfig2),
  vZ: () => (/* binding */ defaultConfig_default),
  Vg: () => (/* binding */ detectType),
  Bf: () => (/* binding */ detectors),
  Zn: () => (/* binding */ directiveRegex),
  ku: () => (/* binding */ evaluate),
  M6: () => (/* binding */ frontMatterRegex),
  Mx: () => (/* binding */ getAccDescription),
  eu: () => (/* binding */ getAccTitle),
  iE: () => (/* binding */ getConfig),
  nV: () => (/* binding */ getConfig2),
  _7: () => (/* binding */ getDiagram),
  cq: () => (/* binding */ getDiagramLoader),
  Kr: () => (/* binding */ getDiagramTitle),
  ZD: () => (/* binding */ getSiteConfig),
  xN: () => (/* binding */ getThemeVariables3),
  Gr: () => (/* binding */ getUrl),
  qO: () => (/* binding */ getUserDefinedConfig),
  l0: () => (/* binding */ hasKatex),
  Vw: () => (/* binding */ lineBreakRegex),
  UO: () => (/* binding */ parseGenericTypes),
  Cq: () => (/* binding */ registerDiagram),
  KO: () => (/* binding */ registerLazyLoadedDiagrams),
  ur: () => (/* binding */ renderKatexSanitized),
  mc: () => (/* binding */ chunk_ABZYJK2D_reset),
  NM: () => (/* binding */ sanitizeDirective),
  oO: () => (/* binding */ sanitizeText),
  uX: () => (/* binding */ sanitizeText3),
  dY: () => (/* binding */ saveConfigFromInitialize),
  U$: () => (/* binding */ setAccDescription),
  GN: () => (/* binding */ setAccTitle),
  v6: () => (/* binding */ setConfig),
  Y4: () => (/* binding */ setConfig2),
  g2: () => (/* binding */ setDiagramTitle),
  Yn: () => (/* binding */ setSiteConfig),
  j7: () => (/* binding */ setupGraphViewbox),
  Rw: () => (/* binding */ setupGraphViewbox2),
  Ee: () => (/* binding */ styles_default),
  _j: () => (/* binding */ themes_default),
  Tb: () => (/* binding */ updateSiteConfig)
});

// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-AGHRB4JF.mjs
var chunk_AGHRB4JF = __webpack_require__(74999);
// EXTERNAL MODULE: ../node_modules/khroma/dist/color/index.js + 4 modules
var dist_color = __webpack_require__(42528);
// EXTERNAL MODULE: ../node_modules/khroma/dist/methods/change.js
var change = __webpack_require__(95323);
;// CONCATENATED MODULE: ../node_modules/khroma/dist/methods/adjust.js
/* IMPORT */


/* MAIN */
const adjust = (color, channels) => {
    const ch = dist_color/* default */.Z.parse(color);
    const changes = {};
    for (const c in channels) {
        if (!channels[c])
            continue;
        changes[c] = ch[c] + channels[c];
    }
    return (0,change/* default */.Z)(color, changes);
};
/* EXPORT */
/* harmony default export */ const methods_adjust = (adjust);

// EXTERNAL MODULE: ../node_modules/khroma/dist/methods/rgba.js
var rgba = __webpack_require__(14728);
;// CONCATENATED MODULE: ../node_modules/khroma/dist/methods/mix.js
/* IMPORT */


/* MAIN */
//SOURCE: https://github.com/sass/dart-sass/blob/7457d2e9e7e623d9844ffd037a070cf32d39c348/lib/src/functions/color.dart#L718-L756
const mix = (color1, color2, weight = 50) => {
    const { r: r1, g: g1, b: b1, a: a1 } = dist_color/* default */.Z.parse(color1);
    const { r: r2, g: g2, b: b2, a: a2 } = dist_color/* default */.Z.parse(color2);
    const weightScale = weight / 100;
    const weightNormalized = (weightScale * 2) - 1;
    const alphaDelta = a1 - a2;
    const weight1combined = ((weightNormalized * alphaDelta) === -1) ? weightNormalized : (weightNormalized + alphaDelta) / (1 + weightNormalized * alphaDelta);
    const weight1 = (weight1combined + 1) / 2;
    const weight2 = 1 - weight1;
    const r = (r1 * weight1) + (r2 * weight2);
    const g = (g1 * weight1) + (g2 * weight2);
    const b = (b1 * weight1) + (b2 * weight2);
    const a = (a1 * weightScale) + (a2 * (1 - weightScale));
    return (0,rgba/* default */.Z)(r, g, b, a);
};
/* EXPORT */
/* harmony default export */ const methods_mix = (mix);

;// CONCATENATED MODULE: ../node_modules/khroma/dist/methods/invert.js
/* IMPORT */


/* MAIN */
const invert = (color, weight = 100) => {
    const inverse = dist_color/* default */.Z.parse(color);
    inverse.r = 255 - inverse.r;
    inverse.g = 255 - inverse.g;
    inverse.b = 255 - inverse.b;
    return methods_mix(inverse, color, weight);
};
/* EXPORT */
/* harmony default export */ const methods_invert = (invert);

// EXTERNAL MODULE: ../node_modules/khroma/dist/methods/darken.js
var darken = __webpack_require__(57838);
// EXTERNAL MODULE: ../node_modules/khroma/dist/methods/lighten.js
var lighten = __webpack_require__(28482);
// EXTERNAL MODULE: ../node_modules/khroma/dist/methods/is_dark.js + 2 modules
var is_dark = __webpack_require__(28186);
// EXTERNAL MODULE: ../node_modules/dompurify/dist/purify.es.mjs
var purify_es = __webpack_require__(4719);
;// CONCATENATED MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-ABZYJK2D.mjs


// src/diagram-api/regexes.ts
var frontMatterRegex = /^-{3}\s*[\n\r](.*?)[\n\r]-{3}\s*[\n\r]+/s;
var directiveRegex = /%{2}{\s*(?:(\w+)\s*:|(\w+))\s*(?:(\w+)|((?:(?!}%{2}).|\r?\n)*))?\s*(?:}%{2})?/gi;
var anyCommentRegex = /\s*%%.*\n/gm;

// src/errors.ts
var UnknownDiagramError = class extends Error {
  static {
    (0,chunk_AGHRB4JF/* __name */.eW)(this, "UnknownDiagramError");
  }
  constructor(message) {
    super(message);
    this.name = "UnknownDiagramError";
  }
};

// src/diagram-api/detectType.ts
var detectors = {};
var detectType = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(text, config2) {
  text = text.replace(frontMatterRegex, "").replace(directiveRegex, "").replace(anyCommentRegex, "\n");
  for (const [key, { detector }] of Object.entries(detectors)) {
    const diagram = detector(text, config2);
    if (diagram) {
      return key;
    }
  }
  throw new UnknownDiagramError(
    `No diagram type detected matching given configuration for text: ${text}`
  );
}, "detectType");
var registerLazyLoadedDiagrams = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((...diagrams2) => {
  for (const { id, detector, loader } of diagrams2) {
    addDetector(id, detector, loader);
  }
}, "registerLazyLoadedDiagrams");
var addDetector = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((key, detector, loader) => {
  if (detectors[key]) {
    chunk_AGHRB4JF/* log */.cM.warn(`Detector with key ${key} already exists. Overwriting.`);
  }
  detectors[key] = { detector, loader };
  chunk_AGHRB4JF/* log */.cM.debug(`Detector with key ${key} added${loader ? " with loader" : ""}`);
}, "addDetector");
var getDiagramLoader = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((key) => {
  return detectors[key].loader;
}, "getDiagramLoader");

// src/assignWithDepth.ts
var assignWithDepth = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((dst, src, { depth = 2, clobber = false } = {}) => {
  const config2 = { depth, clobber };
  if (Array.isArray(src) && !Array.isArray(dst)) {
    src.forEach((s) => assignWithDepth(dst, s, config2));
    return dst;
  } else if (Array.isArray(src) && Array.isArray(dst)) {
    src.forEach((s) => {
      if (!dst.includes(s)) {
        dst.push(s);
      }
    });
    return dst;
  }
  if (dst === void 0 || depth <= 0) {
    if (dst !== void 0 && dst !== null && typeof dst === "object" && typeof src === "object") {
      return Object.assign(dst, src);
    } else {
      return src;
    }
  }
  if (src !== void 0 && typeof dst === "object" && typeof src === "object") {
    Object.keys(src).forEach((key) => {
      if (typeof src[key] === "object" && (dst[key] === void 0 || typeof dst[key] === "object")) {
        if (dst[key] === void 0) {
          dst[key] = Array.isArray(src[key]) ? [] : {};
        }
        dst[key] = assignWithDepth(dst[key], src[key], { depth: depth - 1, clobber });
      } else if (clobber || typeof dst[key] !== "object" && typeof src[key] !== "object") {
        dst[key] = src[key];
      }
    });
  }
  return dst;
}, "assignWithDepth");
var assignWithDepth_default = assignWithDepth;

// src/themes/theme-base.js


// src/themes/erDiagram-oldHardcodedValues.ts
var oldAttributeBackgroundColorOdd = "#ffffff";
var oldAttributeBackgroundColorEven = "#f2f2f2";

// src/themes/theme-helpers.js

var mkBorder = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((col, darkMode) => darkMode ? methods_adjust(col, { s: -40, l: 10 }) : methods_adjust(col, { s: -40, l: -10 }), "mkBorder");

// src/themes/theme-base.js
var Theme = class {
  static {
    (0,chunk_AGHRB4JF/* __name */.eW)(this, "Theme");
  }
  constructor() {
    this.background = "#f4f4f4";
    this.primaryColor = "#fff4dd";
    this.noteBkgColor = "#fff5ad";
    this.noteTextColor = "#333";
    this.THEME_COLOR_LIMIT = 12;
    this.fontFamily = '"trebuchet ms", verdana, arial, sans-serif';
    this.fontSize = "16px";
  }
  updateColors() {
    this.primaryTextColor = this.primaryTextColor || (this.darkMode ? "#eee" : "#333");
    this.secondaryColor = this.secondaryColor || methods_adjust(this.primaryColor, { h: -120 });
    this.tertiaryColor = this.tertiaryColor || methods_adjust(this.primaryColor, { h: 180, l: 5 });
    this.primaryBorderColor = this.primaryBorderColor || mkBorder(this.primaryColor, this.darkMode);
    this.secondaryBorderColor = this.secondaryBorderColor || mkBorder(this.secondaryColor, this.darkMode);
    this.tertiaryBorderColor = this.tertiaryBorderColor || mkBorder(this.tertiaryColor, this.darkMode);
    this.noteBorderColor = this.noteBorderColor || mkBorder(this.noteBkgColor, this.darkMode);
    this.noteBkgColor = this.noteBkgColor || "#fff5ad";
    this.noteTextColor = this.noteTextColor || "#333";
    this.secondaryTextColor = this.secondaryTextColor || methods_invert(this.secondaryColor);
    this.tertiaryTextColor = this.tertiaryTextColor || methods_invert(this.tertiaryColor);
    this.lineColor = this.lineColor || methods_invert(this.background);
    this.arrowheadColor = this.arrowheadColor || methods_invert(this.background);
    this.textColor = this.textColor || this.primaryTextColor;
    this.border2 = this.border2 || this.tertiaryBorderColor;
    this.nodeBkg = this.nodeBkg || this.primaryColor;
    this.mainBkg = this.mainBkg || this.primaryColor;
    this.nodeBorder = this.nodeBorder || this.primaryBorderColor;
    this.clusterBkg = this.clusterBkg || this.tertiaryColor;
    this.clusterBorder = this.clusterBorder || this.tertiaryBorderColor;
    this.defaultLinkColor = this.defaultLinkColor || this.lineColor;
    this.titleColor = this.titleColor || this.tertiaryTextColor;
    this.edgeLabelBackground = this.edgeLabelBackground || (this.darkMode ? (0,darken/* default */.Z)(this.secondaryColor, 30) : this.secondaryColor);
    this.nodeTextColor = this.nodeTextColor || this.primaryTextColor;
    this.actorBorder = this.actorBorder || this.primaryBorderColor;
    this.actorBkg = this.actorBkg || this.mainBkg;
    this.actorTextColor = this.actorTextColor || this.primaryTextColor;
    this.actorLineColor = this.actorLineColor || this.actorBorder;
    this.labelBoxBkgColor = this.labelBoxBkgColor || this.actorBkg;
    this.signalColor = this.signalColor || this.textColor;
    this.signalTextColor = this.signalTextColor || this.textColor;
    this.labelBoxBorderColor = this.labelBoxBorderColor || this.actorBorder;
    this.labelTextColor = this.labelTextColor || this.actorTextColor;
    this.loopTextColor = this.loopTextColor || this.actorTextColor;
    this.activationBorderColor = this.activationBorderColor || (0,darken/* default */.Z)(this.secondaryColor, 10);
    this.activationBkgColor = this.activationBkgColor || this.secondaryColor;
    this.sequenceNumberColor = this.sequenceNumberColor || methods_invert(this.lineColor);
    this.sectionBkgColor = this.sectionBkgColor || this.tertiaryColor;
    this.altSectionBkgColor = this.altSectionBkgColor || "white";
    this.sectionBkgColor = this.sectionBkgColor || this.secondaryColor;
    this.sectionBkgColor2 = this.sectionBkgColor2 || this.primaryColor;
    this.excludeBkgColor = this.excludeBkgColor || "#eeeeee";
    this.taskBorderColor = this.taskBorderColor || this.primaryBorderColor;
    this.taskBkgColor = this.taskBkgColor || this.primaryColor;
    this.activeTaskBorderColor = this.activeTaskBorderColor || this.primaryColor;
    this.activeTaskBkgColor = this.activeTaskBkgColor || (0,lighten/* default */.Z)(this.primaryColor, 23);
    this.gridColor = this.gridColor || "lightgrey";
    this.doneTaskBkgColor = this.doneTaskBkgColor || "lightgrey";
    this.doneTaskBorderColor = this.doneTaskBorderColor || "grey";
    this.critBorderColor = this.critBorderColor || "#ff8888";
    this.critBkgColor = this.critBkgColor || "red";
    this.todayLineColor = this.todayLineColor || "red";
    this.vertLineColor = this.vertLineColor || "navy";
    this.taskTextColor = this.taskTextColor || this.textColor;
    this.taskTextOutsideColor = this.taskTextOutsideColor || this.textColor;
    this.taskTextLightColor = this.taskTextLightColor || this.textColor;
    this.taskTextColor = this.taskTextColor || this.primaryTextColor;
    this.taskTextDarkColor = this.taskTextDarkColor || this.textColor;
    this.taskTextClickableColor = this.taskTextClickableColor || "#003163";
    this.personBorder = this.personBorder || this.primaryBorderColor;
    this.personBkg = this.personBkg || this.mainBkg;
    if (this.darkMode) {
      this.rowOdd = this.rowOdd || (0,darken/* default */.Z)(this.mainBkg, 5) || "#ffffff";
      this.rowEven = this.rowEven || (0,darken/* default */.Z)(this.mainBkg, 10);
    } else {
      this.rowOdd = this.rowOdd || (0,lighten/* default */.Z)(this.mainBkg, 75) || "#ffffff";
      this.rowEven = this.rowEven || (0,lighten/* default */.Z)(this.mainBkg, 5);
    }
    this.transitionColor = this.transitionColor || this.lineColor;
    this.transitionLabelColor = this.transitionLabelColor || this.textColor;
    this.stateLabelColor = this.stateLabelColor || this.stateBkg || this.primaryTextColor;
    this.stateBkg = this.stateBkg || this.mainBkg;
    this.labelBackgroundColor = this.labelBackgroundColor || this.stateBkg;
    this.compositeBackground = this.compositeBackground || this.background || this.tertiaryColor;
    this.altBackground = this.altBackground || this.tertiaryColor;
    this.compositeTitleBackground = this.compositeTitleBackground || this.mainBkg;
    this.compositeBorder = this.compositeBorder || this.nodeBorder;
    this.innerEndBackground = this.nodeBorder;
    this.errorBkgColor = this.errorBkgColor || this.tertiaryColor;
    this.errorTextColor = this.errorTextColor || this.tertiaryTextColor;
    this.transitionColor = this.transitionColor || this.lineColor;
    this.specialStateColor = this.lineColor;
    this.cScale0 = this.cScale0 || this.primaryColor;
    this.cScale1 = this.cScale1 || this.secondaryColor;
    this.cScale2 = this.cScale2 || this.tertiaryColor;
    this.cScale3 = this.cScale3 || methods_adjust(this.primaryColor, { h: 30 });
    this.cScale4 = this.cScale4 || methods_adjust(this.primaryColor, { h: 60 });
    this.cScale5 = this.cScale5 || methods_adjust(this.primaryColor, { h: 90 });
    this.cScale6 = this.cScale6 || methods_adjust(this.primaryColor, { h: 120 });
    this.cScale7 = this.cScale7 || methods_adjust(this.primaryColor, { h: 150 });
    this.cScale8 = this.cScale8 || methods_adjust(this.primaryColor, { h: 210, l: 150 });
    this.cScale9 = this.cScale9 || methods_adjust(this.primaryColor, { h: 270 });
    this.cScale10 = this.cScale10 || methods_adjust(this.primaryColor, { h: 300 });
    this.cScale11 = this.cScale11 || methods_adjust(this.primaryColor, { h: 330 });
    if (this.darkMode) {
      for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
        this["cScale" + i] = (0,darken/* default */.Z)(this["cScale" + i], 75);
      }
    } else {
      for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
        this["cScale" + i] = (0,darken/* default */.Z)(this["cScale" + i], 25);
      }
    }
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["cScaleInv" + i] = this["cScaleInv" + i] || methods_invert(this["cScale" + i]);
    }
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      if (this.darkMode) {
        this["cScalePeer" + i] = this["cScalePeer" + i] || (0,lighten/* default */.Z)(this["cScale" + i], 10);
      } else {
        this["cScalePeer" + i] = this["cScalePeer" + i] || (0,darken/* default */.Z)(this["cScale" + i], 10);
      }
    }
    this.scaleLabelColor = this.scaleLabelColor || this.labelTextColor;
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["cScaleLabel" + i] = this["cScaleLabel" + i] || this.scaleLabelColor;
    }
    const multiplier = this.darkMode ? -4 : -1;
    for (let i = 0; i < 5; i++) {
      this["surface" + i] = this["surface" + i] || methods_adjust(this.mainBkg, { h: 180, s: -15, l: multiplier * (5 + i * 3) });
      this["surfacePeer" + i] = this["surfacePeer" + i] || methods_adjust(this.mainBkg, { h: 180, s: -15, l: multiplier * (8 + i * 3) });
    }
    this.classText = this.classText || this.textColor;
    this.fillType0 = this.fillType0 || this.primaryColor;
    this.fillType1 = this.fillType1 || this.secondaryColor;
    this.fillType2 = this.fillType2 || methods_adjust(this.primaryColor, { h: 64 });
    this.fillType3 = this.fillType3 || methods_adjust(this.secondaryColor, { h: 64 });
    this.fillType4 = this.fillType4 || methods_adjust(this.primaryColor, { h: -64 });
    this.fillType5 = this.fillType5 || methods_adjust(this.secondaryColor, { h: -64 });
    this.fillType6 = this.fillType6 || methods_adjust(this.primaryColor, { h: 128 });
    this.fillType7 = this.fillType7 || methods_adjust(this.secondaryColor, { h: 128 });
    this.pie1 = this.pie1 || this.primaryColor;
    this.pie2 = this.pie2 || this.secondaryColor;
    this.pie3 = this.pie3 || this.tertiaryColor;
    this.pie4 = this.pie4 || methods_adjust(this.primaryColor, { l: -10 });
    this.pie5 = this.pie5 || methods_adjust(this.secondaryColor, { l: -10 });
    this.pie6 = this.pie6 || methods_adjust(this.tertiaryColor, { l: -10 });
    this.pie7 = this.pie7 || methods_adjust(this.primaryColor, { h: 60, l: -10 });
    this.pie8 = this.pie8 || methods_adjust(this.primaryColor, { h: -60, l: -10 });
    this.pie9 = this.pie9 || methods_adjust(this.primaryColor, { h: 120, l: 0 });
    this.pie10 = this.pie10 || methods_adjust(this.primaryColor, { h: 60, l: -20 });
    this.pie11 = this.pie11 || methods_adjust(this.primaryColor, { h: -60, l: -20 });
    this.pie12 = this.pie12 || methods_adjust(this.primaryColor, { h: 120, l: -10 });
    this.pieTitleTextSize = this.pieTitleTextSize || "25px";
    this.pieTitleTextColor = this.pieTitleTextColor || this.taskTextDarkColor;
    this.pieSectionTextSize = this.pieSectionTextSize || "17px";
    this.pieSectionTextColor = this.pieSectionTextColor || this.textColor;
    this.pieLegendTextSize = this.pieLegendTextSize || "17px";
    this.pieLegendTextColor = this.pieLegendTextColor || this.taskTextDarkColor;
    this.pieStrokeColor = this.pieStrokeColor || "black";
    this.pieStrokeWidth = this.pieStrokeWidth || "2px";
    this.pieOuterStrokeWidth = this.pieOuterStrokeWidth || "2px";
    this.pieOuterStrokeColor = this.pieOuterStrokeColor || "black";
    this.pieOpacity = this.pieOpacity || "0.7";
    this.radar = {
      axisColor: this.radar?.axisColor || this.lineColor,
      axisStrokeWidth: this.radar?.axisStrokeWidth || 2,
      axisLabelFontSize: this.radar?.axisLabelFontSize || 12,
      curveOpacity: this.radar?.curveOpacity || 0.5,
      curveStrokeWidth: this.radar?.curveStrokeWidth || 2,
      graticuleColor: this.radar?.graticuleColor || "#DEDEDE",
      graticuleStrokeWidth: this.radar?.graticuleStrokeWidth || 1,
      graticuleOpacity: this.radar?.graticuleOpacity || 0.3,
      legendBoxSize: this.radar?.legendBoxSize || 12,
      legendFontSize: this.radar?.legendFontSize || 12
    };
    this.archEdgeColor = this.archEdgeColor || "#777";
    this.archEdgeArrowColor = this.archEdgeArrowColor || "#777";
    this.archEdgeWidth = this.archEdgeWidth || "3";
    this.archGroupBorderColor = this.archGroupBorderColor || "#000";
    this.archGroupBorderWidth = this.archGroupBorderWidth || "2px";
    this.quadrant1Fill = this.quadrant1Fill || this.primaryColor;
    this.quadrant2Fill = this.quadrant2Fill || methods_adjust(this.primaryColor, { r: 5, g: 5, b: 5 });
    this.quadrant3Fill = this.quadrant3Fill || methods_adjust(this.primaryColor, { r: 10, g: 10, b: 10 });
    this.quadrant4Fill = this.quadrant4Fill || methods_adjust(this.primaryColor, { r: 15, g: 15, b: 15 });
    this.quadrant1TextFill = this.quadrant1TextFill || this.primaryTextColor;
    this.quadrant2TextFill = this.quadrant2TextFill || methods_adjust(this.primaryTextColor, { r: -5, g: -5, b: -5 });
    this.quadrant3TextFill = this.quadrant3TextFill || methods_adjust(this.primaryTextColor, { r: -10, g: -10, b: -10 });
    this.quadrant4TextFill = this.quadrant4TextFill || methods_adjust(this.primaryTextColor, { r: -15, g: -15, b: -15 });
    this.quadrantPointFill = this.quadrantPointFill || (0,is_dark/* default */.Z)(this.quadrant1Fill) ? (0,lighten/* default */.Z)(this.quadrant1Fill) : (0,darken/* default */.Z)(this.quadrant1Fill);
    this.quadrantPointTextFill = this.quadrantPointTextFill || this.primaryTextColor;
    this.quadrantXAxisTextFill = this.quadrantXAxisTextFill || this.primaryTextColor;
    this.quadrantYAxisTextFill = this.quadrantYAxisTextFill || this.primaryTextColor;
    this.quadrantInternalBorderStrokeFill = this.quadrantInternalBorderStrokeFill || this.primaryBorderColor;
    this.quadrantExternalBorderStrokeFill = this.quadrantExternalBorderStrokeFill || this.primaryBorderColor;
    this.quadrantTitleFill = this.quadrantTitleFill || this.primaryTextColor;
    this.xyChart = {
      backgroundColor: this.xyChart?.backgroundColor || this.background,
      titleColor: this.xyChart?.titleColor || this.primaryTextColor,
      xAxisTitleColor: this.xyChart?.xAxisTitleColor || this.primaryTextColor,
      xAxisLabelColor: this.xyChart?.xAxisLabelColor || this.primaryTextColor,
      xAxisTickColor: this.xyChart?.xAxisTickColor || this.primaryTextColor,
      xAxisLineColor: this.xyChart?.xAxisLineColor || this.primaryTextColor,
      yAxisTitleColor: this.xyChart?.yAxisTitleColor || this.primaryTextColor,
      yAxisLabelColor: this.xyChart?.yAxisLabelColor || this.primaryTextColor,
      yAxisTickColor: this.xyChart?.yAxisTickColor || this.primaryTextColor,
      yAxisLineColor: this.xyChart?.yAxisLineColor || this.primaryTextColor,
      plotColorPalette: this.xyChart?.plotColorPalette || "#FFF4DD,#FFD8B1,#FFA07A,#ECEFF1,#D6DBDF,#C3E0A8,#FFB6A4,#FFD74D,#738FA7,#FFFFF0"
    };
    this.requirementBackground = this.requirementBackground || this.primaryColor;
    this.requirementBorderColor = this.requirementBorderColor || this.primaryBorderColor;
    this.requirementBorderSize = this.requirementBorderSize || "1";
    this.requirementTextColor = this.requirementTextColor || this.primaryTextColor;
    this.relationColor = this.relationColor || this.lineColor;
    this.relationLabelBackground = this.relationLabelBackground || (this.darkMode ? (0,darken/* default */.Z)(this.secondaryColor, 30) : this.secondaryColor);
    this.relationLabelColor = this.relationLabelColor || this.actorTextColor;
    this.git0 = this.git0 || this.primaryColor;
    this.git1 = this.git1 || this.secondaryColor;
    this.git2 = this.git2 || this.tertiaryColor;
    this.git3 = this.git3 || methods_adjust(this.primaryColor, { h: -30 });
    this.git4 = this.git4 || methods_adjust(this.primaryColor, { h: -60 });
    this.git5 = this.git5 || methods_adjust(this.primaryColor, { h: -90 });
    this.git6 = this.git6 || methods_adjust(this.primaryColor, { h: 60 });
    this.git7 = this.git7 || methods_adjust(this.primaryColor, { h: 120 });
    if (this.darkMode) {
      this.git0 = (0,lighten/* default */.Z)(this.git0, 25);
      this.git1 = (0,lighten/* default */.Z)(this.git1, 25);
      this.git2 = (0,lighten/* default */.Z)(this.git2, 25);
      this.git3 = (0,lighten/* default */.Z)(this.git3, 25);
      this.git4 = (0,lighten/* default */.Z)(this.git4, 25);
      this.git5 = (0,lighten/* default */.Z)(this.git5, 25);
      this.git6 = (0,lighten/* default */.Z)(this.git6, 25);
      this.git7 = (0,lighten/* default */.Z)(this.git7, 25);
    } else {
      this.git0 = (0,darken/* default */.Z)(this.git0, 25);
      this.git1 = (0,darken/* default */.Z)(this.git1, 25);
      this.git2 = (0,darken/* default */.Z)(this.git2, 25);
      this.git3 = (0,darken/* default */.Z)(this.git3, 25);
      this.git4 = (0,darken/* default */.Z)(this.git4, 25);
      this.git5 = (0,darken/* default */.Z)(this.git5, 25);
      this.git6 = (0,darken/* default */.Z)(this.git6, 25);
      this.git7 = (0,darken/* default */.Z)(this.git7, 25);
    }
    this.gitInv0 = this.gitInv0 || methods_invert(this.git0);
    this.gitInv1 = this.gitInv1 || methods_invert(this.git1);
    this.gitInv2 = this.gitInv2 || methods_invert(this.git2);
    this.gitInv3 = this.gitInv3 || methods_invert(this.git3);
    this.gitInv4 = this.gitInv4 || methods_invert(this.git4);
    this.gitInv5 = this.gitInv5 || methods_invert(this.git5);
    this.gitInv6 = this.gitInv6 || methods_invert(this.git6);
    this.gitInv7 = this.gitInv7 || methods_invert(this.git7);
    this.branchLabelColor = this.branchLabelColor || (this.darkMode ? "black" : this.labelTextColor);
    this.gitBranchLabel0 = this.gitBranchLabel0 || this.branchLabelColor;
    this.gitBranchLabel1 = this.gitBranchLabel1 || this.branchLabelColor;
    this.gitBranchLabel2 = this.gitBranchLabel2 || this.branchLabelColor;
    this.gitBranchLabel3 = this.gitBranchLabel3 || this.branchLabelColor;
    this.gitBranchLabel4 = this.gitBranchLabel4 || this.branchLabelColor;
    this.gitBranchLabel5 = this.gitBranchLabel5 || this.branchLabelColor;
    this.gitBranchLabel6 = this.gitBranchLabel6 || this.branchLabelColor;
    this.gitBranchLabel7 = this.gitBranchLabel7 || this.branchLabelColor;
    this.tagLabelColor = this.tagLabelColor || this.primaryTextColor;
    this.tagLabelBackground = this.tagLabelBackground || this.primaryColor;
    this.tagLabelBorder = this.tagBorder || this.primaryBorderColor;
    this.tagLabelFontSize = this.tagLabelFontSize || "10px";
    this.commitLabelColor = this.commitLabelColor || this.secondaryTextColor;
    this.commitLabelBackground = this.commitLabelBackground || this.secondaryColor;
    this.commitLabelFontSize = this.commitLabelFontSize || "10px";
    this.attributeBackgroundColorOdd = this.attributeBackgroundColorOdd || oldAttributeBackgroundColorOdd;
    this.attributeBackgroundColorEven = this.attributeBackgroundColorEven || oldAttributeBackgroundColorEven;
  }
  calculate(overrides) {
    if (typeof overrides !== "object") {
      this.updateColors();
      return;
    }
    const keys = Object.keys(overrides);
    keys.forEach((k) => {
      this[k] = overrides[k];
    });
    this.updateColors();
    keys.forEach((k) => {
      this[k] = overrides[k];
    });
  }
};
var getThemeVariables = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((userOverrides) => {
  const theme = new Theme();
  theme.calculate(userOverrides);
  return theme;
}, "getThemeVariables");

// src/themes/theme-dark.js

var Theme2 = class {
  static {
    (0,chunk_AGHRB4JF/* __name */.eW)(this, "Theme");
  }
  constructor() {
    this.background = "#333";
    this.primaryColor = "#1f2020";
    this.secondaryColor = (0,lighten/* default */.Z)(this.primaryColor, 16);
    this.tertiaryColor = methods_adjust(this.primaryColor, { h: -160 });
    this.primaryBorderColor = methods_invert(this.background);
    this.secondaryBorderColor = mkBorder(this.secondaryColor, this.darkMode);
    this.tertiaryBorderColor = mkBorder(this.tertiaryColor, this.darkMode);
    this.primaryTextColor = methods_invert(this.primaryColor);
    this.secondaryTextColor = methods_invert(this.secondaryColor);
    this.tertiaryTextColor = methods_invert(this.tertiaryColor);
    this.lineColor = methods_invert(this.background);
    this.textColor = methods_invert(this.background);
    this.mainBkg = "#1f2020";
    this.secondBkg = "calculated";
    this.mainContrastColor = "lightgrey";
    this.darkTextColor = (0,lighten/* default */.Z)(methods_invert("#323D47"), 10);
    this.lineColor = "calculated";
    this.border1 = "#ccc";
    this.border2 = (0,rgba/* default */.Z)(255, 255, 255, 0.25);
    this.arrowheadColor = "calculated";
    this.fontFamily = '"trebuchet ms", verdana, arial, sans-serif';
    this.fontSize = "16px";
    this.labelBackground = "#181818";
    this.textColor = "#ccc";
    this.THEME_COLOR_LIMIT = 12;
    this.nodeBkg = "calculated";
    this.nodeBorder = "calculated";
    this.clusterBkg = "calculated";
    this.clusterBorder = "calculated";
    this.defaultLinkColor = "calculated";
    this.titleColor = "#F9FFFE";
    this.edgeLabelBackground = "calculated";
    this.actorBorder = "calculated";
    this.actorBkg = "calculated";
    this.actorTextColor = "calculated";
    this.actorLineColor = "calculated";
    this.signalColor = "calculated";
    this.signalTextColor = "calculated";
    this.labelBoxBkgColor = "calculated";
    this.labelBoxBorderColor = "calculated";
    this.labelTextColor = "calculated";
    this.loopTextColor = "calculated";
    this.noteBorderColor = "calculated";
    this.noteBkgColor = "#fff5ad";
    this.noteTextColor = "calculated";
    this.activationBorderColor = "calculated";
    this.activationBkgColor = "calculated";
    this.sequenceNumberColor = "black";
    this.sectionBkgColor = (0,darken/* default */.Z)("#EAE8D9", 30);
    this.altSectionBkgColor = "calculated";
    this.sectionBkgColor2 = "#EAE8D9";
    this.excludeBkgColor = (0,darken/* default */.Z)(this.sectionBkgColor, 10);
    this.taskBorderColor = (0,rgba/* default */.Z)(255, 255, 255, 70);
    this.taskBkgColor = "calculated";
    this.taskTextColor = "calculated";
    this.taskTextLightColor = "calculated";
    this.taskTextOutsideColor = "calculated";
    this.taskTextClickableColor = "#003163";
    this.activeTaskBorderColor = (0,rgba/* default */.Z)(255, 255, 255, 50);
    this.activeTaskBkgColor = "#81B1DB";
    this.gridColor = "calculated";
    this.doneTaskBkgColor = "calculated";
    this.doneTaskBorderColor = "grey";
    this.critBorderColor = "#E83737";
    this.critBkgColor = "#E83737";
    this.taskTextDarkColor = "calculated";
    this.todayLineColor = "#DB5757";
    this.vertLineColor = "#00BFFF";
    this.personBorder = this.primaryBorderColor;
    this.personBkg = this.mainBkg;
    this.archEdgeColor = "calculated";
    this.archEdgeArrowColor = "calculated";
    this.archEdgeWidth = "3";
    this.archGroupBorderColor = this.primaryBorderColor;
    this.archGroupBorderWidth = "2px";
    this.rowOdd = this.rowOdd || (0,lighten/* default */.Z)(this.mainBkg, 5) || "#ffffff";
    this.rowEven = this.rowEven || (0,darken/* default */.Z)(this.mainBkg, 10);
    this.labelColor = "calculated";
    this.errorBkgColor = "#a44141";
    this.errorTextColor = "#ddd";
  }
  updateColors() {
    this.secondBkg = (0,lighten/* default */.Z)(this.mainBkg, 16);
    this.lineColor = this.mainContrastColor;
    this.arrowheadColor = this.mainContrastColor;
    this.nodeBkg = this.mainBkg;
    this.nodeBorder = this.border1;
    this.clusterBkg = this.secondBkg;
    this.clusterBorder = this.border2;
    this.defaultLinkColor = this.lineColor;
    this.edgeLabelBackground = (0,lighten/* default */.Z)(this.labelBackground, 25);
    this.actorBorder = this.border1;
    this.actorBkg = this.mainBkg;
    this.actorTextColor = this.mainContrastColor;
    this.actorLineColor = this.actorBorder;
    this.signalColor = this.mainContrastColor;
    this.signalTextColor = this.mainContrastColor;
    this.labelBoxBkgColor = this.actorBkg;
    this.labelBoxBorderColor = this.actorBorder;
    this.labelTextColor = this.mainContrastColor;
    this.loopTextColor = this.mainContrastColor;
    this.noteBorderColor = this.secondaryBorderColor;
    this.noteBkgColor = this.secondBkg;
    this.noteTextColor = this.secondaryTextColor;
    this.activationBorderColor = this.border1;
    this.activationBkgColor = this.secondBkg;
    this.altSectionBkgColor = this.background;
    this.taskBkgColor = (0,lighten/* default */.Z)(this.mainBkg, 23);
    this.taskTextColor = this.darkTextColor;
    this.taskTextLightColor = this.mainContrastColor;
    this.taskTextOutsideColor = this.taskTextLightColor;
    this.gridColor = this.mainContrastColor;
    this.doneTaskBkgColor = this.mainContrastColor;
    this.taskTextDarkColor = this.darkTextColor;
    this.archEdgeColor = this.lineColor;
    this.archEdgeArrowColor = this.lineColor;
    this.transitionColor = this.transitionColor || this.lineColor;
    this.transitionLabelColor = this.transitionLabelColor || this.textColor;
    this.stateLabelColor = this.stateLabelColor || this.stateBkg || this.primaryTextColor;
    this.stateBkg = this.stateBkg || this.mainBkg;
    this.labelBackgroundColor = this.labelBackgroundColor || this.stateBkg;
    this.compositeBackground = this.compositeBackground || this.background || this.tertiaryColor;
    this.altBackground = this.altBackground || "#555";
    this.compositeTitleBackground = this.compositeTitleBackground || this.mainBkg;
    this.compositeBorder = this.compositeBorder || this.nodeBorder;
    this.innerEndBackground = this.primaryBorderColor;
    this.specialStateColor = "#f4f4f4";
    this.errorBkgColor = this.errorBkgColor || this.tertiaryColor;
    this.errorTextColor = this.errorTextColor || this.tertiaryTextColor;
    this.fillType0 = this.primaryColor;
    this.fillType1 = this.secondaryColor;
    this.fillType2 = methods_adjust(this.primaryColor, { h: 64 });
    this.fillType3 = methods_adjust(this.secondaryColor, { h: 64 });
    this.fillType4 = methods_adjust(this.primaryColor, { h: -64 });
    this.fillType5 = methods_adjust(this.secondaryColor, { h: -64 });
    this.fillType6 = methods_adjust(this.primaryColor, { h: 128 });
    this.fillType7 = methods_adjust(this.secondaryColor, { h: 128 });
    this.cScale1 = this.cScale1 || "#0b0000";
    this.cScale2 = this.cScale2 || "#4d1037";
    this.cScale3 = this.cScale3 || "#3f5258";
    this.cScale4 = this.cScale4 || "#4f2f1b";
    this.cScale5 = this.cScale5 || "#6e0a0a";
    this.cScale6 = this.cScale6 || "#3b0048";
    this.cScale7 = this.cScale7 || "#995a01";
    this.cScale8 = this.cScale8 || "#154706";
    this.cScale9 = this.cScale9 || "#161722";
    this.cScale10 = this.cScale10 || "#00296f";
    this.cScale11 = this.cScale11 || "#01629c";
    this.cScale12 = this.cScale12 || "#010029";
    this.cScale0 = this.cScale0 || this.primaryColor;
    this.cScale1 = this.cScale1 || this.secondaryColor;
    this.cScale2 = this.cScale2 || this.tertiaryColor;
    this.cScale3 = this.cScale3 || methods_adjust(this.primaryColor, { h: 30 });
    this.cScale4 = this.cScale4 || methods_adjust(this.primaryColor, { h: 60 });
    this.cScale5 = this.cScale5 || methods_adjust(this.primaryColor, { h: 90 });
    this.cScale6 = this.cScale6 || methods_adjust(this.primaryColor, { h: 120 });
    this.cScale7 = this.cScale7 || methods_adjust(this.primaryColor, { h: 150 });
    this.cScale8 = this.cScale8 || methods_adjust(this.primaryColor, { h: 210 });
    this.cScale9 = this.cScale9 || methods_adjust(this.primaryColor, { h: 270 });
    this.cScale10 = this.cScale10 || methods_adjust(this.primaryColor, { h: 300 });
    this.cScale11 = this.cScale11 || methods_adjust(this.primaryColor, { h: 330 });
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["cScaleInv" + i] = this["cScaleInv" + i] || methods_invert(this["cScale" + i]);
    }
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["cScalePeer" + i] = this["cScalePeer" + i] || (0,lighten/* default */.Z)(this["cScale" + i], 10);
    }
    for (let i = 0; i < 5; i++) {
      this["surface" + i] = this["surface" + i] || methods_adjust(this.mainBkg, { h: 30, s: -30, l: -(-10 + i * 4) });
      this["surfacePeer" + i] = this["surfacePeer" + i] || methods_adjust(this.mainBkg, { h: 30, s: -30, l: -(-7 + i * 4) });
    }
    this.scaleLabelColor = this.scaleLabelColor || (this.darkMode ? "black" : this.labelTextColor);
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["cScaleLabel" + i] = this["cScaleLabel" + i] || this.scaleLabelColor;
    }
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["pie" + i] = this["cScale" + i];
    }
    this.pieTitleTextSize = this.pieTitleTextSize || "25px";
    this.pieTitleTextColor = this.pieTitleTextColor || this.taskTextDarkColor;
    this.pieSectionTextSize = this.pieSectionTextSize || "17px";
    this.pieSectionTextColor = this.pieSectionTextColor || this.textColor;
    this.pieLegendTextSize = this.pieLegendTextSize || "17px";
    this.pieLegendTextColor = this.pieLegendTextColor || this.taskTextDarkColor;
    this.pieStrokeColor = this.pieStrokeColor || "black";
    this.pieStrokeWidth = this.pieStrokeWidth || "2px";
    this.pieOuterStrokeWidth = this.pieOuterStrokeWidth || "2px";
    this.pieOuterStrokeColor = this.pieOuterStrokeColor || "black";
    this.pieOpacity = this.pieOpacity || "0.7";
    this.quadrant1Fill = this.quadrant1Fill || this.primaryColor;
    this.quadrant2Fill = this.quadrant2Fill || methods_adjust(this.primaryColor, { r: 5, g: 5, b: 5 });
    this.quadrant3Fill = this.quadrant3Fill || methods_adjust(this.primaryColor, { r: 10, g: 10, b: 10 });
    this.quadrant4Fill = this.quadrant4Fill || methods_adjust(this.primaryColor, { r: 15, g: 15, b: 15 });
    this.quadrant1TextFill = this.quadrant1TextFill || this.primaryTextColor;
    this.quadrant2TextFill = this.quadrant2TextFill || methods_adjust(this.primaryTextColor, { r: -5, g: -5, b: -5 });
    this.quadrant3TextFill = this.quadrant3TextFill || methods_adjust(this.primaryTextColor, { r: -10, g: -10, b: -10 });
    this.quadrant4TextFill = this.quadrant4TextFill || methods_adjust(this.primaryTextColor, { r: -15, g: -15, b: -15 });
    this.quadrantPointFill = this.quadrantPointFill || (0,is_dark/* default */.Z)(this.quadrant1Fill) ? (0,lighten/* default */.Z)(this.quadrant1Fill) : (0,darken/* default */.Z)(this.quadrant1Fill);
    this.quadrantPointTextFill = this.quadrantPointTextFill || this.primaryTextColor;
    this.quadrantXAxisTextFill = this.quadrantXAxisTextFill || this.primaryTextColor;
    this.quadrantYAxisTextFill = this.quadrantYAxisTextFill || this.primaryTextColor;
    this.quadrantInternalBorderStrokeFill = this.quadrantInternalBorderStrokeFill || this.primaryBorderColor;
    this.quadrantExternalBorderStrokeFill = this.quadrantExternalBorderStrokeFill || this.primaryBorderColor;
    this.quadrantTitleFill = this.quadrantTitleFill || this.primaryTextColor;
    this.xyChart = {
      backgroundColor: this.xyChart?.backgroundColor || this.background,
      titleColor: this.xyChart?.titleColor || this.primaryTextColor,
      xAxisTitleColor: this.xyChart?.xAxisTitleColor || this.primaryTextColor,
      xAxisLabelColor: this.xyChart?.xAxisLabelColor || this.primaryTextColor,
      xAxisTickColor: this.xyChart?.xAxisTickColor || this.primaryTextColor,
      xAxisLineColor: this.xyChart?.xAxisLineColor || this.primaryTextColor,
      yAxisTitleColor: this.xyChart?.yAxisTitleColor || this.primaryTextColor,
      yAxisLabelColor: this.xyChart?.yAxisLabelColor || this.primaryTextColor,
      yAxisTickColor: this.xyChart?.yAxisTickColor || this.primaryTextColor,
      yAxisLineColor: this.xyChart?.yAxisLineColor || this.primaryTextColor,
      plotColorPalette: this.xyChart?.plotColorPalette || "#3498db,#2ecc71,#e74c3c,#f1c40f,#bdc3c7,#ffffff,#34495e,#9b59b6,#1abc9c,#e67e22"
    };
    this.packet = {
      startByteColor: this.primaryTextColor,
      endByteColor: this.primaryTextColor,
      labelColor: this.primaryTextColor,
      titleColor: this.primaryTextColor,
      blockStrokeColor: this.primaryTextColor,
      blockFillColor: this.background
    };
    this.radar = {
      axisColor: this.radar?.axisColor || this.lineColor,
      axisStrokeWidth: this.radar?.axisStrokeWidth || 2,
      axisLabelFontSize: this.radar?.axisLabelFontSize || 12,
      curveOpacity: this.radar?.curveOpacity || 0.5,
      curveStrokeWidth: this.radar?.curveStrokeWidth || 2,
      graticuleColor: this.radar?.graticuleColor || "#DEDEDE",
      graticuleStrokeWidth: this.radar?.graticuleStrokeWidth || 1,
      graticuleOpacity: this.radar?.graticuleOpacity || 0.3,
      legendBoxSize: this.radar?.legendBoxSize || 12,
      legendFontSize: this.radar?.legendFontSize || 12
    };
    this.classText = this.primaryTextColor;
    this.requirementBackground = this.requirementBackground || this.primaryColor;
    this.requirementBorderColor = this.requirementBorderColor || this.primaryBorderColor;
    this.requirementBorderSize = this.requirementBorderSize || "1";
    this.requirementTextColor = this.requirementTextColor || this.primaryTextColor;
    this.relationColor = this.relationColor || this.lineColor;
    this.relationLabelBackground = this.relationLabelBackground || (this.darkMode ? (0,darken/* default */.Z)(this.secondaryColor, 30) : this.secondaryColor);
    this.relationLabelColor = this.relationLabelColor || this.actorTextColor;
    this.git0 = (0,lighten/* default */.Z)(this.secondaryColor, 20);
    this.git1 = (0,lighten/* default */.Z)(this.pie2 || this.secondaryColor, 20);
    this.git2 = (0,lighten/* default */.Z)(this.pie3 || this.tertiaryColor, 20);
    this.git3 = (0,lighten/* default */.Z)(this.pie4 || methods_adjust(this.primaryColor, { h: -30 }), 20);
    this.git4 = (0,lighten/* default */.Z)(this.pie5 || methods_adjust(this.primaryColor, { h: -60 }), 20);
    this.git5 = (0,lighten/* default */.Z)(this.pie6 || methods_adjust(this.primaryColor, { h: -90 }), 10);
    this.git6 = (0,lighten/* default */.Z)(this.pie7 || methods_adjust(this.primaryColor, { h: 60 }), 10);
    this.git7 = (0,lighten/* default */.Z)(this.pie8 || methods_adjust(this.primaryColor, { h: 120 }), 20);
    this.gitInv0 = this.gitInv0 || methods_invert(this.git0);
    this.gitInv1 = this.gitInv1 || methods_invert(this.git1);
    this.gitInv2 = this.gitInv2 || methods_invert(this.git2);
    this.gitInv3 = this.gitInv3 || methods_invert(this.git3);
    this.gitInv4 = this.gitInv4 || methods_invert(this.git4);
    this.gitInv5 = this.gitInv5 || methods_invert(this.git5);
    this.gitInv6 = this.gitInv6 || methods_invert(this.git6);
    this.gitInv7 = this.gitInv7 || methods_invert(this.git7);
    this.gitBranchLabel0 = this.gitBranchLabel0 || methods_invert(this.labelTextColor);
    this.gitBranchLabel1 = this.gitBranchLabel1 || this.labelTextColor;
    this.gitBranchLabel2 = this.gitBranchLabel2 || this.labelTextColor;
    this.gitBranchLabel3 = this.gitBranchLabel3 || methods_invert(this.labelTextColor);
    this.gitBranchLabel4 = this.gitBranchLabel4 || this.labelTextColor;
    this.gitBranchLabel5 = this.gitBranchLabel5 || this.labelTextColor;
    this.gitBranchLabel6 = this.gitBranchLabel6 || this.labelTextColor;
    this.gitBranchLabel7 = this.gitBranchLabel7 || this.labelTextColor;
    this.tagLabelColor = this.tagLabelColor || this.primaryTextColor;
    this.tagLabelBackground = this.tagLabelBackground || this.primaryColor;
    this.tagLabelBorder = this.tagBorder || this.primaryBorderColor;
    this.tagLabelFontSize = this.tagLabelFontSize || "10px";
    this.commitLabelColor = this.commitLabelColor || this.secondaryTextColor;
    this.commitLabelBackground = this.commitLabelBackground || this.secondaryColor;
    this.commitLabelFontSize = this.commitLabelFontSize || "10px";
    this.attributeBackgroundColorOdd = this.attributeBackgroundColorOdd || (0,lighten/* default */.Z)(this.background, 12);
    this.attributeBackgroundColorEven = this.attributeBackgroundColorEven || (0,lighten/* default */.Z)(this.background, 2);
    this.nodeBorder = this.nodeBorder || "#999";
  }
  calculate(overrides) {
    if (typeof overrides !== "object") {
      this.updateColors();
      return;
    }
    const keys = Object.keys(overrides);
    keys.forEach((k) => {
      this[k] = overrides[k];
    });
    this.updateColors();
    keys.forEach((k) => {
      this[k] = overrides[k];
    });
  }
};
var getThemeVariables2 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((userOverrides) => {
  const theme = new Theme2();
  theme.calculate(userOverrides);
  return theme;
}, "getThemeVariables");

// src/themes/theme-default.js

var Theme3 = class {
  static {
    (0,chunk_AGHRB4JF/* __name */.eW)(this, "Theme");
  }
  constructor() {
    this.background = "#f4f4f4";
    this.primaryColor = "#ECECFF";
    this.secondaryColor = methods_adjust(this.primaryColor, { h: 120 });
    this.secondaryColor = "#ffffde";
    this.tertiaryColor = methods_adjust(this.primaryColor, { h: -160 });
    this.primaryBorderColor = mkBorder(this.primaryColor, this.darkMode);
    this.secondaryBorderColor = mkBorder(this.secondaryColor, this.darkMode);
    this.tertiaryBorderColor = mkBorder(this.tertiaryColor, this.darkMode);
    this.primaryTextColor = methods_invert(this.primaryColor);
    this.secondaryTextColor = methods_invert(this.secondaryColor);
    this.tertiaryTextColor = methods_invert(this.tertiaryColor);
    this.lineColor = methods_invert(this.background);
    this.textColor = methods_invert(this.background);
    this.background = "white";
    this.mainBkg = "#ECECFF";
    this.secondBkg = "#ffffde";
    this.lineColor = "#333333";
    this.border1 = "#9370DB";
    this.border2 = "#aaaa33";
    this.arrowheadColor = "#333333";
    this.fontFamily = '"trebuchet ms", verdana, arial, sans-serif';
    this.fontSize = "16px";
    this.labelBackground = "rgba(232,232,232, 0.8)";
    this.textColor = "#333";
    this.THEME_COLOR_LIMIT = 12;
    this.nodeBkg = "calculated";
    this.nodeBorder = "calculated";
    this.clusterBkg = "calculated";
    this.clusterBorder = "calculated";
    this.defaultLinkColor = "calculated";
    this.titleColor = "calculated";
    this.edgeLabelBackground = "calculated";
    this.actorBorder = "calculated";
    this.actorBkg = "calculated";
    this.actorTextColor = "black";
    this.actorLineColor = "calculated";
    this.signalColor = "calculated";
    this.signalTextColor = "calculated";
    this.labelBoxBkgColor = "calculated";
    this.labelBoxBorderColor = "calculated";
    this.labelTextColor = "calculated";
    this.loopTextColor = "calculated";
    this.noteBorderColor = "calculated";
    this.noteBkgColor = "#fff5ad";
    this.noteTextColor = "calculated";
    this.activationBorderColor = "#666";
    this.activationBkgColor = "#f4f4f4";
    this.sequenceNumberColor = "white";
    this.sectionBkgColor = "calculated";
    this.altSectionBkgColor = "calculated";
    this.sectionBkgColor2 = "calculated";
    this.excludeBkgColor = "#eeeeee";
    this.taskBorderColor = "calculated";
    this.taskBkgColor = "calculated";
    this.taskTextLightColor = "calculated";
    this.taskTextColor = this.taskTextLightColor;
    this.taskTextDarkColor = "calculated";
    this.taskTextOutsideColor = this.taskTextDarkColor;
    this.taskTextClickableColor = "calculated";
    this.activeTaskBorderColor = "calculated";
    this.activeTaskBkgColor = "calculated";
    this.gridColor = "calculated";
    this.doneTaskBkgColor = "calculated";
    this.doneTaskBorderColor = "calculated";
    this.critBorderColor = "calculated";
    this.critBkgColor = "calculated";
    this.todayLineColor = "calculated";
    this.vertLineColor = "calculated";
    this.sectionBkgColor = (0,rgba/* default */.Z)(102, 102, 255, 0.49);
    this.altSectionBkgColor = "white";
    this.sectionBkgColor2 = "#fff400";
    this.taskBorderColor = "#534fbc";
    this.taskBkgColor = "#8a90dd";
    this.taskTextLightColor = "white";
    this.taskTextColor = "calculated";
    this.taskTextDarkColor = "black";
    this.taskTextOutsideColor = "calculated";
    this.taskTextClickableColor = "#003163";
    this.activeTaskBorderColor = "#534fbc";
    this.activeTaskBkgColor = "#bfc7ff";
    this.gridColor = "lightgrey";
    this.doneTaskBkgColor = "lightgrey";
    this.doneTaskBorderColor = "grey";
    this.critBorderColor = "#ff8888";
    this.critBkgColor = "red";
    this.todayLineColor = "red";
    this.vertLineColor = "navy";
    this.personBorder = this.primaryBorderColor;
    this.personBkg = this.mainBkg;
    this.archEdgeColor = "calculated";
    this.archEdgeArrowColor = "calculated";
    this.archEdgeWidth = "3";
    this.archGroupBorderColor = this.primaryBorderColor;
    this.archGroupBorderWidth = "2px";
    this.rowOdd = "calculated";
    this.rowEven = "calculated";
    this.labelColor = "black";
    this.errorBkgColor = "#552222";
    this.errorTextColor = "#552222";
    this.updateColors();
  }
  updateColors() {
    this.cScale0 = this.cScale0 || this.primaryColor;
    this.cScale1 = this.cScale1 || this.secondaryColor;
    this.cScale2 = this.cScale2 || this.tertiaryColor;
    this.cScale3 = this.cScale3 || methods_adjust(this.primaryColor, { h: 30 });
    this.cScale4 = this.cScale4 || methods_adjust(this.primaryColor, { h: 60 });
    this.cScale5 = this.cScale5 || methods_adjust(this.primaryColor, { h: 90 });
    this.cScale6 = this.cScale6 || methods_adjust(this.primaryColor, { h: 120 });
    this.cScale7 = this.cScale7 || methods_adjust(this.primaryColor, { h: 150 });
    this.cScale8 = this.cScale8 || methods_adjust(this.primaryColor, { h: 210 });
    this.cScale9 = this.cScale9 || methods_adjust(this.primaryColor, { h: 270 });
    this.cScale10 = this.cScale10 || methods_adjust(this.primaryColor, { h: 300 });
    this.cScale11 = this.cScale11 || methods_adjust(this.primaryColor, { h: 330 });
    this["cScalePeer1"] = this["cScalePeer1"] || (0,darken/* default */.Z)(this.secondaryColor, 45);
    this["cScalePeer2"] = this["cScalePeer2"] || (0,darken/* default */.Z)(this.tertiaryColor, 40);
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["cScale" + i] = (0,darken/* default */.Z)(this["cScale" + i], 10);
      this["cScalePeer" + i] = this["cScalePeer" + i] || (0,darken/* default */.Z)(this["cScale" + i], 25);
    }
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["cScaleInv" + i] = this["cScaleInv" + i] || methods_adjust(this["cScale" + i], { h: 180 });
    }
    for (let i = 0; i < 5; i++) {
      this["surface" + i] = this["surface" + i] || methods_adjust(this.mainBkg, { h: 30, l: -(5 + i * 5) });
      this["surfacePeer" + i] = this["surfacePeer" + i] || methods_adjust(this.mainBkg, { h: 30, l: -(7 + i * 5) });
    }
    this.scaleLabelColor = this.scaleLabelColor !== "calculated" && this.scaleLabelColor ? this.scaleLabelColor : this.labelTextColor;
    if (this.labelTextColor !== "calculated") {
      this.cScaleLabel0 = this.cScaleLabel0 || methods_invert(this.labelTextColor);
      this.cScaleLabel3 = this.cScaleLabel3 || methods_invert(this.labelTextColor);
      for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
        this["cScaleLabel" + i] = this["cScaleLabel" + i] || this.labelTextColor;
      }
    }
    this.nodeBkg = this.mainBkg;
    this.nodeBorder = this.border1;
    this.clusterBkg = this.secondBkg;
    this.clusterBorder = this.border2;
    this.defaultLinkColor = this.lineColor;
    this.titleColor = this.textColor;
    this.edgeLabelBackground = this.labelBackground;
    this.actorBorder = (0,lighten/* default */.Z)(this.border1, 23);
    this.actorBkg = this.mainBkg;
    this.labelBoxBkgColor = this.actorBkg;
    this.signalColor = this.textColor;
    this.signalTextColor = this.textColor;
    this.labelBoxBorderColor = this.actorBorder;
    this.labelTextColor = this.actorTextColor;
    this.loopTextColor = this.actorTextColor;
    this.noteBorderColor = this.border2;
    this.noteTextColor = this.actorTextColor;
    this.actorLineColor = this.actorBorder;
    this.taskTextColor = this.taskTextLightColor;
    this.taskTextOutsideColor = this.taskTextDarkColor;
    this.archEdgeColor = this.lineColor;
    this.archEdgeArrowColor = this.lineColor;
    this.rowOdd = this.rowOdd || (0,lighten/* default */.Z)(this.primaryColor, 75) || "#ffffff";
    this.rowEven = this.rowEven || (0,lighten/* default */.Z)(this.primaryColor, 1);
    this.transitionColor = this.transitionColor || this.lineColor;
    this.transitionLabelColor = this.transitionLabelColor || this.textColor;
    this.stateLabelColor = this.stateLabelColor || this.stateBkg || this.primaryTextColor;
    this.stateBkg = this.stateBkg || this.mainBkg;
    this.labelBackgroundColor = this.labelBackgroundColor || this.stateBkg;
    this.compositeBackground = this.compositeBackground || this.background || this.tertiaryColor;
    this.altBackground = this.altBackground || "#f0f0f0";
    this.compositeTitleBackground = this.compositeTitleBackground || this.mainBkg;
    this.compositeBorder = this.compositeBorder || this.nodeBorder;
    this.innerEndBackground = this.nodeBorder;
    this.specialStateColor = this.lineColor;
    this.errorBkgColor = this.errorBkgColor || this.tertiaryColor;
    this.errorTextColor = this.errorTextColor || this.tertiaryTextColor;
    this.transitionColor = this.transitionColor || this.lineColor;
    this.classText = this.primaryTextColor;
    this.fillType0 = this.primaryColor;
    this.fillType1 = this.secondaryColor;
    this.fillType2 = methods_adjust(this.primaryColor, { h: 64 });
    this.fillType3 = methods_adjust(this.secondaryColor, { h: 64 });
    this.fillType4 = methods_adjust(this.primaryColor, { h: -64 });
    this.fillType5 = methods_adjust(this.secondaryColor, { h: -64 });
    this.fillType6 = methods_adjust(this.primaryColor, { h: 128 });
    this.fillType7 = methods_adjust(this.secondaryColor, { h: 128 });
    this.pie1 = this.pie1 || this.primaryColor;
    this.pie2 = this.pie2 || this.secondaryColor;
    this.pie3 = this.pie3 || methods_adjust(this.tertiaryColor, { l: -40 });
    this.pie4 = this.pie4 || methods_adjust(this.primaryColor, { l: -10 });
    this.pie5 = this.pie5 || methods_adjust(this.secondaryColor, { l: -30 });
    this.pie6 = this.pie6 || methods_adjust(this.tertiaryColor, { l: -20 });
    this.pie7 = this.pie7 || methods_adjust(this.primaryColor, { h: 60, l: -20 });
    this.pie8 = this.pie8 || methods_adjust(this.primaryColor, { h: -60, l: -40 });
    this.pie9 = this.pie9 || methods_adjust(this.primaryColor, { h: 120, l: -40 });
    this.pie10 = this.pie10 || methods_adjust(this.primaryColor, { h: 60, l: -40 });
    this.pie11 = this.pie11 || methods_adjust(this.primaryColor, { h: -90, l: -40 });
    this.pie12 = this.pie12 || methods_adjust(this.primaryColor, { h: 120, l: -30 });
    this.pieTitleTextSize = this.pieTitleTextSize || "25px";
    this.pieTitleTextColor = this.pieTitleTextColor || this.taskTextDarkColor;
    this.pieSectionTextSize = this.pieSectionTextSize || "17px";
    this.pieSectionTextColor = this.pieSectionTextColor || this.textColor;
    this.pieLegendTextSize = this.pieLegendTextSize || "17px";
    this.pieLegendTextColor = this.pieLegendTextColor || this.taskTextDarkColor;
    this.pieStrokeColor = this.pieStrokeColor || "black";
    this.pieStrokeWidth = this.pieStrokeWidth || "2px";
    this.pieOuterStrokeWidth = this.pieOuterStrokeWidth || "2px";
    this.pieOuterStrokeColor = this.pieOuterStrokeColor || "black";
    this.pieOpacity = this.pieOpacity || "0.7";
    this.quadrant1Fill = this.quadrant1Fill || this.primaryColor;
    this.quadrant2Fill = this.quadrant2Fill || methods_adjust(this.primaryColor, { r: 5, g: 5, b: 5 });
    this.quadrant3Fill = this.quadrant3Fill || methods_adjust(this.primaryColor, { r: 10, g: 10, b: 10 });
    this.quadrant4Fill = this.quadrant4Fill || methods_adjust(this.primaryColor, { r: 15, g: 15, b: 15 });
    this.quadrant1TextFill = this.quadrant1TextFill || this.primaryTextColor;
    this.quadrant2TextFill = this.quadrant2TextFill || methods_adjust(this.primaryTextColor, { r: -5, g: -5, b: -5 });
    this.quadrant3TextFill = this.quadrant3TextFill || methods_adjust(this.primaryTextColor, { r: -10, g: -10, b: -10 });
    this.quadrant4TextFill = this.quadrant4TextFill || methods_adjust(this.primaryTextColor, { r: -15, g: -15, b: -15 });
    this.quadrantPointFill = this.quadrantPointFill || (0,is_dark/* default */.Z)(this.quadrant1Fill) ? (0,lighten/* default */.Z)(this.quadrant1Fill) : (0,darken/* default */.Z)(this.quadrant1Fill);
    this.quadrantPointTextFill = this.quadrantPointTextFill || this.primaryTextColor;
    this.quadrantXAxisTextFill = this.quadrantXAxisTextFill || this.primaryTextColor;
    this.quadrantYAxisTextFill = this.quadrantYAxisTextFill || this.primaryTextColor;
    this.quadrantInternalBorderStrokeFill = this.quadrantInternalBorderStrokeFill || this.primaryBorderColor;
    this.quadrantExternalBorderStrokeFill = this.quadrantExternalBorderStrokeFill || this.primaryBorderColor;
    this.quadrantTitleFill = this.quadrantTitleFill || this.primaryTextColor;
    this.radar = {
      axisColor: this.radar?.axisColor || this.lineColor,
      axisStrokeWidth: this.radar?.axisStrokeWidth || 2,
      axisLabelFontSize: this.radar?.axisLabelFontSize || 12,
      curveOpacity: this.radar?.curveOpacity || 0.5,
      curveStrokeWidth: this.radar?.curveStrokeWidth || 2,
      graticuleColor: this.radar?.graticuleColor || "#DEDEDE",
      graticuleStrokeWidth: this.radar?.graticuleStrokeWidth || 1,
      graticuleOpacity: this.radar?.graticuleOpacity || 0.3,
      legendBoxSize: this.radar?.legendBoxSize || 12,
      legendFontSize: this.radar?.legendFontSize || 12
    };
    this.xyChart = {
      backgroundColor: this.xyChart?.backgroundColor || this.background,
      titleColor: this.xyChart?.titleColor || this.primaryTextColor,
      xAxisTitleColor: this.xyChart?.xAxisTitleColor || this.primaryTextColor,
      xAxisLabelColor: this.xyChart?.xAxisLabelColor || this.primaryTextColor,
      xAxisTickColor: this.xyChart?.xAxisTickColor || this.primaryTextColor,
      xAxisLineColor: this.xyChart?.xAxisLineColor || this.primaryTextColor,
      yAxisTitleColor: this.xyChart?.yAxisTitleColor || this.primaryTextColor,
      yAxisLabelColor: this.xyChart?.yAxisLabelColor || this.primaryTextColor,
      yAxisTickColor: this.xyChart?.yAxisTickColor || this.primaryTextColor,
      yAxisLineColor: this.xyChart?.yAxisLineColor || this.primaryTextColor,
      plotColorPalette: this.xyChart?.plotColorPalette || "#ECECFF,#8493A6,#FFC3A0,#DCDDE1,#B8E994,#D1A36F,#C3CDE6,#FFB6C1,#496078,#F8F3E3"
    };
    this.requirementBackground = this.requirementBackground || this.primaryColor;
    this.requirementBorderColor = this.requirementBorderColor || this.primaryBorderColor;
    this.requirementBorderSize = this.requirementBorderSize || "1";
    this.requirementTextColor = this.requirementTextColor || this.primaryTextColor;
    this.relationColor = this.relationColor || this.lineColor;
    this.relationLabelBackground = this.relationLabelBackground || this.labelBackground;
    this.relationLabelColor = this.relationLabelColor || this.actorTextColor;
    this.git0 = this.git0 || this.primaryColor;
    this.git1 = this.git1 || this.secondaryColor;
    this.git2 = this.git2 || this.tertiaryColor;
    this.git3 = this.git3 || methods_adjust(this.primaryColor, { h: -30 });
    this.git4 = this.git4 || methods_adjust(this.primaryColor, { h: -60 });
    this.git5 = this.git5 || methods_adjust(this.primaryColor, { h: -90 });
    this.git6 = this.git6 || methods_adjust(this.primaryColor, { h: 60 });
    this.git7 = this.git7 || methods_adjust(this.primaryColor, { h: 120 });
    if (this.darkMode) {
      this.git0 = (0,lighten/* default */.Z)(this.git0, 25);
      this.git1 = (0,lighten/* default */.Z)(this.git1, 25);
      this.git2 = (0,lighten/* default */.Z)(this.git2, 25);
      this.git3 = (0,lighten/* default */.Z)(this.git3, 25);
      this.git4 = (0,lighten/* default */.Z)(this.git4, 25);
      this.git5 = (0,lighten/* default */.Z)(this.git5, 25);
      this.git6 = (0,lighten/* default */.Z)(this.git6, 25);
      this.git7 = (0,lighten/* default */.Z)(this.git7, 25);
    } else {
      this.git0 = (0,darken/* default */.Z)(this.git0, 25);
      this.git1 = (0,darken/* default */.Z)(this.git1, 25);
      this.git2 = (0,darken/* default */.Z)(this.git2, 25);
      this.git3 = (0,darken/* default */.Z)(this.git3, 25);
      this.git4 = (0,darken/* default */.Z)(this.git4, 25);
      this.git5 = (0,darken/* default */.Z)(this.git5, 25);
      this.git6 = (0,darken/* default */.Z)(this.git6, 25);
      this.git7 = (0,darken/* default */.Z)(this.git7, 25);
    }
    this.gitInv0 = this.gitInv0 || (0,darken/* default */.Z)(methods_invert(this.git0), 25);
    this.gitInv1 = this.gitInv1 || methods_invert(this.git1);
    this.gitInv2 = this.gitInv2 || methods_invert(this.git2);
    this.gitInv3 = this.gitInv3 || methods_invert(this.git3);
    this.gitInv4 = this.gitInv4 || methods_invert(this.git4);
    this.gitInv5 = this.gitInv5 || methods_invert(this.git5);
    this.gitInv6 = this.gitInv6 || methods_invert(this.git6);
    this.gitInv7 = this.gitInv7 || methods_invert(this.git7);
    this.gitBranchLabel0 = this.gitBranchLabel0 || methods_invert(this.labelTextColor);
    this.gitBranchLabel1 = this.gitBranchLabel1 || this.labelTextColor;
    this.gitBranchLabel2 = this.gitBranchLabel2 || this.labelTextColor;
    this.gitBranchLabel3 = this.gitBranchLabel3 || methods_invert(this.labelTextColor);
    this.gitBranchLabel4 = this.gitBranchLabel4 || this.labelTextColor;
    this.gitBranchLabel5 = this.gitBranchLabel5 || this.labelTextColor;
    this.gitBranchLabel6 = this.gitBranchLabel6 || this.labelTextColor;
    this.gitBranchLabel7 = this.gitBranchLabel7 || this.labelTextColor;
    this.tagLabelColor = this.tagLabelColor || this.primaryTextColor;
    this.tagLabelBackground = this.tagLabelBackground || this.primaryColor;
    this.tagLabelBorder = this.tagBorder || this.primaryBorderColor;
    this.tagLabelFontSize = this.tagLabelFontSize || "10px";
    this.commitLabelColor = this.commitLabelColor || this.secondaryTextColor;
    this.commitLabelBackground = this.commitLabelBackground || this.secondaryColor;
    this.commitLabelFontSize = this.commitLabelFontSize || "10px";
    this.attributeBackgroundColorOdd = this.attributeBackgroundColorOdd || oldAttributeBackgroundColorOdd;
    this.attributeBackgroundColorEven = this.attributeBackgroundColorEven || oldAttributeBackgroundColorEven;
  }
  calculate(overrides) {
    Object.keys(this).forEach((k) => {
      if (this[k] === "calculated") {
        this[k] = void 0;
      }
    });
    if (typeof overrides !== "object") {
      this.updateColors();
      return;
    }
    const keys = Object.keys(overrides);
    keys.forEach((k) => {
      this[k] = overrides[k];
    });
    this.updateColors();
    keys.forEach((k) => {
      this[k] = overrides[k];
    });
  }
};
var getThemeVariables3 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((userOverrides) => {
  const theme = new Theme3();
  theme.calculate(userOverrides);
  return theme;
}, "getThemeVariables");

// src/themes/theme-forest.js

var Theme4 = class {
  static {
    (0,chunk_AGHRB4JF/* __name */.eW)(this, "Theme");
  }
  constructor() {
    this.background = "#f4f4f4";
    this.primaryColor = "#cde498";
    this.secondaryColor = "#cdffb2";
    this.background = "white";
    this.mainBkg = "#cde498";
    this.secondBkg = "#cdffb2";
    this.lineColor = "green";
    this.border1 = "#13540c";
    this.border2 = "#6eaa49";
    this.arrowheadColor = "green";
    this.fontFamily = '"trebuchet ms", verdana, arial, sans-serif';
    this.fontSize = "16px";
    this.tertiaryColor = (0,lighten/* default */.Z)("#cde498", 10);
    this.primaryBorderColor = mkBorder(this.primaryColor, this.darkMode);
    this.secondaryBorderColor = mkBorder(this.secondaryColor, this.darkMode);
    this.tertiaryBorderColor = mkBorder(this.tertiaryColor, this.darkMode);
    this.primaryTextColor = methods_invert(this.primaryColor);
    this.secondaryTextColor = methods_invert(this.secondaryColor);
    this.tertiaryTextColor = methods_invert(this.primaryColor);
    this.lineColor = methods_invert(this.background);
    this.textColor = methods_invert(this.background);
    this.THEME_COLOR_LIMIT = 12;
    this.nodeBkg = "calculated";
    this.nodeBorder = "calculated";
    this.clusterBkg = "calculated";
    this.clusterBorder = "calculated";
    this.defaultLinkColor = "calculated";
    this.titleColor = "#333";
    this.edgeLabelBackground = "#e8e8e8";
    this.actorBorder = "calculated";
    this.actorBkg = "calculated";
    this.actorTextColor = "black";
    this.actorLineColor = "calculated";
    this.signalColor = "#333";
    this.signalTextColor = "#333";
    this.labelBoxBkgColor = "calculated";
    this.labelBoxBorderColor = "#326932";
    this.labelTextColor = "calculated";
    this.loopTextColor = "calculated";
    this.noteBorderColor = "calculated";
    this.noteBkgColor = "#fff5ad";
    this.noteTextColor = "calculated";
    this.activationBorderColor = "#666";
    this.activationBkgColor = "#f4f4f4";
    this.sequenceNumberColor = "white";
    this.sectionBkgColor = "#6eaa49";
    this.altSectionBkgColor = "white";
    this.sectionBkgColor2 = "#6eaa49";
    this.excludeBkgColor = "#eeeeee";
    this.taskBorderColor = "calculated";
    this.taskBkgColor = "#487e3a";
    this.taskTextLightColor = "white";
    this.taskTextColor = "calculated";
    this.taskTextDarkColor = "black";
    this.taskTextOutsideColor = "calculated";
    this.taskTextClickableColor = "#003163";
    this.activeTaskBorderColor = "calculated";
    this.activeTaskBkgColor = "calculated";
    this.gridColor = "lightgrey";
    this.doneTaskBkgColor = "lightgrey";
    this.doneTaskBorderColor = "grey";
    this.critBorderColor = "#ff8888";
    this.critBkgColor = "red";
    this.todayLineColor = "red";
    this.vertLineColor = "#00BFFF";
    this.personBorder = this.primaryBorderColor;
    this.personBkg = this.mainBkg;
    this.archEdgeColor = "calculated";
    this.archEdgeArrowColor = "calculated";
    this.archEdgeWidth = "3";
    this.archGroupBorderColor = this.primaryBorderColor;
    this.archGroupBorderWidth = "2px";
    this.labelColor = "black";
    this.errorBkgColor = "#552222";
    this.errorTextColor = "#552222";
  }
  updateColors() {
    this.actorBorder = (0,darken/* default */.Z)(this.mainBkg, 20);
    this.actorBkg = this.mainBkg;
    this.labelBoxBkgColor = this.actorBkg;
    this.labelTextColor = this.actorTextColor;
    this.loopTextColor = this.actorTextColor;
    this.noteBorderColor = this.border2;
    this.noteTextColor = this.actorTextColor;
    this.actorLineColor = this.actorBorder;
    this.cScale0 = this.cScale0 || this.primaryColor;
    this.cScale1 = this.cScale1 || this.secondaryColor;
    this.cScale2 = this.cScale2 || this.tertiaryColor;
    this.cScale3 = this.cScale3 || methods_adjust(this.primaryColor, { h: 30 });
    this.cScale4 = this.cScale4 || methods_adjust(this.primaryColor, { h: 60 });
    this.cScale5 = this.cScale5 || methods_adjust(this.primaryColor, { h: 90 });
    this.cScale6 = this.cScale6 || methods_adjust(this.primaryColor, { h: 120 });
    this.cScale7 = this.cScale7 || methods_adjust(this.primaryColor, { h: 150 });
    this.cScale8 = this.cScale8 || methods_adjust(this.primaryColor, { h: 210 });
    this.cScale9 = this.cScale9 || methods_adjust(this.primaryColor, { h: 270 });
    this.cScale10 = this.cScale10 || methods_adjust(this.primaryColor, { h: 300 });
    this.cScale11 = this.cScale11 || methods_adjust(this.primaryColor, { h: 330 });
    this["cScalePeer1"] = this["cScalePeer1"] || (0,darken/* default */.Z)(this.secondaryColor, 45);
    this["cScalePeer2"] = this["cScalePeer2"] || (0,darken/* default */.Z)(this.tertiaryColor, 40);
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["cScale" + i] = (0,darken/* default */.Z)(this["cScale" + i], 10);
      this["cScalePeer" + i] = this["cScalePeer" + i] || (0,darken/* default */.Z)(this["cScale" + i], 25);
    }
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["cScaleInv" + i] = this["cScaleInv" + i] || methods_adjust(this["cScale" + i], { h: 180 });
    }
    this.scaleLabelColor = this.scaleLabelColor !== "calculated" && this.scaleLabelColor ? this.scaleLabelColor : this.labelTextColor;
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["cScaleLabel" + i] = this["cScaleLabel" + i] || this.scaleLabelColor;
    }
    for (let i = 0; i < 5; i++) {
      this["surface" + i] = this["surface" + i] || methods_adjust(this.mainBkg, { h: 30, s: -30, l: -(5 + i * 5) });
      this["surfacePeer" + i] = this["surfacePeer" + i] || methods_adjust(this.mainBkg, { h: 30, s: -30, l: -(8 + i * 5) });
    }
    this.nodeBkg = this.mainBkg;
    this.nodeBorder = this.border1;
    this.clusterBkg = this.secondBkg;
    this.clusterBorder = this.border2;
    this.defaultLinkColor = this.lineColor;
    this.taskBorderColor = this.border1;
    this.taskTextColor = this.taskTextLightColor;
    this.taskTextOutsideColor = this.taskTextDarkColor;
    this.activeTaskBorderColor = this.taskBorderColor;
    this.activeTaskBkgColor = this.mainBkg;
    this.archEdgeColor = this.lineColor;
    this.archEdgeArrowColor = this.lineColor;
    this.rowOdd = this.rowOdd || (0,lighten/* default */.Z)(this.mainBkg, 75) || "#ffffff";
    this.rowEven = this.rowEven || (0,lighten/* default */.Z)(this.mainBkg, 20);
    this.transitionColor = this.transitionColor || this.lineColor;
    this.transitionLabelColor = this.transitionLabelColor || this.textColor;
    this.stateLabelColor = this.stateLabelColor || this.stateBkg || this.primaryTextColor;
    this.stateBkg = this.stateBkg || this.mainBkg;
    this.labelBackgroundColor = this.labelBackgroundColor || this.stateBkg;
    this.compositeBackground = this.compositeBackground || this.background || this.tertiaryColor;
    this.altBackground = this.altBackground || "#f0f0f0";
    this.compositeTitleBackground = this.compositeTitleBackground || this.mainBkg;
    this.compositeBorder = this.compositeBorder || this.nodeBorder;
    this.innerEndBackground = this.primaryBorderColor;
    this.specialStateColor = this.lineColor;
    this.errorBkgColor = this.errorBkgColor || this.tertiaryColor;
    this.errorTextColor = this.errorTextColor || this.tertiaryTextColor;
    this.transitionColor = this.transitionColor || this.lineColor;
    this.classText = this.primaryTextColor;
    this.fillType0 = this.primaryColor;
    this.fillType1 = this.secondaryColor;
    this.fillType2 = methods_adjust(this.primaryColor, { h: 64 });
    this.fillType3 = methods_adjust(this.secondaryColor, { h: 64 });
    this.fillType4 = methods_adjust(this.primaryColor, { h: -64 });
    this.fillType5 = methods_adjust(this.secondaryColor, { h: -64 });
    this.fillType6 = methods_adjust(this.primaryColor, { h: 128 });
    this.fillType7 = methods_adjust(this.secondaryColor, { h: 128 });
    this.pie1 = this.pie1 || this.primaryColor;
    this.pie2 = this.pie2 || this.secondaryColor;
    this.pie3 = this.pie3 || this.tertiaryColor;
    this.pie4 = this.pie4 || methods_adjust(this.primaryColor, { l: -30 });
    this.pie5 = this.pie5 || methods_adjust(this.secondaryColor, { l: -30 });
    this.pie6 = this.pie6 || methods_adjust(this.tertiaryColor, { h: 40, l: -40 });
    this.pie7 = this.pie7 || methods_adjust(this.primaryColor, { h: 60, l: -10 });
    this.pie8 = this.pie8 || methods_adjust(this.primaryColor, { h: -60, l: -10 });
    this.pie9 = this.pie9 || methods_adjust(this.primaryColor, { h: 120, l: 0 });
    this.pie10 = this.pie10 || methods_adjust(this.primaryColor, { h: 60, l: -50 });
    this.pie11 = this.pie11 || methods_adjust(this.primaryColor, { h: -60, l: -50 });
    this.pie12 = this.pie12 || methods_adjust(this.primaryColor, { h: 120, l: -50 });
    this.pieTitleTextSize = this.pieTitleTextSize || "25px";
    this.pieTitleTextColor = this.pieTitleTextColor || this.taskTextDarkColor;
    this.pieSectionTextSize = this.pieSectionTextSize || "17px";
    this.pieSectionTextColor = this.pieSectionTextColor || this.textColor;
    this.pieLegendTextSize = this.pieLegendTextSize || "17px";
    this.pieLegendTextColor = this.pieLegendTextColor || this.taskTextDarkColor;
    this.pieStrokeColor = this.pieStrokeColor || "black";
    this.pieStrokeWidth = this.pieStrokeWidth || "2px";
    this.pieOuterStrokeWidth = this.pieOuterStrokeWidth || "2px";
    this.pieOuterStrokeColor = this.pieOuterStrokeColor || "black";
    this.pieOpacity = this.pieOpacity || "0.7";
    this.quadrant1Fill = this.quadrant1Fill || this.primaryColor;
    this.quadrant2Fill = this.quadrant2Fill || methods_adjust(this.primaryColor, { r: 5, g: 5, b: 5 });
    this.quadrant3Fill = this.quadrant3Fill || methods_adjust(this.primaryColor, { r: 10, g: 10, b: 10 });
    this.quadrant4Fill = this.quadrant4Fill || methods_adjust(this.primaryColor, { r: 15, g: 15, b: 15 });
    this.quadrant1TextFill = this.quadrant1TextFill || this.primaryTextColor;
    this.quadrant2TextFill = this.quadrant2TextFill || methods_adjust(this.primaryTextColor, { r: -5, g: -5, b: -5 });
    this.quadrant3TextFill = this.quadrant3TextFill || methods_adjust(this.primaryTextColor, { r: -10, g: -10, b: -10 });
    this.quadrant4TextFill = this.quadrant4TextFill || methods_adjust(this.primaryTextColor, { r: -15, g: -15, b: -15 });
    this.quadrantPointFill = this.quadrantPointFill || (0,is_dark/* default */.Z)(this.quadrant1Fill) ? (0,lighten/* default */.Z)(this.quadrant1Fill) : (0,darken/* default */.Z)(this.quadrant1Fill);
    this.quadrantPointTextFill = this.quadrantPointTextFill || this.primaryTextColor;
    this.quadrantXAxisTextFill = this.quadrantXAxisTextFill || this.primaryTextColor;
    this.quadrantYAxisTextFill = this.quadrantYAxisTextFill || this.primaryTextColor;
    this.quadrantInternalBorderStrokeFill = this.quadrantInternalBorderStrokeFill || this.primaryBorderColor;
    this.quadrantExternalBorderStrokeFill = this.quadrantExternalBorderStrokeFill || this.primaryBorderColor;
    this.quadrantTitleFill = this.quadrantTitleFill || this.primaryTextColor;
    this.packet = {
      startByteColor: this.primaryTextColor,
      endByteColor: this.primaryTextColor,
      labelColor: this.primaryTextColor,
      titleColor: this.primaryTextColor,
      blockStrokeColor: this.primaryTextColor,
      blockFillColor: this.mainBkg
    };
    this.radar = {
      axisColor: this.radar?.axisColor || this.lineColor,
      axisStrokeWidth: this.radar?.axisStrokeWidth || 2,
      axisLabelFontSize: this.radar?.axisLabelFontSize || 12,
      curveOpacity: this.radar?.curveOpacity || 0.5,
      curveStrokeWidth: this.radar?.curveStrokeWidth || 2,
      graticuleColor: this.radar?.graticuleColor || "#DEDEDE",
      graticuleStrokeWidth: this.radar?.graticuleStrokeWidth || 1,
      graticuleOpacity: this.radar?.graticuleOpacity || 0.3,
      legendBoxSize: this.radar?.legendBoxSize || 12,
      legendFontSize: this.radar?.legendFontSize || 12
    };
    this.xyChart = {
      backgroundColor: this.xyChart?.backgroundColor || this.background,
      titleColor: this.xyChart?.titleColor || this.primaryTextColor,
      xAxisTitleColor: this.xyChart?.xAxisTitleColor || this.primaryTextColor,
      xAxisLabelColor: this.xyChart?.xAxisLabelColor || this.primaryTextColor,
      xAxisTickColor: this.xyChart?.xAxisTickColor || this.primaryTextColor,
      xAxisLineColor: this.xyChart?.xAxisLineColor || this.primaryTextColor,
      yAxisTitleColor: this.xyChart?.yAxisTitleColor || this.primaryTextColor,
      yAxisLabelColor: this.xyChart?.yAxisLabelColor || this.primaryTextColor,
      yAxisTickColor: this.xyChart?.yAxisTickColor || this.primaryTextColor,
      yAxisLineColor: this.xyChart?.yAxisLineColor || this.primaryTextColor,
      plotColorPalette: this.xyChart?.plotColorPalette || "#CDE498,#FF6B6B,#A0D2DB,#D7BDE2,#F0F0F0,#FFC3A0,#7FD8BE,#FF9A8B,#FAF3E0,#FFF176"
    };
    this.requirementBackground = this.requirementBackground || this.primaryColor;
    this.requirementBorderColor = this.requirementBorderColor || this.primaryBorderColor;
    this.requirementBorderSize = this.requirementBorderSize || "1";
    this.requirementTextColor = this.requirementTextColor || this.primaryTextColor;
    this.relationColor = this.relationColor || this.lineColor;
    this.relationLabelBackground = this.relationLabelBackground || this.edgeLabelBackground;
    this.relationLabelColor = this.relationLabelColor || this.actorTextColor;
    this.git0 = this.git0 || this.primaryColor;
    this.git1 = this.git1 || this.secondaryColor;
    this.git2 = this.git2 || this.tertiaryColor;
    this.git3 = this.git3 || methods_adjust(this.primaryColor, { h: -30 });
    this.git4 = this.git4 || methods_adjust(this.primaryColor, { h: -60 });
    this.git5 = this.git5 || methods_adjust(this.primaryColor, { h: -90 });
    this.git6 = this.git6 || methods_adjust(this.primaryColor, { h: 60 });
    this.git7 = this.git7 || methods_adjust(this.primaryColor, { h: 120 });
    if (this.darkMode) {
      this.git0 = (0,lighten/* default */.Z)(this.git0, 25);
      this.git1 = (0,lighten/* default */.Z)(this.git1, 25);
      this.git2 = (0,lighten/* default */.Z)(this.git2, 25);
      this.git3 = (0,lighten/* default */.Z)(this.git3, 25);
      this.git4 = (0,lighten/* default */.Z)(this.git4, 25);
      this.git5 = (0,lighten/* default */.Z)(this.git5, 25);
      this.git6 = (0,lighten/* default */.Z)(this.git6, 25);
      this.git7 = (0,lighten/* default */.Z)(this.git7, 25);
    } else {
      this.git0 = (0,darken/* default */.Z)(this.git0, 25);
      this.git1 = (0,darken/* default */.Z)(this.git1, 25);
      this.git2 = (0,darken/* default */.Z)(this.git2, 25);
      this.git3 = (0,darken/* default */.Z)(this.git3, 25);
      this.git4 = (0,darken/* default */.Z)(this.git4, 25);
      this.git5 = (0,darken/* default */.Z)(this.git5, 25);
      this.git6 = (0,darken/* default */.Z)(this.git6, 25);
      this.git7 = (0,darken/* default */.Z)(this.git7, 25);
    }
    this.gitInv0 = this.gitInv0 || methods_invert(this.git0);
    this.gitInv1 = this.gitInv1 || methods_invert(this.git1);
    this.gitInv2 = this.gitInv2 || methods_invert(this.git2);
    this.gitInv3 = this.gitInv3 || methods_invert(this.git3);
    this.gitInv4 = this.gitInv4 || methods_invert(this.git4);
    this.gitInv5 = this.gitInv5 || methods_invert(this.git5);
    this.gitInv6 = this.gitInv6 || methods_invert(this.git6);
    this.gitInv7 = this.gitInv7 || methods_invert(this.git7);
    this.gitBranchLabel0 = this.gitBranchLabel0 || methods_invert(this.labelTextColor);
    this.gitBranchLabel1 = this.gitBranchLabel1 || this.labelTextColor;
    this.gitBranchLabel2 = this.gitBranchLabel2 || this.labelTextColor;
    this.gitBranchLabel3 = this.gitBranchLabel3 || methods_invert(this.labelTextColor);
    this.gitBranchLabel4 = this.gitBranchLabel4 || this.labelTextColor;
    this.gitBranchLabel5 = this.gitBranchLabel5 || this.labelTextColor;
    this.gitBranchLabel6 = this.gitBranchLabel6 || this.labelTextColor;
    this.gitBranchLabel7 = this.gitBranchLabel7 || this.labelTextColor;
    this.tagLabelColor = this.tagLabelColor || this.primaryTextColor;
    this.tagLabelBackground = this.tagLabelBackground || this.primaryColor;
    this.tagLabelBorder = this.tagBorder || this.primaryBorderColor;
    this.tagLabelFontSize = this.tagLabelFontSize || "10px";
    this.commitLabelColor = this.commitLabelColor || this.secondaryTextColor;
    this.commitLabelBackground = this.commitLabelBackground || this.secondaryColor;
    this.commitLabelFontSize = this.commitLabelFontSize || "10px";
    this.attributeBackgroundColorOdd = this.attributeBackgroundColorOdd || oldAttributeBackgroundColorOdd;
    this.attributeBackgroundColorEven = this.attributeBackgroundColorEven || oldAttributeBackgroundColorEven;
  }
  calculate(overrides) {
    if (typeof overrides !== "object") {
      this.updateColors();
      return;
    }
    const keys = Object.keys(overrides);
    keys.forEach((k) => {
      this[k] = overrides[k];
    });
    this.updateColors();
    keys.forEach((k) => {
      this[k] = overrides[k];
    });
  }
};
var getThemeVariables4 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((userOverrides) => {
  const theme = new Theme4();
  theme.calculate(userOverrides);
  return theme;
}, "getThemeVariables");

// src/themes/theme-neutral.js

var Theme5 = class {
  static {
    (0,chunk_AGHRB4JF/* __name */.eW)(this, "Theme");
  }
  constructor() {
    this.primaryColor = "#eee";
    this.contrast = "#707070";
    this.secondaryColor = (0,lighten/* default */.Z)(this.contrast, 55);
    this.background = "#ffffff";
    this.tertiaryColor = methods_adjust(this.primaryColor, { h: -160 });
    this.primaryBorderColor = mkBorder(this.primaryColor, this.darkMode);
    this.secondaryBorderColor = mkBorder(this.secondaryColor, this.darkMode);
    this.tertiaryBorderColor = mkBorder(this.tertiaryColor, this.darkMode);
    this.primaryTextColor = methods_invert(this.primaryColor);
    this.secondaryTextColor = methods_invert(this.secondaryColor);
    this.tertiaryTextColor = methods_invert(this.tertiaryColor);
    this.lineColor = methods_invert(this.background);
    this.textColor = methods_invert(this.background);
    this.mainBkg = "#eee";
    this.secondBkg = "calculated";
    this.lineColor = "#666";
    this.border1 = "#999";
    this.border2 = "calculated";
    this.note = "#ffa";
    this.text = "#333";
    this.critical = "#d42";
    this.done = "#bbb";
    this.arrowheadColor = "#333333";
    this.fontFamily = '"trebuchet ms", verdana, arial, sans-serif';
    this.fontSize = "16px";
    this.THEME_COLOR_LIMIT = 12;
    this.nodeBkg = "calculated";
    this.nodeBorder = "calculated";
    this.clusterBkg = "calculated";
    this.clusterBorder = "calculated";
    this.defaultLinkColor = "calculated";
    this.titleColor = "calculated";
    this.edgeLabelBackground = "white";
    this.actorBorder = "calculated";
    this.actorBkg = "calculated";
    this.actorTextColor = "calculated";
    this.actorLineColor = this.actorBorder;
    this.signalColor = "calculated";
    this.signalTextColor = "calculated";
    this.labelBoxBkgColor = "calculated";
    this.labelBoxBorderColor = "calculated";
    this.labelTextColor = "calculated";
    this.loopTextColor = "calculated";
    this.noteBorderColor = "calculated";
    this.noteBkgColor = "calculated";
    this.noteTextColor = "calculated";
    this.activationBorderColor = "#666";
    this.activationBkgColor = "#f4f4f4";
    this.sequenceNumberColor = "white";
    this.sectionBkgColor = "calculated";
    this.altSectionBkgColor = "white";
    this.sectionBkgColor2 = "calculated";
    this.excludeBkgColor = "#eeeeee";
    this.taskBorderColor = "calculated";
    this.taskBkgColor = "calculated";
    this.taskTextLightColor = "white";
    this.taskTextColor = "calculated";
    this.taskTextDarkColor = "calculated";
    this.taskTextOutsideColor = "calculated";
    this.taskTextClickableColor = "#003163";
    this.activeTaskBorderColor = "calculated";
    this.activeTaskBkgColor = "calculated";
    this.gridColor = "calculated";
    this.doneTaskBkgColor = "calculated";
    this.doneTaskBorderColor = "calculated";
    this.critBkgColor = "calculated";
    this.critBorderColor = "calculated";
    this.todayLineColor = "calculated";
    this.vertLineColor = "calculated";
    this.personBorder = this.primaryBorderColor;
    this.personBkg = this.mainBkg;
    this.archEdgeColor = "calculated";
    this.archEdgeArrowColor = "calculated";
    this.archEdgeWidth = "3";
    this.archGroupBorderColor = this.primaryBorderColor;
    this.archGroupBorderWidth = "2px";
    this.rowOdd = this.rowOdd || (0,lighten/* default */.Z)(this.mainBkg, 75) || "#ffffff";
    this.rowEven = this.rowEven || "#f4f4f4";
    this.labelColor = "black";
    this.errorBkgColor = "#552222";
    this.errorTextColor = "#552222";
  }
  updateColors() {
    this.secondBkg = (0,lighten/* default */.Z)(this.contrast, 55);
    this.border2 = this.contrast;
    this.actorBorder = (0,lighten/* default */.Z)(this.border1, 23);
    this.actorBkg = this.mainBkg;
    this.actorTextColor = this.text;
    this.actorLineColor = this.actorBorder;
    this.signalColor = this.text;
    this.signalTextColor = this.text;
    this.labelBoxBkgColor = this.actorBkg;
    this.labelBoxBorderColor = this.actorBorder;
    this.labelTextColor = this.text;
    this.loopTextColor = this.text;
    this.noteBorderColor = "#999";
    this.noteBkgColor = "#666";
    this.noteTextColor = "#fff";
    this.cScale0 = this.cScale0 || "#555";
    this.cScale1 = this.cScale1 || "#F4F4F4";
    this.cScale2 = this.cScale2 || "#555";
    this.cScale3 = this.cScale3 || "#BBB";
    this.cScale4 = this.cScale4 || "#777";
    this.cScale5 = this.cScale5 || "#999";
    this.cScale6 = this.cScale6 || "#DDD";
    this.cScale7 = this.cScale7 || "#FFF";
    this.cScale8 = this.cScale8 || "#DDD";
    this.cScale9 = this.cScale9 || "#BBB";
    this.cScale10 = this.cScale10 || "#999";
    this.cScale11 = this.cScale11 || "#777";
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["cScaleInv" + i] = this["cScaleInv" + i] || methods_invert(this["cScale" + i]);
    }
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      if (this.darkMode) {
        this["cScalePeer" + i] = this["cScalePeer" + i] || (0,lighten/* default */.Z)(this["cScale" + i], 10);
      } else {
        this["cScalePeer" + i] = this["cScalePeer" + i] || (0,darken/* default */.Z)(this["cScale" + i], 10);
      }
    }
    this.scaleLabelColor = this.scaleLabelColor || (this.darkMode ? "black" : this.labelTextColor);
    this.cScaleLabel0 = this.cScaleLabel0 || this.cScale1;
    this.cScaleLabel2 = this.cScaleLabel2 || this.cScale1;
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["cScaleLabel" + i] = this["cScaleLabel" + i] || this.scaleLabelColor;
    }
    for (let i = 0; i < 5; i++) {
      this["surface" + i] = this["surface" + i] || methods_adjust(this.mainBkg, { l: -(5 + i * 5) });
      this["surfacePeer" + i] = this["surfacePeer" + i] || methods_adjust(this.mainBkg, { l: -(8 + i * 5) });
    }
    this.nodeBkg = this.mainBkg;
    this.nodeBorder = this.border1;
    this.clusterBkg = this.secondBkg;
    this.clusterBorder = this.border2;
    this.defaultLinkColor = this.lineColor;
    this.titleColor = this.text;
    this.sectionBkgColor = (0,lighten/* default */.Z)(this.contrast, 30);
    this.sectionBkgColor2 = (0,lighten/* default */.Z)(this.contrast, 30);
    this.taskBorderColor = (0,darken/* default */.Z)(this.contrast, 10);
    this.taskBkgColor = this.contrast;
    this.taskTextColor = this.taskTextLightColor;
    this.taskTextDarkColor = this.text;
    this.taskTextOutsideColor = this.taskTextDarkColor;
    this.activeTaskBorderColor = this.taskBorderColor;
    this.activeTaskBkgColor = this.mainBkg;
    this.gridColor = (0,lighten/* default */.Z)(this.border1, 30);
    this.doneTaskBkgColor = this.done;
    this.doneTaskBorderColor = this.lineColor;
    this.critBkgColor = this.critical;
    this.critBorderColor = (0,darken/* default */.Z)(this.critBkgColor, 10);
    this.todayLineColor = this.critBkgColor;
    this.vertLineColor = this.critBkgColor;
    this.archEdgeColor = this.lineColor;
    this.archEdgeArrowColor = this.lineColor;
    this.transitionColor = this.transitionColor || "#000";
    this.transitionLabelColor = this.transitionLabelColor || this.textColor;
    this.stateLabelColor = this.stateLabelColor || this.stateBkg || this.primaryTextColor;
    this.stateBkg = this.stateBkg || this.mainBkg;
    this.labelBackgroundColor = this.labelBackgroundColor || this.stateBkg;
    this.compositeBackground = this.compositeBackground || this.background || this.tertiaryColor;
    this.altBackground = this.altBackground || "#f4f4f4";
    this.compositeTitleBackground = this.compositeTitleBackground || this.mainBkg;
    this.stateBorder = this.stateBorder || "#000";
    this.innerEndBackground = this.primaryBorderColor;
    this.specialStateColor = "#222";
    this.errorBkgColor = this.errorBkgColor || this.tertiaryColor;
    this.errorTextColor = this.errorTextColor || this.tertiaryTextColor;
    this.classText = this.primaryTextColor;
    this.fillType0 = this.primaryColor;
    this.fillType1 = this.secondaryColor;
    this.fillType2 = methods_adjust(this.primaryColor, { h: 64 });
    this.fillType3 = methods_adjust(this.secondaryColor, { h: 64 });
    this.fillType4 = methods_adjust(this.primaryColor, { h: -64 });
    this.fillType5 = methods_adjust(this.secondaryColor, { h: -64 });
    this.fillType6 = methods_adjust(this.primaryColor, { h: 128 });
    this.fillType7 = methods_adjust(this.secondaryColor, { h: 128 });
    for (let i = 0; i < this.THEME_COLOR_LIMIT; i++) {
      this["pie" + i] = this["cScale" + i];
    }
    this.pie12 = this.pie0;
    this.pieTitleTextSize = this.pieTitleTextSize || "25px";
    this.pieTitleTextColor = this.pieTitleTextColor || this.taskTextDarkColor;
    this.pieSectionTextSize = this.pieSectionTextSize || "17px";
    this.pieSectionTextColor = this.pieSectionTextColor || this.textColor;
    this.pieLegendTextSize = this.pieLegendTextSize || "17px";
    this.pieLegendTextColor = this.pieLegendTextColor || this.taskTextDarkColor;
    this.pieStrokeColor = this.pieStrokeColor || "black";
    this.pieStrokeWidth = this.pieStrokeWidth || "2px";
    this.pieOuterStrokeWidth = this.pieOuterStrokeWidth || "2px";
    this.pieOuterStrokeColor = this.pieOuterStrokeColor || "black";
    this.pieOpacity = this.pieOpacity || "0.7";
    this.quadrant1Fill = this.quadrant1Fill || this.primaryColor;
    this.quadrant2Fill = this.quadrant2Fill || methods_adjust(this.primaryColor, { r: 5, g: 5, b: 5 });
    this.quadrant3Fill = this.quadrant3Fill || methods_adjust(this.primaryColor, { r: 10, g: 10, b: 10 });
    this.quadrant4Fill = this.quadrant4Fill || methods_adjust(this.primaryColor, { r: 15, g: 15, b: 15 });
    this.quadrant1TextFill = this.quadrant1TextFill || this.primaryTextColor;
    this.quadrant2TextFill = this.quadrant2TextFill || methods_adjust(this.primaryTextColor, { r: -5, g: -5, b: -5 });
    this.quadrant3TextFill = this.quadrant3TextFill || methods_adjust(this.primaryTextColor, { r: -10, g: -10, b: -10 });
    this.quadrant4TextFill = this.quadrant4TextFill || methods_adjust(this.primaryTextColor, { r: -15, g: -15, b: -15 });
    this.quadrantPointFill = this.quadrantPointFill || (0,is_dark/* default */.Z)(this.quadrant1Fill) ? (0,lighten/* default */.Z)(this.quadrant1Fill) : (0,darken/* default */.Z)(this.quadrant1Fill);
    this.quadrantPointTextFill = this.quadrantPointTextFill || this.primaryTextColor;
    this.quadrantXAxisTextFill = this.quadrantXAxisTextFill || this.primaryTextColor;
    this.quadrantYAxisTextFill = this.quadrantYAxisTextFill || this.primaryTextColor;
    this.quadrantInternalBorderStrokeFill = this.quadrantInternalBorderStrokeFill || this.primaryBorderColor;
    this.quadrantExternalBorderStrokeFill = this.quadrantExternalBorderStrokeFill || this.primaryBorderColor;
    this.quadrantTitleFill = this.quadrantTitleFill || this.primaryTextColor;
    this.xyChart = {
      backgroundColor: this.xyChart?.backgroundColor || this.background,
      titleColor: this.xyChart?.titleColor || this.primaryTextColor,
      xAxisTitleColor: this.xyChart?.xAxisTitleColor || this.primaryTextColor,
      xAxisLabelColor: this.xyChart?.xAxisLabelColor || this.primaryTextColor,
      xAxisTickColor: this.xyChart?.xAxisTickColor || this.primaryTextColor,
      xAxisLineColor: this.xyChart?.xAxisLineColor || this.primaryTextColor,
      yAxisTitleColor: this.xyChart?.yAxisTitleColor || this.primaryTextColor,
      yAxisLabelColor: this.xyChart?.yAxisLabelColor || this.primaryTextColor,
      yAxisTickColor: this.xyChart?.yAxisTickColor || this.primaryTextColor,
      yAxisLineColor: this.xyChart?.yAxisLineColor || this.primaryTextColor,
      plotColorPalette: this.xyChart?.plotColorPalette || "#EEE,#6BB8E4,#8ACB88,#C7ACD6,#E8DCC2,#FFB2A8,#FFF380,#7E8D91,#FFD8B1,#FAF3E0"
    };
    this.radar = {
      axisColor: this.radar?.axisColor || this.lineColor,
      axisStrokeWidth: this.radar?.axisStrokeWidth || 2,
      axisLabelFontSize: this.radar?.axisLabelFontSize || 12,
      curveOpacity: this.radar?.curveOpacity || 0.5,
      curveStrokeWidth: this.radar?.curveStrokeWidth || 2,
      graticuleColor: this.radar?.graticuleColor || "#DEDEDE",
      graticuleStrokeWidth: this.radar?.graticuleStrokeWidth || 1,
      graticuleOpacity: this.radar?.graticuleOpacity || 0.3,
      legendBoxSize: this.radar?.legendBoxSize || 12,
      legendFontSize: this.radar?.legendFontSize || 12
    };
    this.requirementBackground = this.requirementBackground || this.primaryColor;
    this.requirementBorderColor = this.requirementBorderColor || this.primaryBorderColor;
    this.requirementBorderSize = this.requirementBorderSize || "1";
    this.requirementTextColor = this.requirementTextColor || this.primaryTextColor;
    this.relationColor = this.relationColor || this.lineColor;
    this.relationLabelBackground = this.relationLabelBackground || this.edgeLabelBackground;
    this.relationLabelColor = this.relationLabelColor || this.actorTextColor;
    this.git0 = (0,darken/* default */.Z)(this.pie1, 25) || this.primaryColor;
    this.git1 = this.pie2 || this.secondaryColor;
    this.git2 = this.pie3 || this.tertiaryColor;
    this.git3 = this.pie4 || methods_adjust(this.primaryColor, { h: -30 });
    this.git4 = this.pie5 || methods_adjust(this.primaryColor, { h: -60 });
    this.git5 = this.pie6 || methods_adjust(this.primaryColor, { h: -90 });
    this.git6 = this.pie7 || methods_adjust(this.primaryColor, { h: 60 });
    this.git7 = this.pie8 || methods_adjust(this.primaryColor, { h: 120 });
    this.gitInv0 = this.gitInv0 || methods_invert(this.git0);
    this.gitInv1 = this.gitInv1 || methods_invert(this.git1);
    this.gitInv2 = this.gitInv2 || methods_invert(this.git2);
    this.gitInv3 = this.gitInv3 || methods_invert(this.git3);
    this.gitInv4 = this.gitInv4 || methods_invert(this.git4);
    this.gitInv5 = this.gitInv5 || methods_invert(this.git5);
    this.gitInv6 = this.gitInv6 || methods_invert(this.git6);
    this.gitInv7 = this.gitInv7 || methods_invert(this.git7);
    this.branchLabelColor = this.branchLabelColor || this.labelTextColor;
    this.gitBranchLabel0 = this.branchLabelColor;
    this.gitBranchLabel1 = "white";
    this.gitBranchLabel2 = this.branchLabelColor;
    this.gitBranchLabel3 = "white";
    this.gitBranchLabel4 = this.branchLabelColor;
    this.gitBranchLabel5 = this.branchLabelColor;
    this.gitBranchLabel6 = this.branchLabelColor;
    this.gitBranchLabel7 = this.branchLabelColor;
    this.tagLabelColor = this.tagLabelColor || this.primaryTextColor;
    this.tagLabelBackground = this.tagLabelBackground || this.primaryColor;
    this.tagLabelBorder = this.tagBorder || this.primaryBorderColor;
    this.tagLabelFontSize = this.tagLabelFontSize || "10px";
    this.commitLabelColor = this.commitLabelColor || this.secondaryTextColor;
    this.commitLabelBackground = this.commitLabelBackground || this.secondaryColor;
    this.commitLabelFontSize = this.commitLabelFontSize || "10px";
    this.attributeBackgroundColorOdd = this.attributeBackgroundColorOdd || oldAttributeBackgroundColorOdd;
    this.attributeBackgroundColorEven = this.attributeBackgroundColorEven || oldAttributeBackgroundColorEven;
  }
  calculate(overrides) {
    if (typeof overrides !== "object") {
      this.updateColors();
      return;
    }
    const keys = Object.keys(overrides);
    keys.forEach((k) => {
      this[k] = overrides[k];
    });
    this.updateColors();
    keys.forEach((k) => {
      this[k] = overrides[k];
    });
  }
};
var getThemeVariables5 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((userOverrides) => {
  const theme = new Theme5();
  theme.calculate(userOverrides);
  return theme;
}, "getThemeVariables");

// src/themes/index.js
var themes_default = {
  base: {
    getThemeVariables
  },
  dark: {
    getThemeVariables: getThemeVariables2
  },
  default: {
    getThemeVariables: getThemeVariables3
  },
  forest: {
    getThemeVariables: getThemeVariables4
  },
  neutral: {
    getThemeVariables: getThemeVariables5
  }
};

// src/schemas/config.schema.yaml?only-defaults=true
var config_schema_default = {
  "flowchart": {
    "useMaxWidth": true,
    "titleTopMargin": 25,
    "subGraphTitleMargin": {
      "top": 0,
      "bottom": 0
    },
    "diagramPadding": 8,
    "htmlLabels": true,
    "nodeSpacing": 50,
    "rankSpacing": 50,
    "curve": "basis",
    "padding": 15,
    "defaultRenderer": "dagre-wrapper",
    "wrappingWidth": 200,
    "inheritDir": false
  },
  "sequence": {
    "useMaxWidth": true,
    "hideUnusedParticipants": false,
    "activationWidth": 10,
    "diagramMarginX": 50,
    "diagramMarginY": 10,
    "actorMargin": 50,
    "width": 150,
    "height": 65,
    "boxMargin": 10,
    "boxTextMargin": 5,
    "noteMargin": 10,
    "messageMargin": 35,
    "messageAlign": "center",
    "mirrorActors": true,
    "forceMenus": false,
    "bottomMarginAdj": 1,
    "rightAngles": false,
    "showSequenceNumbers": false,
    "actorFontSize": 14,
    "actorFontFamily": '"Open Sans", sans-serif',
    "actorFontWeight": 400,
    "noteFontSize": 14,
    "noteFontFamily": '"trebuchet ms", verdana, arial, sans-serif',
    "noteFontWeight": 400,
    "noteAlign": "center",
    "messageFontSize": 16,
    "messageFontFamily": '"trebuchet ms", verdana, arial, sans-serif',
    "messageFontWeight": 400,
    "wrap": false,
    "wrapPadding": 10,
    "labelBoxWidth": 50,
    "labelBoxHeight": 20
  },
  "gantt": {
    "useMaxWidth": true,
    "titleTopMargin": 25,
    "barHeight": 20,
    "barGap": 4,
    "topPadding": 50,
    "rightPadding": 75,
    "leftPadding": 75,
    "gridLineStartPadding": 35,
    "fontSize": 11,
    "sectionFontSize": 11,
    "numberSectionStyles": 4,
    "axisFormat": "%Y-%m-%d",
    "topAxis": false,
    "displayMode": "",
    "weekday": "sunday"
  },
  "journey": {
    "useMaxWidth": true,
    "diagramMarginX": 50,
    "diagramMarginY": 10,
    "leftMargin": 150,
    "maxLabelWidth": 360,
    "width": 150,
    "height": 50,
    "boxMargin": 10,
    "boxTextMargin": 5,
    "noteMargin": 10,
    "messageMargin": 35,
    "messageAlign": "center",
    "bottomMarginAdj": 1,
    "rightAngles": false,
    "taskFontSize": 14,
    "taskFontFamily": '"Open Sans", sans-serif',
    "taskMargin": 50,
    "activationWidth": 10,
    "textPlacement": "fo",
    "actorColours": [
      "#8FBC8F",
      "#7CFC00",
      "#00FFFF",
      "#20B2AA",
      "#B0E0E6",
      "#FFFFE0"
    ],
    "sectionFills": [
      "#191970",
      "#8B008B",
      "#4B0082",
      "#2F4F4F",
      "#800000",
      "#8B4513",
      "#00008B"
    ],
    "sectionColours": [
      "#fff"
    ],
    "titleColor": "",
    "titleFontFamily": '"trebuchet ms", verdana, arial, sans-serif',
    "titleFontSize": "4ex"
  },
  "class": {
    "useMaxWidth": true,
    "titleTopMargin": 25,
    "arrowMarkerAbsolute": false,
    "dividerMargin": 10,
    "padding": 5,
    "textHeight": 10,
    "defaultRenderer": "dagre-wrapper",
    "htmlLabels": false,
    "hideEmptyMembersBox": false
  },
  "state": {
    "useMaxWidth": true,
    "titleTopMargin": 25,
    "dividerMargin": 10,
    "sizeUnit": 5,
    "padding": 8,
    "textHeight": 10,
    "titleShift": -15,
    "noteMargin": 10,
    "forkWidth": 70,
    "forkHeight": 7,
    "miniPadding": 2,
    "fontSizeFactor": 5.02,
    "fontSize": 24,
    "labelHeight": 16,
    "edgeLengthFactor": "20",
    "compositTitleSize": 35,
    "radius": 5,
    "defaultRenderer": "dagre-wrapper"
  },
  "er": {
    "useMaxWidth": true,
    "titleTopMargin": 25,
    "diagramPadding": 20,
    "layoutDirection": "TB",
    "minEntityWidth": 100,
    "minEntityHeight": 75,
    "entityPadding": 15,
    "nodeSpacing": 140,
    "rankSpacing": 80,
    "stroke": "gray",
    "fill": "honeydew",
    "fontSize": 12
  },
  "pie": {
    "useMaxWidth": true,
    "textPosition": 0.75
  },
  "quadrantChart": {
    "useMaxWidth": true,
    "chartWidth": 500,
    "chartHeight": 500,
    "titleFontSize": 20,
    "titlePadding": 10,
    "quadrantPadding": 5,
    "xAxisLabelPadding": 5,
    "yAxisLabelPadding": 5,
    "xAxisLabelFontSize": 16,
    "yAxisLabelFontSize": 16,
    "quadrantLabelFontSize": 16,
    "quadrantTextTopPadding": 5,
    "pointTextPadding": 5,
    "pointLabelFontSize": 12,
    "pointRadius": 5,
    "xAxisPosition": "top",
    "yAxisPosition": "left",
    "quadrantInternalBorderStrokeWidth": 1,
    "quadrantExternalBorderStrokeWidth": 2
  },
  "xyChart": {
    "useMaxWidth": true,
    "width": 700,
    "height": 500,
    "titleFontSize": 20,
    "titlePadding": 10,
    "showDataLabel": false,
    "showTitle": true,
    "xAxis": {
      "$ref": "#/$defs/XYChartAxisConfig",
      "showLabel": true,
      "labelFontSize": 14,
      "labelPadding": 5,
      "showTitle": true,
      "titleFontSize": 16,
      "titlePadding": 5,
      "showTick": true,
      "tickLength": 5,
      "tickWidth": 2,
      "showAxisLine": true,
      "axisLineWidth": 2
    },
    "yAxis": {
      "$ref": "#/$defs/XYChartAxisConfig",
      "showLabel": true,
      "labelFontSize": 14,
      "labelPadding": 5,
      "showTitle": true,
      "titleFontSize": 16,
      "titlePadding": 5,
      "showTick": true,
      "tickLength": 5,
      "tickWidth": 2,
      "showAxisLine": true,
      "axisLineWidth": 2
    },
    "chartOrientation": "vertical",
    "plotReservedSpacePercent": 50
  },
  "requirement": {
    "useMaxWidth": true,
    "rect_fill": "#f9f9f9",
    "text_color": "#333",
    "rect_border_size": "0.5px",
    "rect_border_color": "#bbb",
    "rect_min_width": 200,
    "rect_min_height": 200,
    "fontSize": 14,
    "rect_padding": 10,
    "line_height": 20
  },
  "mindmap": {
    "useMaxWidth": true,
    "padding": 10,
    "maxNodeWidth": 200,
    "layoutAlgorithm": "cose-bilkent"
  },
  "kanban": {
    "useMaxWidth": true,
    "padding": 8,
    "sectionWidth": 200,
    "ticketBaseUrl": ""
  },
  "timeline": {
    "useMaxWidth": true,
    "diagramMarginX": 50,
    "diagramMarginY": 10,
    "leftMargin": 150,
    "width": 150,
    "height": 50,
    "boxMargin": 10,
    "boxTextMargin": 5,
    "noteMargin": 10,
    "messageMargin": 35,
    "messageAlign": "center",
    "bottomMarginAdj": 1,
    "rightAngles": false,
    "taskFontSize": 14,
    "taskFontFamily": '"Open Sans", sans-serif',
    "taskMargin": 50,
    "activationWidth": 10,
    "textPlacement": "fo",
    "actorColours": [
      "#8FBC8F",
      "#7CFC00",
      "#00FFFF",
      "#20B2AA",
      "#B0E0E6",
      "#FFFFE0"
    ],
    "sectionFills": [
      "#191970",
      "#8B008B",
      "#4B0082",
      "#2F4F4F",
      "#800000",
      "#8B4513",
      "#00008B"
    ],
    "sectionColours": [
      "#fff"
    ],
    "disableMulticolor": false
  },
  "gitGraph": {
    "useMaxWidth": true,
    "titleTopMargin": 25,
    "diagramPadding": 8,
    "nodeLabel": {
      "width": 75,
      "height": 100,
      "x": -25,
      "y": 0
    },
    "mainBranchName": "main",
    "mainBranchOrder": 0,
    "showCommitLabel": true,
    "showBranches": true,
    "rotateCommitLabel": true,
    "parallelCommits": false,
    "arrowMarkerAbsolute": false
  },
  "c4": {
    "useMaxWidth": true,
    "diagramMarginX": 50,
    "diagramMarginY": 10,
    "c4ShapeMargin": 50,
    "c4ShapePadding": 20,
    "width": 216,
    "height": 60,
    "boxMargin": 10,
    "c4ShapeInRow": 4,
    "nextLinePaddingX": 0,
    "c4BoundaryInRow": 2,
    "personFontSize": 14,
    "personFontFamily": '"Open Sans", sans-serif',
    "personFontWeight": "normal",
    "external_personFontSize": 14,
    "external_personFontFamily": '"Open Sans", sans-serif',
    "external_personFontWeight": "normal",
    "systemFontSize": 14,
    "systemFontFamily": '"Open Sans", sans-serif',
    "systemFontWeight": "normal",
    "external_systemFontSize": 14,
    "external_systemFontFamily": '"Open Sans", sans-serif',
    "external_systemFontWeight": "normal",
    "system_dbFontSize": 14,
    "system_dbFontFamily": '"Open Sans", sans-serif',
    "system_dbFontWeight": "normal",
    "external_system_dbFontSize": 14,
    "external_system_dbFontFamily": '"Open Sans", sans-serif',
    "external_system_dbFontWeight": "normal",
    "system_queueFontSize": 14,
    "system_queueFontFamily": '"Open Sans", sans-serif',
    "system_queueFontWeight": "normal",
    "external_system_queueFontSize": 14,
    "external_system_queueFontFamily": '"Open Sans", sans-serif',
    "external_system_queueFontWeight": "normal",
    "boundaryFontSize": 14,
    "boundaryFontFamily": '"Open Sans", sans-serif',
    "boundaryFontWeight": "normal",
    "messageFontSize": 12,
    "messageFontFamily": '"Open Sans", sans-serif',
    "messageFontWeight": "normal",
    "containerFontSize": 14,
    "containerFontFamily": '"Open Sans", sans-serif',
    "containerFontWeight": "normal",
    "external_containerFontSize": 14,
    "external_containerFontFamily": '"Open Sans", sans-serif',
    "external_containerFontWeight": "normal",
    "container_dbFontSize": 14,
    "container_dbFontFamily": '"Open Sans", sans-serif',
    "container_dbFontWeight": "normal",
    "external_container_dbFontSize": 14,
    "external_container_dbFontFamily": '"Open Sans", sans-serif',
    "external_container_dbFontWeight": "normal",
    "container_queueFontSize": 14,
    "container_queueFontFamily": '"Open Sans", sans-serif',
    "container_queueFontWeight": "normal",
    "external_container_queueFontSize": 14,
    "external_container_queueFontFamily": '"Open Sans", sans-serif',
    "external_container_queueFontWeight": "normal",
    "componentFontSize": 14,
    "componentFontFamily": '"Open Sans", sans-serif',
    "componentFontWeight": "normal",
    "external_componentFontSize": 14,
    "external_componentFontFamily": '"Open Sans", sans-serif',
    "external_componentFontWeight": "normal",
    "component_dbFontSize": 14,
    "component_dbFontFamily": '"Open Sans", sans-serif',
    "component_dbFontWeight": "normal",
    "external_component_dbFontSize": 14,
    "external_component_dbFontFamily": '"Open Sans", sans-serif',
    "external_component_dbFontWeight": "normal",
    "component_queueFontSize": 14,
    "component_queueFontFamily": '"Open Sans", sans-serif',
    "component_queueFontWeight": "normal",
    "external_component_queueFontSize": 14,
    "external_component_queueFontFamily": '"Open Sans", sans-serif',
    "external_component_queueFontWeight": "normal",
    "wrap": true,
    "wrapPadding": 10,
    "person_bg_color": "#08427B",
    "person_border_color": "#073B6F",
    "external_person_bg_color": "#686868",
    "external_person_border_color": "#8A8A8A",
    "system_bg_color": "#1168BD",
    "system_border_color": "#3C7FC0",
    "system_db_bg_color": "#1168BD",
    "system_db_border_color": "#3C7FC0",
    "system_queue_bg_color": "#1168BD",
    "system_queue_border_color": "#3C7FC0",
    "external_system_bg_color": "#999999",
    "external_system_border_color": "#8A8A8A",
    "external_system_db_bg_color": "#999999",
    "external_system_db_border_color": "#8A8A8A",
    "external_system_queue_bg_color": "#999999",
    "external_system_queue_border_color": "#8A8A8A",
    "container_bg_color": "#438DD5",
    "container_border_color": "#3C7FC0",
    "container_db_bg_color": "#438DD5",
    "container_db_border_color": "#3C7FC0",
    "container_queue_bg_color": "#438DD5",
    "container_queue_border_color": "#3C7FC0",
    "external_container_bg_color": "#B3B3B3",
    "external_container_border_color": "#A6A6A6",
    "external_container_db_bg_color": "#B3B3B3",
    "external_container_db_border_color": "#A6A6A6",
    "external_container_queue_bg_color": "#B3B3B3",
    "external_container_queue_border_color": "#A6A6A6",
    "component_bg_color": "#85BBF0",
    "component_border_color": "#78A8D8",
    "component_db_bg_color": "#85BBF0",
    "component_db_border_color": "#78A8D8",
    "component_queue_bg_color": "#85BBF0",
    "component_queue_border_color": "#78A8D8",
    "external_component_bg_color": "#CCCCCC",
    "external_component_border_color": "#BFBFBF",
    "external_component_db_bg_color": "#CCCCCC",
    "external_component_db_border_color": "#BFBFBF",
    "external_component_queue_bg_color": "#CCCCCC",
    "external_component_queue_border_color": "#BFBFBF"
  },
  "sankey": {
    "useMaxWidth": true,
    "width": 600,
    "height": 400,
    "linkColor": "gradient",
    "nodeAlignment": "justify",
    "showValues": true,
    "prefix": "",
    "suffix": ""
  },
  "block": {
    "useMaxWidth": true,
    "padding": 8
  },
  "packet": {
    "useMaxWidth": true,
    "rowHeight": 32,
    "bitWidth": 32,
    "bitsPerRow": 32,
    "showBits": true,
    "paddingX": 5,
    "paddingY": 5
  },
  "architecture": {
    "useMaxWidth": true,
    "padding": 40,
    "iconSize": 80,
    "fontSize": 16
  },
  "radar": {
    "useMaxWidth": true,
    "width": 600,
    "height": 600,
    "marginTop": 50,
    "marginRight": 50,
    "marginBottom": 50,
    "marginLeft": 50,
    "axisScaleFactor": 1,
    "axisLabelFactor": 1.05,
    "curveTension": 0.17
  },
  "theme": "default",
  "look": "classic",
  "handDrawnSeed": 0,
  "layout": "dagre",
  "maxTextSize": 5e4,
  "maxEdges": 500,
  "darkMode": false,
  "fontFamily": '"trebuchet ms", verdana, arial, sans-serif;',
  "logLevel": 5,
  "securityLevel": "strict",
  "startOnLoad": true,
  "arrowMarkerAbsolute": false,
  "secure": [
    "secure",
    "securityLevel",
    "startOnLoad",
    "maxTextSize",
    "suppressErrorRendering",
    "maxEdges"
  ],
  "legacyMathML": false,
  "forceLegacyMathML": false,
  "deterministicIds": false,
  "fontSize": 16,
  "markdownAutoWrap": true,
  "suppressErrorRendering": false
};

// src/defaultConfig.ts
var config = {
  ...config_schema_default,
  // Set, even though they're `undefined` so that `configKeys` finds these keys
  // TODO: Should we replace these with `null` so that they can go in the JSON Schema?
  deterministicIDSeed: void 0,
  elk: {
    // mergeEdges is needed here to be considered
    mergeEdges: false,
    nodePlacementStrategy: "BRANDES_KOEPF",
    forceNodeModelOrder: false,
    considerModelOrder: "NODES_AND_EDGES"
  },
  themeCSS: void 0,
  // add non-JSON default config values
  themeVariables: themes_default.default.getThemeVariables(),
  sequence: {
    ...config_schema_default.sequence,
    messageFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.messageFontFamily,
        fontSize: this.messageFontSize,
        fontWeight: this.messageFontWeight
      };
    }, "messageFont"),
    noteFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.noteFontFamily,
        fontSize: this.noteFontSize,
        fontWeight: this.noteFontWeight
      };
    }, "noteFont"),
    actorFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.actorFontFamily,
        fontSize: this.actorFontSize,
        fontWeight: this.actorFontWeight
      };
    }, "actorFont")
  },
  class: {
    hideEmptyMembersBox: false
  },
  gantt: {
    ...config_schema_default.gantt,
    tickInterval: void 0,
    useWidth: void 0
    // can probably be removed since `configKeys` already includes this
  },
  c4: {
    ...config_schema_default.c4,
    useWidth: void 0,
    personFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.personFontFamily,
        fontSize: this.personFontSize,
        fontWeight: this.personFontWeight
      };
    }, "personFont"),
    flowchart: {
      ...config_schema_default.flowchart,
      inheritDir: false
      // default to legacy behavior
    },
    external_personFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.external_personFontFamily,
        fontSize: this.external_personFontSize,
        fontWeight: this.external_personFontWeight
      };
    }, "external_personFont"),
    systemFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.systemFontFamily,
        fontSize: this.systemFontSize,
        fontWeight: this.systemFontWeight
      };
    }, "systemFont"),
    external_systemFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.external_systemFontFamily,
        fontSize: this.external_systemFontSize,
        fontWeight: this.external_systemFontWeight
      };
    }, "external_systemFont"),
    system_dbFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.system_dbFontFamily,
        fontSize: this.system_dbFontSize,
        fontWeight: this.system_dbFontWeight
      };
    }, "system_dbFont"),
    external_system_dbFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.external_system_dbFontFamily,
        fontSize: this.external_system_dbFontSize,
        fontWeight: this.external_system_dbFontWeight
      };
    }, "external_system_dbFont"),
    system_queueFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.system_queueFontFamily,
        fontSize: this.system_queueFontSize,
        fontWeight: this.system_queueFontWeight
      };
    }, "system_queueFont"),
    external_system_queueFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.external_system_queueFontFamily,
        fontSize: this.external_system_queueFontSize,
        fontWeight: this.external_system_queueFontWeight
      };
    }, "external_system_queueFont"),
    containerFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.containerFontFamily,
        fontSize: this.containerFontSize,
        fontWeight: this.containerFontWeight
      };
    }, "containerFont"),
    external_containerFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.external_containerFontFamily,
        fontSize: this.external_containerFontSize,
        fontWeight: this.external_containerFontWeight
      };
    }, "external_containerFont"),
    container_dbFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.container_dbFontFamily,
        fontSize: this.container_dbFontSize,
        fontWeight: this.container_dbFontWeight
      };
    }, "container_dbFont"),
    external_container_dbFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.external_container_dbFontFamily,
        fontSize: this.external_container_dbFontSize,
        fontWeight: this.external_container_dbFontWeight
      };
    }, "external_container_dbFont"),
    container_queueFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.container_queueFontFamily,
        fontSize: this.container_queueFontSize,
        fontWeight: this.container_queueFontWeight
      };
    }, "container_queueFont"),
    external_container_queueFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.external_container_queueFontFamily,
        fontSize: this.external_container_queueFontSize,
        fontWeight: this.external_container_queueFontWeight
      };
    }, "external_container_queueFont"),
    componentFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.componentFontFamily,
        fontSize: this.componentFontSize,
        fontWeight: this.componentFontWeight
      };
    }, "componentFont"),
    external_componentFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.external_componentFontFamily,
        fontSize: this.external_componentFontSize,
        fontWeight: this.external_componentFontWeight
      };
    }, "external_componentFont"),
    component_dbFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.component_dbFontFamily,
        fontSize: this.component_dbFontSize,
        fontWeight: this.component_dbFontWeight
      };
    }, "component_dbFont"),
    external_component_dbFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.external_component_dbFontFamily,
        fontSize: this.external_component_dbFontSize,
        fontWeight: this.external_component_dbFontWeight
      };
    }, "external_component_dbFont"),
    component_queueFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.component_queueFontFamily,
        fontSize: this.component_queueFontSize,
        fontWeight: this.component_queueFontWeight
      };
    }, "component_queueFont"),
    external_component_queueFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.external_component_queueFontFamily,
        fontSize: this.external_component_queueFontSize,
        fontWeight: this.external_component_queueFontWeight
      };
    }, "external_component_queueFont"),
    boundaryFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.boundaryFontFamily,
        fontSize: this.boundaryFontSize,
        fontWeight: this.boundaryFontWeight
      };
    }, "boundaryFont"),
    messageFont: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
      return {
        fontFamily: this.messageFontFamily,
        fontSize: this.messageFontSize,
        fontWeight: this.messageFontWeight
      };
    }, "messageFont")
  },
  pie: {
    ...config_schema_default.pie,
    useWidth: 984
  },
  xyChart: {
    ...config_schema_default.xyChart,
    useWidth: void 0
  },
  requirement: {
    ...config_schema_default.requirement,
    useWidth: void 0
  },
  packet: {
    ...config_schema_default.packet
  },
  radar: {
    ...config_schema_default.radar
  },
  treemap: {
    useMaxWidth: true,
    padding: 10,
    diagramPadding: 8,
    showValues: true,
    nodeWidth: 100,
    nodeHeight: 40,
    borderWidth: 1,
    valueFontSize: 12,
    labelFontSize: 14,
    valueFormat: ","
  }
};
var keyify = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((obj, prefix = "") => Object.keys(obj).reduce((res, el) => {
  if (Array.isArray(obj[el])) {
    return res;
  } else if (typeof obj[el] === "object" && obj[el] !== null) {
    return [...res, prefix + el, ...keyify(obj[el], "")];
  }
  return [...res, prefix + el];
}, []), "keyify");
var configKeys = new Set(keyify(config, ""));
var defaultConfig_default = config;

// src/utils/sanitizeDirective.ts
var sanitizeDirective = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((args) => {
  chunk_AGHRB4JF/* log */.cM.debug("sanitizeDirective called with", args);
  if (typeof args !== "object" || args == null) {
    return;
  }
  if (Array.isArray(args)) {
    args.forEach((arg) => sanitizeDirective(arg));
    return;
  }
  for (const key of Object.keys(args)) {
    chunk_AGHRB4JF/* log */.cM.debug("Checking key", key);
    if (key.startsWith("__") || key.includes("proto") || key.includes("constr") || !configKeys.has(key) || args[key] == null) {
      chunk_AGHRB4JF/* log */.cM.debug("sanitize deleting key: ", key);
      delete args[key];
      continue;
    }
    if (typeof args[key] === "object") {
      chunk_AGHRB4JF/* log */.cM.debug("sanitizing object", key);
      sanitizeDirective(args[key]);
      continue;
    }
    const cssMatchers = ["themeCSS", "fontFamily", "altFontFamily"];
    for (const cssKey of cssMatchers) {
      if (key.includes(cssKey)) {
        chunk_AGHRB4JF/* log */.cM.debug("sanitizing css option", key);
        args[key] = sanitizeCss(args[key]);
      }
    }
  }
  if (args.themeVariables) {
    for (const k of Object.keys(args.themeVariables)) {
      const val = args.themeVariables[k];
      if (val?.match && !val.match(/^[\d "#%(),.;A-Za-z]+$/)) {
        args.themeVariables[k] = "";
      }
    }
  }
  chunk_AGHRB4JF/* log */.cM.debug("After sanitization", args);
}, "sanitizeDirective");
var sanitizeCss = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((str) => {
  let startCnt = 0;
  let endCnt = 0;
  for (const element of str) {
    if (startCnt < endCnt) {
      return "{ /* ERROR: Unbalanced CSS */ }";
    }
    if (element === "{") {
      startCnt++;
    } else if (element === "}") {
      endCnt++;
    }
  }
  if (startCnt !== endCnt) {
    return "{ /* ERROR: Unbalanced CSS */ }";
  }
  return str;
}, "sanitizeCss");

// src/config.ts
var defaultConfig = Object.freeze(defaultConfig_default);
var siteConfig = assignWithDepth_default({}, defaultConfig);
var configFromInitialize;
var directives = [];
var currentConfig = assignWithDepth_default({}, defaultConfig);
var updateCurrentConfig = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((siteCfg, _directives) => {
  let cfg = assignWithDepth_default({}, siteCfg);
  let sumOfDirectives = {};
  for (const d of _directives) {
    sanitize(d);
    sumOfDirectives = assignWithDepth_default(sumOfDirectives, d);
  }
  cfg = assignWithDepth_default(cfg, sumOfDirectives);
  if (sumOfDirectives.theme && sumOfDirectives.theme in themes_default) {
    const tmpConfigFromInitialize = assignWithDepth_default({}, configFromInitialize);
    const themeVariables = assignWithDepth_default(
      tmpConfigFromInitialize.themeVariables || {},
      sumOfDirectives.themeVariables
    );
    if (cfg.theme && cfg.theme in themes_default) {
      cfg.themeVariables = themes_default[cfg.theme].getThemeVariables(themeVariables);
    }
  }
  currentConfig = cfg;
  checkConfig(currentConfig);
  return currentConfig;
}, "updateCurrentConfig");
var setSiteConfig = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((conf) => {
  siteConfig = assignWithDepth_default({}, defaultConfig);
  siteConfig = assignWithDepth_default(siteConfig, conf);
  if (conf.theme && themes_default[conf.theme]) {
    siteConfig.themeVariables = themes_default[conf.theme].getThemeVariables(conf.themeVariables);
  }
  updateCurrentConfig(siteConfig, directives);
  return siteConfig;
}, "setSiteConfig");
var saveConfigFromInitialize = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((conf) => {
  configFromInitialize = assignWithDepth_default({}, conf);
}, "saveConfigFromInitialize");
var updateSiteConfig = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((conf) => {
  siteConfig = assignWithDepth_default(siteConfig, conf);
  updateCurrentConfig(siteConfig, directives);
  return siteConfig;
}, "updateSiteConfig");
var getSiteConfig = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
  return assignWithDepth_default({}, siteConfig);
}, "getSiteConfig");
var setConfig = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((conf) => {
  checkConfig(conf);
  assignWithDepth_default(currentConfig, conf);
  return getConfig();
}, "setConfig");
var getConfig = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
  return assignWithDepth_default({}, currentConfig);
}, "getConfig");
var sanitize = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((options) => {
  if (!options) {
    return;
  }
  ["secure", ...siteConfig.secure ?? []].forEach((key) => {
    if (Object.hasOwn(options, key)) {
      chunk_AGHRB4JF/* log */.cM.debug(`Denied attempt to modify a secure key ${key}`, options[key]);
      delete options[key];
    }
  });
  Object.keys(options).forEach((key) => {
    if (key.startsWith("__")) {
      delete options[key];
    }
  });
  Object.keys(options).forEach((key) => {
    if (typeof options[key] === "string" && (options[key].includes("<") || options[key].includes(">") || options[key].includes("url(data:"))) {
      delete options[key];
    }
    if (typeof options[key] === "object") {
      sanitize(options[key]);
    }
  });
}, "sanitize");
var addDirective = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((directive) => {
  sanitizeDirective(directive);
  if (directive.fontFamily && !directive.themeVariables?.fontFamily) {
    directive.themeVariables = {
      ...directive.themeVariables,
      fontFamily: directive.fontFamily
    };
  }
  directives.push(directive);
  updateCurrentConfig(siteConfig, directives);
}, "addDirective");
var chunk_ABZYJK2D_reset = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((config2 = siteConfig) => {
  directives = [];
  updateCurrentConfig(config2, directives);
}, "reset");
var ConfigWarning = {
  LAZY_LOAD_DEPRECATED: "The configuration options lazyLoadedDiagrams and loadExternalDiagramsAtStartup are deprecated. Please use registerExternalDiagrams instead."
};
var issuedWarnings = {};
var issueWarning = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((warning) => {
  if (issuedWarnings[warning]) {
    return;
  }
  chunk_AGHRB4JF/* log */.cM.warn(ConfigWarning[warning]);
  issuedWarnings[warning] = true;
}, "issueWarning");
var checkConfig = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((config2) => {
  if (!config2) {
    return;
  }
  if (config2.lazyLoadedDiagrams || config2.loadExternalDiagramsAtStartup) {
    issueWarning("LAZY_LOAD_DEPRECATED");
  }
}, "checkConfig");
var getUserDefinedConfig = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
  let userConfig = {};
  if (configFromInitialize) {
    userConfig = assignWithDepth_default(userConfig, configFromInitialize);
  }
  for (const d of directives) {
    userConfig = assignWithDepth_default(userConfig, d);
  }
  return userConfig;
}, "getUserDefinedConfig");

// src/diagrams/common/common.ts

var lineBreakRegex = /<br\s*\/?>/gi;
var getRows = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((s) => {
  if (!s) {
    return [""];
  }
  const str = breakToPlaceholder(s).replace(/\\n/g, "#br#");
  return str.split("#br#");
}, "getRows");
var setupDompurifyHooksIfNotSetup = /* @__PURE__ */ (() => {
  let setup = false;
  return () => {
    if (!setup) {
      setupDompurifyHooks();
      setup = true;
    }
  };
})();
function setupDompurifyHooks() {
  const TEMPORARY_ATTRIBUTE = "data-temp-href-target";
  purify_es/* default */.Z.addHook("beforeSanitizeAttributes", (node) => {
    if (node.tagName === "A" && node.hasAttribute("target")) {
      node.setAttribute(TEMPORARY_ATTRIBUTE, node.getAttribute("target") ?? "");
    }
  });
  purify_es/* default */.Z.addHook("afterSanitizeAttributes", (node) => {
    if (node.tagName === "A" && node.hasAttribute(TEMPORARY_ATTRIBUTE)) {
      node.setAttribute("target", node.getAttribute(TEMPORARY_ATTRIBUTE) ?? "");
      node.removeAttribute(TEMPORARY_ATTRIBUTE);
      if (node.getAttribute("target") === "_blank") {
        node.setAttribute("rel", "noopener");
      }
    }
  });
}
(0,chunk_AGHRB4JF/* __name */.eW)(setupDompurifyHooks, "setupDompurifyHooks");
var removeScript = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  setupDompurifyHooksIfNotSetup();
  const sanitizedText = purify_es/* default */.Z.sanitize(txt);
  return sanitizedText;
}, "removeScript");
var sanitizeMore = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((text, config2) => {
  if (config2.flowchart?.htmlLabels !== false) {
    const level = config2.securityLevel;
    if (level === "antiscript" || level === "strict") {
      text = removeScript(text);
    } else if (level !== "loose") {
      text = breakToPlaceholder(text);
      text = text.replace(/</g, "&lt;").replace(/>/g, "&gt;");
      text = text.replace(/=/g, "&equals;");
      text = placeholderToBreak(text);
    }
  }
  return text;
}, "sanitizeMore");
var sanitizeText = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((text, config2) => {
  if (!text) {
    return text;
  }
  if (config2.dompurifyConfig) {
    text = purify_es/* default */.Z.sanitize(sanitizeMore(text, config2), config2.dompurifyConfig).toString();
  } else {
    text = purify_es/* default */.Z.sanitize(sanitizeMore(text, config2), {
      FORBID_TAGS: ["style"]
    }).toString();
  }
  return text;
}, "sanitizeText");
var sanitizeTextOrArray = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((a, config2) => {
  if (typeof a === "string") {
    return sanitizeText(a, config2);
  }
  return a.flat().map((x) => sanitizeText(x, config2));
}, "sanitizeTextOrArray");
var hasBreaks = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((text) => {
  return lineBreakRegex.test(text);
}, "hasBreaks");
var splitBreaks = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((text) => {
  return text.split(lineBreakRegex);
}, "splitBreaks");
var placeholderToBreak = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((s) => {
  return s.replace(/#br#/g, "<br/>");
}, "placeholderToBreak");
var breakToPlaceholder = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((s) => {
  return s.replace(lineBreakRegex, "#br#");
}, "breakToPlaceholder");
var getUrl = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((useAbsolute) => {
  let url = "";
  if (useAbsolute) {
    url = window.location.protocol + "//" + window.location.host + window.location.pathname + window.location.search;
    url = CSS.escape(url);
  }
  return url;
}, "getUrl");
var evaluate = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((val) => val === false || ["false", "null", "0"].includes(String(val).trim().toLowerCase()) ? false : true, "evaluate");
var getMax = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(...values) {
  const newValues = values.filter((value) => {
    return !isNaN(value);
  });
  return Math.max(...newValues);
}, "getMax");
var getMin = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(...values) {
  const newValues = values.filter((value) => {
    return !isNaN(value);
  });
  return Math.min(...newValues);
}, "getMin");
var parseGenericTypes = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(input) {
  const inputSets = input.split(/(,)/);
  const output = [];
  for (let i = 0; i < inputSets.length; i++) {
    let thisSet = inputSets[i];
    if (thisSet === "," && i > 0 && i + 1 < inputSets.length) {
      const previousSet = inputSets[i - 1];
      const nextSet = inputSets[i + 1];
      if (shouldCombineSets(previousSet, nextSet)) {
        thisSet = previousSet + "," + nextSet;
        i++;
        output.pop();
      }
    }
    output.push(processSet(thisSet));
  }
  return output.join("");
}, "parseGenericTypes");
var countOccurrence = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((string, substring) => {
  return Math.max(0, string.split(substring).length - 1);
}, "countOccurrence");
var shouldCombineSets = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((previousSet, nextSet) => {
  const prevCount = countOccurrence(previousSet, "~");
  const nextCount = countOccurrence(nextSet, "~");
  return prevCount === 1 && nextCount === 1;
}, "shouldCombineSets");
var processSet = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((input) => {
  const tildeCount = countOccurrence(input, "~");
  let hasStartingTilde = false;
  if (tildeCount <= 1) {
    return input;
  }
  if (tildeCount % 2 !== 0 && input.startsWith("~")) {
    input = input.substring(1);
    hasStartingTilde = true;
  }
  const chars = [...input];
  let first = chars.indexOf("~");
  let last = chars.lastIndexOf("~");
  while (first !== -1 && last !== -1 && first !== last) {
    chars[first] = "<";
    chars[last] = ">";
    first = chars.indexOf("~");
    last = chars.lastIndexOf("~");
  }
  if (hasStartingTilde) {
    chars.unshift("~");
  }
  return chars.join("");
}, "processSet");
var isMathMLSupported = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => window.MathMLElement !== void 0, "isMathMLSupported");
var katexRegex = /\$\$(.*)\$\$/g;
var hasKatex = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((text) => (text.match(katexRegex)?.length ?? 0) > 0, "hasKatex");
var calculateMathMLDimensions = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async (text, config2) => {
  const divElem = document.createElement("div");
  divElem.innerHTML = await renderKatexSanitized(text, config2);
  divElem.id = "katex-temp";
  divElem.style.visibility = "hidden";
  divElem.style.position = "absolute";
  divElem.style.top = "0";
  const body = document.querySelector("body");
  body?.insertAdjacentElement("beforeend", divElem);
  const dim = { width: divElem.clientWidth, height: divElem.clientHeight };
  divElem.remove();
  return dim;
}, "calculateMathMLDimensions");
var renderKatexUnsanitized = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async (text, config2) => {
  if (!hasKatex(text)) {
    return text;
  }
  if (!(isMathMLSupported() || config2.legacyMathML || config2.forceLegacyMathML)) {
    return text.replace(katexRegex, "MathML is unsupported in this environment.");
  }
  if (true) {
    const { default: katex } = await __webpack_require__.e(/* import() */ 100).then(__webpack_require__.bind(__webpack_require__, 60100));
    const outputMode = config2.forceLegacyMathML || !isMathMLSupported() && config2.legacyMathML ? "htmlAndMathml" : "mathml";
    return text.split(lineBreakRegex).map(
      (line) => hasKatex(line) ? `<div style="display: flex; align-items: center; justify-content: center; white-space: nowrap;">${line}</div>` : `<div>${line}</div>`
    ).join("").replace(
      katexRegex,
      (_, c) => katex.renderToString(c, {
        throwOnError: true,
        displayMode: true,
        output: outputMode
      }).replace(/\n/g, " ").replace(/<annotation.*<\/annotation>/g, "")
    );
  }
  return text.replace(
    katexRegex,
    "Katex is not supported in @mermaid-js/tiny. Please use the full mermaid library."
  );
}, "renderKatexUnsanitized");
var renderKatexSanitized = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async (text, config2) => {
  return sanitizeText(await renderKatexUnsanitized(text, config2), config2);
}, "renderKatexSanitized");
var common_default = {
  getRows,
  sanitizeText,
  sanitizeTextOrArray,
  hasBreaks,
  splitBreaks,
  lineBreakRegex,
  removeScript,
  getUrl,
  evaluate,
  getMax,
  getMin
};

// src/setupGraphViewbox.js
var d3Attrs = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(d3Elem, attrs) {
  for (let attr of attrs) {
    d3Elem.attr(attr[0], attr[1]);
  }
}, "d3Attrs");
var calculateSvgSizeAttrs = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(height, width, useMaxWidth) {
  let attrs = /* @__PURE__ */ new Map();
  if (useMaxWidth) {
    attrs.set("width", "100%");
    attrs.set("style", `max-width: ${width}px;`);
  } else {
    attrs.set("height", height);
    attrs.set("width", width);
  }
  return attrs;
}, "calculateSvgSizeAttrs");
var configureSvgSize = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(svgElem, height, width, useMaxWidth) {
  const attrs = calculateSvgSizeAttrs(height, width, useMaxWidth);
  d3Attrs(svgElem, attrs);
}, "configureSvgSize");
var setupGraphViewbox = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(graph, svgElem, padding, useMaxWidth) {
  const svgBounds = svgElem.node().getBBox();
  const sWidth = svgBounds.width;
  const sHeight = svgBounds.height;
  chunk_AGHRB4JF/* log */.cM.info(`SVG bounds: ${sWidth}x${sHeight}`, svgBounds);
  let width = 0;
  let height = 0;
  chunk_AGHRB4JF/* log */.cM.info(`Graph bounds: ${width}x${height}`, graph);
  width = sWidth + padding * 2;
  height = sHeight + padding * 2;
  chunk_AGHRB4JF/* log */.cM.info(`Calculated bounds: ${width}x${height}`);
  configureSvgSize(svgElem, height, width, useMaxWidth);
  const vBox = `${svgBounds.x - padding} ${svgBounds.y - padding} ${svgBounds.width + 2 * padding} ${svgBounds.height + 2 * padding}`;
  svgElem.attr("viewBox", vBox);
}, "setupGraphViewbox");

// src/styles.ts
var themes = {};
var getStyles = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((type, userStyles, options) => {
  let diagramStyles = "";
  if (type in themes && themes[type]) {
    diagramStyles = themes[type](options);
  } else {
    chunk_AGHRB4JF/* log */.cM.warn(`No theme found for ${type}`);
  }
  return ` & {
    font-family: ${options.fontFamily};
    font-size: ${options.fontSize};
    fill: ${options.textColor}
  }
  @keyframes edge-animation-frame {
    from {
      stroke-dashoffset: 0;
    }
  }
  @keyframes dash {
    to {
      stroke-dashoffset: 0;
    }
  }
  & .edge-animation-slow {
    stroke-dasharray: 9,5 !important;
    stroke-dashoffset: 900;
    animation: dash 50s linear infinite;
    stroke-linecap: round;
  }
  & .edge-animation-fast {
    stroke-dasharray: 9,5 !important;
    stroke-dashoffset: 900;
    animation: dash 20s linear infinite;
    stroke-linecap: round;
  }
  /* Classes common for multiple diagrams */

  & .error-icon {
    fill: ${options.errorBkgColor};
  }
  & .error-text {
    fill: ${options.errorTextColor};
    stroke: ${options.errorTextColor};
  }

  & .edge-thickness-normal {
    stroke-width: 1px;
  }
  & .edge-thickness-thick {
    stroke-width: 3.5px
  }
  & .edge-pattern-solid {
    stroke-dasharray: 0;
  }
  & .edge-thickness-invisible {
    stroke-width: 0;
    fill: none;
  }
  & .edge-pattern-dashed{
    stroke-dasharray: 3;
  }
  .edge-pattern-dotted {
    stroke-dasharray: 2;
  }

  & .marker {
    fill: ${options.lineColor};
    stroke: ${options.lineColor};
  }
  & .marker.cross {
    stroke: ${options.lineColor};
  }

  & svg {
    font-family: ${options.fontFamily};
    font-size: ${options.fontSize};
  }
   & p {
    margin: 0
   }

  ${diagramStyles}

  ${userStyles}
`;
}, "getStyles");
var addStylesForDiagram = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((type, diagramTheme) => {
  if (diagramTheme !== void 0) {
    themes[type] = diagramTheme;
  }
}, "addStylesForDiagram");
var styles_default = getStyles;

// src/diagrams/common/commonDb.ts
var commonDb_exports = {};
(0,chunk_AGHRB4JF/* __export */.r2)(commonDb_exports, {
  clear: () => clear,
  getAccDescription: () => getAccDescription,
  getAccTitle: () => getAccTitle,
  getDiagramTitle: () => getDiagramTitle,
  setAccDescription: () => setAccDescription,
  setAccTitle: () => setAccTitle,
  setDiagramTitle: () => setDiagramTitle
});
var accTitle = "";
var diagramTitle = "";
var accDescription = "";
var sanitizeText2 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => sanitizeText(txt, getConfig()), "sanitizeText");
var clear = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
  accTitle = "";
  accDescription = "";
  diagramTitle = "";
}, "clear");
var setAccTitle = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  accTitle = sanitizeText2(txt).replace(/^\s+/g, "");
}, "setAccTitle");
var getAccTitle = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => accTitle, "getAccTitle");
var setAccDescription = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  accDescription = sanitizeText2(txt).replace(/\n\s+/g, "\n");
}, "setAccDescription");
var getAccDescription = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => accDescription, "getAccDescription");
var setDiagramTitle = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  diagramTitle = sanitizeText2(txt);
}, "setDiagramTitle");
var getDiagramTitle = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => diagramTitle, "getDiagramTitle");

// src/diagram-api/diagramAPI.ts
var log2 = chunk_AGHRB4JF/* log */.cM;
var setLogLevel2 = chunk_AGHRB4JF/* setLogLevel */.Ub;
var getConfig2 = getConfig;
var setConfig2 = setConfig;
var defaultConfig2 = defaultConfig;
var sanitizeText3 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((text) => sanitizeText(text, getConfig2()), "sanitizeText");
var setupGraphViewbox2 = setupGraphViewbox;
var getCommonDb = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
  return commonDb_exports;
}, "getCommonDb");
var diagrams = {};
var registerDiagram = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((id, diagram, detector) => {
  if (diagrams[id]) {
    log2.warn(`Diagram with id ${id} already registered. Overwriting.`);
  }
  diagrams[id] = diagram;
  if (detector) {
    addDetector(id, detector);
  }
  addStylesForDiagram(id, diagram.styles);
  diagram.injectUtils?.(
    log2,
    setLogLevel2,
    getConfig2,
    sanitizeText3,
    setupGraphViewbox2,
    getCommonDb(),
    () => {
    }
  );
}, "registerDiagram");
var getDiagram = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((name) => {
  if (name in diagrams) {
    return diagrams[name];
  }
  throw new DiagramNotFoundError(name);
}, "getDiagram");
var DiagramNotFoundError = class extends Error {
  static {
    (0,chunk_AGHRB4JF/* __name */.eW)(this, "DiagramNotFoundError");
  }
  constructor(name) {
    super(`Diagram ${name} not found.`);
  }
};




/***/ }),

/***/ 74999:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Ub: () => (/* binding */ setLogLevel),
/* harmony export */   cM: () => (/* binding */ log),
/* harmony export */   eW: () => (/* binding */ __name),
/* harmony export */   r2: () => (/* binding */ __export)
/* harmony export */ });
/* harmony import */ var dayjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(59395);
/* harmony import */ var dayjs__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(dayjs__WEBPACK_IMPORTED_MODULE_0__);
var __defProp = Object.defineProperty;
var __name = (target, value) => __defProp(target, "name", { value, configurable: true });
var __export = (target, all) => {
  for (var name in all)
    __defProp(target, name, { get: all[name], enumerable: true });
};

// src/logger.ts

var LEVELS = {
  trace: 0,
  debug: 1,
  info: 2,
  warn: 3,
  error: 4,
  fatal: 5
};
var log = {
  trace: /* @__PURE__ */ __name((..._args) => {
  }, "trace"),
  debug: /* @__PURE__ */ __name((..._args) => {
  }, "debug"),
  info: /* @__PURE__ */ __name((..._args) => {
  }, "info"),
  warn: /* @__PURE__ */ __name((..._args) => {
  }, "warn"),
  error: /* @__PURE__ */ __name((..._args) => {
  }, "error"),
  fatal: /* @__PURE__ */ __name((..._args) => {
  }, "fatal")
};
var setLogLevel = /* @__PURE__ */ __name(function(level = "fatal") {
  let numericLevel = LEVELS.fatal;
  if (typeof level === "string") {
    if (level.toLowerCase() in LEVELS) {
      numericLevel = LEVELS[level];
    }
  } else if (typeof level === "number") {
    numericLevel = level;
  }
  log.trace = () => {
  };
  log.debug = () => {
  };
  log.info = () => {
  };
  log.warn = () => {
  };
  log.error = () => {
  };
  log.fatal = () => {
  };
  if (numericLevel <= LEVELS.fatal) {
    log.fatal = console.error ? console.error.bind(console, format("FATAL"), "color: orange") : console.log.bind(console, "\x1B[35m", format("FATAL"));
  }
  if (numericLevel <= LEVELS.error) {
    log.error = console.error ? console.error.bind(console, format("ERROR"), "color: orange") : console.log.bind(console, "\x1B[31m", format("ERROR"));
  }
  if (numericLevel <= LEVELS.warn) {
    log.warn = console.warn ? console.warn.bind(console, format("WARN"), "color: orange") : console.log.bind(console, `\x1B[33m`, format("WARN"));
  }
  if (numericLevel <= LEVELS.info) {
    log.info = console.info ? console.info.bind(console, format("INFO"), "color: lightblue") : console.log.bind(console, "\x1B[34m", format("INFO"));
  }
  if (numericLevel <= LEVELS.debug) {
    log.debug = console.debug ? console.debug.bind(console, format("DEBUG"), "color: lightgreen") : console.log.bind(console, "\x1B[32m", format("DEBUG"));
  }
  if (numericLevel <= LEVELS.trace) {
    log.trace = console.debug ? console.debug.bind(console, format("TRACE"), "color: lightgreen") : console.log.bind(console, "\x1B[32m", format("TRACE"));
  }
}, "setLogLevel");
var format = /* @__PURE__ */ __name((level) => {
  const time = dayjs__WEBPACK_IMPORTED_MODULE_0___default()().format("ss.SSS");
  return `%c${time} : ${level} : `;
}, "format");




/***/ }),

/***/ 37756:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Fh: () => (/* binding */ isLabelStyle),
/* harmony export */   TA: () => (/* binding */ solidStateFill),
/* harmony export */   UG: () => (/* binding */ styles2String),
/* harmony export */   _q: () => (/* binding */ userNodeOverrides),
/* harmony export */   dT: () => (/* binding */ compileStyles)
/* harmony export */ });
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(74999);



// src/rendering-util/rendering-elements/shapes/handDrawnShapeStyles.ts
var solidStateFill = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((color) => {
  const { handDrawnSeed } = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getConfig2 */ .nV)();
  return {
    fill: color,
    hachureAngle: 120,
    // angle of hachure,
    hachureGap: 4,
    fillWeight: 2,
    roughness: 0.7,
    stroke: color,
    seed: handDrawnSeed
  };
}, "solidStateFill");
var compileStyles = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((node) => {
  const stylesMap = styles2Map([
    ...node.cssCompiledStyles || [],
    ...node.cssStyles || [],
    ...node.labelStyle || []
  ]);
  return { stylesMap, stylesArray: [...stylesMap] };
}, "compileStyles");
var styles2Map = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((styles) => {
  const styleMap = /* @__PURE__ */ new Map();
  styles.forEach((style) => {
    const [key, value] = style.split(":");
    styleMap.set(key.trim(), value?.trim());
  });
  return styleMap;
}, "styles2Map");
var isLabelStyle = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((key) => {
  return key === "color" || key === "font-size" || key === "font-family" || key === "font-weight" || key === "font-style" || key === "text-decoration" || key === "text-align" || key === "text-transform" || key === "line-height" || key === "letter-spacing" || key === "word-spacing" || key === "text-shadow" || key === "text-overflow" || key === "white-space" || key === "word-wrap" || key === "word-break" || key === "overflow-wrap" || key === "hyphens";
}, "isLabelStyle");
var styles2String = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((node) => {
  const { stylesArray } = compileStyles(node);
  const labelStyles = [];
  const nodeStyles = [];
  const borderStyles = [];
  const backgroundStyles = [];
  stylesArray.forEach((style) => {
    const key = style[0];
    if (isLabelStyle(key)) {
      labelStyles.push(style.join(":") + " !important");
    } else {
      nodeStyles.push(style.join(":") + " !important");
      if (key.includes("stroke")) {
        borderStyles.push(style.join(":") + " !important");
      }
      if (key === "fill") {
        backgroundStyles.push(style.join(":") + " !important");
      }
    }
  });
  return {
    labelStyles: labelStyles.join(";"),
    nodeStyles: nodeStyles.join(";"),
    stylesArray,
    borderStyles,
    backgroundStyles
  };
}, "styles2String");
var userNodeOverrides = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((node, options) => {
  const { themeVariables, handDrawnSeed } = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getConfig2 */ .nV)();
  const { nodeBorder, mainBkg } = themeVariables;
  const { stylesMap } = compileStyles(node);
  const result = Object.assign(
    {
      roughness: 0.7,
      fill: stylesMap.get("fill") || mainBkg,
      fillStyle: "hachure",
      // solid fill
      fillWeight: 4,
      hachureGap: 5.2,
      stroke: stylesMap.get("stroke") || nodeBorder,
      seed: handDrawnSeed,
      strokeWidth: stylesMap.get("stroke-width")?.replace("px", "") || 1.3,
      fillLineDash: [0, 0],
      strokeLineDash: getStrokeDashArray(stylesMap.get("stroke-dasharray"))
    },
    options
  );
  return result;
}, "userNodeOverrides");
var getStrokeDashArray = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((strokeDasharrayStyle) => {
  if (!strokeDasharrayStyle) {
    return [0, 0];
  }
  const dashArray = strokeDasharrayStyle.trim().split(/\s+/).map(Number);
  if (dashArray.length === 1) {
    const val = isNaN(dashArray[0]) ? 0 : dashArray[0];
    return [val, val];
  }
  const first = isNaN(dashArray[0]) ? 0 : dashArray[0];
  const second = isNaN(dashArray[1]) ? 0 : dashArray[1];
  return [first, second];
}, "getStrokeDashArray");




/***/ }),

/***/ 22462:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   L: () => (/* binding */ getSubGraphTitleMargins)
/* harmony export */ });
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(74999);


// src/utils/subGraphTitleMargins.ts
var getSubGraphTitleMargins = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(({
  flowchart
}) => {
  const subGraphTitleTopMargin = flowchart?.subGraphTitleMargin?.top ?? 0;
  const subGraphTitleBottomMargin = flowchart?.subGraphTitleMargin?.bottom ?? 0;
  const subGraphTitleTotalMargin = subGraphTitleTopMargin + subGraphTitleBottomMargin;
  return {
    subGraphTitleTopMargin,
    subGraphTitleBottomMargin,
    subGraphTitleTotalMargin
  };
}, "getSubGraphTitleMargins");




/***/ }),

/***/ 98154:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   P: () => (/* binding */ selectSvgElement)
/* harmony export */ });
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(74999);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(35321);



// src/rendering-util/selectSvgElement.ts

var selectSvgElement = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((id) => {
  const { securityLevel } = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getConfig2 */ .nV)();
  let root = (0,d3__WEBPACK_IMPORTED_MODULE_2__/* .select */ .Ys)("body");
  if (securityLevel === "sandbox") {
    const sandboxElement = (0,d3__WEBPACK_IMPORTED_MODULE_2__/* .select */ .Ys)(`#i${id}`);
    const doc = sandboxElement.node()?.contentDocument ?? document;
    root = (0,d3__WEBPACK_IMPORTED_MODULE_2__/* .select */ .Ys)(doc.body);
  }
  const svg = root.select(`#${id}`);
  return svg;
}, "selectSvgElement");




/***/ }),

/***/ 67894:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Os: () => (/* binding */ markerOffsets),
/* harmony export */   oZ: () => (/* binding */ markerOffsets2),
/* harmony export */   om: () => (/* binding */ getLineFunctionsWithOffset)
/* harmony export */ });
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(74999);


// src/utils/lineWithOffset.ts
var markerOffsets = {
  aggregation: 17.25,
  extension: 17.25,
  composition: 17.25,
  dependency: 6,
  lollipop: 13.5,
  arrow_point: 4
  //arrow_cross: 24,
};
var markerOffsets2 = {
  arrow_point: 9,
  arrow_cross: 12.5,
  arrow_circle: 12.5
};
function calculateDeltaAndAngle(point1, point2) {
  if (point1 === void 0 || point2 === void 0) {
    return { angle: 0, deltaX: 0, deltaY: 0 };
  }
  point1 = pointTransformer(point1);
  point2 = pointTransformer(point2);
  const [x1, y1] = [point1.x, point1.y];
  const [x2, y2] = [point2.x, point2.y];
  const deltaX = x2 - x1;
  const deltaY = y2 - y1;
  return { angle: Math.atan(deltaY / deltaX), deltaX, deltaY };
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(calculateDeltaAndAngle, "calculateDeltaAndAngle");
var pointTransformer = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((data) => {
  if (Array.isArray(data)) {
    return { x: data[0], y: data[1] };
  }
  return data;
}, "pointTransformer");
var getLineFunctionsWithOffset = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)((edge) => {
  return {
    x: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function(d, i, data) {
      let offset = 0;
      const DIRECTION = pointTransformer(data[0]).x < pointTransformer(data[data.length - 1]).x ? "left" : "right";
      if (i === 0 && Object.hasOwn(markerOffsets, edge.arrowTypeStart)) {
        const { angle, deltaX } = calculateDeltaAndAngle(data[0], data[1]);
        offset = markerOffsets[edge.arrowTypeStart] * Math.cos(angle) * (deltaX >= 0 ? 1 : -1);
      } else if (i === data.length - 1 && Object.hasOwn(markerOffsets, edge.arrowTypeEnd)) {
        const { angle, deltaX } = calculateDeltaAndAngle(
          data[data.length - 1],
          data[data.length - 2]
        );
        offset = markerOffsets[edge.arrowTypeEnd] * Math.cos(angle) * (deltaX >= 0 ? 1 : -1);
      }
      const differenceToEnd = Math.abs(
        pointTransformer(d).x - pointTransformer(data[data.length - 1]).x
      );
      const differenceInYEnd = Math.abs(
        pointTransformer(d).y - pointTransformer(data[data.length - 1]).y
      );
      const differenceToStart = Math.abs(pointTransformer(d).x - pointTransformer(data[0]).x);
      const differenceInYStart = Math.abs(pointTransformer(d).y - pointTransformer(data[0]).y);
      const startMarkerHeight = markerOffsets[edge.arrowTypeStart];
      const endMarkerHeight = markerOffsets[edge.arrowTypeEnd];
      const extraRoom = 1;
      if (differenceToEnd < endMarkerHeight && differenceToEnd > 0 && differenceInYEnd < endMarkerHeight) {
        let adjustment = endMarkerHeight + extraRoom - differenceToEnd;
        adjustment *= DIRECTION === "right" ? -1 : 1;
        offset -= adjustment;
      }
      if (differenceToStart < startMarkerHeight && differenceToStart > 0 && differenceInYStart < startMarkerHeight) {
        let adjustment = startMarkerHeight + extraRoom - differenceToStart;
        adjustment *= DIRECTION === "right" ? -1 : 1;
        offset += adjustment;
      }
      return pointTransformer(d).x + offset;
    }, "x"),
    y: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function(d, i, data) {
      let offset = 0;
      const DIRECTION = pointTransformer(data[0]).y < pointTransformer(data[data.length - 1]).y ? "down" : "up";
      if (i === 0 && Object.hasOwn(markerOffsets, edge.arrowTypeStart)) {
        const { angle, deltaY } = calculateDeltaAndAngle(data[0], data[1]);
        offset = markerOffsets[edge.arrowTypeStart] * Math.abs(Math.sin(angle)) * (deltaY >= 0 ? 1 : -1);
      } else if (i === data.length - 1 && Object.hasOwn(markerOffsets, edge.arrowTypeEnd)) {
        const { angle, deltaY } = calculateDeltaAndAngle(
          data[data.length - 1],
          data[data.length - 2]
        );
        offset = markerOffsets[edge.arrowTypeEnd] * Math.abs(Math.sin(angle)) * (deltaY >= 0 ? 1 : -1);
      }
      const differenceToEnd = Math.abs(
        pointTransformer(d).y - pointTransformer(data[data.length - 1]).y
      );
      const differenceInXEnd = Math.abs(
        pointTransformer(d).x - pointTransformer(data[data.length - 1]).x
      );
      const differenceToStart = Math.abs(pointTransformer(d).y - pointTransformer(data[0]).y);
      const differenceInXStart = Math.abs(pointTransformer(d).x - pointTransformer(data[0]).x);
      const startMarkerHeight = markerOffsets[edge.arrowTypeStart];
      const endMarkerHeight = markerOffsets[edge.arrowTypeEnd];
      const extraRoom = 1;
      if (differenceToEnd < endMarkerHeight && differenceToEnd > 0 && differenceInXEnd < endMarkerHeight) {
        let adjustment = endMarkerHeight + extraRoom - differenceToEnd;
        adjustment *= DIRECTION === "up" ? -1 : 1;
        offset -= adjustment;
      }
      if (differenceToStart < startMarkerHeight && differenceToStart > 0 && differenceInXStart < startMarkerHeight) {
        let adjustment = startMarkerHeight + extraRoom - differenceToStart;
        adjustment *= DIRECTION === "up" ? -1 : 1;
        offset += adjustment;
      }
      return pointTransformer(d).y + offset;
    }, "y")
  };
}, "getLineFunctionsWithOffset");
if (void 0) {
  const { it, expect, describe } = void 0;
  describe("calculateDeltaAndAngle", () => {
    it("should calculate the angle and deltas between two points", () => {
      expect(calculateDeltaAndAngle([0, 0], [0, 1])).toStrictEqual({
        angle: 1.5707963267948966,
        deltaX: 0,
        deltaY: 1
      });
      expect(calculateDeltaAndAngle([1, 0], [0, -1])).toStrictEqual({
        angle: 0.7853981633974483,
        deltaX: -1,
        deltaY: -1
      });
      expect(calculateDeltaAndAngle({ x: 1, y: 0 }, [0, -1])).toStrictEqual({
        angle: 0.7853981633974483,
        deltaX: -1,
        deltaY: -1
      });
      expect(calculateDeltaAndAngle({ x: 1, y: 0 }, { x: 1, y: 0 })).toStrictEqual({
        angle: NaN,
        deltaX: 0,
        deltaY: 0
      });
    });
    it("should calculate the angle and deltas if one point in undefined", () => {
      expect(calculateDeltaAndAngle(void 0, [0, 1])).toStrictEqual({
        angle: 0,
        deltaX: 0,
        deltaY: 0
      });
      expect(calculateDeltaAndAngle([0, 1], void 0)).toStrictEqual({
        angle: 0,
        deltaX: 0,
        deltaY: 0
      });
    });
  });
}




/***/ }),

/***/ 2111:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  QA: () => (/* binding */ computeDimensionOfText),
  rw: () => (/* binding */ createText),
  s4: () => (/* binding */ getIconSVG),
  ef: () => (/* binding */ registerIconPacks),
  EY: () => (/* binding */ replaceIconSubstring),
  cN: () => (/* binding */ unknownIcon)
});

// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-S3R3BYOJ.mjs
var chunk_S3R3BYOJ = __webpack_require__(17175);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-ABZYJK2D.mjs + 3 modules
var chunk_ABZYJK2D = __webpack_require__(61805);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-AGHRB4JF.mjs
var chunk_AGHRB4JF = __webpack_require__(74999);
;// CONCATENATED MODULE: ../node_modules/@iconify/utils/lib/icon/name.js
/**
* Expression to test part of icon name.
*
* Used when loading icons from Iconify API due to project naming convension.
* Ignored when using custom icon sets - convension does not apply.
*/
const matchIconName = /^[a-z0-9]+(-[a-z0-9]+)*$/;
/**
* Convert string icon name to IconifyIconName object.
*/
const stringToIcon = (value, validate, allowSimpleName, provider = "") => {
	const colonSeparated = value.split(":");
	if (value.slice(0, 1) === "@") {
		if (colonSeparated.length < 2 || colonSeparated.length > 3) return null;
		provider = colonSeparated.shift().slice(1);
	}
	if (colonSeparated.length > 3 || !colonSeparated.length) return null;
	if (colonSeparated.length > 1) {
		const name$1 = colonSeparated.pop();
		const prefix = colonSeparated.pop();
		const result = {
			provider: colonSeparated.length > 0 ? colonSeparated[0] : provider,
			prefix,
			name: name$1
		};
		return validate && !validateIconName(result) ? null : result;
	}
	const name = colonSeparated[0];
	const dashSeparated = name.split("-");
	if (dashSeparated.length > 1) {
		const result = {
			provider,
			prefix: dashSeparated.shift(),
			name: dashSeparated.join("-")
		};
		return validate && !validateIconName(result) ? null : result;
	}
	if (allowSimpleName && provider === "") {
		const result = {
			provider,
			prefix: "",
			name
		};
		return validate && !validateIconName(result, allowSimpleName) ? null : result;
	}
	return null;
};
/**
* Check if icon is valid.
*
* This function is not part of stringToIcon because validation is not needed for most code.
*/
const validateIconName = (icon, allowSimpleName) => {
	if (!icon) return false;
	return !!((allowSimpleName && icon.prefix === "" || !!icon.prefix) && !!icon.name);
};


;// CONCATENATED MODULE: ../node_modules/@iconify/utils/lib/icon/defaults.js
/**
* Default values for dimensions
*/
const defaultIconDimensions = Object.freeze({
	left: 0,
	top: 0,
	width: 16,
	height: 16
});
/**
* Default values for transformations
*/
const defaultIconTransformations = Object.freeze({
	rotate: 0,
	vFlip: false,
	hFlip: false
});
/**
* Default values for all optional IconifyIcon properties
*/
const defaultIconProps = Object.freeze({
	...defaultIconDimensions,
	...defaultIconTransformations
});
/**
* Default values for all properties used in ExtendedIconifyIcon
*/
const defaultExtendedIconProps = Object.freeze({
	...defaultIconProps,
	body: "",
	hidden: false
});


;// CONCATENATED MODULE: ../node_modules/@iconify/utils/lib/icon/transformations.js
/**
* Merge transformations
*/
function mergeIconTransformations(obj1, obj2) {
	const result = {};
	if (!obj1.hFlip !== !obj2.hFlip) result.hFlip = true;
	if (!obj1.vFlip !== !obj2.vFlip) result.vFlip = true;
	const rotate = ((obj1.rotate || 0) + (obj2.rotate || 0)) % 4;
	if (rotate) result.rotate = rotate;
	return result;
}


;// CONCATENATED MODULE: ../node_modules/@iconify/utils/lib/icon/merge.js



/**
* Merge icon and alias
*
* Can also be used to merge default values and icon
*/
function mergeIconData(parent, child) {
	const result = mergeIconTransformations(parent, child);
	for (const key in defaultExtendedIconProps) if (key in defaultIconTransformations) {
		if (key in parent && !(key in result)) result[key] = defaultIconTransformations[key];
	} else if (key in child) result[key] = child[key];
	else if (key in parent) result[key] = parent[key];
	return result;
}


;// CONCATENATED MODULE: ../node_modules/@iconify/utils/lib/icon-set/tree.js
/**
* Resolve icon set icons
*
* Returns parent icon for each icon
*/
function getIconsTree(data, names) {
	const icons = data.icons;
	const aliases = data.aliases || Object.create(null);
	const resolved = Object.create(null);
	function resolve(name) {
		if (icons[name]) return resolved[name] = [];
		if (!(name in resolved)) {
			resolved[name] = null;
			const parent = aliases[name] && aliases[name].parent;
			const value = parent && resolve(parent);
			if (value) resolved[name] = [parent].concat(value);
		}
		return resolved[name];
	}
	(names || Object.keys(icons).concat(Object.keys(aliases))).forEach(resolve);
	return resolved;
}


;// CONCATENATED MODULE: ../node_modules/@iconify/utils/lib/icon-set/get-icon.js



/**
* Get icon data, using prepared aliases tree
*/
function internalGetIconData(data, name, tree) {
	const icons = data.icons;
	const aliases = data.aliases || Object.create(null);
	let currentProps = {};
	function parse(name$1) {
		currentProps = mergeIconData(icons[name$1] || aliases[name$1], currentProps);
	}
	parse(name);
	tree.forEach(parse);
	return mergeIconData(data, currentProps);
}
/**
* Get data for icon
*/
function getIconData(data, name) {
	if (data.icons[name]) return internalGetIconData(data, name, []);
	const tree = getIconsTree(data, [name])[name];
	return tree ? internalGetIconData(data, name, tree) : null;
}


;// CONCATENATED MODULE: ../node_modules/@iconify/utils/lib/customisations/defaults.js


/**
* Default icon customisations values
*/
const defaultIconSizeCustomisations = Object.freeze({
	width: null,
	height: null
});
const defaultIconCustomisations = Object.freeze({
	...defaultIconSizeCustomisations,
	...defaultIconTransformations
});


;// CONCATENATED MODULE: ../node_modules/@iconify/utils/lib/svg/size.js
/**
* Regular expressions for calculating dimensions
*/
const unitsSplit = /(-?[0-9.]*[0-9]+[0-9.]*)/g;
const unitsTest = /^-?[0-9.]*[0-9]+[0-9.]*$/g;
function calculateSize(size, ratio, precision) {
	if (ratio === 1) return size;
	precision = precision || 100;
	if (typeof size === "number") return Math.ceil(size * ratio * precision) / precision;
	if (typeof size !== "string") return size;
	const oldParts = size.split(unitsSplit);
	if (oldParts === null || !oldParts.length) return size;
	const newParts = [];
	let code = oldParts.shift();
	let isNumber = unitsTest.test(code);
	while (true) {
		if (isNumber) {
			const num = parseFloat(code);
			if (isNaN(num)) newParts.push(code);
			else newParts.push(Math.ceil(num * ratio * precision) / precision);
		} else newParts.push(code);
		code = oldParts.shift();
		if (code === void 0) return newParts.join("");
		isNumber = !isNumber;
	}
}


;// CONCATENATED MODULE: ../node_modules/@iconify/utils/lib/svg/defs.js
function splitSVGDefs(content, tag = "defs") {
	let defs = "";
	const index = content.indexOf("<" + tag);
	while (index >= 0) {
		const start = content.indexOf(">", index);
		const end = content.indexOf("</" + tag);
		if (start === -1 || end === -1) break;
		const endEnd = content.indexOf(">", end);
		if (endEnd === -1) break;
		defs += content.slice(start + 1, end).trim();
		content = content.slice(0, index).trim() + content.slice(endEnd + 1);
	}
	return {
		defs,
		content
	};
}
/**
* Merge defs and content
*/
function mergeDefsAndContent(defs, content) {
	return defs ? "<defs>" + defs + "</defs>" + content : content;
}
/**
* Wrap SVG content, without wrapping definitions
*/
function wrapSVGContent(body, start, end) {
	const split = splitSVGDefs(body);
	return mergeDefsAndContent(split.defs, start + split.content + end);
}


;// CONCATENATED MODULE: ../node_modules/@iconify/utils/lib/svg/build.js





/**
* Check if value should be unset. Allows multiple keywords
*/
const isUnsetKeyword = (value) => value === "unset" || value === "undefined" || value === "none";
/**
* Get SVG attributes and content from icon + customisations
*
* Does not generate style to make it compatible with frameworks that use objects for style, such as React.
* Instead, it generates 'inline' value. If true, rendering engine should add verticalAlign: -0.125em to icon.
*
* Customisations should be normalised by platform specific parser.
* Result should be converted to <svg> by platform specific parser.
* Use replaceIDs to generate unique IDs for body.
*/
function iconToSVG(icon, customisations) {
	const fullIcon = {
		...defaultIconProps,
		...icon
	};
	const fullCustomisations = {
		...defaultIconCustomisations,
		...customisations
	};
	const box = {
		left: fullIcon.left,
		top: fullIcon.top,
		width: fullIcon.width,
		height: fullIcon.height
	};
	let body = fullIcon.body;
	[fullIcon, fullCustomisations].forEach((props) => {
		const transformations = [];
		const hFlip = props.hFlip;
		const vFlip = props.vFlip;
		let rotation = props.rotate;
		if (hFlip) if (vFlip) rotation += 2;
		else {
			transformations.push("translate(" + (box.width + box.left).toString() + " " + (0 - box.top).toString() + ")");
			transformations.push("scale(-1 1)");
			box.top = box.left = 0;
		}
		else if (vFlip) {
			transformations.push("translate(" + (0 - box.left).toString() + " " + (box.height + box.top).toString() + ")");
			transformations.push("scale(1 -1)");
			box.top = box.left = 0;
		}
		let tempValue;
		if (rotation < 0) rotation -= Math.floor(rotation / 4) * 4;
		rotation = rotation % 4;
		switch (rotation) {
			case 1:
				tempValue = box.height / 2 + box.top;
				transformations.unshift("rotate(90 " + tempValue.toString() + " " + tempValue.toString() + ")");
				break;
			case 2:
				transformations.unshift("rotate(180 " + (box.width / 2 + box.left).toString() + " " + (box.height / 2 + box.top).toString() + ")");
				break;
			case 3:
				tempValue = box.width / 2 + box.left;
				transformations.unshift("rotate(-90 " + tempValue.toString() + " " + tempValue.toString() + ")");
				break;
		}
		if (rotation % 2 === 1) {
			if (box.left !== box.top) {
				tempValue = box.left;
				box.left = box.top;
				box.top = tempValue;
			}
			if (box.width !== box.height) {
				tempValue = box.width;
				box.width = box.height;
				box.height = tempValue;
			}
		}
		if (transformations.length) body = wrapSVGContent(body, "<g transform=\"" + transformations.join(" ") + "\">", "</g>");
	});
	const customisationsWidth = fullCustomisations.width;
	const customisationsHeight = fullCustomisations.height;
	const boxWidth = box.width;
	const boxHeight = box.height;
	let width;
	let height;
	if (customisationsWidth === null) {
		height = customisationsHeight === null ? "1em" : customisationsHeight === "auto" ? boxHeight : customisationsHeight;
		width = calculateSize(height, boxWidth / boxHeight);
	} else {
		width = customisationsWidth === "auto" ? boxWidth : customisationsWidth;
		height = customisationsHeight === null ? calculateSize(width, boxHeight / boxWidth) : customisationsHeight === "auto" ? boxHeight : customisationsHeight;
	}
	const attributes = {};
	const setAttr = (prop, value) => {
		if (!isUnsetKeyword(value)) attributes[prop] = value.toString();
	};
	setAttr("width", width);
	setAttr("height", height);
	const viewBox = [
		box.left,
		box.top,
		boxWidth,
		boxHeight
	];
	attributes.viewBox = viewBox.join(" ");
	return {
		attributes,
		viewBox,
		body
	};
}


;// CONCATENATED MODULE: ../node_modules/@iconify/utils/lib/svg/html.js
/**
* Generate <svg>
*/
function iconToHTML(body, attributes) {
	let renderAttribsHTML = body.indexOf("xlink:") === -1 ? "" : " xmlns:xlink=\"http://www.w3.org/1999/xlink\"";
	for (const attr in attributes) renderAttribsHTML += " " + attr + "=\"" + attributes[attr] + "\"";
	return "<svg xmlns=\"http://www.w3.org/2000/svg\"" + renderAttribsHTML + ">" + body + "</svg>";
}


;// CONCATENATED MODULE: ../node_modules/@iconify/utils/lib/svg/id.js
/**
* IDs usage:
*
* id="{id}"
* xlink:href="#{id}"
* url(#{id})
*
* From SVG animations:
*
* begin="0;{id}.end"
* begin="{id}.end"
* begin="{id}.click"
*/
/**
* Regular expression for finding ids
*/
const regex = /\sid="(\S+)"/g;
/**
* New random-ish prefix for ids
*
* Do not use dash, it cannot be used in SVG 2 animations
*/
const randomPrefix = "IconifyId" + Date.now().toString(16) + (Math.random() * 16777216 | 0).toString(16);
/**
* Counter for ids, increasing with every replacement
*/
let counter = 0;
/**
* Replace IDs in SVG output with unique IDs
*/
function replaceIDs(body, prefix = randomPrefix) {
	const ids = [];
	let match;
	while (match = regex.exec(body)) ids.push(match[1]);
	if (!ids.length) return body;
	const suffix = "suffix" + (Math.random() * 16777216 | Date.now()).toString(16);
	ids.forEach((id) => {
		const newID = typeof prefix === "function" ? prefix(id) : prefix + (counter++).toString();
		const escapedID = id.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
		body = body.replace(new RegExp("([#;\"])(" + escapedID + ")([\")]|\\.[a-z])", "g"), "$1" + newID + suffix + "$3");
	});
	body = body.replace(new RegExp(suffix, "g"), "");
	return body;
}


// EXTERNAL MODULE: ../node_modules/d3/src/index.js + 99 modules
var src = __webpack_require__(35321);
// EXTERNAL MODULE: consume shared module (default) marked@^16.2.1 (strict) (fallback: ../node_modules/marked/lib/marked.esm.js)
var marked_esm_js_ = __webpack_require__(95486);
// EXTERNAL MODULE: ../node_modules/ts-dedent/esm/index.js
var esm = __webpack_require__(11464);
;// CONCATENATED MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-JA3XYJ7Z.mjs




// src/rendering-util/icons.ts

var unknownIcon = {
  body: '<g><rect width="80" height="80" style="fill: #087ebf; stroke-width: 0px;"/><text transform="translate(21.16 64.67)" style="fill: #fff; font-family: ArialMT, Arial; font-size: 67.75px;"><tspan x="0" y="0">?</tspan></text></g>',
  height: 80,
  width: 80
};
var iconsStore = /* @__PURE__ */ new Map();
var loaderStore = /* @__PURE__ */ new Map();
var registerIconPacks = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((iconLoaders) => {
  for (const iconLoader of iconLoaders) {
    if (!iconLoader.name) {
      throw new Error(
        'Invalid icon loader. Must have a "name" property with non-empty string value.'
      );
    }
    chunk_AGHRB4JF/* log */.cM.debug("Registering icon pack:", iconLoader.name);
    if ("loader" in iconLoader) {
      loaderStore.set(iconLoader.name, iconLoader.loader);
    } else if ("icons" in iconLoader) {
      iconsStore.set(iconLoader.name, iconLoader.icons);
    } else {
      chunk_AGHRB4JF/* log */.cM.error("Invalid icon loader:", iconLoader);
      throw new Error('Invalid icon loader. Must have either "icons" or "loader" property.');
    }
  }
}, "registerIconPacks");
var getRegisteredIconData = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async (iconName, fallbackPrefix) => {
  const data = stringToIcon(iconName, true, fallbackPrefix !== void 0);
  if (!data) {
    throw new Error(`Invalid icon name: ${iconName}`);
  }
  const prefix = data.prefix || fallbackPrefix;
  if (!prefix) {
    throw new Error(`Icon name must contain a prefix: ${iconName}`);
  }
  let icons = iconsStore.get(prefix);
  if (!icons) {
    const loader = loaderStore.get(prefix);
    if (!loader) {
      throw new Error(`Icon set not found: ${data.prefix}`);
    }
    try {
      const loaded = await loader();
      icons = { ...loaded, prefix };
      iconsStore.set(prefix, icons);
    } catch (e) {
      chunk_AGHRB4JF/* log */.cM.error(e);
      throw new Error(`Failed to load icon set: ${data.prefix}`);
    }
  }
  const iconData = getIconData(icons, data.name);
  if (!iconData) {
    throw new Error(`Icon not found: ${iconName}`);
  }
  return iconData;
}, "getRegisteredIconData");
var isIconAvailable = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async (iconName) => {
  try {
    await getRegisteredIconData(iconName);
    return true;
  } catch {
    return false;
  }
}, "isIconAvailable");
var getIconSVG = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async (iconName, customisations, extraAttributes) => {
  let iconData;
  try {
    iconData = await getRegisteredIconData(iconName, customisations?.fallbackPrefix);
  } catch (e) {
    chunk_AGHRB4JF/* log */.cM.error(e);
    iconData = unknownIcon;
  }
  const renderData = iconToSVG(iconData, customisations);
  const svg = iconToHTML(replaceIDs(renderData.body), {
    ...renderData.attributes,
    ...extraAttributes
  });
  return (0,chunk_ABZYJK2D/* sanitizeText */.oO)(svg, (0,chunk_ABZYJK2D/* getConfig */.iE)());
}, "getIconSVG");

// src/rendering-util/createText.ts


// src/rendering-util/handle-markdown-text.ts


function preprocessMarkdown(markdown, { markdownAutoWrap }) {
  const withoutBR = markdown.replace(/<br\/>/g, "\n");
  const withoutMultipleNewlines = withoutBR.replace(/\n{2,}/g, "\n");
  const withoutExtraSpaces = (0,esm/* dedent */.Z)(withoutMultipleNewlines);
  if (markdownAutoWrap === false) {
    return withoutExtraSpaces.replace(/ /g, "&nbsp;");
  }
  return withoutExtraSpaces;
}
(0,chunk_AGHRB4JF/* __name */.eW)(preprocessMarkdown, "preprocessMarkdown");
function markdownToLines(markdown, config = {}) {
  const preprocessedMarkdown = preprocessMarkdown(markdown, config);
  const nodes = marked_esm_js_.marked.lexer(preprocessedMarkdown);
  const lines = [[]];
  let currentLine = 0;
  function processNode(node, parentType = "normal") {
    if (node.type === "text") {
      const textLines = node.text.split("\n");
      textLines.forEach((textLine, index) => {
        if (index !== 0) {
          currentLine++;
          lines.push([]);
        }
        textLine.split(" ").forEach((word) => {
          word = word.replace(/&#39;/g, `'`);
          if (word) {
            lines[currentLine].push({ content: word, type: parentType });
          }
        });
      });
    } else if (node.type === "strong" || node.type === "em") {
      node.tokens.forEach((contentNode) => {
        processNode(contentNode, node.type);
      });
    } else if (node.type === "html") {
      lines[currentLine].push({ content: node.text, type: "normal" });
    }
  }
  (0,chunk_AGHRB4JF/* __name */.eW)(processNode, "processNode");
  nodes.forEach((treeNode) => {
    if (treeNode.type === "paragraph") {
      treeNode.tokens?.forEach((contentNode) => {
        processNode(contentNode);
      });
    } else if (treeNode.type === "html") {
      lines[currentLine].push({ content: treeNode.text, type: "normal" });
    } else {
      lines[currentLine].push({ content: treeNode.raw, type: "normal" });
    }
  });
  return lines;
}
(0,chunk_AGHRB4JF/* __name */.eW)(markdownToLines, "markdownToLines");
function markdownToHTML(markdown, { markdownAutoWrap } = {}) {
  const nodes = marked_esm_js_.marked.lexer(markdown);
  function output(node) {
    if (node.type === "text") {
      if (markdownAutoWrap === false) {
        return node.text.replace(/\n */g, "<br/>").replace(/ /g, "&nbsp;");
      }
      return node.text.replace(/\n */g, "<br/>");
    } else if (node.type === "strong") {
      return `<strong>${node.tokens?.map(output).join("")}</strong>`;
    } else if (node.type === "em") {
      return `<em>${node.tokens?.map(output).join("")}</em>`;
    } else if (node.type === "paragraph") {
      return `<p>${node.tokens?.map(output).join("")}</p>`;
    } else if (node.type === "space") {
      return "";
    } else if (node.type === "html") {
      return `${node.text}`;
    } else if (node.type === "escape") {
      return node.text;
    }
    chunk_AGHRB4JF/* log */.cM.warn(`Unsupported markdown: ${node.type}`);
    return node.raw;
  }
  (0,chunk_AGHRB4JF/* __name */.eW)(output, "output");
  return nodes.map(output).join("");
}
(0,chunk_AGHRB4JF/* __name */.eW)(markdownToHTML, "markdownToHTML");

// src/rendering-util/splitText.ts
function splitTextToChars(text) {
  if (Intl.Segmenter) {
    return [...new Intl.Segmenter().segment(text)].map((s) => s.segment);
  }
  return [...text];
}
(0,chunk_AGHRB4JF/* __name */.eW)(splitTextToChars, "splitTextToChars");
function splitWordToFitWidth(checkFit, word) {
  const characters = splitTextToChars(word.content);
  return splitWordToFitWidthRecursion(checkFit, [], characters, word.type);
}
(0,chunk_AGHRB4JF/* __name */.eW)(splitWordToFitWidth, "splitWordToFitWidth");
function splitWordToFitWidthRecursion(checkFit, usedChars, remainingChars, type) {
  if (remainingChars.length === 0) {
    return [
      { content: usedChars.join(""), type },
      { content: "", type }
    ];
  }
  const [nextChar, ...rest] = remainingChars;
  const newWord = [...usedChars, nextChar];
  if (checkFit([{ content: newWord.join(""), type }])) {
    return splitWordToFitWidthRecursion(checkFit, newWord, rest, type);
  }
  if (usedChars.length === 0 && nextChar) {
    usedChars.push(nextChar);
    remainingChars.shift();
  }
  return [
    { content: usedChars.join(""), type },
    { content: remainingChars.join(""), type }
  ];
}
(0,chunk_AGHRB4JF/* __name */.eW)(splitWordToFitWidthRecursion, "splitWordToFitWidthRecursion");
function splitLineToFitWidth(line, checkFit) {
  if (line.some(({ content }) => content.includes("\n"))) {
    throw new Error("splitLineToFitWidth does not support newlines in the line");
  }
  return splitLineToFitWidthRecursion(line, checkFit);
}
(0,chunk_AGHRB4JF/* __name */.eW)(splitLineToFitWidth, "splitLineToFitWidth");
function splitLineToFitWidthRecursion(words, checkFit, lines = [], newLine = []) {
  if (words.length === 0) {
    if (newLine.length > 0) {
      lines.push(newLine);
    }
    return lines.length > 0 ? lines : [];
  }
  let joiner = "";
  if (words[0].content === " ") {
    joiner = " ";
    words.shift();
  }
  const nextWord = words.shift() ?? { content: " ", type: "normal" };
  const lineWithNextWord = [...newLine];
  if (joiner !== "") {
    lineWithNextWord.push({ content: joiner, type: "normal" });
  }
  lineWithNextWord.push(nextWord);
  if (checkFit(lineWithNextWord)) {
    return splitLineToFitWidthRecursion(words, checkFit, lines, lineWithNextWord);
  }
  if (newLine.length > 0) {
    lines.push(newLine);
    words.unshift(nextWord);
  } else if (nextWord.content) {
    const [line, rest] = splitWordToFitWidth(checkFit, nextWord);
    lines.push([line]);
    if (rest.content) {
      words.unshift(rest);
    }
  }
  return splitLineToFitWidthRecursion(words, checkFit, lines);
}
(0,chunk_AGHRB4JF/* __name */.eW)(splitLineToFitWidthRecursion, "splitLineToFitWidthRecursion");

// src/rendering-util/createText.ts
function applyStyle(dom, styleFn) {
  if (styleFn) {
    dom.attr("style", styleFn);
  }
}
(0,chunk_AGHRB4JF/* __name */.eW)(applyStyle, "applyStyle");
async function addHtmlSpan(element, node, width, classes, addBackground = false, config = (0,chunk_ABZYJK2D/* getConfig */.iE)()) {
  const fo = element.append("foreignObject");
  fo.attr("width", `${10 * width}px`);
  fo.attr("height", `${10 * width}px`);
  const div = fo.append("xhtml:div");
  const sanitizedLabel = (0,chunk_ABZYJK2D/* hasKatex */.l0)(node.label) ? await (0,chunk_ABZYJK2D/* renderKatexSanitized */.ur)(node.label.replace(chunk_ABZYJK2D/* common_default */.SY.lineBreakRegex, "\n"), config) : (0,chunk_ABZYJK2D/* sanitizeText */.oO)(node.label, config);
  const labelClass = node.isNode ? "nodeLabel" : "edgeLabel";
  const span = div.append("span");
  span.html(sanitizedLabel);
  applyStyle(span, node.labelStyle);
  span.attr("class", `${labelClass} ${classes}`);
  applyStyle(div, node.labelStyle);
  div.style("display", "table-cell");
  div.style("white-space", "nowrap");
  div.style("line-height", "1.5");
  div.style("max-width", width + "px");
  div.style("text-align", "center");
  div.attr("xmlns", "http://www.w3.org/1999/xhtml");
  if (addBackground) {
    div.attr("class", "labelBkg");
  }
  let bbox = div.node().getBoundingClientRect();
  if (bbox.width === width) {
    div.style("display", "table");
    div.style("white-space", "break-spaces");
    div.style("width", width + "px");
    bbox = div.node().getBoundingClientRect();
  }
  return fo.node();
}
(0,chunk_AGHRB4JF/* __name */.eW)(addHtmlSpan, "addHtmlSpan");
function createTspan(textElement, lineIndex, lineHeight) {
  return textElement.append("tspan").attr("class", "text-outer-tspan").attr("x", 0).attr("y", lineIndex * lineHeight - 0.1 + "em").attr("dy", lineHeight + "em");
}
(0,chunk_AGHRB4JF/* __name */.eW)(createTspan, "createTspan");
function computeWidthOfText(parentNode, lineHeight, line) {
  const testElement = parentNode.append("text");
  const testSpan = createTspan(testElement, 1, lineHeight);
  updateTextContentAndStyles(testSpan, line);
  const textLength = testSpan.node().getComputedTextLength();
  testElement.remove();
  return textLength;
}
(0,chunk_AGHRB4JF/* __name */.eW)(computeWidthOfText, "computeWidthOfText");
function computeDimensionOfText(parentNode, lineHeight, text) {
  const testElement = parentNode.append("text");
  const testSpan = createTspan(testElement, 1, lineHeight);
  updateTextContentAndStyles(testSpan, [{ content: text, type: "normal" }]);
  const textDimension = testSpan.node()?.getBoundingClientRect();
  if (textDimension) {
    testElement.remove();
  }
  return textDimension;
}
(0,chunk_AGHRB4JF/* __name */.eW)(computeDimensionOfText, "computeDimensionOfText");
function createFormattedText(width, g, structuredText, addBackground = false) {
  const lineHeight = 1.1;
  const labelGroup = g.append("g");
  const bkg = labelGroup.insert("rect").attr("class", "background").attr("style", "stroke: none");
  const textElement = labelGroup.append("text").attr("y", "-10.1");
  let lineIndex = 0;
  for (const line of structuredText) {
    const checkWidth = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((line2) => computeWidthOfText(labelGroup, lineHeight, line2) <= width, "checkWidth");
    const linesUnderWidth = checkWidth(line) ? [line] : splitLineToFitWidth(line, checkWidth);
    for (const preparedLine of linesUnderWidth) {
      const tspan = createTspan(textElement, lineIndex, lineHeight);
      updateTextContentAndStyles(tspan, preparedLine);
      lineIndex++;
    }
  }
  if (addBackground) {
    const bbox = textElement.node().getBBox();
    const padding = 2;
    bkg.attr("x", bbox.x - padding).attr("y", bbox.y - padding).attr("width", bbox.width + 2 * padding).attr("height", bbox.height + 2 * padding);
    return labelGroup.node();
  } else {
    return textElement.node();
  }
}
(0,chunk_AGHRB4JF/* __name */.eW)(createFormattedText, "createFormattedText");
function updateTextContentAndStyles(tspan, wrappedLine) {
  tspan.text("");
  wrappedLine.forEach((word, index) => {
    const innerTspan = tspan.append("tspan").attr("font-style", word.type === "em" ? "italic" : "normal").attr("class", "text-inner-tspan").attr("font-weight", word.type === "strong" ? "bold" : "normal");
    if (index === 0) {
      innerTspan.text(word.content);
    } else {
      innerTspan.text(" " + word.content);
    }
  });
}
(0,chunk_AGHRB4JF/* __name */.eW)(updateTextContentAndStyles, "updateTextContentAndStyles");
async function replaceIconSubstring(text, config = {}) {
  const pendingReplacements = [];
  text.replace(/(fa[bklrs]?):fa-([\w-]+)/g, (fullMatch, prefix, iconName) => {
    pendingReplacements.push(
      (async () => {
        const registeredIconName = `${prefix}:${iconName}`;
        if (await isIconAvailable(registeredIconName)) {
          return await getIconSVG(registeredIconName, void 0, { class: "label-icon" });
        } else {
          return `<i class='${(0,chunk_ABZYJK2D/* sanitizeText */.oO)(fullMatch, config).replace(":", " ")}'></i>`;
        }
      })()
    );
    return fullMatch;
  });
  const replacements = await Promise.all(pendingReplacements);
  return text.replace(/(fa[bklrs]?):fa-([\w-]+)/g, () => replacements.shift() ?? "");
}
(0,chunk_AGHRB4JF/* __name */.eW)(replaceIconSubstring, "replaceIconSubstring");
var createText = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async (el, text = "", {
  style = "",
  isTitle = false,
  classes = "",
  useHtmlLabels = true,
  isNode = true,
  width = 200,
  addSvgBackground = false
} = {}, config) => {
  chunk_AGHRB4JF/* log */.cM.debug(
    "XYZ createText",
    text,
    style,
    isTitle,
    classes,
    useHtmlLabels,
    isNode,
    "addSvgBackground: ",
    addSvgBackground
  );
  if (useHtmlLabels) {
    const htmlText = markdownToHTML(text, config);
    const decodedReplacedText = await replaceIconSubstring((0,chunk_S3R3BYOJ/* decodeEntities */.SH)(htmlText), config);
    const inputForKatex = text.replace(/\\\\/g, "\\");
    const node = {
      isNode,
      label: (0,chunk_ABZYJK2D/* hasKatex */.l0)(text) ? inputForKatex : decodedReplacedText,
      labelStyle: style.replace("fill:", "color:")
    };
    const vertexNode = await addHtmlSpan(el, node, width, classes, addSvgBackground, config);
    return vertexNode;
  } else {
    const sanitizeBR = text.replace(/<br\s*\/?>/g, "<br/>");
    const structuredText = markdownToLines(sanitizeBR.replace("<br>", "<br/>"), config);
    const svgLabel = createFormattedText(
      width,
      el,
      structuredText,
      text ? addSvgBackground : false
    );
    if (isNode) {
      if (/stroke:/.exec(style)) {
        style = style.replace("stroke:", "lineColor:");
      }
      const nodeLabelTextStyle = style.replace(/stroke:[^;]+;?/g, "").replace(/stroke-width:[^;]+;?/g, "").replace(/fill:[^;]+;?/g, "").replace(/color:/g, "fill:");
      (0,src/* select */.Ys)(svgLabel).attr("style", nodeLabelTextStyle);
    } else {
      const edgeLabelRectStyle = style.replace(/stroke:[^;]+;?/g, "").replace(/stroke-width:[^;]+;?/g, "").replace(/fill:[^;]+;?/g, "").replace(/background:/g, "fill:");
      (0,src/* select */.Ys)(svgLabel).select("rect").attr("style", edgeLabelRectStyle.replace(/background:/g, "fill:"));
      const edgeLabelTextStyle = style.replace(/stroke:[^;]+;?/g, "").replace(/stroke-width:[^;]+;?/g, "").replace(/fill:[^;]+;?/g, "").replace(/color:/g, "fill:");
      (0,src/* select */.Ys)(svgLabel).select("text").attr("style", edgeLabelTextStyle);
    }
    return svgLabel;
  }
}, "createText");




/***/ }),

/***/ 26212:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   C1: () => (/* binding */ labelHelper),
/* harmony export */   Lf: () => (/* binding */ insertNode),
/* harmony export */   XO: () => (/* binding */ createLabel_default),
/* harmony export */   Yn: () => (/* binding */ setNodeElem),
/* harmony export */   ZH: () => (/* binding */ clear),
/* harmony export */   aH: () => (/* binding */ positionNode),
/* harmony export */   dW: () => (/* binding */ isValidShape),
/* harmony export */   gU: () => (/* binding */ clear2),
/* harmony export */   jr: () => (/* binding */ updateNodeBounds),
/* harmony export */   us: () => (/* binding */ insertCluster)
/* harmony export */ });
/* harmony import */ var _chunk_CVBHYZKI_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(22462);
/* harmony import */ var _chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(37756);
/* harmony import */ var _chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(2111);
/* harmony import */ var _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(17175);
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(74999);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(35321);
/* harmony import */ var roughjs__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(36366);







// src/rendering-util/rendering-elements/shapes/util.ts

var labelHelper = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(async (parent, node, _classes) => {
  let cssClasses;
  const useHtmlLabels = node.useHtmlLabels || (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)()?.htmlLabels);
  if (!_classes) {
    cssClasses = "node default";
  } else {
    cssClasses = _classes;
  }
  const shapeSvg = parent.insert("g").attr("class", cssClasses).attr("id", node.domId || node.id);
  const labelEl = shapeSvg.insert("g").attr("class", "label").attr("style", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(node.labelStyle));
  let label;
  if (node.label === void 0) {
    label = "";
  } else {
    label = typeof node.label === "string" ? node.label : node.label[0];
  }
  const text2 = await (0,_chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_2__/* .createText */ .rw)(labelEl, (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .sanitizeText */ .oO)((0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .decodeEntities */ .SH)(label), (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)()), {
    useHtmlLabels,
    width: node.width || (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)().flowchart?.wrappingWidth,
    // @ts-expect-error -- This is currently not used. Should this be `classes` instead?
    cssClasses: "markdown-node-label",
    style: node.labelStyle,
    addSvgBackground: !!node.icon || !!node.img
  });
  let bbox = text2.getBBox();
  const halfPadding = (node?.padding ?? 0) / 2;
  if (useHtmlLabels) {
    const div = text2.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(text2);
    const images = div.getElementsByTagName("img");
    if (images) {
      const noImgText = label.replace(/<img[^>]*>/g, "").trim() === "";
      await Promise.all(
        [...images].map(
          (img) => new Promise((res) => {
            function setupImage() {
              img.style.display = "flex";
              img.style.flexDirection = "column";
              if (noImgText) {
                const bodyFontSize = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)().fontSize ? (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)().fontSize : window.getComputedStyle(document.body).fontSize;
                const enlargingFactor = 5;
                const [parsedBodyFontSize = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .defaultConfig_default */ .vZ.fontSize] = (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .parseFontSize */ .VG)(bodyFontSize);
                const width = parsedBodyFontSize * enlargingFactor + "px";
                img.style.minWidth = width;
                img.style.maxWidth = width;
              } else {
                img.style.width = "100%";
              }
              res(img);
            }
            (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(setupImage, "setupImage");
            setTimeout(() => {
              if (img.complete) {
                setupImage();
              }
            });
            img.addEventListener("error", setupImage);
            img.addEventListener("load", setupImage);
          })
        )
      );
    }
    bbox = div.getBoundingClientRect();
    dv.attr("width", bbox.width);
    dv.attr("height", bbox.height);
  }
  if (useHtmlLabels) {
    labelEl.attr("transform", "translate(" + -bbox.width / 2 + ", " + -bbox.height / 2 + ")");
  } else {
    labelEl.attr("transform", "translate(0, " + -bbox.height / 2 + ")");
  }
  if (node.centerLabel) {
    labelEl.attr("transform", "translate(" + -bbox.width / 2 + ", " + -bbox.height / 2 + ")");
  }
  labelEl.insert("rect", ":first-child");
  return { shapeSvg, bbox, halfPadding, label: labelEl };
}, "labelHelper");
var insertLabel = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(async (parent, label, options) => {
  const useHtmlLabels = options.useHtmlLabels || (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)()?.flowchart?.htmlLabels);
  const labelEl = parent.insert("g").attr("class", "label").attr("style", options.labelStyle || "");
  const text2 = await (0,_chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_2__/* .createText */ .rw)(labelEl, (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .sanitizeText */ .oO)((0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .decodeEntities */ .SH)(label), (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)()), {
    useHtmlLabels,
    width: options.width || (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)()?.flowchart?.wrappingWidth,
    style: options.labelStyle,
    addSvgBackground: !!options.icon || !!options.img
  });
  let bbox = text2.getBBox();
  const halfPadding = options.padding / 2;
  if ((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)()?.flowchart?.htmlLabels)) {
    const div = text2.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(text2);
    bbox = div.getBoundingClientRect();
    dv.attr("width", bbox.width);
    dv.attr("height", bbox.height);
  }
  if (useHtmlLabels) {
    labelEl.attr("transform", "translate(" + -bbox.width / 2 + ", " + -bbox.height / 2 + ")");
  } else {
    labelEl.attr("transform", "translate(0, " + -bbox.height / 2 + ")");
  }
  if (options.centerLabel) {
    labelEl.attr("transform", "translate(" + -bbox.width / 2 + ", " + -bbox.height / 2 + ")");
  }
  labelEl.insert("rect", ":first-child");
  return { shapeSvg: parent, bbox, halfPadding, label: labelEl };
}, "insertLabel");
var updateNodeBounds = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((node, element) => {
  const bbox = element.node().getBBox();
  node.width = bbox.width;
  node.height = bbox.height;
}, "updateNodeBounds");
var getNodeClasses = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((node, extra) => (node.look === "handDrawn" ? "rough-node" : "node") + " " + node.cssClasses + " " + (extra || ""), "getNodeClasses");
function createPathFromPoints(points) {
  const pointStrings = points.map((p, i) => `${i === 0 ? "M" : "L"}${p.x},${p.y}`);
  pointStrings.push("Z");
  return pointStrings.join(" ");
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(createPathFromPoints, "createPathFromPoints");
function generateFullSineWavePoints(x1, y1, x2, y2, amplitude, numCycles) {
  const points = [];
  const steps = 50;
  const deltaX = x2 - x1;
  const deltaY = y2 - y1;
  const cycleLength = deltaX / numCycles;
  const frequency = 2 * Math.PI / cycleLength;
  const midY = y1 + deltaY / 2;
  for (let i = 0; i <= steps; i++) {
    const t = i / steps;
    const x = x1 + t * deltaX;
    const y = midY + amplitude * Math.sin(frequency * (x - x1));
    points.push({ x, y });
  }
  return points;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(generateFullSineWavePoints, "generateFullSineWavePoints");
function generateCirclePoints(centerX, centerY, radius, numPoints, startAngle, endAngle) {
  const points = [];
  const startAngleRad = startAngle * Math.PI / 180;
  const endAngleRad = endAngle * Math.PI / 180;
  const angleRange = endAngleRad - startAngleRad;
  const angleStep = angleRange / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    const angle = startAngleRad + i * angleStep;
    const x = centerX + radius * Math.cos(angle);
    const y = centerY + radius * Math.sin(angle);
    points.push({ x: -x, y: -y });
  }
  return points;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(generateCirclePoints, "generateCirclePoints");

// src/rendering-util/rendering-elements/clusters.js



// src/rendering-util/rendering-elements/intersect/intersect-rect.js
var intersectRect = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((node, point) => {
  var x = node.x;
  var y = node.y;
  var dx = point.x - x;
  var dy = point.y - y;
  var w = node.width / 2;
  var h = node.height / 2;
  var sx, sy;
  if (Math.abs(dy) * w > Math.abs(dx) * h) {
    if (dy < 0) {
      h = -h;
    }
    sx = dy === 0 ? 0 : h * dx / dy;
    sy = h;
  } else {
    if (dx < 0) {
      w = -w;
    }
    sx = w;
    sy = dx === 0 ? 0 : w * dy / dx;
  }
  return { x: x + sx, y: y + sy };
}, "intersectRect");
var intersect_rect_default = intersectRect;

// src/rendering-util/rendering-elements/createLabel.js

function applyStyle(dom, styleFn) {
  if (styleFn) {
    dom.attr("style", styleFn);
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(applyStyle, "applyStyle");
async function addHtmlLabel(node) {
  const fo = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(document.createElementNS("http://www.w3.org/2000/svg", "foreignObject"));
  const div = fo.append("xhtml:div");
  const config = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)();
  let label = node.label;
  if (node.label && (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .hasKatex */ .l0)(node.label)) {
    label = await (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .renderKatexSanitized */ .ur)(node.label.replace(_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .common_default */ .SY.lineBreakRegex, "\n"), config);
  }
  const labelClass = node.isNode ? "nodeLabel" : "edgeLabel";
  const labelSpan = '<span class="' + labelClass + '" ' + (node.labelStyle ? 'style="' + node.labelStyle + '"' : "") + // codeql [js/html-constructed-from-input] : false positive
  ">" + label + "</span>";
  div.html((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .sanitizeText */ .oO)(labelSpan, config));
  applyStyle(div, node.labelStyle);
  div.style("display", "inline-block");
  div.style("padding-right", "1px");
  div.style("white-space", "nowrap");
  div.attr("xmlns", "http://www.w3.org/1999/xhtml");
  return fo.node();
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(addHtmlLabel, "addHtmlLabel");
var createLabel = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(async (_vertexText, style, isTitle, isNode) => {
  let vertexText = _vertexText || "";
  if (typeof vertexText === "object") {
    vertexText = vertexText[0];
  }
  if ((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)().flowchart.htmlLabels)) {
    vertexText = vertexText.replace(/\\n|\n/g, "<br />");
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("vertexText" + vertexText);
    const node = {
      isNode,
      label: (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .decodeEntities */ .SH)(vertexText).replace(
        /fa[blrs]?:fa-[\w-]+/g,
        (s) => `<i class='${s.replace(":", " ")}'></i>`
      ),
      labelStyle: style ? style.replace("fill:", "color:") : style
    };
    let vertexNode = await addHtmlLabel(node);
    return vertexNode;
  } else {
    const svgLabel = document.createElementNS("http://www.w3.org/2000/svg", "text");
    svgLabel.setAttribute("style", style.replace("color:", "fill:"));
    let rows = [];
    if (typeof vertexText === "string") {
      rows = vertexText.split(/\\n|\n|<br\s*\/?>/gi);
    } else if (Array.isArray(vertexText)) {
      rows = vertexText;
    } else {
      rows = [];
    }
    for (const row of rows) {
      const tspan = document.createElementNS("http://www.w3.org/2000/svg", "tspan");
      tspan.setAttributeNS("http://www.w3.org/XML/1998/namespace", "xml:space", "preserve");
      tspan.setAttribute("dy", "1em");
      tspan.setAttribute("x", "0");
      if (isTitle) {
        tspan.setAttribute("class", "title-row");
      } else {
        tspan.setAttribute("class", "row");
      }
      tspan.textContent = row.trim();
      svgLabel.appendChild(tspan);
    }
    return svgLabel;
  }
}, "createLabel");
var createLabel_default = createLabel;

// src/rendering-util/rendering-elements/shapes/roundedRectPath.ts
var createRoundedRectPathD = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((x, y, totalWidth, totalHeight, radius) => [
  "M",
  x + radius,
  y,
  // Move to the first point
  "H",
  x + totalWidth - radius,
  // Draw horizontal line to the beginning of the right corner
  "A",
  radius,
  radius,
  0,
  0,
  1,
  x + totalWidth,
  y + radius,
  // Draw arc to the right top corner
  "V",
  y + totalHeight - radius,
  // Draw vertical line down to the beginning of the right bottom corner
  "A",
  radius,
  radius,
  0,
  0,
  1,
  x + totalWidth - radius,
  y + totalHeight,
  // Draw arc to the right bottom corner
  "H",
  x + radius,
  // Draw horizontal line to the beginning of the left bottom corner
  "A",
  radius,
  radius,
  0,
  0,
  1,
  x,
  y + totalHeight - radius,
  // Draw arc to the left bottom corner
  "V",
  y + radius,
  // Draw vertical line up to the beginning of the left top corner
  "A",
  radius,
  radius,
  0,
  0,
  1,
  x + radius,
  y,
  // Draw arc to the left top corner
  "Z"
  // Close the path
].join(" "), "createRoundedRectPathD");

// src/rendering-util/rendering-elements/clusters.js
var rect = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(async (parent, node) => {
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("Creating subgraph rect for ", node.id, node);
  const siteConfig = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)();
  const { themeVariables, handDrawnSeed } = siteConfig;
  const { clusterBkg, clusterBorder } = themeVariables;
  const { labelStyles, nodeStyles, borderStyles, backgroundStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  const shapeSvg = parent.insert("g").attr("class", "cluster " + node.cssClasses).attr("id", node.id).attr("data-look", node.look);
  const useHtmlLabels = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)(siteConfig.flowchart.htmlLabels);
  const labelEl = shapeSvg.insert("g").attr("class", "cluster-label ");
  const text2 = await (0,_chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_2__/* .createText */ .rw)(labelEl, node.label, {
    style: node.labelStyle,
    useHtmlLabels,
    isNode: true
  });
  let bbox = text2.getBBox();
  if ((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)(siteConfig.flowchart.htmlLabels)) {
    const div = text2.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(text2);
    bbox = div.getBoundingClientRect();
    dv.attr("width", bbox.width);
    dv.attr("height", bbox.height);
  }
  const width = node.width <= bbox.width + node.padding ? bbox.width + node.padding : node.width;
  if (node.width <= bbox.width + node.padding) {
    node.diff = (width - node.width) / 2 - node.padding;
  } else {
    node.diff = -node.padding;
  }
  const height = node.height;
  const x = node.x - width / 2;
  const y = node.y - height / 2;
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.trace("Data ", node, JSON.stringify(node));
  let rect2;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {
      roughness: 0.7,
      fill: clusterBkg,
      // fill: 'red',
      stroke: clusterBorder,
      fillWeight: 3,
      seed: handDrawnSeed
    });
    const roughNode = rc.path(createRoundedRectPathD(x, y, width, height, 0), options);
    rect2 = shapeSvg.insert(() => {
      _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.debug("Rough node insert CXC", roughNode);
      return roughNode;
    }, ":first-child");
    rect2.select("path:nth-child(2)").attr("style", borderStyles.join(";"));
    rect2.select("path").attr("style", backgroundStyles.join(";").replace("fill", "stroke"));
  } else {
    rect2 = shapeSvg.insert("rect", ":first-child");
    rect2.attr("style", nodeStyles).attr("rx", node.rx).attr("ry", node.ry).attr("x", x).attr("y", y).attr("width", width).attr("height", height);
  }
  const { subGraphTitleTopMargin } = (0,_chunk_CVBHYZKI_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getSubGraphTitleMargins */ .L)(siteConfig);
  labelEl.attr(
    "transform",
    // This puts the label on top of the box instead of inside it
    `translate(${node.x - bbox.width / 2}, ${node.y - node.height / 2 + subGraphTitleTopMargin})`
  );
  if (labelStyles) {
    const span = labelEl.select("span");
    if (span) {
      span.attr("style", labelStyles);
    }
  }
  const rectBox = rect2.node().getBBox();
  node.offsetX = 0;
  node.width = rectBox.width;
  node.height = rectBox.height;
  node.offsetY = bbox.height - node.padding / 2;
  node.intersect = function(point) {
    return intersect_rect_default(node, point);
  };
  return { cluster: shapeSvg, labelBBox: bbox };
}, "rect");
var noteGroup = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((parent, node) => {
  const shapeSvg = parent.insert("g").attr("class", "note-cluster").attr("id", node.id);
  const rect2 = shapeSvg.insert("rect", ":first-child");
  const padding = 0 * node.padding;
  const halfPadding = padding / 2;
  rect2.attr("rx", node.rx).attr("ry", node.ry).attr("x", node.x - node.width / 2 - halfPadding).attr("y", node.y - node.height / 2 - halfPadding).attr("width", node.width + padding).attr("height", node.height + padding).attr("fill", "none");
  const rectBox = rect2.node().getBBox();
  node.width = rectBox.width;
  node.height = rectBox.height;
  node.intersect = function(point) {
    return intersect_rect_default(node, point);
  };
  return { cluster: shapeSvg, labelBBox: { width: 0, height: 0 } };
}, "noteGroup");
var roundedWithTitle = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(async (parent, node) => {
  const siteConfig = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)();
  const { themeVariables, handDrawnSeed } = siteConfig;
  const { altBackground, compositeBackground, compositeTitleBackground, nodeBorder } = themeVariables;
  const shapeSvg = parent.insert("g").attr("class", node.cssClasses).attr("id", node.id).attr("data-id", node.id).attr("data-look", node.look);
  const outerRectG = shapeSvg.insert("g", ":first-child");
  const label = shapeSvg.insert("g").attr("class", "cluster-label");
  let innerRect = shapeSvg.append("rect");
  const text2 = label.node().appendChild(await createLabel_default(node.label, node.labelStyle, void 0, true));
  let bbox = text2.getBBox();
  if ((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)(siteConfig.flowchart.htmlLabels)) {
    const div = text2.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(text2);
    bbox = div.getBoundingClientRect();
    dv.attr("width", bbox.width);
    dv.attr("height", bbox.height);
  }
  const padding = 0 * node.padding;
  const halfPadding = padding / 2;
  const width = (node.width <= bbox.width + node.padding ? bbox.width + node.padding : node.width) + padding;
  if (node.width <= bbox.width + node.padding) {
    node.diff = (width - node.width) / 2 - node.padding;
  } else {
    node.diff = -node.padding;
  }
  const height = node.height + padding;
  const innerHeight = node.height + padding - bbox.height - 6;
  const x = node.x - width / 2;
  const y = node.y - height / 2;
  node.width = width;
  const innerY = node.y - node.height / 2 - halfPadding + bbox.height + 2;
  let rect2;
  if (node.look === "handDrawn") {
    const isAlt = node.cssClasses.includes("statediagram-cluster-alt");
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const roughOuterNode = node.rx || node.ry ? rc.path(createRoundedRectPathD(x, y, width, height, 10), {
      roughness: 0.7,
      fill: compositeTitleBackground,
      fillStyle: "solid",
      stroke: nodeBorder,
      seed: handDrawnSeed
    }) : rc.rectangle(x, y, width, height, { seed: handDrawnSeed });
    rect2 = shapeSvg.insert(() => roughOuterNode, ":first-child");
    const roughInnerNode = rc.rectangle(x, innerY, width, innerHeight, {
      fill: isAlt ? altBackground : compositeBackground,
      fillStyle: isAlt ? "hachure" : "solid",
      stroke: nodeBorder,
      seed: handDrawnSeed
    });
    rect2 = shapeSvg.insert(() => roughOuterNode, ":first-child");
    innerRect = shapeSvg.insert(() => roughInnerNode);
  } else {
    rect2 = outerRectG.insert("rect", ":first-child");
    const outerRectClass = "outer";
    rect2.attr("class", outerRectClass).attr("x", x).attr("y", y).attr("width", width).attr("height", height).attr("data-look", node.look);
    innerRect.attr("class", "inner").attr("x", x).attr("y", innerY).attr("width", width).attr("height", innerHeight);
  }
  label.attr(
    "transform",
    `translate(${node.x - bbox.width / 2}, ${y + 1 - ((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)(siteConfig.flowchart.htmlLabels) ? 0 : 3)})`
  );
  const rectBox = rect2.node().getBBox();
  node.height = rectBox.height;
  node.offsetX = 0;
  node.offsetY = bbox.height - node.padding / 2;
  node.labelBBox = bbox;
  node.intersect = function(point) {
    return intersect_rect_default(node, point);
  };
  return { cluster: shapeSvg, labelBBox: bbox };
}, "roundedWithTitle");
var kanbanSection = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(async (parent, node) => {
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("Creating subgraph rect for ", node.id, node);
  const siteConfig = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)();
  const { themeVariables, handDrawnSeed } = siteConfig;
  const { clusterBkg, clusterBorder } = themeVariables;
  const { labelStyles, nodeStyles, borderStyles, backgroundStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  const shapeSvg = parent.insert("g").attr("class", "cluster " + node.cssClasses).attr("id", node.id).attr("data-look", node.look);
  const useHtmlLabels = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)(siteConfig.flowchart.htmlLabels);
  const labelEl = shapeSvg.insert("g").attr("class", "cluster-label ");
  const text2 = await (0,_chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_2__/* .createText */ .rw)(labelEl, node.label, {
    style: node.labelStyle,
    useHtmlLabels,
    isNode: true,
    width: node.width
  });
  let bbox = text2.getBBox();
  if ((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)(siteConfig.flowchart.htmlLabels)) {
    const div = text2.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(text2);
    bbox = div.getBoundingClientRect();
    dv.attr("width", bbox.width);
    dv.attr("height", bbox.height);
  }
  const width = node.width <= bbox.width + node.padding ? bbox.width + node.padding : node.width;
  if (node.width <= bbox.width + node.padding) {
    node.diff = (width - node.width) / 2 - node.padding;
  } else {
    node.diff = -node.padding;
  }
  const height = node.height;
  const x = node.x - width / 2;
  const y = node.y - height / 2;
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.trace("Data ", node, JSON.stringify(node));
  let rect2;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {
      roughness: 0.7,
      fill: clusterBkg,
      // fill: 'red',
      stroke: clusterBorder,
      fillWeight: 4,
      seed: handDrawnSeed
    });
    const roughNode = rc.path(createRoundedRectPathD(x, y, width, height, node.rx), options);
    rect2 = shapeSvg.insert(() => {
      _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.debug("Rough node insert CXC", roughNode);
      return roughNode;
    }, ":first-child");
    rect2.select("path:nth-child(2)").attr("style", borderStyles.join(";"));
    rect2.select("path").attr("style", backgroundStyles.join(";").replace("fill", "stroke"));
  } else {
    rect2 = shapeSvg.insert("rect", ":first-child");
    rect2.attr("style", nodeStyles).attr("rx", node.rx).attr("ry", node.ry).attr("x", x).attr("y", y).attr("width", width).attr("height", height);
  }
  const { subGraphTitleTopMargin } = (0,_chunk_CVBHYZKI_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getSubGraphTitleMargins */ .L)(siteConfig);
  labelEl.attr(
    "transform",
    // This puts the label on top of the box instead of inside it
    `translate(${node.x - bbox.width / 2}, ${node.y - node.height / 2 + subGraphTitleTopMargin})`
  );
  if (labelStyles) {
    const span = labelEl.select("span");
    if (span) {
      span.attr("style", labelStyles);
    }
  }
  const rectBox = rect2.node().getBBox();
  node.offsetX = 0;
  node.width = rectBox.width;
  node.height = rectBox.height;
  node.offsetY = bbox.height - node.padding / 2;
  node.intersect = function(point) {
    return intersect_rect_default(node, point);
  };
  return { cluster: shapeSvg, labelBBox: bbox };
}, "kanbanSection");
var divider = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((parent, node) => {
  const siteConfig = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)();
  const { themeVariables, handDrawnSeed } = siteConfig;
  const { nodeBorder } = themeVariables;
  const shapeSvg = parent.insert("g").attr("class", node.cssClasses).attr("id", node.id).attr("data-look", node.look);
  const outerRectG = shapeSvg.insert("g", ":first-child");
  const padding = 0 * node.padding;
  const width = node.width + padding;
  node.diff = -node.padding;
  const height = node.height + padding;
  const x = node.x - width / 2;
  const y = node.y - height / 2;
  node.width = width;
  let rect2;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const roughOuterNode = rc.rectangle(x, y, width, height, {
      fill: "lightgrey",
      roughness: 0.5,
      strokeLineDash: [5],
      stroke: nodeBorder,
      seed: handDrawnSeed
    });
    rect2 = shapeSvg.insert(() => roughOuterNode, ":first-child");
  } else {
    rect2 = outerRectG.insert("rect", ":first-child");
    const outerRectClass = "divider";
    rect2.attr("class", outerRectClass).attr("x", x).attr("y", y).attr("width", width).attr("height", height).attr("data-look", node.look);
  }
  const rectBox = rect2.node().getBBox();
  node.height = rectBox.height;
  node.offsetX = 0;
  node.offsetY = 0;
  node.intersect = function(point) {
    return intersect_rect_default(node, point);
  };
  return { cluster: shapeSvg, labelBBox: {} };
}, "divider");
var squareRect = rect;
var shapes = {
  rect,
  squareRect,
  roundedWithTitle,
  noteGroup,
  divider,
  kanbanSection
};
var clusterElems = /* @__PURE__ */ new Map();
var insertCluster = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(async (elem, node) => {
  const shape = node.shape || "rect";
  const cluster = await shapes[shape](elem, node);
  clusterElems.set(node.id, cluster);
  return cluster;
}, "insertCluster");
var clear = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(() => {
  clusterElems = /* @__PURE__ */ new Map();
}, "clear");

// src/rendering-util/rendering-elements/intersect/intersect-node.js
function intersectNode(node, point) {
  return node.intersect(point);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(intersectNode, "intersectNode");
var intersect_node_default = intersectNode;

// src/rendering-util/rendering-elements/intersect/intersect-ellipse.js
function intersectEllipse(node, rx, ry, point) {
  var cx = node.x;
  var cy = node.y;
  var px = cx - point.x;
  var py = cy - point.y;
  var det = Math.sqrt(rx * rx * py * py + ry * ry * px * px);
  var dx = Math.abs(rx * ry * px / det);
  if (point.x < cx) {
    dx = -dx;
  }
  var dy = Math.abs(rx * ry * py / det);
  if (point.y < cy) {
    dy = -dy;
  }
  return { x: cx + dx, y: cy + dy };
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(intersectEllipse, "intersectEllipse");
var intersect_ellipse_default = intersectEllipse;

// src/rendering-util/rendering-elements/intersect/intersect-circle.js
function intersectCircle(node, rx, point) {
  return intersect_ellipse_default(node, rx, rx, point);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(intersectCircle, "intersectCircle");
var intersect_circle_default = intersectCircle;

// src/rendering-util/rendering-elements/intersect/intersect-line.js
function intersectLine(p1, p2, q1, q2) {
  {
    const a1 = p2.y - p1.y;
    const b1 = p1.x - p2.x;
    const c1 = p2.x * p1.y - p1.x * p2.y;
    const r3 = a1 * q1.x + b1 * q1.y + c1;
    const r4 = a1 * q2.x + b1 * q2.y + c1;
    const epsilon = 1e-6;
    if (r3 !== 0 && r4 !== 0 && sameSign(r3, r4)) {
      return;
    }
    const a2 = q2.y - q1.y;
    const b2 = q1.x - q2.x;
    const c2 = q2.x * q1.y - q1.x * q2.y;
    const r1 = a2 * p1.x + b2 * p1.y + c2;
    const r2 = a2 * p2.x + b2 * p2.y + c2;
    if (Math.abs(r1) < epsilon && Math.abs(r2) < epsilon && sameSign(r1, r2)) {
      return;
    }
    const denom = a1 * b2 - a2 * b1;
    if (denom === 0) {
      return;
    }
    const offset = Math.abs(denom / 2);
    let num = b1 * c2 - b2 * c1;
    const x = num < 0 ? (num - offset) / denom : (num + offset) / denom;
    num = a2 * c1 - a1 * c2;
    const y = num < 0 ? (num - offset) / denom : (num + offset) / denom;
    return { x, y };
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(intersectLine, "intersectLine");
function sameSign(r1, r2) {
  return r1 * r2 > 0;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(sameSign, "sameSign");
var intersect_line_default = intersectLine;

// src/rendering-util/rendering-elements/intersect/intersect-polygon.js
function intersectPolygon(node, polyPoints, point) {
  let x1 = node.x;
  let y1 = node.y;
  let intersections = [];
  let minX = Number.POSITIVE_INFINITY;
  let minY = Number.POSITIVE_INFINITY;
  if (typeof polyPoints.forEach === "function") {
    polyPoints.forEach(function(entry) {
      minX = Math.min(minX, entry.x);
      minY = Math.min(minY, entry.y);
    });
  } else {
    minX = Math.min(minX, polyPoints.x);
    minY = Math.min(minY, polyPoints.y);
  }
  let left = x1 - node.width / 2 - minX;
  let top = y1 - node.height / 2 - minY;
  for (let i = 0; i < polyPoints.length; i++) {
    let p1 = polyPoints[i];
    let p2 = polyPoints[i < polyPoints.length - 1 ? i + 1 : 0];
    let intersect = intersect_line_default(
      node,
      point,
      { x: left + p1.x, y: top + p1.y },
      { x: left + p2.x, y: top + p2.y }
    );
    if (intersect) {
      intersections.push(intersect);
    }
  }
  if (!intersections.length) {
    return node;
  }
  if (intersections.length > 1) {
    intersections.sort(function(p, q) {
      let pdx = p.x - point.x;
      let pdy = p.y - point.y;
      let distp = Math.sqrt(pdx * pdx + pdy * pdy);
      let qdx = q.x - point.x;
      let qdy = q.y - point.y;
      let distq = Math.sqrt(qdx * qdx + qdy * qdy);
      return distp < distq ? -1 : distp === distq ? 0 : 1;
    });
  }
  return intersections[0];
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(intersectPolygon, "intersectPolygon");
var intersect_polygon_default = intersectPolygon;

// src/rendering-util/rendering-elements/intersect/index.js
var intersect_default = {
  node: intersect_node_default,
  circle: intersect_circle_default,
  ellipse: intersect_ellipse_default,
  polygon: intersect_polygon_default,
  rect: intersect_rect_default
};

// src/rendering-util/rendering-elements/shapes/anchor.ts

function anchor(parent, node) {
  const { labelStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const classes = getNodeClasses(node);
  let cssClasses = classes;
  if (!classes) {
    cssClasses = "anchor";
  }
  const shapeSvg = parent.insert("g").attr("class", cssClasses).attr("id", node.domId || node.id);
  const radius = 1;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, { fill: "black", stroke: "none", fillStyle: "solid" });
  if (node.look !== "handDrawn") {
    options.roughness = 0;
  }
  const roughNode = rc.circle(0, 0, radius * 2, options);
  const circleElem = shapeSvg.insert(() => roughNode, ":first-child");
  circleElem.attr("class", "anchor").attr("style", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(cssStyles));
  updateNodeBounds(node, circleElem);
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("Circle intersect", node, radius, point);
    return intersect_default.circle(node, radius, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(anchor, "anchor");

// src/rendering-util/rendering-elements/shapes/bowTieRect.ts

function generateArcPoints(x1, y1, x2, y2, rx, ry, clockwise) {
  const numPoints = 20;
  const midX = (x1 + x2) / 2;
  const midY = (y1 + y2) / 2;
  const angle = Math.atan2(y2 - y1, x2 - x1);
  const dx = (x2 - x1) / 2;
  const dy = (y2 - y1) / 2;
  const transformedX = dx / rx;
  const transformedY = dy / ry;
  const distance = Math.sqrt(transformedX ** 2 + transformedY ** 2);
  if (distance > 1) {
    throw new Error("The given radii are too small to create an arc between the points.");
  }
  const scaledCenterDistance = Math.sqrt(1 - distance ** 2);
  const centerX = midX + scaledCenterDistance * ry * Math.sin(angle) * (clockwise ? -1 : 1);
  const centerY = midY - scaledCenterDistance * rx * Math.cos(angle) * (clockwise ? -1 : 1);
  const startAngle = Math.atan2((y1 - centerY) / ry, (x1 - centerX) / rx);
  const endAngle = Math.atan2((y2 - centerY) / ry, (x2 - centerX) / rx);
  let angleRange = endAngle - startAngle;
  if (clockwise && angleRange < 0) {
    angleRange += 2 * Math.PI;
  }
  if (!clockwise && angleRange > 0) {
    angleRange -= 2 * Math.PI;
  }
  const points = [];
  for (let i = 0; i < numPoints; i++) {
    const t = i / (numPoints - 1);
    const angle2 = startAngle + t * angleRange;
    const x = centerX + rx * Math.cos(angle2);
    const y = centerY + ry * Math.sin(angle2);
    points.push({ x, y });
  }
  return points;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(generateArcPoints, "generateArcPoints");
async function bowTieRect(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const w = bbox.width + node.padding + 20;
  const h = bbox.height + node.padding;
  const ry = h / 2;
  const rx = ry / (2.5 + h / 50);
  const { cssStyles } = node;
  const points = [
    { x: w / 2, y: -h / 2 },
    { x: -w / 2, y: -h / 2 },
    ...generateArcPoints(-w / 2, -h / 2, -w / 2, h / 2, rx, ry, false),
    { x: w / 2, y: h / 2 },
    ...generateArcPoints(w / 2, h / 2, w / 2, -h / 2, rx, ry, true)
  ];
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const bowTieRectPath = createPathFromPoints(points);
  const bowTieRectShapePath = rc.path(bowTieRectPath, options);
  const bowTieRectShape = shapeSvg.insert(() => bowTieRectShapePath, ":first-child");
  bowTieRectShape.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    bowTieRectShape.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    bowTieRectShape.selectAll("path").attr("style", nodeStyles);
  }
  bowTieRectShape.attr("transform", `translate(${rx / 2}, 0)`);
  updateNodeBounds(node, bowTieRectShape);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(bowTieRect, "bowTieRect");

// src/rendering-util/rendering-elements/shapes/card.ts


// src/rendering-util/rendering-elements/shapes/insertPolygonShape.ts
function insertPolygonShape(parent, w, h, points) {
  return parent.insert("polygon", ":first-child").attr(
    "points",
    points.map(function(d) {
      return d.x + "," + d.y;
    }).join(" ")
  ).attr("class", "label-container").attr("transform", "translate(" + -w / 2 + "," + h / 2 + ")");
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(insertPolygonShape, "insertPolygonShape");

// src/rendering-util/rendering-elements/shapes/card.ts
async function card(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const h = bbox.height + node.padding;
  const padding = 12;
  const w = bbox.width + node.padding + padding;
  const left = 0;
  const right = w;
  const top = -h;
  const bottom = 0;
  const points = [
    { x: left + padding, y: top },
    { x: right, y: top },
    { x: right, y: bottom },
    { x: left, y: bottom },
    { x: left, y: top + padding },
    { x: left + padding, y: top }
  ];
  let polygon;
  const { cssStyles } = node;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
    const pathData = createPathFromPoints(points);
    const roughNode = rc.path(pathData, options);
    polygon = shapeSvg.insert(() => roughNode, ":first-child").attr("transform", `translate(${-w / 2}, ${h / 2})`);
    if (cssStyles) {
      polygon.attr("style", cssStyles);
    }
  } else {
    polygon = insertPolygonShape(shapeSvg, w, h, points);
  }
  if (nodeStyles) {
    polygon.attr("style", nodeStyles);
  }
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    return intersect_default.polygon(node, points, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(card, "card");

// src/rendering-util/rendering-elements/shapes/choice.ts

function choice(parent, node) {
  const { nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.label = "";
  const shapeSvg = parent.insert("g").attr("class", getNodeClasses(node)).attr("id", node.domId ?? node.id);
  const { cssStyles } = node;
  const s = Math.max(28, node.width ?? 0);
  const points = [
    { x: 0, y: s / 2 },
    { x: s / 2, y: 0 },
    { x: 0, y: -s / 2 },
    { x: -s / 2, y: 0 }
  ];
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const choicePath = createPathFromPoints(points);
  const roughNode = rc.path(choicePath, options);
  const choiceShape = shapeSvg.insert(() => roughNode, ":first-child");
  if (cssStyles && node.look !== "handDrawn") {
    choiceShape.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    choiceShape.selectAll("path").attr("style", nodeStyles);
  }
  node.width = 28;
  node.height = 28;
  node.intersect = function(point) {
    return intersect_default.polygon(node, points, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(choice, "choice");

// src/rendering-util/rendering-elements/shapes/circle.ts

async function circle(parent, node, options) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, halfPadding } = await labelHelper(parent, node, getNodeClasses(node));
  const padding = options?.padding ?? halfPadding;
  const radius = bbox.width / 2 + padding;
  let circleElem;
  const { cssStyles } = node;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options2 = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
    const roughNode = rc.circle(0, 0, radius * 2, options2);
    circleElem = shapeSvg.insert(() => roughNode, ":first-child");
    circleElem.attr("class", "basic label-container").attr("style", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(cssStyles));
  } else {
    circleElem = shapeSvg.insert("circle", ":first-child").attr("class", "basic label-container").attr("style", nodeStyles).attr("r", radius).attr("cx", 0).attr("cy", 0);
  }
  updateNodeBounds(node, circleElem);
  node.calcIntersect = function(bounds, point) {
    const radius2 = bounds.width / 2;
    return intersect_default.circle(bounds, radius2, point);
  };
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("Circle intersect", node, radius, point);
    return intersect_default.circle(node, radius, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(circle, "circle");

// src/rendering-util/rendering-elements/shapes/crossedCircle.ts

function createLine(r) {
  const xAxis45 = Math.cos(Math.PI / 4);
  const yAxis45 = Math.sin(Math.PI / 4);
  const lineLength = r * 2;
  const pointQ1 = { x: lineLength / 2 * xAxis45, y: lineLength / 2 * yAxis45 };
  const pointQ2 = { x: -(lineLength / 2) * xAxis45, y: lineLength / 2 * yAxis45 };
  const pointQ3 = { x: -(lineLength / 2) * xAxis45, y: -(lineLength / 2) * yAxis45 };
  const pointQ4 = { x: lineLength / 2 * xAxis45, y: -(lineLength / 2) * yAxis45 };
  return `M ${pointQ2.x},${pointQ2.y} L ${pointQ4.x},${pointQ4.y}
                   M ${pointQ1.x},${pointQ1.y} L ${pointQ3.x},${pointQ3.y}`;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(createLine, "createLine");
function crossedCircle(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  node.label = "";
  const shapeSvg = parent.insert("g").attr("class", getNodeClasses(node)).attr("id", node.domId ?? node.id);
  const radius = Math.max(30, node?.width ?? 0);
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const circleNode = rc.circle(0, 0, radius * 2, options);
  const linePath = createLine(radius);
  const lineNode = rc.path(linePath, options);
  const crossedCircle2 = shapeSvg.insert(() => circleNode, ":first-child");
  crossedCircle2.insert(() => lineNode);
  if (cssStyles && node.look !== "handDrawn") {
    crossedCircle2.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    crossedCircle2.selectAll("path").attr("style", nodeStyles);
  }
  updateNodeBounds(node, crossedCircle2);
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("crossedCircle intersect", node, { radius, point });
    const pos = intersect_default.circle(node, radius, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(crossedCircle, "crossedCircle");

// src/rendering-util/rendering-elements/shapes/curlyBraceLeft.ts

function generateCirclePoints2(centerX, centerY, radius, numPoints = 100, startAngle = 0, endAngle = 180) {
  const points = [];
  const startAngleRad = startAngle * Math.PI / 180;
  const endAngleRad = endAngle * Math.PI / 180;
  const angleRange = endAngleRad - startAngleRad;
  const angleStep = angleRange / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    const angle = startAngleRad + i * angleStep;
    const x = centerX + radius * Math.cos(angle);
    const y = centerY + radius * Math.sin(angle);
    points.push({ x: -x, y: -y });
  }
  return points;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(generateCirclePoints2, "generateCirclePoints");
async function curlyBraceLeft(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = bbox.width + (node.padding ?? 0);
  const h = bbox.height + (node.padding ?? 0);
  const radius = Math.max(5, h * 0.1);
  const { cssStyles } = node;
  const points = [
    ...generateCirclePoints2(w / 2, -h / 2, radius, 30, -90, 0),
    { x: -w / 2 - radius, y: radius },
    ...generateCirclePoints2(w / 2 + radius * 2, -radius, radius, 20, -180, -270),
    ...generateCirclePoints2(w / 2 + radius * 2, radius, radius, 20, -90, -180),
    { x: -w / 2 - radius, y: -h / 2 },
    ...generateCirclePoints2(w / 2, h / 2, radius, 20, 0, 90)
  ];
  const rectPoints = [
    { x: w / 2, y: -h / 2 - radius },
    { x: -w / 2, y: -h / 2 - radius },
    ...generateCirclePoints2(w / 2, -h / 2, radius, 20, -90, 0),
    { x: -w / 2 - radius, y: -radius },
    ...generateCirclePoints2(w / 2 + w * 0.1, -radius, radius, 20, -180, -270),
    ...generateCirclePoints2(w / 2 + w * 0.1, radius, radius, 20, -90, -180),
    { x: -w / 2 - radius, y: h / 2 },
    ...generateCirclePoints2(w / 2, h / 2, radius, 20, 0, 90),
    { x: -w / 2, y: h / 2 + radius },
    { x: w / 2, y: h / 2 + radius }
  ];
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, { fill: "none" });
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const curlyBraceLeftPath = createPathFromPoints(points);
  const newCurlyBracePath = curlyBraceLeftPath.replace("Z", "");
  const curlyBraceLeftNode = rc.path(newCurlyBracePath, options);
  const rectPath = createPathFromPoints(rectPoints);
  const rectShape = rc.path(rectPath, { ...options });
  const curlyBraceLeftShape = shapeSvg.insert("g", ":first-child");
  curlyBraceLeftShape.insert(() => rectShape, ":first-child").attr("stroke-opacity", 0);
  curlyBraceLeftShape.insert(() => curlyBraceLeftNode, ":first-child");
  curlyBraceLeftShape.attr("class", "text");
  if (cssStyles && node.look !== "handDrawn") {
    curlyBraceLeftShape.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    curlyBraceLeftShape.selectAll("path").attr("style", nodeStyles);
  }
  curlyBraceLeftShape.attr("transform", `translate(${radius}, 0)`);
  label.attr(
    "transform",
    `translate(${-w / 2 + radius - (bbox.x - (bbox.left ?? 0))},${-h / 2 + (node.padding ?? 0) / 2 - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, curlyBraceLeftShape);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, rectPoints, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(curlyBraceLeft, "curlyBraceLeft");

// src/rendering-util/rendering-elements/shapes/curlyBraceRight.ts

function generateCirclePoints3(centerX, centerY, radius, numPoints = 100, startAngle = 0, endAngle = 180) {
  const points = [];
  const startAngleRad = startAngle * Math.PI / 180;
  const endAngleRad = endAngle * Math.PI / 180;
  const angleRange = endAngleRad - startAngleRad;
  const angleStep = angleRange / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    const angle = startAngleRad + i * angleStep;
    const x = centerX + radius * Math.cos(angle);
    const y = centerY + radius * Math.sin(angle);
    points.push({ x, y });
  }
  return points;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(generateCirclePoints3, "generateCirclePoints");
async function curlyBraceRight(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = bbox.width + (node.padding ?? 0);
  const h = bbox.height + (node.padding ?? 0);
  const radius = Math.max(5, h * 0.1);
  const { cssStyles } = node;
  const points = [
    ...generateCirclePoints3(w / 2, -h / 2, radius, 20, -90, 0),
    { x: w / 2 + radius, y: -radius },
    ...generateCirclePoints3(w / 2 + radius * 2, -radius, radius, 20, -180, -270),
    ...generateCirclePoints3(w / 2 + radius * 2, radius, radius, 20, -90, -180),
    { x: w / 2 + radius, y: h / 2 },
    ...generateCirclePoints3(w / 2, h / 2, radius, 20, 0, 90)
  ];
  const rectPoints = [
    { x: -w / 2, y: -h / 2 - radius },
    { x: w / 2, y: -h / 2 - radius },
    ...generateCirclePoints3(w / 2, -h / 2, radius, 20, -90, 0),
    { x: w / 2 + radius, y: -radius },
    ...generateCirclePoints3(w / 2 + radius * 2, -radius, radius, 20, -180, -270),
    ...generateCirclePoints3(w / 2 + radius * 2, radius, radius, 20, -90, -180),
    { x: w / 2 + radius, y: h / 2 },
    ...generateCirclePoints3(w / 2, h / 2, radius, 20, 0, 90),
    { x: w / 2, y: h / 2 + radius },
    { x: -w / 2, y: h / 2 + radius }
  ];
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, { fill: "none" });
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const curlyBraceRightPath = createPathFromPoints(points);
  const newCurlyBracePath = curlyBraceRightPath.replace("Z", "");
  const curlyBraceRightNode = rc.path(newCurlyBracePath, options);
  const rectPath = createPathFromPoints(rectPoints);
  const rectShape = rc.path(rectPath, { ...options });
  const curlyBraceRightShape = shapeSvg.insert("g", ":first-child");
  curlyBraceRightShape.insert(() => rectShape, ":first-child").attr("stroke-opacity", 0);
  curlyBraceRightShape.insert(() => curlyBraceRightNode, ":first-child");
  curlyBraceRightShape.attr("class", "text");
  if (cssStyles && node.look !== "handDrawn") {
    curlyBraceRightShape.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    curlyBraceRightShape.selectAll("path").attr("style", nodeStyles);
  }
  curlyBraceRightShape.attr("transform", `translate(${-radius}, 0)`);
  label.attr(
    "transform",
    `translate(${-w / 2 + (node.padding ?? 0) / 2 - (bbox.x - (bbox.left ?? 0))},${-h / 2 + (node.padding ?? 0) / 2 - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, curlyBraceRightShape);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, rectPoints, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(curlyBraceRight, "curlyBraceRight");

// src/rendering-util/rendering-elements/shapes/curlyBraces.ts

function generateCirclePoints4(centerX, centerY, radius, numPoints = 100, startAngle = 0, endAngle = 180) {
  const points = [];
  const startAngleRad = startAngle * Math.PI / 180;
  const endAngleRad = endAngle * Math.PI / 180;
  const angleRange = endAngleRad - startAngleRad;
  const angleStep = angleRange / (numPoints - 1);
  for (let i = 0; i < numPoints; i++) {
    const angle = startAngleRad + i * angleStep;
    const x = centerX + radius * Math.cos(angle);
    const y = centerY + radius * Math.sin(angle);
    points.push({ x: -x, y: -y });
  }
  return points;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(generateCirclePoints4, "generateCirclePoints");
async function curlyBraces(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = bbox.width + (node.padding ?? 0);
  const h = bbox.height + (node.padding ?? 0);
  const radius = Math.max(5, h * 0.1);
  const { cssStyles } = node;
  const leftCurlyBracePoints = [
    ...generateCirclePoints4(w / 2, -h / 2, radius, 30, -90, 0),
    { x: -w / 2 - radius, y: radius },
    ...generateCirclePoints4(w / 2 + radius * 2, -radius, radius, 20, -180, -270),
    ...generateCirclePoints4(w / 2 + radius * 2, radius, radius, 20, -90, -180),
    { x: -w / 2 - radius, y: -h / 2 },
    ...generateCirclePoints4(w / 2, h / 2, radius, 20, 0, 90)
  ];
  const rightCurlyBracePoints = [
    ...generateCirclePoints4(-w / 2 + radius + radius / 2, -h / 2, radius, 20, -90, -180),
    { x: w / 2 - radius / 2, y: radius },
    ...generateCirclePoints4(-w / 2 - radius / 2, -radius, radius, 20, 0, 90),
    ...generateCirclePoints4(-w / 2 - radius / 2, radius, radius, 20, -90, 0),
    { x: w / 2 - radius / 2, y: -radius },
    ...generateCirclePoints4(-w / 2 + radius + radius / 2, h / 2, radius, 30, -180, -270)
  ];
  const rectPoints = [
    { x: w / 2, y: -h / 2 - radius },
    { x: -w / 2, y: -h / 2 - radius },
    ...generateCirclePoints4(w / 2, -h / 2, radius, 20, -90, 0),
    { x: -w / 2 - radius, y: -radius },
    ...generateCirclePoints4(w / 2 + radius * 2, -radius, radius, 20, -180, -270),
    ...generateCirclePoints4(w / 2 + radius * 2, radius, radius, 20, -90, -180),
    { x: -w / 2 - radius, y: h / 2 },
    ...generateCirclePoints4(w / 2, h / 2, radius, 20, 0, 90),
    { x: -w / 2, y: h / 2 + radius },
    { x: w / 2 - radius - radius / 2, y: h / 2 + radius },
    ...generateCirclePoints4(-w / 2 + radius + radius / 2, -h / 2, radius, 20, -90, -180),
    { x: w / 2 - radius / 2, y: radius },
    ...generateCirclePoints4(-w / 2 - radius / 2, -radius, radius, 20, 0, 90),
    ...generateCirclePoints4(-w / 2 - radius / 2, radius, radius, 20, -90, 0),
    { x: w / 2 - radius / 2, y: -radius },
    ...generateCirclePoints4(-w / 2 + radius + radius / 2, h / 2, radius, 30, -180, -270)
  ];
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, { fill: "none" });
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const leftCurlyBracePath = createPathFromPoints(leftCurlyBracePoints);
  const newLeftCurlyBracePath = leftCurlyBracePath.replace("Z", "");
  const leftCurlyBraceNode = rc.path(newLeftCurlyBracePath, options);
  const rightCurlyBracePath = createPathFromPoints(rightCurlyBracePoints);
  const newRightCurlyBracePath = rightCurlyBracePath.replace("Z", "");
  const rightCurlyBraceNode = rc.path(newRightCurlyBracePath, options);
  const rectPath = createPathFromPoints(rectPoints);
  const rectShape = rc.path(rectPath, { ...options });
  const curlyBracesShape = shapeSvg.insert("g", ":first-child");
  curlyBracesShape.insert(() => rectShape, ":first-child").attr("stroke-opacity", 0);
  curlyBracesShape.insert(() => leftCurlyBraceNode, ":first-child");
  curlyBracesShape.insert(() => rightCurlyBraceNode, ":first-child");
  curlyBracesShape.attr("class", "text");
  if (cssStyles && node.look !== "handDrawn") {
    curlyBracesShape.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    curlyBracesShape.selectAll("path").attr("style", nodeStyles);
  }
  curlyBracesShape.attr("transform", `translate(${radius - radius / 4}, 0)`);
  label.attr(
    "transform",
    `translate(${-w / 2 + (node.padding ?? 0) / 2 - (bbox.x - (bbox.left ?? 0))},${-h / 2 + (node.padding ?? 0) / 2 - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, curlyBracesShape);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, rectPoints, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(curlyBraces, "curlyBraces");

// src/rendering-util/rendering-elements/shapes/curvedTrapezoid.ts

async function curvedTrapezoid(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const minWidth = 80, minHeight = 20;
  const w = Math.max(minWidth, (bbox.width + (node.padding ?? 0) * 2) * 1.25, node?.width ?? 0);
  const h = Math.max(minHeight, bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const radius = h / 2;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const totalWidth = w, totalHeight = h;
  const rw = totalWidth - radius;
  const tw = totalHeight / 4;
  const points = [
    { x: rw, y: 0 },
    { x: tw, y: 0 },
    { x: 0, y: totalHeight / 2 },
    { x: tw, y: totalHeight },
    { x: rw, y: totalHeight },
    ...generateCirclePoints(-rw, -totalHeight / 2, radius, 50, 270, 90)
  ];
  const pathData = createPathFromPoints(points);
  const shapeNode = rc.path(pathData, options);
  const polygon = shapeSvg.insert(() => shapeNode, ":first-child");
  polygon.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", nodeStyles);
  }
  polygon.attr("transform", `translate(${-w / 2}, ${-h / 2})`);
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(curvedTrapezoid, "curvedTrapezoid");

// src/rendering-util/rendering-elements/shapes/cylinder.ts

var createCylinderPathD = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((x, y, width, height, rx, ry) => {
  return [
    `M${x},${y + ry}`,
    `a${rx},${ry} 0,0,0 ${width},0`,
    `a${rx},${ry} 0,0,0 ${-width},0`,
    `l0,${height}`,
    `a${rx},${ry} 0,0,0 ${width},0`,
    `l0,${-height}`
  ].join(" ");
}, "createCylinderPathD");
var createOuterCylinderPathD = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((x, y, width, height, rx, ry) => {
  return [
    `M${x},${y + ry}`,
    `M${x + width},${y + ry}`,
    `a${rx},${ry} 0,0,0 ${-width},0`,
    `l0,${height}`,
    `a${rx},${ry} 0,0,0 ${width},0`,
    `l0,${-height}`
  ].join(" ");
}, "createOuterCylinderPathD");
var createInnerCylinderPathD = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((x, y, width, height, rx, ry) => {
  return [`M${x - width / 2},${-height / 2}`, `a${rx},${ry} 0,0,0 ${width},0`].join(" ");
}, "createInnerCylinderPathD");
async function cylinder(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + node.padding, node.width ?? 0);
  const rx = w / 2;
  const ry = rx / (2.5 + w / 50);
  const h = Math.max(bbox.height + ry + node.padding, node.height ?? 0);
  let cylinder2;
  const { cssStyles } = node;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const outerPathData = createOuterCylinderPathD(0, 0, w, h, rx, ry);
    const innerPathData = createInnerCylinderPathD(0, ry, w, h, rx, ry);
    const outerNode = rc.path(outerPathData, (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {}));
    const innerLine = rc.path(innerPathData, (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, { fill: "none" }));
    cylinder2 = shapeSvg.insert(() => innerLine, ":first-child");
    cylinder2 = shapeSvg.insert(() => outerNode, ":first-child");
    cylinder2.attr("class", "basic label-container");
    if (cssStyles) {
      cylinder2.attr("style", cssStyles);
    }
  } else {
    const pathData = createCylinderPathD(0, 0, w, h, rx, ry);
    cylinder2 = shapeSvg.insert("path", ":first-child").attr("d", pathData).attr("class", "basic label-container").attr("style", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(cssStyles)).attr("style", nodeStyles);
  }
  cylinder2.attr("label-offset-y", ry);
  cylinder2.attr("transform", `translate(${-w / 2}, ${-(h / 2 + ry)})`);
  updateNodeBounds(node, cylinder2);
  label.attr(
    "transform",
    `translate(${-(bbox.width / 2) - (bbox.x - (bbox.left ?? 0))}, ${-(bbox.height / 2) + (node.padding ?? 0) / 1.5 - (bbox.y - (bbox.top ?? 0))})`
  );
  node.intersect = function(point) {
    const pos = intersect_default.rect(node, point);
    const x = pos.x - (node.x ?? 0);
    if (rx != 0 && (Math.abs(x) < (node.width ?? 0) / 2 || Math.abs(x) == (node.width ?? 0) / 2 && Math.abs(pos.y - (node.y ?? 0)) > (node.height ?? 0) / 2 - ry)) {
      let y = ry * ry * (1 - x * x / (rx * rx));
      if (y > 0) {
        y = Math.sqrt(y);
      }
      y = ry - y;
      if (point.y - (node.y ?? 0) > 0) {
        y = -y;
      }
      pos.y += y;
    }
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(cylinder, "cylinder");

// src/rendering-util/rendering-elements/shapes/dividedRect.ts

async function dividedRectangle(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = bbox.width + node.padding;
  const h = bbox.height + node.padding;
  const rectOffset = h * 0.2;
  const x = -w / 2;
  const y = -h / 2 - rectOffset / 2;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const pts = [
    { x, y: y + rectOffset },
    { x: -x, y: y + rectOffset },
    { x: -x, y: -y },
    { x, y: -y },
    { x, y },
    { x: -x, y },
    { x: -x, y: y + rectOffset }
  ];
  const poly = rc.polygon(
    pts.map((p) => [p.x, p.y]),
    options
  );
  const polygon = shapeSvg.insert(() => poly, ":first-child");
  polygon.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    polygon.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    polygon.selectAll("path").attr("style", nodeStyles);
  }
  label.attr(
    "transform",
    `translate(${x + (node.padding ?? 0) / 2 - (bbox.x - (bbox.left ?? 0))}, ${y + rectOffset + (node.padding ?? 0) / 2 - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    const pos = intersect_default.rect(node, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(dividedRectangle, "dividedRectangle");

// src/rendering-util/rendering-elements/shapes/doubleCircle.ts

async function doublecircle(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, halfPadding } = await labelHelper(parent, node, getNodeClasses(node));
  const gap = 5;
  const outerRadius = bbox.width / 2 + halfPadding + gap;
  const innerRadius = bbox.width / 2 + halfPadding;
  let circleGroup;
  const { cssStyles } = node;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const outerOptions = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, { roughness: 0.2, strokeWidth: 2.5 });
    const innerOptions = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, { roughness: 0.2, strokeWidth: 1.5 });
    const outerRoughNode = rc.circle(0, 0, outerRadius * 2, outerOptions);
    const innerRoughNode = rc.circle(0, 0, innerRadius * 2, innerOptions);
    circleGroup = shapeSvg.insert("g", ":first-child");
    circleGroup.attr("class", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(node.cssClasses)).attr("style", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(cssStyles));
    circleGroup.node()?.appendChild(outerRoughNode);
    circleGroup.node()?.appendChild(innerRoughNode);
  } else {
    circleGroup = shapeSvg.insert("g", ":first-child");
    const outerCircle = circleGroup.insert("circle", ":first-child");
    const innerCircle = circleGroup.insert("circle");
    circleGroup.attr("class", "basic label-container").attr("style", nodeStyles);
    outerCircle.attr("class", "outer-circle").attr("style", nodeStyles).attr("r", outerRadius).attr("cx", 0).attr("cy", 0);
    innerCircle.attr("class", "inner-circle").attr("style", nodeStyles).attr("r", innerRadius).attr("cx", 0).attr("cy", 0);
  }
  updateNodeBounds(node, circleGroup);
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("DoubleCircle intersect", node, outerRadius, point);
    return intersect_default.circle(node, outerRadius, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(doublecircle, "doublecircle");

// src/rendering-util/rendering-elements/shapes/filledCircle.ts

function filledCircle(parent, node, { config: { themeVariables } }) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.label = "";
  node.labelStyle = labelStyles;
  const shapeSvg = parent.insert("g").attr("class", getNodeClasses(node)).attr("id", node.domId ?? node.id);
  const radius = 7;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const { nodeBorder } = themeVariables;
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, { fillStyle: "solid" });
  if (node.look !== "handDrawn") {
    options.roughness = 0;
  }
  const circleNode = rc.circle(0, 0, radius * 2, options);
  const filledCircle2 = shapeSvg.insert(() => circleNode, ":first-child");
  filledCircle2.selectAll("path").attr("style", `fill: ${nodeBorder} !important;`);
  if (cssStyles && cssStyles.length > 0 && node.look !== "handDrawn") {
    filledCircle2.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    filledCircle2.selectAll("path").attr("style", nodeStyles);
  }
  updateNodeBounds(node, filledCircle2);
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("filledCircle intersect", node, { radius, point });
    const pos = intersect_default.circle(node, radius, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(filledCircle, "filledCircle");

// src/rendering-util/rendering-elements/shapes/flippedTriangle.ts

async function flippedTriangle(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = bbox.width + (node.padding ?? 0);
  const h = w + bbox.height;
  const tw = w + bbox.height;
  const points = [
    { x: 0, y: -h },
    { x: tw, y: -h },
    { x: tw / 2, y: 0 }
  ];
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const pathData = createPathFromPoints(points);
  const roughNode = rc.path(pathData, options);
  const flippedTriangle2 = shapeSvg.insert(() => roughNode, ":first-child").attr("transform", `translate(${-h / 2}, ${h / 2})`);
  if (cssStyles && node.look !== "handDrawn") {
    flippedTriangle2.selectChildren("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    flippedTriangle2.selectChildren("path").attr("style", nodeStyles);
  }
  node.width = w;
  node.height = h;
  updateNodeBounds(node, flippedTriangle2);
  label.attr(
    "transform",
    `translate(${-bbox.width / 2 - (bbox.x - (bbox.left ?? 0))}, ${-h / 2 + (node.padding ?? 0) / 2 + (bbox.y - (bbox.top ?? 0))})`
  );
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("Triangle intersect", node, points, point);
    return intersect_default.polygon(node, points, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(flippedTriangle, "flippedTriangle");

// src/rendering-util/rendering-elements/shapes/forkJoin.ts

function forkJoin(parent, node, { dir, config: { state: state2, themeVariables } }) {
  const { nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.label = "";
  const shapeSvg = parent.insert("g").attr("class", getNodeClasses(node)).attr("id", node.domId ?? node.id);
  const { cssStyles } = node;
  let width = Math.max(70, node?.width ?? 0);
  let height = Math.max(10, node?.height ?? 0);
  if (dir === "LR") {
    width = Math.max(10, node?.width ?? 0);
    height = Math.max(70, node?.height ?? 0);
  }
  const x = -1 * width / 2;
  const y = -1 * height / 2;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {
    stroke: themeVariables.lineColor,
    fill: themeVariables.lineColor
  });
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const roughNode = rc.rectangle(x, y, width, height, options);
  const shape = shapeSvg.insert(() => roughNode, ":first-child");
  if (cssStyles && node.look !== "handDrawn") {
    shape.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    shape.selectAll("path").attr("style", nodeStyles);
  }
  updateNodeBounds(node, shape);
  const padding = state2?.padding ?? 0;
  if (node.width && node.height) {
    node.width += padding / 2 || 0;
    node.height += padding / 2 || 0;
  }
  node.intersect = function(point) {
    return intersect_default.rect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(forkJoin, "forkJoin");

// src/rendering-util/rendering-elements/shapes/halfRoundedRectangle.ts

async function halfRoundedRectangle(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const minWidth = 80, minHeight = 50;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(minWidth, bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const h = Math.max(minHeight, bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const radius = h / 2;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const points = [
    { x: -w / 2, y: -h / 2 },
    { x: w / 2 - radius, y: -h / 2 },
    ...generateCirclePoints(-w / 2 + radius, 0, radius, 50, 90, 270),
    { x: w / 2 - radius, y: h / 2 },
    { x: -w / 2, y: h / 2 }
  ];
  const pathData = createPathFromPoints(points);
  const shapeNode = rc.path(pathData, options);
  const polygon = shapeSvg.insert(() => shapeNode, ":first-child");
  polygon.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", nodeStyles);
  }
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("Pill intersect", node, { radius, point });
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(halfRoundedRectangle, "halfRoundedRectangle");

// src/rendering-util/rendering-elements/shapes/hexagon.ts

async function hexagon(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const h = bbox.height + (node.padding ?? 0);
  const w = bbox.width + (node.padding ?? 0) * 2.5;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  let halfWidth = w / 2;
  const m = halfWidth / 6;
  halfWidth = halfWidth + m;
  const halfHeight = h / 2;
  const fixedLength = halfHeight / 2;
  const deducedWidth = halfWidth - fixedLength;
  const points = [
    { x: -deducedWidth, y: -halfHeight },
    { x: 0, y: -halfHeight },
    { x: deducedWidth, y: -halfHeight },
    { x: halfWidth, y: 0 },
    { x: deducedWidth, y: halfHeight },
    { x: 0, y: halfHeight },
    { x: -deducedWidth, y: halfHeight },
    { x: -halfWidth, y: 0 }
  ];
  const pathData = createPathFromPoints(points);
  const shapeNode = rc.path(pathData, options);
  const polygon = shapeSvg.insert(() => shapeNode, ":first-child");
  polygon.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", nodeStyles);
  }
  node.width = w;
  node.height = h;
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    return intersect_default.polygon(node, points, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(hexagon, "hexagon");

// src/rendering-util/rendering-elements/shapes/hourglass.ts

async function hourglass(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.label = "";
  node.labelStyle = labelStyles;
  const { shapeSvg } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(30, node?.width ?? 0);
  const h = Math.max(30, node?.height ?? 0);
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const points = [
    { x: 0, y: 0 },
    { x: w, y: 0 },
    { x: 0, y: h },
    { x: w, y: h }
  ];
  const pathData = createPathFromPoints(points);
  const shapeNode = rc.path(pathData, options);
  const polygon = shapeSvg.insert(() => shapeNode, ":first-child");
  polygon.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", nodeStyles);
  }
  polygon.attr("transform", `translate(${-w / 2}, ${-h / 2})`);
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("Pill intersect", node, { points });
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(hourglass, "hourglass");

// src/rendering-util/rendering-elements/shapes/icon.ts

async function icon(parent, node, { config: { themeVariables, flowchart } }) {
  const { labelStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const assetHeight = node.assetHeight ?? 48;
  const assetWidth = node.assetWidth ?? 48;
  const iconSize = Math.max(assetHeight, assetWidth);
  const defaultWidth = flowchart?.wrappingWidth;
  node.width = Math.max(iconSize, defaultWidth ?? 0);
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, "icon-shape default");
  const topLabel = node.pos === "t";
  const height = iconSize;
  const width = iconSize;
  const { nodeBorder } = themeVariables;
  const { stylesMap } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .compileStyles */ .dT)(node);
  const x = -width / 2;
  const y = -height / 2;
  const labelPadding = node.label ? 8 : 0;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, { stroke: "none", fill: "none" });
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const iconNode = rc.rectangle(x, y, width, height, options);
  const outerWidth = Math.max(width, bbox.width);
  const outerHeight = height + bbox.height + labelPadding;
  const outerNode = rc.rectangle(-outerWidth / 2, -outerHeight / 2, outerWidth, outerHeight, {
    ...options,
    fill: "transparent",
    stroke: "none"
  });
  const iconShape = shapeSvg.insert(() => iconNode, ":first-child");
  const outerShape = shapeSvg.insert(() => outerNode);
  if (node.icon) {
    const iconElem = shapeSvg.append("g");
    iconElem.html(
      `<g>${await (0,_chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getIconSVG */ .s4)(node.icon, {
        height: iconSize,
        width: iconSize,
        fallbackPrefix: ""
      })}</g>`
    );
    const iconBBox = iconElem.node().getBBox();
    const iconWidth = iconBBox.width;
    const iconHeight = iconBBox.height;
    const iconX = iconBBox.x;
    const iconY = iconBBox.y;
    iconElem.attr(
      "transform",
      `translate(${-iconWidth / 2 - iconX},${topLabel ? bbox.height / 2 + labelPadding / 2 - iconHeight / 2 - iconY : -bbox.height / 2 - labelPadding / 2 - iconHeight / 2 - iconY})`
    );
    iconElem.attr("style", `color: ${stylesMap.get("stroke") ?? nodeBorder};`);
  }
  label.attr(
    "transform",
    `translate(${-bbox.width / 2 - (bbox.x - (bbox.left ?? 0))},${topLabel ? -outerHeight / 2 : outerHeight / 2 - bbox.height})`
  );
  iconShape.attr(
    "transform",
    `translate(${0},${topLabel ? bbox.height / 2 + labelPadding / 2 : -bbox.height / 2 - labelPadding / 2})`
  );
  updateNodeBounds(node, outerShape);
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("iconSquare intersect", node, point);
    if (!node.label) {
      return intersect_default.rect(node, point);
    }
    const dx = node.x ?? 0;
    const dy = node.y ?? 0;
    const nodeHeight = node.height ?? 0;
    let points = [];
    if (topLabel) {
      points = [
        { x: dx - bbox.width / 2, y: dy - nodeHeight / 2 },
        { x: dx + bbox.width / 2, y: dy - nodeHeight / 2 },
        { x: dx + bbox.width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding },
        { x: dx + width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding },
        { x: dx + width / 2, y: dy + nodeHeight / 2 },
        { x: dx - width / 2, y: dy + nodeHeight / 2 },
        { x: dx - width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding },
        { x: dx - bbox.width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding }
      ];
    } else {
      points = [
        { x: dx - width / 2, y: dy - nodeHeight / 2 },
        { x: dx + width / 2, y: dy - nodeHeight / 2 },
        { x: dx + width / 2, y: dy - nodeHeight / 2 + height },
        { x: dx + bbox.width / 2, y: dy - nodeHeight / 2 + height },
        { x: dx + bbox.width / 2 / 2, y: dy + nodeHeight / 2 },
        { x: dx - bbox.width / 2, y: dy + nodeHeight / 2 },
        { x: dx - bbox.width / 2, y: dy - nodeHeight / 2 + height },
        { x: dx - width / 2, y: dy - nodeHeight / 2 + height }
      ];
    }
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(icon, "icon");

// src/rendering-util/rendering-elements/shapes/iconCircle.ts

async function iconCircle(parent, node, { config: { themeVariables, flowchart } }) {
  const { labelStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const assetHeight = node.assetHeight ?? 48;
  const assetWidth = node.assetWidth ?? 48;
  const iconSize = Math.max(assetHeight, assetWidth);
  const defaultWidth = flowchart?.wrappingWidth;
  node.width = Math.max(iconSize, defaultWidth ?? 0);
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, "icon-shape default");
  const padding = 20;
  const labelPadding = node.label ? 8 : 0;
  const topLabel = node.pos === "t";
  const { nodeBorder, mainBkg } = themeVariables;
  const { stylesMap } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .compileStyles */ .dT)(node);
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const fill = stylesMap.get("fill");
  options.stroke = fill ?? mainBkg;
  const iconElem = shapeSvg.append("g");
  if (node.icon) {
    iconElem.html(
      `<g>${await (0,_chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getIconSVG */ .s4)(node.icon, {
        height: iconSize,
        width: iconSize,
        fallbackPrefix: ""
      })}</g>`
    );
  }
  const iconBBox = iconElem.node().getBBox();
  const iconWidth = iconBBox.width;
  const iconHeight = iconBBox.height;
  const iconX = iconBBox.x;
  const iconY = iconBBox.y;
  const diameter = Math.max(iconWidth, iconHeight) * Math.SQRT2 + padding * 2;
  const iconNode = rc.circle(0, 0, diameter, options);
  const outerWidth = Math.max(diameter, bbox.width);
  const outerHeight = diameter + bbox.height + labelPadding;
  const outerNode = rc.rectangle(-outerWidth / 2, -outerHeight / 2, outerWidth, outerHeight, {
    ...options,
    fill: "transparent",
    stroke: "none"
  });
  const iconShape = shapeSvg.insert(() => iconNode, ":first-child");
  const outerShape = shapeSvg.insert(() => outerNode);
  iconElem.attr(
    "transform",
    `translate(${-iconWidth / 2 - iconX},${topLabel ? bbox.height / 2 + labelPadding / 2 - iconHeight / 2 - iconY : -bbox.height / 2 - labelPadding / 2 - iconHeight / 2 - iconY})`
  );
  iconElem.attr("style", `color: ${stylesMap.get("stroke") ?? nodeBorder};`);
  label.attr(
    "transform",
    `translate(${-bbox.width / 2 - (bbox.x - (bbox.left ?? 0))},${topLabel ? -outerHeight / 2 : outerHeight / 2 - bbox.height})`
  );
  iconShape.attr(
    "transform",
    `translate(${0},${topLabel ? bbox.height / 2 + labelPadding / 2 : -bbox.height / 2 - labelPadding / 2})`
  );
  updateNodeBounds(node, outerShape);
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("iconSquare intersect", node, point);
    const pos = intersect_default.rect(node, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(iconCircle, "iconCircle");

// src/rendering-util/rendering-elements/shapes/iconRounded.ts

async function iconRounded(parent, node, { config: { themeVariables, flowchart } }) {
  const { labelStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const assetHeight = node.assetHeight ?? 48;
  const assetWidth = node.assetWidth ?? 48;
  const iconSize = Math.max(assetHeight, assetWidth);
  const defaultWidth = flowchart?.wrappingWidth;
  node.width = Math.max(iconSize, defaultWidth ?? 0);
  const { shapeSvg, bbox, halfPadding, label } = await labelHelper(
    parent,
    node,
    "icon-shape default"
  );
  const topLabel = node.pos === "t";
  const height = iconSize + halfPadding * 2;
  const width = iconSize + halfPadding * 2;
  const { nodeBorder, mainBkg } = themeVariables;
  const { stylesMap } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .compileStyles */ .dT)(node);
  const x = -width / 2;
  const y = -height / 2;
  const labelPadding = node.label ? 8 : 0;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const fill = stylesMap.get("fill");
  options.stroke = fill ?? mainBkg;
  const iconNode = rc.path(createRoundedRectPathD(x, y, width, height, 5), options);
  const outerWidth = Math.max(width, bbox.width);
  const outerHeight = height + bbox.height + labelPadding;
  const outerNode = rc.rectangle(-outerWidth / 2, -outerHeight / 2, outerWidth, outerHeight, {
    ...options,
    fill: "transparent",
    stroke: "none"
  });
  const iconShape = shapeSvg.insert(() => iconNode, ":first-child").attr("class", "icon-shape2");
  const outerShape = shapeSvg.insert(() => outerNode);
  if (node.icon) {
    const iconElem = shapeSvg.append("g");
    iconElem.html(
      `<g>${await (0,_chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getIconSVG */ .s4)(node.icon, {
        height: iconSize,
        width: iconSize,
        fallbackPrefix: ""
      })}</g>`
    );
    const iconBBox = iconElem.node().getBBox();
    const iconWidth = iconBBox.width;
    const iconHeight = iconBBox.height;
    const iconX = iconBBox.x;
    const iconY = iconBBox.y;
    iconElem.attr(
      "transform",
      `translate(${-iconWidth / 2 - iconX},${topLabel ? bbox.height / 2 + labelPadding / 2 - iconHeight / 2 - iconY : -bbox.height / 2 - labelPadding / 2 - iconHeight / 2 - iconY})`
    );
    iconElem.attr("style", `color: ${stylesMap.get("stroke") ?? nodeBorder};`);
  }
  label.attr(
    "transform",
    `translate(${-bbox.width / 2 - (bbox.x - (bbox.left ?? 0))},${topLabel ? -outerHeight / 2 : outerHeight / 2 - bbox.height})`
  );
  iconShape.attr(
    "transform",
    `translate(${0},${topLabel ? bbox.height / 2 + labelPadding / 2 : -bbox.height / 2 - labelPadding / 2})`
  );
  updateNodeBounds(node, outerShape);
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("iconSquare intersect", node, point);
    if (!node.label) {
      return intersect_default.rect(node, point);
    }
    const dx = node.x ?? 0;
    const dy = node.y ?? 0;
    const nodeHeight = node.height ?? 0;
    let points = [];
    if (topLabel) {
      points = [
        { x: dx - bbox.width / 2, y: dy - nodeHeight / 2 },
        { x: dx + bbox.width / 2, y: dy - nodeHeight / 2 },
        { x: dx + bbox.width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding },
        { x: dx + width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding },
        { x: dx + width / 2, y: dy + nodeHeight / 2 },
        { x: dx - width / 2, y: dy + nodeHeight / 2 },
        { x: dx - width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding },
        { x: dx - bbox.width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding }
      ];
    } else {
      points = [
        { x: dx - width / 2, y: dy - nodeHeight / 2 },
        { x: dx + width / 2, y: dy - nodeHeight / 2 },
        { x: dx + width / 2, y: dy - nodeHeight / 2 + height },
        { x: dx + bbox.width / 2, y: dy - nodeHeight / 2 + height },
        { x: dx + bbox.width / 2 / 2, y: dy + nodeHeight / 2 },
        { x: dx - bbox.width / 2, y: dy + nodeHeight / 2 },
        { x: dx - bbox.width / 2, y: dy - nodeHeight / 2 + height },
        { x: dx - width / 2, y: dy - nodeHeight / 2 + height }
      ];
    }
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(iconRounded, "iconRounded");

// src/rendering-util/rendering-elements/shapes/iconSquare.ts

async function iconSquare(parent, node, { config: { themeVariables, flowchart } }) {
  const { labelStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const assetHeight = node.assetHeight ?? 48;
  const assetWidth = node.assetWidth ?? 48;
  const iconSize = Math.max(assetHeight, assetWidth);
  const defaultWidth = flowchart?.wrappingWidth;
  node.width = Math.max(iconSize, defaultWidth ?? 0);
  const { shapeSvg, bbox, halfPadding, label } = await labelHelper(
    parent,
    node,
    "icon-shape default"
  );
  const topLabel = node.pos === "t";
  const height = iconSize + halfPadding * 2;
  const width = iconSize + halfPadding * 2;
  const { nodeBorder, mainBkg } = themeVariables;
  const { stylesMap } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .compileStyles */ .dT)(node);
  const x = -width / 2;
  const y = -height / 2;
  const labelPadding = node.label ? 8 : 0;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const fill = stylesMap.get("fill");
  options.stroke = fill ?? mainBkg;
  const iconNode = rc.path(createRoundedRectPathD(x, y, width, height, 0.1), options);
  const outerWidth = Math.max(width, bbox.width);
  const outerHeight = height + bbox.height + labelPadding;
  const outerNode = rc.rectangle(-outerWidth / 2, -outerHeight / 2, outerWidth, outerHeight, {
    ...options,
    fill: "transparent",
    stroke: "none"
  });
  const iconShape = shapeSvg.insert(() => iconNode, ":first-child");
  const outerShape = shapeSvg.insert(() => outerNode);
  if (node.icon) {
    const iconElem = shapeSvg.append("g");
    iconElem.html(
      `<g>${await (0,_chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getIconSVG */ .s4)(node.icon, {
        height: iconSize,
        width: iconSize,
        fallbackPrefix: ""
      })}</g>`
    );
    const iconBBox = iconElem.node().getBBox();
    const iconWidth = iconBBox.width;
    const iconHeight = iconBBox.height;
    const iconX = iconBBox.x;
    const iconY = iconBBox.y;
    iconElem.attr(
      "transform",
      `translate(${-iconWidth / 2 - iconX},${topLabel ? bbox.height / 2 + labelPadding / 2 - iconHeight / 2 - iconY : -bbox.height / 2 - labelPadding / 2 - iconHeight / 2 - iconY})`
    );
    iconElem.attr("style", `color: ${stylesMap.get("stroke") ?? nodeBorder};`);
  }
  label.attr(
    "transform",
    `translate(${-bbox.width / 2 - (bbox.x - (bbox.left ?? 0))},${topLabel ? -outerHeight / 2 : outerHeight / 2 - bbox.height})`
  );
  iconShape.attr(
    "transform",
    `translate(${0},${topLabel ? bbox.height / 2 + labelPadding / 2 : -bbox.height / 2 - labelPadding / 2})`
  );
  updateNodeBounds(node, outerShape);
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("iconSquare intersect", node, point);
    if (!node.label) {
      return intersect_default.rect(node, point);
    }
    const dx = node.x ?? 0;
    const dy = node.y ?? 0;
    const nodeHeight = node.height ?? 0;
    let points = [];
    if (topLabel) {
      points = [
        { x: dx - bbox.width / 2, y: dy - nodeHeight / 2 },
        { x: dx + bbox.width / 2, y: dy - nodeHeight / 2 },
        { x: dx + bbox.width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding },
        { x: dx + width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding },
        { x: dx + width / 2, y: dy + nodeHeight / 2 },
        { x: dx - width / 2, y: dy + nodeHeight / 2 },
        { x: dx - width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding },
        { x: dx - bbox.width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding }
      ];
    } else {
      points = [
        { x: dx - width / 2, y: dy - nodeHeight / 2 },
        { x: dx + width / 2, y: dy - nodeHeight / 2 },
        { x: dx + width / 2, y: dy - nodeHeight / 2 + height },
        { x: dx + bbox.width / 2, y: dy - nodeHeight / 2 + height },
        { x: dx + bbox.width / 2 / 2, y: dy + nodeHeight / 2 },
        { x: dx - bbox.width / 2, y: dy + nodeHeight / 2 },
        { x: dx - bbox.width / 2, y: dy - nodeHeight / 2 + height },
        { x: dx - width / 2, y: dy - nodeHeight / 2 + height }
      ];
    }
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(iconSquare, "iconSquare");

// src/rendering-util/rendering-elements/shapes/imageSquare.ts

async function imageSquare(parent, node, { config: { flowchart } }) {
  const img = new Image();
  img.src = node?.img ?? "";
  await img.decode();
  const imageNaturalWidth = Number(img.naturalWidth.toString().replace("px", ""));
  const imageNaturalHeight = Number(img.naturalHeight.toString().replace("px", ""));
  node.imageAspectRatio = imageNaturalWidth / imageNaturalHeight;
  const { labelStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const defaultWidth = flowchart?.wrappingWidth;
  node.defaultWidth = flowchart?.wrappingWidth;
  const imageRawWidth = Math.max(
    node.label ? defaultWidth ?? 0 : 0,
    node?.assetWidth ?? imageNaturalWidth
  );
  const imageWidth = node.constraint === "on" ? node?.assetHeight ? node.assetHeight * node.imageAspectRatio : imageRawWidth : imageRawWidth;
  const imageHeight = node.constraint === "on" ? imageWidth / node.imageAspectRatio : node?.assetHeight ?? imageNaturalHeight;
  node.width = Math.max(imageWidth, defaultWidth ?? 0);
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, "image-shape default");
  const topLabel = node.pos === "t";
  const x = -imageWidth / 2;
  const y = -imageHeight / 2;
  const labelPadding = node.label ? 8 : 0;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const imageNode = rc.rectangle(x, y, imageWidth, imageHeight, options);
  const outerWidth = Math.max(imageWidth, bbox.width);
  const outerHeight = imageHeight + bbox.height + labelPadding;
  const outerNode = rc.rectangle(-outerWidth / 2, -outerHeight / 2, outerWidth, outerHeight, {
    ...options,
    fill: "none",
    stroke: "none"
  });
  const iconShape = shapeSvg.insert(() => imageNode, ":first-child");
  const outerShape = shapeSvg.insert(() => outerNode);
  if (node.img) {
    const image = shapeSvg.append("image");
    image.attr("href", node.img);
    image.attr("width", imageWidth);
    image.attr("height", imageHeight);
    image.attr("preserveAspectRatio", "none");
    image.attr(
      "transform",
      `translate(${-imageWidth / 2},${topLabel ? outerHeight / 2 - imageHeight : -outerHeight / 2})`
    );
  }
  label.attr(
    "transform",
    `translate(${-bbox.width / 2 - (bbox.x - (bbox.left ?? 0))},${topLabel ? -imageHeight / 2 - bbox.height / 2 - labelPadding / 2 : imageHeight / 2 - bbox.height / 2 + labelPadding / 2})`
  );
  iconShape.attr(
    "transform",
    `translate(${0},${topLabel ? bbox.height / 2 + labelPadding / 2 : -bbox.height / 2 - labelPadding / 2})`
  );
  updateNodeBounds(node, outerShape);
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("iconSquare intersect", node, point);
    if (!node.label) {
      return intersect_default.rect(node, point);
    }
    const dx = node.x ?? 0;
    const dy = node.y ?? 0;
    const nodeHeight = node.height ?? 0;
    let points = [];
    if (topLabel) {
      points = [
        { x: dx - bbox.width / 2, y: dy - nodeHeight / 2 },
        { x: dx + bbox.width / 2, y: dy - nodeHeight / 2 },
        { x: dx + bbox.width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding },
        { x: dx + imageWidth / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding },
        { x: dx + imageWidth / 2, y: dy + nodeHeight / 2 },
        { x: dx - imageWidth / 2, y: dy + nodeHeight / 2 },
        { x: dx - imageWidth / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding },
        { x: dx - bbox.width / 2, y: dy - nodeHeight / 2 + bbox.height + labelPadding }
      ];
    } else {
      points = [
        { x: dx - imageWidth / 2, y: dy - nodeHeight / 2 },
        { x: dx + imageWidth / 2, y: dy - nodeHeight / 2 },
        { x: dx + imageWidth / 2, y: dy - nodeHeight / 2 + imageHeight },
        { x: dx + bbox.width / 2, y: dy - nodeHeight / 2 + imageHeight },
        { x: dx + bbox.width / 2 / 2, y: dy + nodeHeight / 2 },
        { x: dx - bbox.width / 2, y: dy + nodeHeight / 2 },
        { x: dx - bbox.width / 2, y: dy - nodeHeight / 2 + imageHeight },
        { x: dx - imageWidth / 2, y: dy - nodeHeight / 2 + imageHeight }
      ];
    }
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(imageSquare, "imageSquare");

// src/rendering-util/rendering-elements/shapes/invertedTrapezoid.ts

async function inv_trapezoid(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const h = Math.max(bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const points = [
    { x: 0, y: 0 },
    { x: w, y: 0 },
    { x: w + 3 * h / 6, y: -h },
    { x: -3 * h / 6, y: -h }
  ];
  let polygon;
  const { cssStyles } = node;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
    const pathData = createPathFromPoints(points);
    const roughNode = rc.path(pathData, options);
    polygon = shapeSvg.insert(() => roughNode, ":first-child").attr("transform", `translate(${-w / 2}, ${h / 2})`);
    if (cssStyles) {
      polygon.attr("style", cssStyles);
    }
  } else {
    polygon = insertPolygonShape(shapeSvg, w, h, points);
  }
  if (nodeStyles) {
    polygon.attr("style", nodeStyles);
  }
  node.width = w;
  node.height = h;
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    return intersect_default.polygon(node, points, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(inv_trapezoid, "inv_trapezoid");

// src/rendering-util/rendering-elements/shapes/drawRect.ts

async function drawRect(parent, node, options) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const totalWidth = Math.max(bbox.width + options.labelPaddingX * 2, node?.width || 0);
  const totalHeight = Math.max(bbox.height + options.labelPaddingY * 2, node?.height || 0);
  const x = -totalWidth / 2;
  const y = -totalHeight / 2;
  let rect2;
  let { rx, ry } = node;
  const { cssStyles } = node;
  if (options?.rx && options.ry) {
    rx = options.rx;
    ry = options.ry;
  }
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options2 = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
    const roughNode = rx || ry ? rc.path(createRoundedRectPathD(x, y, totalWidth, totalHeight, rx || 0), options2) : rc.rectangle(x, y, totalWidth, totalHeight, options2);
    rect2 = shapeSvg.insert(() => roughNode, ":first-child");
    rect2.attr("class", "basic label-container").attr("style", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(cssStyles));
  } else {
    rect2 = shapeSvg.insert("rect", ":first-child");
    rect2.attr("class", "basic label-container").attr("style", nodeStyles).attr("rx", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(rx)).attr("ry", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(ry)).attr("x", x).attr("y", y).attr("width", totalWidth).attr("height", totalHeight);
  }
  updateNodeBounds(node, rect2);
  node.calcIntersect = function(bounds, point) {
    return intersect_default.rect(bounds, point);
  };
  node.intersect = function(point) {
    return intersect_default.rect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(drawRect, "drawRect");

// src/rendering-util/rendering-elements/shapes/labelRect.ts
async function labelRect(parent, node) {
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, "label");
  const rect2 = shapeSvg.insert("rect", ":first-child");
  const totalWidth = 0.1;
  const totalHeight = 0.1;
  rect2.attr("width", totalWidth).attr("height", totalHeight);
  shapeSvg.attr("class", "label edgeLabel");
  label.attr(
    "transform",
    `translate(${-(bbox.width / 2) - (bbox.x - (bbox.left ?? 0))}, ${-(bbox.height / 2) - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, rect2);
  node.intersect = function(point) {
    return intersect_default.rect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(labelRect, "labelRect");

// src/rendering-util/rendering-elements/shapes/leanLeft.ts

async function lean_left(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + (node.padding ?? 0), node?.width ?? 0);
  const h = Math.max(bbox.height + (node.padding ?? 0), node?.height ?? 0);
  const points = [
    { x: 0, y: 0 },
    { x: w + 3 * h / 6, y: 0 },
    { x: w, y: -h },
    { x: -(3 * h) / 6, y: -h }
  ];
  let polygon;
  const { cssStyles } = node;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
    const pathData = createPathFromPoints(points);
    const roughNode = rc.path(pathData, options);
    polygon = shapeSvg.insert(() => roughNode, ":first-child").attr("transform", `translate(${-w / 2}, ${h / 2})`);
    if (cssStyles) {
      polygon.attr("style", cssStyles);
    }
  } else {
    polygon = insertPolygonShape(shapeSvg, w, h, points);
  }
  if (nodeStyles) {
    polygon.attr("style", nodeStyles);
  }
  node.width = w;
  node.height = h;
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    return intersect_default.polygon(node, points, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(lean_left, "lean_left");

// src/rendering-util/rendering-elements/shapes/leanRight.ts

async function lean_right(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + (node.padding ?? 0), node?.width ?? 0);
  const h = Math.max(bbox.height + (node.padding ?? 0), node?.height ?? 0);
  const points = [
    { x: -3 * h / 6, y: 0 },
    { x: w, y: 0 },
    { x: w + 3 * h / 6, y: -h },
    { x: 0, y: -h }
  ];
  let polygon;
  const { cssStyles } = node;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
    const pathData = createPathFromPoints(points);
    const roughNode = rc.path(pathData, options);
    polygon = shapeSvg.insert(() => roughNode, ":first-child").attr("transform", `translate(${-w / 2}, ${h / 2})`);
    if (cssStyles) {
      polygon.attr("style", cssStyles);
    }
  } else {
    polygon = insertPolygonShape(shapeSvg, w, h, points);
  }
  if (nodeStyles) {
    polygon.attr("style", nodeStyles);
  }
  node.width = w;
  node.height = h;
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    return intersect_default.polygon(node, points, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(lean_right, "lean_right");

// src/rendering-util/rendering-elements/shapes/lightningBolt.ts

function lightningBolt(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.label = "";
  node.labelStyle = labelStyles;
  const shapeSvg = parent.insert("g").attr("class", getNodeClasses(node)).attr("id", node.domId ?? node.id);
  const { cssStyles } = node;
  const width = Math.max(35, node?.width ?? 0);
  const height = Math.max(35, node?.height ?? 0);
  const gap = 7;
  const points = [
    { x: width, y: 0 },
    { x: 0, y: height + gap / 2 },
    { x: width - 2 * gap, y: height + gap / 2 },
    { x: 0, y: 2 * height },
    { x: width, y: height - gap / 2 },
    { x: 2 * gap, y: height - gap / 2 }
  ];
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const linePath = createPathFromPoints(points);
  const lineNode = rc.path(linePath, options);
  const lightningBolt2 = shapeSvg.insert(() => lineNode, ":first-child");
  if (cssStyles && node.look !== "handDrawn") {
    lightningBolt2.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    lightningBolt2.selectAll("path").attr("style", nodeStyles);
  }
  lightningBolt2.attr("transform", `translate(-${width / 2},${-height})`);
  updateNodeBounds(node, lightningBolt2);
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("lightningBolt intersect", node, point);
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(lightningBolt, "lightningBolt");

// src/rendering-util/rendering-elements/shapes/linedCylinder.ts

var createCylinderPathD2 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((x, y, width, height, rx, ry, outerOffset) => {
  return [
    `M${x},${y + ry}`,
    `a${rx},${ry} 0,0,0 ${width},0`,
    `a${rx},${ry} 0,0,0 ${-width},0`,
    `l0,${height}`,
    `a${rx},${ry} 0,0,0 ${width},0`,
    `l0,${-height}`,
    `M${x},${y + ry + outerOffset}`,
    `a${rx},${ry} 0,0,0 ${width},0`
  ].join(" ");
}, "createCylinderPathD");
var createOuterCylinderPathD2 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((x, y, width, height, rx, ry, outerOffset) => {
  return [
    `M${x},${y + ry}`,
    `M${x + width},${y + ry}`,
    `a${rx},${ry} 0,0,0 ${-width},0`,
    `l0,${height}`,
    `a${rx},${ry} 0,0,0 ${width},0`,
    `l0,${-height}`,
    `M${x},${y + ry + outerOffset}`,
    `a${rx},${ry} 0,0,0 ${width},0`
  ].join(" ");
}, "createOuterCylinderPathD");
var createInnerCylinderPathD2 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((x, y, width, height, rx, ry) => {
  return [`M${x - width / 2},${-height / 2}`, `a${rx},${ry} 0,0,0 ${width},0`].join(" ");
}, "createInnerCylinderPathD");
async function linedCylinder(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + (node.padding ?? 0), node.width ?? 0);
  const rx = w / 2;
  const ry = rx / (2.5 + w / 50);
  const h = Math.max(bbox.height + ry + (node.padding ?? 0), node.height ?? 0);
  const outerOffset = h * 0.1;
  let cylinder2;
  const { cssStyles } = node;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const outerPathData = createOuterCylinderPathD2(0, 0, w, h, rx, ry, outerOffset);
    const innerPathData = createInnerCylinderPathD2(0, ry, w, h, rx, ry);
    const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
    const outerNode = rc.path(outerPathData, options);
    const innerLine = rc.path(innerPathData, options);
    const innerLineEl = shapeSvg.insert(() => innerLine, ":first-child");
    innerLineEl.attr("class", "line");
    cylinder2 = shapeSvg.insert(() => outerNode, ":first-child");
    cylinder2.attr("class", "basic label-container");
    if (cssStyles) {
      cylinder2.attr("style", cssStyles);
    }
  } else {
    const pathData = createCylinderPathD2(0, 0, w, h, rx, ry, outerOffset);
    cylinder2 = shapeSvg.insert("path", ":first-child").attr("d", pathData).attr("class", "basic label-container").attr("style", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(cssStyles)).attr("style", nodeStyles);
  }
  cylinder2.attr("label-offset-y", ry);
  cylinder2.attr("transform", `translate(${-w / 2}, ${-(h / 2 + ry)})`);
  updateNodeBounds(node, cylinder2);
  label.attr(
    "transform",
    `translate(${-(bbox.width / 2) - (bbox.x - (bbox.left ?? 0))}, ${-(bbox.height / 2) + ry - (bbox.y - (bbox.top ?? 0))})`
  );
  node.intersect = function(point) {
    const pos = intersect_default.rect(node, point);
    const x = pos.x - (node.x ?? 0);
    if (rx != 0 && (Math.abs(x) < (node.width ?? 0) / 2 || Math.abs(x) == (node.width ?? 0) / 2 && Math.abs(pos.y - (node.y ?? 0)) > (node.height ?? 0) / 2 - ry)) {
      let y = ry * ry * (1 - x * x / (rx * rx));
      if (y > 0) {
        y = Math.sqrt(y);
      }
      y = ry - y;
      if (point.y - (node.y ?? 0) > 0) {
        y = -y;
      }
      pos.y += y;
    }
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(linedCylinder, "linedCylinder");

// src/rendering-util/rendering-elements/shapes/linedWaveEdgedRect.ts

async function linedWaveEdgedRect(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const h = Math.max(bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const waveAmplitude = h / 4;
  const finalH = h + waveAmplitude;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const points = [
    { x: -w / 2 - w / 2 * 0.1, y: -finalH / 2 },
    { x: -w / 2 - w / 2 * 0.1, y: finalH / 2 },
    ...generateFullSineWavePoints(
      -w / 2 - w / 2 * 0.1,
      finalH / 2,
      w / 2 + w / 2 * 0.1,
      finalH / 2,
      waveAmplitude,
      0.8
    ),
    { x: w / 2 + w / 2 * 0.1, y: -finalH / 2 },
    { x: -w / 2 - w / 2 * 0.1, y: -finalH / 2 },
    { x: -w / 2, y: -finalH / 2 },
    { x: -w / 2, y: finalH / 2 * 1.1 },
    { x: -w / 2, y: -finalH / 2 }
  ];
  const poly = rc.polygon(
    points.map((p) => [p.x, p.y]),
    options
  );
  const waveEdgeRect = shapeSvg.insert(() => poly, ":first-child");
  waveEdgeRect.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    waveEdgeRect.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    waveEdgeRect.selectAll("path").attr("style", nodeStyles);
  }
  waveEdgeRect.attr("transform", `translate(0,${-waveAmplitude / 2})`);
  label.attr(
    "transform",
    `translate(${-w / 2 + (node.padding ?? 0) + w / 2 * 0.1 / 2 - (bbox.x - (bbox.left ?? 0))},${-h / 2 + (node.padding ?? 0) - waveAmplitude / 2 - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, waveEdgeRect);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(linedWaveEdgedRect, "linedWaveEdgedRect");

// src/rendering-util/rendering-elements/shapes/multiRect.ts

async function multiRect(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const h = Math.max(bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const rectOffset = 5;
  const x = -w / 2;
  const y = -h / 2;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  const outerPathPoints = [
    { x: x - rectOffset, y: y + rectOffset },
    { x: x - rectOffset, y: y + h + rectOffset },
    { x: x + w - rectOffset, y: y + h + rectOffset },
    { x: x + w - rectOffset, y: y + h },
    { x: x + w, y: y + h },
    { x: x + w, y: y + h - rectOffset },
    { x: x + w + rectOffset, y: y + h - rectOffset },
    { x: x + w + rectOffset, y: y - rectOffset },
    { x: x + rectOffset, y: y - rectOffset },
    { x: x + rectOffset, y },
    { x, y },
    { x, y: y + rectOffset }
  ];
  const innerPathPoints = [
    { x, y: y + rectOffset },
    { x: x + w - rectOffset, y: y + rectOffset },
    { x: x + w - rectOffset, y: y + h },
    { x: x + w, y: y + h },
    { x: x + w, y },
    { x, y }
  ];
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const outerPath = createPathFromPoints(outerPathPoints);
  const outerNode = rc.path(outerPath, options);
  const innerPath = createPathFromPoints(innerPathPoints);
  const innerNode = rc.path(innerPath, { ...options, fill: "none" });
  const multiRect2 = shapeSvg.insert(() => innerNode, ":first-child");
  multiRect2.insert(() => outerNode, ":first-child");
  multiRect2.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    multiRect2.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    multiRect2.selectAll("path").attr("style", nodeStyles);
  }
  label.attr(
    "transform",
    `translate(${-(bbox.width / 2) - rectOffset - (bbox.x - (bbox.left ?? 0))}, ${-(bbox.height / 2) + rectOffset - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, multiRect2);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, outerPathPoints, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(multiRect, "multiRect");

// src/rendering-util/rendering-elements/shapes/multiWaveEdgedRectangle.ts

async function multiWaveEdgedRectangle(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const h = Math.max(bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const waveAmplitude = h / 4;
  const finalH = h + waveAmplitude;
  const x = -w / 2;
  const y = -finalH / 2;
  const rectOffset = 5;
  const { cssStyles } = node;
  const wavePoints = generateFullSineWavePoints(
    x - rectOffset,
    y + finalH + rectOffset,
    x + w - rectOffset,
    y + finalH + rectOffset,
    waveAmplitude,
    0.8
  );
  const lastWavePoint = wavePoints?.[wavePoints.length - 1];
  const outerPathPoints = [
    { x: x - rectOffset, y: y + rectOffset },
    { x: x - rectOffset, y: y + finalH + rectOffset },
    ...wavePoints,
    { x: x + w - rectOffset, y: lastWavePoint.y - rectOffset },
    { x: x + w, y: lastWavePoint.y - rectOffset },
    { x: x + w, y: lastWavePoint.y - 2 * rectOffset },
    { x: x + w + rectOffset, y: lastWavePoint.y - 2 * rectOffset },
    { x: x + w + rectOffset, y: y - rectOffset },
    { x: x + rectOffset, y: y - rectOffset },
    { x: x + rectOffset, y },
    { x, y },
    { x, y: y + rectOffset }
  ];
  const innerPathPoints = [
    { x, y: y + rectOffset },
    { x: x + w - rectOffset, y: y + rectOffset },
    { x: x + w - rectOffset, y: lastWavePoint.y - rectOffset },
    { x: x + w, y: lastWavePoint.y - rectOffset },
    { x: x + w, y },
    { x, y }
  ];
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const outerPath = createPathFromPoints(outerPathPoints);
  const outerNode = rc.path(outerPath, options);
  const innerPath = createPathFromPoints(innerPathPoints);
  const innerNode = rc.path(innerPath, options);
  const shape = shapeSvg.insert(() => outerNode, ":first-child");
  shape.insert(() => innerNode);
  shape.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    shape.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    shape.selectAll("path").attr("style", nodeStyles);
  }
  shape.attr("transform", `translate(0,${-waveAmplitude / 2})`);
  label.attr(
    "transform",
    `translate(${-(bbox.width / 2) - rectOffset - (bbox.x - (bbox.left ?? 0))}, ${-(bbox.height / 2) + rectOffset - waveAmplitude / 2 - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, shape);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, outerPathPoints, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(multiWaveEdgedRectangle, "multiWaveEdgedRectangle");

// src/rendering-util/rendering-elements/shapes/note.ts

async function note(parent, node, { config: { themeVariables } }) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const useHtmlLabels = node.useHtmlLabels || (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig */ .iE)().flowchart?.htmlLabels !== false;
  if (!useHtmlLabels) {
    node.centerLabel = true;
  }
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const totalWidth = Math.max(bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const totalHeight = Math.max(bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const x = -totalWidth / 2;
  const y = -totalHeight / 2;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {
    fill: themeVariables.noteBkgColor,
    stroke: themeVariables.noteBorderColor
  });
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const noteShapeNode = rc.rectangle(x, y, totalWidth, totalHeight, options);
  const rect2 = shapeSvg.insert(() => noteShapeNode, ":first-child");
  rect2.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    rect2.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    rect2.selectAll("path").attr("style", nodeStyles);
  }
  label.attr(
    "transform",
    `translate(${-bbox.width / 2 - (bbox.x - (bbox.left ?? 0))}, ${-(bbox.height / 2) - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, rect2);
  node.intersect = function(point) {
    return intersect_default.rect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(note, "note");

// src/rendering-util/rendering-elements/shapes/question.ts

var createDecisionBoxPathD = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((x, y, size) => {
  return [
    `M${x + size / 2},${y}`,
    `L${x + size},${y - size / 2}`,
    `L${x + size / 2},${y - size}`,
    `L${x},${y - size / 2}`,
    "Z"
  ].join(" ");
}, "createDecisionBoxPathD");
async function question(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const w = bbox.width + node.padding;
  const h = bbox.height + node.padding;
  const s = w + h;
  const adjustment = 0.5;
  const points = [
    { x: s / 2, y: 0 },
    { x: s, y: -s / 2 },
    { x: s / 2, y: -s },
    { x: 0, y: -s / 2 }
  ];
  let polygon;
  const { cssStyles } = node;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
    const pathData = createDecisionBoxPathD(0, 0, s);
    const roughNode = rc.path(pathData, options);
    polygon = shapeSvg.insert(() => roughNode, ":first-child").attr("transform", `translate(${-s / 2 + adjustment}, ${s / 2})`);
    if (cssStyles) {
      polygon.attr("style", cssStyles);
    }
  } else {
    polygon = insertPolygonShape(shapeSvg, s, s, points);
    polygon.attr("transform", `translate(${-s / 2 + adjustment}, ${s / 2})`);
  }
  if (nodeStyles) {
    polygon.attr("style", nodeStyles);
  }
  updateNodeBounds(node, polygon);
  node.calcIntersect = function(bounds, point) {
    const s2 = bounds.width;
    const points2 = [
      { x: s2 / 2, y: 0 },
      { x: s2, y: -s2 / 2 },
      { x: s2 / 2, y: -s2 },
      { x: 0, y: -s2 / 2 }
    ];
    const res = intersect_default.polygon(bounds, points2, point);
    return { x: res.x - 0.5, y: res.y - 0.5 };
  };
  node.intersect = function(point) {
    return this.calcIntersect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(question, "question");

// src/rendering-util/rendering-elements/shapes/rectLeftInvArrow.ts

async function rect_left_inv_arrow(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + (node.padding ?? 0), node?.width ?? 0);
  const h = Math.max(bbox.height + (node.padding ?? 0), node?.height ?? 0);
  const x = -w / 2;
  const y = -h / 2;
  const notch = y / 2;
  const points = [
    { x: x + notch, y },
    { x, y: 0 },
    { x: x + notch, y: -y },
    { x: -x, y: -y },
    { x: -x, y }
  ];
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const pathData = createPathFromPoints(points);
  const roughNode = rc.path(pathData, options);
  const polygon = shapeSvg.insert(() => roughNode, ":first-child");
  polygon.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    polygon.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    polygon.selectAll("path").attr("style", nodeStyles);
  }
  polygon.attr("transform", `translate(${-notch / 2},0)`);
  label.attr(
    "transform",
    `translate(${-notch / 2 - bbox.width / 2 - (bbox.x - (bbox.left ?? 0))}, ${-(bbox.height / 2) - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    return intersect_default.polygon(node, points, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(rect_left_inv_arrow, "rect_left_inv_arrow");

// src/rendering-util/rendering-elements/shapes/rectWithTitle.ts


async function rectWithTitle(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  let classes;
  if (!node.cssClasses) {
    classes = "node default";
  } else {
    classes = "node " + node.cssClasses;
  }
  const shapeSvg = parent.insert("g").attr("class", classes).attr("id", node.domId || node.id);
  const g = shapeSvg.insert("g");
  const label = shapeSvg.insert("g").attr("class", "label").attr("style", nodeStyles);
  const description = node.description;
  const title = node.label;
  const text2 = label.node().appendChild(await createLabel_default(title, node.labelStyle, true, true));
  let bbox = { width: 0, height: 0 };
  if ((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)()?.flowchart?.htmlLabels)) {
    const div2 = text2.children[0];
    const dv2 = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(text2);
    bbox = div2.getBoundingClientRect();
    dv2.attr("width", bbox.width);
    dv2.attr("height", bbox.height);
  }
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("Text 2", description);
  const textRows = description || [];
  const titleBox = text2.getBBox();
  const descr = label.node().appendChild(
    await createLabel_default(
      textRows.join ? textRows.join("<br/>") : textRows,
      node.labelStyle,
      true,
      true
    )
  );
  const div = descr.children[0];
  const dv = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(descr);
  bbox = div.getBoundingClientRect();
  dv.attr("width", bbox.width);
  dv.attr("height", bbox.height);
  const halfPadding = (node.padding || 0) / 2;
  (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(descr).attr(
    "transform",
    "translate( " + (bbox.width > titleBox.width ? 0 : (titleBox.width - bbox.width) / 2) + ", " + (titleBox.height + halfPadding + 5) + ")"
  );
  (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(text2).attr(
    "transform",
    "translate( " + (bbox.width < titleBox.width ? 0 : -(titleBox.width - bbox.width) / 2) + ", 0)"
  );
  bbox = label.node().getBBox();
  label.attr(
    "transform",
    "translate(" + -bbox.width / 2 + ", " + (-bbox.height / 2 - halfPadding + 3) + ")"
  );
  const totalWidth = bbox.width + (node.padding || 0);
  const totalHeight = bbox.height + (node.padding || 0);
  const x = -bbox.width / 2 - halfPadding;
  const y = -bbox.height / 2 - halfPadding;
  let rect2;
  let innerLine;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
    const roughNode = rc.path(
      createRoundedRectPathD(x, y, totalWidth, totalHeight, node.rx || 0),
      options
    );
    const roughLine = rc.line(
      -bbox.width / 2 - halfPadding,
      -bbox.height / 2 - halfPadding + titleBox.height + halfPadding,
      bbox.width / 2 + halfPadding,
      -bbox.height / 2 - halfPadding + titleBox.height + halfPadding,
      options
    );
    innerLine = shapeSvg.insert(() => {
      _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.debug("Rough node insert CXC", roughNode);
      return roughLine;
    }, ":first-child");
    rect2 = shapeSvg.insert(() => {
      _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.debug("Rough node insert CXC", roughNode);
      return roughNode;
    }, ":first-child");
  } else {
    rect2 = g.insert("rect", ":first-child");
    innerLine = g.insert("line");
    rect2.attr("class", "outer title-state").attr("style", nodeStyles).attr("x", -bbox.width / 2 - halfPadding).attr("y", -bbox.height / 2 - halfPadding).attr("width", bbox.width + (node.padding || 0)).attr("height", bbox.height + (node.padding || 0));
    innerLine.attr("class", "divider").attr("x1", -bbox.width / 2 - halfPadding).attr("x2", bbox.width / 2 + halfPadding).attr("y1", -bbox.height / 2 - halfPadding + titleBox.height + halfPadding).attr("y2", -bbox.height / 2 - halfPadding + titleBox.height + halfPadding);
  }
  updateNodeBounds(node, rect2);
  node.intersect = function(point) {
    return intersect_default.rect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(rectWithTitle, "rectWithTitle");

// src/rendering-util/rendering-elements/shapes/roundedRect.ts

function generateArcPoints2(x1, y1, x2, y2, rx, ry, clockwise) {
  const numPoints = 20;
  const midX = (x1 + x2) / 2;
  const midY = (y1 + y2) / 2;
  const angle = Math.atan2(y2 - y1, x2 - x1);
  const dx = (x2 - x1) / 2;
  const dy = (y2 - y1) / 2;
  const transformedX = dx / rx;
  const transformedY = dy / ry;
  const distance = Math.sqrt(transformedX ** 2 + transformedY ** 2);
  if (distance > 1) {
    throw new Error("The given radii are too small to create an arc between the points.");
  }
  const scaledCenterDistance = Math.sqrt(1 - distance ** 2);
  const centerX = midX + scaledCenterDistance * ry * Math.sin(angle) * (clockwise ? -1 : 1);
  const centerY = midY - scaledCenterDistance * rx * Math.cos(angle) * (clockwise ? -1 : 1);
  const startAngle = Math.atan2((y1 - centerY) / ry, (x1 - centerX) / rx);
  const endAngle = Math.atan2((y2 - centerY) / ry, (x2 - centerX) / rx);
  let angleRange = endAngle - startAngle;
  if (clockwise && angleRange < 0) {
    angleRange += 2 * Math.PI;
  }
  if (!clockwise && angleRange > 0) {
    angleRange -= 2 * Math.PI;
  }
  const points = [];
  for (let i = 0; i < numPoints; i++) {
    const t = i / (numPoints - 1);
    const angle2 = startAngle + t * angleRange;
    const x = centerX + rx * Math.cos(angle2);
    const y = centerY + ry * Math.sin(angle2);
    points.push({ x, y });
  }
  return points;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(generateArcPoints2, "generateArcPoints");
async function roundedRect(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const labelPaddingX = node?.padding ?? 0;
  const labelPaddingY = node?.padding ?? 0;
  const w = (node?.width ? node?.width : bbox.width) + labelPaddingX * 2;
  const h = (node?.height ? node?.height : bbox.height) + labelPaddingY * 2;
  const radius = node.radius || 5;
  const taper = node.taper || 5;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.stroke) {
    options.stroke = node.stroke;
  }
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const points = [
    // Top edge (left to right)
    { x: -w / 2 + taper, y: -h / 2 },
    // Top-left corner start (1)
    { x: w / 2 - taper, y: -h / 2 },
    // Top-right corner start (2)
    ...generateArcPoints2(w / 2 - taper, -h / 2, w / 2, -h / 2 + taper, radius, radius, true),
    // Top-left arc (2 to 3)
    // Right edge (top to bottom)
    { x: w / 2, y: -h / 2 + taper },
    // Top-right taper point (3)
    { x: w / 2, y: h / 2 - taper },
    // Bottom-right taper point (4)
    ...generateArcPoints2(w / 2, h / 2 - taper, w / 2 - taper, h / 2, radius, radius, true),
    // Top-left arc (4 to 5)
    // Bottom edge (right to left)
    { x: w / 2 - taper, y: h / 2 },
    // Bottom-right corner start (5)
    { x: -w / 2 + taper, y: h / 2 },
    // Bottom-left corner start (6)
    ...generateArcPoints2(-w / 2 + taper, h / 2, -w / 2, h / 2 - taper, radius, radius, true),
    // Top-left arc (4 to 5)
    // Left edge (bottom to top)
    { x: -w / 2, y: h / 2 - taper },
    // Bottom-left taper point (7)
    { x: -w / 2, y: -h / 2 + taper },
    // Top-left taper point (8)
    ...generateArcPoints2(-w / 2, -h / 2 + taper, -w / 2 + taper, -h / 2, radius, radius, true)
    // Top-left arc (4 to 5)
  ];
  const pathData = createPathFromPoints(points);
  const shapeNode = rc.path(pathData, options);
  const polygon = shapeSvg.insert(() => shapeNode, ":first-child");
  polygon.attr("class", "basic label-container outer-path");
  if (cssStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", nodeStyles);
  }
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(roundedRect, "roundedRect");

// src/rendering-util/rendering-elements/shapes/shadedProcess.ts

async function shadedProcess(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const halfPadding = node?.padding ?? 0;
  const w = Math.max(bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const h = Math.max(bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const x = -bbox.width / 2 - halfPadding;
  const y = -bbox.height / 2 - halfPadding;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const points = [
    { x, y },
    { x: x + w + 8, y },
    { x: x + w + 8, y: y + h },
    { x: x - 8, y: y + h },
    { x: x - 8, y },
    { x, y },
    { x, y: y + h }
  ];
  const roughNode = rc.polygon(
    points.map((p) => [p.x, p.y]),
    options
  );
  const rect2 = shapeSvg.insert(() => roughNode, ":first-child");
  rect2.attr("class", "basic label-container").attr("style", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(cssStyles));
  if (nodeStyles && node.look !== "handDrawn") {
    rect2.selectAll("path").attr("style", nodeStyles);
  }
  if (cssStyles && node.look !== "handDrawn") {
    rect2.selectAll("path").attr("style", nodeStyles);
  }
  label.attr(
    "transform",
    `translate(${-w / 2 + 4 + (node.padding ?? 0) - (bbox.x - (bbox.left ?? 0))},${-h / 2 + (node.padding ?? 0) - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, rect2);
  node.intersect = function(point) {
    return intersect_default.rect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(shadedProcess, "shadedProcess");

// src/rendering-util/rendering-elements/shapes/slopedRect.ts

async function slopedRect(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const h = Math.max(bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const x = -w / 2;
  const y = -h / 2;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const points = [
    { x, y },
    { x, y: y + h },
    { x: x + w, y: y + h },
    { x: x + w, y: y - h / 2 }
  ];
  const pathData = createPathFromPoints(points);
  const shapeNode = rc.path(pathData, options);
  const polygon = shapeSvg.insert(() => shapeNode, ":first-child");
  polygon.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", nodeStyles);
  }
  polygon.attr("transform", `translate(0, ${h / 4})`);
  label.attr(
    "transform",
    `translate(${-w / 2 + (node.padding ?? 0) - (bbox.x - (bbox.left ?? 0))}, ${-h / 4 + (node.padding ?? 0) - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(slopedRect, "slopedRect");

// src/rendering-util/rendering-elements/shapes/squareRect.ts
async function squareRect2(parent, node) {
  const options = {
    rx: 0,
    ry: 0,
    classes: "",
    labelPaddingX: node.labelPaddingX ?? (node?.padding || 0) * 2,
    labelPaddingY: (node?.padding || 0) * 1
  };
  return drawRect(parent, node, options);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(squareRect2, "squareRect");

// src/rendering-util/rendering-elements/shapes/stadium.ts

async function stadium(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const h = bbox.height + node.padding;
  const w = bbox.width + h / 4 + node.padding;
  const radius = h / 2;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const points = [
    { x: -w / 2 + radius, y: -h / 2 },
    { x: w / 2 - radius, y: -h / 2 },
    ...generateCirclePoints(-w / 2 + radius, 0, radius, 50, 90, 270),
    { x: w / 2 - radius, y: h / 2 },
    ...generateCirclePoints(w / 2 - radius, 0, radius, 50, 270, 450)
  ];
  const pathData = createPathFromPoints(points);
  const shapeNode = rc.path(pathData, options);
  const polygon = shapeSvg.insert(() => shapeNode, ":first-child");
  polygon.attr("class", "basic label-container outer-path");
  if (cssStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", nodeStyles);
  }
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(stadium, "stadium");

// src/rendering-util/rendering-elements/shapes/state.ts
async function state(parent, node) {
  const options = {
    rx: 5,
    ry: 5,
    classes: "flowchart-node"
  };
  return drawRect(parent, node, options);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(state, "state");

// src/rendering-util/rendering-elements/shapes/stateEnd.ts

function stateEnd(parent, node, { config: { themeVariables } }) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { cssStyles } = node;
  const { lineColor, stateBorder, nodeBorder } = themeVariables;
  const shapeSvg = parent.insert("g").attr("class", "node default").attr("id", node.domId || node.id);
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const roughNode = rc.circle(0, 0, 14, {
    ...options,
    stroke: lineColor,
    strokeWidth: 2
  });
  const innerFill = stateBorder ?? nodeBorder;
  const roughInnerNode = rc.circle(0, 0, 5, {
    ...options,
    fill: innerFill,
    stroke: innerFill,
    strokeWidth: 2,
    fillStyle: "solid"
  });
  const circle2 = shapeSvg.insert(() => roughNode, ":first-child");
  circle2.insert(() => roughInnerNode);
  if (cssStyles) {
    circle2.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles) {
    circle2.selectAll("path").attr("style", nodeStyles);
  }
  updateNodeBounds(node, circle2);
  node.intersect = function(point) {
    return intersect_default.circle(node, 7, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(stateEnd, "stateEnd");

// src/rendering-util/rendering-elements/shapes/stateStart.ts

function stateStart(parent, node, { config: { themeVariables } }) {
  const { lineColor } = themeVariables;
  const shapeSvg = parent.insert("g").attr("class", "node default").attr("id", node.domId || node.id);
  let circle2;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const roughNode = rc.circle(0, 0, 14, (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .solidStateFill */ .TA)(lineColor));
    circle2 = shapeSvg.insert(() => roughNode);
    circle2.attr("class", "state-start").attr("r", 7).attr("width", 14).attr("height", 14);
  } else {
    circle2 = shapeSvg.insert("circle", ":first-child");
    circle2.attr("class", "state-start").attr("r", 7).attr("width", 14).attr("height", 14);
  }
  updateNodeBounds(node, circle2);
  node.intersect = function(point) {
    return intersect_default.circle(node, 7, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(stateStart, "stateStart");

// src/rendering-util/rendering-elements/shapes/subroutine.ts

async function subroutine(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const halfPadding = (node?.padding || 0) / 2;
  const w = bbox.width + node.padding;
  const h = bbox.height + node.padding;
  const x = -bbox.width / 2 - halfPadding;
  const y = -bbox.height / 2 - halfPadding;
  const points = [
    { x: 0, y: 0 },
    { x: w, y: 0 },
    { x: w, y: -h },
    { x: 0, y: -h },
    { x: 0, y: 0 },
    { x: -8, y: 0 },
    { x: w + 8, y: 0 },
    { x: w + 8, y: -h },
    { x: -8, y: -h },
    { x: -8, y: 0 }
  ];
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
    const roughNode = rc.rectangle(x - 8, y, w + 16, h, options);
    const l1 = rc.line(x, y, x, y + h, options);
    const l2 = rc.line(x + w, y, x + w, y + h, options);
    shapeSvg.insert(() => l1, ":first-child");
    shapeSvg.insert(() => l2, ":first-child");
    const rect2 = shapeSvg.insert(() => roughNode, ":first-child");
    const { cssStyles } = node;
    rect2.attr("class", "basic label-container").attr("style", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(cssStyles));
    updateNodeBounds(node, rect2);
  } else {
    const el = insertPolygonShape(shapeSvg, w, h, points);
    if (nodeStyles) {
      el.attr("style", nodeStyles);
    }
    updateNodeBounds(node, el);
  }
  node.intersect = function(point) {
    return intersect_default.polygon(node, points, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(subroutine, "subroutine");

// src/rendering-util/rendering-elements/shapes/taggedRect.ts

async function taggedRect(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const h = Math.max(bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const x = -w / 2;
  const y = -h / 2;
  const tagWidth = 0.2 * h;
  const tagHeight = 0.2 * h;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  const rectPoints = [
    { x: x - tagWidth / 2, y },
    { x: x + w + tagWidth / 2, y },
    { x: x + w + tagWidth / 2, y: y + h },
    { x: x - tagWidth / 2, y: y + h }
  ];
  const tagPoints = [
    { x: x + w - tagWidth / 2, y: y + h },
    { x: x + w + tagWidth / 2, y: y + h },
    { x: x + w + tagWidth / 2, y: y + h - tagHeight }
  ];
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const rectPath = createPathFromPoints(rectPoints);
  const rectNode = rc.path(rectPath, options);
  const tagPath = createPathFromPoints(tagPoints);
  const tagNode = rc.path(tagPath, { ...options, fillStyle: "solid" });
  const taggedRect2 = shapeSvg.insert(() => tagNode, ":first-child");
  taggedRect2.insert(() => rectNode, ":first-child");
  taggedRect2.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    taggedRect2.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    taggedRect2.selectAll("path").attr("style", nodeStyles);
  }
  updateNodeBounds(node, taggedRect2);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, rectPoints, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(taggedRect, "taggedRect");

// src/rendering-util/rendering-elements/shapes/taggedWaveEdgedRectangle.ts

async function taggedWaveEdgedRectangle(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const h = Math.max(bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const waveAmplitude = h / 4;
  const tagWidth = 0.2 * w;
  const tagHeight = 0.2 * h;
  const finalH = h + waveAmplitude;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const points = [
    { x: -w / 2 - w / 2 * 0.1, y: finalH / 2 },
    ...generateFullSineWavePoints(
      -w / 2 - w / 2 * 0.1,
      finalH / 2,
      w / 2 + w / 2 * 0.1,
      finalH / 2,
      waveAmplitude,
      0.8
    ),
    { x: w / 2 + w / 2 * 0.1, y: -finalH / 2 },
    { x: -w / 2 - w / 2 * 0.1, y: -finalH / 2 }
  ];
  const x = -w / 2 + w / 2 * 0.1;
  const y = -finalH / 2 - tagHeight * 0.4;
  const tagPoints = [
    { x: x + w - tagWidth, y: (y + h) * 1.4 },
    { x: x + w, y: y + h - tagHeight },
    { x: x + w, y: (y + h) * 0.9 },
    ...generateFullSineWavePoints(
      x + w,
      (y + h) * 1.3,
      x + w - tagWidth,
      (y + h) * 1.5,
      -h * 0.03,
      0.5
    )
  ];
  const waveEdgeRectPath = createPathFromPoints(points);
  const waveEdgeRectNode = rc.path(waveEdgeRectPath, options);
  const taggedWaveEdgeRectPath = createPathFromPoints(tagPoints);
  const taggedWaveEdgeRectNode = rc.path(taggedWaveEdgeRectPath, {
    ...options,
    fillStyle: "solid"
  });
  const waveEdgeRect = shapeSvg.insert(() => taggedWaveEdgeRectNode, ":first-child");
  waveEdgeRect.insert(() => waveEdgeRectNode, ":first-child");
  waveEdgeRect.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    waveEdgeRect.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    waveEdgeRect.selectAll("path").attr("style", nodeStyles);
  }
  waveEdgeRect.attr("transform", `translate(0,${-waveAmplitude / 2})`);
  label.attr(
    "transform",
    `translate(${-w / 2 + (node.padding ?? 0) - (bbox.x - (bbox.left ?? 0))},${-h / 2 + (node.padding ?? 0) - waveAmplitude / 2 - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, waveEdgeRect);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(taggedWaveEdgedRectangle, "taggedWaveEdgedRectangle");

// src/rendering-util/rendering-elements/shapes/text.ts
async function text(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const totalWidth = Math.max(bbox.width + node.padding, node?.width || 0);
  const totalHeight = Math.max(bbox.height + node.padding, node?.height || 0);
  const x = -totalWidth / 2;
  const y = -totalHeight / 2;
  const rect2 = shapeSvg.insert("rect", ":first-child");
  rect2.attr("class", "text").attr("style", nodeStyles).attr("rx", 0).attr("ry", 0).attr("x", x).attr("y", y).attr("width", totalWidth).attr("height", totalHeight);
  updateNodeBounds(node, rect2);
  node.intersect = function(point) {
    return intersect_default.rect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(text, "text");

// src/rendering-util/rendering-elements/shapes/tiltedCylinder.ts

var createCylinderPathD3 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((x, y, width, height, rx, ry) => {
  return `M${x},${y}
    a${rx},${ry} 0,0,1 ${0},${-height}
    l${width},${0}
    a${rx},${ry} 0,0,1 ${0},${height}
    M${width},${-height}
    a${rx},${ry} 0,0,0 ${0},${height}
    l${-width},${0}`;
}, "createCylinderPathD");
var createOuterCylinderPathD3 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((x, y, width, height, rx, ry) => {
  return [
    `M${x},${y}`,
    `M${x + width},${y}`,
    `a${rx},${ry} 0,0,0 ${0},${-height}`,
    `l${-width},0`,
    `a${rx},${ry} 0,0,0 ${0},${height}`,
    `l${width},0`
  ].join(" ");
}, "createOuterCylinderPathD");
var createInnerCylinderPathD3 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((x, y, width, height, rx, ry) => {
  return [`M${x + width / 2},${-height / 2}`, `a${rx},${ry} 0,0,0 0,${height}`].join(" ");
}, "createInnerCylinderPathD");
async function tiltedCylinder(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label, halfPadding } = await labelHelper(
    parent,
    node,
    getNodeClasses(node)
  );
  const labelPadding = node.look === "neo" ? halfPadding * 2 : halfPadding;
  const h = bbox.height + labelPadding;
  const ry = h / 2;
  const rx = ry / (2.5 + h / 50);
  const w = bbox.width + rx + labelPadding;
  const { cssStyles } = node;
  let cylinder2;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const outerPathData = createOuterCylinderPathD3(0, 0, w, h, rx, ry);
    const innerPathData = createInnerCylinderPathD3(0, 0, w, h, rx, ry);
    const outerNode = rc.path(outerPathData, (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {}));
    const innerLine = rc.path(innerPathData, (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, { fill: "none" }));
    cylinder2 = shapeSvg.insert(() => innerLine, ":first-child");
    cylinder2 = shapeSvg.insert(() => outerNode, ":first-child");
    cylinder2.attr("class", "basic label-container");
    if (cssStyles) {
      cylinder2.attr("style", cssStyles);
    }
  } else {
    const pathData = createCylinderPathD3(0, 0, w, h, rx, ry);
    cylinder2 = shapeSvg.insert("path", ":first-child").attr("d", pathData).attr("class", "basic label-container").attr("style", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(cssStyles)).attr("style", nodeStyles);
    cylinder2.attr("class", "basic label-container");
    if (cssStyles) {
      cylinder2.selectAll("path").attr("style", cssStyles);
    }
    if (nodeStyles) {
      cylinder2.selectAll("path").attr("style", nodeStyles);
    }
  }
  cylinder2.attr("label-offset-x", rx);
  cylinder2.attr("transform", `translate(${-w / 2}, ${h / 2} )`);
  label.attr(
    "transform",
    `translate(${-(bbox.width / 2) - rx - (bbox.x - (bbox.left ?? 0))}, ${-(bbox.height / 2) - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, cylinder2);
  node.intersect = function(point) {
    const pos = intersect_default.rect(node, point);
    const y = pos.y - (node.y ?? 0);
    if (ry != 0 && (Math.abs(y) < (node.height ?? 0) / 2 || Math.abs(y) == (node.height ?? 0) / 2 && Math.abs(pos.x - (node.x ?? 0)) > (node.width ?? 0) / 2 - rx)) {
      let x = rx * rx * (1 - y * y / (ry * ry));
      if (x != 0) {
        x = Math.sqrt(Math.abs(x));
      }
      x = rx - x;
      if (point.x - (node.x ?? 0) > 0) {
        x = -x;
      }
      pos.x += x;
    }
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(tiltedCylinder, "tiltedCylinder");

// src/rendering-util/rendering-elements/shapes/trapezoid.ts

async function trapezoid(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const w = bbox.width + node.padding;
  const h = bbox.height + node.padding;
  const points = [
    { x: -3 * h / 6, y: 0 },
    { x: w + 3 * h / 6, y: 0 },
    { x: w, y: -h },
    { x: 0, y: -h }
  ];
  let polygon;
  const { cssStyles } = node;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
    const pathData = createPathFromPoints(points);
    const roughNode = rc.path(pathData, options);
    polygon = shapeSvg.insert(() => roughNode, ":first-child").attr("transform", `translate(${-w / 2}, ${h / 2})`);
    if (cssStyles) {
      polygon.attr("style", cssStyles);
    }
  } else {
    polygon = insertPolygonShape(shapeSvg, w, h, points);
  }
  if (nodeStyles) {
    polygon.attr("style", nodeStyles);
  }
  node.width = w;
  node.height = h;
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    return intersect_default.polygon(node, points, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(trapezoid, "trapezoid");

// src/rendering-util/rendering-elements/shapes/trapezoidalPentagon.ts

async function trapezoidalPentagon(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const minWidth = 60, minHeight = 20;
  const w = Math.max(minWidth, bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const h = Math.max(minHeight, bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const points = [
    { x: -w / 2 * 0.8, y: -h / 2 },
    { x: w / 2 * 0.8, y: -h / 2 },
    { x: w / 2, y: -h / 2 * 0.6 },
    { x: w / 2, y: h / 2 },
    { x: -w / 2, y: h / 2 },
    { x: -w / 2, y: -h / 2 * 0.6 }
  ];
  const pathData = createPathFromPoints(points);
  const shapeNode = rc.path(pathData, options);
  const polygon = shapeSvg.insert(() => shapeNode, ":first-child");
  polygon.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", nodeStyles);
  }
  updateNodeBounds(node, polygon);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(trapezoidalPentagon, "trapezoidalPentagon");

// src/rendering-util/rendering-elements/shapes/triangle.ts

async function triangle(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const useHtmlLabels = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)().flowchart?.htmlLabels);
  const w = bbox.width + (node.padding ?? 0);
  const h = w + bbox.height;
  const tw = w + bbox.height;
  const points = [
    { x: 0, y: 0 },
    { x: tw, y: 0 },
    { x: tw / 2, y: -h }
  ];
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const pathData = createPathFromPoints(points);
  const roughNode = rc.path(pathData, options);
  const polygon = shapeSvg.insert(() => roughNode, ":first-child").attr("transform", `translate(${-h / 2}, ${h / 2})`);
  if (cssStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    polygon.selectChildren("path").attr("style", nodeStyles);
  }
  node.width = w;
  node.height = h;
  updateNodeBounds(node, polygon);
  label.attr(
    "transform",
    `translate(${-bbox.width / 2 - (bbox.x - (bbox.left ?? 0))}, ${h / 2 - (bbox.height + (node.padding ?? 0) / (useHtmlLabels ? 2 : 1) - (bbox.y - (bbox.top ?? 0)))})`
  );
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("Triangle intersect", node, points, point);
    return intersect_default.polygon(node, points, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(triangle, "triangle");

// src/rendering-util/rendering-elements/shapes/waveEdgedRectangle.ts

async function waveEdgedRectangle(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const h = Math.max(bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const waveAmplitude = h / 8;
  const finalH = h + waveAmplitude;
  const { cssStyles } = node;
  const minWidth = 70;
  const widthDif = minWidth - w;
  const extraW = widthDif > 0 ? widthDif / 2 : 0;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const points = [
    { x: -w / 2 - extraW, y: finalH / 2 },
    ...generateFullSineWavePoints(
      -w / 2 - extraW,
      finalH / 2,
      w / 2 + extraW,
      finalH / 2,
      waveAmplitude,
      0.8
    ),
    { x: w / 2 + extraW, y: -finalH / 2 },
    { x: -w / 2 - extraW, y: -finalH / 2 }
  ];
  const waveEdgeRectPath = createPathFromPoints(points);
  const waveEdgeRectNode = rc.path(waveEdgeRectPath, options);
  const waveEdgeRect = shapeSvg.insert(() => waveEdgeRectNode, ":first-child");
  waveEdgeRect.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    waveEdgeRect.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    waveEdgeRect.selectAll("path").attr("style", nodeStyles);
  }
  waveEdgeRect.attr("transform", `translate(0,${-waveAmplitude / 2})`);
  label.attr(
    "transform",
    `translate(${-w / 2 + (node.padding ?? 0) - (bbox.x - (bbox.left ?? 0))},${-h / 2 + (node.padding ?? 0) - waveAmplitude - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, waveEdgeRect);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(waveEdgedRectangle, "waveEdgedRectangle");

// src/rendering-util/rendering-elements/shapes/waveRectangle.ts

async function waveRectangle(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox } = await labelHelper(parent, node, getNodeClasses(node));
  const minWidth = 100;
  const minHeight = 50;
  const baseWidth = Math.max(bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const baseHeight = Math.max(bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const aspectRatio = baseWidth / baseHeight;
  let w = baseWidth;
  let h = baseHeight;
  if (w > h * aspectRatio) {
    h = w / aspectRatio;
  } else {
    w = h * aspectRatio;
  }
  w = Math.max(w, minWidth);
  h = Math.max(h, minHeight);
  const waveAmplitude = Math.min(h * 0.2, h / 4);
  const finalH = h + waveAmplitude * 2;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const points = [
    { x: -w / 2, y: finalH / 2 },
    ...generateFullSineWavePoints(-w / 2, finalH / 2, w / 2, finalH / 2, waveAmplitude, 1),
    { x: w / 2, y: -finalH / 2 },
    ...generateFullSineWavePoints(w / 2, -finalH / 2, -w / 2, -finalH / 2, waveAmplitude, -1)
  ];
  const waveRectPath = createPathFromPoints(points);
  const waveRectNode = rc.path(waveRectPath, options);
  const waveRect = shapeSvg.insert(() => waveRectNode, ":first-child");
  waveRect.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    waveRect.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    waveRect.selectAll("path").attr("style", nodeStyles);
  }
  updateNodeBounds(node, waveRect);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, points, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(waveRectangle, "waveRectangle");

// src/rendering-util/rendering-elements/shapes/windowPane.ts

async function windowPane(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, label } = await labelHelper(parent, node, getNodeClasses(node));
  const w = Math.max(bbox.width + (node.padding ?? 0) * 2, node?.width ?? 0);
  const h = Math.max(bbox.height + (node.padding ?? 0) * 2, node?.height ?? 0);
  const rectOffset = 5;
  const x = -w / 2;
  const y = -h / 2;
  const { cssStyles } = node;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  const outerPathPoints = [
    { x: x - rectOffset, y: y - rectOffset },
    { x: x - rectOffset, y: y + h },
    { x: x + w, y: y + h },
    { x: x + w, y: y - rectOffset }
  ];
  const path = `M${x - rectOffset},${y - rectOffset} L${x + w},${y - rectOffset} L${x + w},${y + h} L${x - rectOffset},${y + h} L${x - rectOffset},${y - rectOffset}
                M${x - rectOffset},${y} L${x + w},${y}
                M${x},${y - rectOffset} L${x},${y + h}`;
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const no = rc.path(path, options);
  const windowPane2 = shapeSvg.insert(() => no, ":first-child");
  windowPane2.attr("transform", `translate(${rectOffset / 2}, ${rectOffset / 2})`);
  windowPane2.attr("class", "basic label-container");
  if (cssStyles && node.look !== "handDrawn") {
    windowPane2.selectAll("path").attr("style", cssStyles);
  }
  if (nodeStyles && node.look !== "handDrawn") {
    windowPane2.selectAll("path").attr("style", nodeStyles);
  }
  label.attr(
    "transform",
    `translate(${-(bbox.width / 2) + rectOffset / 2 - (bbox.x - (bbox.left ?? 0))}, ${-(bbox.height / 2) + rectOffset / 2 - (bbox.y - (bbox.top ?? 0))})`
  );
  updateNodeBounds(node, windowPane2);
  node.intersect = function(point) {
    const pos = intersect_default.polygon(node, outerPathPoints, point);
    return pos;
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(windowPane, "windowPane");

// src/rendering-util/rendering-elements/shapes/erBox.ts


async function erBox(parent, node) {
  const entityNode = node;
  if (entityNode.alias) {
    node.label = entityNode.alias;
  }
  if (node.look === "handDrawn") {
    const { themeVariables: themeVariables2 } = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig */ .iE)();
    const { background } = themeVariables2;
    const backgroundNode = {
      ...node,
      id: node.id + "-background",
      look: "default",
      cssStyles: ["stroke: none", `fill: ${background}`]
    };
    await erBox(parent, backgroundNode);
  }
  const config = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig */ .iE)();
  node.useHtmlLabels = config.htmlLabels;
  let PADDING = config.er?.diagramPadding ?? 10;
  let TEXT_PADDING = config.er?.entityPadding ?? 6;
  const { cssStyles } = node;
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  if (entityNode.attributes.length === 0 && node.label) {
    const options2 = {
      rx: 0,
      ry: 0,
      labelPaddingX: PADDING,
      labelPaddingY: PADDING * 1.5,
      classes: ""
    };
    if ((0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .calculateTextWidth */ .Cq)(node.label, config) + options2.labelPaddingX * 2 < config.er.minEntityWidth) {
      node.width = config.er.minEntityWidth;
    }
    const shapeSvg2 = await drawRect(parent, node, options2);
    if (!(0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)(config.htmlLabels)) {
      const textElement = shapeSvg2.select("text");
      const bbox = textElement.node()?.getBBox();
      textElement.attr("transform", `translate(${-bbox.width / 2}, 0)`);
    }
    return shapeSvg2;
  }
  if (!config.htmlLabels) {
    PADDING *= 1.25;
    TEXT_PADDING *= 1.25;
  }
  let cssClasses = getNodeClasses(node);
  if (!cssClasses) {
    cssClasses = "node default";
  }
  const shapeSvg = parent.insert("g").attr("class", cssClasses).attr("id", node.domId || node.id);
  const nameBBox = await addText(shapeSvg, node.label ?? "", config, 0, 0, ["name"], labelStyles);
  nameBBox.height += TEXT_PADDING;
  let yOffset = 0;
  const yOffsets = [];
  const rows = [];
  let maxTypeWidth = 0;
  let maxNameWidth = 0;
  let maxKeysWidth = 0;
  let maxCommentWidth = 0;
  let keysPresent = true;
  let commentPresent = true;
  for (const attribute of entityNode.attributes) {
    const typeBBox = await addText(
      shapeSvg,
      attribute.type,
      config,
      0,
      yOffset,
      ["attribute-type"],
      labelStyles
    );
    maxTypeWidth = Math.max(maxTypeWidth, typeBBox.width + PADDING);
    const nameBBox2 = await addText(
      shapeSvg,
      attribute.name,
      config,
      0,
      yOffset,
      ["attribute-name"],
      labelStyles
    );
    maxNameWidth = Math.max(maxNameWidth, nameBBox2.width + PADDING);
    const keysBBox = await addText(
      shapeSvg,
      attribute.keys.join(),
      config,
      0,
      yOffset,
      ["attribute-keys"],
      labelStyles
    );
    maxKeysWidth = Math.max(maxKeysWidth, keysBBox.width + PADDING);
    const commentBBox = await addText(
      shapeSvg,
      attribute.comment,
      config,
      0,
      yOffset,
      ["attribute-comment"],
      labelStyles
    );
    maxCommentWidth = Math.max(maxCommentWidth, commentBBox.width + PADDING);
    const rowHeight = Math.max(typeBBox.height, nameBBox2.height, keysBBox.height, commentBBox.height) + TEXT_PADDING;
    rows.push({ yOffset, rowHeight });
    yOffset += rowHeight;
  }
  let totalWidthSections = 4;
  if (maxKeysWidth <= PADDING) {
    keysPresent = false;
    maxKeysWidth = 0;
    totalWidthSections--;
  }
  if (maxCommentWidth <= PADDING) {
    commentPresent = false;
    maxCommentWidth = 0;
    totalWidthSections--;
  }
  const shapeBBox = shapeSvg.node().getBBox();
  if (nameBBox.width + PADDING * 2 - (maxTypeWidth + maxNameWidth + maxKeysWidth + maxCommentWidth) > 0) {
    const difference = nameBBox.width + PADDING * 2 - (maxTypeWidth + maxNameWidth + maxKeysWidth + maxCommentWidth);
    maxTypeWidth += difference / totalWidthSections;
    maxNameWidth += difference / totalWidthSections;
    if (maxKeysWidth > 0) {
      maxKeysWidth += difference / totalWidthSections;
    }
    if (maxCommentWidth > 0) {
      maxCommentWidth += difference / totalWidthSections;
    }
  }
  const maxWidth = maxTypeWidth + maxNameWidth + maxKeysWidth + maxCommentWidth;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  let totalShapeBBoxHeight = 0;
  if (rows.length > 0) {
    totalShapeBBoxHeight = rows.reduce((sum, row) => sum + (row?.rowHeight ?? 0), 0);
  }
  const w = Math.max(shapeBBox.width + PADDING * 2, node?.width || 0, maxWidth);
  const h = Math.max((totalShapeBBoxHeight ?? 0) + nameBBox.height, node?.height || 0);
  const x = -w / 2;
  const y = -h / 2;
  shapeSvg.selectAll("g:not(:first-child)").each((_, i, nodes) => {
    const text2 = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(nodes[i]);
    const transform = text2.attr("transform");
    let translateX = 0;
    let translateY = 0;
    if (transform) {
      const regex = RegExp(/translate\(([^,]+),([^)]+)\)/);
      const translate = regex.exec(transform);
      if (translate) {
        translateX = parseFloat(translate[1]);
        translateY = parseFloat(translate[2]);
        if (text2.attr("class").includes("attribute-name")) {
          translateX += maxTypeWidth;
        } else if (text2.attr("class").includes("attribute-keys")) {
          translateX += maxTypeWidth + maxNameWidth;
        } else if (text2.attr("class").includes("attribute-comment")) {
          translateX += maxTypeWidth + maxNameWidth + maxKeysWidth;
        }
      }
    }
    text2.attr(
      "transform",
      `translate(${x + PADDING / 2 + translateX}, ${translateY + y + nameBBox.height + TEXT_PADDING / 2})`
    );
  });
  shapeSvg.select(".name").attr("transform", "translate(" + -nameBBox.width / 2 + ", " + (y + TEXT_PADDING / 2) + ")");
  const roughRect = rc.rectangle(x, y, w, h, options);
  const rect2 = shapeSvg.insert(() => roughRect, ":first-child").attr("style", cssStyles.join(""));
  const { themeVariables } = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig */ .iE)();
  const { rowEven, rowOdd, nodeBorder } = themeVariables;
  yOffsets.push(0);
  for (const [i, row] of rows.entries()) {
    const contentRowIndex = i + 1;
    const isEven = contentRowIndex % 2 === 0 && row.yOffset !== 0;
    const roughRect2 = rc.rectangle(x, nameBBox.height + y + row?.yOffset, w, row?.rowHeight, {
      ...options,
      fill: isEven ? rowEven : rowOdd,
      stroke: nodeBorder
    });
    shapeSvg.insert(() => roughRect2, "g.label").attr("style", cssStyles.join("")).attr("class", `row-rect-${isEven ? "even" : "odd"}`);
  }
  let roughLine = rc.line(x, nameBBox.height + y, w + x, nameBBox.height + y, options);
  shapeSvg.insert(() => roughLine).attr("class", "divider");
  roughLine = rc.line(maxTypeWidth + x, nameBBox.height + y, maxTypeWidth + x, h + y, options);
  shapeSvg.insert(() => roughLine).attr("class", "divider");
  if (keysPresent) {
    roughLine = rc.line(
      maxTypeWidth + maxNameWidth + x,
      nameBBox.height + y,
      maxTypeWidth + maxNameWidth + x,
      h + y,
      options
    );
    shapeSvg.insert(() => roughLine).attr("class", "divider");
  }
  if (commentPresent) {
    roughLine = rc.line(
      maxTypeWidth + maxNameWidth + maxKeysWidth + x,
      nameBBox.height + y,
      maxTypeWidth + maxNameWidth + maxKeysWidth + x,
      h + y,
      options
    );
    shapeSvg.insert(() => roughLine).attr("class", "divider");
  }
  for (const yOffset2 of yOffsets) {
    roughLine = rc.line(
      x,
      nameBBox.height + y + yOffset2,
      w + x,
      nameBBox.height + y + yOffset2,
      options
    );
    shapeSvg.insert(() => roughLine).attr("class", "divider");
  }
  updateNodeBounds(node, rect2);
  if (nodeStyles && node.look !== "handDrawn") {
    const allStyle = nodeStyles.split(";");
    const strokeStyles = allStyle?.filter((e) => {
      return e.includes("stroke");
    })?.map((s) => `${s}`).join("; ");
    shapeSvg.selectAll("path").attr("style", strokeStyles ?? "");
    shapeSvg.selectAll(".row-rect-even path").attr("style", nodeStyles);
  }
  node.intersect = function(point) {
    return intersect_default.rect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(erBox, "erBox");
async function addText(shapeSvg, labelText, config, translateX = 0, translateY = 0, classes = [], style = "") {
  const label = shapeSvg.insert("g").attr("class", `label ${classes.join(" ")}`).attr("transform", `translate(${translateX}, ${translateY})`).attr("style", style);
  if (labelText !== (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .parseGenericTypes */ .UO)(labelText)) {
    labelText = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .parseGenericTypes */ .UO)(labelText);
    labelText = labelText.replaceAll("<", "&lt;").replaceAll(">", "&gt;");
  }
  const text2 = label.node().appendChild(
    await (0,_chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_2__/* .createText */ .rw)(
      label,
      labelText,
      {
        width: (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .calculateTextWidth */ .Cq)(labelText, config) + 100,
        style,
        useHtmlLabels: config.htmlLabels
      },
      config
    )
  );
  if (labelText.includes("&lt;") || labelText.includes("&gt;")) {
    let child = text2.children[0];
    child.textContent = child.textContent.replaceAll("&lt;", "<").replaceAll("&gt;", ">");
    while (child.childNodes[0]) {
      child = child.childNodes[0];
      child.textContent = child.textContent.replaceAll("&lt;", "<").replaceAll("&gt;", ">");
    }
  }
  let bbox = text2.getBBox();
  if ((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)(config.htmlLabels)) {
    const div = text2.children[0];
    div.style.textAlign = "start";
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(text2);
    bbox = div.getBoundingClientRect();
    dv.attr("width", bbox.width);
    dv.attr("height", bbox.height);
  }
  return bbox;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(addText, "addText");

// src/rendering-util/rendering-elements/shapes/classBox.ts



// src/diagrams/class/shapeUtil.ts

async function textHelper(parent, node, config, useHtmlLabels, GAP = config.class.padding ?? 12) {
  const TEXT_PADDING = !useHtmlLabels ? 3 : 0;
  const shapeSvg = parent.insert("g").attr("class", getNodeClasses(node)).attr("id", node.domId || node.id);
  let annotationGroup = null;
  let labelGroup = null;
  let membersGroup = null;
  let methodsGroup = null;
  let annotationGroupHeight = 0;
  let labelGroupHeight = 0;
  let membersGroupHeight = 0;
  annotationGroup = shapeSvg.insert("g").attr("class", "annotation-group text");
  if (node.annotations.length > 0) {
    const annotation = node.annotations[0];
    await addText2(annotationGroup, { text: `\xAB${annotation}\xBB` }, 0);
    const annotationGroupBBox = annotationGroup.node().getBBox();
    annotationGroupHeight = annotationGroupBBox.height;
  }
  labelGroup = shapeSvg.insert("g").attr("class", "label-group text");
  await addText2(labelGroup, node, 0, ["font-weight: bolder"]);
  const labelGroupBBox = labelGroup.node().getBBox();
  labelGroupHeight = labelGroupBBox.height;
  membersGroup = shapeSvg.insert("g").attr("class", "members-group text");
  let yOffset = 0;
  for (const member of node.members) {
    const height = await addText2(membersGroup, member, yOffset, [member.parseClassifier()]);
    yOffset += height + TEXT_PADDING;
  }
  membersGroupHeight = membersGroup.node().getBBox().height;
  if (membersGroupHeight <= 0) {
    membersGroupHeight = GAP / 2;
  }
  methodsGroup = shapeSvg.insert("g").attr("class", "methods-group text");
  let methodsYOffset = 0;
  for (const method of node.methods) {
    const height = await addText2(methodsGroup, method, methodsYOffset, [method.parseClassifier()]);
    methodsYOffset += height + TEXT_PADDING;
  }
  let bbox = shapeSvg.node().getBBox();
  if (annotationGroup !== null) {
    const annotationGroupBBox = annotationGroup.node().getBBox();
    annotationGroup.attr("transform", `translate(${-annotationGroupBBox.width / 2})`);
  }
  labelGroup.attr("transform", `translate(${-labelGroupBBox.width / 2}, ${annotationGroupHeight})`);
  bbox = shapeSvg.node().getBBox();
  membersGroup.attr(
    "transform",
    `translate(${0}, ${annotationGroupHeight + labelGroupHeight + GAP * 2})`
  );
  bbox = shapeSvg.node().getBBox();
  methodsGroup.attr(
    "transform",
    `translate(${0}, ${annotationGroupHeight + labelGroupHeight + (membersGroupHeight ? membersGroupHeight + GAP * 4 : GAP * 2)})`
  );
  bbox = shapeSvg.node().getBBox();
  return { shapeSvg, bbox };
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(textHelper, "textHelper");
async function addText2(parentGroup, node, yOffset, styles = []) {
  const textEl = parentGroup.insert("g").attr("class", "label").attr("style", styles.join("; "));
  const config = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig */ .iE)();
  let useHtmlLabels = "useHtmlLabels" in node ? node.useHtmlLabels : (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)(config.htmlLabels) ?? true;
  let textContent = "";
  if ("text" in node) {
    textContent = node.text;
  } else {
    textContent = node.label;
  }
  if (!useHtmlLabels && textContent.startsWith("\\")) {
    textContent = textContent.substring(1);
  }
  if ((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .hasKatex */ .l0)(textContent)) {
    useHtmlLabels = true;
  }
  const text2 = await (0,_chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_2__/* .createText */ .rw)(
    textEl,
    (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .sanitizeText2 */ .uX)((0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .decodeEntities */ .SH)(textContent)),
    {
      width: (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .calculateTextWidth */ .Cq)(textContent, config) + 50,
      // Add room for error when splitting text into multiple lines
      classes: "markdown-node-label",
      useHtmlLabels
    },
    config
  );
  let bbox;
  let numberOfLines = 1;
  if (!useHtmlLabels) {
    if (styles.includes("font-weight: bolder")) {
      (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(text2).selectAll("tspan").attr("font-weight", "");
    }
    numberOfLines = text2.children.length;
    const textChild = text2.children[0];
    if (text2.textContent === "" || text2.textContent.includes("&gt")) {
      textChild.textContent = textContent[0] + textContent.substring(1).replaceAll("&gt;", ">").replaceAll("&lt;", "<").trim();
      const preserveSpace = textContent[1] === " ";
      if (preserveSpace) {
        textChild.textContent = textChild.textContent[0] + " " + textChild.textContent.substring(1);
      }
    }
    if (textChild.textContent === "undefined") {
      textChild.textContent = "";
    }
    bbox = text2.getBBox();
  } else {
    const div = text2.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(text2);
    numberOfLines = div.innerHTML.split("<br>").length;
    if (div.innerHTML.includes("</math>")) {
      numberOfLines += div.innerHTML.split("<mrow>").length - 1;
    }
    const images = div.getElementsByTagName("img");
    if (images) {
      const noImgText = textContent.replace(/<img[^>]*>/g, "").trim() === "";
      await Promise.all(
        [...images].map(
          (img) => new Promise((res) => {
            function setupImage() {
              img.style.display = "flex";
              img.style.flexDirection = "column";
              if (noImgText) {
                const bodyFontSize = config.fontSize?.toString() ?? window.getComputedStyle(document.body).fontSize;
                const enlargingFactor = 5;
                const width = parseInt(bodyFontSize, 10) * enlargingFactor + "px";
                img.style.minWidth = width;
                img.style.maxWidth = width;
              } else {
                img.style.width = "100%";
              }
              res(img);
            }
            (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(setupImage, "setupImage");
            setTimeout(() => {
              if (img.complete) {
                setupImage();
              }
            });
            img.addEventListener("error", setupImage);
            img.addEventListener("load", setupImage);
          })
        )
      );
    }
    bbox = div.getBoundingClientRect();
    dv.attr("width", bbox.width);
    dv.attr("height", bbox.height);
  }
  textEl.attr("transform", "translate(0," + (-bbox.height / (2 * numberOfLines) + yOffset) + ")");
  return bbox.height;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(addText2, "addText");

// src/rendering-util/rendering-elements/shapes/classBox.ts
async function classBox(parent, node) {
  const config = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)();
  const PADDING = config.class.padding ?? 12;
  const GAP = PADDING;
  const useHtmlLabels = node.useHtmlLabels ?? (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .evaluate */ .ku)(config.htmlLabels) ?? true;
  const classNode = node;
  classNode.annotations = classNode.annotations ?? [];
  classNode.members = classNode.members ?? [];
  classNode.methods = classNode.methods ?? [];
  const { shapeSvg, bbox } = await textHelper(parent, node, config, useHtmlLabels, GAP);
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  node.cssStyles = classNode.styles || "";
  const styles = classNode.styles?.join(";") || nodeStyles || "";
  if (!node.cssStyles) {
    node.cssStyles = styles.replaceAll("!important", "").split(";");
  }
  const renderExtraBox = classNode.members.length === 0 && classNode.methods.length === 0 && !config.class?.hideEmptyMembersBox;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const w = bbox.width;
  let h = bbox.height;
  if (classNode.members.length === 0 && classNode.methods.length === 0) {
    h += GAP;
  } else if (classNode.members.length > 0 && classNode.methods.length === 0) {
    h += GAP * 2;
  }
  const x = -w / 2;
  const y = -h / 2;
  const roughRect = rc.rectangle(
    x - PADDING,
    y - PADDING - (renderExtraBox ? PADDING : classNode.members.length === 0 && classNode.methods.length === 0 ? -PADDING / 2 : 0),
    w + 2 * PADDING,
    h + 2 * PADDING + (renderExtraBox ? PADDING * 2 : classNode.members.length === 0 && classNode.methods.length === 0 ? -PADDING : 0),
    options
  );
  const rect2 = shapeSvg.insert(() => roughRect, ":first-child");
  rect2.attr("class", "basic label-container");
  const rectBBox = rect2.node().getBBox();
  shapeSvg.selectAll(".text").each((_, i, nodes) => {
    const text2 = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(nodes[i]);
    const transform = text2.attr("transform");
    let translateY = 0;
    if (transform) {
      const regex = RegExp(/translate\(([^,]+),([^)]+)\)/);
      const translate = regex.exec(transform);
      if (translate) {
        translateY = parseFloat(translate[2]);
      }
    }
    let newTranslateY = translateY + y + PADDING - (renderExtraBox ? PADDING : classNode.members.length === 0 && classNode.methods.length === 0 ? -PADDING / 2 : 0);
    if (!useHtmlLabels) {
      newTranslateY -= 4;
    }
    let newTranslateX = x;
    if (text2.attr("class").includes("label-group") || text2.attr("class").includes("annotation-group")) {
      newTranslateX = -text2.node()?.getBBox().width / 2 || 0;
      shapeSvg.selectAll("text").each(function(_2, i2, nodes2) {
        if (window.getComputedStyle(nodes2[i2]).textAnchor === "middle") {
          newTranslateX = 0;
        }
      });
    }
    text2.attr("transform", `translate(${newTranslateX}, ${newTranslateY})`);
  });
  const annotationGroupHeight = shapeSvg.select(".annotation-group").node().getBBox().height - (renderExtraBox ? PADDING / 2 : 0) || 0;
  const labelGroupHeight = shapeSvg.select(".label-group").node().getBBox().height - (renderExtraBox ? PADDING / 2 : 0) || 0;
  const membersGroupHeight = shapeSvg.select(".members-group").node().getBBox().height - (renderExtraBox ? PADDING / 2 : 0) || 0;
  if (classNode.members.length > 0 || classNode.methods.length > 0 || renderExtraBox) {
    const roughLine = rc.line(
      rectBBox.x,
      annotationGroupHeight + labelGroupHeight + y + PADDING,
      rectBBox.x + rectBBox.width,
      annotationGroupHeight + labelGroupHeight + y + PADDING,
      options
    );
    const line = shapeSvg.insert(() => roughLine);
    line.attr("class", "divider").attr("style", styles);
  }
  if (renderExtraBox || classNode.members.length > 0 || classNode.methods.length > 0) {
    const roughLine = rc.line(
      rectBBox.x,
      annotationGroupHeight + labelGroupHeight + membersGroupHeight + y + GAP * 2 + PADDING,
      rectBBox.x + rectBBox.width,
      annotationGroupHeight + labelGroupHeight + membersGroupHeight + y + PADDING + GAP * 2,
      options
    );
    const line = shapeSvg.insert(() => roughLine);
    line.attr("class", "divider").attr("style", styles);
  }
  if (classNode.look !== "handDrawn") {
    shapeSvg.selectAll("path").attr("style", styles);
  }
  rect2.select(":nth-child(2)").attr("style", styles);
  shapeSvg.selectAll(".divider").select("path").attr("style", styles);
  if (node.labelStyle) {
    shapeSvg.selectAll("span").attr("style", node.labelStyle);
  } else {
    shapeSvg.selectAll("span").attr("style", styles);
  }
  if (!useHtmlLabels) {
    const colorRegex = RegExp(/color\s*:\s*([^;]*)/);
    const match = colorRegex.exec(styles);
    if (match) {
      const colorStyle = match[0].replace("color", "fill");
      shapeSvg.selectAll("tspan").attr("style", colorStyle);
    } else if (labelStyles) {
      const match2 = colorRegex.exec(labelStyles);
      if (match2) {
        const colorStyle = match2[0].replace("color", "fill");
        shapeSvg.selectAll("tspan").attr("style", colorStyle);
      }
    }
  }
  updateNodeBounds(node, rect2);
  node.intersect = function(point) {
    return intersect_default.rect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(classBox, "classBox");

// src/rendering-util/rendering-elements/shapes/requirementBox.ts


async function requirementBox(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const requirementNode = node;
  const elementNode = node;
  const padding = 20;
  const gap = 20;
  const isRequirementNode = "verifyMethod" in node;
  const classes = getNodeClasses(node);
  const shapeSvg = parent.insert("g").attr("class", classes).attr("id", node.domId ?? node.id);
  let typeHeight;
  if (isRequirementNode) {
    typeHeight = await addText3(
      shapeSvg,
      `&lt;&lt;${requirementNode.type}&gt;&gt;`,
      0,
      node.labelStyle
    );
  } else {
    typeHeight = await addText3(shapeSvg, "&lt;&lt;Element&gt;&gt;", 0, node.labelStyle);
  }
  let accumulativeHeight = typeHeight;
  const nameHeight = await addText3(
    shapeSvg,
    requirementNode.name,
    accumulativeHeight,
    node.labelStyle + "; font-weight: bold;"
  );
  accumulativeHeight += nameHeight + gap;
  if (isRequirementNode) {
    const idHeight = await addText3(
      shapeSvg,
      `${requirementNode.requirementId ? `ID: ${requirementNode.requirementId}` : ""}`,
      accumulativeHeight,
      node.labelStyle
    );
    accumulativeHeight += idHeight;
    const textHeight = await addText3(
      shapeSvg,
      `${requirementNode.text ? `Text: ${requirementNode.text}` : ""}`,
      accumulativeHeight,
      node.labelStyle
    );
    accumulativeHeight += textHeight;
    const riskHeight = await addText3(
      shapeSvg,
      `${requirementNode.risk ? `Risk: ${requirementNode.risk}` : ""}`,
      accumulativeHeight,
      node.labelStyle
    );
    accumulativeHeight += riskHeight;
    await addText3(
      shapeSvg,
      `${requirementNode.verifyMethod ? `Verification: ${requirementNode.verifyMethod}` : ""}`,
      accumulativeHeight,
      node.labelStyle
    );
  } else {
    const typeHeight2 = await addText3(
      shapeSvg,
      `${elementNode.type ? `Type: ${elementNode.type}` : ""}`,
      accumulativeHeight,
      node.labelStyle
    );
    accumulativeHeight += typeHeight2;
    await addText3(
      shapeSvg,
      `${elementNode.docRef ? `Doc Ref: ${elementNode.docRef}` : ""}`,
      accumulativeHeight,
      node.labelStyle
    );
  }
  const totalWidth = (shapeSvg.node()?.getBBox().width ?? 200) + padding;
  const totalHeight = (shapeSvg.node()?.getBBox().height ?? 200) + padding;
  const x = -totalWidth / 2;
  const y = -totalHeight / 2;
  const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
  const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
  if (node.look !== "handDrawn") {
    options.roughness = 0;
    options.fillStyle = "solid";
  }
  const roughRect = rc.rectangle(x, y, totalWidth, totalHeight, options);
  const rect2 = shapeSvg.insert(() => roughRect, ":first-child");
  rect2.attr("class", "basic label-container").attr("style", nodeStyles);
  shapeSvg.selectAll(".label").each((_, i, nodes) => {
    const text2 = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(nodes[i]);
    const transform = text2.attr("transform");
    let translateX = 0;
    let translateY = 0;
    if (transform) {
      const regex = RegExp(/translate\(([^,]+),([^)]+)\)/);
      const translate = regex.exec(transform);
      if (translate) {
        translateX = parseFloat(translate[1]);
        translateY = parseFloat(translate[2]);
      }
    }
    const newTranslateY = translateY - totalHeight / 2;
    let newTranslateX = x + padding / 2;
    if (i === 0 || i === 1) {
      newTranslateX = translateX;
    }
    text2.attr("transform", `translate(${newTranslateX}, ${newTranslateY + padding})`);
  });
  if (accumulativeHeight > typeHeight + nameHeight + gap) {
    const roughLine = rc.line(
      x,
      y + typeHeight + nameHeight + gap,
      x + totalWidth,
      y + typeHeight + nameHeight + gap,
      options
    );
    const dividerLine = shapeSvg.insert(() => roughLine);
    dividerLine.attr("style", nodeStyles);
  }
  updateNodeBounds(node, rect2);
  node.intersect = function(point) {
    return intersect_default.rect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(requirementBox, "requirementBox");
async function addText3(parentGroup, inputText, yOffset, style = "") {
  if (inputText === "") {
    return 0;
  }
  const textEl = parentGroup.insert("g").attr("class", "label").attr("style", style);
  const config = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .getConfig2 */ .nV)();
  const useHtmlLabels = config.htmlLabels ?? true;
  const text2 = await (0,_chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_2__/* .createText */ .rw)(
    textEl,
    (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_4__/* .sanitizeText2 */ .uX)((0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .decodeEntities */ .SH)(inputText)),
    {
      width: (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .calculateTextWidth */ .Cq)(inputText, config) + 50,
      // Add room for error when splitting text into multiple lines
      classes: "markdown-node-label",
      useHtmlLabels,
      style
    },
    config
  );
  let bbox;
  if (!useHtmlLabels) {
    const textChild = text2.children[0];
    for (const child of textChild.children) {
      child.textContent = child.textContent.replaceAll("&gt;", ">").replaceAll("&lt;", "<");
      if (style) {
        child.setAttribute("style", style);
      }
    }
    bbox = text2.getBBox();
    bbox.height += 6;
  } else {
    const div = text2.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_6__/* .select */ .Ys)(text2);
    bbox = div.getBoundingClientRect();
    dv.attr("width", bbox.width);
    dv.attr("height", bbox.height);
  }
  textEl.attr("transform", `translate(${-bbox.width / 2},${-bbox.height / 2 + yOffset})`);
  return bbox.height;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(addText3, "addText");

// src/rendering-util/rendering-elements/shapes/kanbanItem.ts

var colorFromPriority = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((priority) => {
  switch (priority) {
    case "Very High":
      return "red";
    case "High":
      return "orange";
    case "Medium":
      return null;
    // no stroke
    case "Low":
      return "blue";
    case "Very Low":
      return "lightblue";
  }
}, "colorFromPriority");
async function kanbanItem(parent, kanbanNode, { config }) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(kanbanNode);
  kanbanNode.labelStyle = labelStyles || "";
  const labelPaddingX = 10;
  const orgWidth = kanbanNode.width;
  kanbanNode.width = (kanbanNode.width ?? 200) - 10;
  const {
    shapeSvg,
    bbox,
    label: labelElTitle
  } = await labelHelper(parent, kanbanNode, getNodeClasses(kanbanNode));
  const padding = kanbanNode.padding || 10;
  let ticketUrl = "";
  let link;
  if ("ticket" in kanbanNode && kanbanNode.ticket && config?.kanban?.ticketBaseUrl) {
    ticketUrl = config?.kanban?.ticketBaseUrl.replace("#TICKET#", kanbanNode.ticket);
    link = shapeSvg.insert("svg:a", ":first-child").attr("class", "kanban-ticket-link").attr("xlink:href", ticketUrl).attr("target", "_blank");
  }
  const options = {
    useHtmlLabels: kanbanNode.useHtmlLabels,
    labelStyle: kanbanNode.labelStyle || "",
    width: kanbanNode.width,
    img: kanbanNode.img,
    padding: kanbanNode.padding || 8,
    centerLabel: false
  };
  let labelEl, bbox2;
  if (link) {
    ({ label: labelEl, bbox: bbox2 } = await insertLabel(
      link,
      "ticket" in kanbanNode && kanbanNode.ticket || "",
      options
    ));
  } else {
    ({ label: labelEl, bbox: bbox2 } = await insertLabel(
      shapeSvg,
      "ticket" in kanbanNode && kanbanNode.ticket || "",
      options
    ));
  }
  const { label: labelElAssigned, bbox: bboxAssigned } = await insertLabel(
    shapeSvg,
    "assigned" in kanbanNode && kanbanNode.assigned || "",
    options
  );
  kanbanNode.width = orgWidth;
  const labelPaddingY = 10;
  const totalWidth = kanbanNode?.width || 0;
  const heightAdj = Math.max(bbox2.height, bboxAssigned.height) / 2;
  const totalHeight = Math.max(bbox.height + labelPaddingY * 2, kanbanNode?.height || 0) + heightAdj;
  const x = -totalWidth / 2;
  const y = -totalHeight / 2;
  labelElTitle.attr(
    "transform",
    "translate(" + (padding - totalWidth / 2) + ", " + (-heightAdj - bbox.height / 2) + ")"
  );
  labelEl.attr(
    "transform",
    "translate(" + (padding - totalWidth / 2) + ", " + (-heightAdj + bbox.height / 2) + ")"
  );
  labelElAssigned.attr(
    "transform",
    "translate(" + (padding + totalWidth / 2 - bboxAssigned.width - 2 * labelPaddingX) + ", " + (-heightAdj + bbox.height / 2) + ")"
  );
  let rect2;
  const { rx, ry } = kanbanNode;
  const { cssStyles } = kanbanNode;
  if (kanbanNode.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options2 = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(kanbanNode, {});
    const roughNode = rx || ry ? rc.path(createRoundedRectPathD(x, y, totalWidth, totalHeight, rx || 0), options2) : rc.rectangle(x, y, totalWidth, totalHeight, options2);
    rect2 = shapeSvg.insert(() => roughNode, ":first-child");
    rect2.attr("class", "basic label-container").attr("style", cssStyles ? cssStyles : null);
  } else {
    rect2 = shapeSvg.insert("rect", ":first-child");
    rect2.attr("class", "basic label-container __APA__").attr("style", nodeStyles).attr("rx", rx ?? 5).attr("ry", ry ?? 5).attr("x", x).attr("y", y).attr("width", totalWidth).attr("height", totalHeight);
    const priority = "priority" in kanbanNode && kanbanNode.priority;
    if (priority) {
      const line = shapeSvg.append("line");
      const lineX = x + 2;
      const y1 = y + Math.floor((rx ?? 0) / 2);
      const y2 = y + totalHeight - Math.floor((rx ?? 0) / 2);
      line.attr("x1", lineX).attr("y1", y1).attr("x2", lineX).attr("y2", y2).attr("stroke-width", "4").attr("stroke", colorFromPriority(priority));
    }
  }
  updateNodeBounds(kanbanNode, rect2);
  kanbanNode.height = totalHeight;
  kanbanNode.intersect = function(point) {
    return intersect_default.rect(kanbanNode, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(kanbanItem, "kanbanItem");

// src/rendering-util/rendering-elements/shapes/bang.ts

async function bang(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, halfPadding, label } = await labelHelper(
    parent,
    node,
    getNodeClasses(node)
  );
  const w = bbox.width + 10 * halfPadding;
  const h = bbox.height + 8 * halfPadding;
  const r = 0.15 * w;
  const { cssStyles } = node;
  const minWidth = bbox.width + 20;
  const minHeight = bbox.height + 20;
  const effectiveWidth = Math.max(w, minWidth);
  const effectiveHeight = Math.max(h, minHeight);
  label.attr("transform", `translate(${-bbox.width / 2}, ${-bbox.height / 2})`);
  let bangElem;
  const path = `M0 0 
    a${r},${r} 1 0,0 ${effectiveWidth * 0.25},${-1 * effectiveHeight * 0.1}
    a${r},${r} 1 0,0 ${effectiveWidth * 0.25},${0}
    a${r},${r} 1 0,0 ${effectiveWidth * 0.25},${0}
    a${r},${r} 1 0,0 ${effectiveWidth * 0.25},${effectiveHeight * 0.1}

    a${r},${r} 1 0,0 ${effectiveWidth * 0.15},${effectiveHeight * 0.33}
    a${r * 0.8},${r * 0.8} 1 0,0 0,${effectiveHeight * 0.34}
    a${r},${r} 1 0,0 ${-1 * effectiveWidth * 0.15},${effectiveHeight * 0.33}

    a${r},${r} 1 0,0 ${-1 * effectiveWidth * 0.25},${effectiveHeight * 0.15}
    a${r},${r} 1 0,0 ${-1 * effectiveWidth * 0.25},0
    a${r},${r} 1 0,0 ${-1 * effectiveWidth * 0.25},0
    a${r},${r} 1 0,0 ${-1 * effectiveWidth * 0.25},${-1 * effectiveHeight * 0.15}

    a${r},${r} 1 0,0 ${-1 * effectiveWidth * 0.1},${-1 * effectiveHeight * 0.33}
    a${r * 0.8},${r * 0.8} 1 0,0 0,${-1 * effectiveHeight * 0.34}
    a${r},${r} 1 0,0 ${effectiveWidth * 0.1},${-1 * effectiveHeight * 0.33}
  H0 V0 Z`;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
    const roughNode = rc.path(path, options);
    bangElem = shapeSvg.insert(() => roughNode, ":first-child");
    bangElem.attr("class", "basic label-container").attr("style", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(cssStyles));
  } else {
    bangElem = shapeSvg.insert("path", ":first-child").attr("class", "basic label-container").attr("style", nodeStyles).attr("d", path);
  }
  bangElem.attr("transform", `translate(${-effectiveWidth / 2}, ${-effectiveHeight / 2})`);
  updateNodeBounds(node, bangElem);
  node.calcIntersect = function(bounds, point) {
    return intersect_default.rect(bounds, point);
  };
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("Bang intersect", node, point);
    return intersect_default.rect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(bang, "bang");

// src/rendering-util/rendering-elements/shapes/cloud.ts

async function cloud(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, halfPadding, label } = await labelHelper(
    parent,
    node,
    getNodeClasses(node)
  );
  const w = bbox.width + 2 * halfPadding;
  const h = bbox.height + 2 * halfPadding;
  const r1 = 0.15 * w;
  const r2 = 0.25 * w;
  const r3 = 0.35 * w;
  const r4 = 0.2 * w;
  const { cssStyles } = node;
  let cloudElem;
  const path = `M0 0 
    a${r1},${r1} 0 0,1 ${w * 0.25},${-1 * w * 0.1}
    a${r3},${r3} 1 0,1 ${w * 0.4},${-1 * w * 0.1}
    a${r2},${r2} 1 0,1 ${w * 0.35},${w * 0.2}

    a${r1},${r1} 1 0,1 ${w * 0.15},${h * 0.35}
    a${r4},${r4} 1 0,1 ${-1 * w * 0.15},${h * 0.65}

    a${r2},${r1} 1 0,1 ${-1 * w * 0.25},${w * 0.15}
    a${r3},${r3} 1 0,1 ${-1 * w * 0.5},0
    a${r1},${r1} 1 0,1 ${-1 * w * 0.25},${-1 * w * 0.15}

    a${r1},${r1} 1 0,1 ${-1 * w * 0.1},${-1 * h * 0.35}
    a${r4},${r4} 1 0,1 ${w * 0.1},${-1 * h * 0.65}
  H0 V0 Z`;
  if (node.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_7__/* ["default"] */ .Z.svg(shapeSvg);
    const options = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .userNodeOverrides */ ._q)(node, {});
    const roughNode = rc.path(path, options);
    cloudElem = shapeSvg.insert(() => roughNode, ":first-child");
    cloudElem.attr("class", "basic label-container").attr("style", (0,_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_3__/* .handleUndefinedAttr */ .R7)(cssStyles));
  } else {
    cloudElem = shapeSvg.insert("path", ":first-child").attr("class", "basic label-container").attr("style", nodeStyles).attr("d", path);
  }
  label.attr("transform", `translate(${-bbox.width / 2}, ${-bbox.height / 2})`);
  cloudElem.attr("transform", `translate(${-w / 2}, ${-h / 2})`);
  updateNodeBounds(node, cloudElem);
  node.calcIntersect = function(bounds, point) {
    return intersect_default.rect(bounds, point);
  };
  node.intersect = function(point) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.info("Cloud intersect", node, point);
    return intersect_default.rect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(cloud, "cloud");

// src/rendering-util/rendering-elements/shapes/defaultMindmapNode.ts
async function defaultMindmapNode(parent, node) {
  const { labelStyles, nodeStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_1__/* .styles2String */ .UG)(node);
  node.labelStyle = labelStyles;
  const { shapeSvg, bbox, halfPadding, label } = await labelHelper(
    parent,
    node,
    getNodeClasses(node)
  );
  const w = bbox.width + 8 * halfPadding;
  const h = bbox.height + 2 * halfPadding;
  const rd = 5;
  const rectPath = `
    M${-w / 2} ${h / 2 - rd}
    v${-h + 2 * rd}
    q0,-${rd} ${rd},-${rd}
    h${w - 2 * rd}
    q${rd},0 ${rd},${rd}
    v${h - 2 * rd}
    q0,${rd} -${rd},${rd}
    h${-w + 2 * rd}
    q-${rd},0 -${rd},-${rd}
    Z
  `;
  const bg = shapeSvg.append("path").attr("id", "node-" + node.id).attr("class", "node-bkg node-" + node.type).attr("style", nodeStyles).attr("d", rectPath);
  shapeSvg.append("line").attr("class", "node-line-").attr("x1", -w / 2).attr("y1", h / 2).attr("x2", w / 2).attr("y2", h / 2);
  label.attr("transform", `translate(${-bbox.width / 2}, ${-bbox.height / 2})`);
  shapeSvg.append(() => label.node());
  updateNodeBounds(node, bg);
  node.calcIntersect = function(bounds, point) {
    return intersect_default.rect(bounds, point);
  };
  node.intersect = function(point) {
    return intersect_default.rect(node, point);
  };
  return shapeSvg;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(defaultMindmapNode, "defaultMindmapNode");

// src/rendering-util/rendering-elements/shapes/mindmapCircle.ts
async function mindmapCircle(parent, node) {
  const options = {
    padding: node.padding ?? 0
  };
  return circle(parent, node, options);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(mindmapCircle, "mindmapCircle");

// src/rendering-util/rendering-elements/shapes.ts
var shapesDefs = [
  {
    semanticName: "Process",
    name: "Rectangle",
    shortName: "rect",
    description: "Standard process shape",
    aliases: ["proc", "process", "rectangle"],
    internalAliases: ["squareRect"],
    handler: squareRect2
  },
  {
    semanticName: "Event",
    name: "Rounded Rectangle",
    shortName: "rounded",
    description: "Represents an event",
    aliases: ["event"],
    internalAliases: ["roundedRect"],
    handler: roundedRect
  },
  {
    semanticName: "Terminal Point",
    name: "Stadium",
    shortName: "stadium",
    description: "Terminal point",
    aliases: ["terminal", "pill"],
    handler: stadium
  },
  {
    semanticName: "Subprocess",
    name: "Framed Rectangle",
    shortName: "fr-rect",
    description: "Subprocess",
    aliases: ["subprocess", "subproc", "framed-rectangle", "subroutine"],
    handler: subroutine
  },
  {
    semanticName: "Database",
    name: "Cylinder",
    shortName: "cyl",
    description: "Database storage",
    aliases: ["db", "database", "cylinder"],
    handler: cylinder
  },
  {
    semanticName: "Start",
    name: "Circle",
    shortName: "circle",
    description: "Starting point",
    aliases: ["circ"],
    handler: circle
  },
  {
    semanticName: "Bang",
    name: "Bang",
    shortName: "bang",
    description: "Bang",
    aliases: ["bang"],
    handler: bang
  },
  {
    semanticName: "Cloud",
    name: "Cloud",
    shortName: "cloud",
    description: "cloud",
    aliases: ["cloud"],
    handler: cloud
  },
  {
    semanticName: "Decision",
    name: "Diamond",
    shortName: "diam",
    description: "Decision-making step",
    aliases: ["decision", "diamond", "question"],
    handler: question
  },
  {
    semanticName: "Prepare Conditional",
    name: "Hexagon",
    shortName: "hex",
    description: "Preparation or condition step",
    aliases: ["hexagon", "prepare"],
    handler: hexagon
  },
  {
    semanticName: "Data Input/Output",
    name: "Lean Right",
    shortName: "lean-r",
    description: "Represents input or output",
    aliases: ["lean-right", "in-out"],
    internalAliases: ["lean_right"],
    handler: lean_right
  },
  {
    semanticName: "Data Input/Output",
    name: "Lean Left",
    shortName: "lean-l",
    description: "Represents output or input",
    aliases: ["lean-left", "out-in"],
    internalAliases: ["lean_left"],
    handler: lean_left
  },
  {
    semanticName: "Priority Action",
    name: "Trapezoid Base Bottom",
    shortName: "trap-b",
    description: "Priority action",
    aliases: ["priority", "trapezoid-bottom", "trapezoid"],
    handler: trapezoid
  },
  {
    semanticName: "Manual Operation",
    name: "Trapezoid Base Top",
    shortName: "trap-t",
    description: "Represents a manual task",
    aliases: ["manual", "trapezoid-top", "inv-trapezoid"],
    internalAliases: ["inv_trapezoid"],
    handler: inv_trapezoid
  },
  {
    semanticName: "Stop",
    name: "Double Circle",
    shortName: "dbl-circ",
    description: "Represents a stop point",
    aliases: ["double-circle"],
    internalAliases: ["doublecircle"],
    handler: doublecircle
  },
  {
    semanticName: "Text Block",
    name: "Text Block",
    shortName: "text",
    description: "Text block",
    handler: text
  },
  {
    semanticName: "Card",
    name: "Notched Rectangle",
    shortName: "notch-rect",
    description: "Represents a card",
    aliases: ["card", "notched-rectangle"],
    handler: card
  },
  {
    semanticName: "Lined/Shaded Process",
    name: "Lined Rectangle",
    shortName: "lin-rect",
    description: "Lined process shape",
    aliases: ["lined-rectangle", "lined-process", "lin-proc", "shaded-process"],
    handler: shadedProcess
  },
  {
    semanticName: "Start",
    name: "Small Circle",
    shortName: "sm-circ",
    description: "Small starting point",
    aliases: ["start", "small-circle"],
    internalAliases: ["stateStart"],
    handler: stateStart
  },
  {
    semanticName: "Stop",
    name: "Framed Circle",
    shortName: "fr-circ",
    description: "Stop point",
    aliases: ["stop", "framed-circle"],
    internalAliases: ["stateEnd"],
    handler: stateEnd
  },
  {
    semanticName: "Fork/Join",
    name: "Filled Rectangle",
    shortName: "fork",
    description: "Fork or join in process flow",
    aliases: ["join"],
    internalAliases: ["forkJoin"],
    handler: forkJoin
  },
  {
    semanticName: "Collate",
    name: "Hourglass",
    shortName: "hourglass",
    description: "Represents a collate operation",
    aliases: ["hourglass", "collate"],
    handler: hourglass
  },
  {
    semanticName: "Comment",
    name: "Curly Brace",
    shortName: "brace",
    description: "Adds a comment",
    aliases: ["comment", "brace-l"],
    handler: curlyBraceLeft
  },
  {
    semanticName: "Comment Right",
    name: "Curly Brace",
    shortName: "brace-r",
    description: "Adds a comment",
    handler: curlyBraceRight
  },
  {
    semanticName: "Comment with braces on both sides",
    name: "Curly Braces",
    shortName: "braces",
    description: "Adds a comment",
    handler: curlyBraces
  },
  {
    semanticName: "Com Link",
    name: "Lightning Bolt",
    shortName: "bolt",
    description: "Communication link",
    aliases: ["com-link", "lightning-bolt"],
    handler: lightningBolt
  },
  {
    semanticName: "Document",
    name: "Document",
    shortName: "doc",
    description: "Represents a document",
    aliases: ["doc", "document"],
    handler: waveEdgedRectangle
  },
  {
    semanticName: "Delay",
    name: "Half-Rounded Rectangle",
    shortName: "delay",
    description: "Represents a delay",
    aliases: ["half-rounded-rectangle"],
    handler: halfRoundedRectangle
  },
  {
    semanticName: "Direct Access Storage",
    name: "Horizontal Cylinder",
    shortName: "h-cyl",
    description: "Direct access storage",
    aliases: ["das", "horizontal-cylinder"],
    handler: tiltedCylinder
  },
  {
    semanticName: "Disk Storage",
    name: "Lined Cylinder",
    shortName: "lin-cyl",
    description: "Disk storage",
    aliases: ["disk", "lined-cylinder"],
    handler: linedCylinder
  },
  {
    semanticName: "Display",
    name: "Curved Trapezoid",
    shortName: "curv-trap",
    description: "Represents a display",
    aliases: ["curved-trapezoid", "display"],
    handler: curvedTrapezoid
  },
  {
    semanticName: "Divided Process",
    name: "Divided Rectangle",
    shortName: "div-rect",
    description: "Divided process shape",
    aliases: ["div-proc", "divided-rectangle", "divided-process"],
    handler: dividedRectangle
  },
  {
    semanticName: "Extract",
    name: "Triangle",
    shortName: "tri",
    description: "Extraction process",
    aliases: ["extract", "triangle"],
    handler: triangle
  },
  {
    semanticName: "Internal Storage",
    name: "Window Pane",
    shortName: "win-pane",
    description: "Internal storage",
    aliases: ["internal-storage", "window-pane"],
    handler: windowPane
  },
  {
    semanticName: "Junction",
    name: "Filled Circle",
    shortName: "f-circ",
    description: "Junction point",
    aliases: ["junction", "filled-circle"],
    handler: filledCircle
  },
  {
    semanticName: "Loop Limit",
    name: "Trapezoidal Pentagon",
    shortName: "notch-pent",
    description: "Loop limit step",
    aliases: ["loop-limit", "notched-pentagon"],
    handler: trapezoidalPentagon
  },
  {
    semanticName: "Manual File",
    name: "Flipped Triangle",
    shortName: "flip-tri",
    description: "Manual file operation",
    aliases: ["manual-file", "flipped-triangle"],
    handler: flippedTriangle
  },
  {
    semanticName: "Manual Input",
    name: "Sloped Rectangle",
    shortName: "sl-rect",
    description: "Manual input step",
    aliases: ["manual-input", "sloped-rectangle"],
    handler: slopedRect
  },
  {
    semanticName: "Multi-Document",
    name: "Stacked Document",
    shortName: "docs",
    description: "Multiple documents",
    aliases: ["documents", "st-doc", "stacked-document"],
    handler: multiWaveEdgedRectangle
  },
  {
    semanticName: "Multi-Process",
    name: "Stacked Rectangle",
    shortName: "st-rect",
    description: "Multiple processes",
    aliases: ["procs", "processes", "stacked-rectangle"],
    handler: multiRect
  },
  {
    semanticName: "Stored Data",
    name: "Bow Tie Rectangle",
    shortName: "bow-rect",
    description: "Stored data",
    aliases: ["stored-data", "bow-tie-rectangle"],
    handler: bowTieRect
  },
  {
    semanticName: "Summary",
    name: "Crossed Circle",
    shortName: "cross-circ",
    description: "Summary",
    aliases: ["summary", "crossed-circle"],
    handler: crossedCircle
  },
  {
    semanticName: "Tagged Document",
    name: "Tagged Document",
    shortName: "tag-doc",
    description: "Tagged document",
    aliases: ["tag-doc", "tagged-document"],
    handler: taggedWaveEdgedRectangle
  },
  {
    semanticName: "Tagged Process",
    name: "Tagged Rectangle",
    shortName: "tag-rect",
    description: "Tagged process",
    aliases: ["tagged-rectangle", "tag-proc", "tagged-process"],
    handler: taggedRect
  },
  {
    semanticName: "Paper Tape",
    name: "Flag",
    shortName: "flag",
    description: "Paper tape",
    aliases: ["paper-tape"],
    handler: waveRectangle
  },
  {
    semanticName: "Odd",
    name: "Odd",
    shortName: "odd",
    description: "Odd shape",
    internalAliases: ["rect_left_inv_arrow"],
    handler: rect_left_inv_arrow
  },
  {
    semanticName: "Lined Document",
    name: "Lined Document",
    shortName: "lin-doc",
    description: "Lined document",
    aliases: ["lined-document"],
    handler: linedWaveEdgedRect
  }
];
var generateShapeMap = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(() => {
  const undocumentedShapes = {
    // States
    state,
    choice,
    note,
    // Rectangles
    rectWithTitle,
    labelRect,
    // Icons
    iconSquare,
    iconCircle,
    icon,
    iconRounded,
    imageSquare,
    anchor,
    // Kanban diagram
    kanbanItem,
    //Mindmap diagram
    mindmapCircle,
    defaultMindmapNode,
    // class diagram
    classBox,
    // er diagram
    erBox,
    // Requirement diagram
    requirementBox
  };
  const entries = [
    ...Object.entries(undocumentedShapes),
    ...shapesDefs.flatMap((shape) => {
      const aliases = [
        shape.shortName,
        ..."aliases" in shape ? shape.aliases : [],
        ..."internalAliases" in shape ? shape.internalAliases : []
      ];
      return aliases.map((alias) => [alias, shape.handler]);
    })
  ];
  return Object.fromEntries(entries);
}, "generateShapeMap");
var shapes2 = generateShapeMap();
function isValidShape(shape) {
  return shape in shapes2;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(isValidShape, "isValidShape");

// src/rendering-util/rendering-elements/nodes.ts
var nodeElems = /* @__PURE__ */ new Map();
async function insertNode(elem, node, renderOptions) {
  let newEl;
  let el;
  if (node.shape === "rect") {
    if (node.rx && node.ry) {
      node.shape = "roundedRect";
    } else {
      node.shape = "squareRect";
    }
  }
  const shapeHandler = node.shape ? shapes2[node.shape] : void 0;
  if (!shapeHandler) {
    throw new Error(`No such shape: ${node.shape}. Please check your syntax.`);
  }
  if (node.link) {
    let target;
    if (renderOptions.config.securityLevel === "sandbox") {
      target = "_top";
    } else if (node.linkTarget) {
      target = node.linkTarget || "_blank";
    }
    newEl = elem.insert("svg:a").attr("xlink:href", node.link).attr("target", target ?? null);
    el = await shapeHandler(newEl, node, renderOptions);
  } else {
    el = await shapeHandler(elem, node, renderOptions);
    newEl = el;
  }
  if (node.tooltip) {
    el.attr("title", node.tooltip);
  }
  nodeElems.set(node.id, newEl);
  if (node.haveCallback) {
    newEl.attr("class", newEl.attr("class") + " clickable");
  }
  return newEl;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(insertNode, "insertNode");
var setNodeElem = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((elem, node) => {
  nodeElems.set(node.id, elem);
}, "setNodeElem");
var clear2 = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)(() => {
  nodeElems.clear();
}, "clear");
var positionNode = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .__name */ .eW)((node) => {
  const el = nodeElems.get(node.id);
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_5__/* .log */ .cM.trace(
    "Transforming node",
    node.diff,
    node,
    "translate(" + (node.x - node.width / 2 - 5) + ", " + node.width / 2 + ")"
  );
  const padding = 8;
  const diff = node.diff || 0;
  if (node.clusterNode) {
    el.attr(
      "transform",
      "translate(" + (node.x + diff - node.width / 2) + ", " + (node.y - node.height / 2 - padding) + ")"
    );
  } else {
    el.attr("transform", "translate(" + node.x + ", " + node.y + ")");
  }
  return diff;
}, "positionNode");




/***/ }),

/***/ 71262:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   A: () => (/* binding */ JSON_SCHEMA),
/* harmony export */   z: () => (/* binding */ load)
/* harmony export */ });
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(74999);


// ../../node_modules/.pnpm/js-yaml@4.1.0/node_modules/js-yaml/dist/js-yaml.mjs
function isNothing(subject) {
  return typeof subject === "undefined" || subject === null;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isNothing, "isNothing");
function isObject(subject) {
  return typeof subject === "object" && subject !== null;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isObject, "isObject");
function toArray(sequence) {
  if (Array.isArray(sequence)) return sequence;
  else if (isNothing(sequence)) return [];
  return [sequence];
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(toArray, "toArray");
function extend(target, source) {
  var index, length, key, sourceKeys;
  if (source) {
    sourceKeys = Object.keys(source);
    for (index = 0, length = sourceKeys.length; index < length; index += 1) {
      key = sourceKeys[index];
      target[key] = source[key];
    }
  }
  return target;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(extend, "extend");
function repeat(string, count) {
  var result = "", cycle;
  for (cycle = 0; cycle < count; cycle += 1) {
    result += string;
  }
  return result;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(repeat, "repeat");
function isNegativeZero(number) {
  return number === 0 && Number.NEGATIVE_INFINITY === 1 / number;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isNegativeZero, "isNegativeZero");
var isNothing_1 = isNothing;
var isObject_1 = isObject;
var toArray_1 = toArray;
var repeat_1 = repeat;
var isNegativeZero_1 = isNegativeZero;
var extend_1 = extend;
var common = {
  isNothing: isNothing_1,
  isObject: isObject_1,
  toArray: toArray_1,
  repeat: repeat_1,
  isNegativeZero: isNegativeZero_1,
  extend: extend_1
};
function formatError(exception2, compact) {
  var where = "", message = exception2.reason || "(unknown reason)";
  if (!exception2.mark) return message;
  if (exception2.mark.name) {
    where += 'in "' + exception2.mark.name + '" ';
  }
  where += "(" + (exception2.mark.line + 1) + ":" + (exception2.mark.column + 1) + ")";
  if (!compact && exception2.mark.snippet) {
    where += "\n\n" + exception2.mark.snippet;
  }
  return message + " " + where;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(formatError, "formatError");
function YAMLException$1(reason, mark) {
  Error.call(this);
  this.name = "YAMLException";
  this.reason = reason;
  this.mark = mark;
  this.message = formatError(this, false);
  if (Error.captureStackTrace) {
    Error.captureStackTrace(this, this.constructor);
  } else {
    this.stack = new Error().stack || "";
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(YAMLException$1, "YAMLException$1");
YAMLException$1.prototype = Object.create(Error.prototype);
YAMLException$1.prototype.constructor = YAMLException$1;
YAMLException$1.prototype.toString = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function toString(compact) {
  return this.name + ": " + formatError(this, compact);
}, "toString");
var exception = YAMLException$1;
function getLine(buffer, lineStart, lineEnd, position, maxLineLength) {
  var head = "";
  var tail = "";
  var maxHalfLength = Math.floor(maxLineLength / 2) - 1;
  if (position - lineStart > maxHalfLength) {
    head = " ... ";
    lineStart = position - maxHalfLength + head.length;
  }
  if (lineEnd - position > maxHalfLength) {
    tail = " ...";
    lineEnd = position + maxHalfLength - tail.length;
  }
  return {
    str: head + buffer.slice(lineStart, lineEnd).replace(/\t/g, "\u2192") + tail,
    pos: position - lineStart + head.length
    // relative position
  };
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(getLine, "getLine");
function padStart(string, max) {
  return common.repeat(" ", max - string.length) + string;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(padStart, "padStart");
function makeSnippet(mark, options) {
  options = Object.create(options || null);
  if (!mark.buffer) return null;
  if (!options.maxLength) options.maxLength = 79;
  if (typeof options.indent !== "number") options.indent = 1;
  if (typeof options.linesBefore !== "number") options.linesBefore = 3;
  if (typeof options.linesAfter !== "number") options.linesAfter = 2;
  var re = /\r?\n|\r|\0/g;
  var lineStarts = [0];
  var lineEnds = [];
  var match;
  var foundLineNo = -1;
  while (match = re.exec(mark.buffer)) {
    lineEnds.push(match.index);
    lineStarts.push(match.index + match[0].length);
    if (mark.position <= match.index && foundLineNo < 0) {
      foundLineNo = lineStarts.length - 2;
    }
  }
  if (foundLineNo < 0) foundLineNo = lineStarts.length - 1;
  var result = "", i, line;
  var lineNoLength = Math.min(mark.line + options.linesAfter, lineEnds.length).toString().length;
  var maxLineLength = options.maxLength - (options.indent + lineNoLength + 3);
  for (i = 1; i <= options.linesBefore; i++) {
    if (foundLineNo - i < 0) break;
    line = getLine(
      mark.buffer,
      lineStarts[foundLineNo - i],
      lineEnds[foundLineNo - i],
      mark.position - (lineStarts[foundLineNo] - lineStarts[foundLineNo - i]),
      maxLineLength
    );
    result = common.repeat(" ", options.indent) + padStart((mark.line - i + 1).toString(), lineNoLength) + " | " + line.str + "\n" + result;
  }
  line = getLine(mark.buffer, lineStarts[foundLineNo], lineEnds[foundLineNo], mark.position, maxLineLength);
  result += common.repeat(" ", options.indent) + padStart((mark.line + 1).toString(), lineNoLength) + " | " + line.str + "\n";
  result += common.repeat("-", options.indent + lineNoLength + 3 + line.pos) + "^\n";
  for (i = 1; i <= options.linesAfter; i++) {
    if (foundLineNo + i >= lineEnds.length) break;
    line = getLine(
      mark.buffer,
      lineStarts[foundLineNo + i],
      lineEnds[foundLineNo + i],
      mark.position - (lineStarts[foundLineNo] - lineStarts[foundLineNo + i]),
      maxLineLength
    );
    result += common.repeat(" ", options.indent) + padStart((mark.line + i + 1).toString(), lineNoLength) + " | " + line.str + "\n";
  }
  return result.replace(/\n$/, "");
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(makeSnippet, "makeSnippet");
var snippet = makeSnippet;
var TYPE_CONSTRUCTOR_OPTIONS = [
  "kind",
  "multi",
  "resolve",
  "construct",
  "instanceOf",
  "predicate",
  "represent",
  "representName",
  "defaultStyle",
  "styleAliases"
];
var YAML_NODE_KINDS = [
  "scalar",
  "sequence",
  "mapping"
];
function compileStyleAliases(map2) {
  var result = {};
  if (map2 !== null) {
    Object.keys(map2).forEach(function(style) {
      map2[style].forEach(function(alias) {
        result[String(alias)] = style;
      });
    });
  }
  return result;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(compileStyleAliases, "compileStyleAliases");
function Type$1(tag, options) {
  options = options || {};
  Object.keys(options).forEach(function(name) {
    if (TYPE_CONSTRUCTOR_OPTIONS.indexOf(name) === -1) {
      throw new exception('Unknown option "' + name + '" is met in definition of "' + tag + '" YAML type.');
    }
  });
  this.options = options;
  this.tag = tag;
  this.kind = options["kind"] || null;
  this.resolve = options["resolve"] || function() {
    return true;
  };
  this.construct = options["construct"] || function(data) {
    return data;
  };
  this.instanceOf = options["instanceOf"] || null;
  this.predicate = options["predicate"] || null;
  this.represent = options["represent"] || null;
  this.representName = options["representName"] || null;
  this.defaultStyle = options["defaultStyle"] || null;
  this.multi = options["multi"] || false;
  this.styleAliases = compileStyleAliases(options["styleAliases"] || null);
  if (YAML_NODE_KINDS.indexOf(this.kind) === -1) {
    throw new exception('Unknown kind "' + this.kind + '" is specified for "' + tag + '" YAML type.');
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(Type$1, "Type$1");
var type = Type$1;
function compileList(schema2, name) {
  var result = [];
  schema2[name].forEach(function(currentType) {
    var newIndex = result.length;
    result.forEach(function(previousType, previousIndex) {
      if (previousType.tag === currentType.tag && previousType.kind === currentType.kind && previousType.multi === currentType.multi) {
        newIndex = previousIndex;
      }
    });
    result[newIndex] = currentType;
  });
  return result;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(compileList, "compileList");
function compileMap() {
  var result = {
    scalar: {},
    sequence: {},
    mapping: {},
    fallback: {},
    multi: {
      scalar: [],
      sequence: [],
      mapping: [],
      fallback: []
    }
  }, index, length;
  function collectType(type2) {
    if (type2.multi) {
      result.multi[type2.kind].push(type2);
      result.multi["fallback"].push(type2);
    } else {
      result[type2.kind][type2.tag] = result["fallback"][type2.tag] = type2;
    }
  }
  (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(collectType, "collectType");
  for (index = 0, length = arguments.length; index < length; index += 1) {
    arguments[index].forEach(collectType);
  }
  return result;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(compileMap, "compileMap");
function Schema$1(definition) {
  return this.extend(definition);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(Schema$1, "Schema$1");
Schema$1.prototype.extend = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function extend2(definition) {
  var implicit = [];
  var explicit = [];
  if (definition instanceof type) {
    explicit.push(definition);
  } else if (Array.isArray(definition)) {
    explicit = explicit.concat(definition);
  } else if (definition && (Array.isArray(definition.implicit) || Array.isArray(definition.explicit))) {
    if (definition.implicit) implicit = implicit.concat(definition.implicit);
    if (definition.explicit) explicit = explicit.concat(definition.explicit);
  } else {
    throw new exception("Schema.extend argument should be a Type, [ Type ], or a schema definition ({ implicit: [...], explicit: [...] })");
  }
  implicit.forEach(function(type$1) {
    if (!(type$1 instanceof type)) {
      throw new exception("Specified list of YAML types (or a single Type object) contains a non-Type object.");
    }
    if (type$1.loadKind && type$1.loadKind !== "scalar") {
      throw new exception("There is a non-scalar type in the implicit list of a schema. Implicit resolving of such types is not supported.");
    }
    if (type$1.multi) {
      throw new exception("There is a multi type in the implicit list of a schema. Multi tags can only be listed as explicit.");
    }
  });
  explicit.forEach(function(type$1) {
    if (!(type$1 instanceof type)) {
      throw new exception("Specified list of YAML types (or a single Type object) contains a non-Type object.");
    }
  });
  var result = Object.create(Schema$1.prototype);
  result.implicit = (this.implicit || []).concat(implicit);
  result.explicit = (this.explicit || []).concat(explicit);
  result.compiledImplicit = compileList(result, "implicit");
  result.compiledExplicit = compileList(result, "explicit");
  result.compiledTypeMap = compileMap(result.compiledImplicit, result.compiledExplicit);
  return result;
}, "extend");
var schema = Schema$1;
var str = new type("tag:yaml.org,2002:str", {
  kind: "scalar",
  construct: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function(data) {
    return data !== null ? data : "";
  }, "construct")
});
var seq = new type("tag:yaml.org,2002:seq", {
  kind: "sequence",
  construct: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function(data) {
    return data !== null ? data : [];
  }, "construct")
});
var map = new type("tag:yaml.org,2002:map", {
  kind: "mapping",
  construct: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function(data) {
    return data !== null ? data : {};
  }, "construct")
});
var failsafe = new schema({
  explicit: [
    str,
    seq,
    map
  ]
});
function resolveYamlNull(data) {
  if (data === null) return true;
  var max = data.length;
  return max === 1 && data === "~" || max === 4 && (data === "null" || data === "Null" || data === "NULL");
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(resolveYamlNull, "resolveYamlNull");
function constructYamlNull() {
  return null;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(constructYamlNull, "constructYamlNull");
function isNull(object) {
  return object === null;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isNull, "isNull");
var _null = new type("tag:yaml.org,2002:null", {
  kind: "scalar",
  resolve: resolveYamlNull,
  construct: constructYamlNull,
  predicate: isNull,
  represent: {
    canonical: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function() {
      return "~";
    }, "canonical"),
    lowercase: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function() {
      return "null";
    }, "lowercase"),
    uppercase: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function() {
      return "NULL";
    }, "uppercase"),
    camelcase: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function() {
      return "Null";
    }, "camelcase"),
    empty: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function() {
      return "";
    }, "empty")
  },
  defaultStyle: "lowercase"
});
function resolveYamlBoolean(data) {
  if (data === null) return false;
  var max = data.length;
  return max === 4 && (data === "true" || data === "True" || data === "TRUE") || max === 5 && (data === "false" || data === "False" || data === "FALSE");
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(resolveYamlBoolean, "resolveYamlBoolean");
function constructYamlBoolean(data) {
  return data === "true" || data === "True" || data === "TRUE";
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(constructYamlBoolean, "constructYamlBoolean");
function isBoolean(object) {
  return Object.prototype.toString.call(object) === "[object Boolean]";
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isBoolean, "isBoolean");
var bool = new type("tag:yaml.org,2002:bool", {
  kind: "scalar",
  resolve: resolveYamlBoolean,
  construct: constructYamlBoolean,
  predicate: isBoolean,
  represent: {
    lowercase: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function(object) {
      return object ? "true" : "false";
    }, "lowercase"),
    uppercase: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function(object) {
      return object ? "TRUE" : "FALSE";
    }, "uppercase"),
    camelcase: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function(object) {
      return object ? "True" : "False";
    }, "camelcase")
  },
  defaultStyle: "lowercase"
});
function isHexCode(c) {
  return 48 <= c && c <= 57 || 65 <= c && c <= 70 || 97 <= c && c <= 102;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isHexCode, "isHexCode");
function isOctCode(c) {
  return 48 <= c && c <= 55;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isOctCode, "isOctCode");
function isDecCode(c) {
  return 48 <= c && c <= 57;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isDecCode, "isDecCode");
function resolveYamlInteger(data) {
  if (data === null) return false;
  var max = data.length, index = 0, hasDigits = false, ch;
  if (!max) return false;
  ch = data[index];
  if (ch === "-" || ch === "+") {
    ch = data[++index];
  }
  if (ch === "0") {
    if (index + 1 === max) return true;
    ch = data[++index];
    if (ch === "b") {
      index++;
      for (; index < max; index++) {
        ch = data[index];
        if (ch === "_") continue;
        if (ch !== "0" && ch !== "1") return false;
        hasDigits = true;
      }
      return hasDigits && ch !== "_";
    }
    if (ch === "x") {
      index++;
      for (; index < max; index++) {
        ch = data[index];
        if (ch === "_") continue;
        if (!isHexCode(data.charCodeAt(index))) return false;
        hasDigits = true;
      }
      return hasDigits && ch !== "_";
    }
    if (ch === "o") {
      index++;
      for (; index < max; index++) {
        ch = data[index];
        if (ch === "_") continue;
        if (!isOctCode(data.charCodeAt(index))) return false;
        hasDigits = true;
      }
      return hasDigits && ch !== "_";
    }
  }
  if (ch === "_") return false;
  for (; index < max; index++) {
    ch = data[index];
    if (ch === "_") continue;
    if (!isDecCode(data.charCodeAt(index))) {
      return false;
    }
    hasDigits = true;
  }
  if (!hasDigits || ch === "_") return false;
  return true;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(resolveYamlInteger, "resolveYamlInteger");
function constructYamlInteger(data) {
  var value = data, sign = 1, ch;
  if (value.indexOf("_") !== -1) {
    value = value.replace(/_/g, "");
  }
  ch = value[0];
  if (ch === "-" || ch === "+") {
    if (ch === "-") sign = -1;
    value = value.slice(1);
    ch = value[0];
  }
  if (value === "0") return 0;
  if (ch === "0") {
    if (value[1] === "b") return sign * parseInt(value.slice(2), 2);
    if (value[1] === "x") return sign * parseInt(value.slice(2), 16);
    if (value[1] === "o") return sign * parseInt(value.slice(2), 8);
  }
  return sign * parseInt(value, 10);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(constructYamlInteger, "constructYamlInteger");
function isInteger(object) {
  return Object.prototype.toString.call(object) === "[object Number]" && (object % 1 === 0 && !common.isNegativeZero(object));
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isInteger, "isInteger");
var int = new type("tag:yaml.org,2002:int", {
  kind: "scalar",
  resolve: resolveYamlInteger,
  construct: constructYamlInteger,
  predicate: isInteger,
  represent: {
    binary: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function(obj) {
      return obj >= 0 ? "0b" + obj.toString(2) : "-0b" + obj.toString(2).slice(1);
    }, "binary"),
    octal: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function(obj) {
      return obj >= 0 ? "0o" + obj.toString(8) : "-0o" + obj.toString(8).slice(1);
    }, "octal"),
    decimal: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function(obj) {
      return obj.toString(10);
    }, "decimal"),
    /* eslint-disable max-len */
    hexadecimal: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function(obj) {
      return obj >= 0 ? "0x" + obj.toString(16).toUpperCase() : "-0x" + obj.toString(16).toUpperCase().slice(1);
    }, "hexadecimal")
  },
  defaultStyle: "decimal",
  styleAliases: {
    binary: [2, "bin"],
    octal: [8, "oct"],
    decimal: [10, "dec"],
    hexadecimal: [16, "hex"]
  }
});
var YAML_FLOAT_PATTERN = new RegExp(
  // 2.5e4, 2.5 and integers
  "^(?:[-+]?(?:[0-9][0-9_]*)(?:\\.[0-9_]*)?(?:[eE][-+]?[0-9]+)?|\\.[0-9_]+(?:[eE][-+]?[0-9]+)?|[-+]?\\.(?:inf|Inf|INF)|\\.(?:nan|NaN|NAN))$"
);
function resolveYamlFloat(data) {
  if (data === null) return false;
  if (!YAML_FLOAT_PATTERN.test(data) || // Quick hack to not allow integers end with `_`
  // Probably should update regexp & check speed
  data[data.length - 1] === "_") {
    return false;
  }
  return true;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(resolveYamlFloat, "resolveYamlFloat");
function constructYamlFloat(data) {
  var value, sign;
  value = data.replace(/_/g, "").toLowerCase();
  sign = value[0] === "-" ? -1 : 1;
  if ("+-".indexOf(value[0]) >= 0) {
    value = value.slice(1);
  }
  if (value === ".inf") {
    return sign === 1 ? Number.POSITIVE_INFINITY : Number.NEGATIVE_INFINITY;
  } else if (value === ".nan") {
    return NaN;
  }
  return sign * parseFloat(value, 10);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(constructYamlFloat, "constructYamlFloat");
var SCIENTIFIC_WITHOUT_DOT = /^[-+]?[0-9]+e/;
function representYamlFloat(object, style) {
  var res;
  if (isNaN(object)) {
    switch (style) {
      case "lowercase":
        return ".nan";
      case "uppercase":
        return ".NAN";
      case "camelcase":
        return ".NaN";
    }
  } else if (Number.POSITIVE_INFINITY === object) {
    switch (style) {
      case "lowercase":
        return ".inf";
      case "uppercase":
        return ".INF";
      case "camelcase":
        return ".Inf";
    }
  } else if (Number.NEGATIVE_INFINITY === object) {
    switch (style) {
      case "lowercase":
        return "-.inf";
      case "uppercase":
        return "-.INF";
      case "camelcase":
        return "-.Inf";
    }
  } else if (common.isNegativeZero(object)) {
    return "-0.0";
  }
  res = object.toString(10);
  return SCIENTIFIC_WITHOUT_DOT.test(res) ? res.replace("e", ".e") : res;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(representYamlFloat, "representYamlFloat");
function isFloat(object) {
  return Object.prototype.toString.call(object) === "[object Number]" && (object % 1 !== 0 || common.isNegativeZero(object));
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isFloat, "isFloat");
var float = new type("tag:yaml.org,2002:float", {
  kind: "scalar",
  resolve: resolveYamlFloat,
  construct: constructYamlFloat,
  predicate: isFloat,
  represent: representYamlFloat,
  defaultStyle: "lowercase"
});
var json = failsafe.extend({
  implicit: [
    _null,
    bool,
    int,
    float
  ]
});
var core = json;
var YAML_DATE_REGEXP = new RegExp(
  "^([0-9][0-9][0-9][0-9])-([0-9][0-9])-([0-9][0-9])$"
);
var YAML_TIMESTAMP_REGEXP = new RegExp(
  "^([0-9][0-9][0-9][0-9])-([0-9][0-9]?)-([0-9][0-9]?)(?:[Tt]|[ \\t]+)([0-9][0-9]?):([0-9][0-9]):([0-9][0-9])(?:\\.([0-9]*))?(?:[ \\t]*(Z|([-+])([0-9][0-9]?)(?::([0-9][0-9]))?))?$"
);
function resolveYamlTimestamp(data) {
  if (data === null) return false;
  if (YAML_DATE_REGEXP.exec(data) !== null) return true;
  if (YAML_TIMESTAMP_REGEXP.exec(data) !== null) return true;
  return false;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(resolveYamlTimestamp, "resolveYamlTimestamp");
function constructYamlTimestamp(data) {
  var match, year, month, day, hour, minute, second, fraction = 0, delta = null, tz_hour, tz_minute, date;
  match = YAML_DATE_REGEXP.exec(data);
  if (match === null) match = YAML_TIMESTAMP_REGEXP.exec(data);
  if (match === null) throw new Error("Date resolve error");
  year = +match[1];
  month = +match[2] - 1;
  day = +match[3];
  if (!match[4]) {
    return new Date(Date.UTC(year, month, day));
  }
  hour = +match[4];
  minute = +match[5];
  second = +match[6];
  if (match[7]) {
    fraction = match[7].slice(0, 3);
    while (fraction.length < 3) {
      fraction += "0";
    }
    fraction = +fraction;
  }
  if (match[9]) {
    tz_hour = +match[10];
    tz_minute = +(match[11] || 0);
    delta = (tz_hour * 60 + tz_minute) * 6e4;
    if (match[9] === "-") delta = -delta;
  }
  date = new Date(Date.UTC(year, month, day, hour, minute, second, fraction));
  if (delta) date.setTime(date.getTime() - delta);
  return date;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(constructYamlTimestamp, "constructYamlTimestamp");
function representYamlTimestamp(object) {
  return object.toISOString();
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(representYamlTimestamp, "representYamlTimestamp");
var timestamp = new type("tag:yaml.org,2002:timestamp", {
  kind: "scalar",
  resolve: resolveYamlTimestamp,
  construct: constructYamlTimestamp,
  instanceOf: Date,
  represent: representYamlTimestamp
});
function resolveYamlMerge(data) {
  return data === "<<" || data === null;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(resolveYamlMerge, "resolveYamlMerge");
var merge = new type("tag:yaml.org,2002:merge", {
  kind: "scalar",
  resolve: resolveYamlMerge
});
var BASE64_MAP = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=\n\r";
function resolveYamlBinary(data) {
  if (data === null) return false;
  var code, idx, bitlen = 0, max = data.length, map2 = BASE64_MAP;
  for (idx = 0; idx < max; idx++) {
    code = map2.indexOf(data.charAt(idx));
    if (code > 64) continue;
    if (code < 0) return false;
    bitlen += 6;
  }
  return bitlen % 8 === 0;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(resolveYamlBinary, "resolveYamlBinary");
function constructYamlBinary(data) {
  var idx, tailbits, input = data.replace(/[\r\n=]/g, ""), max = input.length, map2 = BASE64_MAP, bits = 0, result = [];
  for (idx = 0; idx < max; idx++) {
    if (idx % 4 === 0 && idx) {
      result.push(bits >> 16 & 255);
      result.push(bits >> 8 & 255);
      result.push(bits & 255);
    }
    bits = bits << 6 | map2.indexOf(input.charAt(idx));
  }
  tailbits = max % 4 * 6;
  if (tailbits === 0) {
    result.push(bits >> 16 & 255);
    result.push(bits >> 8 & 255);
    result.push(bits & 255);
  } else if (tailbits === 18) {
    result.push(bits >> 10 & 255);
    result.push(bits >> 2 & 255);
  } else if (tailbits === 12) {
    result.push(bits >> 4 & 255);
  }
  return new Uint8Array(result);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(constructYamlBinary, "constructYamlBinary");
function representYamlBinary(object) {
  var result = "", bits = 0, idx, tail, max = object.length, map2 = BASE64_MAP;
  for (idx = 0; idx < max; idx++) {
    if (idx % 3 === 0 && idx) {
      result += map2[bits >> 18 & 63];
      result += map2[bits >> 12 & 63];
      result += map2[bits >> 6 & 63];
      result += map2[bits & 63];
    }
    bits = (bits << 8) + object[idx];
  }
  tail = max % 3;
  if (tail === 0) {
    result += map2[bits >> 18 & 63];
    result += map2[bits >> 12 & 63];
    result += map2[bits >> 6 & 63];
    result += map2[bits & 63];
  } else if (tail === 2) {
    result += map2[bits >> 10 & 63];
    result += map2[bits >> 4 & 63];
    result += map2[bits << 2 & 63];
    result += map2[64];
  } else if (tail === 1) {
    result += map2[bits >> 2 & 63];
    result += map2[bits << 4 & 63];
    result += map2[64];
    result += map2[64];
  }
  return result;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(representYamlBinary, "representYamlBinary");
function isBinary(obj) {
  return Object.prototype.toString.call(obj) === "[object Uint8Array]";
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isBinary, "isBinary");
var binary = new type("tag:yaml.org,2002:binary", {
  kind: "scalar",
  resolve: resolveYamlBinary,
  construct: constructYamlBinary,
  predicate: isBinary,
  represent: representYamlBinary
});
var _hasOwnProperty$3 = Object.prototype.hasOwnProperty;
var _toString$2 = Object.prototype.toString;
function resolveYamlOmap(data) {
  if (data === null) return true;
  var objectKeys = [], index, length, pair, pairKey, pairHasKey, object = data;
  for (index = 0, length = object.length; index < length; index += 1) {
    pair = object[index];
    pairHasKey = false;
    if (_toString$2.call(pair) !== "[object Object]") return false;
    for (pairKey in pair) {
      if (_hasOwnProperty$3.call(pair, pairKey)) {
        if (!pairHasKey) pairHasKey = true;
        else return false;
      }
    }
    if (!pairHasKey) return false;
    if (objectKeys.indexOf(pairKey) === -1) objectKeys.push(pairKey);
    else return false;
  }
  return true;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(resolveYamlOmap, "resolveYamlOmap");
function constructYamlOmap(data) {
  return data !== null ? data : [];
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(constructYamlOmap, "constructYamlOmap");
var omap = new type("tag:yaml.org,2002:omap", {
  kind: "sequence",
  resolve: resolveYamlOmap,
  construct: constructYamlOmap
});
var _toString$1 = Object.prototype.toString;
function resolveYamlPairs(data) {
  if (data === null) return true;
  var index, length, pair, keys, result, object = data;
  result = new Array(object.length);
  for (index = 0, length = object.length; index < length; index += 1) {
    pair = object[index];
    if (_toString$1.call(pair) !== "[object Object]") return false;
    keys = Object.keys(pair);
    if (keys.length !== 1) return false;
    result[index] = [keys[0], pair[keys[0]]];
  }
  return true;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(resolveYamlPairs, "resolveYamlPairs");
function constructYamlPairs(data) {
  if (data === null) return [];
  var index, length, pair, keys, result, object = data;
  result = new Array(object.length);
  for (index = 0, length = object.length; index < length; index += 1) {
    pair = object[index];
    keys = Object.keys(pair);
    result[index] = [keys[0], pair[keys[0]]];
  }
  return result;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(constructYamlPairs, "constructYamlPairs");
var pairs = new type("tag:yaml.org,2002:pairs", {
  kind: "sequence",
  resolve: resolveYamlPairs,
  construct: constructYamlPairs
});
var _hasOwnProperty$2 = Object.prototype.hasOwnProperty;
function resolveYamlSet(data) {
  if (data === null) return true;
  var key, object = data;
  for (key in object) {
    if (_hasOwnProperty$2.call(object, key)) {
      if (object[key] !== null) return false;
    }
  }
  return true;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(resolveYamlSet, "resolveYamlSet");
function constructYamlSet(data) {
  return data !== null ? data : {};
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(constructYamlSet, "constructYamlSet");
var set = new type("tag:yaml.org,2002:set", {
  kind: "mapping",
  resolve: resolveYamlSet,
  construct: constructYamlSet
});
var _default = core.extend({
  implicit: [
    timestamp,
    merge
  ],
  explicit: [
    binary,
    omap,
    pairs,
    set
  ]
});
var _hasOwnProperty$1 = Object.prototype.hasOwnProperty;
var CONTEXT_FLOW_IN = 1;
var CONTEXT_FLOW_OUT = 2;
var CONTEXT_BLOCK_IN = 3;
var CONTEXT_BLOCK_OUT = 4;
var CHOMPING_CLIP = 1;
var CHOMPING_STRIP = 2;
var CHOMPING_KEEP = 3;
var PATTERN_NON_PRINTABLE = /[\x00-\x08\x0B\x0C\x0E-\x1F\x7F-\x84\x86-\x9F\uFFFE\uFFFF]|[\uD800-\uDBFF](?![\uDC00-\uDFFF])|(?:[^\uD800-\uDBFF]|^)[\uDC00-\uDFFF]/;
var PATTERN_NON_ASCII_LINE_BREAKS = /[\x85\u2028\u2029]/;
var PATTERN_FLOW_INDICATORS = /[,\[\]\{\}]/;
var PATTERN_TAG_HANDLE = /^(?:!|!!|![a-z\-]+!)$/i;
var PATTERN_TAG_URI = /^(?:!|[^,\[\]\{\}])(?:%[0-9a-f]{2}|[0-9a-z\-#;\/\?:@&=\+\$,_\.!~\*'\(\)\[\]])*$/i;
function _class(obj) {
  return Object.prototype.toString.call(obj);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(_class, "_class");
function is_EOL(c) {
  return c === 10 || c === 13;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(is_EOL, "is_EOL");
function is_WHITE_SPACE(c) {
  return c === 9 || c === 32;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(is_WHITE_SPACE, "is_WHITE_SPACE");
function is_WS_OR_EOL(c) {
  return c === 9 || c === 32 || c === 10 || c === 13;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(is_WS_OR_EOL, "is_WS_OR_EOL");
function is_FLOW_INDICATOR(c) {
  return c === 44 || c === 91 || c === 93 || c === 123 || c === 125;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(is_FLOW_INDICATOR, "is_FLOW_INDICATOR");
function fromHexCode(c) {
  var lc;
  if (48 <= c && c <= 57) {
    return c - 48;
  }
  lc = c | 32;
  if (97 <= lc && lc <= 102) {
    return lc - 97 + 10;
  }
  return -1;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(fromHexCode, "fromHexCode");
function escapedHexLen(c) {
  if (c === 120) {
    return 2;
  }
  if (c === 117) {
    return 4;
  }
  if (c === 85) {
    return 8;
  }
  return 0;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(escapedHexLen, "escapedHexLen");
function fromDecimalCode(c) {
  if (48 <= c && c <= 57) {
    return c - 48;
  }
  return -1;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(fromDecimalCode, "fromDecimalCode");
function simpleEscapeSequence(c) {
  return c === 48 ? "\0" : c === 97 ? "\x07" : c === 98 ? "\b" : c === 116 ? "	" : c === 9 ? "	" : c === 110 ? "\n" : c === 118 ? "\v" : c === 102 ? "\f" : c === 114 ? "\r" : c === 101 ? "\x1B" : c === 32 ? " " : c === 34 ? '"' : c === 47 ? "/" : c === 92 ? "\\" : c === 78 ? "\x85" : c === 95 ? "\xA0" : c === 76 ? "\u2028" : c === 80 ? "\u2029" : "";
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(simpleEscapeSequence, "simpleEscapeSequence");
function charFromCodepoint(c) {
  if (c <= 65535) {
    return String.fromCharCode(c);
  }
  return String.fromCharCode(
    (c - 65536 >> 10) + 55296,
    (c - 65536 & 1023) + 56320
  );
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(charFromCodepoint, "charFromCodepoint");
var simpleEscapeCheck = new Array(256);
var simpleEscapeMap = new Array(256);
for (i = 0; i < 256; i++) {
  simpleEscapeCheck[i] = simpleEscapeSequence(i) ? 1 : 0;
  simpleEscapeMap[i] = simpleEscapeSequence(i);
}
var i;
function State$1(input, options) {
  this.input = input;
  this.filename = options["filename"] || null;
  this.schema = options["schema"] || _default;
  this.onWarning = options["onWarning"] || null;
  this.legacy = options["legacy"] || false;
  this.json = options["json"] || false;
  this.listener = options["listener"] || null;
  this.implicitTypes = this.schema.compiledImplicit;
  this.typeMap = this.schema.compiledTypeMap;
  this.length = input.length;
  this.position = 0;
  this.line = 0;
  this.lineStart = 0;
  this.lineIndent = 0;
  this.firstTabInLine = -1;
  this.documents = [];
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(State$1, "State$1");
function generateError(state, message) {
  var mark = {
    name: state.filename,
    buffer: state.input.slice(0, -1),
    // omit trailing \0
    position: state.position,
    line: state.line,
    column: state.position - state.lineStart
  };
  mark.snippet = snippet(mark);
  return new exception(message, mark);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(generateError, "generateError");
function throwError(state, message) {
  throw generateError(state, message);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(throwError, "throwError");
function throwWarning(state, message) {
  if (state.onWarning) {
    state.onWarning.call(null, generateError(state, message));
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(throwWarning, "throwWarning");
var directiveHandlers = {
  YAML: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function handleYamlDirective(state, name, args) {
    var match, major, minor;
    if (state.version !== null) {
      throwError(state, "duplication of %YAML directive");
    }
    if (args.length !== 1) {
      throwError(state, "YAML directive accepts exactly one argument");
    }
    match = /^([0-9]+)\.([0-9]+)$/.exec(args[0]);
    if (match === null) {
      throwError(state, "ill-formed argument of the YAML directive");
    }
    major = parseInt(match[1], 10);
    minor = parseInt(match[2], 10);
    if (major !== 1) {
      throwError(state, "unacceptable YAML version of the document");
    }
    state.version = args[0];
    state.checkLineBreaks = minor < 2;
    if (minor !== 1 && minor !== 2) {
      throwWarning(state, "unsupported YAML version of the document");
    }
  }, "handleYamlDirective"),
  TAG: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(function handleTagDirective(state, name, args) {
    var handle, prefix;
    if (args.length !== 2) {
      throwError(state, "TAG directive accepts exactly two arguments");
    }
    handle = args[0];
    prefix = args[1];
    if (!PATTERN_TAG_HANDLE.test(handle)) {
      throwError(state, "ill-formed tag handle (first argument) of the TAG directive");
    }
    if (_hasOwnProperty$1.call(state.tagMap, handle)) {
      throwError(state, 'there is a previously declared suffix for "' + handle + '" tag handle');
    }
    if (!PATTERN_TAG_URI.test(prefix)) {
      throwError(state, "ill-formed tag prefix (second argument) of the TAG directive");
    }
    try {
      prefix = decodeURIComponent(prefix);
    } catch (err) {
      throwError(state, "tag prefix is malformed: " + prefix);
    }
    state.tagMap[handle] = prefix;
  }, "handleTagDirective")
};
function captureSegment(state, start, end, checkJson) {
  var _position, _length, _character, _result;
  if (start < end) {
    _result = state.input.slice(start, end);
    if (checkJson) {
      for (_position = 0, _length = _result.length; _position < _length; _position += 1) {
        _character = _result.charCodeAt(_position);
        if (!(_character === 9 || 32 <= _character && _character <= 1114111)) {
          throwError(state, "expected valid JSON character");
        }
      }
    } else if (PATTERN_NON_PRINTABLE.test(_result)) {
      throwError(state, "the stream contains non-printable characters");
    }
    state.result += _result;
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(captureSegment, "captureSegment");
function mergeMappings(state, destination, source, overridableKeys) {
  var sourceKeys, key, index, quantity;
  if (!common.isObject(source)) {
    throwError(state, "cannot merge mappings; the provided source object is unacceptable");
  }
  sourceKeys = Object.keys(source);
  for (index = 0, quantity = sourceKeys.length; index < quantity; index += 1) {
    key = sourceKeys[index];
    if (!_hasOwnProperty$1.call(destination, key)) {
      destination[key] = source[key];
      overridableKeys[key] = true;
    }
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(mergeMappings, "mergeMappings");
function storeMappingPair(state, _result, overridableKeys, keyTag, keyNode, valueNode, startLine, startLineStart, startPos) {
  var index, quantity;
  if (Array.isArray(keyNode)) {
    keyNode = Array.prototype.slice.call(keyNode);
    for (index = 0, quantity = keyNode.length; index < quantity; index += 1) {
      if (Array.isArray(keyNode[index])) {
        throwError(state, "nested arrays are not supported inside keys");
      }
      if (typeof keyNode === "object" && _class(keyNode[index]) === "[object Object]") {
        keyNode[index] = "[object Object]";
      }
    }
  }
  if (typeof keyNode === "object" && _class(keyNode) === "[object Object]") {
    keyNode = "[object Object]";
  }
  keyNode = String(keyNode);
  if (_result === null) {
    _result = {};
  }
  if (keyTag === "tag:yaml.org,2002:merge") {
    if (Array.isArray(valueNode)) {
      for (index = 0, quantity = valueNode.length; index < quantity; index += 1) {
        mergeMappings(state, _result, valueNode[index], overridableKeys);
      }
    } else {
      mergeMappings(state, _result, valueNode, overridableKeys);
    }
  } else {
    if (!state.json && !_hasOwnProperty$1.call(overridableKeys, keyNode) && _hasOwnProperty$1.call(_result, keyNode)) {
      state.line = startLine || state.line;
      state.lineStart = startLineStart || state.lineStart;
      state.position = startPos || state.position;
      throwError(state, "duplicated mapping key");
    }
    if (keyNode === "__proto__") {
      Object.defineProperty(_result, keyNode, {
        configurable: true,
        enumerable: true,
        writable: true,
        value: valueNode
      });
    } else {
      _result[keyNode] = valueNode;
    }
    delete overridableKeys[keyNode];
  }
  return _result;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(storeMappingPair, "storeMappingPair");
function readLineBreak(state) {
  var ch;
  ch = state.input.charCodeAt(state.position);
  if (ch === 10) {
    state.position++;
  } else if (ch === 13) {
    state.position++;
    if (state.input.charCodeAt(state.position) === 10) {
      state.position++;
    }
  } else {
    throwError(state, "a line break is expected");
  }
  state.line += 1;
  state.lineStart = state.position;
  state.firstTabInLine = -1;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(readLineBreak, "readLineBreak");
function skipSeparationSpace(state, allowComments, checkIndent) {
  var lineBreaks = 0, ch = state.input.charCodeAt(state.position);
  while (ch !== 0) {
    while (is_WHITE_SPACE(ch)) {
      if (ch === 9 && state.firstTabInLine === -1) {
        state.firstTabInLine = state.position;
      }
      ch = state.input.charCodeAt(++state.position);
    }
    if (allowComments && ch === 35) {
      do {
        ch = state.input.charCodeAt(++state.position);
      } while (ch !== 10 && ch !== 13 && ch !== 0);
    }
    if (is_EOL(ch)) {
      readLineBreak(state);
      ch = state.input.charCodeAt(state.position);
      lineBreaks++;
      state.lineIndent = 0;
      while (ch === 32) {
        state.lineIndent++;
        ch = state.input.charCodeAt(++state.position);
      }
    } else {
      break;
    }
  }
  if (checkIndent !== -1 && lineBreaks !== 0 && state.lineIndent < checkIndent) {
    throwWarning(state, "deficient indentation");
  }
  return lineBreaks;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(skipSeparationSpace, "skipSeparationSpace");
function testDocumentSeparator(state) {
  var _position = state.position, ch;
  ch = state.input.charCodeAt(_position);
  if ((ch === 45 || ch === 46) && ch === state.input.charCodeAt(_position + 1) && ch === state.input.charCodeAt(_position + 2)) {
    _position += 3;
    ch = state.input.charCodeAt(_position);
    if (ch === 0 || is_WS_OR_EOL(ch)) {
      return true;
    }
  }
  return false;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(testDocumentSeparator, "testDocumentSeparator");
function writeFoldedLines(state, count) {
  if (count === 1) {
    state.result += " ";
  } else if (count > 1) {
    state.result += common.repeat("\n", count - 1);
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(writeFoldedLines, "writeFoldedLines");
function readPlainScalar(state, nodeIndent, withinFlowCollection) {
  var preceding, following, captureStart, captureEnd, hasPendingContent, _line, _lineStart, _lineIndent, _kind = state.kind, _result = state.result, ch;
  ch = state.input.charCodeAt(state.position);
  if (is_WS_OR_EOL(ch) || is_FLOW_INDICATOR(ch) || ch === 35 || ch === 38 || ch === 42 || ch === 33 || ch === 124 || ch === 62 || ch === 39 || ch === 34 || ch === 37 || ch === 64 || ch === 96) {
    return false;
  }
  if (ch === 63 || ch === 45) {
    following = state.input.charCodeAt(state.position + 1);
    if (is_WS_OR_EOL(following) || withinFlowCollection && is_FLOW_INDICATOR(following)) {
      return false;
    }
  }
  state.kind = "scalar";
  state.result = "";
  captureStart = captureEnd = state.position;
  hasPendingContent = false;
  while (ch !== 0) {
    if (ch === 58) {
      following = state.input.charCodeAt(state.position + 1);
      if (is_WS_OR_EOL(following) || withinFlowCollection && is_FLOW_INDICATOR(following)) {
        break;
      }
    } else if (ch === 35) {
      preceding = state.input.charCodeAt(state.position - 1);
      if (is_WS_OR_EOL(preceding)) {
        break;
      }
    } else if (state.position === state.lineStart && testDocumentSeparator(state) || withinFlowCollection && is_FLOW_INDICATOR(ch)) {
      break;
    } else if (is_EOL(ch)) {
      _line = state.line;
      _lineStart = state.lineStart;
      _lineIndent = state.lineIndent;
      skipSeparationSpace(state, false, -1);
      if (state.lineIndent >= nodeIndent) {
        hasPendingContent = true;
        ch = state.input.charCodeAt(state.position);
        continue;
      } else {
        state.position = captureEnd;
        state.line = _line;
        state.lineStart = _lineStart;
        state.lineIndent = _lineIndent;
        break;
      }
    }
    if (hasPendingContent) {
      captureSegment(state, captureStart, captureEnd, false);
      writeFoldedLines(state, state.line - _line);
      captureStart = captureEnd = state.position;
      hasPendingContent = false;
    }
    if (!is_WHITE_SPACE(ch)) {
      captureEnd = state.position + 1;
    }
    ch = state.input.charCodeAt(++state.position);
  }
  captureSegment(state, captureStart, captureEnd, false);
  if (state.result) {
    return true;
  }
  state.kind = _kind;
  state.result = _result;
  return false;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(readPlainScalar, "readPlainScalar");
function readSingleQuotedScalar(state, nodeIndent) {
  var ch, captureStart, captureEnd;
  ch = state.input.charCodeAt(state.position);
  if (ch !== 39) {
    return false;
  }
  state.kind = "scalar";
  state.result = "";
  state.position++;
  captureStart = captureEnd = state.position;
  while ((ch = state.input.charCodeAt(state.position)) !== 0) {
    if (ch === 39) {
      captureSegment(state, captureStart, state.position, true);
      ch = state.input.charCodeAt(++state.position);
      if (ch === 39) {
        captureStart = state.position;
        state.position++;
        captureEnd = state.position;
      } else {
        return true;
      }
    } else if (is_EOL(ch)) {
      captureSegment(state, captureStart, captureEnd, true);
      writeFoldedLines(state, skipSeparationSpace(state, false, nodeIndent));
      captureStart = captureEnd = state.position;
    } else if (state.position === state.lineStart && testDocumentSeparator(state)) {
      throwError(state, "unexpected end of the document within a single quoted scalar");
    } else {
      state.position++;
      captureEnd = state.position;
    }
  }
  throwError(state, "unexpected end of the stream within a single quoted scalar");
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(readSingleQuotedScalar, "readSingleQuotedScalar");
function readDoubleQuotedScalar(state, nodeIndent) {
  var captureStart, captureEnd, hexLength, hexResult, tmp, ch;
  ch = state.input.charCodeAt(state.position);
  if (ch !== 34) {
    return false;
  }
  state.kind = "scalar";
  state.result = "";
  state.position++;
  captureStart = captureEnd = state.position;
  while ((ch = state.input.charCodeAt(state.position)) !== 0) {
    if (ch === 34) {
      captureSegment(state, captureStart, state.position, true);
      state.position++;
      return true;
    } else if (ch === 92) {
      captureSegment(state, captureStart, state.position, true);
      ch = state.input.charCodeAt(++state.position);
      if (is_EOL(ch)) {
        skipSeparationSpace(state, false, nodeIndent);
      } else if (ch < 256 && simpleEscapeCheck[ch]) {
        state.result += simpleEscapeMap[ch];
        state.position++;
      } else if ((tmp = escapedHexLen(ch)) > 0) {
        hexLength = tmp;
        hexResult = 0;
        for (; hexLength > 0; hexLength--) {
          ch = state.input.charCodeAt(++state.position);
          if ((tmp = fromHexCode(ch)) >= 0) {
            hexResult = (hexResult << 4) + tmp;
          } else {
            throwError(state, "expected hexadecimal character");
          }
        }
        state.result += charFromCodepoint(hexResult);
        state.position++;
      } else {
        throwError(state, "unknown escape sequence");
      }
      captureStart = captureEnd = state.position;
    } else if (is_EOL(ch)) {
      captureSegment(state, captureStart, captureEnd, true);
      writeFoldedLines(state, skipSeparationSpace(state, false, nodeIndent));
      captureStart = captureEnd = state.position;
    } else if (state.position === state.lineStart && testDocumentSeparator(state)) {
      throwError(state, "unexpected end of the document within a double quoted scalar");
    } else {
      state.position++;
      captureEnd = state.position;
    }
  }
  throwError(state, "unexpected end of the stream within a double quoted scalar");
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(readDoubleQuotedScalar, "readDoubleQuotedScalar");
function readFlowCollection(state, nodeIndent) {
  var readNext = true, _line, _lineStart, _pos, _tag = state.tag, _result, _anchor = state.anchor, following, terminator, isPair, isExplicitPair, isMapping, overridableKeys = /* @__PURE__ */ Object.create(null), keyNode, keyTag, valueNode, ch;
  ch = state.input.charCodeAt(state.position);
  if (ch === 91) {
    terminator = 93;
    isMapping = false;
    _result = [];
  } else if (ch === 123) {
    terminator = 125;
    isMapping = true;
    _result = {};
  } else {
    return false;
  }
  if (state.anchor !== null) {
    state.anchorMap[state.anchor] = _result;
  }
  ch = state.input.charCodeAt(++state.position);
  while (ch !== 0) {
    skipSeparationSpace(state, true, nodeIndent);
    ch = state.input.charCodeAt(state.position);
    if (ch === terminator) {
      state.position++;
      state.tag = _tag;
      state.anchor = _anchor;
      state.kind = isMapping ? "mapping" : "sequence";
      state.result = _result;
      return true;
    } else if (!readNext) {
      throwError(state, "missed comma between flow collection entries");
    } else if (ch === 44) {
      throwError(state, "expected the node content, but found ','");
    }
    keyTag = keyNode = valueNode = null;
    isPair = isExplicitPair = false;
    if (ch === 63) {
      following = state.input.charCodeAt(state.position + 1);
      if (is_WS_OR_EOL(following)) {
        isPair = isExplicitPair = true;
        state.position++;
        skipSeparationSpace(state, true, nodeIndent);
      }
    }
    _line = state.line;
    _lineStart = state.lineStart;
    _pos = state.position;
    composeNode(state, nodeIndent, CONTEXT_FLOW_IN, false, true);
    keyTag = state.tag;
    keyNode = state.result;
    skipSeparationSpace(state, true, nodeIndent);
    ch = state.input.charCodeAt(state.position);
    if ((isExplicitPair || state.line === _line) && ch === 58) {
      isPair = true;
      ch = state.input.charCodeAt(++state.position);
      skipSeparationSpace(state, true, nodeIndent);
      composeNode(state, nodeIndent, CONTEXT_FLOW_IN, false, true);
      valueNode = state.result;
    }
    if (isMapping) {
      storeMappingPair(state, _result, overridableKeys, keyTag, keyNode, valueNode, _line, _lineStart, _pos);
    } else if (isPair) {
      _result.push(storeMappingPair(state, null, overridableKeys, keyTag, keyNode, valueNode, _line, _lineStart, _pos));
    } else {
      _result.push(keyNode);
    }
    skipSeparationSpace(state, true, nodeIndent);
    ch = state.input.charCodeAt(state.position);
    if (ch === 44) {
      readNext = true;
      ch = state.input.charCodeAt(++state.position);
    } else {
      readNext = false;
    }
  }
  throwError(state, "unexpected end of the stream within a flow collection");
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(readFlowCollection, "readFlowCollection");
function readBlockScalar(state, nodeIndent) {
  var captureStart, folding, chomping = CHOMPING_CLIP, didReadContent = false, detectedIndent = false, textIndent = nodeIndent, emptyLines = 0, atMoreIndented = false, tmp, ch;
  ch = state.input.charCodeAt(state.position);
  if (ch === 124) {
    folding = false;
  } else if (ch === 62) {
    folding = true;
  } else {
    return false;
  }
  state.kind = "scalar";
  state.result = "";
  while (ch !== 0) {
    ch = state.input.charCodeAt(++state.position);
    if (ch === 43 || ch === 45) {
      if (CHOMPING_CLIP === chomping) {
        chomping = ch === 43 ? CHOMPING_KEEP : CHOMPING_STRIP;
      } else {
        throwError(state, "repeat of a chomping mode identifier");
      }
    } else if ((tmp = fromDecimalCode(ch)) >= 0) {
      if (tmp === 0) {
        throwError(state, "bad explicit indentation width of a block scalar; it cannot be less than one");
      } else if (!detectedIndent) {
        textIndent = nodeIndent + tmp - 1;
        detectedIndent = true;
      } else {
        throwError(state, "repeat of an indentation width identifier");
      }
    } else {
      break;
    }
  }
  if (is_WHITE_SPACE(ch)) {
    do {
      ch = state.input.charCodeAt(++state.position);
    } while (is_WHITE_SPACE(ch));
    if (ch === 35) {
      do {
        ch = state.input.charCodeAt(++state.position);
      } while (!is_EOL(ch) && ch !== 0);
    }
  }
  while (ch !== 0) {
    readLineBreak(state);
    state.lineIndent = 0;
    ch = state.input.charCodeAt(state.position);
    while ((!detectedIndent || state.lineIndent < textIndent) && ch === 32) {
      state.lineIndent++;
      ch = state.input.charCodeAt(++state.position);
    }
    if (!detectedIndent && state.lineIndent > textIndent) {
      textIndent = state.lineIndent;
    }
    if (is_EOL(ch)) {
      emptyLines++;
      continue;
    }
    if (state.lineIndent < textIndent) {
      if (chomping === CHOMPING_KEEP) {
        state.result += common.repeat("\n", didReadContent ? 1 + emptyLines : emptyLines);
      } else if (chomping === CHOMPING_CLIP) {
        if (didReadContent) {
          state.result += "\n";
        }
      }
      break;
    }
    if (folding) {
      if (is_WHITE_SPACE(ch)) {
        atMoreIndented = true;
        state.result += common.repeat("\n", didReadContent ? 1 + emptyLines : emptyLines);
      } else if (atMoreIndented) {
        atMoreIndented = false;
        state.result += common.repeat("\n", emptyLines + 1);
      } else if (emptyLines === 0) {
        if (didReadContent) {
          state.result += " ";
        }
      } else {
        state.result += common.repeat("\n", emptyLines);
      }
    } else {
      state.result += common.repeat("\n", didReadContent ? 1 + emptyLines : emptyLines);
    }
    didReadContent = true;
    detectedIndent = true;
    emptyLines = 0;
    captureStart = state.position;
    while (!is_EOL(ch) && ch !== 0) {
      ch = state.input.charCodeAt(++state.position);
    }
    captureSegment(state, captureStart, state.position, false);
  }
  return true;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(readBlockScalar, "readBlockScalar");
function readBlockSequence(state, nodeIndent) {
  var _line, _tag = state.tag, _anchor = state.anchor, _result = [], following, detected = false, ch;
  if (state.firstTabInLine !== -1) return false;
  if (state.anchor !== null) {
    state.anchorMap[state.anchor] = _result;
  }
  ch = state.input.charCodeAt(state.position);
  while (ch !== 0) {
    if (state.firstTabInLine !== -1) {
      state.position = state.firstTabInLine;
      throwError(state, "tab characters must not be used in indentation");
    }
    if (ch !== 45) {
      break;
    }
    following = state.input.charCodeAt(state.position + 1);
    if (!is_WS_OR_EOL(following)) {
      break;
    }
    detected = true;
    state.position++;
    if (skipSeparationSpace(state, true, -1)) {
      if (state.lineIndent <= nodeIndent) {
        _result.push(null);
        ch = state.input.charCodeAt(state.position);
        continue;
      }
    }
    _line = state.line;
    composeNode(state, nodeIndent, CONTEXT_BLOCK_IN, false, true);
    _result.push(state.result);
    skipSeparationSpace(state, true, -1);
    ch = state.input.charCodeAt(state.position);
    if ((state.line === _line || state.lineIndent > nodeIndent) && ch !== 0) {
      throwError(state, "bad indentation of a sequence entry");
    } else if (state.lineIndent < nodeIndent) {
      break;
    }
  }
  if (detected) {
    state.tag = _tag;
    state.anchor = _anchor;
    state.kind = "sequence";
    state.result = _result;
    return true;
  }
  return false;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(readBlockSequence, "readBlockSequence");
function readBlockMapping(state, nodeIndent, flowIndent) {
  var following, allowCompact, _line, _keyLine, _keyLineStart, _keyPos, _tag = state.tag, _anchor = state.anchor, _result = {}, overridableKeys = /* @__PURE__ */ Object.create(null), keyTag = null, keyNode = null, valueNode = null, atExplicitKey = false, detected = false, ch;
  if (state.firstTabInLine !== -1) return false;
  if (state.anchor !== null) {
    state.anchorMap[state.anchor] = _result;
  }
  ch = state.input.charCodeAt(state.position);
  while (ch !== 0) {
    if (!atExplicitKey && state.firstTabInLine !== -1) {
      state.position = state.firstTabInLine;
      throwError(state, "tab characters must not be used in indentation");
    }
    following = state.input.charCodeAt(state.position + 1);
    _line = state.line;
    if ((ch === 63 || ch === 58) && is_WS_OR_EOL(following)) {
      if (ch === 63) {
        if (atExplicitKey) {
          storeMappingPair(state, _result, overridableKeys, keyTag, keyNode, null, _keyLine, _keyLineStart, _keyPos);
          keyTag = keyNode = valueNode = null;
        }
        detected = true;
        atExplicitKey = true;
        allowCompact = true;
      } else if (atExplicitKey) {
        atExplicitKey = false;
        allowCompact = true;
      } else {
        throwError(state, "incomplete explicit mapping pair; a key node is missed; or followed by a non-tabulated empty line");
      }
      state.position += 1;
      ch = following;
    } else {
      _keyLine = state.line;
      _keyLineStart = state.lineStart;
      _keyPos = state.position;
      if (!composeNode(state, flowIndent, CONTEXT_FLOW_OUT, false, true)) {
        break;
      }
      if (state.line === _line) {
        ch = state.input.charCodeAt(state.position);
        while (is_WHITE_SPACE(ch)) {
          ch = state.input.charCodeAt(++state.position);
        }
        if (ch === 58) {
          ch = state.input.charCodeAt(++state.position);
          if (!is_WS_OR_EOL(ch)) {
            throwError(state, "a whitespace character is expected after the key-value separator within a block mapping");
          }
          if (atExplicitKey) {
            storeMappingPair(state, _result, overridableKeys, keyTag, keyNode, null, _keyLine, _keyLineStart, _keyPos);
            keyTag = keyNode = valueNode = null;
          }
          detected = true;
          atExplicitKey = false;
          allowCompact = false;
          keyTag = state.tag;
          keyNode = state.result;
        } else if (detected) {
          throwError(state, "can not read an implicit mapping pair; a colon is missed");
        } else {
          state.tag = _tag;
          state.anchor = _anchor;
          return true;
        }
      } else if (detected) {
        throwError(state, "can not read a block mapping entry; a multiline key may not be an implicit key");
      } else {
        state.tag = _tag;
        state.anchor = _anchor;
        return true;
      }
    }
    if (state.line === _line || state.lineIndent > nodeIndent) {
      if (atExplicitKey) {
        _keyLine = state.line;
        _keyLineStart = state.lineStart;
        _keyPos = state.position;
      }
      if (composeNode(state, nodeIndent, CONTEXT_BLOCK_OUT, true, allowCompact)) {
        if (atExplicitKey) {
          keyNode = state.result;
        } else {
          valueNode = state.result;
        }
      }
      if (!atExplicitKey) {
        storeMappingPair(state, _result, overridableKeys, keyTag, keyNode, valueNode, _keyLine, _keyLineStart, _keyPos);
        keyTag = keyNode = valueNode = null;
      }
      skipSeparationSpace(state, true, -1);
      ch = state.input.charCodeAt(state.position);
    }
    if ((state.line === _line || state.lineIndent > nodeIndent) && ch !== 0) {
      throwError(state, "bad indentation of a mapping entry");
    } else if (state.lineIndent < nodeIndent) {
      break;
    }
  }
  if (atExplicitKey) {
    storeMappingPair(state, _result, overridableKeys, keyTag, keyNode, null, _keyLine, _keyLineStart, _keyPos);
  }
  if (detected) {
    state.tag = _tag;
    state.anchor = _anchor;
    state.kind = "mapping";
    state.result = _result;
  }
  return detected;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(readBlockMapping, "readBlockMapping");
function readTagProperty(state) {
  var _position, isVerbatim = false, isNamed = false, tagHandle, tagName, ch;
  ch = state.input.charCodeAt(state.position);
  if (ch !== 33) return false;
  if (state.tag !== null) {
    throwError(state, "duplication of a tag property");
  }
  ch = state.input.charCodeAt(++state.position);
  if (ch === 60) {
    isVerbatim = true;
    ch = state.input.charCodeAt(++state.position);
  } else if (ch === 33) {
    isNamed = true;
    tagHandle = "!!";
    ch = state.input.charCodeAt(++state.position);
  } else {
    tagHandle = "!";
  }
  _position = state.position;
  if (isVerbatim) {
    do {
      ch = state.input.charCodeAt(++state.position);
    } while (ch !== 0 && ch !== 62);
    if (state.position < state.length) {
      tagName = state.input.slice(_position, state.position);
      ch = state.input.charCodeAt(++state.position);
    } else {
      throwError(state, "unexpected end of the stream within a verbatim tag");
    }
  } else {
    while (ch !== 0 && !is_WS_OR_EOL(ch)) {
      if (ch === 33) {
        if (!isNamed) {
          tagHandle = state.input.slice(_position - 1, state.position + 1);
          if (!PATTERN_TAG_HANDLE.test(tagHandle)) {
            throwError(state, "named tag handle cannot contain such characters");
          }
          isNamed = true;
          _position = state.position + 1;
        } else {
          throwError(state, "tag suffix cannot contain exclamation marks");
        }
      }
      ch = state.input.charCodeAt(++state.position);
    }
    tagName = state.input.slice(_position, state.position);
    if (PATTERN_FLOW_INDICATORS.test(tagName)) {
      throwError(state, "tag suffix cannot contain flow indicator characters");
    }
  }
  if (tagName && !PATTERN_TAG_URI.test(tagName)) {
    throwError(state, "tag name cannot contain such characters: " + tagName);
  }
  try {
    tagName = decodeURIComponent(tagName);
  } catch (err) {
    throwError(state, "tag name is malformed: " + tagName);
  }
  if (isVerbatim) {
    state.tag = tagName;
  } else if (_hasOwnProperty$1.call(state.tagMap, tagHandle)) {
    state.tag = state.tagMap[tagHandle] + tagName;
  } else if (tagHandle === "!") {
    state.tag = "!" + tagName;
  } else if (tagHandle === "!!") {
    state.tag = "tag:yaml.org,2002:" + tagName;
  } else {
    throwError(state, 'undeclared tag handle "' + tagHandle + '"');
  }
  return true;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(readTagProperty, "readTagProperty");
function readAnchorProperty(state) {
  var _position, ch;
  ch = state.input.charCodeAt(state.position);
  if (ch !== 38) return false;
  if (state.anchor !== null) {
    throwError(state, "duplication of an anchor property");
  }
  ch = state.input.charCodeAt(++state.position);
  _position = state.position;
  while (ch !== 0 && !is_WS_OR_EOL(ch) && !is_FLOW_INDICATOR(ch)) {
    ch = state.input.charCodeAt(++state.position);
  }
  if (state.position === _position) {
    throwError(state, "name of an anchor node must contain at least one character");
  }
  state.anchor = state.input.slice(_position, state.position);
  return true;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(readAnchorProperty, "readAnchorProperty");
function readAlias(state) {
  var _position, alias, ch;
  ch = state.input.charCodeAt(state.position);
  if (ch !== 42) return false;
  ch = state.input.charCodeAt(++state.position);
  _position = state.position;
  while (ch !== 0 && !is_WS_OR_EOL(ch) && !is_FLOW_INDICATOR(ch)) {
    ch = state.input.charCodeAt(++state.position);
  }
  if (state.position === _position) {
    throwError(state, "name of an alias node must contain at least one character");
  }
  alias = state.input.slice(_position, state.position);
  if (!_hasOwnProperty$1.call(state.anchorMap, alias)) {
    throwError(state, 'unidentified alias "' + alias + '"');
  }
  state.result = state.anchorMap[alias];
  skipSeparationSpace(state, true, -1);
  return true;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(readAlias, "readAlias");
function composeNode(state, parentIndent, nodeContext, allowToSeek, allowCompact) {
  var allowBlockStyles, allowBlockScalars, allowBlockCollections, indentStatus = 1, atNewLine = false, hasContent = false, typeIndex, typeQuantity, typeList, type2, flowIndent, blockIndent;
  if (state.listener !== null) {
    state.listener("open", state);
  }
  state.tag = null;
  state.anchor = null;
  state.kind = null;
  state.result = null;
  allowBlockStyles = allowBlockScalars = allowBlockCollections = CONTEXT_BLOCK_OUT === nodeContext || CONTEXT_BLOCK_IN === nodeContext;
  if (allowToSeek) {
    if (skipSeparationSpace(state, true, -1)) {
      atNewLine = true;
      if (state.lineIndent > parentIndent) {
        indentStatus = 1;
      } else if (state.lineIndent === parentIndent) {
        indentStatus = 0;
      } else if (state.lineIndent < parentIndent) {
        indentStatus = -1;
      }
    }
  }
  if (indentStatus === 1) {
    while (readTagProperty(state) || readAnchorProperty(state)) {
      if (skipSeparationSpace(state, true, -1)) {
        atNewLine = true;
        allowBlockCollections = allowBlockStyles;
        if (state.lineIndent > parentIndent) {
          indentStatus = 1;
        } else if (state.lineIndent === parentIndent) {
          indentStatus = 0;
        } else if (state.lineIndent < parentIndent) {
          indentStatus = -1;
        }
      } else {
        allowBlockCollections = false;
      }
    }
  }
  if (allowBlockCollections) {
    allowBlockCollections = atNewLine || allowCompact;
  }
  if (indentStatus === 1 || CONTEXT_BLOCK_OUT === nodeContext) {
    if (CONTEXT_FLOW_IN === nodeContext || CONTEXT_FLOW_OUT === nodeContext) {
      flowIndent = parentIndent;
    } else {
      flowIndent = parentIndent + 1;
    }
    blockIndent = state.position - state.lineStart;
    if (indentStatus === 1) {
      if (allowBlockCollections && (readBlockSequence(state, blockIndent) || readBlockMapping(state, blockIndent, flowIndent)) || readFlowCollection(state, flowIndent)) {
        hasContent = true;
      } else {
        if (allowBlockScalars && readBlockScalar(state, flowIndent) || readSingleQuotedScalar(state, flowIndent) || readDoubleQuotedScalar(state, flowIndent)) {
          hasContent = true;
        } else if (readAlias(state)) {
          hasContent = true;
          if (state.tag !== null || state.anchor !== null) {
            throwError(state, "alias node should not have any properties");
          }
        } else if (readPlainScalar(state, flowIndent, CONTEXT_FLOW_IN === nodeContext)) {
          hasContent = true;
          if (state.tag === null) {
            state.tag = "?";
          }
        }
        if (state.anchor !== null) {
          state.anchorMap[state.anchor] = state.result;
        }
      }
    } else if (indentStatus === 0) {
      hasContent = allowBlockCollections && readBlockSequence(state, blockIndent);
    }
  }
  if (state.tag === null) {
    if (state.anchor !== null) {
      state.anchorMap[state.anchor] = state.result;
    }
  } else if (state.tag === "?") {
    if (state.result !== null && state.kind !== "scalar") {
      throwError(state, 'unacceptable node kind for !<?> tag; it should be "scalar", not "' + state.kind + '"');
    }
    for (typeIndex = 0, typeQuantity = state.implicitTypes.length; typeIndex < typeQuantity; typeIndex += 1) {
      type2 = state.implicitTypes[typeIndex];
      if (type2.resolve(state.result)) {
        state.result = type2.construct(state.result);
        state.tag = type2.tag;
        if (state.anchor !== null) {
          state.anchorMap[state.anchor] = state.result;
        }
        break;
      }
    }
  } else if (state.tag !== "!") {
    if (_hasOwnProperty$1.call(state.typeMap[state.kind || "fallback"], state.tag)) {
      type2 = state.typeMap[state.kind || "fallback"][state.tag];
    } else {
      type2 = null;
      typeList = state.typeMap.multi[state.kind || "fallback"];
      for (typeIndex = 0, typeQuantity = typeList.length; typeIndex < typeQuantity; typeIndex += 1) {
        if (state.tag.slice(0, typeList[typeIndex].tag.length) === typeList[typeIndex].tag) {
          type2 = typeList[typeIndex];
          break;
        }
      }
    }
    if (!type2) {
      throwError(state, "unknown tag !<" + state.tag + ">");
    }
    if (state.result !== null && type2.kind !== state.kind) {
      throwError(state, "unacceptable node kind for !<" + state.tag + '> tag; it should be "' + type2.kind + '", not "' + state.kind + '"');
    }
    if (!type2.resolve(state.result, state.tag)) {
      throwError(state, "cannot resolve a node with !<" + state.tag + "> explicit tag");
    } else {
      state.result = type2.construct(state.result, state.tag);
      if (state.anchor !== null) {
        state.anchorMap[state.anchor] = state.result;
      }
    }
  }
  if (state.listener !== null) {
    state.listener("close", state);
  }
  return state.tag !== null || state.anchor !== null || hasContent;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(composeNode, "composeNode");
function readDocument(state) {
  var documentStart = state.position, _position, directiveName, directiveArgs, hasDirectives = false, ch;
  state.version = null;
  state.checkLineBreaks = state.legacy;
  state.tagMap = /* @__PURE__ */ Object.create(null);
  state.anchorMap = /* @__PURE__ */ Object.create(null);
  while ((ch = state.input.charCodeAt(state.position)) !== 0) {
    skipSeparationSpace(state, true, -1);
    ch = state.input.charCodeAt(state.position);
    if (state.lineIndent > 0 || ch !== 37) {
      break;
    }
    hasDirectives = true;
    ch = state.input.charCodeAt(++state.position);
    _position = state.position;
    while (ch !== 0 && !is_WS_OR_EOL(ch)) {
      ch = state.input.charCodeAt(++state.position);
    }
    directiveName = state.input.slice(_position, state.position);
    directiveArgs = [];
    if (directiveName.length < 1) {
      throwError(state, "directive name must not be less than one character in length");
    }
    while (ch !== 0) {
      while (is_WHITE_SPACE(ch)) {
        ch = state.input.charCodeAt(++state.position);
      }
      if (ch === 35) {
        do {
          ch = state.input.charCodeAt(++state.position);
        } while (ch !== 0 && !is_EOL(ch));
        break;
      }
      if (is_EOL(ch)) break;
      _position = state.position;
      while (ch !== 0 && !is_WS_OR_EOL(ch)) {
        ch = state.input.charCodeAt(++state.position);
      }
      directiveArgs.push(state.input.slice(_position, state.position));
    }
    if (ch !== 0) readLineBreak(state);
    if (_hasOwnProperty$1.call(directiveHandlers, directiveName)) {
      directiveHandlers[directiveName](state, directiveName, directiveArgs);
    } else {
      throwWarning(state, 'unknown document directive "' + directiveName + '"');
    }
  }
  skipSeparationSpace(state, true, -1);
  if (state.lineIndent === 0 && state.input.charCodeAt(state.position) === 45 && state.input.charCodeAt(state.position + 1) === 45 && state.input.charCodeAt(state.position + 2) === 45) {
    state.position += 3;
    skipSeparationSpace(state, true, -1);
  } else if (hasDirectives) {
    throwError(state, "directives end mark is expected");
  }
  composeNode(state, state.lineIndent - 1, CONTEXT_BLOCK_OUT, false, true);
  skipSeparationSpace(state, true, -1);
  if (state.checkLineBreaks && PATTERN_NON_ASCII_LINE_BREAKS.test(state.input.slice(documentStart, state.position))) {
    throwWarning(state, "non-ASCII line breaks are interpreted as content");
  }
  state.documents.push(state.result);
  if (state.position === state.lineStart && testDocumentSeparator(state)) {
    if (state.input.charCodeAt(state.position) === 46) {
      state.position += 3;
      skipSeparationSpace(state, true, -1);
    }
    return;
  }
  if (state.position < state.length - 1) {
    throwError(state, "end of the stream or a document separator is expected");
  } else {
    return;
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(readDocument, "readDocument");
function loadDocuments(input, options) {
  input = String(input);
  options = options || {};
  if (input.length !== 0) {
    if (input.charCodeAt(input.length - 1) !== 10 && input.charCodeAt(input.length - 1) !== 13) {
      input += "\n";
    }
    if (input.charCodeAt(0) === 65279) {
      input = input.slice(1);
    }
  }
  var state = new State$1(input, options);
  var nullpos = input.indexOf("\0");
  if (nullpos !== -1) {
    state.position = nullpos;
    throwError(state, "null byte is not allowed in input");
  }
  state.input += "\0";
  while (state.input.charCodeAt(state.position) === 32) {
    state.lineIndent += 1;
    state.position += 1;
  }
  while (state.position < state.length - 1) {
    readDocument(state);
  }
  return state.documents;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(loadDocuments, "loadDocuments");
function loadAll$1(input, iterator, options) {
  if (iterator !== null && typeof iterator === "object" && typeof options === "undefined") {
    options = iterator;
    iterator = null;
  }
  var documents = loadDocuments(input, options);
  if (typeof iterator !== "function") {
    return documents;
  }
  for (var index = 0, length = documents.length; index < length; index += 1) {
    iterator(documents[index]);
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(loadAll$1, "loadAll$1");
function load$1(input, options) {
  var documents = loadDocuments(input, options);
  if (documents.length === 0) {
    return void 0;
  } else if (documents.length === 1) {
    return documents[0];
  }
  throw new exception("expected a single document in the stream, but found more");
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(load$1, "load$1");
var loadAll_1 = loadAll$1;
var load_1 = load$1;
var loader = {
  loadAll: loadAll_1,
  load: load_1
};
var _toString = Object.prototype.toString;
var _hasOwnProperty = Object.prototype.hasOwnProperty;
var CHAR_BOM = 65279;
var CHAR_TAB = 9;
var CHAR_LINE_FEED = 10;
var CHAR_CARRIAGE_RETURN = 13;
var CHAR_SPACE = 32;
var CHAR_EXCLAMATION = 33;
var CHAR_DOUBLE_QUOTE = 34;
var CHAR_SHARP = 35;
var CHAR_PERCENT = 37;
var CHAR_AMPERSAND = 38;
var CHAR_SINGLE_QUOTE = 39;
var CHAR_ASTERISK = 42;
var CHAR_COMMA = 44;
var CHAR_MINUS = 45;
var CHAR_COLON = 58;
var CHAR_EQUALS = 61;
var CHAR_GREATER_THAN = 62;
var CHAR_QUESTION = 63;
var CHAR_COMMERCIAL_AT = 64;
var CHAR_LEFT_SQUARE_BRACKET = 91;
var CHAR_RIGHT_SQUARE_BRACKET = 93;
var CHAR_GRAVE_ACCENT = 96;
var CHAR_LEFT_CURLY_BRACKET = 123;
var CHAR_VERTICAL_LINE = 124;
var CHAR_RIGHT_CURLY_BRACKET = 125;
var ESCAPE_SEQUENCES = {};
ESCAPE_SEQUENCES[0] = "\\0";
ESCAPE_SEQUENCES[7] = "\\a";
ESCAPE_SEQUENCES[8] = "\\b";
ESCAPE_SEQUENCES[9] = "\\t";
ESCAPE_SEQUENCES[10] = "\\n";
ESCAPE_SEQUENCES[11] = "\\v";
ESCAPE_SEQUENCES[12] = "\\f";
ESCAPE_SEQUENCES[13] = "\\r";
ESCAPE_SEQUENCES[27] = "\\e";
ESCAPE_SEQUENCES[34] = '\\"';
ESCAPE_SEQUENCES[92] = "\\\\";
ESCAPE_SEQUENCES[133] = "\\N";
ESCAPE_SEQUENCES[160] = "\\_";
ESCAPE_SEQUENCES[8232] = "\\L";
ESCAPE_SEQUENCES[8233] = "\\P";
var DEPRECATED_BOOLEANS_SYNTAX = [
  "y",
  "Y",
  "yes",
  "Yes",
  "YES",
  "on",
  "On",
  "ON",
  "n",
  "N",
  "no",
  "No",
  "NO",
  "off",
  "Off",
  "OFF"
];
var DEPRECATED_BASE60_SYNTAX = /^[-+]?[0-9_]+(?::[0-9_]+)+(?:\.[0-9_]*)?$/;
function compileStyleMap(schema2, map2) {
  var result, keys, index, length, tag, style, type2;
  if (map2 === null) return {};
  result = {};
  keys = Object.keys(map2);
  for (index = 0, length = keys.length; index < length; index += 1) {
    tag = keys[index];
    style = String(map2[tag]);
    if (tag.slice(0, 2) === "!!") {
      tag = "tag:yaml.org,2002:" + tag.slice(2);
    }
    type2 = schema2.compiledTypeMap["fallback"][tag];
    if (type2 && _hasOwnProperty.call(type2.styleAliases, style)) {
      style = type2.styleAliases[style];
    }
    result[tag] = style;
  }
  return result;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(compileStyleMap, "compileStyleMap");
function encodeHex(character) {
  var string, handle, length;
  string = character.toString(16).toUpperCase();
  if (character <= 255) {
    handle = "x";
    length = 2;
  } else if (character <= 65535) {
    handle = "u";
    length = 4;
  } else if (character <= 4294967295) {
    handle = "U";
    length = 8;
  } else {
    throw new exception("code point within a string may not be greater than 0xFFFFFFFF");
  }
  return "\\" + handle + common.repeat("0", length - string.length) + string;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(encodeHex, "encodeHex");
var QUOTING_TYPE_SINGLE = 1;
var QUOTING_TYPE_DOUBLE = 2;
function State(options) {
  this.schema = options["schema"] || _default;
  this.indent = Math.max(1, options["indent"] || 2);
  this.noArrayIndent = options["noArrayIndent"] || false;
  this.skipInvalid = options["skipInvalid"] || false;
  this.flowLevel = common.isNothing(options["flowLevel"]) ? -1 : options["flowLevel"];
  this.styleMap = compileStyleMap(this.schema, options["styles"] || null);
  this.sortKeys = options["sortKeys"] || false;
  this.lineWidth = options["lineWidth"] || 80;
  this.noRefs = options["noRefs"] || false;
  this.noCompatMode = options["noCompatMode"] || false;
  this.condenseFlow = options["condenseFlow"] || false;
  this.quotingType = options["quotingType"] === '"' ? QUOTING_TYPE_DOUBLE : QUOTING_TYPE_SINGLE;
  this.forceQuotes = options["forceQuotes"] || false;
  this.replacer = typeof options["replacer"] === "function" ? options["replacer"] : null;
  this.implicitTypes = this.schema.compiledImplicit;
  this.explicitTypes = this.schema.compiledExplicit;
  this.tag = null;
  this.result = "";
  this.duplicates = [];
  this.usedDuplicates = null;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(State, "State");
function indentString(string, spaces) {
  var ind = common.repeat(" ", spaces), position = 0, next = -1, result = "", line, length = string.length;
  while (position < length) {
    next = string.indexOf("\n", position);
    if (next === -1) {
      line = string.slice(position);
      position = length;
    } else {
      line = string.slice(position, next + 1);
      position = next + 1;
    }
    if (line.length && line !== "\n") result += ind;
    result += line;
  }
  return result;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(indentString, "indentString");
function generateNextLine(state, level) {
  return "\n" + common.repeat(" ", state.indent * level);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(generateNextLine, "generateNextLine");
function testImplicitResolving(state, str2) {
  var index, length, type2;
  for (index = 0, length = state.implicitTypes.length; index < length; index += 1) {
    type2 = state.implicitTypes[index];
    if (type2.resolve(str2)) {
      return true;
    }
  }
  return false;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(testImplicitResolving, "testImplicitResolving");
function isWhitespace(c) {
  return c === CHAR_SPACE || c === CHAR_TAB;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isWhitespace, "isWhitespace");
function isPrintable(c) {
  return 32 <= c && c <= 126 || 161 <= c && c <= 55295 && c !== 8232 && c !== 8233 || 57344 <= c && c <= 65533 && c !== CHAR_BOM || 65536 <= c && c <= 1114111;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isPrintable, "isPrintable");
function isNsCharOrWhitespace(c) {
  return isPrintable(c) && c !== CHAR_BOM && c !== CHAR_CARRIAGE_RETURN && c !== CHAR_LINE_FEED;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isNsCharOrWhitespace, "isNsCharOrWhitespace");
function isPlainSafe(c, prev, inblock) {
  var cIsNsCharOrWhitespace = isNsCharOrWhitespace(c);
  var cIsNsChar = cIsNsCharOrWhitespace && !isWhitespace(c);
  return (
    // ns-plain-safe
    (inblock ? (
      // c = flow-in
      cIsNsCharOrWhitespace
    ) : cIsNsCharOrWhitespace && c !== CHAR_COMMA && c !== CHAR_LEFT_SQUARE_BRACKET && c !== CHAR_RIGHT_SQUARE_BRACKET && c !== CHAR_LEFT_CURLY_BRACKET && c !== CHAR_RIGHT_CURLY_BRACKET) && c !== CHAR_SHARP && !(prev === CHAR_COLON && !cIsNsChar) || isNsCharOrWhitespace(prev) && !isWhitespace(prev) && c === CHAR_SHARP || prev === CHAR_COLON && cIsNsChar
  );
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isPlainSafe, "isPlainSafe");
function isPlainSafeFirst(c) {
  return isPrintable(c) && c !== CHAR_BOM && !isWhitespace(c) && c !== CHAR_MINUS && c !== CHAR_QUESTION && c !== CHAR_COLON && c !== CHAR_COMMA && c !== CHAR_LEFT_SQUARE_BRACKET && c !== CHAR_RIGHT_SQUARE_BRACKET && c !== CHAR_LEFT_CURLY_BRACKET && c !== CHAR_RIGHT_CURLY_BRACKET && c !== CHAR_SHARP && c !== CHAR_AMPERSAND && c !== CHAR_ASTERISK && c !== CHAR_EXCLAMATION && c !== CHAR_VERTICAL_LINE && c !== CHAR_EQUALS && c !== CHAR_GREATER_THAN && c !== CHAR_SINGLE_QUOTE && c !== CHAR_DOUBLE_QUOTE && c !== CHAR_PERCENT && c !== CHAR_COMMERCIAL_AT && c !== CHAR_GRAVE_ACCENT;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isPlainSafeFirst, "isPlainSafeFirst");
function isPlainSafeLast(c) {
  return !isWhitespace(c) && c !== CHAR_COLON;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(isPlainSafeLast, "isPlainSafeLast");
function codePointAt(string, pos) {
  var first = string.charCodeAt(pos), second;
  if (first >= 55296 && first <= 56319 && pos + 1 < string.length) {
    second = string.charCodeAt(pos + 1);
    if (second >= 56320 && second <= 57343) {
      return (first - 55296) * 1024 + second - 56320 + 65536;
    }
  }
  return first;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(codePointAt, "codePointAt");
function needIndentIndicator(string) {
  var leadingSpaceRe = /^\n* /;
  return leadingSpaceRe.test(string);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(needIndentIndicator, "needIndentIndicator");
var STYLE_PLAIN = 1;
var STYLE_SINGLE = 2;
var STYLE_LITERAL = 3;
var STYLE_FOLDED = 4;
var STYLE_DOUBLE = 5;
function chooseScalarStyle(string, singleLineOnly, indentPerLevel, lineWidth, testAmbiguousType, quotingType, forceQuotes, inblock) {
  var i;
  var char = 0;
  var prevChar = null;
  var hasLineBreak = false;
  var hasFoldableLine = false;
  var shouldTrackWidth = lineWidth !== -1;
  var previousLineBreak = -1;
  var plain = isPlainSafeFirst(codePointAt(string, 0)) && isPlainSafeLast(codePointAt(string, string.length - 1));
  if (singleLineOnly || forceQuotes) {
    for (i = 0; i < string.length; char >= 65536 ? i += 2 : i++) {
      char = codePointAt(string, i);
      if (!isPrintable(char)) {
        return STYLE_DOUBLE;
      }
      plain = plain && isPlainSafe(char, prevChar, inblock);
      prevChar = char;
    }
  } else {
    for (i = 0; i < string.length; char >= 65536 ? i += 2 : i++) {
      char = codePointAt(string, i);
      if (char === CHAR_LINE_FEED) {
        hasLineBreak = true;
        if (shouldTrackWidth) {
          hasFoldableLine = hasFoldableLine || // Foldable line = too long, and not more-indented.
          i - previousLineBreak - 1 > lineWidth && string[previousLineBreak + 1] !== " ";
          previousLineBreak = i;
        }
      } else if (!isPrintable(char)) {
        return STYLE_DOUBLE;
      }
      plain = plain && isPlainSafe(char, prevChar, inblock);
      prevChar = char;
    }
    hasFoldableLine = hasFoldableLine || shouldTrackWidth && (i - previousLineBreak - 1 > lineWidth && string[previousLineBreak + 1] !== " ");
  }
  if (!hasLineBreak && !hasFoldableLine) {
    if (plain && !forceQuotes && !testAmbiguousType(string)) {
      return STYLE_PLAIN;
    }
    return quotingType === QUOTING_TYPE_DOUBLE ? STYLE_DOUBLE : STYLE_SINGLE;
  }
  if (indentPerLevel > 9 && needIndentIndicator(string)) {
    return STYLE_DOUBLE;
  }
  if (!forceQuotes) {
    return hasFoldableLine ? STYLE_FOLDED : STYLE_LITERAL;
  }
  return quotingType === QUOTING_TYPE_DOUBLE ? STYLE_DOUBLE : STYLE_SINGLE;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(chooseScalarStyle, "chooseScalarStyle");
function writeScalar(state, string, level, iskey, inblock) {
  state.dump = (function() {
    if (string.length === 0) {
      return state.quotingType === QUOTING_TYPE_DOUBLE ? '""' : "''";
    }
    if (!state.noCompatMode) {
      if (DEPRECATED_BOOLEANS_SYNTAX.indexOf(string) !== -1 || DEPRECATED_BASE60_SYNTAX.test(string)) {
        return state.quotingType === QUOTING_TYPE_DOUBLE ? '"' + string + '"' : "'" + string + "'";
      }
    }
    var indent = state.indent * Math.max(1, level);
    var lineWidth = state.lineWidth === -1 ? -1 : Math.max(Math.min(state.lineWidth, 40), state.lineWidth - indent);
    var singleLineOnly = iskey || state.flowLevel > -1 && level >= state.flowLevel;
    function testAmbiguity(string2) {
      return testImplicitResolving(state, string2);
    }
    (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(testAmbiguity, "testAmbiguity");
    switch (chooseScalarStyle(
      string,
      singleLineOnly,
      state.indent,
      lineWidth,
      testAmbiguity,
      state.quotingType,
      state.forceQuotes && !iskey,
      inblock
    )) {
      case STYLE_PLAIN:
        return string;
      case STYLE_SINGLE:
        return "'" + string.replace(/'/g, "''") + "'";
      case STYLE_LITERAL:
        return "|" + blockHeader(string, state.indent) + dropEndingNewline(indentString(string, indent));
      case STYLE_FOLDED:
        return ">" + blockHeader(string, state.indent) + dropEndingNewline(indentString(foldString(string, lineWidth), indent));
      case STYLE_DOUBLE:
        return '"' + escapeString(string) + '"';
      default:
        throw new exception("impossible error: invalid scalar style");
    }
  })();
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(writeScalar, "writeScalar");
function blockHeader(string, indentPerLevel) {
  var indentIndicator = needIndentIndicator(string) ? String(indentPerLevel) : "";
  var clip = string[string.length - 1] === "\n";
  var keep = clip && (string[string.length - 2] === "\n" || string === "\n");
  var chomp = keep ? "+" : clip ? "" : "-";
  return indentIndicator + chomp + "\n";
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(blockHeader, "blockHeader");
function dropEndingNewline(string) {
  return string[string.length - 1] === "\n" ? string.slice(0, -1) : string;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(dropEndingNewline, "dropEndingNewline");
function foldString(string, width) {
  var lineRe = /(\n+)([^\n]*)/g;
  var result = (function() {
    var nextLF = string.indexOf("\n");
    nextLF = nextLF !== -1 ? nextLF : string.length;
    lineRe.lastIndex = nextLF;
    return foldLine(string.slice(0, nextLF), width);
  })();
  var prevMoreIndented = string[0] === "\n" || string[0] === " ";
  var moreIndented;
  var match;
  while (match = lineRe.exec(string)) {
    var prefix = match[1], line = match[2];
    moreIndented = line[0] === " ";
    result += prefix + (!prevMoreIndented && !moreIndented && line !== "" ? "\n" : "") + foldLine(line, width);
    prevMoreIndented = moreIndented;
  }
  return result;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(foldString, "foldString");
function foldLine(line, width) {
  if (line === "" || line[0] === " ") return line;
  var breakRe = / [^ ]/g;
  var match;
  var start = 0, end, curr = 0, next = 0;
  var result = "";
  while (match = breakRe.exec(line)) {
    next = match.index;
    if (next - start > width) {
      end = curr > start ? curr : next;
      result += "\n" + line.slice(start, end);
      start = end + 1;
    }
    curr = next;
  }
  result += "\n";
  if (line.length - start > width && curr > start) {
    result += line.slice(start, curr) + "\n" + line.slice(curr + 1);
  } else {
    result += line.slice(start);
  }
  return result.slice(1);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(foldLine, "foldLine");
function escapeString(string) {
  var result = "";
  var char = 0;
  var escapeSeq;
  for (var i = 0; i < string.length; char >= 65536 ? i += 2 : i++) {
    char = codePointAt(string, i);
    escapeSeq = ESCAPE_SEQUENCES[char];
    if (!escapeSeq && isPrintable(char)) {
      result += string[i];
      if (char >= 65536) result += string[i + 1];
    } else {
      result += escapeSeq || encodeHex(char);
    }
  }
  return result;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(escapeString, "escapeString");
function writeFlowSequence(state, level, object) {
  var _result = "", _tag = state.tag, index, length, value;
  for (index = 0, length = object.length; index < length; index += 1) {
    value = object[index];
    if (state.replacer) {
      value = state.replacer.call(object, String(index), value);
    }
    if (writeNode(state, level, value, false, false) || typeof value === "undefined" && writeNode(state, level, null, false, false)) {
      if (_result !== "") _result += "," + (!state.condenseFlow ? " " : "");
      _result += state.dump;
    }
  }
  state.tag = _tag;
  state.dump = "[" + _result + "]";
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(writeFlowSequence, "writeFlowSequence");
function writeBlockSequence(state, level, object, compact) {
  var _result = "", _tag = state.tag, index, length, value;
  for (index = 0, length = object.length; index < length; index += 1) {
    value = object[index];
    if (state.replacer) {
      value = state.replacer.call(object, String(index), value);
    }
    if (writeNode(state, level + 1, value, true, true, false, true) || typeof value === "undefined" && writeNode(state, level + 1, null, true, true, false, true)) {
      if (!compact || _result !== "") {
        _result += generateNextLine(state, level);
      }
      if (state.dump && CHAR_LINE_FEED === state.dump.charCodeAt(0)) {
        _result += "-";
      } else {
        _result += "- ";
      }
      _result += state.dump;
    }
  }
  state.tag = _tag;
  state.dump = _result || "[]";
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(writeBlockSequence, "writeBlockSequence");
function writeFlowMapping(state, level, object) {
  var _result = "", _tag = state.tag, objectKeyList = Object.keys(object), index, length, objectKey, objectValue, pairBuffer;
  for (index = 0, length = objectKeyList.length; index < length; index += 1) {
    pairBuffer = "";
    if (_result !== "") pairBuffer += ", ";
    if (state.condenseFlow) pairBuffer += '"';
    objectKey = objectKeyList[index];
    objectValue = object[objectKey];
    if (state.replacer) {
      objectValue = state.replacer.call(object, objectKey, objectValue);
    }
    if (!writeNode(state, level, objectKey, false, false)) {
      continue;
    }
    if (state.dump.length > 1024) pairBuffer += "? ";
    pairBuffer += state.dump + (state.condenseFlow ? '"' : "") + ":" + (state.condenseFlow ? "" : " ");
    if (!writeNode(state, level, objectValue, false, false)) {
      continue;
    }
    pairBuffer += state.dump;
    _result += pairBuffer;
  }
  state.tag = _tag;
  state.dump = "{" + _result + "}";
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(writeFlowMapping, "writeFlowMapping");
function writeBlockMapping(state, level, object, compact) {
  var _result = "", _tag = state.tag, objectKeyList = Object.keys(object), index, length, objectKey, objectValue, explicitPair, pairBuffer;
  if (state.sortKeys === true) {
    objectKeyList.sort();
  } else if (typeof state.sortKeys === "function") {
    objectKeyList.sort(state.sortKeys);
  } else if (state.sortKeys) {
    throw new exception("sortKeys must be a boolean or a function");
  }
  for (index = 0, length = objectKeyList.length; index < length; index += 1) {
    pairBuffer = "";
    if (!compact || _result !== "") {
      pairBuffer += generateNextLine(state, level);
    }
    objectKey = objectKeyList[index];
    objectValue = object[objectKey];
    if (state.replacer) {
      objectValue = state.replacer.call(object, objectKey, objectValue);
    }
    if (!writeNode(state, level + 1, objectKey, true, true, true)) {
      continue;
    }
    explicitPair = state.tag !== null && state.tag !== "?" || state.dump && state.dump.length > 1024;
    if (explicitPair) {
      if (state.dump && CHAR_LINE_FEED === state.dump.charCodeAt(0)) {
        pairBuffer += "?";
      } else {
        pairBuffer += "? ";
      }
    }
    pairBuffer += state.dump;
    if (explicitPair) {
      pairBuffer += generateNextLine(state, level);
    }
    if (!writeNode(state, level + 1, objectValue, true, explicitPair)) {
      continue;
    }
    if (state.dump && CHAR_LINE_FEED === state.dump.charCodeAt(0)) {
      pairBuffer += ":";
    } else {
      pairBuffer += ": ";
    }
    pairBuffer += state.dump;
    _result += pairBuffer;
  }
  state.tag = _tag;
  state.dump = _result || "{}";
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(writeBlockMapping, "writeBlockMapping");
function detectType(state, object, explicit) {
  var _result, typeList, index, length, type2, style;
  typeList = explicit ? state.explicitTypes : state.implicitTypes;
  for (index = 0, length = typeList.length; index < length; index += 1) {
    type2 = typeList[index];
    if ((type2.instanceOf || type2.predicate) && (!type2.instanceOf || typeof object === "object" && object instanceof type2.instanceOf) && (!type2.predicate || type2.predicate(object))) {
      if (explicit) {
        if (type2.multi && type2.representName) {
          state.tag = type2.representName(object);
        } else {
          state.tag = type2.tag;
        }
      } else {
        state.tag = "?";
      }
      if (type2.represent) {
        style = state.styleMap[type2.tag] || type2.defaultStyle;
        if (_toString.call(type2.represent) === "[object Function]") {
          _result = type2.represent(object, style);
        } else if (_hasOwnProperty.call(type2.represent, style)) {
          _result = type2.represent[style](object, style);
        } else {
          throw new exception("!<" + type2.tag + '> tag resolver accepts not "' + style + '" style');
        }
        state.dump = _result;
      }
      return true;
    }
  }
  return false;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(detectType, "detectType");
function writeNode(state, level, object, block, compact, iskey, isblockseq) {
  state.tag = null;
  state.dump = object;
  if (!detectType(state, object, false)) {
    detectType(state, object, true);
  }
  var type2 = _toString.call(state.dump);
  var inblock = block;
  var tagStr;
  if (block) {
    block = state.flowLevel < 0 || state.flowLevel > level;
  }
  var objectOrArray = type2 === "[object Object]" || type2 === "[object Array]", duplicateIndex, duplicate;
  if (objectOrArray) {
    duplicateIndex = state.duplicates.indexOf(object);
    duplicate = duplicateIndex !== -1;
  }
  if (state.tag !== null && state.tag !== "?" || duplicate || state.indent !== 2 && level > 0) {
    compact = false;
  }
  if (duplicate && state.usedDuplicates[duplicateIndex]) {
    state.dump = "*ref_" + duplicateIndex;
  } else {
    if (objectOrArray && duplicate && !state.usedDuplicates[duplicateIndex]) {
      state.usedDuplicates[duplicateIndex] = true;
    }
    if (type2 === "[object Object]") {
      if (block && Object.keys(state.dump).length !== 0) {
        writeBlockMapping(state, level, state.dump, compact);
        if (duplicate) {
          state.dump = "&ref_" + duplicateIndex + state.dump;
        }
      } else {
        writeFlowMapping(state, level, state.dump);
        if (duplicate) {
          state.dump = "&ref_" + duplicateIndex + " " + state.dump;
        }
      }
    } else if (type2 === "[object Array]") {
      if (block && state.dump.length !== 0) {
        if (state.noArrayIndent && !isblockseq && level > 0) {
          writeBlockSequence(state, level - 1, state.dump, compact);
        } else {
          writeBlockSequence(state, level, state.dump, compact);
        }
        if (duplicate) {
          state.dump = "&ref_" + duplicateIndex + state.dump;
        }
      } else {
        writeFlowSequence(state, level, state.dump);
        if (duplicate) {
          state.dump = "&ref_" + duplicateIndex + " " + state.dump;
        }
      }
    } else if (type2 === "[object String]") {
      if (state.tag !== "?") {
        writeScalar(state, state.dump, level, iskey, inblock);
      }
    } else if (type2 === "[object Undefined]") {
      return false;
    } else {
      if (state.skipInvalid) return false;
      throw new exception("unacceptable kind of an object to dump " + type2);
    }
    if (state.tag !== null && state.tag !== "?") {
      tagStr = encodeURI(
        state.tag[0] === "!" ? state.tag.slice(1) : state.tag
      ).replace(/!/g, "%21");
      if (state.tag[0] === "!") {
        tagStr = "!" + tagStr;
      } else if (tagStr.slice(0, 18) === "tag:yaml.org,2002:") {
        tagStr = "!!" + tagStr.slice(18);
      } else {
        tagStr = "!<" + tagStr + ">";
      }
      state.dump = tagStr + " " + state.dump;
    }
  }
  return true;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(writeNode, "writeNode");
function getDuplicateReferences(object, state) {
  var objects = [], duplicatesIndexes = [], index, length;
  inspectNode(object, objects, duplicatesIndexes);
  for (index = 0, length = duplicatesIndexes.length; index < length; index += 1) {
    state.duplicates.push(objects[duplicatesIndexes[index]]);
  }
  state.usedDuplicates = new Array(length);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(getDuplicateReferences, "getDuplicateReferences");
function inspectNode(object, objects, duplicatesIndexes) {
  var objectKeyList, index, length;
  if (object !== null && typeof object === "object") {
    index = objects.indexOf(object);
    if (index !== -1) {
      if (duplicatesIndexes.indexOf(index) === -1) {
        duplicatesIndexes.push(index);
      }
    } else {
      objects.push(object);
      if (Array.isArray(object)) {
        for (index = 0, length = object.length; index < length; index += 1) {
          inspectNode(object[index], objects, duplicatesIndexes);
        }
      } else {
        objectKeyList = Object.keys(object);
        for (index = 0, length = objectKeyList.length; index < length; index += 1) {
          inspectNode(object[objectKeyList[index]], objects, duplicatesIndexes);
        }
      }
    }
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(inspectNode, "inspectNode");
function dump$1(input, options) {
  options = options || {};
  var state = new State(options);
  if (!state.noRefs) getDuplicateReferences(input, state);
  var value = input;
  if (state.replacer) {
    value = state.replacer.call({ "": value }, "", value);
  }
  if (writeNode(state, 0, value, true, true)) return state.dump + "\n";
  return "";
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(dump$1, "dump$1");
var dump_1 = dump$1;
var dumper = {
  dump: dump_1
};
function renamed(from, to) {
  return function() {
    throw new Error("Function yaml." + from + " is removed in js-yaml 4. Use yaml." + to + " instead, which is now safe by default.");
  };
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(renamed, "renamed");
var JSON_SCHEMA = json;
var load = loader.load;
var loadAll = loader.loadAll;
var dump = dumper.dump;
var safeLoad = renamed("safeLoad", "load");
var safeLoadAll = renamed("safeLoadAll", "loadAll");
var safeDump = renamed("safeDump", "dump");


/*! Bundled license information:

js-yaml/dist/js-yaml.mjs:
  (*! js-yaml 4.1.0 https://github.com/nodeca/js-yaml @license MIT *)
*/


/***/ }),

/***/ 54160:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   _b: () => (/* binding */ getRegisteredLayoutAlgorithm),
/* harmony export */   jM: () => (/* binding */ registerLayoutLoaders),
/* harmony export */   sY: () => (/* binding */ render)
/* harmony export */ });
/* harmony import */ var _chunk_QXUST7PY_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(55195);
/* harmony import */ var _chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(26212);
/* harmony import */ var _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(17175);
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(74999);






// src/internals.ts
var internalHelpers = {
  common: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .common_default */ .SY,
  getConfig: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_3__/* .getConfig */ .iE,
  insertCluster: _chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_1__/* .insertCluster */ .us,
  insertEdge: _chunk_QXUST7PY_mjs__WEBPACK_IMPORTED_MODULE_0__/* .insertEdge */ .QP,
  insertEdgeLabel: _chunk_QXUST7PY_mjs__WEBPACK_IMPORTED_MODULE_0__/* .insertEdgeLabel */ .I_,
  insertMarkers: _chunk_QXUST7PY_mjs__WEBPACK_IMPORTED_MODULE_0__/* .markers_default */ .DQ,
  insertNode: _chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_1__/* .insertNode */ .Lf,
  interpolateToCurve: _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_2__/* .interpolateToCurve */ .le,
  labelHelper: _chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_1__/* .labelHelper */ .C1,
  log: _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .log */ .cM,
  positionEdgeLabel: _chunk_QXUST7PY_mjs__WEBPACK_IMPORTED_MODULE_0__/* .positionEdgeLabel */ .Jj
};

// src/rendering-util/render.ts
var layoutAlgorithms = {};
var registerLayoutLoaders = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)((loaders) => {
  for (const loader of loaders) {
    layoutAlgorithms[loader.name] = loader;
  }
}, "registerLayoutLoaders");
var registerDefaultLayoutLoaders = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)(() => {
  registerLayoutLoaders([
    {
      name: "dagre",
      loader: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)(async () => await Promise.all(/* import() */[__webpack_require__.e(2228), __webpack_require__.e(4276), __webpack_require__.e(7239)]).then(__webpack_require__.bind(__webpack_require__, 57239)), "loader")
    },
    ... true ? [
      {
        name: "cose-bilkent",
        loader: /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)(async () => await Promise.all(/* import() */[__webpack_require__.e(3207), __webpack_require__.e(4047)]).then(__webpack_require__.bind(__webpack_require__, 34047)), "loader")
      }
    ] : 0
  ]);
}, "registerDefaultLayoutLoaders");
registerDefaultLayoutLoaders();
var render = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)(async (data4Layout, svg) => {
  if (!(data4Layout.layoutAlgorithm in layoutAlgorithms)) {
    throw new Error(`Unknown layout algorithm: ${data4Layout.layoutAlgorithm}`);
  }
  const layoutDefinition = layoutAlgorithms[data4Layout.layoutAlgorithm];
  const layoutRenderer = await layoutDefinition.loader();
  return layoutRenderer.render(data4Layout, svg, internalHelpers, {
    algorithm: layoutDefinition.algorithm
  });
}, "render");
var getRegisteredLayoutAlgorithm = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .__name */ .eW)((algorithm = "", { fallback = "dagre" } = {}) => {
  if (algorithm in layoutAlgorithms) {
    return algorithm;
  }
  if (fallback in layoutAlgorithms) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_4__/* .log */ .cM.warn(`Layout algorithm ${algorithm} is not registered. Using ${fallback} as fallback.`);
    return fallback;
  }
  throw new Error(`Both layout algorithms ${algorithm} and ${fallback} are not registered.`);
}, "getRegisteredLayoutAlgorithm");




/***/ }),

/***/ 55195:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   DQ: () => (/* binding */ markers_default),
/* harmony export */   I_: () => (/* binding */ insertEdgeLabel),
/* harmony export */   Jj: () => (/* binding */ positionEdgeLabel),
/* harmony export */   QP: () => (/* binding */ insertEdge),
/* harmony export */   ZH: () => (/* binding */ clear)
/* harmony export */ });
/* harmony import */ var _chunk_HN2XXSSU_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(67894);
/* harmony import */ var _chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(26212);
/* harmony import */ var _chunk_CVBHYZKI_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(22462);
/* harmony import */ var _chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(37756);
/* harmony import */ var _chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(2111);
/* harmony import */ var _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(17175);
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(74999);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(35321);
/* harmony import */ var roughjs__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(36366);









// src/rendering-util/rendering-elements/edges.js



// src/rendering-util/rendering-elements/edgeMarker.ts
var addEdgeMarkers = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((svgPath, edge, url, id, diagramType, strokeColor) => {
  if (edge.arrowTypeStart) {
    addEdgeMarker(svgPath, "start", edge.arrowTypeStart, url, id, diagramType, strokeColor);
  }
  if (edge.arrowTypeEnd) {
    addEdgeMarker(svgPath, "end", edge.arrowTypeEnd, url, id, diagramType, strokeColor);
  }
}, "addEdgeMarkers");
var arrowTypesMap = {
  arrow_cross: { type: "cross", fill: false },
  arrow_point: { type: "point", fill: true },
  arrow_barb: { type: "barb", fill: true },
  arrow_circle: { type: "circle", fill: false },
  aggregation: { type: "aggregation", fill: false },
  extension: { type: "extension", fill: false },
  composition: { type: "composition", fill: true },
  dependency: { type: "dependency", fill: true },
  lollipop: { type: "lollipop", fill: false },
  only_one: { type: "onlyOne", fill: false },
  zero_or_one: { type: "zeroOrOne", fill: false },
  one_or_more: { type: "oneOrMore", fill: false },
  zero_or_more: { type: "zeroOrMore", fill: false },
  requirement_arrow: { type: "requirement_arrow", fill: false },
  requirement_contains: { type: "requirement_contains", fill: false }
};
var addEdgeMarker = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((svgPath, position, arrowType, url, id, diagramType, strokeColor) => {
  const arrowTypeInfo = arrowTypesMap[arrowType];
  if (!arrowTypeInfo) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.warn(`Unknown arrow type: ${arrowType}`);
    return;
  }
  const endMarkerType = arrowTypeInfo.type;
  const suffix = position === "start" ? "Start" : "End";
  const originalMarkerId = `${id}_${diagramType}-${endMarkerType}${suffix}`;
  if (strokeColor && strokeColor.trim() !== "") {
    const colorId = strokeColor.replace(/[^\dA-Za-z]/g, "_");
    const coloredMarkerId = `${originalMarkerId}_${colorId}`;
    if (!document.getElementById(coloredMarkerId)) {
      const originalMarker = document.getElementById(originalMarkerId);
      if (originalMarker) {
        const coloredMarker = originalMarker.cloneNode(true);
        coloredMarker.id = coloredMarkerId;
        const paths = coloredMarker.querySelectorAll("path, circle, line");
        paths.forEach((path) => {
          path.setAttribute("stroke", strokeColor);
          if (arrowTypeInfo.fill) {
            path.setAttribute("fill", strokeColor);
          }
        });
        originalMarker.parentNode?.appendChild(coloredMarker);
      }
    }
    svgPath.attr(`marker-${position}`, `url(${url}#${coloredMarkerId})`);
  } else {
    svgPath.attr(`marker-${position}`, `url(${url}#${originalMarkerId})`);
  }
}, "addEdgeMarker");

// src/rendering-util/rendering-elements/edges.js
var edgeLabels = /* @__PURE__ */ new Map();
var terminalLabels = /* @__PURE__ */ new Map();
var clear = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)(() => {
  edgeLabels.clear();
  terminalLabels.clear();
}, "clear");
var getLabelStyles = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((styleArray) => {
  let styles = styleArray ? styleArray.reduce((acc, style) => acc + ";" + style, "") : "";
  return styles;
}, "getLabelStyles");
var insertEdgeLabel = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)(async (elem, edge) => {
  let useHtmlLabels = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_6__/* .evaluate */ .ku)((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_6__/* .getConfig2 */ .nV)().flowchart.htmlLabels);
  const { labelStyles } = (0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_3__/* .styles2String */ .UG)(edge);
  edge.labelStyle = labelStyles;
  const labelElement = await (0,_chunk_JA3XYJ7Z_mjs__WEBPACK_IMPORTED_MODULE_4__/* .createText */ .rw)(elem, edge.label, {
    style: edge.labelStyle,
    useHtmlLabels,
    addSvgBackground: true,
    isNode: false
  });
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.info("abc82", edge, edge.labelType);
  const edgeLabel = elem.insert("g").attr("class", "edgeLabel");
  const label = edgeLabel.insert("g").attr("class", "label").attr("data-id", edge.id);
  label.node().appendChild(labelElement);
  let bbox = labelElement.getBBox();
  if (useHtmlLabels) {
    const div = labelElement.children[0];
    const dv = (0,d3__WEBPACK_IMPORTED_MODULE_8__/* .select */ .Ys)(labelElement);
    bbox = div.getBoundingClientRect();
    dv.attr("width", bbox.width);
    dv.attr("height", bbox.height);
  }
  label.attr("transform", "translate(" + -bbox.width / 2 + ", " + -bbox.height / 2 + ")");
  edgeLabels.set(edge.id, edgeLabel);
  edge.width = bbox.width;
  edge.height = bbox.height;
  let fo;
  if (edge.startLabelLeft) {
    const startLabelElement = await (0,_chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_1__/* .createLabel_default */ .XO)(
      edge.startLabelLeft,
      getLabelStyles(edge.labelStyle)
    );
    const startEdgeLabelLeft = elem.insert("g").attr("class", "edgeTerminals");
    const inner = startEdgeLabelLeft.insert("g").attr("class", "inner");
    fo = inner.node().appendChild(startLabelElement);
    const slBox = startLabelElement.getBBox();
    inner.attr("transform", "translate(" + -slBox.width / 2 + ", " + -slBox.height / 2 + ")");
    if (!terminalLabels.get(edge.id)) {
      terminalLabels.set(edge.id, {});
    }
    terminalLabels.get(edge.id).startLeft = startEdgeLabelLeft;
    setTerminalWidth(fo, edge.startLabelLeft);
  }
  if (edge.startLabelRight) {
    const startLabelElement = await (0,_chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_1__/* .createLabel_default */ .XO)(
      edge.startLabelRight,
      getLabelStyles(edge.labelStyle)
    );
    const startEdgeLabelRight = elem.insert("g").attr("class", "edgeTerminals");
    const inner = startEdgeLabelRight.insert("g").attr("class", "inner");
    fo = startEdgeLabelRight.node().appendChild(startLabelElement);
    inner.node().appendChild(startLabelElement);
    const slBox = startLabelElement.getBBox();
    inner.attr("transform", "translate(" + -slBox.width / 2 + ", " + -slBox.height / 2 + ")");
    if (!terminalLabels.get(edge.id)) {
      terminalLabels.set(edge.id, {});
    }
    terminalLabels.get(edge.id).startRight = startEdgeLabelRight;
    setTerminalWidth(fo, edge.startLabelRight);
  }
  if (edge.endLabelLeft) {
    const endLabelElement = await (0,_chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_1__/* .createLabel_default */ .XO)(edge.endLabelLeft, getLabelStyles(edge.labelStyle));
    const endEdgeLabelLeft = elem.insert("g").attr("class", "edgeTerminals");
    const inner = endEdgeLabelLeft.insert("g").attr("class", "inner");
    fo = inner.node().appendChild(endLabelElement);
    const slBox = endLabelElement.getBBox();
    inner.attr("transform", "translate(" + -slBox.width / 2 + ", " + -slBox.height / 2 + ")");
    endEdgeLabelLeft.node().appendChild(endLabelElement);
    if (!terminalLabels.get(edge.id)) {
      terminalLabels.set(edge.id, {});
    }
    terminalLabels.get(edge.id).endLeft = endEdgeLabelLeft;
    setTerminalWidth(fo, edge.endLabelLeft);
  }
  if (edge.endLabelRight) {
    const endLabelElement = await (0,_chunk_JZLCHNYA_mjs__WEBPACK_IMPORTED_MODULE_1__/* .createLabel_default */ .XO)(edge.endLabelRight, getLabelStyles(edge.labelStyle));
    const endEdgeLabelRight = elem.insert("g").attr("class", "edgeTerminals");
    const inner = endEdgeLabelRight.insert("g").attr("class", "inner");
    fo = inner.node().appendChild(endLabelElement);
    const slBox = endLabelElement.getBBox();
    inner.attr("transform", "translate(" + -slBox.width / 2 + ", " + -slBox.height / 2 + ")");
    endEdgeLabelRight.node().appendChild(endLabelElement);
    if (!terminalLabels.get(edge.id)) {
      terminalLabels.set(edge.id, {});
    }
    terminalLabels.get(edge.id).endRight = endEdgeLabelRight;
    setTerminalWidth(fo, edge.endLabelRight);
  }
  return labelElement;
}, "insertEdgeLabel");
function setTerminalWidth(fo, value) {
  if ((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_6__/* .getConfig2 */ .nV)().flowchart.htmlLabels && fo) {
    fo.style.width = value.length * 9 + "px";
    fo.style.height = "12px";
  }
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)(setTerminalWidth, "setTerminalWidth");
var positionEdgeLabel = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((edge, paths) => {
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.debug("Moving label abc88 ", edge.id, edge.label, edgeLabels.get(edge.id), paths);
  let path = paths.updatedPath ? paths.updatedPath : paths.originalPath;
  const siteConfig = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_6__/* .getConfig2 */ .nV)();
  const { subGraphTitleTotalMargin } = (0,_chunk_CVBHYZKI_mjs__WEBPACK_IMPORTED_MODULE_2__/* .getSubGraphTitleMargins */ .L)(siteConfig);
  if (edge.label) {
    const el = edgeLabels.get(edge.id);
    let x = edge.x;
    let y = edge.y;
    if (path) {
      const pos = _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_5__/* .utils_default */ .w8.calcLabelPosition(path);
      _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.debug(
        "Moving label " + edge.label + " from (",
        x,
        ",",
        y,
        ") to (",
        pos.x,
        ",",
        pos.y,
        ") abc88"
      );
      if (paths.updatedPath) {
        x = pos.x;
        y = pos.y;
      }
    }
    el.attr("transform", `translate(${x}, ${y + subGraphTitleTotalMargin / 2})`);
  }
  if (edge.startLabelLeft) {
    const el = terminalLabels.get(edge.id).startLeft;
    let x = edge.x;
    let y = edge.y;
    if (path) {
      const pos = _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_5__/* .utils_default */ .w8.calcTerminalLabelPosition(edge.arrowTypeStart ? 10 : 0, "start_left", path);
      x = pos.x;
      y = pos.y;
    }
    el.attr("transform", `translate(${x}, ${y})`);
  }
  if (edge.startLabelRight) {
    const el = terminalLabels.get(edge.id).startRight;
    let x = edge.x;
    let y = edge.y;
    if (path) {
      const pos = _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_5__/* .utils_default */ .w8.calcTerminalLabelPosition(
        edge.arrowTypeStart ? 10 : 0,
        "start_right",
        path
      );
      x = pos.x;
      y = pos.y;
    }
    el.attr("transform", `translate(${x}, ${y})`);
  }
  if (edge.endLabelLeft) {
    const el = terminalLabels.get(edge.id).endLeft;
    let x = edge.x;
    let y = edge.y;
    if (path) {
      const pos = _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_5__/* .utils_default */ .w8.calcTerminalLabelPosition(edge.arrowTypeEnd ? 10 : 0, "end_left", path);
      x = pos.x;
      y = pos.y;
    }
    el.attr("transform", `translate(${x}, ${y})`);
  }
  if (edge.endLabelRight) {
    const el = terminalLabels.get(edge.id).endRight;
    let x = edge.x;
    let y = edge.y;
    if (path) {
      const pos = _chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_5__/* .utils_default */ .w8.calcTerminalLabelPosition(edge.arrowTypeEnd ? 10 : 0, "end_right", path);
      x = pos.x;
      y = pos.y;
    }
    el.attr("transform", `translate(${x}, ${y})`);
  }
}, "positionEdgeLabel");
var outsideNode = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((node, point2) => {
  const x = node.x;
  const y = node.y;
  const dx = Math.abs(point2.x - x);
  const dy = Math.abs(point2.y - y);
  const w = node.width / 2;
  const h = node.height / 2;
  return dx >= w || dy >= h;
}, "outsideNode");
var intersection = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((node, outsidePoint, insidePoint) => {
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.debug(`intersection calc abc89:
  outsidePoint: ${JSON.stringify(outsidePoint)}
  insidePoint : ${JSON.stringify(insidePoint)}
  node        : x:${node.x} y:${node.y} w:${node.width} h:${node.height}`);
  const x = node.x;
  const y = node.y;
  const dx = Math.abs(x - insidePoint.x);
  const w = node.width / 2;
  let r = insidePoint.x < outsidePoint.x ? w - dx : w + dx;
  const h = node.height / 2;
  const Q = Math.abs(outsidePoint.y - insidePoint.y);
  const R = Math.abs(outsidePoint.x - insidePoint.x);
  if (Math.abs(y - outsidePoint.y) * w > Math.abs(x - outsidePoint.x) * h) {
    let q = insidePoint.y < outsidePoint.y ? outsidePoint.y - h - y : y - h - outsidePoint.y;
    r = R * q / Q;
    const res = {
      x: insidePoint.x < outsidePoint.x ? insidePoint.x + r : insidePoint.x - R + r,
      y: insidePoint.y < outsidePoint.y ? insidePoint.y + Q - q : insidePoint.y - Q + q
    };
    if (r === 0) {
      res.x = outsidePoint.x;
      res.y = outsidePoint.y;
    }
    if (R === 0) {
      res.x = outsidePoint.x;
    }
    if (Q === 0) {
      res.y = outsidePoint.y;
    }
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.debug(`abc89 top/bottom calc, Q ${Q}, q ${q}, R ${R}, r ${r}`, res);
    return res;
  } else {
    if (insidePoint.x < outsidePoint.x) {
      r = outsidePoint.x - w - x;
    } else {
      r = x - w - outsidePoint.x;
    }
    let q = Q * r / R;
    let _x = insidePoint.x < outsidePoint.x ? insidePoint.x + R - r : insidePoint.x - R + r;
    let _y = insidePoint.y < outsidePoint.y ? insidePoint.y + q : insidePoint.y - q;
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.debug(`sides calc abc89, Q ${Q}, q ${q}, R ${R}, r ${r}`, { _x, _y });
    if (r === 0) {
      _x = outsidePoint.x;
      _y = outsidePoint.y;
    }
    if (R === 0) {
      _x = outsidePoint.x;
    }
    if (Q === 0) {
      _y = outsidePoint.y;
    }
    return { x: _x, y: _y };
  }
}, "intersection");
var cutPathAtIntersect = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((_points, boundaryNode) => {
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.warn("abc88 cutPathAtIntersect", _points, boundaryNode);
  let points = [];
  let lastPointOutside = _points[0];
  let isInside = false;
  _points.forEach((point2) => {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.info("abc88 checking point", point2, boundaryNode);
    if (!outsideNode(boundaryNode, point2) && !isInside) {
      const inter = intersection(boundaryNode, lastPointOutside, point2);
      _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.debug("abc88 inside", point2, lastPointOutside, inter);
      _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.debug("abc88 intersection", inter, boundaryNode);
      let pointPresent = false;
      points.forEach((p) => {
        pointPresent = pointPresent || p.x === inter.x && p.y === inter.y;
      });
      if (!points.some((e) => e.x === inter.x && e.y === inter.y)) {
        points.push(inter);
      } else {
        _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.warn("abc88 no intersect", inter, points);
      }
      isInside = true;
    } else {
      _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.warn("abc88 outside", point2, lastPointOutside);
      lastPointOutside = point2;
      if (!isInside) {
        points.push(point2);
      }
    }
  });
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.debug("returning points", points);
  return points;
}, "cutPathAtIntersect");
function extractCornerPoints(points) {
  const cornerPoints = [];
  const cornerPointPositions = [];
  for (let i = 1; i < points.length - 1; i++) {
    const prev = points[i - 1];
    const curr = points[i];
    const next = points[i + 1];
    if (prev.x === curr.x && curr.y === next.y && Math.abs(curr.x - next.x) > 5 && Math.abs(curr.y - prev.y) > 5) {
      cornerPoints.push(curr);
      cornerPointPositions.push(i);
    } else if (prev.y === curr.y && curr.x === next.x && Math.abs(curr.x - prev.x) > 5 && Math.abs(curr.y - next.y) > 5) {
      cornerPoints.push(curr);
      cornerPointPositions.push(i);
    }
  }
  return { cornerPoints, cornerPointPositions };
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)(extractCornerPoints, "extractCornerPoints");
var findAdjacentPoint = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)(function(pointA, pointB, distance) {
  const xDiff = pointB.x - pointA.x;
  const yDiff = pointB.y - pointA.y;
  const length = Math.sqrt(xDiff * xDiff + yDiff * yDiff);
  const ratio = distance / length;
  return { x: pointB.x - ratio * xDiff, y: pointB.y - ratio * yDiff };
}, "findAdjacentPoint");
var fixCorners = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)(function(lineData) {
  const { cornerPointPositions } = extractCornerPoints(lineData);
  const newLineData = [];
  for (let i = 0; i < lineData.length; i++) {
    if (cornerPointPositions.includes(i)) {
      const prevPoint = lineData[i - 1];
      const nextPoint = lineData[i + 1];
      const cornerPoint = lineData[i];
      const newPrevPoint = findAdjacentPoint(prevPoint, cornerPoint, 5);
      const newNextPoint = findAdjacentPoint(nextPoint, cornerPoint, 5);
      const xDiff = newNextPoint.x - newPrevPoint.x;
      const yDiff = newNextPoint.y - newPrevPoint.y;
      newLineData.push(newPrevPoint);
      const a = Math.sqrt(2) * 2;
      let newCornerPoint = { x: cornerPoint.x, y: cornerPoint.y };
      if (Math.abs(nextPoint.x - prevPoint.x) > 10 && Math.abs(nextPoint.y - prevPoint.y) >= 10) {
        _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.debug(
          "Corner point fixing",
          Math.abs(nextPoint.x - prevPoint.x),
          Math.abs(nextPoint.y - prevPoint.y)
        );
        const r = 5;
        if (cornerPoint.x === newPrevPoint.x) {
          newCornerPoint = {
            x: xDiff < 0 ? newPrevPoint.x - r + a : newPrevPoint.x + r - a,
            y: yDiff < 0 ? newPrevPoint.y - a : newPrevPoint.y + a
          };
        } else {
          newCornerPoint = {
            x: xDiff < 0 ? newPrevPoint.x - a : newPrevPoint.x + a,
            y: yDiff < 0 ? newPrevPoint.y - r + a : newPrevPoint.y + r - a
          };
        }
      } else {
        _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.debug(
          "Corner point skipping fixing",
          Math.abs(nextPoint.x - prevPoint.x),
          Math.abs(nextPoint.y - prevPoint.y)
        );
      }
      newLineData.push(newCornerPoint, newNextPoint);
    } else {
      newLineData.push(lineData[i]);
    }
  }
  return newLineData;
}, "fixCorners");
var generateDashArray = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((len, oValueS, oValueE) => {
  const middleLength = len - oValueS - oValueE;
  const dashLength = 2;
  const gapLength = 2;
  const dashGapPairLength = dashLength + gapLength;
  const numberOfPairs = Math.floor(middleLength / dashGapPairLength);
  const middlePattern = Array(numberOfPairs).fill(`${dashLength} ${gapLength}`).join(" ");
  const dashArray = `0 ${oValueS} ${middlePattern} ${oValueE}`;
  return dashArray;
}, "generateDashArray");
var insertEdge = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)(function(elem, edge, clusterDb, diagramType, startNode, endNode, id, skipIntersect = false) {
  const { handDrawnSeed } = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_6__/* .getConfig2 */ .nV)();
  let points = edge.points;
  let pointsHasChanged = false;
  const tail = startNode;
  var head = endNode;
  const edgeClassStyles = [];
  for (const key in edge.cssCompiledStyles) {
    if ((0,_chunk_ATLVNIR6_mjs__WEBPACK_IMPORTED_MODULE_3__/* .isLabelStyle */ .Fh)(key)) {
      continue;
    }
    edgeClassStyles.push(edge.cssCompiledStyles[key]);
  }
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.debug("UIO intersect check", edge.points, head.x, tail.x);
  if (head.intersect && tail.intersect && !skipIntersect) {
    points = points.slice(1, edge.points.length - 1);
    points.unshift(tail.intersect(points[0]));
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.debug(
      "Last point UIO",
      edge.start,
      "-->",
      edge.end,
      points[points.length - 1],
      head,
      head.intersect(points[points.length - 1])
    );
    points.push(head.intersect(points[points.length - 1]));
  }
  const pointsStr = btoa(JSON.stringify(points));
  if (edge.toCluster) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.info("to cluster abc88", clusterDb.get(edge.toCluster));
    points = cutPathAtIntersect(edge.points, clusterDb.get(edge.toCluster).node);
    pointsHasChanged = true;
  }
  if (edge.fromCluster) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.debug(
      "from cluster abc88",
      clusterDb.get(edge.fromCluster),
      JSON.stringify(points, null, 2)
    );
    points = cutPathAtIntersect(points.reverse(), clusterDb.get(edge.fromCluster).node).reverse();
    pointsHasChanged = true;
  }
  let lineData = points.filter((p) => !Number.isNaN(p.y));
  lineData = fixCorners(lineData);
  let curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveBasis */ .$0Z;
  curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveLinear */ .c_6;
  switch (edge.curve) {
    case "linear":
      curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveLinear */ .c_6;
      break;
    case "basis":
      curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveBasis */ .$0Z;
      break;
    case "cardinal":
      curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveCardinal */ .YY7;
      break;
    case "bumpX":
      curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveBumpX */ .qpX;
      break;
    case "bumpY":
      curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveBumpY */ .u93;
      break;
    case "catmullRom":
      curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveCatmullRom */ .zgE;
      break;
    case "monotoneX":
      curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveMonotoneX */ .FdL;
      break;
    case "monotoneY":
      curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveMonotoneY */ .ak_;
      break;
    case "natural":
      curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveNatural */ .SxZ;
      break;
    case "step":
      curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveStep */ .eA_;
      break;
    case "stepAfter":
      curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveStepAfter */ .jsv;
      break;
    case "stepBefore":
      curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveStepBefore */ .iJ;
      break;
    default:
      curve = d3__WEBPACK_IMPORTED_MODULE_8__/* .curveBasis */ .$0Z;
  }
  const { x, y } = (0,_chunk_HN2XXSSU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getLineFunctionsWithOffset */ .om)(edge);
  const lineFunction = (0,d3__WEBPACK_IMPORTED_MODULE_8__/* .line */ .jvg)().x(x).y(y).curve(curve);
  let strokeClasses;
  switch (edge.thickness) {
    case "normal":
      strokeClasses = "edge-thickness-normal";
      break;
    case "thick":
      strokeClasses = "edge-thickness-thick";
      break;
    case "invisible":
      strokeClasses = "edge-thickness-invisible";
      break;
    default:
      strokeClasses = "edge-thickness-normal";
  }
  switch (edge.pattern) {
    case "solid":
      strokeClasses += " edge-pattern-solid";
      break;
    case "dotted":
      strokeClasses += " edge-pattern-dotted";
      break;
    case "dashed":
      strokeClasses += " edge-pattern-dashed";
      break;
    default:
      strokeClasses += " edge-pattern-solid";
  }
  let svgPath;
  let linePath = edge.curve === "rounded" ? generateRoundedPath(applyMarkerOffsetsToPoints(lineData, edge), 5) : lineFunction(lineData);
  const edgeStyles = Array.isArray(edge.style) ? edge.style : [edge.style];
  let strokeColor = edgeStyles.find((style) => style?.startsWith("stroke:"));
  let animatedEdge = false;
  if (edge.look === "handDrawn") {
    const rc = roughjs__WEBPACK_IMPORTED_MODULE_9__/* ["default"] */ .Z.svg(elem);
    Object.assign([], lineData);
    const svgPathNode = rc.path(linePath, {
      roughness: 0.3,
      seed: handDrawnSeed
    });
    strokeClasses += " transition";
    svgPath = (0,d3__WEBPACK_IMPORTED_MODULE_8__/* .select */ .Ys)(svgPathNode).select("path").attr("id", edge.id).attr("class", " " + strokeClasses + (edge.classes ? " " + edge.classes : "")).attr("style", edgeStyles ? edgeStyles.reduce((acc, style) => acc + ";" + style, "") : "");
    let d = svgPath.attr("d");
    svgPath.attr("d", d);
    elem.node().appendChild(svgPath.node());
  } else {
    const stylesFromClasses = edgeClassStyles.join(";");
    const styles = edgeStyles ? edgeStyles.reduce((acc, style) => acc + style + ";", "") : "";
    let animationClass = "";
    if (edge.animate) {
      animationClass = " edge-animation-fast";
    }
    if (edge.animation) {
      animationClass = " edge-animation-" + edge.animation;
    }
    const pathStyle = (stylesFromClasses ? stylesFromClasses + ";" + styles + ";" : styles) + ";" + (edgeStyles ? edgeStyles.reduce((acc, style) => acc + ";" + style, "") : "");
    svgPath = elem.append("path").attr("d", linePath).attr("id", edge.id).attr(
      "class",
      " " + strokeClasses + (edge.classes ? " " + edge.classes : "") + (animationClass ?? "")
    ).attr("style", pathStyle);
    strokeColor = pathStyle.match(/stroke:([^;]+)/)?.[1];
    animatedEdge = edge.animate === true || !!edge.animation || stylesFromClasses.includes("animation");
    const pathNode = svgPath.node();
    const len = typeof pathNode.getTotalLength === "function" ? pathNode.getTotalLength() : 0;
    const oValueS = _chunk_HN2XXSSU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .markerOffsets2 */ .oZ[edge.arrowTypeStart] || 0;
    const oValueE = _chunk_HN2XXSSU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .markerOffsets2 */ .oZ[edge.arrowTypeEnd] || 0;
    if (edge.look === "neo" && !animatedEdge) {
      const dashArray = edge.pattern === "dotted" || edge.pattern === "dashed" ? generateDashArray(len, oValueS, oValueE) : `0 ${oValueS} ${len - oValueS - oValueE} ${oValueE}`;
      const mOffset = `stroke-dasharray: ${dashArray}; stroke-dashoffset: 0;`;
      svgPath.attr("style", mOffset + svgPath.attr("style"));
    }
  }
  svgPath.attr("data-edge", true);
  svgPath.attr("data-et", "edge");
  svgPath.attr("data-id", edge.id);
  svgPath.attr("data-points", pointsStr);
  if (edge.showPoints) {
    lineData.forEach((point3) => {
      elem.append("circle").style("stroke", "red").style("fill", "red").attr("r", 1).attr("cx", point3.x).attr("cy", point3.y);
    });
  }
  let url = "";
  if ((0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_6__/* .getConfig2 */ .nV)().flowchart.arrowMarkerAbsolute || (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_6__/* .getConfig2 */ .nV)().state.arrowMarkerAbsolute) {
    url = window.location.protocol + "//" + window.location.host + window.location.pathname + window.location.search;
    url = url.replace(/\(/g, "\\(").replace(/\)/g, "\\)");
  }
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.info("arrowTypeStart", edge.arrowTypeStart);
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.info("arrowTypeEnd", edge.arrowTypeEnd);
  addEdgeMarkers(svgPath, edge, url, id, diagramType, strokeColor);
  const midIndex = Math.floor(points.length / 2);
  const point2 = points[midIndex];
  if (!_chunk_S3R3BYOJ_mjs__WEBPACK_IMPORTED_MODULE_5__/* .utils_default */ .w8.isLabelCoordinateInPath(point2, svgPath.attr("d"))) {
    pointsHasChanged = true;
  }
  let paths = {};
  if (pointsHasChanged) {
    paths.updatedPath = points;
  }
  paths.originalPath = edge.points;
  return paths;
}, "insertEdge");
function generateRoundedPath(points, radius) {
  if (points.length < 2) {
    return "";
  }
  let path = "";
  const size = points.length;
  const epsilon = 1e-5;
  for (let i = 0; i < size; i++) {
    const currPoint = points[i];
    const prevPoint = points[i - 1];
    const nextPoint = points[i + 1];
    if (i === 0) {
      path += `M${currPoint.x},${currPoint.y}`;
    } else if (i === size - 1) {
      path += `L${currPoint.x},${currPoint.y}`;
    } else {
      const dx1 = currPoint.x - prevPoint.x;
      const dy1 = currPoint.y - prevPoint.y;
      const dx2 = nextPoint.x - currPoint.x;
      const dy2 = nextPoint.y - currPoint.y;
      const len1 = Math.hypot(dx1, dy1);
      const len2 = Math.hypot(dx2, dy2);
      if (len1 < epsilon || len2 < epsilon) {
        path += `L${currPoint.x},${currPoint.y}`;
        continue;
      }
      const nx1 = dx1 / len1;
      const ny1 = dy1 / len1;
      const nx2 = dx2 / len2;
      const ny2 = dy2 / len2;
      const dot = nx1 * nx2 + ny1 * ny2;
      const clampedDot = Math.max(-1, Math.min(1, dot));
      const angle = Math.acos(clampedDot);
      if (angle < epsilon || Math.abs(Math.PI - angle) < epsilon) {
        path += `L${currPoint.x},${currPoint.y}`;
        continue;
      }
      const cutLen = Math.min(radius / Math.sin(angle / 2), len1 / 2, len2 / 2);
      const startX = currPoint.x - nx1 * cutLen;
      const startY = currPoint.y - ny1 * cutLen;
      const endX = currPoint.x + nx2 * cutLen;
      const endY = currPoint.y + ny2 * cutLen;
      path += `L${startX},${startY}`;
      path += `Q${currPoint.x},${currPoint.y} ${endX},${endY}`;
    }
  }
  return path;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)(generateRoundedPath, "generateRoundedPath");
function calculateDeltaAndAngle(point1, point2) {
  if (!point1 || !point2) {
    return { angle: 0, deltaX: 0, deltaY: 0 };
  }
  const deltaX = point2.x - point1.x;
  const deltaY = point2.y - point1.y;
  const angle = Math.atan2(deltaY, deltaX);
  return { angle, deltaX, deltaY };
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)(calculateDeltaAndAngle, "calculateDeltaAndAngle");
function applyMarkerOffsetsToPoints(points, edge) {
  const newPoints = points.map((point2) => ({ ...point2 }));
  if (points.length >= 2 && _chunk_HN2XXSSU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .markerOffsets */ .Os[edge.arrowTypeStart]) {
    const offsetValue = _chunk_HN2XXSSU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .markerOffsets */ .Os[edge.arrowTypeStart];
    const point1 = points[0];
    const point2 = points[1];
    const { angle } = calculateDeltaAndAngle(point1, point2);
    const offsetX = offsetValue * Math.cos(angle);
    const offsetY = offsetValue * Math.sin(angle);
    newPoints[0].x = point1.x + offsetX;
    newPoints[0].y = point1.y + offsetY;
  }
  const n = points.length;
  if (n >= 2 && _chunk_HN2XXSSU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .markerOffsets */ .Os[edge.arrowTypeEnd]) {
    const offsetValue = _chunk_HN2XXSSU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .markerOffsets */ .Os[edge.arrowTypeEnd];
    const point1 = points[n - 1];
    const point2 = points[n - 2];
    const { angle } = calculateDeltaAndAngle(point2, point1);
    const offsetX = offsetValue * Math.cos(angle);
    const offsetY = offsetValue * Math.sin(angle);
    newPoints[n - 1].x = point1.x - offsetX;
    newPoints[n - 1].y = point1.y - offsetY;
  }
  return newPoints;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)(applyMarkerOffsetsToPoints, "applyMarkerOffsetsToPoints");

// src/rendering-util/rendering-elements/markers.js
var insertMarkers = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, markerArray, type, id) => {
  markerArray.forEach((markerName) => {
    markers[markerName](elem, type, id);
  });
}, "insertMarkers");
var extension = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .log */ .cM.trace("Making markers for ", id);
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-extensionStart").attr("class", "marker extension " + type).attr("refX", 18).attr("refY", 7).attr("markerWidth", 190).attr("markerHeight", 240).attr("orient", "auto").append("path").attr("d", "M 1,7 L18,13 V 1 Z");
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-extensionEnd").attr("class", "marker extension " + type).attr("refX", 1).attr("refY", 7).attr("markerWidth", 20).attr("markerHeight", 28).attr("orient", "auto").append("path").attr("d", "M 1,1 V 13 L18,7 Z");
}, "extension");
var composition = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-compositionStart").attr("class", "marker composition " + type).attr("refX", 18).attr("refY", 7).attr("markerWidth", 190).attr("markerHeight", 240).attr("orient", "auto").append("path").attr("d", "M 18,7 L9,13 L1,7 L9,1 Z");
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-compositionEnd").attr("class", "marker composition " + type).attr("refX", 1).attr("refY", 7).attr("markerWidth", 20).attr("markerHeight", 28).attr("orient", "auto").append("path").attr("d", "M 18,7 L9,13 L1,7 L9,1 Z");
}, "composition");
var aggregation = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-aggregationStart").attr("class", "marker aggregation " + type).attr("refX", 18).attr("refY", 7).attr("markerWidth", 190).attr("markerHeight", 240).attr("orient", "auto").append("path").attr("d", "M 18,7 L9,13 L1,7 L9,1 Z");
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-aggregationEnd").attr("class", "marker aggregation " + type).attr("refX", 1).attr("refY", 7).attr("markerWidth", 20).attr("markerHeight", 28).attr("orient", "auto").append("path").attr("d", "M 18,7 L9,13 L1,7 L9,1 Z");
}, "aggregation");
var dependency = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-dependencyStart").attr("class", "marker dependency " + type).attr("refX", 6).attr("refY", 7).attr("markerWidth", 190).attr("markerHeight", 240).attr("orient", "auto").append("path").attr("d", "M 5,7 L9,13 L1,7 L9,1 Z");
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-dependencyEnd").attr("class", "marker dependency " + type).attr("refX", 13).attr("refY", 7).attr("markerWidth", 20).attr("markerHeight", 28).attr("orient", "auto").append("path").attr("d", "M 18,7 L9,13 L14,7 L9,1 Z");
}, "dependency");
var lollipop = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-lollipopStart").attr("class", "marker lollipop " + type).attr("refX", 13).attr("refY", 7).attr("markerWidth", 190).attr("markerHeight", 240).attr("orient", "auto").append("circle").attr("stroke", "black").attr("fill", "transparent").attr("cx", 7).attr("cy", 7).attr("r", 6);
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-lollipopEnd").attr("class", "marker lollipop " + type).attr("refX", 1).attr("refY", 7).attr("markerWidth", 190).attr("markerHeight", 240).attr("orient", "auto").append("circle").attr("stroke", "black").attr("fill", "transparent").attr("cx", 7).attr("cy", 7).attr("r", 6);
}, "lollipop");
var point = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  elem.append("marker").attr("id", id + "_" + type + "-pointEnd").attr("class", "marker " + type).attr("viewBox", "0 0 10 10").attr("refX", 5).attr("refY", 5).attr("markerUnits", "userSpaceOnUse").attr("markerWidth", 8).attr("markerHeight", 8).attr("orient", "auto").append("path").attr("d", "M 0 0 L 10 5 L 0 10 z").attr("class", "arrowMarkerPath").style("stroke-width", 1).style("stroke-dasharray", "1,0");
  elem.append("marker").attr("id", id + "_" + type + "-pointStart").attr("class", "marker " + type).attr("viewBox", "0 0 10 10").attr("refX", 4.5).attr("refY", 5).attr("markerUnits", "userSpaceOnUse").attr("markerWidth", 8).attr("markerHeight", 8).attr("orient", "auto").append("path").attr("d", "M 0 5 L 10 10 L 10 0 z").attr("class", "arrowMarkerPath").style("stroke-width", 1).style("stroke-dasharray", "1,0");
}, "point");
var circle = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  elem.append("marker").attr("id", id + "_" + type + "-circleEnd").attr("class", "marker " + type).attr("viewBox", "0 0 10 10").attr("refX", 11).attr("refY", 5).attr("markerUnits", "userSpaceOnUse").attr("markerWidth", 11).attr("markerHeight", 11).attr("orient", "auto").append("circle").attr("cx", "5").attr("cy", "5").attr("r", "5").attr("class", "arrowMarkerPath").style("stroke-width", 1).style("stroke-dasharray", "1,0");
  elem.append("marker").attr("id", id + "_" + type + "-circleStart").attr("class", "marker " + type).attr("viewBox", "0 0 10 10").attr("refX", -1).attr("refY", 5).attr("markerUnits", "userSpaceOnUse").attr("markerWidth", 11).attr("markerHeight", 11).attr("orient", "auto").append("circle").attr("cx", "5").attr("cy", "5").attr("r", "5").attr("class", "arrowMarkerPath").style("stroke-width", 1).style("stroke-dasharray", "1,0");
}, "circle");
var cross = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  elem.append("marker").attr("id", id + "_" + type + "-crossEnd").attr("class", "marker cross " + type).attr("viewBox", "0 0 11 11").attr("refX", 12).attr("refY", 5.2).attr("markerUnits", "userSpaceOnUse").attr("markerWidth", 11).attr("markerHeight", 11).attr("orient", "auto").append("path").attr("d", "M 1,1 l 9,9 M 10,1 l -9,9").attr("class", "arrowMarkerPath").style("stroke-width", 2).style("stroke-dasharray", "1,0");
  elem.append("marker").attr("id", id + "_" + type + "-crossStart").attr("class", "marker cross " + type).attr("viewBox", "0 0 11 11").attr("refX", -1).attr("refY", 5.2).attr("markerUnits", "userSpaceOnUse").attr("markerWidth", 11).attr("markerHeight", 11).attr("orient", "auto").append("path").attr("d", "M 1,1 l 9,9 M 10,1 l -9,9").attr("class", "arrowMarkerPath").style("stroke-width", 2).style("stroke-dasharray", "1,0");
}, "cross");
var barb = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-barbEnd").attr("refX", 19).attr("refY", 7).attr("markerWidth", 20).attr("markerHeight", 14).attr("markerUnits", "userSpaceOnUse").attr("orient", "auto").append("path").attr("d", "M 19,7 L9,13 L14,7 L9,1 Z");
}, "barb");
var only_one = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-onlyOneStart").attr("class", "marker onlyOne " + type).attr("refX", 0).attr("refY", 9).attr("markerWidth", 18).attr("markerHeight", 18).attr("orient", "auto").append("path").attr("d", "M9,0 L9,18 M15,0 L15,18");
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-onlyOneEnd").attr("class", "marker onlyOne " + type).attr("refX", 18).attr("refY", 9).attr("markerWidth", 18).attr("markerHeight", 18).attr("orient", "auto").append("path").attr("d", "M3,0 L3,18 M9,0 L9,18");
}, "only_one");
var zero_or_one = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  const startMarker = elem.append("defs").append("marker").attr("id", id + "_" + type + "-zeroOrOneStart").attr("class", "marker zeroOrOne " + type).attr("refX", 0).attr("refY", 9).attr("markerWidth", 30).attr("markerHeight", 18).attr("orient", "auto");
  startMarker.append("circle").attr("fill", "white").attr("cx", 21).attr("cy", 9).attr("r", 6);
  startMarker.append("path").attr("d", "M9,0 L9,18");
  const endMarker = elem.append("defs").append("marker").attr("id", id + "_" + type + "-zeroOrOneEnd").attr("class", "marker zeroOrOne " + type).attr("refX", 30).attr("refY", 9).attr("markerWidth", 30).attr("markerHeight", 18).attr("orient", "auto");
  endMarker.append("circle").attr("fill", "white").attr("cx", 9).attr("cy", 9).attr("r", 6);
  endMarker.append("path").attr("d", "M21,0 L21,18");
}, "zero_or_one");
var one_or_more = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-oneOrMoreStart").attr("class", "marker oneOrMore " + type).attr("refX", 18).attr("refY", 18).attr("markerWidth", 45).attr("markerHeight", 36).attr("orient", "auto").append("path").attr("d", "M0,18 Q 18,0 36,18 Q 18,36 0,18 M42,9 L42,27");
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-oneOrMoreEnd").attr("class", "marker oneOrMore " + type).attr("refX", 27).attr("refY", 18).attr("markerWidth", 45).attr("markerHeight", 36).attr("orient", "auto").append("path").attr("d", "M3,9 L3,27 M9,18 Q27,0 45,18 Q27,36 9,18");
}, "one_or_more");
var zero_or_more = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  const startMarker = elem.append("defs").append("marker").attr("id", id + "_" + type + "-zeroOrMoreStart").attr("class", "marker zeroOrMore " + type).attr("refX", 18).attr("refY", 18).attr("markerWidth", 57).attr("markerHeight", 36).attr("orient", "auto");
  startMarker.append("circle").attr("fill", "white").attr("cx", 48).attr("cy", 18).attr("r", 6);
  startMarker.append("path").attr("d", "M0,18 Q18,0 36,18 Q18,36 0,18");
  const endMarker = elem.append("defs").append("marker").attr("id", id + "_" + type + "-zeroOrMoreEnd").attr("class", "marker zeroOrMore " + type).attr("refX", 39).attr("refY", 18).attr("markerWidth", 57).attr("markerHeight", 36).attr("orient", "auto");
  endMarker.append("circle").attr("fill", "white").attr("cx", 9).attr("cy", 18).attr("r", 6);
  endMarker.append("path").attr("d", "M21,18 Q39,0 57,18 Q39,36 21,18");
}, "zero_or_more");
var requirement_arrow = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  elem.append("defs").append("marker").attr("id", id + "_" + type + "-requirement_arrowEnd").attr("refX", 20).attr("refY", 10).attr("markerWidth", 20).attr("markerHeight", 20).attr("orient", "auto").append("path").attr(
    "d",
    `M0,0
      L20,10
      M20,10
      L0,20`
  );
}, "requirement_arrow");
var requirement_contains = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_7__/* .__name */ .eW)((elem, type, id) => {
  const containsNode = elem.append("defs").append("marker").attr("id", id + "_" + type + "-requirement_containsStart").attr("refX", 0).attr("refY", 10).attr("markerWidth", 20).attr("markerHeight", 20).attr("orient", "auto").append("g");
  containsNode.append("circle").attr("cx", 10).attr("cy", 10).attr("r", 9).attr("fill", "none");
  containsNode.append("line").attr("x1", 1).attr("x2", 19).attr("y1", 10).attr("y2", 10);
  containsNode.append("line").attr("y1", 1).attr("y2", 19).attr("x1", 10).attr("x2", 10);
}, "requirement_contains");
var markers = {
  extension,
  composition,
  aggregation,
  dependency,
  lollipop,
  point,
  circle,
  cross,
  barb,
  only_one,
  zero_or_one,
  one_or_more,
  zero_or_more,
  requirement_arrow,
  requirement_contains
};
var markers_default = insertMarkers;




/***/ }),

/***/ 17175:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   $m: () => (/* binding */ ZERO_WIDTH_SPACE),
/* harmony export */   Cq: () => (/* binding */ calculateTextWidth),
/* harmony export */   Ln: () => (/* binding */ getEdgeId),
/* harmony export */   MX: () => (/* binding */ random),
/* harmony export */   Ox: () => (/* binding */ generateId),
/* harmony export */   R7: () => (/* binding */ handleUndefinedAttr),
/* harmony export */   Rb: () => (/* binding */ cleanAndMerge),
/* harmony export */   SH: () => (/* binding */ decodeEntities),
/* harmony export */   VG: () => (/* binding */ parseFontSize),
/* harmony export */   Vy: () => (/* binding */ encodeEntities),
/* harmony export */   X4: () => (/* binding */ wrapLabel),
/* harmony export */   XD: () => (/* binding */ calculateTextHeight),
/* harmony export */   bZ: () => (/* binding */ isDetailedError),
/* harmony export */   be: () => (/* binding */ getStylesFromArray),
/* harmony export */   le: () => (/* binding */ interpolateToCurve),
/* harmony export */   tf: () => (/* binding */ removeDirectives),
/* harmony export */   w8: () => (/* binding */ utils_default)
/* harmony export */ });
/* harmony import */ var _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(61805);
/* harmony import */ var _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(74999);
/* harmony import */ var _braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(7608);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(35321);
/* harmony import */ var lodash_es_memoize_js__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(86861);
/* harmony import */ var lodash_es_merge_js__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(20895);



// src/utils.ts




var ZERO_WIDTH_SPACE = "\u200B";
var d3CurveTypes = {
  curveBasis: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveBasis */ .$0Z,
  curveBasisClosed: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveBasisClosed */ .Dts,
  curveBasisOpen: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveBasisOpen */ .WQY,
  curveBumpX: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveBumpX */ .qpX,
  curveBumpY: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveBumpY */ .u93,
  curveBundle: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveBundle */ .tFB,
  curveCardinalClosed: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveCardinalClosed */ .OvA,
  curveCardinalOpen: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveCardinalOpen */ .dCK,
  curveCardinal: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveCardinal */ .YY7,
  curveCatmullRomClosed: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveCatmullRomClosed */ .fGX,
  curveCatmullRomOpen: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveCatmullRomOpen */ .$m7,
  curveCatmullRom: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveCatmullRom */ .zgE,
  curveLinear: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveLinear */ .c_6,
  curveLinearClosed: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveLinearClosed */ .fxm,
  curveMonotoneX: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveMonotoneX */ .FdL,
  curveMonotoneY: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveMonotoneY */ .ak_,
  curveNatural: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveNatural */ .SxZ,
  curveStep: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveStep */ .eA_,
  curveStepAfter: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveStepAfter */ .jsv,
  curveStepBefore: d3__WEBPACK_IMPORTED_MODULE_3__/* .curveStepBefore */ .iJ
};
var directiveWithoutOpen = /\s*(?:(\w+)(?=:):|(\w+))\s*(?:(\w+)|((?:(?!}%{2}).|\r?\n)*))?\s*(?:}%{2})?/gi;
var detectInit = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(function(text, config) {
  const inits = detectDirective(text, /(?:init\b)|(?:initialize\b)/);
  let results = {};
  if (Array.isArray(inits)) {
    const args = inits.map((init) => init.args);
    (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .sanitizeDirective */ .NM)(args);
    results = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .assignWithDepth_default */ .Yc)(results, [...args]);
  } else {
    results = inits.args;
  }
  if (!results) {
    return;
  }
  let type = (0,_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .detectType */ .Vg)(text, config);
  const prop = "config";
  if (results[prop] !== void 0) {
    if (type === "flowchart-v2") {
      type = "flowchart";
    }
    results[type] = results[prop];
    delete results[prop];
  }
  return results;
}, "detectInit");
var detectDirective = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(function(text, type = null) {
  try {
    const commentWithoutDirectives = new RegExp(
      `[%]{2}(?![{]${directiveWithoutOpen.source})(?=[}][%]{2}).*
`,
      "ig"
    );
    text = text.trim().replace(commentWithoutDirectives, "").replace(/'/gm, '"');
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .log */ .cM.debug(
      `Detecting diagram directive${type !== null ? " type:" + type : ""} based on the text:${text}`
    );
    let match;
    const result = [];
    while ((match = _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .directiveRegex */ .Zn.exec(text)) !== null) {
      if (match.index === _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .directiveRegex */ .Zn.lastIndex) {
        _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .directiveRegex */ .Zn.lastIndex++;
      }
      if (match && !type || type && match[1]?.match(type) || type && match[2]?.match(type)) {
        const type2 = match[1] ? match[1] : match[2];
        const args = match[3] ? match[3].trim() : match[4] ? JSON.parse(match[4].trim()) : null;
        result.push({ type: type2, args });
      }
    }
    if (result.length === 0) {
      return { type: text, args: null };
    }
    return result.length === 1 ? result[0] : result;
  } catch (error) {
    _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .log */ .cM.error(
      `ERROR: ${error.message} - Unable to parse directive type: '${type}' based on the text: '${text}'`
    );
    return { type: void 0, args: null };
  }
}, "detectDirective");
var removeDirectives = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(function(text) {
  return text.replace(_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .directiveRegex */ .Zn, "");
}, "removeDirectives");
var isSubstringInArray = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(function(str, arr) {
  for (const [i, element] of arr.entries()) {
    if (element.match(str)) {
      return i;
    }
  }
  return -1;
}, "isSubstringInArray");
function interpolateToCurve(interpolate, defaultCurve) {
  if (!interpolate) {
    return defaultCurve;
  }
  const curveName = `curve${interpolate.charAt(0).toUpperCase() + interpolate.slice(1)}`;
  return d3CurveTypes[curveName] ?? defaultCurve;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(interpolateToCurve, "interpolateToCurve");
function formatUrl(linkStr, config) {
  const url = linkStr.trim();
  if (!url) {
    return void 0;
  }
  if (config.securityLevel !== "loose") {
    return (0,_braintree_sanitize_url__WEBPACK_IMPORTED_MODULE_2__/* .sanitizeUrl */ .N)(url);
  }
  return url;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(formatUrl, "formatUrl");
var runFunc = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((functionName, ...params) => {
  const arrPaths = functionName.split(".");
  const len = arrPaths.length - 1;
  const fnName = arrPaths[len];
  let obj = window;
  for (let i = 0; i < len; i++) {
    obj = obj[arrPaths[i]];
    if (!obj) {
      _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .log */ .cM.error(`Function name: ${functionName} not found in window`);
      return;
    }
  }
  obj[fnName](...params);
}, "runFunc");
function distance(p1, p2) {
  if (!p1 || !p2) {
    return 0;
  }
  return Math.sqrt(Math.pow(p2.x - p1.x, 2) + Math.pow(p2.y - p1.y, 2));
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(distance, "distance");
function traverseEdge(points) {
  let prevPoint;
  let totalDistance = 0;
  points.forEach((point) => {
    totalDistance += distance(point, prevPoint);
    prevPoint = point;
  });
  const remainingDistance = totalDistance / 2;
  return calculatePoint(points, remainingDistance);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(traverseEdge, "traverseEdge");
function calcLabelPosition(points) {
  if (points.length === 1) {
    return points[0];
  }
  return traverseEdge(points);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(calcLabelPosition, "calcLabelPosition");
var roundNumber = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((num, precision = 2) => {
  const factor = Math.pow(10, precision);
  return Math.round(num * factor) / factor;
}, "roundNumber");
var calculatePoint = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((points, distanceToTraverse) => {
  let prevPoint = void 0;
  let remainingDistance = distanceToTraverse;
  for (const point of points) {
    if (prevPoint) {
      const vectorDistance = distance(point, prevPoint);
      if (vectorDistance === 0) {
        return prevPoint;
      }
      if (vectorDistance < remainingDistance) {
        remainingDistance -= vectorDistance;
      } else {
        const distanceRatio = remainingDistance / vectorDistance;
        if (distanceRatio <= 0) {
          return prevPoint;
        }
        if (distanceRatio >= 1) {
          return { x: point.x, y: point.y };
        }
        if (distanceRatio > 0 && distanceRatio < 1) {
          return {
            x: roundNumber((1 - distanceRatio) * prevPoint.x + distanceRatio * point.x, 5),
            y: roundNumber((1 - distanceRatio) * prevPoint.y + distanceRatio * point.y, 5)
          };
        }
      }
    }
    prevPoint = point;
  }
  throw new Error("Could not find a suitable point for the given distance");
}, "calculatePoint");
var calcCardinalityPosition = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((isRelationTypePresent, points, initialPosition) => {
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .log */ .cM.info(`our points ${JSON.stringify(points)}`);
  if (points[0] !== initialPosition) {
    points = points.reverse();
  }
  const distanceToCardinalityPoint = 25;
  const center = calculatePoint(points, distanceToCardinalityPoint);
  const d = isRelationTypePresent ? 10 : 5;
  const angle = Math.atan2(points[0].y - center.y, points[0].x - center.x);
  const cardinalityPosition = { x: 0, y: 0 };
  cardinalityPosition.x = Math.sin(angle) * d + (points[0].x + center.x) / 2;
  cardinalityPosition.y = -Math.cos(angle) * d + (points[0].y + center.y) / 2;
  return cardinalityPosition;
}, "calcCardinalityPosition");
function calcTerminalLabelPosition(terminalMarkerSize, position, _points) {
  const points = structuredClone(_points);
  _chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .log */ .cM.info("our points", points);
  if (position !== "start_left" && position !== "start_right") {
    points.reverse();
  }
  const distanceToCardinalityPoint = 25 + terminalMarkerSize;
  const center = calculatePoint(points, distanceToCardinalityPoint);
  const d = 10 + terminalMarkerSize * 0.5;
  const angle = Math.atan2(points[0].y - center.y, points[0].x - center.x);
  const cardinalityPosition = { x: 0, y: 0 };
  if (position === "start_left") {
    cardinalityPosition.x = Math.sin(angle + Math.PI) * d + (points[0].x + center.x) / 2;
    cardinalityPosition.y = -Math.cos(angle + Math.PI) * d + (points[0].y + center.y) / 2;
  } else if (position === "end_right") {
    cardinalityPosition.x = Math.sin(angle - Math.PI) * d + (points[0].x + center.x) / 2 - 5;
    cardinalityPosition.y = -Math.cos(angle - Math.PI) * d + (points[0].y + center.y) / 2 - 5;
  } else if (position === "end_left") {
    cardinalityPosition.x = Math.sin(angle) * d + (points[0].x + center.x) / 2 - 5;
    cardinalityPosition.y = -Math.cos(angle) * d + (points[0].y + center.y) / 2 - 5;
  } else {
    cardinalityPosition.x = Math.sin(angle) * d + (points[0].x + center.x) / 2;
    cardinalityPosition.y = -Math.cos(angle) * d + (points[0].y + center.y) / 2;
  }
  return cardinalityPosition;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(calcTerminalLabelPosition, "calcTerminalLabelPosition");
function getStylesFromArray(arr) {
  let style = "";
  let labelStyle = "";
  for (const element of arr) {
    if (element !== void 0) {
      if (element.startsWith("color:") || element.startsWith("text-align:")) {
        labelStyle = labelStyle + element + ";";
      } else {
        style = style + element + ";";
      }
    }
  }
  return { style, labelStyle };
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(getStylesFromArray, "getStylesFromArray");
var cnt = 0;
var generateId = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(() => {
  cnt++;
  return "id-" + Math.random().toString(36).substr(2, 12) + "-" + cnt;
}, "generateId");
function makeRandomHex(length) {
  let result = "";
  const characters = "0123456789abcdef";
  const charactersLength = characters.length;
  for (let i = 0; i < length; i++) {
    result += characters.charAt(Math.floor(Math.random() * charactersLength));
  }
  return result;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(makeRandomHex, "makeRandomHex");
var random = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((options) => {
  return makeRandomHex(options.length);
}, "random");
var getTextObj = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(function() {
  return {
    x: 0,
    y: 0,
    fill: void 0,
    anchor: "start",
    style: "#666",
    width: 100,
    height: 100,
    textMargin: 0,
    rx: 0,
    ry: 0,
    valign: void 0,
    text: ""
  };
}, "getTextObj");
var drawSimpleText = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(function(elem, textData) {
  const nText = textData.text.replace(_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .common_default */ .SY.lineBreakRegex, " ");
  const [, _fontSizePx] = parseFontSize(textData.fontSize);
  const textElem = elem.append("text");
  textElem.attr("x", textData.x);
  textElem.attr("y", textData.y);
  textElem.style("text-anchor", textData.anchor);
  textElem.style("font-family", textData.fontFamily);
  textElem.style("font-size", _fontSizePx);
  textElem.style("font-weight", textData.fontWeight);
  textElem.attr("fill", textData.fill);
  if (textData.class !== void 0) {
    textElem.attr("class", textData.class);
  }
  const span = textElem.append("tspan");
  span.attr("x", textData.x + textData.textMargin * 2);
  span.attr("fill", textData.fill);
  span.text(nText);
  return textElem;
}, "drawSimpleText");
var wrapLabel = (0,lodash_es_memoize_js__WEBPACK_IMPORTED_MODULE_4__/* ["default"] */ .Z)(
  (label, maxWidth, config) => {
    if (!label) {
      return label;
    }
    config = Object.assign(
      { fontSize: 12, fontWeight: 400, fontFamily: "Arial", joinWith: "<br/>" },
      config
    );
    if (_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .common_default */ .SY.lineBreakRegex.test(label)) {
      return label;
    }
    const words = label.split(" ").filter(Boolean);
    const completedLines = [];
    let nextLine = "";
    words.forEach((word, index) => {
      const wordLength = calculateTextWidth(`${word} `, config);
      const nextLineLength = calculateTextWidth(nextLine, config);
      if (wordLength > maxWidth) {
        const { hyphenatedStrings, remainingWord } = breakString(word, maxWidth, "-", config);
        completedLines.push(nextLine, ...hyphenatedStrings);
        nextLine = remainingWord;
      } else if (nextLineLength + wordLength >= maxWidth) {
        completedLines.push(nextLine);
        nextLine = word;
      } else {
        nextLine = [nextLine, word].filter(Boolean).join(" ");
      }
      const currentWord = index + 1;
      const isLastWord = currentWord === words.length;
      if (isLastWord) {
        completedLines.push(nextLine);
      }
    });
    return completedLines.filter((line) => line !== "").join(config.joinWith);
  },
  (label, maxWidth, config) => `${label}${maxWidth}${config.fontSize}${config.fontWeight}${config.fontFamily}${config.joinWith}`
);
var breakString = (0,lodash_es_memoize_js__WEBPACK_IMPORTED_MODULE_4__/* ["default"] */ .Z)(
  (word, maxWidth, hyphenCharacter = "-", config) => {
    config = Object.assign(
      { fontSize: 12, fontWeight: 400, fontFamily: "Arial", margin: 0 },
      config
    );
    const characters = [...word];
    const lines = [];
    let currentLine = "";
    characters.forEach((character, index) => {
      const nextLine = `${currentLine}${character}`;
      const lineWidth = calculateTextWidth(nextLine, config);
      if (lineWidth >= maxWidth) {
        const currentCharacter = index + 1;
        const isLastLine = characters.length === currentCharacter;
        const hyphenatedNextLine = `${nextLine}${hyphenCharacter}`;
        lines.push(isLastLine ? nextLine : hyphenatedNextLine);
        currentLine = "";
      } else {
        currentLine = nextLine;
      }
    });
    return { hyphenatedStrings: lines, remainingWord: currentLine };
  },
  (word, maxWidth, hyphenCharacter = "-", config) => `${word}${maxWidth}${hyphenCharacter}${config.fontSize}${config.fontWeight}${config.fontFamily}`
);
function calculateTextHeight(text, config) {
  return calculateTextDimensions(text, config).height;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(calculateTextHeight, "calculateTextHeight");
function calculateTextWidth(text, config) {
  return calculateTextDimensions(text, config).width;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(calculateTextWidth, "calculateTextWidth");
var calculateTextDimensions = (0,lodash_es_memoize_js__WEBPACK_IMPORTED_MODULE_4__/* ["default"] */ .Z)(
  (text, config) => {
    const { fontSize = 12, fontFamily = "Arial", fontWeight = 400 } = config;
    if (!text) {
      return { width: 0, height: 0 };
    }
    const [, _fontSizePx] = parseFontSize(fontSize);
    const fontFamilies = ["sans-serif", fontFamily];
    const lines = text.split(_chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .common_default */ .SY.lineBreakRegex);
    const dims = [];
    const body = (0,d3__WEBPACK_IMPORTED_MODULE_3__/* .select */ .Ys)("body");
    if (!body.remove) {
      return { width: 0, height: 0, lineHeight: 0 };
    }
    const g = body.append("svg");
    for (const fontFamily2 of fontFamilies) {
      let cHeight = 0;
      const dim = { width: 0, height: 0, lineHeight: 0 };
      for (const line of lines) {
        const textObj = getTextObj();
        textObj.text = line || ZERO_WIDTH_SPACE;
        const textElem = drawSimpleText(g, textObj).style("font-size", _fontSizePx).style("font-weight", fontWeight).style("font-family", fontFamily2);
        const bBox = (textElem._groups || textElem)[0][0].getBBox();
        if (bBox.width === 0 && bBox.height === 0) {
          throw new Error("svg element not in render tree");
        }
        dim.width = Math.round(Math.max(dim.width, bBox.width));
        cHeight = Math.round(bBox.height);
        dim.height += cHeight;
        dim.lineHeight = Math.round(Math.max(dim.lineHeight, cHeight));
      }
      dims.push(dim);
    }
    g.remove();
    const index = isNaN(dims[1].height) || isNaN(dims[1].width) || isNaN(dims[1].lineHeight) || dims[0].height > dims[1].height && dims[0].width > dims[1].width && dims[0].lineHeight > dims[1].lineHeight ? 0 : 1;
    return dims[index];
  },
  (text, config) => `${text}${config.fontSize}${config.fontWeight}${config.fontFamily}`
);
var InitIDGenerator = class {
  constructor(deterministic = false, seed) {
    this.count = 0;
    this.count = seed ? seed.length : 0;
    this.next = deterministic ? () => this.count++ : () => Date.now();
  }
  static {
    (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(this, "InitIDGenerator");
  }
};
var decoder;
var entityDecode = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(function(html) {
  decoder = decoder || document.createElement("div");
  html = escape(html).replace(/%26/g, "&").replace(/%23/g, "#").replace(/%3B/g, ";");
  decoder.innerHTML = html;
  return unescape(decoder.textContent);
}, "entityDecode");
function isDetailedError(error) {
  return "str" in error;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(isDetailedError, "isDetailedError");
var insertTitle = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((parent, cssClass, titleTopMargin, title) => {
  if (!title) {
    return;
  }
  const bounds = parent.node()?.getBBox();
  if (!bounds) {
    return;
  }
  parent.append("text").text(title).attr("text-anchor", "middle").attr("x", bounds.x + bounds.width / 2).attr("y", -titleTopMargin).attr("class", cssClass);
}, "insertTitle");
var parseFontSize = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((fontSize) => {
  if (typeof fontSize === "number") {
    return [fontSize, fontSize + "px"];
  }
  const fontSizeNumber = parseInt(fontSize ?? "", 10);
  if (Number.isNaN(fontSizeNumber)) {
    return [void 0, void 0];
  } else if (fontSize === String(fontSizeNumber)) {
    return [fontSizeNumber, fontSize + "px"];
  } else {
    return [fontSizeNumber, fontSize];
  }
}, "parseFontSize");
function cleanAndMerge(defaultData, data) {
  return (0,lodash_es_merge_js__WEBPACK_IMPORTED_MODULE_5__/* ["default"] */ .Z)({}, defaultData, data);
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(cleanAndMerge, "cleanAndMerge");
var utils_default = {
  assignWithDepth: _chunk_ABZYJK2D_mjs__WEBPACK_IMPORTED_MODULE_0__/* .assignWithDepth_default */ .Yc,
  wrapLabel,
  calculateTextHeight,
  calculateTextWidth,
  calculateTextDimensions,
  cleanAndMerge,
  detectInit,
  detectDirective,
  isSubstringInArray,
  interpolateToCurve,
  calcLabelPosition,
  calcCardinalityPosition,
  calcTerminalLabelPosition,
  formatUrl,
  getStylesFromArray,
  generateId,
  random,
  runFunc,
  entityDecode,
  insertTitle,
  isLabelCoordinateInPath,
  parseFontSize,
  InitIDGenerator
};
var encodeEntities = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(function(text) {
  let txt = text;
  txt = txt.replace(/style.*:\S*#.*;/g, function(s) {
    return s.substring(0, s.length - 1);
  });
  txt = txt.replace(/classDef.*:\S*#.*;/g, function(s) {
    return s.substring(0, s.length - 1);
  });
  txt = txt.replace(/#\w+;/g, function(s) {
    const innerTxt = s.substring(1, s.length - 1);
    const isInt = /^\+?\d+$/.test(innerTxt);
    if (isInt) {
      return "\uFB02\xB0\xB0" + innerTxt + "\xB6\xDF";
    } else {
      return "\uFB02\xB0" + innerTxt + "\xB6\xDF";
    }
  });
  return txt;
}, "encodeEntities");
var decodeEntities = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(function(text) {
  return text.replace(/ﬂ°°/g, "&#").replace(/ﬂ°/g, "&").replace(/¶ß/g, ";");
}, "decodeEntities");
var getEdgeId = /* @__PURE__ */ (0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)((from, to, {
  counter = 0,
  prefix,
  suffix
}, id) => {
  if (id) {
    return id;
  }
  return `${prefix ? `${prefix}_` : ""}${from}_${to}_${counter}${suffix ? `_${suffix}` : ""}`;
}, "getEdgeId");
function handleUndefinedAttr(attrValue) {
  return attrValue ?? null;
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(handleUndefinedAttr, "handleUndefinedAttr");
function isLabelCoordinateInPath(point, dAttr) {
  const roundedX = Math.round(point.x);
  const roundedY = Math.round(point.y);
  const sanitizedD = dAttr.replace(
    /(\d+\.\d+)/g,
    (match) => Math.round(parseFloat(match)).toString()
  );
  return sanitizedD.includes(roundedX.toString()) || sanitizedD.includes(roundedY.toString());
}
(0,_chunk_AGHRB4JF_mjs__WEBPACK_IMPORTED_MODULE_1__/* .__name */ .eW)(isLabelCoordinateInPath, "isLabelCoordinateInPath");




/***/ }),

/***/ 94406:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
// ESM COMPAT FLAG
__webpack_require__.r(__webpack_exports__);

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  "default": () => (/* binding */ mermaid_default)
});

// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-6MN3ZHY7.mjs
var chunk_6MN3ZHY7 = __webpack_require__(50692);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-EXTU4WIE.mjs
var chunk_EXTU4WIE = __webpack_require__(98154);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-MI3HLSF2.mjs
var chunk_MI3HLSF2 = __webpack_require__(71262);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-N4CR4FBY.mjs
var chunk_N4CR4FBY = __webpack_require__(54160);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-QXUST7PY.mjs
var chunk_QXUST7PY = __webpack_require__(55195);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-HN2XXSSU.mjs
var chunk_HN2XXSSU = __webpack_require__(67894);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-JZLCHNYA.mjs
var chunk_JZLCHNYA = __webpack_require__(26212);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-CVBHYZKI.mjs
var chunk_CVBHYZKI = __webpack_require__(22462);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-ATLVNIR6.mjs
var chunk_ATLVNIR6 = __webpack_require__(37756);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-JA3XYJ7Z.mjs + 12 modules
var chunk_JA3XYJ7Z = __webpack_require__(2111);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-S3R3BYOJ.mjs
var chunk_S3R3BYOJ = __webpack_require__(17175);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-ABZYJK2D.mjs + 3 modules
var chunk_ABZYJK2D = __webpack_require__(61805);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-AGHRB4JF.mjs
var chunk_AGHRB4JF = __webpack_require__(74999);
// EXTERNAL MODULE: ../node_modules/ts-dedent/esm/index.js
var esm = __webpack_require__(11464);
// EXTERNAL MODULE: ../node_modules/d3/src/index.js + 99 modules
var src = __webpack_require__(35321);
;// CONCATENATED MODULE: ../node_modules/stylis/src/Enum.js
var MS = '-ms-'
var MOZ = '-moz-'
var WEBKIT = '-webkit-'

var COMMENT = 'comm'
var RULESET = 'rule'
var DECLARATION = 'decl'

var PAGE = '@page'
var MEDIA = '@media'
var IMPORT = '@import'
var CHARSET = '@charset'
var VIEWPORT = '@viewport'
var SUPPORTS = '@supports'
var DOCUMENT = '@document'
var NAMESPACE = '@namespace'
var KEYFRAMES = '@keyframes'
var FONT_FACE = '@font-face'
var COUNTER_STYLE = '@counter-style'
var FONT_FEATURE_VALUES = '@font-feature-values'
var LAYER = '@layer'
var SCOPE = '@scope'

;// CONCATENATED MODULE: ../node_modules/stylis/src/Utility.js
/**
 * @param {number}
 * @return {number}
 */
var abs = Math.abs

/**
 * @param {number}
 * @return {string}
 */
var Utility_from = String.fromCharCode

/**
 * @param {object}
 * @return {object}
 */
var Utility_assign = Object.assign

/**
 * @param {string} value
 * @param {number} length
 * @return {number}
 */
function hash (value, length) {
	return charat(value, 0) ^ 45 ? (((((((length << 2) ^ charat(value, 0)) << 2) ^ charat(value, 1)) << 2) ^ charat(value, 2)) << 2) ^ charat(value, 3) : 0
}

/**
 * @param {string} value
 * @return {string}
 */
function trim (value) {
	return value.trim()
}

/**
 * @param {string} value
 * @param {RegExp} pattern
 * @return {string?}
 */
function match (value, pattern) {
	return (value = pattern.exec(value)) ? value[0] : value
}

/**
 * @param {string} value
 * @param {(string|RegExp)} pattern
 * @param {string} replacement
 * @return {string}
 */
function replace (value, pattern, replacement) {
	return value.replace(pattern, replacement)
}

/**
 * @param {string} value
 * @param {string} search
 * @param {number} position
 * @return {number}
 */
function indexof (value, search, position) {
	return value.indexOf(search, position)
}

/**
 * @param {string} value
 * @param {number} index
 * @return {number}
 */
function charat (value, index) {
	return value.charCodeAt(index) | 0
}

/**
 * @param {string} value
 * @param {number} begin
 * @param {number} end
 * @return {string}
 */
function substr (value, begin, end) {
	return value.slice(begin, end)
}

/**
 * @param {string} value
 * @return {number}
 */
function strlen (value) {
	return value.length
}

/**
 * @param {any[]} value
 * @return {number}
 */
function sizeof (value) {
	return value.length
}

/**
 * @param {any} value
 * @param {any[]} array
 * @return {any}
 */
function Utility_append (value, array) {
	return array.push(value), value
}

/**
 * @param {string[]} array
 * @param {function} callback
 * @return {string}
 */
function combine (array, callback) {
	return array.map(callback).join('')
}

/**
 * @param {string[]} array
 * @param {RegExp} pattern
 * @return {string[]}
 */
function filter (array, pattern) {
	return array.filter(function (value) { return !match(value, pattern) })
}

;// CONCATENATED MODULE: ../node_modules/stylis/src/Serializer.js



/**
 * @param {object[]} children
 * @param {function} callback
 * @return {string}
 */
function serialize (children, callback) {
	var output = ''

	for (var i = 0; i < children.length; i++)
		output += callback(children[i], i, children, callback) || ''

	return output
}

/**
 * @param {object} element
 * @param {number} index
 * @param {object[]} children
 * @param {function} callback
 * @return {string}
 */
function stringify (element, index, children, callback) {
	switch (element.type) {
		case LAYER: if (element.children.length) break
		case IMPORT: case NAMESPACE: case DECLARATION: return element.return = element.return || element.value
		case COMMENT: return ''
		case KEYFRAMES: return element.return = element.value + '{' + serialize(element.children, callback) + '}'
		case RULESET: if (!strlen(element.value = element.props.join(','))) return ''
	}

	return strlen(children = serialize(element.children, callback)) ? element.return = element.value + '{' + children + '}' : ''
}

;// CONCATENATED MODULE: ../node_modules/stylis/src/Tokenizer.js


var line = 1
var column = 1
var Tokenizer_length = 0
var position = 0
var character = 0
var characters = ''

/**
 * @param {string} value
 * @param {object | null} root
 * @param {object | null} parent
 * @param {string} type
 * @param {string[] | string} props
 * @param {object[] | string} children
 * @param {object[]} siblings
 * @param {number} length
 */
function node (value, root, parent, type, props, children, length, siblings) {
	return {value: value, root: root, parent: parent, type: type, props: props, children: children, line: line, column: column, length: length, return: '', siblings: siblings}
}

/**
 * @param {object} root
 * @param {object} props
 * @return {object}
 */
function copy (root, props) {
	return assign(node('', null, null, '', null, null, 0, root.siblings), root, {length: -root.length}, props)
}

/**
 * @param {object} root
 */
function lift (root) {
	while (root.root)
		root = copy(root.root, {children: [root]})

	append(root, root.siblings)
}

/**
 * @return {number}
 */
function Tokenizer_char () {
	return character
}

/**
 * @return {number}
 */
function prev () {
	character = position > 0 ? charat(characters, --position) : 0

	if (column--, character === 10)
		column = 1, line--

	return character
}

/**
 * @return {number}
 */
function next () {
	character = position < Tokenizer_length ? charat(characters, position++) : 0

	if (column++, character === 10)
		column = 1, line++

	return character
}

/**
 * @return {number}
 */
function peek () {
	return charat(characters, position)
}

/**
 * @return {number}
 */
function caret () {
	return position
}

/**
 * @param {number} begin
 * @param {number} end
 * @return {string}
 */
function slice (begin, end) {
	return substr(characters, begin, end)
}

/**
 * @param {number} type
 * @return {number}
 */
function token (type) {
	switch (type) {
		// \0 \t \n \r \s whitespace token
		case 0: case 9: case 10: case 13: case 32:
			return 5
		// ! + , / > @ ~ isolate token
		case 33: case 43: case 44: case 47: case 62: case 64: case 126:
		// ; { } breakpoint token
		case 59: case 123: case 125:
			return 4
		// : accompanied token
		case 58:
			return 3
		// " ' ( [ opening delimit token
		case 34: case 39: case 40: case 91:
			return 2
		// ) ] closing delimit token
		case 41: case 93:
			return 1
	}

	return 0
}

/**
 * @param {string} value
 * @return {any[]}
 */
function alloc (value) {
	return line = column = 1, Tokenizer_length = strlen(characters = value), position = 0, []
}

/**
 * @param {any} value
 * @return {any}
 */
function dealloc (value) {
	return characters = '', value
}

/**
 * @param {number} type
 * @return {string}
 */
function delimit (type) {
	return trim(slice(position - 1, delimiter(type === 91 ? type + 2 : type === 40 ? type + 1 : type)))
}

/**
 * @param {string} value
 * @return {string[]}
 */
function tokenize (value) {
	return dealloc(tokenizer(alloc(value)))
}

/**
 * @param {number} type
 * @return {string}
 */
function whitespace (type) {
	while (character = peek())
		if (character < 33)
			next()
		else
			break

	return token(type) > 2 || token(character) > 3 ? '' : ' '
}

/**
 * @param {string[]} children
 * @return {string[]}
 */
function tokenizer (children) {
	while (next())
		switch (token(character)) {
			case 0: append(identifier(position - 1), children)
				break
			case 2: append(delimit(character), children)
				break
			default: append(from(character), children)
		}

	return children
}

/**
 * @param {number} index
 * @param {number} count
 * @return {string}
 */
function escaping (index, count) {
	while (--count && next())
		// not 0-9 A-F a-f
		if (character < 48 || character > 102 || (character > 57 && character < 65) || (character > 70 && character < 97))
			break

	return slice(index, caret() + (count < 6 && peek() == 32 && next() == 32))
}

/**
 * @param {number} type
 * @return {number}
 */
function delimiter (type) {
	while (next())
		switch (character) {
			// ] ) " '
			case type:
				return position
			// " '
			case 34: case 39:
				if (type !== 34 && type !== 39)
					delimiter(character)
				break
			// (
			case 40:
				if (type === 41)
					delimiter(type)
				break
			// \
			case 92:
				next()
				break
		}

	return position
}

/**
 * @param {number} type
 * @param {number} index
 * @return {number}
 */
function commenter (type, index) {
	while (next())
		// //
		if (type + character === 47 + 10)
			break
		// /*
		else if (type + character === 42 + 42 && peek() === 47)
			break

	return '/*' + slice(index, position - 1) + '*' + Utility_from(type === 47 ? type : next())
}

/**
 * @param {number} index
 * @return {string}
 */
function identifier (index) {
	while (!token(peek()))
		next()

	return slice(index, position)
}

;// CONCATENATED MODULE: ../node_modules/stylis/src/Parser.js




/**
 * @param {string} value
 * @return {object[]}
 */
function compile (value) {
	return dealloc(parse('', null, null, null, [''], value = alloc(value), 0, [0], value))
}

/**
 * @param {string} value
 * @param {object} root
 * @param {object?} parent
 * @param {string[]} rule
 * @param {string[]} rules
 * @param {string[]} rulesets
 * @param {number[]} pseudo
 * @param {number[]} points
 * @param {string[]} declarations
 * @return {object}
 */
function parse (value, root, parent, rule, rules, rulesets, pseudo, points, declarations) {
	var index = 0
	var offset = 0
	var length = pseudo
	var atrule = 0
	var property = 0
	var previous = 0
	var variable = 1
	var scanning = 1
	var ampersand = 1
	var character = 0
	var type = ''
	var props = rules
	var children = rulesets
	var reference = rule
	var characters = type

	while (scanning)
		switch (previous = character, character = next()) {
			// (
			case 40:
				if (previous != 108 && charat(characters, length - 1) == 58) {
					if (indexof(characters += replace(delimit(character), '&', '&\f'), '&\f', abs(index ? points[index - 1] : 0)) != -1)
						ampersand = -1
					break
				}
			// " ' [
			case 34: case 39: case 91:
				characters += delimit(character)
				break
			// \t \n \r \s
			case 9: case 10: case 13: case 32:
				characters += whitespace(previous)
				break
			// \
			case 92:
				characters += escaping(caret() - 1, 7)
				continue
			// /
			case 47:
				switch (peek()) {
					case 42: case 47:
						Utility_append(comment(commenter(next(), caret()), root, parent, declarations), declarations)
						if ((token(previous || 1) == 5 || token(peek() || 1) == 5) && strlen(characters) && substr(characters, -1, void 0) !== ' ') characters += ' '
						break
					default:
						characters += '/'
				}
				break
			// {
			case 123 * variable:
				points[index++] = strlen(characters) * ampersand
			// } ; \0
			case 125 * variable: case 59: case 0:
				switch (character) {
					// \0 }
					case 0: case 125: scanning = 0
					// ;
					case 59 + offset: if (ampersand == -1) characters = replace(characters, /\f/g, '')
						if (property > 0 && (strlen(characters) - length || (variable === 0 && previous === 47)))
							Utility_append(property > 32 ? declaration(characters + ';', rule, parent, length - 1, declarations) : declaration(replace(characters, ' ', '') + ';', rule, parent, length - 2, declarations), declarations)
						break
					// @ ;
					case 59: characters += ';'
					// { rule/at-rule
					default:
						Utility_append(reference = ruleset(characters, root, parent, index, offset, rules, points, type, props = [], children = [], length, rulesets), rulesets)

						if (character === 123)
							if (offset === 0)
								parse(characters, root, reference, reference, props, rulesets, length, points, children)
							else {
								switch (atrule) {
									// c(ontainer)
									case 99:
										if (charat(characters, 3) === 110) break
									// l(ayer)
									case 108:
										if (charat(characters, 2) === 97) break
									default:
										offset = 0
									// d(ocument) m(edia) s(upports)
									case 100: case 109: case 115:
								}
								if (offset) parse(value, reference, reference, rule && Utility_append(ruleset(value, reference, reference, 0, 0, rules, points, type, rules, props = [], length, children), children), rules, children, length, points, rule ? props : children)
								else parse(characters, reference, reference, reference, [''], children, 0, points, children)
							}
				}

				index = offset = property = 0, variable = ampersand = 1, type = characters = '', length = pseudo
				break
			// :
			case 58:
				length = 1 + strlen(characters), property = previous
			default:
				if (variable < 1)
					if (character == 123)
						--variable
					else if (character == 125 && variable++ == 0 && prev() == 125)
						continue

				switch (characters += Utility_from(character), character * variable) {
					// &
					case 38:
						ampersand = offset > 0 ? 1 : (characters += '\f', -1)
						break
					// ,
					case 44:
						points[index++] = (strlen(characters) - 1) * ampersand, ampersand = 1
						break
					// @
					case 64:
						// -
						if (peek() === 45)
							characters += delimit(next())

						atrule = peek(), offset = length = strlen(type = characters += identifier(caret())), character++
						break
					// -
					case 45:
						if (previous === 45 && strlen(characters) == 2)
							variable = 0
				}
		}

	return rulesets
}

/**
 * @param {string} value
 * @param {object} root
 * @param {object?} parent
 * @param {number} index
 * @param {number} offset
 * @param {string[]} rules
 * @param {number[]} points
 * @param {string} type
 * @param {string[]} props
 * @param {string[]} children
 * @param {number} length
 * @param {object[]} siblings
 * @return {object}
 */
function ruleset (value, root, parent, index, offset, rules, points, type, props, children, length, siblings) {
	var post = offset - 1
	var rule = offset === 0 ? rules : ['']
	var size = sizeof(rule)

	for (var i = 0, j = 0, k = 0; i < index; ++i)
		for (var x = 0, y = substr(value, post + 1, post = abs(j = points[i])), z = value; x < size; ++x)
			if (z = trim(j > 0 ? rule[x] + ' ' + y : replace(y, /&\f/g, rule[x])))
				props[k++] = z

	return node(value, root, parent, offset === 0 ? RULESET : type, props, children, length, siblings)
}

/**
 * @param {number} value
 * @param {object} root
 * @param {object?} parent
 * @param {object[]} siblings
 * @return {object}
 */
function comment (value, root, parent, siblings) {
	return node(value, root, parent, COMMENT, Utility_from(Tokenizer_char()), substr(value, 2, -2), 0, siblings)
}

/**
 * @param {string} value
 * @param {object} root
 * @param {object?} parent
 * @param {number} length
 * @param {object[]} siblings
 * @return {object}
 */
function declaration (value, root, parent, length, siblings) {
	return node(value, root, parent, DECLARATION, substr(value, 0, length), substr(value, length + 1, -1), length, siblings)
}

// EXTERNAL MODULE: ../node_modules/dompurify/dist/purify.es.mjs
var purify_es = __webpack_require__(4719);
// EXTERNAL MODULE: ../node_modules/lodash-es/isEmpty.js
var isEmpty = __webpack_require__(66400);
;// CONCATENATED MODULE: ../node_modules/mermaid/dist/mermaid.core.mjs














// src/mermaid.ts


// src/diagrams/c4/c4Detector.ts
var id = "c4";
var detector = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*C4Context|C4Container|C4Component|C4Dynamic|C4Deployment/.test(txt);
}, "detector");
var loader = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 2140).then(__webpack_require__.bind(__webpack_require__, 52140));
  return { id, diagram: diagram2 };
}, "loader");
var mermaid_core_plugin = {
  id,
  detector,
  loader
};
var c4Detector_default = mermaid_core_plugin;

// src/diagrams/flowchart/flowDetector.ts
var id2 = "flowchart";
var detector2 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt, config) => {
  if (config?.flowchart?.defaultRenderer === "dagre-wrapper" || config?.flowchart?.defaultRenderer === "elk") {
    return false;
  }
  return /^\s*graph/.test(txt);
}, "detector");
var loader2 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 652).then(__webpack_require__.bind(__webpack_require__, 652));
  return { id: id2, diagram: diagram2 };
}, "loader");
var plugin2 = {
  id: id2,
  detector: detector2,
  loader: loader2
};
var flowDetector_default = plugin2;

// src/diagrams/flowchart/flowDetector-v2.ts
var id3 = "flowchart-v2";
var detector3 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt, config) => {
  if (config?.flowchart?.defaultRenderer === "dagre-d3") {
    return false;
  }
  if (config?.flowchart?.defaultRenderer === "elk") {
    config.layout = "elk";
  }
  if (/^\s*graph/.test(txt) && config?.flowchart?.defaultRenderer === "dagre-wrapper") {
    return true;
  }
  return /^\s*flowchart/.test(txt);
}, "detector");
var loader3 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 652).then(__webpack_require__.bind(__webpack_require__, 652));
  return { id: id3, diagram: diagram2 };
}, "loader");
var plugin3 = {
  id: id3,
  detector: detector3,
  loader: loader3
};
var flowDetector_v2_default = plugin3;

// src/diagrams/er/erDetector.ts
var id4 = "er";
var detector4 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*erDiagram/.test(txt);
}, "detector");
var loader4 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 7988).then(__webpack_require__.bind(__webpack_require__, 57988));
  return { id: id4, diagram: diagram2 };
}, "loader");
var plugin4 = {
  id: id4,
  detector: detector4,
  loader: loader4
};
var erDetector_default = plugin4;

// src/diagrams/git/gitGraphDetector.ts
var id5 = "gitGraph";
var detector5 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*gitGraph/.test(txt);
}, "detector");
var loader5 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await Promise.all(/* import() */[__webpack_require__.e(2228), __webpack_require__.e(3197), __webpack_require__.e(3591)]).then(__webpack_require__.bind(__webpack_require__, 33591));
  return { id: id5, diagram: diagram2 };
}, "loader");
var plugin5 = {
  id: id5,
  detector: detector5,
  loader: loader5
};
var gitGraphDetector_default = plugin5;

// src/diagrams/gantt/ganttDetector.ts
var id6 = "gantt";
var detector6 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*gantt/.test(txt);
}, "detector");
var loader6 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 1830).then(__webpack_require__.bind(__webpack_require__, 11830));
  return { id: id6, diagram: diagram2 };
}, "loader");
var plugin6 = {
  id: id6,
  detector: detector6,
  loader: loader6
};
var ganttDetector_default = plugin6;

// src/diagrams/info/infoDetector.ts
var id7 = "info";
var detector7 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*info/.test(txt);
}, "detector");
var loader7 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await Promise.all(/* import() */[__webpack_require__.e(2228), __webpack_require__.e(3197), __webpack_require__.e(9373)]).then(__webpack_require__.bind(__webpack_require__, 99373));
  return { id: id7, diagram: diagram2 };
}, "loader");
var info = {
  id: id7,
  detector: detector7,
  loader: loader7
};

// src/diagrams/pie/pieDetector.ts
var id8 = "pie";
var detector8 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*pie/.test(txt);
}, "detector");
var loader8 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await Promise.all(/* import() */[__webpack_require__.e(2228), __webpack_require__.e(3197), __webpack_require__.e(1616)]).then(__webpack_require__.bind(__webpack_require__, 1616));
  return { id: id8, diagram: diagram2 };
}, "loader");
var pie = {
  id: id8,
  detector: detector8,
  loader: loader8
};

// src/diagrams/quadrant-chart/quadrantDetector.ts
var id9 = "quadrantChart";
var detector9 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*quadrantChart/.test(txt);
}, "detector");
var loader9 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 1679).then(__webpack_require__.bind(__webpack_require__, 61679));
  return { id: id9, diagram: diagram2 };
}, "loader");
var plugin7 = {
  id: id9,
  detector: detector9,
  loader: loader9
};
var quadrantDetector_default = plugin7;

// src/diagrams/xychart/xychartDetector.ts
var id10 = "xychart";
var detector10 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*xychart(-beta)?/.test(txt);
}, "detector");
var loader10 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 5468).then(__webpack_require__.bind(__webpack_require__, 95468));
  return { id: id10, diagram: diagram2 };
}, "loader");
var plugin8 = {
  id: id10,
  detector: detector10,
  loader: loader10
};
var xychartDetector_default = plugin8;

// src/diagrams/requirement/requirementDetector.ts
var id11 = "requirement";
var detector11 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*requirement(Diagram)?/.test(txt);
}, "detector");
var loader11 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 9848).then(__webpack_require__.bind(__webpack_require__, 89848));
  return { id: id11, diagram: diagram2 };
}, "loader");
var plugin9 = {
  id: id11,
  detector: detector11,
  loader: loader11
};
var requirementDetector_default = plugin9;

// src/diagrams/sequence/sequenceDetector.ts
var id12 = "sequence";
var detector12 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*sequenceDiagram/.test(txt);
}, "detector");
var loader12 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 9966).then(__webpack_require__.bind(__webpack_require__, 39966));
  return { id: id12, diagram: diagram2 };
}, "loader");
var plugin10 = {
  id: id12,
  detector: detector12,
  loader: loader12
};
var sequenceDetector_default = plugin10;

// src/diagrams/class/classDetector.ts
var id13 = "class";
var detector13 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt, config) => {
  if (config?.class?.defaultRenderer === "dagre-wrapper") {
    return false;
  }
  return /^\s*classDiagram/.test(txt);
}, "detector");
var loader13 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await Promise.all(/* import() */[__webpack_require__.e(4682), __webpack_require__.e(4058)]).then(__webpack_require__.bind(__webpack_require__, 36364));
  return { id: id13, diagram: diagram2 };
}, "loader");
var plugin11 = {
  id: id13,
  detector: detector13,
  loader: loader13
};
var classDetector_default = plugin11;

// src/diagrams/class/classDetector-V2.ts
var id14 = "classDiagram";
var detector14 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt, config) => {
  if (/^\s*classDiagram/.test(txt) && config?.class?.defaultRenderer === "dagre-wrapper") {
    return true;
  }
  return /^\s*classDiagram-v2/.test(txt);
}, "detector");
var loader14 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await Promise.all(/* import() */[__webpack_require__.e(4682), __webpack_require__.e(4135)]).then(__webpack_require__.bind(__webpack_require__, 74135));
  return { id: id14, diagram: diagram2 };
}, "loader");
var plugin12 = {
  id: id14,
  detector: detector14,
  loader: loader14
};
var classDetector_V2_default = plugin12;

// src/diagrams/state/stateDetector.ts
var id15 = "state";
var detector15 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt, config) => {
  if (config?.state?.defaultRenderer === "dagre-wrapper") {
    return false;
  }
  return /^\s*stateDiagram/.test(txt);
}, "detector");
var loader15 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await Promise.all(/* import() */[__webpack_require__.e(2228), __webpack_require__.e(4276), __webpack_require__.e(5726), __webpack_require__.e(8140)]).then(__webpack_require__.bind(__webpack_require__, 98140));
  return { id: id15, diagram: diagram2 };
}, "loader");
var plugin13 = {
  id: id15,
  detector: detector15,
  loader: loader15
};
var stateDetector_default = plugin13;

// src/diagrams/state/stateDetector-V2.ts
var id16 = "stateDiagram";
var detector16 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt, config) => {
  if (/^\s*stateDiagram-v2/.test(txt)) {
    return true;
  }
  if (/^\s*stateDiagram/.test(txt) && config?.state?.defaultRenderer === "dagre-wrapper") {
    return true;
  }
  return false;
}, "detector");
var loader16 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await Promise.all(/* import() */[__webpack_require__.e(5726), __webpack_require__.e(4062)]).then(__webpack_require__.bind(__webpack_require__, 54062));
  return { id: id16, diagram: diagram2 };
}, "loader");
var plugin14 = {
  id: id16,
  detector: detector16,
  loader: loader16
};
var stateDetector_V2_default = plugin14;

// src/diagrams/user-journey/journeyDetector.ts
var id17 = "journey";
var detector17 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*journey/.test(txt);
}, "detector");
var loader17 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 1225).then(__webpack_require__.bind(__webpack_require__, 34885));
  return { id: id17, diagram: diagram2 };
}, "loader");
var plugin15 = {
  id: id17,
  detector: detector17,
  loader: loader17
};
var journeyDetector_default = plugin15;

// src/diagrams/error/errorRenderer.ts
var draw = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((_text, id28, version) => {
  chunk_AGHRB4JF/* log */.cM.debug("rendering svg for syntax error\n");
  const svg = (0,chunk_EXTU4WIE/* selectSvgElement */.P)(id28);
  const g = svg.append("g");
  svg.attr("viewBox", "0 0 2412 512");
  (0,chunk_ABZYJK2D/* configureSvgSize */.v2)(svg, 100, 512, true);
  g.append("path").attr("class", "error-icon").attr(
    "d",
    "m411.313,123.313c6.25-6.25 6.25-16.375 0-22.625s-16.375-6.25-22.625,0l-32,32-9.375,9.375-20.688-20.688c-12.484-12.5-32.766-12.5-45.25,0l-16,16c-1.261,1.261-2.304,2.648-3.31,4.051-21.739-8.561-45.324-13.426-70.065-13.426-105.867,0-192,86.133-192,192s86.133,192 192,192 192-86.133 192-192c0-24.741-4.864-48.327-13.426-70.065 1.402-1.007 2.79-2.049 4.051-3.31l16-16c12.5-12.492 12.5-32.758 0-45.25l-20.688-20.688 9.375-9.375 32.001-31.999zm-219.313,100.687c-52.938,0-96,43.063-96,96 0,8.836-7.164,16-16,16s-16-7.164-16-16c0-70.578 57.422-128 128-128 8.836,0 16,7.164 16,16s-7.164,16-16,16z"
  );
  g.append("path").attr("class", "error-icon").attr(
    "d",
    "m459.02,148.98c-6.25-6.25-16.375-6.25-22.625,0s-6.25,16.375 0,22.625l16,16c3.125,3.125 7.219,4.688 11.313,4.688 4.094,0 8.188-1.563 11.313-4.688 6.25-6.25 6.25-16.375 0-22.625l-16.001-16z"
  );
  g.append("path").attr("class", "error-icon").attr(
    "d",
    "m340.395,75.605c3.125,3.125 7.219,4.688 11.313,4.688 4.094,0 8.188-1.563 11.313-4.688 6.25-6.25 6.25-16.375 0-22.625l-16-16c-6.25-6.25-16.375-6.25-22.625,0s-6.25,16.375 0,22.625l15.999,16z"
  );
  g.append("path").attr("class", "error-icon").attr(
    "d",
    "m400,64c8.844,0 16-7.164 16-16v-32c0-8.836-7.156-16-16-16-8.844,0-16,7.164-16,16v32c0,8.836 7.156,16 16,16z"
  );
  g.append("path").attr("class", "error-icon").attr(
    "d",
    "m496,96.586h-32c-8.844,0-16,7.164-16,16 0,8.836 7.156,16 16,16h32c8.844,0 16-7.164 16-16 0-8.836-7.156-16-16-16z"
  );
  g.append("path").attr("class", "error-icon").attr(
    "d",
    "m436.98,75.605c3.125,3.125 7.219,4.688 11.313,4.688 4.094,0 8.188-1.563 11.313-4.688l32-32c6.25-6.25 6.25-16.375 0-22.625s-16.375-6.25-22.625,0l-32,32c-6.251,6.25-6.251,16.375-0.001,22.625z"
  );
  g.append("text").attr("class", "error-text").attr("x", 1440).attr("y", 250).attr("font-size", "150px").style("text-anchor", "middle").text("Syntax error in text");
  g.append("text").attr("class", "error-text").attr("x", 1250).attr("y", 400).attr("font-size", "100px").style("text-anchor", "middle").text(`mermaid version ${version}`);
}, "draw");
var renderer = { draw };
var errorRenderer_default = renderer;

// src/diagrams/error/errorDiagram.ts
var diagram = {
  db: {},
  renderer,
  parser: {
    parse: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
      return;
    }, "parse")
  }
};
var errorDiagram_default = diagram;

// src/diagrams/flowchart/elk/detector.ts
var id18 = "flowchart-elk";
var detector18 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt, config = {}) => {
  if (
    // If diagram explicitly states flowchart-elk
    /^\s*flowchart-elk/.test(txt) || // If a flowchart/graph diagram has their default renderer set to elk
    /^\s*(flowchart|graph)/.test(txt) && config?.flowchart?.defaultRenderer === "elk"
  ) {
    config.layout = "elk";
    return true;
  }
  return false;
}, "detector");
var loader18 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 652).then(__webpack_require__.bind(__webpack_require__, 652));
  return { id: id18, diagram: diagram2 };
}, "loader");
var plugin16 = {
  id: id18,
  detector: detector18,
  loader: loader18
};
var detector_default = plugin16;

// src/diagrams/timeline/detector.ts
var id19 = "timeline";
var detector19 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*timeline/.test(txt);
}, "detector");
var loader19 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 6657).then(__webpack_require__.bind(__webpack_require__, 66657));
  return { id: id19, diagram: diagram2 };
}, "loader");
var plugin17 = {
  id: id19,
  detector: detector19,
  loader: loader19
};
var detector_default2 = plugin17;

// src/diagrams/mindmap/detector.ts
var id20 = "mindmap";
var detector20 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*mindmap/.test(txt);
}, "detector");
var loader20 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 3522).then(__webpack_require__.bind(__webpack_require__, 3522));
  return { id: id20, diagram: diagram2 };
}, "loader");
var plugin18 = {
  id: id20,
  detector: detector20,
  loader: loader20
};
var detector_default3 = plugin18;

// src/diagrams/kanban/detector.ts
var id21 = "kanban";
var detector21 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*kanban/.test(txt);
}, "detector");
var loader21 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 369).then(__webpack_require__.bind(__webpack_require__, 30369));
  return { id: id21, diagram: diagram2 };
}, "loader");
var plugin19 = {
  id: id21,
  detector: detector21,
  loader: loader21
};
var detector_default4 = plugin19;

// src/diagrams/sankey/sankeyDetector.ts
var id22 = "sankey";
var detector22 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*sankey(-beta)?/.test(txt);
}, "detector");
var loader22 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await __webpack_require__.e(/* import() */ 9448).then(__webpack_require__.bind(__webpack_require__, 59448));
  return { id: id22, diagram: diagram2 };
}, "loader");
var plugin20 = {
  id: id22,
  detector: detector22,
  loader: loader22
};
var sankeyDetector_default = plugin20;

// src/diagrams/packet/detector.ts
var id23 = "packet";
var detector23 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*packet(-beta)?/.test(txt);
}, "detector");
var loader23 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await Promise.all(/* import() */[__webpack_require__.e(2228), __webpack_require__.e(3197), __webpack_require__.e(5530)]).then(__webpack_require__.bind(__webpack_require__, 25530));
  return { id: id23, diagram: diagram2 };
}, "loader");
var packet = {
  id: id23,
  detector: detector23,
  loader: loader23
};

// src/diagrams/radar/detector.ts
var id24 = "radar";
var detector24 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*radar-beta/.test(txt);
}, "detector");
var loader24 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await Promise.all(/* import() */[__webpack_require__.e(2228), __webpack_require__.e(3197), __webpack_require__.e(7776)]).then(__webpack_require__.bind(__webpack_require__, 57776));
  return { id: id24, diagram: diagram2 };
}, "loader");
var radar = {
  id: id24,
  detector: detector24,
  loader: loader24
};

// src/diagrams/block/blockDetector.ts
var id25 = "block";
var detector25 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*block(-beta)?/.test(txt);
}, "detector");
var loader25 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await Promise.all(/* import() */[__webpack_require__.e(2228), __webpack_require__.e(342)]).then(__webpack_require__.bind(__webpack_require__, 70342));
  return { id: id25, diagram: diagram2 };
}, "loader");
var plugin21 = {
  id: id25,
  detector: detector25,
  loader: loader25
};
var blockDetector_default = plugin21;

// src/diagrams/architecture/architectureDetector.ts
var id26 = "architecture";
var detector26 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*architecture/.test(txt);
}, "detector");
var loader26 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await Promise.all(/* import() */[__webpack_require__.e(2228), __webpack_require__.e(3197), __webpack_require__.e(3207), __webpack_require__.e(5634)]).then(__webpack_require__.bind(__webpack_require__, 75634));
  return { id: id26, diagram: diagram2 };
}, "loader");
var architecture = {
  id: id26,
  detector: detector26,
  loader: loader26
};
var architectureDetector_default = architecture;

// src/diagrams/treemap/detector.ts
var id27 = "treemap";
var detector27 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((txt) => {
  return /^\s*treemap/.test(txt);
}, "detector");
var loader27 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  const { diagram: diagram2 } = await Promise.all(/* import() */[__webpack_require__.e(2228), __webpack_require__.e(3197), __webpack_require__.e(5643)]).then(__webpack_require__.bind(__webpack_require__, 95643));
  return { id: id27, diagram: diagram2 };
}, "loader");
var treemap = {
  id: id27,
  detector: detector27,
  loader: loader27
};

// src/diagram-api/diagram-orchestration.ts
var hasLoadedDiagrams = false;
var addDiagrams = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
  if (hasLoadedDiagrams) {
    return;
  }
  hasLoadedDiagrams = true;
  (0,chunk_ABZYJK2D/* registerDiagram */.Cq)("error", errorDiagram_default, (text) => {
    return text.toLowerCase().trim() === "error";
  });
  (0,chunk_ABZYJK2D/* registerDiagram */.Cq)(
    "---",
    // --- diagram type may appear if YAML front-matter is not parsed correctly
    {
      db: {
        clear: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
        }, "clear")
      },
      styles: {},
      // should never be used
      renderer: {
        draw: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
        }, "draw")
      },
      parser: {
        parse: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
          throw new Error(
            "Diagrams beginning with --- are not valid. If you were trying to use a YAML front-matter, please ensure that you've correctly opened and closed the YAML front-matter with un-indented `---` blocks"
          );
        }, "parse")
      },
      init: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => null, "init")
      // no op
    },
    (text) => {
      return text.toLowerCase().trimStart().startsWith("---");
    }
  );
  if (true) {
    (0,chunk_ABZYJK2D/* registerLazyLoadedDiagrams */.KO)(detector_default, detector_default3, architectureDetector_default);
  }
  (0,chunk_ABZYJK2D/* registerLazyLoadedDiagrams */.KO)(
    c4Detector_default,
    detector_default4,
    classDetector_V2_default,
    classDetector_default,
    erDetector_default,
    ganttDetector_default,
    info,
    pie,
    requirementDetector_default,
    sequenceDetector_default,
    flowDetector_v2_default,
    flowDetector_default,
    detector_default2,
    gitGraphDetector_default,
    stateDetector_V2_default,
    stateDetector_default,
    journeyDetector_default,
    quadrantDetector_default,
    sankeyDetector_default,
    packet,
    xychartDetector_default,
    blockDetector_default,
    radar,
    treemap
  );
}, "addDiagrams");

// src/diagram-api/loadDiagram.ts
var loadRegisteredDiagrams = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  chunk_AGHRB4JF/* log */.cM.debug(`Loading registered diagrams`);
  const results = await Promise.allSettled(
    Object.entries(chunk_ABZYJK2D/* detectors */.Bf).map(async ([key, { detector: detector28, loader: loader28 }]) => {
      if (!loader28) {
        return;
      }
      try {
        (0,chunk_ABZYJK2D/* getDiagram */._7)(key);
      } catch {
        try {
          const { diagram: diagram2, id: id28 } = await loader28();
          (0,chunk_ABZYJK2D/* registerDiagram */.Cq)(id28, diagram2, detector28);
        } catch (err) {
          chunk_AGHRB4JF/* log */.cM.error(`Failed to load external diagram with key ${key}. Removing from detectors.`);
          delete chunk_ABZYJK2D/* detectors */.Bf[key];
          throw err;
        }
      }
    })
  );
  const failed = results.filter((result) => result.status === "rejected");
  if (failed.length > 0) {
    chunk_AGHRB4JF/* log */.cM.error(`Failed to load ${failed.length} external diagrams`);
    for (const res of failed) {
      chunk_AGHRB4JF/* log */.cM.error(res);
    }
    throw new Error(`Failed to load ${failed.length} external diagrams`);
  }
}, "loadRegisteredDiagrams");

// src/mermaidAPI.ts





// src/accessibility.ts
var SVG_ROLE = "graphics-document document";
function setA11yDiagramInfo(svg, diagramType) {
  svg.attr("role", SVG_ROLE);
  if (diagramType !== "") {
    svg.attr("aria-roledescription", diagramType);
  }
}
(0,chunk_AGHRB4JF/* __name */.eW)(setA11yDiagramInfo, "setA11yDiagramInfo");
function addSVGa11yTitleDescription(svg, a11yTitle, a11yDesc, baseId) {
  if (svg.insert === void 0) {
    return;
  }
  if (a11yDesc) {
    const descId = `chart-desc-${baseId}`;
    svg.attr("aria-describedby", descId);
    svg.insert("desc", ":first-child").attr("id", descId).text(a11yDesc);
  }
  if (a11yTitle) {
    const titleId = `chart-title-${baseId}`;
    svg.attr("aria-labelledby", titleId);
    svg.insert("title", ":first-child").attr("id", titleId).text(a11yTitle);
  }
}
(0,chunk_AGHRB4JF/* __name */.eW)(addSVGa11yTitleDescription, "addSVGa11yTitleDescription");

// src/Diagram.ts
var Diagram = class _Diagram {
  constructor(type, text, db, parser, renderer2) {
    this.type = type;
    this.text = text;
    this.db = db;
    this.parser = parser;
    this.renderer = renderer2;
  }
  static {
    (0,chunk_AGHRB4JF/* __name */.eW)(this, "Diagram");
  }
  static async fromText(text, metadata = {}) {
    const config = (0,chunk_ABZYJK2D/* getConfig */.iE)();
    const type = (0,chunk_ABZYJK2D/* detectType */.Vg)(text, config);
    text = (0,chunk_S3R3BYOJ/* encodeEntities */.Vy)(text) + "\n";
    try {
      (0,chunk_ABZYJK2D/* getDiagram */._7)(type);
    } catch {
      const loader28 = (0,chunk_ABZYJK2D/* getDiagramLoader */.cq)(type);
      if (!loader28) {
        throw new chunk_ABZYJK2D/* UnknownDiagramError */.cj(`Diagram ${type} not found.`);
      }
      const { id: id28, diagram: diagram2 } = await loader28();
      (0,chunk_ABZYJK2D/* registerDiagram */.Cq)(id28, diagram2);
    }
    const { db, parser, renderer: renderer2, init: init2 } = (0,chunk_ABZYJK2D/* getDiagram */._7)(type);
    if (parser.parser) {
      parser.parser.yy = db;
    }
    db.clear?.();
    init2?.(config);
    if (metadata.title) {
      db.setDiagramTitle?.(metadata.title);
    }
    await parser.parse(text);
    return new _Diagram(type, text, db, parser, renderer2);
  }
  async render(id28, version) {
    await this.renderer.draw(this.text, id28, version, this);
  }
  getParser() {
    return this.parser;
  }
  getType() {
    return this.type;
  }
};

// src/interactionDb.ts
var interactionFunctions = [];
var attachFunctions = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
  interactionFunctions.forEach((f) => {
    f();
  });
  interactionFunctions = [];
}, "attachFunctions");

// src/diagram-api/comments.ts
var cleanupComments = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((text) => {
  return text.replace(/^\s*%%(?!{)[^\n]+\n?/gm, "").trimStart();
}, "cleanupComments");

// src/diagram-api/frontmatter.ts
function extractFrontMatter(text) {
  const matches = text.match(chunk_ABZYJK2D/* frontMatterRegex */.M6);
  if (!matches) {
    return {
      text,
      metadata: {}
    };
  }
  let parsed = (0,chunk_MI3HLSF2/* load */.z)(matches[1], {
    // To support config, we need JSON schema.
    // https://www.yaml.org/spec/1.2/spec.html#id2803231
    schema: chunk_MI3HLSF2/* JSON_SCHEMA */.A
  }) ?? {};
  parsed = typeof parsed === "object" && !Array.isArray(parsed) ? parsed : {};
  const metadata = {};
  if (parsed.displayMode) {
    metadata.displayMode = parsed.displayMode.toString();
  }
  if (parsed.title) {
    metadata.title = parsed.title.toString();
  }
  if (parsed.config) {
    metadata.config = parsed.config;
  }
  return {
    text: text.slice(matches[0].length),
    metadata
  };
}
(0,chunk_AGHRB4JF/* __name */.eW)(extractFrontMatter, "extractFrontMatter");

// src/preprocess.ts
var cleanupText = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((code) => {
  return code.replace(/\r\n?/g, "\n").replace(
    /<(\w+)([^>]*)>/g,
    (match, tag, attributes) => "<" + tag + attributes.replace(/="([^"]*)"/g, "='$1'") + ">"
  );
}, "cleanupText");
var processFrontmatter = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((code) => {
  const { text, metadata } = extractFrontMatter(code);
  const { displayMode, title, config = {} } = metadata;
  if (displayMode) {
    if (!config.gantt) {
      config.gantt = {};
    }
    config.gantt.displayMode = displayMode;
  }
  return { title, config, text };
}, "processFrontmatter");
var processDirectives = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((code) => {
  const initDirective = chunk_S3R3BYOJ/* utils_default */.w8.detectInit(code) ?? {};
  const wrapDirectives = chunk_S3R3BYOJ/* utils_default */.w8.detectDirective(code, "wrap");
  if (Array.isArray(wrapDirectives)) {
    initDirective.wrap = wrapDirectives.some(({ type }) => type === "wrap");
  } else if (wrapDirectives?.type === "wrap") {
    initDirective.wrap = true;
  }
  return {
    text: (0,chunk_S3R3BYOJ/* removeDirectives */.tf)(code),
    directive: initDirective
  };
}, "processDirectives");
function preprocessDiagram(code) {
  const cleanedCode = cleanupText(code);
  const frontMatterResult = processFrontmatter(cleanedCode);
  const directiveResult = processDirectives(frontMatterResult.text);
  const config = (0,chunk_S3R3BYOJ/* cleanAndMerge */.Rb)(frontMatterResult.config, directiveResult.directive);
  code = cleanupComments(directiveResult.text);
  return {
    code,
    title: frontMatterResult.title,
    config
  };
}
(0,chunk_AGHRB4JF/* __name */.eW)(preprocessDiagram, "preprocessDiagram");

// src/utils/base64.ts
function toBase64(str) {
  const utf8Bytes = new TextEncoder().encode(str);
  const utf8Str = Array.from(utf8Bytes, (byte) => String.fromCodePoint(byte)).join("");
  return btoa(utf8Str);
}
(0,chunk_AGHRB4JF/* __name */.eW)(toBase64, "toBase64");

// src/mermaidAPI.ts
var MAX_TEXTLENGTH = 5e4;
var MAX_TEXTLENGTH_EXCEEDED_MSG = "graph TB;a[Maximum text size in diagram exceeded];style a fill:#faa";
var SECURITY_LVL_SANDBOX = "sandbox";
var SECURITY_LVL_LOOSE = "loose";
var XMLNS_SVG_STD = "http://www.w3.org/2000/svg";
var XMLNS_XLINK_STD = "http://www.w3.org/1999/xlink";
var XMLNS_XHTML_STD = "http://www.w3.org/1999/xhtml";
var IFRAME_WIDTH = "100%";
var IFRAME_HEIGHT = "100%";
var IFRAME_STYLES = "border:0;margin:0;";
var IFRAME_BODY_STYLE = "margin:0";
var IFRAME_SANDBOX_OPTS = "allow-top-navigation-by-user-activation allow-popups";
var IFRAME_NOT_SUPPORTED_MSG = 'The "iframe" tag is not supported by your browser.';
var DOMPURIFY_TAGS = ["foreignobject"];
var DOMPURIFY_ATTR = ["dominant-baseline"];
function processAndSetConfigs(text) {
  const processed = preprocessDiagram(text);
  (0,chunk_ABZYJK2D/* reset */.mc)();
  (0,chunk_ABZYJK2D/* addDirective */.XV)(processed.config ?? {});
  return processed;
}
(0,chunk_AGHRB4JF/* __name */.eW)(processAndSetConfigs, "processAndSetConfigs");
async function mermaid_core_parse(text, parseOptions) {
  addDiagrams();
  try {
    const { code, config } = processAndSetConfigs(text);
    const diagram2 = await getDiagramFromText(code);
    return { diagramType: diagram2.type, config };
  } catch (error) {
    if (parseOptions?.suppressErrors) {
      return false;
    }
    throw error;
  }
}
(0,chunk_AGHRB4JF/* __name */.eW)(mermaid_core_parse, "parse");
var cssImportantStyles = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((cssClass, element, cssClasses = []) => {
  return `
.${cssClass} ${element} { ${cssClasses.join(" !important; ")} !important; }`;
}, "cssImportantStyles");
var createCssStyles = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((config, classDefs = /* @__PURE__ */ new Map()) => {
  let cssStyles = "";
  if (config.themeCSS !== void 0) {
    cssStyles += `
${config.themeCSS}`;
  }
  if (config.fontFamily !== void 0) {
    cssStyles += `
:root { --mermaid-font-family: ${config.fontFamily}}`;
  }
  if (config.altFontFamily !== void 0) {
    cssStyles += `
:root { --mermaid-alt-font-family: ${config.altFontFamily}}`;
  }
  if (classDefs instanceof Map) {
    const htmlLabels = config.htmlLabels ?? config.flowchart?.htmlLabels;
    const cssHtmlElements = ["> *", "span"];
    const cssShapeElements = ["rect", "polygon", "ellipse", "circle", "path"];
    const cssElements = htmlLabels ? cssHtmlElements : cssShapeElements;
    classDefs.forEach((styleClassDef) => {
      if (!(0,isEmpty/* default */.Z)(styleClassDef.styles)) {
        cssElements.forEach((cssElement) => {
          cssStyles += cssImportantStyles(styleClassDef.id, cssElement, styleClassDef.styles);
        });
      }
      if (!(0,isEmpty/* default */.Z)(styleClassDef.textStyles)) {
        cssStyles += cssImportantStyles(
          styleClassDef.id,
          "tspan",
          (styleClassDef?.textStyles || []).map((s) => s.replace("color", "fill"))
        );
      }
    });
  }
  return cssStyles;
}, "createCssStyles");
var createUserStyles = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((config, graphType, classDefs, svgId) => {
  const userCSSstyles = createCssStyles(config, classDefs);
  const allStyles = (0,chunk_ABZYJK2D/* styles_default */.Ee)(graphType, userCSSstyles, config.themeVariables);
  return serialize(compile(`${svgId}{${allStyles}}`), stringify);
}, "createUserStyles");
var cleanUpSvgCode = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((svgCode = "", inSandboxMode, useArrowMarkerUrls) => {
  let cleanedUpSvg = svgCode;
  if (!useArrowMarkerUrls && !inSandboxMode) {
    cleanedUpSvg = cleanedUpSvg.replace(
      /marker-end="url\([\d+./:=?A-Za-z-]*?#/g,
      'marker-end="url(#'
    );
  }
  cleanedUpSvg = (0,chunk_S3R3BYOJ/* decodeEntities */.SH)(cleanedUpSvg);
  cleanedUpSvg = cleanedUpSvg.replace(/<br>/g, "<br/>");
  return cleanedUpSvg;
}, "cleanUpSvgCode");
var putIntoIFrame = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((svgCode = "", svgElement) => {
  const height = svgElement?.viewBox?.baseVal?.height ? svgElement.viewBox.baseVal.height + "px" : IFRAME_HEIGHT;
  const base64encodedSrc = toBase64(`<body style="${IFRAME_BODY_STYLE}">${svgCode}</body>`);
  return `<iframe style="width:${IFRAME_WIDTH};height:${height};${IFRAME_STYLES}" src="data:text/html;charset=UTF-8;base64,${base64encodedSrc}" sandbox="${IFRAME_SANDBOX_OPTS}">
  ${IFRAME_NOT_SUPPORTED_MSG}
</iframe>`;
}, "putIntoIFrame");
var appendDivSvgG = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((parentRoot, id28, enclosingDivId, divStyle, svgXlink) => {
  const enclosingDiv = parentRoot.append("div");
  enclosingDiv.attr("id", enclosingDivId);
  if (divStyle) {
    enclosingDiv.attr("style", divStyle);
  }
  const svgNode = enclosingDiv.append("svg").attr("id", id28).attr("width", "100%").attr("xmlns", XMLNS_SVG_STD);
  if (svgXlink) {
    svgNode.attr("xmlns:xlink", svgXlink);
  }
  svgNode.append("g");
  return parentRoot;
}, "appendDivSvgG");
function sandboxedIframe(parentNode, iFrameId) {
  return parentNode.append("iframe").attr("id", iFrameId).attr("style", "width: 100%; height: 100%;").attr("sandbox", "");
}
(0,chunk_AGHRB4JF/* __name */.eW)(sandboxedIframe, "sandboxedIframe");
var removeExistingElements = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((doc, id28, divId, iFrameId) => {
  doc.getElementById(id28)?.remove();
  doc.getElementById(divId)?.remove();
  doc.getElementById(iFrameId)?.remove();
}, "removeExistingElements");
var render = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async function(id28, text, svgContainingElement) {
  addDiagrams();
  const processed = processAndSetConfigs(text);
  text = processed.code;
  const config = (0,chunk_ABZYJK2D/* getConfig */.iE)();
  chunk_AGHRB4JF/* log */.cM.debug(config);
  if (text.length > (config?.maxTextSize ?? MAX_TEXTLENGTH)) {
    text = MAX_TEXTLENGTH_EXCEEDED_MSG;
  }
  const idSelector = "#" + id28;
  const iFrameID = "i" + id28;
  const iFrameID_selector = "#" + iFrameID;
  const enclosingDivID = "d" + id28;
  const enclosingDivID_selector = "#" + enclosingDivID;
  const removeTempElements = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
    const tmpElementSelector = isSandboxed ? iFrameID_selector : enclosingDivID_selector;
    const node = (0,src/* select */.Ys)(tmpElementSelector).node();
    if (node && "remove" in node) {
      node.remove();
    }
  }, "removeTempElements");
  let root = (0,src/* select */.Ys)("body");
  const isSandboxed = config.securityLevel === SECURITY_LVL_SANDBOX;
  const isLooseSecurityLevel = config.securityLevel === SECURITY_LVL_LOOSE;
  const fontFamily = config.fontFamily;
  if (svgContainingElement !== void 0) {
    if (svgContainingElement) {
      svgContainingElement.innerHTML = "";
    }
    if (isSandboxed) {
      const iframe = sandboxedIframe((0,src/* select */.Ys)(svgContainingElement), iFrameID);
      root = (0,src/* select */.Ys)(iframe.nodes()[0].contentDocument.body);
      root.node().style.margin = 0;
    } else {
      root = (0,src/* select */.Ys)(svgContainingElement);
    }
    appendDivSvgG(root, id28, enclosingDivID, `font-family: ${fontFamily}`, XMLNS_XLINK_STD);
  } else {
    removeExistingElements(document, id28, enclosingDivID, iFrameID);
    if (isSandboxed) {
      const iframe = sandboxedIframe((0,src/* select */.Ys)("body"), iFrameID);
      root = (0,src/* select */.Ys)(iframe.nodes()[0].contentDocument.body);
      root.node().style.margin = 0;
    } else {
      root = (0,src/* select */.Ys)("body");
    }
    appendDivSvgG(root, id28, enclosingDivID);
  }
  let diag;
  let parseEncounteredException;
  try {
    diag = await Diagram.fromText(text, { title: processed.title });
  } catch (error) {
    if (config.suppressErrorRendering) {
      removeTempElements();
      throw error;
    }
    diag = await Diagram.fromText("error");
    parseEncounteredException = error;
  }
  const element = root.select(enclosingDivID_selector).node();
  const diagramType = diag.type;
  const svg = element.firstChild;
  const firstChild = svg.firstChild;
  const diagramClassDefs = diag.renderer.getClasses?.(text, diag);
  const rules = createUserStyles(config, diagramType, diagramClassDefs, idSelector);
  const style1 = document.createElement("style");
  style1.innerHTML = rules;
  svg.insertBefore(style1, firstChild);
  try {
    await diag.renderer.draw(text, id28, chunk_6MN3ZHY7/* package_default */.X.version, diag);
  } catch (e) {
    if (config.suppressErrorRendering) {
      removeTempElements();
    } else {
      errorRenderer_default.draw(text, id28, chunk_6MN3ZHY7/* package_default */.X.version);
    }
    throw e;
  }
  const svgNode = root.select(`${enclosingDivID_selector} svg`);
  const a11yTitle = diag.db.getAccTitle?.();
  const a11yDescr = diag.db.getAccDescription?.();
  addA11yInfo(diagramType, svgNode, a11yTitle, a11yDescr);
  root.select(`[id="${id28}"]`).selectAll("foreignobject > *").attr("xmlns", XMLNS_XHTML_STD);
  let svgCode = root.select(enclosingDivID_selector).node().innerHTML;
  chunk_AGHRB4JF/* log */.cM.debug("config.arrowMarkerAbsolute", config.arrowMarkerAbsolute);
  svgCode = cleanUpSvgCode(svgCode, isSandboxed, (0,chunk_ABZYJK2D/* evaluate */.ku)(config.arrowMarkerAbsolute));
  if (isSandboxed) {
    const svgEl = root.select(enclosingDivID_selector + " svg").node();
    svgCode = putIntoIFrame(svgCode, svgEl);
  } else if (!isLooseSecurityLevel) {
    svgCode = purify_es/* default */.Z.sanitize(svgCode, {
      ADD_TAGS: DOMPURIFY_TAGS,
      ADD_ATTR: DOMPURIFY_ATTR,
      HTML_INTEGRATION_POINTS: { foreignobject: true }
    });
  }
  attachFunctions();
  if (parseEncounteredException) {
    throw parseEncounteredException;
  }
  removeTempElements();
  return {
    diagramType,
    svg: svgCode,
    bindFunctions: diag.db.bindFunctions
  };
}, "render");
function initialize(userOptions = {}) {
  const options = (0,chunk_ABZYJK2D/* assignWithDepth_default */.Yc)({}, userOptions);
  if (options?.fontFamily && !options.themeVariables?.fontFamily) {
    if (!options.themeVariables) {
      options.themeVariables = {};
    }
    options.themeVariables.fontFamily = options.fontFamily;
  }
  (0,chunk_ABZYJK2D/* saveConfigFromInitialize */.dY)(options);
  if (options?.theme && options.theme in chunk_ABZYJK2D/* themes_default */._j) {
    options.themeVariables = chunk_ABZYJK2D/* themes_default */._j[options.theme].getThemeVariables(
      options.themeVariables
    );
  } else if (options) {
    options.themeVariables = chunk_ABZYJK2D/* themes_default */._j.default.getThemeVariables(options.themeVariables);
  }
  const config = typeof options === "object" ? (0,chunk_ABZYJK2D/* setSiteConfig */.Yn)(options) : (0,chunk_ABZYJK2D/* getSiteConfig */.ZD)();
  (0,chunk_AGHRB4JF/* setLogLevel */.Ub)(config.logLevel);
  addDiagrams();
}
(0,chunk_AGHRB4JF/* __name */.eW)(initialize, "initialize");
var getDiagramFromText = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((text, metadata = {}) => {
  const { code } = preprocessDiagram(text);
  return Diagram.fromText(code, metadata);
}, "getDiagramFromText");
function addA11yInfo(diagramType, svgNode, a11yTitle, a11yDescr) {
  setA11yDiagramInfo(svgNode, diagramType);
  addSVGa11yTitleDescription(svgNode, a11yTitle, a11yDescr, svgNode.attr("id"));
}
(0,chunk_AGHRB4JF/* __name */.eW)(addA11yInfo, "addA11yInfo");
var mermaidAPI = Object.freeze({
  render,
  parse: mermaid_core_parse,
  getDiagramFromText,
  initialize,
  getConfig: chunk_ABZYJK2D/* getConfig */.iE,
  setConfig: chunk_ABZYJK2D/* setConfig */.v6,
  getSiteConfig: chunk_ABZYJK2D/* getSiteConfig */.ZD,
  updateSiteConfig: chunk_ABZYJK2D/* updateSiteConfig */.Tb,
  reset: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
    (0,chunk_ABZYJK2D/* reset */.mc)();
  }, "reset"),
  globalReset: /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
    (0,chunk_ABZYJK2D/* reset */.mc)(chunk_ABZYJK2D/* defaultConfig */.u_);
  }, "globalReset"),
  defaultConfig: chunk_ABZYJK2D/* defaultConfig */.u_
});
(0,chunk_AGHRB4JF/* setLogLevel */.Ub)((0,chunk_ABZYJK2D/* getConfig */.iE)().logLevel);
(0,chunk_ABZYJK2D/* reset */.mc)((0,chunk_ABZYJK2D/* getConfig */.iE)());

// src/mermaid.ts
var handleError = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((error, errors, parseError) => {
  chunk_AGHRB4JF/* log */.cM.warn(error);
  if ((0,chunk_S3R3BYOJ/* isDetailedError */.bZ)(error)) {
    if (parseError) {
      parseError(error.str, error.hash);
    }
    errors.push({ ...error, message: error.str, error });
  } else {
    if (parseError) {
      parseError(error);
    }
    if (error instanceof Error) {
      errors.push({
        str: error.message,
        message: error.message,
        hash: error.name,
        error
      });
    }
  }
}, "handleError");
var run = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async function(options = {
  querySelector: ".mermaid"
}) {
  try {
    await runThrowsErrors(options);
  } catch (e) {
    if ((0,chunk_S3R3BYOJ/* isDetailedError */.bZ)(e)) {
      chunk_AGHRB4JF/* log */.cM.error(e.str);
    }
    if (mermaid.parseError) {
      mermaid.parseError(e);
    }
    if (!options.suppressErrors) {
      chunk_AGHRB4JF/* log */.cM.error("Use the suppressErrors option to suppress these errors");
      throw e;
    }
  }
}, "run");
var runThrowsErrors = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async function({ postRenderCallback, querySelector, nodes } = {
  querySelector: ".mermaid"
}) {
  const conf = mermaidAPI.getConfig();
  chunk_AGHRB4JF/* log */.cM.debug(`${!postRenderCallback ? "No " : ""}Callback function found`);
  let nodesToProcess;
  if (nodes) {
    nodesToProcess = nodes;
  } else if (querySelector) {
    nodesToProcess = document.querySelectorAll(querySelector);
  } else {
    throw new Error("Nodes and querySelector are both undefined");
  }
  chunk_AGHRB4JF/* log */.cM.debug(`Found ${nodesToProcess.length} diagrams`);
  if (conf?.startOnLoad !== void 0) {
    chunk_AGHRB4JF/* log */.cM.debug("Start On Load: " + conf?.startOnLoad);
    mermaidAPI.updateSiteConfig({ startOnLoad: conf?.startOnLoad });
  }
  const idGenerator = new chunk_S3R3BYOJ/* utils_default */.w8.InitIDGenerator(conf.deterministicIds, conf.deterministicIDSeed);
  let txt;
  const errors = [];
  for (const element of Array.from(nodesToProcess)) {
    chunk_AGHRB4JF/* log */.cM.info("Rendering diagram: " + element.id);
    if (element.getAttribute("data-processed")) {
      continue;
    }
    element.setAttribute("data-processed", "true");
    const id28 = `mermaid-${idGenerator.next()}`;
    txt = element.innerHTML;
    txt = (0,esm/* dedent */.Z)(chunk_S3R3BYOJ/* utils_default */.w8.entityDecode(txt)).trim().replace(/<br\s*\/?>/gi, "<br/>");
    const init2 = chunk_S3R3BYOJ/* utils_default */.w8.detectInit(txt);
    if (init2) {
      chunk_AGHRB4JF/* log */.cM.debug("Detected early reinit: ", init2);
    }
    try {
      const { svg, bindFunctions } = await render2(id28, txt, element);
      element.innerHTML = svg;
      if (postRenderCallback) {
        await postRenderCallback(id28);
      }
      if (bindFunctions) {
        bindFunctions(element);
      }
    } catch (error) {
      handleError(error, errors, mermaid.parseError);
    }
  }
  if (errors.length > 0) {
    throw errors[0];
  }
}, "runThrowsErrors");
var initialize2 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(config) {
  mermaidAPI.initialize(config);
}, "initialize");
var init = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async function(config, nodes, callback) {
  chunk_AGHRB4JF/* log */.cM.warn("mermaid.init is deprecated. Please use run instead.");
  if (config) {
    initialize2(config);
  }
  const runOptions = { postRenderCallback: callback, querySelector: ".mermaid" };
  if (typeof nodes === "string") {
    runOptions.querySelector = nodes;
  } else if (nodes) {
    if (nodes instanceof HTMLElement) {
      runOptions.nodes = [nodes];
    } else {
      runOptions.nodes = nodes;
    }
  }
  await run(runOptions);
}, "init");
var registerExternalDiagrams = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async (diagrams, {
  lazyLoad = true
} = {}) => {
  addDiagrams();
  (0,chunk_ABZYJK2D/* registerLazyLoadedDiagrams */.KO)(...diagrams);
  if (lazyLoad === false) {
    await loadRegisteredDiagrams();
  }
}, "registerExternalDiagrams");
var contentLoaded = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function() {
  if (mermaid.startOnLoad) {
    const { startOnLoad } = mermaidAPI.getConfig();
    if (startOnLoad) {
      mermaid.run().catch((err) => chunk_AGHRB4JF/* log */.cM.error("Mermaid failed to initialize", err));
    }
  }
}, "contentLoaded");
if (typeof document !== "undefined") {
  window.addEventListener("load", contentLoaded, false);
}
var setParseErrorHandler = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(function(parseErrorHandler) {
  mermaid.parseError = parseErrorHandler;
}, "setParseErrorHandler");
var executionQueue = [];
var executionQueueRunning = false;
var executeQueue = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async () => {
  if (executionQueueRunning) {
    return;
  }
  executionQueueRunning = true;
  while (executionQueue.length > 0) {
    const f = executionQueue.shift();
    if (f) {
      try {
        await f();
      } catch (e) {
        chunk_AGHRB4JF/* log */.cM.error("Error executing queue", e);
      }
    }
  }
  executionQueueRunning = false;
}, "executeQueue");
var parse2 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(async (text, parseOptions) => {
  return new Promise((resolve, reject) => {
    const performCall = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => new Promise((res, rej) => {
      mermaidAPI.parse(text, parseOptions).then(
        (r) => {
          res(r);
          resolve(r);
        },
        (e) => {
          chunk_AGHRB4JF/* log */.cM.error("Error parsing", e);
          mermaid.parseError?.(e);
          rej(e);
          reject(e);
        }
      );
    }), "performCall");
    executionQueue.push(performCall);
    executeQueue().catch(reject);
  });
}, "parse");
var render2 = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)((id28, text, container) => {
  return new Promise((resolve, reject) => {
    const performCall = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => new Promise((res, rej) => {
      mermaidAPI.render(id28, text, container).then(
        (r) => {
          res(r);
          resolve(r);
        },
        (e) => {
          chunk_AGHRB4JF/* log */.cM.error("Error parsing", e);
          mermaid.parseError?.(e);
          rej(e);
          reject(e);
        }
      );
    }), "performCall");
    executionQueue.push(performCall);
    executeQueue().catch(reject);
  });
}, "render");
var getRegisteredDiagramsMetadata = /* @__PURE__ */ (0,chunk_AGHRB4JF/* __name */.eW)(() => {
  return Object.keys(chunk_ABZYJK2D/* detectors */.Bf).map((id28) => ({
    id: id28
  }));
}, "getRegisteredDiagramsMetadata");
var mermaid = {
  startOnLoad: true,
  mermaidAPI,
  parse: parse2,
  render: render2,
  init,
  run,
  registerExternalDiagrams,
  registerLayoutLoaders: chunk_N4CR4FBY/* registerLayoutLoaders */.jM,
  initialize: initialize2,
  parseError: void 0,
  contentLoaded,
  setParseErrorHandler,
  detectType: chunk_ABZYJK2D/* detectType */.Vg,
  registerIconPacks: chunk_JA3XYJ7Z/* registerIconPacks */.ef,
  getRegisteredDiagramsMetadata
};
var mermaid_default = mermaid;

/*! Check if previously processed */
/*!
 * Wait for document loaded before starting the execution
 */


/***/ }),

/***/ 59395:
/***/ (function(module) {

!function(t,e){ true?module.exports=e():0}(this,(function(){"use strict";var t=1e3,e=6e4,n=36e5,r="millisecond",i="second",s="minute",u="hour",a="day",o="week",c="month",f="quarter",h="year",d="date",l="Invalid Date",$=/^(\d{4})[-/]?(\d{1,2})?[-/]?(\d{0,2})[Tt\s]*(\d{1,2})?:?(\d{1,2})?:?(\d{1,2})?[.:]?(\d+)?$/,y=/\[([^\]]+)]|Y{1,4}|M{1,4}|D{1,2}|d{1,4}|H{1,2}|h{1,2}|a|A|m{1,2}|s{1,2}|Z{1,2}|SSS/g,M={name:"en",weekdays:"Sunday_Monday_Tuesday_Wednesday_Thursday_Friday_Saturday".split("_"),months:"January_February_March_April_May_June_July_August_September_October_November_December".split("_"),ordinal:function(t){var e=["th","st","nd","rd"],n=t%100;return"["+t+(e[(n-20)%10]||e[n]||e[0])+"]"}},m=function(t,e,n){var r=String(t);return!r||r.length>=e?t:""+Array(e+1-r.length).join(n)+t},v={s:m,z:function(t){var e=-t.utcOffset(),n=Math.abs(e),r=Math.floor(n/60),i=n%60;return(e<=0?"+":"-")+m(r,2,"0")+":"+m(i,2,"0")},m:function t(e,n){if(e.date()<n.date())return-t(n,e);var r=12*(n.year()-e.year())+(n.month()-e.month()),i=e.clone().add(r,c),s=n-i<0,u=e.clone().add(r+(s?-1:1),c);return+(-(r+(n-i)/(s?i-u:u-i))||0)},a:function(t){return t<0?Math.ceil(t)||0:Math.floor(t)},p:function(t){return{M:c,y:h,w:o,d:a,D:d,h:u,m:s,s:i,ms:r,Q:f}[t]||String(t||"").toLowerCase().replace(/s$/,"")},u:function(t){return void 0===t}},g="en",D={};D[g]=M;var p="$isDayjsObject",S=function(t){return t instanceof _||!(!t||!t[p])},w=function t(e,n,r){var i;if(!e)return g;if("string"==typeof e){var s=e.toLowerCase();D[s]&&(i=s),n&&(D[s]=n,i=s);var u=e.split("-");if(!i&&u.length>1)return t(u[0])}else{var a=e.name;D[a]=e,i=a}return!r&&i&&(g=i),i||!r&&g},O=function(t,e){if(S(t))return t.clone();var n="object"==typeof e?e:{};return n.date=t,n.args=arguments,new _(n)},b=v;b.l=w,b.i=S,b.w=function(t,e){return O(t,{locale:e.$L,utc:e.$u,x:e.$x,$offset:e.$offset})};var _=function(){function M(t){this.$L=w(t.locale,null,!0),this.parse(t),this.$x=this.$x||t.x||{},this[p]=!0}var m=M.prototype;return m.parse=function(t){this.$d=function(t){var e=t.date,n=t.utc;if(null===e)return new Date(NaN);if(b.u(e))return new Date;if(e instanceof Date)return new Date(e);if("string"==typeof e&&!/Z$/i.test(e)){var r=e.match($);if(r){var i=r[2]-1||0,s=(r[7]||"0").substring(0,3);return n?new Date(Date.UTC(r[1],i,r[3]||1,r[4]||0,r[5]||0,r[6]||0,s)):new Date(r[1],i,r[3]||1,r[4]||0,r[5]||0,r[6]||0,s)}}return new Date(e)}(t),this.init()},m.init=function(){var t=this.$d;this.$y=t.getFullYear(),this.$M=t.getMonth(),this.$D=t.getDate(),this.$W=t.getDay(),this.$H=t.getHours(),this.$m=t.getMinutes(),this.$s=t.getSeconds(),this.$ms=t.getMilliseconds()},m.$utils=function(){return b},m.isValid=function(){return!(this.$d.toString()===l)},m.isSame=function(t,e){var n=O(t);return this.startOf(e)<=n&&n<=this.endOf(e)},m.isAfter=function(t,e){return O(t)<this.startOf(e)},m.isBefore=function(t,e){return this.endOf(e)<O(t)},m.$g=function(t,e,n){return b.u(t)?this[e]:this.set(n,t)},m.unix=function(){return Math.floor(this.valueOf()/1e3)},m.valueOf=function(){return this.$d.getTime()},m.startOf=function(t,e){var n=this,r=!!b.u(e)||e,f=b.p(t),l=function(t,e){var i=b.w(n.$u?Date.UTC(n.$y,e,t):new Date(n.$y,e,t),n);return r?i:i.endOf(a)},$=function(t,e){return b.w(n.toDate()[t].apply(n.toDate("s"),(r?[0,0,0,0]:[23,59,59,999]).slice(e)),n)},y=this.$W,M=this.$M,m=this.$D,v="set"+(this.$u?"UTC":"");switch(f){case h:return r?l(1,0):l(31,11);case c:return r?l(1,M):l(0,M+1);case o:var g=this.$locale().weekStart||0,D=(y<g?y+7:y)-g;return l(r?m-D:m+(6-D),M);case a:case d:return $(v+"Hours",0);case u:return $(v+"Minutes",1);case s:return $(v+"Seconds",2);case i:return $(v+"Milliseconds",3);default:return this.clone()}},m.endOf=function(t){return this.startOf(t,!1)},m.$set=function(t,e){var n,o=b.p(t),f="set"+(this.$u?"UTC":""),l=(n={},n[a]=f+"Date",n[d]=f+"Date",n[c]=f+"Month",n[h]=f+"FullYear",n[u]=f+"Hours",n[s]=f+"Minutes",n[i]=f+"Seconds",n[r]=f+"Milliseconds",n)[o],$=o===a?this.$D+(e-this.$W):e;if(o===c||o===h){var y=this.clone().set(d,1);y.$d[l]($),y.init(),this.$d=y.set(d,Math.min(this.$D,y.daysInMonth())).$d}else l&&this.$d[l]($);return this.init(),this},m.set=function(t,e){return this.clone().$set(t,e)},m.get=function(t){return this[b.p(t)]()},m.add=function(r,f){var d,l=this;r=Number(r);var $=b.p(f),y=function(t){var e=O(l);return b.w(e.date(e.date()+Math.round(t*r)),l)};if($===c)return this.set(c,this.$M+r);if($===h)return this.set(h,this.$y+r);if($===a)return y(1);if($===o)return y(7);var M=(d={},d[s]=e,d[u]=n,d[i]=t,d)[$]||1,m=this.$d.getTime()+r*M;return b.w(m,this)},m.subtract=function(t,e){return this.add(-1*t,e)},m.format=function(t){var e=this,n=this.$locale();if(!this.isValid())return n.invalidDate||l;var r=t||"YYYY-MM-DDTHH:mm:ssZ",i=b.z(this),s=this.$H,u=this.$m,a=this.$M,o=n.weekdays,c=n.months,f=n.meridiem,h=function(t,n,i,s){return t&&(t[n]||t(e,r))||i[n].slice(0,s)},d=function(t){return b.s(s%12||12,t,"0")},$=f||function(t,e,n){var r=t<12?"AM":"PM";return n?r.toLowerCase():r};return r.replace(y,(function(t,r){return r||function(t){switch(t){case"YY":return String(e.$y).slice(-2);case"YYYY":return b.s(e.$y,4,"0");case"M":return a+1;case"MM":return b.s(a+1,2,"0");case"MMM":return h(n.monthsShort,a,c,3);case"MMMM":return h(c,a);case"D":return e.$D;case"DD":return b.s(e.$D,2,"0");case"d":return String(e.$W);case"dd":return h(n.weekdaysMin,e.$W,o,2);case"ddd":return h(n.weekdaysShort,e.$W,o,3);case"dddd":return o[e.$W];case"H":return String(s);case"HH":return b.s(s,2,"0");case"h":return d(1);case"hh":return d(2);case"a":return $(s,u,!0);case"A":return $(s,u,!1);case"m":return String(u);case"mm":return b.s(u,2,"0");case"s":return String(e.$s);case"ss":return b.s(e.$s,2,"0");case"SSS":return b.s(e.$ms,3,"0");case"Z":return i}return null}(t)||i.replace(":","")}))},m.utcOffset=function(){return 15*-Math.round(this.$d.getTimezoneOffset()/15)},m.diff=function(r,d,l){var $,y=this,M=b.p(d),m=O(r),v=(m.utcOffset()-this.utcOffset())*e,g=this-m,D=function(){return b.m(y,m)};switch(M){case h:$=D()/12;break;case c:$=D();break;case f:$=D()/3;break;case o:$=(g-v)/6048e5;break;case a:$=(g-v)/864e5;break;case u:$=g/n;break;case s:$=g/e;break;case i:$=g/t;break;default:$=g}return l?$:b.a($)},m.daysInMonth=function(){return this.endOf(c).$D},m.$locale=function(){return D[this.$L]},m.locale=function(t,e){if(!t)return this.$L;var n=this.clone(),r=w(t,e,!0);return r&&(n.$L=r),n},m.clone=function(){return b.w(this.$d,this)},m.toDate=function(){return new Date(this.valueOf())},m.toJSON=function(){return this.isValid()?this.toISOString():null},m.toISOString=function(){return this.$d.toISOString()},m.toString=function(){return this.$d.toUTCString()},M}(),k=_.prototype;return O.prototype=k,[["$ms",r],["$s",i],["$m",s],["$H",u],["$W",a],["$M",c],["$y",h],["$D",d]].forEach((function(t){k[t[1]]=function(e){return this.$g(e,t[0],t[1])}})),O.extend=function(t,e){return t.$i||(t(e,_,O),t.$i=!0),O},O.locale=w,O.isDayjs=S,O.unix=function(t){return O(1e3*t)},O.en=D[g],O.Ls=D,O.p={},O}));

/***/ }),

/***/ 36366:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (/* binding */ at)
/* harmony export */ });
function t(t,e,s){if(t&&t.length){const[n,o]=e,a=Math.PI/180*s,h=Math.cos(a),r=Math.sin(a);for(const e of t){const[t,s]=e;e[0]=(t-n)*h-(s-o)*r+n,e[1]=(t-n)*r+(s-o)*h+o}}}function e(t,e){return t[0]===e[0]&&t[1]===e[1]}function s(s,n,o,a=1){const h=o,r=Math.max(n,.1),i=s[0]&&s[0][0]&&"number"==typeof s[0][0]?[s]:s,c=[0,0];if(h)for(const e of i)t(e,c,h);const l=function(t,s,n){const o=[];for(const s of t){const t=[...s];e(t[0],t[t.length-1])||t.push([t[0][0],t[0][1]]),t.length>2&&o.push(t)}const a=[];s=Math.max(s,.1);const h=[];for(const t of o)for(let e=0;e<t.length-1;e++){const s=t[e],n=t[e+1];if(s[1]!==n[1]){const t=Math.min(s[1],n[1]);h.push({ymin:t,ymax:Math.max(s[1],n[1]),x:t===s[1]?s[0]:n[0],islope:(n[0]-s[0])/(n[1]-s[1])})}}if(h.sort(((t,e)=>t.ymin<e.ymin?-1:t.ymin>e.ymin?1:t.x<e.x?-1:t.x>e.x?1:t.ymax===e.ymax?0:(t.ymax-e.ymax)/Math.abs(t.ymax-e.ymax))),!h.length)return a;let r=[],i=h[0].ymin,c=0;for(;r.length||h.length;){if(h.length){let t=-1;for(let e=0;e<h.length&&!(h[e].ymin>i);e++)t=e;h.splice(0,t+1).forEach((t=>{r.push({s:i,edge:t})}))}if(r=r.filter((t=>!(t.edge.ymax<=i))),r.sort(((t,e)=>t.edge.x===e.edge.x?0:(t.edge.x-e.edge.x)/Math.abs(t.edge.x-e.edge.x))),(1!==n||c%s==0)&&r.length>1)for(let t=0;t<r.length;t+=2){const e=t+1;if(e>=r.length)break;const s=r[t].edge,n=r[e].edge;a.push([[Math.round(s.x),i],[Math.round(n.x),i]])}i+=n,r.forEach((t=>{t.edge.x=t.edge.x+n*t.edge.islope})),c++}return a}(i,r,a);if(h){for(const e of i)t(e,c,-h);!function(e,s,n){const o=[];e.forEach((t=>o.push(...t))),t(o,s,n)}(l,c,-h)}return l}function n(t,e){var n;const o=e.hachureAngle+90;let a=e.hachureGap;a<0&&(a=4*e.strokeWidth),a=Math.round(Math.max(a,.1));let h=1;return e.roughness>=1&&((null===(n=e.randomizer)||void 0===n?void 0:n.next())||Math.random())>.7&&(h=a),s(t,a,o,h||1)}class o{constructor(t){this.helper=t}fillPolygons(t,e){return this._fillPolygons(t,e)}_fillPolygons(t,e){const s=n(t,e);return{type:"fillSketch",ops:this.renderLines(s,e)}}renderLines(t,e){const s=[];for(const n of t)s.push(...this.helper.doubleLineOps(n[0][0],n[0][1],n[1][0],n[1][1],e));return s}}function a(t){const e=t[0],s=t[1];return Math.sqrt(Math.pow(e[0]-s[0],2)+Math.pow(e[1]-s[1],2))}class h extends o{fillPolygons(t,e){let s=e.hachureGap;s<0&&(s=4*e.strokeWidth),s=Math.max(s,.1);const o=n(t,Object.assign({},e,{hachureGap:s})),h=Math.PI/180*e.hachureAngle,r=[],i=.5*s*Math.cos(h),c=.5*s*Math.sin(h);for(const[t,e]of o)a([t,e])&&r.push([[t[0]-i,t[1]+c],[...e]],[[t[0]+i,t[1]-c],[...e]]);return{type:"fillSketch",ops:this.renderLines(r,e)}}}class r extends o{fillPolygons(t,e){const s=this._fillPolygons(t,e),n=Object.assign({},e,{hachureAngle:e.hachureAngle+90}),o=this._fillPolygons(t,n);return s.ops=s.ops.concat(o.ops),s}}class i{constructor(t){this.helper=t}fillPolygons(t,e){const s=n(t,e=Object.assign({},e,{hachureAngle:0}));return this.dotsOnLines(s,e)}dotsOnLines(t,e){const s=[];let n=e.hachureGap;n<0&&(n=4*e.strokeWidth),n=Math.max(n,.1);let o=e.fillWeight;o<0&&(o=e.strokeWidth/2);const h=n/4;for(const r of t){const t=a(r),i=t/n,c=Math.ceil(i)-1,l=t-c*n,u=(r[0][0]+r[1][0])/2-n/4,p=Math.min(r[0][1],r[1][1]);for(let t=0;t<c;t++){const a=p+l+t*n,r=u-h+2*Math.random()*h,i=a-h+2*Math.random()*h,c=this.helper.ellipse(r,i,o,o,e);s.push(...c.ops)}}return{type:"fillSketch",ops:s}}}class c{constructor(t){this.helper=t}fillPolygons(t,e){const s=n(t,e);return{type:"fillSketch",ops:this.dashedLine(s,e)}}dashedLine(t,e){const s=e.dashOffset<0?e.hachureGap<0?4*e.strokeWidth:e.hachureGap:e.dashOffset,n=e.dashGap<0?e.hachureGap<0?4*e.strokeWidth:e.hachureGap:e.dashGap,o=[];return t.forEach((t=>{const h=a(t),r=Math.floor(h/(s+n)),i=(h+n-r*(s+n))/2;let c=t[0],l=t[1];c[0]>l[0]&&(c=t[1],l=t[0]);const u=Math.atan((l[1]-c[1])/(l[0]-c[0]));for(let t=0;t<r;t++){const a=t*(s+n),h=a+s,r=[c[0]+a*Math.cos(u)+i*Math.cos(u),c[1]+a*Math.sin(u)+i*Math.sin(u)],l=[c[0]+h*Math.cos(u)+i*Math.cos(u),c[1]+h*Math.sin(u)+i*Math.sin(u)];o.push(...this.helper.doubleLineOps(r[0],r[1],l[0],l[1],e))}})),o}}class l{constructor(t){this.helper=t}fillPolygons(t,e){const s=e.hachureGap<0?4*e.strokeWidth:e.hachureGap,o=e.zigzagOffset<0?s:e.zigzagOffset,a=n(t,e=Object.assign({},e,{hachureGap:s+o}));return{type:"fillSketch",ops:this.zigzagLines(a,o,e)}}zigzagLines(t,e,s){const n=[];return t.forEach((t=>{const o=a(t),h=Math.round(o/(2*e));let r=t[0],i=t[1];r[0]>i[0]&&(r=t[1],i=t[0]);const c=Math.atan((i[1]-r[1])/(i[0]-r[0]));for(let t=0;t<h;t++){const o=2*t*e,a=2*(t+1)*e,h=Math.sqrt(2*Math.pow(e,2)),i=[r[0]+o*Math.cos(c),r[1]+o*Math.sin(c)],l=[r[0]+a*Math.cos(c),r[1]+a*Math.sin(c)],u=[i[0]+h*Math.cos(c+Math.PI/4),i[1]+h*Math.sin(c+Math.PI/4)];n.push(...this.helper.doubleLineOps(i[0],i[1],u[0],u[1],s),...this.helper.doubleLineOps(u[0],u[1],l[0],l[1],s))}})),n}}const u={};class p{constructor(t){this.seed=t}next(){return this.seed?(2**31-1&(this.seed=Math.imul(48271,this.seed)))/2**31:Math.random()}}const f=0,d=1,g=2,M={A:7,a:7,C:6,c:6,H:1,h:1,L:2,l:2,M:2,m:2,Q:4,q:4,S:4,s:4,T:2,t:2,V:1,v:1,Z:0,z:0};function k(t,e){return t.type===e}function b(t){const e=[],s=function(t){const e=new Array;for(;""!==t;)if(t.match(/^([ \t\r\n,]+)/))t=t.substr(RegExp.$1.length);else if(t.match(/^([aAcChHlLmMqQsStTvVzZ])/))e[e.length]={type:f,text:RegExp.$1},t=t.substr(RegExp.$1.length);else{if(!t.match(/^(([-+]?[0-9]+(\.[0-9]*)?|[-+]?\.[0-9]+)([eE][-+]?[0-9]+)?)/))return[];e[e.length]={type:d,text:`${parseFloat(RegExp.$1)}`},t=t.substr(RegExp.$1.length)}return e[e.length]={type:g,text:""},e}(t);let n="BOD",o=0,a=s[o];for(;!k(a,g);){let h=0;const r=[];if("BOD"===n){if("M"!==a.text&&"m"!==a.text)return b("M0,0"+t);o++,h=M[a.text],n=a.text}else k(a,d)?h=M[n]:(o++,h=M[a.text],n=a.text);if(!(o+h<s.length))throw new Error("Path data ended short");for(let t=o;t<o+h;t++){const e=s[t];if(!k(e,d))throw new Error("Param not a number: "+n+","+e.text);r[r.length]=+e.text}if("number"!=typeof M[n])throw new Error("Bad segment: "+n);{const t={key:n,data:r};e.push(t),o+=h,a=s[o],"M"===n&&(n="L"),"m"===n&&(n="l")}}return e}function y(t){let e=0,s=0,n=0,o=0;const a=[];for(const{key:h,data:r}of t)switch(h){case"M":a.push({key:"M",data:[...r]}),[e,s]=r,[n,o]=r;break;case"m":e+=r[0],s+=r[1],a.push({key:"M",data:[e,s]}),n=e,o=s;break;case"L":a.push({key:"L",data:[...r]}),[e,s]=r;break;case"l":e+=r[0],s+=r[1],a.push({key:"L",data:[e,s]});break;case"C":a.push({key:"C",data:[...r]}),e=r[4],s=r[5];break;case"c":{const t=r.map(((t,n)=>n%2?t+s:t+e));a.push({key:"C",data:t}),e=t[4],s=t[5];break}case"Q":a.push({key:"Q",data:[...r]}),e=r[2],s=r[3];break;case"q":{const t=r.map(((t,n)=>n%2?t+s:t+e));a.push({key:"Q",data:t}),e=t[2],s=t[3];break}case"A":a.push({key:"A",data:[...r]}),e=r[5],s=r[6];break;case"a":e+=r[5],s+=r[6],a.push({key:"A",data:[r[0],r[1],r[2],r[3],r[4],e,s]});break;case"H":a.push({key:"H",data:[...r]}),e=r[0];break;case"h":e+=r[0],a.push({key:"H",data:[e]});break;case"V":a.push({key:"V",data:[...r]}),s=r[0];break;case"v":s+=r[0],a.push({key:"V",data:[s]});break;case"S":a.push({key:"S",data:[...r]}),e=r[2],s=r[3];break;case"s":{const t=r.map(((t,n)=>n%2?t+s:t+e));a.push({key:"S",data:t}),e=t[2],s=t[3];break}case"T":a.push({key:"T",data:[...r]}),e=r[0],s=r[1];break;case"t":e+=r[0],s+=r[1],a.push({key:"T",data:[e,s]});break;case"Z":case"z":a.push({key:"Z",data:[]}),e=n,s=o}return a}function m(t){const e=[];let s="",n=0,o=0,a=0,h=0,r=0,i=0;for(const{key:c,data:l}of t){switch(c){case"M":e.push({key:"M",data:[...l]}),[n,o]=l,[a,h]=l;break;case"C":e.push({key:"C",data:[...l]}),n=l[4],o=l[5],r=l[2],i=l[3];break;case"L":e.push({key:"L",data:[...l]}),[n,o]=l;break;case"H":n=l[0],e.push({key:"L",data:[n,o]});break;case"V":o=l[0],e.push({key:"L",data:[n,o]});break;case"S":{let t=0,a=0;"C"===s||"S"===s?(t=n+(n-r),a=o+(o-i)):(t=n,a=o),e.push({key:"C",data:[t,a,...l]}),r=l[0],i=l[1],n=l[2],o=l[3];break}case"T":{const[t,a]=l;let h=0,c=0;"Q"===s||"T"===s?(h=n+(n-r),c=o+(o-i)):(h=n,c=o);const u=n+2*(h-n)/3,p=o+2*(c-o)/3,f=t+2*(h-t)/3,d=a+2*(c-a)/3;e.push({key:"C",data:[u,p,f,d,t,a]}),r=h,i=c,n=t,o=a;break}case"Q":{const[t,s,a,h]=l,c=n+2*(t-n)/3,u=o+2*(s-o)/3,p=a+2*(t-a)/3,f=h+2*(s-h)/3;e.push({key:"C",data:[c,u,p,f,a,h]}),r=t,i=s,n=a,o=h;break}case"A":{const t=Math.abs(l[0]),s=Math.abs(l[1]),a=l[2],h=l[3],r=l[4],i=l[5],c=l[6];if(0===t||0===s)e.push({key:"C",data:[n,o,i,c,i,c]}),n=i,o=c;else if(n!==i||o!==c){x(n,o,i,c,t,s,a,h,r).forEach((function(t){e.push({key:"C",data:t})})),n=i,o=c}break}case"Z":e.push({key:"Z",data:[]}),n=a,o=h}s=c}return e}function w(t,e,s){return[t*Math.cos(s)-e*Math.sin(s),t*Math.sin(s)+e*Math.cos(s)]}function x(t,e,s,n,o,a,h,r,i,c){const l=(u=h,Math.PI*u/180);var u;let p=[],f=0,d=0,g=0,M=0;if(c)[f,d,g,M]=c;else{[t,e]=w(t,e,-l),[s,n]=w(s,n,-l);const h=(t-s)/2,c=(e-n)/2;let u=h*h/(o*o)+c*c/(a*a);u>1&&(u=Math.sqrt(u),o*=u,a*=u);const p=o*o,k=a*a,b=p*k-p*c*c-k*h*h,y=p*c*c+k*h*h,m=(r===i?-1:1)*Math.sqrt(Math.abs(b/y));g=m*o*c/a+(t+s)/2,M=m*-a*h/o+(e+n)/2,f=Math.asin(parseFloat(((e-M)/a).toFixed(9))),d=Math.asin(parseFloat(((n-M)/a).toFixed(9))),t<g&&(f=Math.PI-f),s<g&&(d=Math.PI-d),f<0&&(f=2*Math.PI+f),d<0&&(d=2*Math.PI+d),i&&f>d&&(f-=2*Math.PI),!i&&d>f&&(d-=2*Math.PI)}let k=d-f;if(Math.abs(k)>120*Math.PI/180){const t=d,e=s,r=n;d=i&&d>f?f+120*Math.PI/180*1:f+120*Math.PI/180*-1,p=x(s=g+o*Math.cos(d),n=M+a*Math.sin(d),e,r,o,a,h,0,i,[d,t,g,M])}k=d-f;const b=Math.cos(f),y=Math.sin(f),m=Math.cos(d),P=Math.sin(d),v=Math.tan(k/4),S=4/3*o*v,O=4/3*a*v,L=[t,e],T=[t+S*y,e-O*b],D=[s+S*P,n-O*m],A=[s,n];if(T[0]=2*L[0]-T[0],T[1]=2*L[1]-T[1],c)return[T,D,A].concat(p);{p=[T,D,A].concat(p);const t=[];for(let e=0;e<p.length;e+=3){const s=w(p[e][0],p[e][1],l),n=w(p[e+1][0],p[e+1][1],l),o=w(p[e+2][0],p[e+2][1],l);t.push([s[0],s[1],n[0],n[1],o[0],o[1]])}return t}}const P={randOffset:function(t,e){return G(t,e)},randOffsetWithRange:function(t,e,s){return E(t,e,s)},ellipse:function(t,e,s,n,o){const a=T(s,n,o);return D(t,e,o,a).opset},doubleLineOps:function(t,e,s,n,o){return $(t,e,s,n,o,!0)}};function v(t,e,s,n,o){return{type:"path",ops:$(t,e,s,n,o)}}function S(t,e,s){const n=(t||[]).length;if(n>2){const o=[];for(let e=0;e<n-1;e++)o.push(...$(t[e][0],t[e][1],t[e+1][0],t[e+1][1],s));return e&&o.push(...$(t[n-1][0],t[n-1][1],t[0][0],t[0][1],s)),{type:"path",ops:o}}return 2===n?v(t[0][0],t[0][1],t[1][0],t[1][1],s):{type:"path",ops:[]}}function O(t,e,s,n,o){return function(t,e){return S(t,!0,e)}([[t,e],[t+s,e],[t+s,e+n],[t,e+n]],o)}function L(t,e){if(t.length){const s="number"==typeof t[0][0]?[t]:t,n=j(s[0],1*(1+.2*e.roughness),e),o=e.disableMultiStroke?[]:j(s[0],1.5*(1+.22*e.roughness),z(e));for(let t=1;t<s.length;t++){const a=s[t];if(a.length){const t=j(a,1*(1+.2*e.roughness),e),s=e.disableMultiStroke?[]:j(a,1.5*(1+.22*e.roughness),z(e));for(const e of t)"move"!==e.op&&n.push(e);for(const t of s)"move"!==t.op&&o.push(t)}}return{type:"path",ops:n.concat(o)}}return{type:"path",ops:[]}}function T(t,e,s){const n=Math.sqrt(2*Math.PI*Math.sqrt((Math.pow(t/2,2)+Math.pow(e/2,2))/2)),o=Math.ceil(Math.max(s.curveStepCount,s.curveStepCount/Math.sqrt(200)*n)),a=2*Math.PI/o;let h=Math.abs(t/2),r=Math.abs(e/2);const i=1-s.curveFitting;return h+=G(h*i,s),r+=G(r*i,s),{increment:a,rx:h,ry:r}}function D(t,e,s,n){const[o,a]=F(n.increment,t,e,n.rx,n.ry,1,n.increment*E(.1,E(.4,1,s),s),s);let h=q(o,null,s);if(!s.disableMultiStroke&&0!==s.roughness){const[o]=F(n.increment,t,e,n.rx,n.ry,1.5,0,s),a=q(o,null,s);h=h.concat(a)}return{estimatedPoints:a,opset:{type:"path",ops:h}}}function A(t,e,s,n,o,a,h,r,i){const c=t,l=e;let u=Math.abs(s/2),p=Math.abs(n/2);u+=G(.01*u,i),p+=G(.01*p,i);let f=o,d=a;for(;f<0;)f+=2*Math.PI,d+=2*Math.PI;d-f>2*Math.PI&&(f=0,d=2*Math.PI);const g=2*Math.PI/i.curveStepCount,M=Math.min(g/2,(d-f)/2),k=V(M,c,l,u,p,f,d,1,i);if(!i.disableMultiStroke){const t=V(M,c,l,u,p,f,d,1.5,i);k.push(...t)}return h&&(r?k.push(...$(c,l,c+u*Math.cos(f),l+p*Math.sin(f),i),...$(c,l,c+u*Math.cos(d),l+p*Math.sin(d),i)):k.push({op:"lineTo",data:[c,l]},{op:"lineTo",data:[c+u*Math.cos(f),l+p*Math.sin(f)]})),{type:"path",ops:k}}function _(t,e){const s=m(y(b(t))),n=[];let o=[0,0],a=[0,0];for(const{key:t,data:h}of s)switch(t){case"M":a=[h[0],h[1]],o=[h[0],h[1]];break;case"L":n.push(...$(a[0],a[1],h[0],h[1],e)),a=[h[0],h[1]];break;case"C":{const[t,s,o,r,i,c]=h;n.push(...Z(t,s,o,r,i,c,a,e)),a=[i,c];break}case"Z":n.push(...$(a[0],a[1],o[0],o[1],e)),a=[o[0],o[1]]}return{type:"path",ops:n}}function I(t,e){const s=[];for(const n of t)if(n.length){const t=e.maxRandomnessOffset||0,o=n.length;if(o>2){s.push({op:"move",data:[n[0][0]+G(t,e),n[0][1]+G(t,e)]});for(let a=1;a<o;a++)s.push({op:"lineTo",data:[n[a][0]+G(t,e),n[a][1]+G(t,e)]})}}return{type:"fillPath",ops:s}}function C(t,e){return function(t,e){let s=t.fillStyle||"hachure";if(!u[s])switch(s){case"zigzag":u[s]||(u[s]=new h(e));break;case"cross-hatch":u[s]||(u[s]=new r(e));break;case"dots":u[s]||(u[s]=new i(e));break;case"dashed":u[s]||(u[s]=new c(e));break;case"zigzag-line":u[s]||(u[s]=new l(e));break;default:s="hachure",u[s]||(u[s]=new o(e))}return u[s]}(e,P).fillPolygons(t,e)}function z(t){const e=Object.assign({},t);return e.randomizer=void 0,t.seed&&(e.seed=t.seed+1),e}function W(t){return t.randomizer||(t.randomizer=new p(t.seed||0)),t.randomizer.next()}function E(t,e,s,n=1){return s.roughness*n*(W(s)*(e-t)+t)}function G(t,e,s=1){return E(-t,t,e,s)}function $(t,e,s,n,o,a=!1){const h=a?o.disableMultiStrokeFill:o.disableMultiStroke,r=R(t,e,s,n,o,!0,!1);if(h)return r;const i=R(t,e,s,n,o,!0,!0);return r.concat(i)}function R(t,e,s,n,o,a,h){const r=Math.pow(t-s,2)+Math.pow(e-n,2),i=Math.sqrt(r);let c=1;c=i<200?1:i>500?.4:-.0016668*i+1.233334;let l=o.maxRandomnessOffset||0;l*l*100>r&&(l=i/10);const u=l/2,p=.2+.2*W(o);let f=o.bowing*o.maxRandomnessOffset*(n-e)/200,d=o.bowing*o.maxRandomnessOffset*(t-s)/200;f=G(f,o,c),d=G(d,o,c);const g=[],M=()=>G(u,o,c),k=()=>G(l,o,c),b=o.preserveVertices;return a&&(h?g.push({op:"move",data:[t+(b?0:M()),e+(b?0:M())]}):g.push({op:"move",data:[t+(b?0:G(l,o,c)),e+(b?0:G(l,o,c))]})),h?g.push({op:"bcurveTo",data:[f+t+(s-t)*p+M(),d+e+(n-e)*p+M(),f+t+2*(s-t)*p+M(),d+e+2*(n-e)*p+M(),s+(b?0:M()),n+(b?0:M())]}):g.push({op:"bcurveTo",data:[f+t+(s-t)*p+k(),d+e+(n-e)*p+k(),f+t+2*(s-t)*p+k(),d+e+2*(n-e)*p+k(),s+(b?0:k()),n+(b?0:k())]}),g}function j(t,e,s){if(!t.length)return[];const n=[];n.push([t[0][0]+G(e,s),t[0][1]+G(e,s)]),n.push([t[0][0]+G(e,s),t[0][1]+G(e,s)]);for(let o=1;o<t.length;o++)n.push([t[o][0]+G(e,s),t[o][1]+G(e,s)]),o===t.length-1&&n.push([t[o][0]+G(e,s),t[o][1]+G(e,s)]);return q(n,null,s)}function q(t,e,s){const n=t.length,o=[];if(n>3){const a=[],h=1-s.curveTightness;o.push({op:"move",data:[t[1][0],t[1][1]]});for(let e=1;e+2<n;e++){const s=t[e];a[0]=[s[0],s[1]],a[1]=[s[0]+(h*t[e+1][0]-h*t[e-1][0])/6,s[1]+(h*t[e+1][1]-h*t[e-1][1])/6],a[2]=[t[e+1][0]+(h*t[e][0]-h*t[e+2][0])/6,t[e+1][1]+(h*t[e][1]-h*t[e+2][1])/6],a[3]=[t[e+1][0],t[e+1][1]],o.push({op:"bcurveTo",data:[a[1][0],a[1][1],a[2][0],a[2][1],a[3][0],a[3][1]]})}if(e&&2===e.length){const t=s.maxRandomnessOffset;o.push({op:"lineTo",data:[e[0]+G(t,s),e[1]+G(t,s)]})}}else 3===n?(o.push({op:"move",data:[t[1][0],t[1][1]]}),o.push({op:"bcurveTo",data:[t[1][0],t[1][1],t[2][0],t[2][1],t[2][0],t[2][1]]})):2===n&&o.push(...R(t[0][0],t[0][1],t[1][0],t[1][1],s,!0,!0));return o}function F(t,e,s,n,o,a,h,r){const i=[],c=[];if(0===r.roughness){t/=4,c.push([e+n*Math.cos(-t),s+o*Math.sin(-t)]);for(let a=0;a<=2*Math.PI;a+=t){const t=[e+n*Math.cos(a),s+o*Math.sin(a)];i.push(t),c.push(t)}c.push([e+n*Math.cos(0),s+o*Math.sin(0)]),c.push([e+n*Math.cos(t),s+o*Math.sin(t)])}else{const l=G(.5,r)-Math.PI/2;c.push([G(a,r)+e+.9*n*Math.cos(l-t),G(a,r)+s+.9*o*Math.sin(l-t)]);const u=2*Math.PI+l-.01;for(let h=l;h<u;h+=t){const t=[G(a,r)+e+n*Math.cos(h),G(a,r)+s+o*Math.sin(h)];i.push(t),c.push(t)}c.push([G(a,r)+e+n*Math.cos(l+2*Math.PI+.5*h),G(a,r)+s+o*Math.sin(l+2*Math.PI+.5*h)]),c.push([G(a,r)+e+.98*n*Math.cos(l+h),G(a,r)+s+.98*o*Math.sin(l+h)]),c.push([G(a,r)+e+.9*n*Math.cos(l+.5*h),G(a,r)+s+.9*o*Math.sin(l+.5*h)])}return[c,i]}function V(t,e,s,n,o,a,h,r,i){const c=a+G(.1,i),l=[];l.push([G(r,i)+e+.9*n*Math.cos(c-t),G(r,i)+s+.9*o*Math.sin(c-t)]);for(let a=c;a<=h;a+=t)l.push([G(r,i)+e+n*Math.cos(a),G(r,i)+s+o*Math.sin(a)]);return l.push([e+n*Math.cos(h),s+o*Math.sin(h)]),l.push([e+n*Math.cos(h),s+o*Math.sin(h)]),q(l,null,i)}function Z(t,e,s,n,o,a,h,r){const i=[],c=[r.maxRandomnessOffset||1,(r.maxRandomnessOffset||1)+.3];let l=[0,0];const u=r.disableMultiStroke?1:2,p=r.preserveVertices;for(let f=0;f<u;f++)0===f?i.push({op:"move",data:[h[0],h[1]]}):i.push({op:"move",data:[h[0]+(p?0:G(c[0],r)),h[1]+(p?0:G(c[0],r))]}),l=p?[o,a]:[o+G(c[f],r),a+G(c[f],r)],i.push({op:"bcurveTo",data:[t+G(c[f],r),e+G(c[f],r),s+G(c[f],r),n+G(c[f],r),l[0],l[1]]});return i}function Q(t){return[...t]}function H(t,e=0){const s=t.length;if(s<3)throw new Error("A curve must have at least three points.");const n=[];if(3===s)n.push(Q(t[0]),Q(t[1]),Q(t[2]),Q(t[2]));else{const s=[];s.push(t[0],t[0]);for(let e=1;e<t.length;e++)s.push(t[e]),e===t.length-1&&s.push(t[e]);const o=[],a=1-e;n.push(Q(s[0]));for(let t=1;t+2<s.length;t++){const e=s[t];o[0]=[e[0],e[1]],o[1]=[e[0]+(a*s[t+1][0]-a*s[t-1][0])/6,e[1]+(a*s[t+1][1]-a*s[t-1][1])/6],o[2]=[s[t+1][0]+(a*s[t][0]-a*s[t+2][0])/6,s[t+1][1]+(a*s[t][1]-a*s[t+2][1])/6],o[3]=[s[t+1][0],s[t+1][1]],n.push(o[1],o[2],o[3])}}return n}function N(t,e){return Math.pow(t[0]-e[0],2)+Math.pow(t[1]-e[1],2)}function B(t,e,s){const n=N(e,s);if(0===n)return N(t,e);let o=((t[0]-e[0])*(s[0]-e[0])+(t[1]-e[1])*(s[1]-e[1]))/n;return o=Math.max(0,Math.min(1,o)),N(t,J(e,s,o))}function J(t,e,s){return[t[0]+(e[0]-t[0])*s,t[1]+(e[1]-t[1])*s]}function K(t,e,s,n){const o=n||[];if(function(t,e){const s=t[e+0],n=t[e+1],o=t[e+2],a=t[e+3];let h=3*n[0]-2*s[0]-a[0];h*=h;let r=3*n[1]-2*s[1]-a[1];r*=r;let i=3*o[0]-2*a[0]-s[0];i*=i;let c=3*o[1]-2*a[1]-s[1];return c*=c,h<i&&(h=i),r<c&&(r=c),h+r}(t,e)<s){const s=t[e+0];if(o.length){(a=o[o.length-1],h=s,Math.sqrt(N(a,h)))>1&&o.push(s)}else o.push(s);o.push(t[e+3])}else{const n=.5,a=t[e+0],h=t[e+1],r=t[e+2],i=t[e+3],c=J(a,h,n),l=J(h,r,n),u=J(r,i,n),p=J(c,l,n),f=J(l,u,n),d=J(p,f,n);K([a,c,p,d],0,s,o),K([d,f,u,i],0,s,o)}var a,h;return o}function U(t,e){return X(t,0,t.length,e)}function X(t,e,s,n,o){const a=o||[],h=t[e],r=t[s-1];let i=0,c=1;for(let n=e+1;n<s-1;++n){const e=B(t[n],h,r);e>i&&(i=e,c=n)}return Math.sqrt(i)>n?(X(t,e,c+1,n,a),X(t,c,s,n,a)):(a.length||a.push(h),a.push(r)),a}function Y(t,e=.15,s){const n=[],o=(t.length-1)/3;for(let s=0;s<o;s++){K(t,3*s,e,n)}return s&&s>0?X(n,0,n.length,s):n}const tt="none";class et{constructor(t){this.defaultOptions={maxRandomnessOffset:2,roughness:1,bowing:1,stroke:"#000",strokeWidth:1,curveTightness:0,curveFitting:.95,curveStepCount:9,fillStyle:"hachure",fillWeight:-1,hachureAngle:-41,hachureGap:-1,dashOffset:-1,dashGap:-1,zigzagOffset:-1,seed:0,disableMultiStroke:!1,disableMultiStrokeFill:!1,preserveVertices:!1,fillShapeRoughnessGain:.8},this.config=t||{},this.config.options&&(this.defaultOptions=this._o(this.config.options))}static newSeed(){return Math.floor(Math.random()*2**31)}_o(t){return t?Object.assign({},this.defaultOptions,t):this.defaultOptions}_d(t,e,s){return{shape:t,sets:e||[],options:s||this.defaultOptions}}line(t,e,s,n,o){const a=this._o(o);return this._d("line",[v(t,e,s,n,a)],a)}rectangle(t,e,s,n,o){const a=this._o(o),h=[],r=O(t,e,s,n,a);if(a.fill){const o=[[t,e],[t+s,e],[t+s,e+n],[t,e+n]];"solid"===a.fillStyle?h.push(I([o],a)):h.push(C([o],a))}return a.stroke!==tt&&h.push(r),this._d("rectangle",h,a)}ellipse(t,e,s,n,o){const a=this._o(o),h=[],r=T(s,n,a),i=D(t,e,a,r);if(a.fill)if("solid"===a.fillStyle){const s=D(t,e,a,r).opset;s.type="fillPath",h.push(s)}else h.push(C([i.estimatedPoints],a));return a.stroke!==tt&&h.push(i.opset),this._d("ellipse",h,a)}circle(t,e,s,n){const o=this.ellipse(t,e,s,s,n);return o.shape="circle",o}linearPath(t,e){const s=this._o(e);return this._d("linearPath",[S(t,!1,s)],s)}arc(t,e,s,n,o,a,h=!1,r){const i=this._o(r),c=[],l=A(t,e,s,n,o,a,h,!0,i);if(h&&i.fill)if("solid"===i.fillStyle){const h=Object.assign({},i);h.disableMultiStroke=!0;const r=A(t,e,s,n,o,a,!0,!1,h);r.type="fillPath",c.push(r)}else c.push(function(t,e,s,n,o,a,h){const r=t,i=e;let c=Math.abs(s/2),l=Math.abs(n/2);c+=G(.01*c,h),l+=G(.01*l,h);let u=o,p=a;for(;u<0;)u+=2*Math.PI,p+=2*Math.PI;p-u>2*Math.PI&&(u=0,p=2*Math.PI);const f=(p-u)/h.curveStepCount,d=[];for(let t=u;t<=p;t+=f)d.push([r+c*Math.cos(t),i+l*Math.sin(t)]);return d.push([r+c*Math.cos(p),i+l*Math.sin(p)]),d.push([r,i]),C([d],h)}(t,e,s,n,o,a,i));return i.stroke!==tt&&c.push(l),this._d("arc",c,i)}curve(t,e){const s=this._o(e),n=[],o=L(t,s);if(s.fill&&s.fill!==tt)if("solid"===s.fillStyle){const e=L(t,Object.assign(Object.assign({},s),{disableMultiStroke:!0,roughness:s.roughness?s.roughness+s.fillShapeRoughnessGain:0}));n.push({type:"fillPath",ops:this._mergedShape(e.ops)})}else{const e=[],o=t;if(o.length){const t="number"==typeof o[0][0]?[o]:o;for(const n of t)n.length<3?e.push(...n):3===n.length?e.push(...Y(H([n[0],n[0],n[1],n[2]]),10,(1+s.roughness)/2)):e.push(...Y(H(n),10,(1+s.roughness)/2))}e.length&&n.push(C([e],s))}return s.stroke!==tt&&n.push(o),this._d("curve",n,s)}polygon(t,e){const s=this._o(e),n=[],o=S(t,!0,s);return s.fill&&("solid"===s.fillStyle?n.push(I([t],s)):n.push(C([t],s))),s.stroke!==tt&&n.push(o),this._d("polygon",n,s)}path(t,e){const s=this._o(e),n=[];if(!t)return this._d("path",n,s);t=(t||"").replace(/\n/g," ").replace(/(-\s)/g,"-").replace("/(ss)/g"," ");const o=s.fill&&"transparent"!==s.fill&&s.fill!==tt,a=s.stroke!==tt,h=!!(s.simplification&&s.simplification<1),r=function(t,e,s){const n=m(y(b(t))),o=[];let a=[],h=[0,0],r=[];const i=()=>{r.length>=4&&a.push(...Y(r,e)),r=[]},c=()=>{i(),a.length&&(o.push(a),a=[])};for(const{key:t,data:e}of n)switch(t){case"M":c(),h=[e[0],e[1]],a.push(h);break;case"L":i(),a.push([e[0],e[1]]);break;case"C":if(!r.length){const t=a.length?a[a.length-1]:h;r.push([t[0],t[1]])}r.push([e[0],e[1]]),r.push([e[2],e[3]]),r.push([e[4],e[5]]);break;case"Z":i(),a.push([h[0],h[1]])}if(c(),!s)return o;const l=[];for(const t of o){const e=U(t,s);e.length&&l.push(e)}return l}(t,1,h?4-4*(s.simplification||1):(1+s.roughness)/2),i=_(t,s);if(o)if("solid"===s.fillStyle)if(1===r.length){const e=_(t,Object.assign(Object.assign({},s),{disableMultiStroke:!0,roughness:s.roughness?s.roughness+s.fillShapeRoughnessGain:0}));n.push({type:"fillPath",ops:this._mergedShape(e.ops)})}else n.push(I(r,s));else n.push(C(r,s));return a&&(h?r.forEach((t=>{n.push(S(t,!1,s))})):n.push(i)),this._d("path",n,s)}opsToPath(t,e){let s="";for(const n of t.ops){const t="number"==typeof e&&e>=0?n.data.map((t=>+t.toFixed(e))):n.data;switch(n.op){case"move":s+=`M${t[0]} ${t[1]} `;break;case"bcurveTo":s+=`C${t[0]} ${t[1]}, ${t[2]} ${t[3]}, ${t[4]} ${t[5]} `;break;case"lineTo":s+=`L${t[0]} ${t[1]} `}}return s.trim()}toPaths(t){const e=t.sets||[],s=t.options||this.defaultOptions,n=[];for(const t of e){let e=null;switch(t.type){case"path":e={d:this.opsToPath(t),stroke:s.stroke,strokeWidth:s.strokeWidth,fill:tt};break;case"fillPath":e={d:this.opsToPath(t),stroke:tt,strokeWidth:0,fill:s.fill||tt};break;case"fillSketch":e=this.fillSketch(t,s)}e&&n.push(e)}return n}fillSketch(t,e){let s=e.fillWeight;return s<0&&(s=e.strokeWidth/2),{d:this.opsToPath(t),stroke:e.fill||tt,strokeWidth:s,fill:tt}}_mergedShape(t){return t.filter(((t,e)=>0===e||"move"!==t.op))}}class st{constructor(t,e){this.canvas=t,this.ctx=this.canvas.getContext("2d"),this.gen=new et(e)}draw(t){const e=t.sets||[],s=t.options||this.getDefaultOptions(),n=this.ctx,o=t.options.fixedDecimalPlaceDigits;for(const a of e)switch(a.type){case"path":n.save(),n.strokeStyle="none"===s.stroke?"transparent":s.stroke,n.lineWidth=s.strokeWidth,s.strokeLineDash&&n.setLineDash(s.strokeLineDash),s.strokeLineDashOffset&&(n.lineDashOffset=s.strokeLineDashOffset),this._drawToContext(n,a,o),n.restore();break;case"fillPath":{n.save(),n.fillStyle=s.fill||"";const e="curve"===t.shape||"polygon"===t.shape||"path"===t.shape?"evenodd":"nonzero";this._drawToContext(n,a,o,e),n.restore();break}case"fillSketch":this.fillSketch(n,a,s)}}fillSketch(t,e,s){let n=s.fillWeight;n<0&&(n=s.strokeWidth/2),t.save(),s.fillLineDash&&t.setLineDash(s.fillLineDash),s.fillLineDashOffset&&(t.lineDashOffset=s.fillLineDashOffset),t.strokeStyle=s.fill||"",t.lineWidth=n,this._drawToContext(t,e,s.fixedDecimalPlaceDigits),t.restore()}_drawToContext(t,e,s,n="nonzero"){t.beginPath();for(const n of e.ops){const e="number"==typeof s&&s>=0?n.data.map((t=>+t.toFixed(s))):n.data;switch(n.op){case"move":t.moveTo(e[0],e[1]);break;case"bcurveTo":t.bezierCurveTo(e[0],e[1],e[2],e[3],e[4],e[5]);break;case"lineTo":t.lineTo(e[0],e[1])}}"fillPath"===e.type?t.fill(n):t.stroke()}get generator(){return this.gen}getDefaultOptions(){return this.gen.defaultOptions}line(t,e,s,n,o){const a=this.gen.line(t,e,s,n,o);return this.draw(a),a}rectangle(t,e,s,n,o){const a=this.gen.rectangle(t,e,s,n,o);return this.draw(a),a}ellipse(t,e,s,n,o){const a=this.gen.ellipse(t,e,s,n,o);return this.draw(a),a}circle(t,e,s,n){const o=this.gen.circle(t,e,s,n);return this.draw(o),o}linearPath(t,e){const s=this.gen.linearPath(t,e);return this.draw(s),s}polygon(t,e){const s=this.gen.polygon(t,e);return this.draw(s),s}arc(t,e,s,n,o,a,h=!1,r){const i=this.gen.arc(t,e,s,n,o,a,h,r);return this.draw(i),i}curve(t,e){const s=this.gen.curve(t,e);return this.draw(s),s}path(t,e){const s=this.gen.path(t,e);return this.draw(s),s}}const nt="http://www.w3.org/2000/svg";class ot{constructor(t,e){this.svg=t,this.gen=new et(e)}draw(t){const e=t.sets||[],s=t.options||this.getDefaultOptions(),n=this.svg.ownerDocument||window.document,o=n.createElementNS(nt,"g"),a=t.options.fixedDecimalPlaceDigits;for(const h of e){let e=null;switch(h.type){case"path":e=n.createElementNS(nt,"path"),e.setAttribute("d",this.opsToPath(h,a)),e.setAttribute("stroke",s.stroke),e.setAttribute("stroke-width",s.strokeWidth+""),e.setAttribute("fill","none"),s.strokeLineDash&&e.setAttribute("stroke-dasharray",s.strokeLineDash.join(" ").trim()),s.strokeLineDashOffset&&e.setAttribute("stroke-dashoffset",`${s.strokeLineDashOffset}`);break;case"fillPath":e=n.createElementNS(nt,"path"),e.setAttribute("d",this.opsToPath(h,a)),e.setAttribute("stroke","none"),e.setAttribute("stroke-width","0"),e.setAttribute("fill",s.fill||""),"curve"!==t.shape&&"polygon"!==t.shape||e.setAttribute("fill-rule","evenodd");break;case"fillSketch":e=this.fillSketch(n,h,s)}e&&o.appendChild(e)}return o}fillSketch(t,e,s){let n=s.fillWeight;n<0&&(n=s.strokeWidth/2);const o=t.createElementNS(nt,"path");return o.setAttribute("d",this.opsToPath(e,s.fixedDecimalPlaceDigits)),o.setAttribute("stroke",s.fill||""),o.setAttribute("stroke-width",n+""),o.setAttribute("fill","none"),s.fillLineDash&&o.setAttribute("stroke-dasharray",s.fillLineDash.join(" ").trim()),s.fillLineDashOffset&&o.setAttribute("stroke-dashoffset",`${s.fillLineDashOffset}`),o}get generator(){return this.gen}getDefaultOptions(){return this.gen.defaultOptions}opsToPath(t,e){return this.gen.opsToPath(t,e)}line(t,e,s,n,o){const a=this.gen.line(t,e,s,n,o);return this.draw(a)}rectangle(t,e,s,n,o){const a=this.gen.rectangle(t,e,s,n,o);return this.draw(a)}ellipse(t,e,s,n,o){const a=this.gen.ellipse(t,e,s,n,o);return this.draw(a)}circle(t,e,s,n){const o=this.gen.circle(t,e,s,n);return this.draw(o)}linearPath(t,e){const s=this.gen.linearPath(t,e);return this.draw(s)}polygon(t,e){const s=this.gen.polygon(t,e);return this.draw(s)}arc(t,e,s,n,o,a,h=!1,r){const i=this.gen.arc(t,e,s,n,o,a,h,r);return this.draw(i)}curve(t,e){const s=this.gen.curve(t,e);return this.draw(s)}path(t,e){const s=this.gen.path(t,e);return this.draw(s)}}var at={canvas:(t,e)=>new st(t,e),svg:(t,e)=>new ot(t,e),generator:t=>new et(t),newSeed:()=>et.newSeed()};


/***/ }),

/***/ 11464:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (/* binding */ dedent)
/* harmony export */ });
function dedent(templ) {
    var values = [];
    for (var _i = 1; _i < arguments.length; _i++) {
        values[_i - 1] = arguments[_i];
    }
    var strings = Array.from(typeof templ === 'string' ? [templ] : templ);
    strings[strings.length - 1] = strings[strings.length - 1].replace(/\r?\n([\t ]*)$/, '');
    var indentLengths = strings.reduce(function (arr, str) {
        var matches = str.match(/\n([\t ]+|(?!\s).)/g);
        if (matches) {
            return arr.concat(matches.map(function (match) { var _a, _b; return (_b = (_a = match.match(/[\t ]/g)) === null || _a === void 0 ? void 0 : _a.length) !== null && _b !== void 0 ? _b : 0; }));
        }
        return arr;
    }, []);
    if (indentLengths.length) {
        var pattern_1 = new RegExp("\n[\t ]{" + Math.min.apply(Math, indentLengths) + "}", 'g');
        strings = strings.map(function (str) { return str.replace(pattern_1, '\n'); });
    }
    strings[0] = strings[0].replace(/^\r?\n/, '');
    var string = strings[0];
    values.forEach(function (value, i) {
        var endentations = string.match(/(?:^|\n)( *)$/);
        var endentation = endentations ? endentations[1] : '';
        var indentedValue = value;
        if (typeof value === 'string' && value.includes('\n')) {
            indentedValue = String(value)
                .split('\n')
                .map(function (str, i) {
                return i === 0 ? str : "" + endentation + str;
            })
                .join('\n');
        }
        string += indentedValue + strings[i + 1];
    });
    return string;
}
/* unused harmony default export */ var __WEBPACK_DEFAULT_EXPORT__ = ((/* unused pure expression or super */ null && (dedent)));
//# sourceMappingURL=index.js.map

/***/ })

}]);
//# sourceMappingURL=4406.b99d8af86330961ec7e2.js.map?v=b99d8af86330961ec7e2
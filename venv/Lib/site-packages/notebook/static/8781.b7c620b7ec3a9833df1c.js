(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[8781],{

/***/ 67417:
/***/ (() => {



/***/ }),

/***/ 60880:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
__webpack_require__.r(__webpack_exports__);
/* harmony import */ var _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(79761);
/* harmony import */ var _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_0__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.

// Inspired by: https://github.com/jupyterlab/jupyterlab/blob/master/dev_mode/index.js



__webpack_require__(72761);
__webpack_require__(67417);

function loadScript(url) {
  return new Promise((resolve, reject) => {
    const newScript = document.createElement('script');
    newScript.onerror = reject;
    newScript.onload = resolve;
    newScript.async = true;
    document.head.appendChild(newScript);
    newScript.src = url;
  });
}
async function loadComponent(url, scope) {
  await loadScript(url);

  // From MIT-licensed https://github.com/module-federation/module-federation-examples/blob/af043acd6be1718ee195b2511adf6011fba4233c/advanced-api/dynamic-remotes/app1/src/App.js#L6-L12
  // eslint-disable-next-line no-undef
  await __webpack_require__.I('default');
  const container = window._JUPYTERLAB[scope];
  // Initialize the container, it may provide shared modules and may need ours
  // eslint-disable-next-line no-undef
  await container.init(__webpack_require__.S.default);
}

async function createModule(scope, module) {
  try {
    const factory = await window._JUPYTERLAB[scope].get(module);
    const instance = factory();
    instance.__scope__ = scope;
    return instance;
  } catch (e) {
    console.warn(
      `Failed to create module: package: ${scope}; module: ${module}`
    );
    throw e;
  }
}

/**
 * The main function
 */
async function main() {
  const mimeExtensionsMods = [
    __webpack_require__(91555),
    __webpack_require__(21599),
    __webpack_require__(31708),
    __webpack_require__(68552),
  ];
  const mimeExtensions = await Promise.all(mimeExtensionsMods);

  // Load the base plugins available on all pages
  let baseMods = [
      __webpack_require__(7249),
  __webpack_require__(93630),
  __webpack_require__(66968),
  __webpack_require__(58215),
  __webpack_require__(28027),
  __webpack_require__(14459),
  __webpack_require__(22414),
  
      (__webpack_require__(63917)["default"].filter)(({id}) => [
       '@jupyterlab/application-extension:commands',
'@jupyterlab/application-extension:context-menu',
'@jupyterlab/application-extension:faviconbusy',
'@jupyterlab/application-extension:router',
'@jupyterlab/application-extension:top-bar',
'@jupyterlab/application-extension:top-spacer',
      ].includes(id)),
      
      (__webpack_require__(14736)["default"].filter)(({id}) => [
       '@jupyterlab/apputils-extension:palette',
'@jupyterlab/apputils-extension:notification',
'@jupyterlab/apputils-extension:sanitizer',
'@jupyterlab/apputils-extension:sessionDialogs',
'@jupyterlab/apputils-extension:settings',
'@jupyterlab/apputils-extension:state',
'@jupyterlab/apputils-extension:themes',
'@jupyterlab/apputils-extension:themes-palette-menu',
'@jupyterlab/apputils-extension:toolbar-registry',
'@jupyterlab/apputils-extension:utilityCommands',
      ].includes(id)),
      __webpack_require__(14914),
  
      (__webpack_require__(34351)["default"].filter)(({id}) => [
       '@jupyterlab/completer-extension:base-service',
'@jupyterlab/completer-extension:inline-completer',
'@jupyterlab/completer-extension:inline-completer-factory',
'@jupyterlab/completer-extension:inline-history',
'@jupyterlab/completer-extension:manager',
      ].includes(id)),
      
      (__webpack_require__(69610)["default"].filter)(({id}) => [
       '@jupyterlab/console-extension:cell-executor',
'@jupyterlab/console-extension:completer',
'@jupyterlab/console-extension:factory',
'@jupyterlab/console-extension:foreign',
'@jupyterlab/console-extension:tracker',
      ].includes(id)),
      __webpack_require__(68889),
  
      (__webpack_require__(61047)["default"].filter)(({id}) => [
       '@jupyterlab/docmanager-extension:plugin',
'@jupyterlab/docmanager-extension:download',
'@jupyterlab/docmanager-extension:contexts',
'@jupyterlab/docmanager-extension:manager',
      ].includes(id)),
      
      (__webpack_require__(53019)["default"].filter)(({id}) => [
       '@jupyterlab/documentsearch-extension:plugin',
      ].includes(id)),
      
      (__webpack_require__(52182)["default"].filter)(({id}) => [
       '@jupyterlab/filebrowser-extension:factory',
'@jupyterlab/filebrowser-extension:default-file-browser',
      ].includes(id)),
      
      (__webpack_require__(91122)["default"].filter)(({id}) => [
       '@jupyterlab/fileeditor-extension:plugin',
      ].includes(id)),
      
      (__webpack_require__(66106)["default"].filter)(({id}) => [
       '@jupyterlab/help-extension:resources',
      ].includes(id)),
      __webpack_require__(64778),
  __webpack_require__(12300),
  __webpack_require__(21495),
  
      (__webpack_require__(69010)["default"].filter)(({id}) => [
       '@jupyterlab/mainmenu-extension:plugin',
      ].includes(id)),
      __webpack_require__(38225),
  __webpack_require__(37924),
  __webpack_require__(12634),
  
      (__webpack_require__(7925)["default"].filter)(({id}) => [
       '@jupyterlab/notebook-extension:cell-executor',
'@jupyterlab/notebook-extension:code-console',
'@jupyterlab/notebook-extension:export',
'@jupyterlab/notebook-extension:factory',
'@jupyterlab/notebook-extension:tracker',
'@jupyterlab/notebook-extension:widget-factory',
      ].includes(id)),
      __webpack_require__(19484),
  __webpack_require__(12293),
  __webpack_require__(36791),
  __webpack_require__(9074),
  __webpack_require__(40103),
  __webpack_require__(52880),
  __webpack_require__(5437),
  __webpack_require__(74469),
  __webpack_require__(11326),
  
  ];

  const page = `/${_jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_0__.PageConfig.getOption('notebookPage')}`;
  switch (page) {
    // list all the other plugins grouped by page
    case '/tree': {
      baseMods = baseMods.concat([
        __webpack_require__(27787),
  
      (__webpack_require__(52182)["default"].filter)(({id}) => [
       '@jupyterlab/filebrowser-extension:browser',
'@jupyterlab/filebrowser-extension:download',
'@jupyterlab/filebrowser-extension:file-upload-status',
'@jupyterlab/filebrowser-extension:open-with',
'@jupyterlab/filebrowser-extension:search',
'@jupyterlab/filebrowser-extension:share-file',
      ].includes(id)),
      __webpack_require__(23003),
  
      (__webpack_require__(78352)["default"].filter)(({id}) => [
       '@jupyterlab/running-extension:plugin',
      ].includes(id)),
      __webpack_require__(47102),
  
      ]);
      break;
    }
    // list all the other plugins grouped by page
    case '/notebooks': {
      baseMods = baseMods.concat([
        __webpack_require__(1286),
  __webpack_require__(61413),
  
      (__webpack_require__(28634)["default"].filter)(({id}) => [
       '@jupyterlab/debugger-extension:config',
'@jupyterlab/debugger-extension:main',
'@jupyterlab/debugger-extension:notebooks',
'@jupyterlab/debugger-extension:service',
'@jupyterlab/debugger-extension:sidebar',
'@jupyterlab/debugger-extension:sources',
      ].includes(id)),
      __webpack_require__(51693),
  
      (__webpack_require__(7925)["default"].filter)(({id}) => [
       '@jupyterlab/notebook-extension:active-cell-tool',
'@jupyterlab/notebook-extension:completer',
'@jupyterlab/notebook-extension:copy-output',
'@jupyterlab/notebook-extension:metadata-editor',
'@jupyterlab/notebook-extension:search',
'@jupyterlab/notebook-extension:toc',
'@jupyterlab/notebook-extension:tools',
'@jupyterlab/notebook-extension:update-raw-mimetype',
      ].includes(id)),
      
      (__webpack_require__(5947)["default"].filter)(({id}) => [
       '@jupyterlab/toc-extension:registry',
'@jupyterlab/toc-extension:tracker',
      ].includes(id)),
      
      (__webpack_require__(74412)["default"].filter)(({id}) => [
       '@jupyterlab/tooltip-extension:manager',
'@jupyterlab/tooltip-extension:notebooks',
      ].includes(id)),
      
      ]);
      break;
    }
    // list all the other plugins grouped by page
    case '/consoles': {
      baseMods = baseMods.concat([
        
      (__webpack_require__(74412)["default"].filter)(({id}) => [
       '@jupyterlab/tooltip-extension:manager',
'@jupyterlab/tooltip-extension:consoles',
      ].includes(id)),
      
      ]);
      break;
    }
    // list all the other plugins grouped by page
    case '/edit': {
      baseMods = baseMods.concat([
        
      (__webpack_require__(91122)["default"].filter)(({id}) => [
       '@jupyterlab/fileeditor-extension:completer',
'@jupyterlab/fileeditor-extension:search',
      ].includes(id)),
      __webpack_require__(57676),
  
      ]);
      break;
    }
  }

  // populate the list of disabled extensions
  const disabled = [];
  const availablePlugins = [];

  /**
   * Iterate over active plugins in an extension.
   *
   * #### Notes
   * This also populates the disabled
   */
  function* activePlugins(extension) {
    // Handle commonjs or es2015 modules
    let exports;
    if (Object.prototype.hasOwnProperty.call(extension, '__esModule')) {
      exports = extension.default;
    } else {
      // CommonJS exports.
      exports = extension;
    }

    let plugins = Array.isArray(exports) ? exports : [exports];
    for (let plugin of plugins) {
      const isDisabled = _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_0__.PageConfig.Extension.isDisabled(plugin.id);
      availablePlugins.push({
        id: plugin.id,
        description: plugin.description,
        requires: plugin.requires ?? [],
        optional: plugin.optional ?? [],
        provides: plugin.provides ?? null,
        autoStart: plugin.autoStart,
        enabled: !isDisabled,
        extension: extension.__scope__
      });
      if (isDisabled) {
        disabled.push(plugin.id);
        continue;
      }
      yield plugin;
    }
  }

  const extension_data = JSON.parse(
    _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_0__.PageConfig.getOption('federated_extensions')
  );

  const mods = [];
  const federatedExtensionPromises = [];
  const federatedMimeExtensionPromises = [];
  const federatedStylePromises = [];

  const extensions = await Promise.allSettled(
    extension_data.map(async data => {
      await loadComponent(
        `${_jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_0__.URLExt.join(
          _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_0__.PageConfig.getOption('fullLabextensionsUrl'),
          data.name,
          data.load
        )}`,
        data.name
      );
      return data;
    })
  );

  extensions.forEach(p => {
    if (p.status === 'rejected') {
      // There was an error loading the component
      console.error(p.reason);
      return;
    }

    const data = p.value;
    if (data.extension) {
      federatedExtensionPromises.push(createModule(data.name, data.extension));
    }
    if (data.mimeExtension) {
      federatedMimeExtensionPromises.push(
        createModule(data.name, data.mimeExtension)
      );
    }
    if (data.style && !_jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_0__.PageConfig.Extension.isDisabled(data.name)) {
      federatedStylePromises.push(createModule(data.name, data.style));
    }
  });

  // Add the base frontend extensions
  const baseFrontendMods = await Promise.all(baseMods);
  baseFrontendMods.forEach(p => {
    for (let plugin of activePlugins(p)) {
      mods.push(plugin);
    }
  });

  // Add the federated extensions.
  const federatedExtensions = await Promise.allSettled(
    federatedExtensionPromises
  );
  federatedExtensions.forEach(p => {
    if (p.status === 'fulfilled') {
      for (let plugin of activePlugins(p.value)) {
        mods.push(plugin);
      }
    } else {
      console.error(p.reason);
    }
  });

  // Add the federated mime extensions.
  const federatedMimeExtensions = await Promise.allSettled(
    federatedMimeExtensionPromises
  );
  federatedMimeExtensions.forEach(p => {
    if (p.status === 'fulfilled') {
      for (let plugin of activePlugins(p.value)) {
        mimeExtensions.push(plugin);
      }
    } else {
      console.error(p.reason);
    }
  });

  // Load all federated component styles and log errors for any that do not
  (await Promise.allSettled(federatedStylePromises))
    .filter(({ status }) => status === 'rejected')
    .forEach(({ reason }) => {
      console.error(reason);
    });

  // Set the list of base notebook multi-page plugins so the app is aware of all
  // its built-in plugins even if they are not loaded on the current page.
  // For example this is useful so the Settings Editor can list the debugger
  // plugin even if the debugger is only loaded on the notebook page.
  _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_0__.PageConfig.setOption('allPlugins', '{"/":{"@jupyter-notebook/application-extension":true,"@jupyter-notebook/console-extension":true,"@jupyter-notebook/docmanager-extension":true,"@jupyter-notebook/documentsearch-extension":true,"@jupyter-notebook/help-extension":true,"@jupyter-notebook/notebook-extension":true,"@jupyter-notebook/terminal-extension":true,"@jupyterlab/application-extension":["@jupyterlab/application-extension:commands","@jupyterlab/application-extension:context-menu","@jupyterlab/application-extension:faviconbusy","@jupyterlab/application-extension:router","@jupyterlab/application-extension:top-bar","@jupyterlab/application-extension:top-spacer"],"@jupyterlab/apputils-extension":["@jupyterlab/apputils-extension:palette","@jupyterlab/apputils-extension:notification","@jupyterlab/apputils-extension:sanitizer","@jupyterlab/apputils-extension:sessionDialogs","@jupyterlab/apputils-extension:settings","@jupyterlab/apputils-extension:state","@jupyterlab/apputils-extension:themes","@jupyterlab/apputils-extension:themes-palette-menu","@jupyterlab/apputils-extension:toolbar-registry","@jupyterlab/apputils-extension:utilityCommands"],"@jupyterlab/codemirror-extension":true,"@jupyterlab/completer-extension":["@jupyterlab/completer-extension:base-service","@jupyterlab/completer-extension:inline-completer","@jupyterlab/completer-extension:inline-completer-factory","@jupyterlab/completer-extension:inline-history","@jupyterlab/completer-extension:manager"],"@jupyterlab/console-extension":["@jupyterlab/console-extension:cell-executor","@jupyterlab/console-extension:completer","@jupyterlab/console-extension:factory","@jupyterlab/console-extension:foreign","@jupyterlab/console-extension:tracker"],"@jupyterlab/csvviewer-extension":true,"@jupyterlab/docmanager-extension":["@jupyterlab/docmanager-extension:plugin","@jupyterlab/docmanager-extension:download","@jupyterlab/docmanager-extension:contexts","@jupyterlab/docmanager-extension:manager"],"@jupyterlab/documentsearch-extension":["@jupyterlab/documentsearch-extension:plugin"],"@jupyterlab/filebrowser-extension":["@jupyterlab/filebrowser-extension:factory","@jupyterlab/filebrowser-extension:default-file-browser"],"@jupyterlab/fileeditor-extension":["@jupyterlab/fileeditor-extension:plugin"],"@jupyterlab/help-extension":["@jupyterlab/help-extension:resources"],"@jupyterlab/htmlviewer-extension":true,"@jupyterlab/imageviewer-extension":true,"@jupyterlab/lsp-extension":true,"@jupyterlab/mainmenu-extension":["@jupyterlab/mainmenu-extension:plugin"],"@jupyterlab/markedparser-extension":true,"@jupyterlab/mathjax-extension":true,"@jupyterlab/mermaid-extension":true,"@jupyterlab/notebook-extension":["@jupyterlab/notebook-extension:cell-executor","@jupyterlab/notebook-extension:code-console","@jupyterlab/notebook-extension:export","@jupyterlab/notebook-extension:factory","@jupyterlab/notebook-extension:tracker","@jupyterlab/notebook-extension:widget-factory"],"@jupyterlab/pluginmanager-extension":true,"@jupyterlab/shortcuts-extension":true,"@jupyterlab/terminal-extension":true,"@jupyterlab/theme-light-extension":true,"@jupyterlab/theme-dark-extension":true,"@jupyterlab/theme-dark-high-contrast-extension":true,"@jupyterlab/translation-extension":true,"@jupyterlab/ui-components-extension":true,"@jupyterlab/hub-extension":true},"/tree":{"@jupyterlab/extensionmanager-extension":true,"@jupyterlab/filebrowser-extension":["@jupyterlab/filebrowser-extension:browser","@jupyterlab/filebrowser-extension:download","@jupyterlab/filebrowser-extension:file-upload-status","@jupyterlab/filebrowser-extension:open-with","@jupyterlab/filebrowser-extension:search","@jupyterlab/filebrowser-extension:share-file"],"@jupyter-notebook/tree-extension":true,"@jupyterlab/running-extension":["@jupyterlab/running-extension:plugin"],"@jupyterlab/settingeditor-extension":true},"/notebooks":{"@jupyterlab/celltags-extension":true,"@jupyterlab/cell-toolbar-extension":true,"@jupyterlab/debugger-extension":["@jupyterlab/debugger-extension:config","@jupyterlab/debugger-extension:main","@jupyterlab/debugger-extension:notebooks","@jupyterlab/debugger-extension:service","@jupyterlab/debugger-extension:sidebar","@jupyterlab/debugger-extension:sources"],"@jupyterlab/metadataform-extension":true,"@jupyterlab/notebook-extension":["@jupyterlab/notebook-extension:active-cell-tool","@jupyterlab/notebook-extension:completer","@jupyterlab/notebook-extension:copy-output","@jupyterlab/notebook-extension:metadata-editor","@jupyterlab/notebook-extension:search","@jupyterlab/notebook-extension:toc","@jupyterlab/notebook-extension:tools","@jupyterlab/notebook-extension:update-raw-mimetype"],"@jupyterlab/toc-extension":["@jupyterlab/toc-extension:registry","@jupyterlab/toc-extension:tracker"],"@jupyterlab/tooltip-extension":["@jupyterlab/tooltip-extension:manager","@jupyterlab/tooltip-extension:notebooks"]},"/consoles":{"@jupyterlab/tooltip-extension":["@jupyterlab/tooltip-extension:manager","@jupyterlab/tooltip-extension:consoles"]},"/edit":{"@jupyterlab/fileeditor-extension":["@jupyterlab/fileeditor-extension:completer","@jupyterlab/fileeditor-extension:search"],"@jupyterlab/markdownviewer-extension":true}}');

  const NotebookApp = (__webpack_require__(8864).NotebookApp);
  const app = new NotebookApp({ mimeExtensions, availablePlugins });

  app.registerPluginModules(mods);

  // Expose global app instance when in dev mode or when toggled explicitly.
  const exposeAppInBrowser =
    (_jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_0__.PageConfig.getOption('exposeAppInBrowser') || '').toLowerCase() === 'true';

  if (exposeAppInBrowser) {
    window.jupyterapp = app;
  }

  await app.start();
}

window.addEventListener('load', main);


/***/ }),

/***/ 72761:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
// ESM COMPAT FLAG
__webpack_require__.r(__webpack_exports__);

// EXTERNAL MODULE: ../node_modules/@jupyterlab/application/style/index.js + 1 modules
var style = __webpack_require__(67663);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/mainmenu/style/index.js
var mainmenu_style = __webpack_require__(70022);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/ui-components/style/index.js + 1 modules
var ui_components_style = __webpack_require__(26238);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/style-loader/dist/runtime/injectStylesIntoStyleTag.js
var injectStylesIntoStyleTag = __webpack_require__(94830);
var injectStylesIntoStyleTag_default = /*#__PURE__*/__webpack_require__.n(injectStylesIntoStyleTag);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/style-loader/dist/runtime/styleDomAPI.js
var styleDomAPI = __webpack_require__(80592);
var styleDomAPI_default = /*#__PURE__*/__webpack_require__.n(styleDomAPI);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/style-loader/dist/runtime/insertBySelector.js
var insertBySelector = __webpack_require__(99763);
var insertBySelector_default = /*#__PURE__*/__webpack_require__.n(insertBySelector);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/style-loader/dist/runtime/setAttributesWithoutAttributes.js
var setAttributesWithoutAttributes = __webpack_require__(28915);
var setAttributesWithoutAttributes_default = /*#__PURE__*/__webpack_require__.n(setAttributesWithoutAttributes);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/style-loader/dist/runtime/insertStyleElement.js
var insertStyleElement = __webpack_require__(80366);
var insertStyleElement_default = /*#__PURE__*/__webpack_require__.n(insertStyleElement);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/style-loader/dist/runtime/styleTagTransform.js
var styleTagTransform = __webpack_require__(17352);
var styleTagTransform_default = /*#__PURE__*/__webpack_require__.n(styleTagTransform);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/css-loader/dist/cjs.js!../packages/application/style/base.css
var base = __webpack_require__(26633);
;// CONCATENATED MODULE: ../packages/application/style/base.css

      
      
      
      
      
      
      
      
      

var options = {};

options.styleTagTransform = (styleTagTransform_default());
options.setAttributes = (setAttributesWithoutAttributes_default());

      options.insert = insertBySelector_default().bind(null, "head");
    
options.domAPI = (styleDomAPI_default());
options.insertStyleElement = (insertStyleElement_default());

var update = injectStylesIntoStyleTag_default()(base/* default */.Z, options);




       /* harmony default export */ const style_base = (base/* default */.Z && base/* default */.Z.locals ? base/* default */.Z.locals : undefined);

// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/css-loader/dist/cjs.js!../packages/application/style/sidepanel.css
var sidepanel = __webpack_require__(93862);
;// CONCATENATED MODULE: ../packages/application/style/sidepanel.css

      
      
      
      
      
      
      
      
      

var sidepanel_options = {};

sidepanel_options.styleTagTransform = (styleTagTransform_default());
sidepanel_options.setAttributes = (setAttributesWithoutAttributes_default());

      sidepanel_options.insert = insertBySelector_default().bind(null, "head");
    
sidepanel_options.domAPI = (styleDomAPI_default());
sidepanel_options.insertStyleElement = (insertStyleElement_default());

var sidepanel_update = injectStylesIntoStyleTag_default()(sidepanel/* default */.Z, sidepanel_options);




       /* harmony default export */ const style_sidepanel = (sidepanel/* default */.Z && sidepanel/* default */.Z.locals ? sidepanel/* default */.Z.locals : undefined);

;// CONCATENATED MODULE: ../packages/application/style/index.js
/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/








// EXTERNAL MODULE: ../node_modules/@lumino/widgets/style/index.js + 1 modules
var widgets_style = __webpack_require__(20959);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/css-loader/dist/cjs.js!../packages/application-extension/style/base.css
var application_extension_style_base = __webpack_require__(58027);
;// CONCATENATED MODULE: ../packages/application-extension/style/base.css

      
      
      
      
      
      
      
      
      

var base_options = {};

base_options.styleTagTransform = (styleTagTransform_default());
base_options.setAttributes = (setAttributesWithoutAttributes_default());

      base_options.insert = insertBySelector_default().bind(null, "head");
    
base_options.domAPI = (styleDomAPI_default());
base_options.insertStyleElement = (insertStyleElement_default());

var base_update = injectStylesIntoStyleTag_default()(application_extension_style_base/* default */.Z, base_options);




       /* harmony default export */ const packages_application_extension_style_base = (application_extension_style_base/* default */.Z && application_extension_style_base/* default */.Z.locals ? application_extension_style_base/* default */.Z.locals : undefined);

;// CONCATENATED MODULE: ../packages/application-extension/style/index.js





// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/css-loader/dist/cjs.js!../packages/console-extension/style/base.css
var console_extension_style_base = __webpack_require__(99776);
;// CONCATENATED MODULE: ../packages/console-extension/style/base.css

      
      
      
      
      
      
      
      
      

var style_base_options = {};

style_base_options.styleTagTransform = (styleTagTransform_default());
style_base_options.setAttributes = (setAttributesWithoutAttributes_default());

      style_base_options.insert = insertBySelector_default().bind(null, "head");
    
style_base_options.domAPI = (styleDomAPI_default());
style_base_options.insertStyleElement = (insertStyleElement_default());

var style_base_update = injectStylesIntoStyleTag_default()(console_extension_style_base/* default */.Z, style_base_options);




       /* harmony default export */ const packages_console_extension_style_base = (console_extension_style_base/* default */.Z && console_extension_style_base/* default */.Z.locals ? console_extension_style_base/* default */.Z.locals : undefined);

;// CONCATENATED MODULE: ../packages/console-extension/style/index.js


// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/css-loader/dist/cjs.js!../packages/docmanager-extension/style/base.css
var docmanager_extension_style_base = __webpack_require__(64399);
;// CONCATENATED MODULE: ../packages/docmanager-extension/style/base.css

      
      
      
      
      
      
      
      
      

var docmanager_extension_style_base_options = {};

docmanager_extension_style_base_options.styleTagTransform = (styleTagTransform_default());
docmanager_extension_style_base_options.setAttributes = (setAttributesWithoutAttributes_default());

      docmanager_extension_style_base_options.insert = insertBySelector_default().bind(null, "head");
    
docmanager_extension_style_base_options.domAPI = (styleDomAPI_default());
docmanager_extension_style_base_options.insertStyleElement = (insertStyleElement_default());

var docmanager_extension_style_base_update = injectStylesIntoStyleTag_default()(docmanager_extension_style_base/* default */.Z, docmanager_extension_style_base_options);




       /* harmony default export */ const packages_docmanager_extension_style_base = (docmanager_extension_style_base/* default */.Z && docmanager_extension_style_base/* default */.Z.locals ? docmanager_extension_style_base/* default */.Z.locals : undefined);

;// CONCATENATED MODULE: ../packages/docmanager-extension/style/index.js


// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/css-loader/dist/cjs.js!../packages/documentsearch-extension/style/base.css
var documentsearch_extension_style_base = __webpack_require__(11112);
;// CONCATENATED MODULE: ../packages/documentsearch-extension/style/base.css

      
      
      
      
      
      
      
      
      

var documentsearch_extension_style_base_options = {};

documentsearch_extension_style_base_options.styleTagTransform = (styleTagTransform_default());
documentsearch_extension_style_base_options.setAttributes = (setAttributesWithoutAttributes_default());

      documentsearch_extension_style_base_options.insert = insertBySelector_default().bind(null, "head");
    
documentsearch_extension_style_base_options.domAPI = (styleDomAPI_default());
documentsearch_extension_style_base_options.insertStyleElement = (insertStyleElement_default());

var documentsearch_extension_style_base_update = injectStylesIntoStyleTag_default()(documentsearch_extension_style_base/* default */.Z, documentsearch_extension_style_base_options);




       /* harmony default export */ const packages_documentsearch_extension_style_base = (documentsearch_extension_style_base/* default */.Z && documentsearch_extension_style_base/* default */.Z.locals ? documentsearch_extension_style_base/* default */.Z.locals : undefined);

;// CONCATENATED MODULE: ../packages/documentsearch-extension/style/index.js


// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/css-loader/dist/cjs.js!../packages/help-extension/style/base.css
var help_extension_style_base = __webpack_require__(49482);
;// CONCATENATED MODULE: ../packages/help-extension/style/base.css

      
      
      
      
      
      
      
      
      

var help_extension_style_base_options = {};

help_extension_style_base_options.styleTagTransform = (styleTagTransform_default());
help_extension_style_base_options.setAttributes = (setAttributesWithoutAttributes_default());

      help_extension_style_base_options.insert = insertBySelector_default().bind(null, "head");
    
help_extension_style_base_options.domAPI = (styleDomAPI_default());
help_extension_style_base_options.insertStyleElement = (insertStyleElement_default());

var help_extension_style_base_update = injectStylesIntoStyleTag_default()(help_extension_style_base/* default */.Z, help_extension_style_base_options);




       /* harmony default export */ const packages_help_extension_style_base = (help_extension_style_base/* default */.Z && help_extension_style_base/* default */.Z.locals ? help_extension_style_base/* default */.Z.locals : undefined);

;// CONCATENATED MODULE: ../packages/help-extension/style/index.js


// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/css-loader/dist/cjs.js!../packages/notebook-extension/style/base.css
var notebook_extension_style_base = __webpack_require__(71564);
;// CONCATENATED MODULE: ../packages/notebook-extension/style/base.css

      
      
      
      
      
      
      
      
      

var notebook_extension_style_base_options = {};

notebook_extension_style_base_options.styleTagTransform = (styleTagTransform_default());
notebook_extension_style_base_options.setAttributes = (setAttributesWithoutAttributes_default());

      notebook_extension_style_base_options.insert = insertBySelector_default().bind(null, "head");
    
notebook_extension_style_base_options.domAPI = (styleDomAPI_default());
notebook_extension_style_base_options.insertStyleElement = (insertStyleElement_default());

var notebook_extension_style_base_update = injectStylesIntoStyleTag_default()(notebook_extension_style_base/* default */.Z, notebook_extension_style_base_options);




       /* harmony default export */ const packages_notebook_extension_style_base = (notebook_extension_style_base/* default */.Z && notebook_extension_style_base/* default */.Z.locals ? notebook_extension_style_base/* default */.Z.locals : undefined);

;// CONCATENATED MODULE: ../packages/notebook-extension/style/index.js


// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/css-loader/dist/cjs.js!../packages/terminal-extension/style/base.css
var terminal_extension_style_base = __webpack_require__(8173);
;// CONCATENATED MODULE: ../packages/terminal-extension/style/base.css

      
      
      
      
      
      
      
      
      

var terminal_extension_style_base_options = {};

terminal_extension_style_base_options.styleTagTransform = (styleTagTransform_default());
terminal_extension_style_base_options.setAttributes = (setAttributesWithoutAttributes_default());

      terminal_extension_style_base_options.insert = insertBySelector_default().bind(null, "head");
    
terminal_extension_style_base_options.domAPI = (styleDomAPI_default());
terminal_extension_style_base_options.insertStyleElement = (insertStyleElement_default());

var terminal_extension_style_base_update = injectStylesIntoStyleTag_default()(terminal_extension_style_base/* default */.Z, terminal_extension_style_base_options);




       /* harmony default export */ const packages_terminal_extension_style_base = (terminal_extension_style_base/* default */.Z && terminal_extension_style_base/* default */.Z.locals ? terminal_extension_style_base/* default */.Z.locals : undefined);

;// CONCATENATED MODULE: ../packages/terminal-extension/style/index.js


// EXTERNAL MODULE: ../node_modules/@jupyterlab/filebrowser/style/index.js + 1 modules
var filebrowser_style = __webpack_require__(40360);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/css-loader/dist/cjs.js!../packages/tree/style/base.css
var tree_style_base = __webpack_require__(30546);
;// CONCATENATED MODULE: ../packages/tree/style/base.css

      
      
      
      
      
      
      
      
      

var tree_style_base_options = {};

tree_style_base_options.styleTagTransform = (styleTagTransform_default());
tree_style_base_options.setAttributes = (setAttributesWithoutAttributes_default());

      tree_style_base_options.insert = insertBySelector_default().bind(null, "head");
    
tree_style_base_options.domAPI = (styleDomAPI_default());
tree_style_base_options.insertStyleElement = (insertStyleElement_default());

var tree_style_base_update = injectStylesIntoStyleTag_default()(tree_style_base/* default */.Z, tree_style_base_options);




       /* harmony default export */ const packages_tree_style_base = (tree_style_base/* default */.Z && tree_style_base/* default */.Z.locals ? tree_style_base/* default */.Z.locals : undefined);

;// CONCATENATED MODULE: ../packages/tree/style/index.js




// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/css-loader/dist/cjs.js!../packages/tree-extension/style/base.css
var tree_extension_style_base = __webpack_require__(94623);
;// CONCATENATED MODULE: ../packages/tree-extension/style/base.css

      
      
      
      
      
      
      
      
      

var tree_extension_style_base_options = {};

tree_extension_style_base_options.styleTagTransform = (styleTagTransform_default());
tree_extension_style_base_options.setAttributes = (setAttributesWithoutAttributes_default());

      tree_extension_style_base_options.insert = insertBySelector_default().bind(null, "head");
    
tree_extension_style_base_options.domAPI = (styleDomAPI_default());
tree_extension_style_base_options.insertStyleElement = (insertStyleElement_default());

var tree_extension_style_base_update = injectStylesIntoStyleTag_default()(tree_extension_style_base/* default */.Z, tree_extension_style_base_options);




       /* harmony default export */ const packages_tree_extension_style_base = (tree_extension_style_base/* default */.Z && tree_extension_style_base/* default */.Z.locals ? tree_extension_style_base/* default */.Z.locals : undefined);

;// CONCATENATED MODULE: ../packages/tree-extension/style/index.js






// EXTERNAL MODULE: ../node_modules/@jupyterlab/builder/node_modules/css-loader/dist/cjs.js!../packages/ui-components/style/base.css
var ui_components_style_base = __webpack_require__(99402);
;// CONCATENATED MODULE: ../packages/ui-components/style/base.css

      
      
      
      
      
      
      
      
      

var ui_components_style_base_options = {};

ui_components_style_base_options.styleTagTransform = (styleTagTransform_default());
ui_components_style_base_options.setAttributes = (setAttributesWithoutAttributes_default());

      ui_components_style_base_options.insert = insertBySelector_default().bind(null, "head");
    
ui_components_style_base_options.domAPI = (styleDomAPI_default());
ui_components_style_base_options.insertStyleElement = (insertStyleElement_default());

var ui_components_style_base_update = injectStylesIntoStyleTag_default()(ui_components_style_base/* default */.Z, ui_components_style_base_options);




       /* harmony default export */ const packages_ui_components_style_base = (ui_components_style_base/* default */.Z && ui_components_style_base/* default */.Z.locals ? ui_components_style_base/* default */.Z.locals : undefined);

;// CONCATENATED MODULE: ../packages/ui-components/style/index.js
/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/



// EXTERNAL MODULE: ../node_modules/@jupyterlab/application-extension/style/index.js + 1 modules
var application_extension_style = __webpack_require__(32723);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/apputils-extension/style/index.js + 1 modules
var apputils_extension_style = __webpack_require__(82649);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/cell-toolbar-extension/style/index.js + 2 modules
var cell_toolbar_extension_style = __webpack_require__(15491);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/celltags-extension/style/index.js + 1 modules
var celltags_extension_style = __webpack_require__(92645);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/codemirror-extension/style/index.js
var codemirror_extension_style = __webpack_require__(18934);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/completer-extension/style/index.js
var completer_extension_style = __webpack_require__(24017);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/console-extension/style/index.js
var console_extension_style = __webpack_require__(72867);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/csvviewer-extension/style/index.js + 2 modules
var csvviewer_extension_style = __webpack_require__(93751);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/debugger-extension/style/index.js + 2 modules
var debugger_extension_style = __webpack_require__(75617);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/docmanager-extension/style/index.js
var docmanager_extension_style = __webpack_require__(91532);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/documentsearch-extension/style/index.js
var documentsearch_extension_style = __webpack_require__(20603);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/extensionmanager-extension/style/index.js + 2 modules
var extensionmanager_extension_style = __webpack_require__(17066);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/filebrowser-extension/style/index.js + 1 modules
var filebrowser_extension_style = __webpack_require__(39380);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/fileeditor-extension/style/index.js
var fileeditor_extension_style = __webpack_require__(9755);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/help-extension/style/index.js + 1 modules
var help_extension_style = __webpack_require__(82164);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/htmlviewer-extension/style/index.js + 2 modules
var htmlviewer_extension_style = __webpack_require__(94938);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/hub-extension/style/index.js
var hub_extension_style = __webpack_require__(53555);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/imageviewer-extension/style/index.js + 2 modules
var imageviewer_extension_style = __webpack_require__(86258);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/javascript-extension/style/index.js + 1 modules
var javascript_extension_style = __webpack_require__(11575);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/json-extension/style/index.js + 1 modules
var json_extension_style = __webpack_require__(10887);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/lsp-extension/style/index.js + 1 modules
var lsp_extension_style = __webpack_require__(26449);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/mainmenu-extension/style/index.js
var mainmenu_extension_style = __webpack_require__(3727);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/markdownviewer-extension/style/index.js + 2 modules
var markdownviewer_extension_style = __webpack_require__(97635);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/markedparser-extension/style/index.js + 1 modules
var markedparser_extension_style = __webpack_require__(16856);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/mathjax-extension/style/index.js + 1 modules
var mathjax_extension_style = __webpack_require__(46165);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/mermaid-extension/style/index.js + 1 modules
var mermaid_extension_style = __webpack_require__(8604);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/metadataform-extension/style/index.js
var metadataform_extension_style = __webpack_require__(26053);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/notebook-extension/style/index.js
var notebook_extension_style = __webpack_require__(84221);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/pdf-extension/style/index.js + 1 modules
var pdf_extension_style = __webpack_require__(53927);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/pluginmanager-extension/style/index.js
var pluginmanager_extension_style = __webpack_require__(67074);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/running-extension/style/index.js
var running_extension_style = __webpack_require__(82799);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/settingeditor-extension/style/index.js + 4 modules
var settingeditor_extension_style = __webpack_require__(58068);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/shortcuts-extension/style/index.js + 1 modules
var shortcuts_extension_style = __webpack_require__(41346);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/terminal-extension/style/index.js + 2 modules
var terminal_extension_style = __webpack_require__(47639);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/toc-extension/style/index.js + 1 modules
var toc_extension_style = __webpack_require__(65315);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/tooltip-extension/style/index.js + 2 modules
var tooltip_extension_style = __webpack_require__(87635);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/translation-extension/style/index.js
var translation_extension_style = __webpack_require__(37609);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/ui-components-extension/style/index.js
var ui_components_extension_style = __webpack_require__(49733);
// EXTERNAL MODULE: ../node_modules/@jupyterlab/vega5-extension/style/index.js + 1 modules
var vega5_extension_style = __webpack_require__(40745);
;// CONCATENATED MODULE: ./build/style.js
/* This is a generated file of CSS imports */
/* It was generated by @jupyterlab/builder in Build.ensureAssets() */




















































/***/ }),

/***/ 58027:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(82971);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11368);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__);
// Imports


var ___CSS_LOADER_EXPORT___ = _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default()((_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default()));
// Module
___CSS_LOADER_EXPORT___.push([module.id, `/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-NotebookSpacer {
  flex-grow: 1;
  flex-shrink: 1;
}

.jp-MainAreaWidget {
  height: 100%;
}

.jp-Toolbar > .jp-Toolbar-item {
  height: unset;
}

#jp-UserMenu {
  flex: 0 0 auto;
  display: flex;
  text-align: center;
  margin-top: 8px;
}

.jp-MimeDocument .jp-RenderedJSON {
  background: var(--jp-layout-color0);
}

/* Hide the stub toolbar that appears above terminals and documents */

.jp-MainAreaWidget > .jp-Toolbar-micro {
  display: none;
}

#jp-NotebookLogo {
  /* bring logo to the front so it is selectable by tab*/
  z-index: 10;
}

/* Hide the notification status item */
.jp-Notification-Status {
  display: none;
}
`, "",{"version":3,"sources":["webpack://./../packages/application-extension/style/base.css"],"names":[],"mappings":"AAAA;;;;8EAI8E;;AAE9E;EACE,YAAY;EACZ,cAAc;AAChB;;AAEA;EACE,YAAY;AACd;;AAEA;EACE,aAAa;AACf;;AAEA;EACE,cAAc;EACd,aAAa;EACb,kBAAkB;EAClB,eAAe;AACjB;;AAEA;EACE,mCAAmC;AACrC;;AAEA,qEAAqE;;AAErE;EACE,aAAa;AACf;;AAEA;EACE,sDAAsD;EACtD,WAAW;AACb;;AAEA,sCAAsC;AACtC;EACE,aAAa;AACf","sourcesContent":["/*-----------------------------------------------------------------------------\n| Copyright (c) Jupyter Development Team.\n|\n| Distributed under the terms of the Modified BSD License.\n|----------------------------------------------------------------------------*/\n\n.jp-NotebookSpacer {\n  flex-grow: 1;\n  flex-shrink: 1;\n}\n\n.jp-MainAreaWidget {\n  height: 100%;\n}\n\n.jp-Toolbar > .jp-Toolbar-item {\n  height: unset;\n}\n\n#jp-UserMenu {\n  flex: 0 0 auto;\n  display: flex;\n  text-align: center;\n  margin-top: 8px;\n}\n\n.jp-MimeDocument .jp-RenderedJSON {\n  background: var(--jp-layout-color0);\n}\n\n/* Hide the stub toolbar that appears above terminals and documents */\n\n.jp-MainAreaWidget > .jp-Toolbar-micro {\n  display: none;\n}\n\n#jp-NotebookLogo {\n  /* bring logo to the front so it is selectable by tab*/\n  z-index: 10;\n}\n\n/* Hide the notification status item */\n.jp-Notification-Status {\n  display: none;\n}\n"],"sourceRoot":""}]);
// Exports
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (___CSS_LOADER_EXPORT___);


/***/ }),

/***/ 26633:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(82971);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11368);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__);
// Imports


var ___CSS_LOADER_EXPORT___ = _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default()((_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default()));
// Module
___CSS_LOADER_EXPORT___.push([module.id, `/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

:root {
  --jp-private-topbar-height: 28px;
  /* Override the layout-2 color for the dark theme */
  --md-grey-800: #323232;
  --jp-notebook-max-width: 1200px;
}

/*
  Override the default background
  See https://github.com/jupyterlab/jupyterlab/pull/16519 for more information
*/
body.jp-ThemedContainer {
  margin: 0;
  padding: 0;
  background: var(--jp-layout-color2);
}

#main.jp-ThemedContainer {
  position: absolute;
  top: 0;
  left: 0;
  right: 0;
  bottom: 0;
  background: var(--jp-layout-color2);
}

#top-panel-wrapper {
  min-height: calc(1.5 * var(--jp-private-topbar-height));
  border-bottom: var(--jp-border-width) solid var(--jp-border-color0);
  background: var(--jp-layout-color1);
}

#top-panel {
  display: flex;
  min-height: calc(1.5 * var(--jp-private-topbar-height));
  padding-left: 5px;
  padding-right: 5px;
  margin-left: auto;
  margin-right: auto;
  max-width: 1200px;
}

#menu-panel-wrapper {
  min-height: var(--jp-private-topbar-height);
  background: var(--jp-layout-color1);
  border-bottom: var(--jp-border-width) solid var(--jp-border-color0);
  box-shadow: var(--jp-elevation-z1);
}

#menu-panel {
  display: flex;
  min-height: var(--jp-private-topbar-height);
  background: var(--jp-layout-color1);
  padding-left: 5px;
  padding-right: 5px;
  margin-left: auto;
  margin-right: auto;
  max-width: var(--jp-notebook-max-width);
}

#main-panel {
  margin-left: auto;
  margin-right: auto;
  max-width: var(--jp-notebook-max-width);
}

#spacer-widget-top {
  min-height: 16px;
}

/* Only edit pages should have a bottom space */

body[data-notebook='edit'] #spacer-widget-bottom {
  min-height: 16px;
}

/* Special case notebooks as document oriented pages */

[data-notebook]:not(body[data-notebook='notebooks']) #main-panel {
  box-shadow: var(--jp-elevation-z4);
}

.jp-TreePanel > .lm-TabPanel-stackedPanel {
  box-shadow: var(--jp-elevation-z4);
}

body[data-notebook='notebooks'] #main-panel {
  margin-left: unset;
  margin-right: unset;
  max-width: unset;
}

body[data-notebook='notebooks'] #spacer-widget-top {
  min-height: unset;
}

#main-panel > .jp-TreePanel {
  padding: 0px 5px;
}

@media only screen and (max-width: 760px) {
  #main-panel > .jp-TreePanel {
    margin: 0px -5px;
  }
}
`, "",{"version":3,"sources":["webpack://./../packages/application/style/base.css"],"names":[],"mappings":"AAAA;;;8EAG8E;;AAE9E;EACE,gCAAgC;EAChC,mDAAmD;EACnD,sBAAsB;EACtB,+BAA+B;AACjC;;AAEA;;;CAGC;AACD;EACE,SAAS;EACT,UAAU;EACV,mCAAmC;AACrC;;AAEA;EACE,kBAAkB;EAClB,MAAM;EACN,OAAO;EACP,QAAQ;EACR,SAAS;EACT,mCAAmC;AACrC;;AAEA;EACE,uDAAuD;EACvD,mEAAmE;EACnE,mCAAmC;AACrC;;AAEA;EACE,aAAa;EACb,uDAAuD;EACvD,iBAAiB;EACjB,kBAAkB;EAClB,iBAAiB;EACjB,kBAAkB;EAClB,iBAAiB;AACnB;;AAEA;EACE,2CAA2C;EAC3C,mCAAmC;EACnC,mEAAmE;EACnE,kCAAkC;AACpC;;AAEA;EACE,aAAa;EACb,2CAA2C;EAC3C,mCAAmC;EACnC,iBAAiB;EACjB,kBAAkB;EAClB,iBAAiB;EACjB,kBAAkB;EAClB,uCAAuC;AACzC;;AAEA;EACE,iBAAiB;EACjB,kBAAkB;EAClB,uCAAuC;AACzC;;AAEA;EACE,gBAAgB;AAClB;;AAEA,+CAA+C;;AAE/C;EACE,gBAAgB;AAClB;;AAEA,sDAAsD;;AAEtD;EACE,kCAAkC;AACpC;;AAEA;EACE,kCAAkC;AACpC;;AAEA;EACE,kBAAkB;EAClB,mBAAmB;EACnB,gBAAgB;AAClB;;AAEA;EACE,iBAAiB;AACnB;;AAEA;EACE,gBAAgB;AAClB;;AAEA;EACE;IACE,gBAAgB;EAClB;AACF","sourcesContent":["/*-----------------------------------------------------------------------------\n| Copyright (c) Jupyter Development Team.\n| Distributed under the terms of the Modified BSD License.\n|----------------------------------------------------------------------------*/\n\n:root {\n  --jp-private-topbar-height: 28px;\n  /* Override the layout-2 color for the dark theme */\n  --md-grey-800: #323232;\n  --jp-notebook-max-width: 1200px;\n}\n\n/*\n  Override the default background\n  See https://github.com/jupyterlab/jupyterlab/pull/16519 for more information\n*/\nbody.jp-ThemedContainer {\n  margin: 0;\n  padding: 0;\n  background: var(--jp-layout-color2);\n}\n\n#main.jp-ThemedContainer {\n  position: absolute;\n  top: 0;\n  left: 0;\n  right: 0;\n  bottom: 0;\n  background: var(--jp-layout-color2);\n}\n\n#top-panel-wrapper {\n  min-height: calc(1.5 * var(--jp-private-topbar-height));\n  border-bottom: var(--jp-border-width) solid var(--jp-border-color0);\n  background: var(--jp-layout-color1);\n}\n\n#top-panel {\n  display: flex;\n  min-height: calc(1.5 * var(--jp-private-topbar-height));\n  padding-left: 5px;\n  padding-right: 5px;\n  margin-left: auto;\n  margin-right: auto;\n  max-width: 1200px;\n}\n\n#menu-panel-wrapper {\n  min-height: var(--jp-private-topbar-height);\n  background: var(--jp-layout-color1);\n  border-bottom: var(--jp-border-width) solid var(--jp-border-color0);\n  box-shadow: var(--jp-elevation-z1);\n}\n\n#menu-panel {\n  display: flex;\n  min-height: var(--jp-private-topbar-height);\n  background: var(--jp-layout-color1);\n  padding-left: 5px;\n  padding-right: 5px;\n  margin-left: auto;\n  margin-right: auto;\n  max-width: var(--jp-notebook-max-width);\n}\n\n#main-panel {\n  margin-left: auto;\n  margin-right: auto;\n  max-width: var(--jp-notebook-max-width);\n}\n\n#spacer-widget-top {\n  min-height: 16px;\n}\n\n/* Only edit pages should have a bottom space */\n\nbody[data-notebook='edit'] #spacer-widget-bottom {\n  min-height: 16px;\n}\n\n/* Special case notebooks as document oriented pages */\n\n[data-notebook]:not(body[data-notebook='notebooks']) #main-panel {\n  box-shadow: var(--jp-elevation-z4);\n}\n\n.jp-TreePanel > .lm-TabPanel-stackedPanel {\n  box-shadow: var(--jp-elevation-z4);\n}\n\nbody[data-notebook='notebooks'] #main-panel {\n  margin-left: unset;\n  margin-right: unset;\n  max-width: unset;\n}\n\nbody[data-notebook='notebooks'] #spacer-widget-top {\n  min-height: unset;\n}\n\n#main-panel > .jp-TreePanel {\n  padding: 0px 5px;\n}\n\n@media only screen and (max-width: 760px) {\n  #main-panel > .jp-TreePanel {\n    margin: 0px -5px;\n  }\n}\n"],"sourceRoot":""}]);
// Exports
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (___CSS_LOADER_EXPORT___);


/***/ }),

/***/ 93862:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(82971);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11368);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__);
// Imports


var ___CSS_LOADER_EXPORT___ = _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default()((_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default()));
// Module
___CSS_LOADER_EXPORT___.push([module.id, `/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|
| Adapted from JupyterLab's packages/application/style/sidepanel.css.
|----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
| Variables
|----------------------------------------------------------------------------*/

:root {
  --jp-private-sidebar-tab-width: 32px;
}

/*-----------------------------------------------------------------------------
| SideBar
|----------------------------------------------------------------------------*/

/* Stack panels */

#jp-right-stack,
#jp-left-stack {
  display: flex;
  flex-direction: column;
  min-width: var(--jp-sidebar-min-width);
}

#jp-left-stack .jp-SidePanel-collapse,
#jp-right-stack .jp-SidePanel-collapse {
  display: flex;
  flex: 0 0 auto;
  min-height: 0;
  padding: 0;
}

#jp-left-stack .jp-SidePanel-collapse {
  justify-content: right;
}

#jp-right-stack .jp-SidePanel-collapse {
  justify-content: left;
}

#jp-left-stack .lm-StackedPanel,
#jp-right-stack .lm-StackedPanel {
  flex: 1 1 auto;
}
`, "",{"version":3,"sources":["webpack://./../packages/application/style/sidepanel.css"],"names":[],"mappings":"AAAA;;;;;8EAK8E;;AAE9E;;8EAE8E;;AAE9E;EACE,oCAAoC;AACtC;;AAEA;;8EAE8E;;AAE9E,iBAAiB;;AAEjB;;EAEE,aAAa;EACb,sBAAsB;EACtB,sCAAsC;AACxC;;AAEA;;EAEE,aAAa;EACb,cAAc;EACd,aAAa;EACb,UAAU;AACZ;;AAEA;EACE,sBAAsB;AACxB;;AAEA;EACE,qBAAqB;AACvB;;AAEA;;EAEE,cAAc;AAChB","sourcesContent":["/*-----------------------------------------------------------------------------\n| Copyright (c) Jupyter Development Team.\n| Distributed under the terms of the Modified BSD License.\n|\n| Adapted from JupyterLab's packages/application/style/sidepanel.css.\n|----------------------------------------------------------------------------*/\n\n/*-----------------------------------------------------------------------------\n| Variables\n|----------------------------------------------------------------------------*/\n\n:root {\n  --jp-private-sidebar-tab-width: 32px;\n}\n\n/*-----------------------------------------------------------------------------\n| SideBar\n|----------------------------------------------------------------------------*/\n\n/* Stack panels */\n\n#jp-right-stack,\n#jp-left-stack {\n  display: flex;\n  flex-direction: column;\n  min-width: var(--jp-sidebar-min-width);\n}\n\n#jp-left-stack .jp-SidePanel-collapse,\n#jp-right-stack .jp-SidePanel-collapse {\n  display: flex;\n  flex: 0 0 auto;\n  min-height: 0;\n  padding: 0;\n}\n\n#jp-left-stack .jp-SidePanel-collapse {\n  justify-content: right;\n}\n\n#jp-right-stack .jp-SidePanel-collapse {\n  justify-content: left;\n}\n\n#jp-left-stack .lm-StackedPanel,\n#jp-right-stack .lm-StackedPanel {\n  flex: 1 1 auto;\n}\n"],"sourceRoot":""}]);
// Exports
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (___CSS_LOADER_EXPORT___);


/***/ }),

/***/ 99776:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(82971);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11368);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__);
// Imports


var ___CSS_LOADER_EXPORT___ = _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default()((_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default()));
// Module
___CSS_LOADER_EXPORT___.push([module.id, ``, "",{"version":3,"sources":[],"names":[],"mappings":"","sourceRoot":""}]);
// Exports
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (___CSS_LOADER_EXPORT___);


/***/ }),

/***/ 64399:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(82971);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11368);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__);
// Imports


var ___CSS_LOADER_EXPORT___ = _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default()((_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default()));
// Module
___CSS_LOADER_EXPORT___.push([module.id, `.jp-Document {
  height: 100%;
}
`, "",{"version":3,"sources":["webpack://./../packages/docmanager-extension/style/base.css"],"names":[],"mappings":"AAAA;EACE,YAAY;AACd","sourcesContent":[".jp-Document {\n  height: 100%;\n}\n"],"sourceRoot":""}]);
// Exports
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (___CSS_LOADER_EXPORT___);


/***/ }),

/***/ 11112:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(82971);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11368);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__);
// Imports


var ___CSS_LOADER_EXPORT___ = _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default()((_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default()));
// Module
___CSS_LOADER_EXPORT___.push([module.id, ``, "",{"version":3,"sources":[],"names":[],"mappings":"","sourceRoot":""}]);
// Exports
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (___CSS_LOADER_EXPORT___);


/***/ }),

/***/ 49482:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(82971);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11368);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__);
// Imports


var ___CSS_LOADER_EXPORT___ = _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default()((_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default()));
// Module
___CSS_LOADER_EXPORT___.push([module.id, `.jp-AboutNotebook .jp-Dialog-header {
  justify-content: center;
  padding: 0;
}

.jp-AboutNotebook-header {
  display: flex;
  flex-direction: row;
  align-items: center;
  padding: var(--jp-flat-button-padding);
}

.jp-AboutNotebook-header-text {
  margin-left: 16px;
}

.jp-AboutNotebook-version {
  color: var(--jp-ui-font-color1);
  font-size: var(--jp-ui-font-size1);
  padding-bottom: 30px;
  font-weight: 400;
  letter-spacing: 0.4px;
  line-height: 1.12;
  min-width: 360px;
  text-align: center;
}

.jp-AboutNotebook-body {
  display: flex;
  font-size: var(--jp-ui-font-size2);
  padding: var(--jp-flat-button-padding);
  color: var(--jp-ui-font-color1);
  text-align: center;
  flex-direction: column;
  min-width: 360px;
  overflow: hidden;
}

.jp-AboutNotebook-about-body pre {
  white-space: pre-wrap;
}

.jp-AboutNotebook-about-externalLinks {
  display: flex;
  flex-direction: column;
  justify-content: flex-start;
  align-items: flex-start;
  color: var(--jp-warn-color0);
}

.jp-AboutNotebook-about-copyright {
  padding-top: 25px;
}
`, "",{"version":3,"sources":["webpack://./../packages/help-extension/style/base.css"],"names":[],"mappings":"AAAA;EACE,uBAAuB;EACvB,UAAU;AACZ;;AAEA;EACE,aAAa;EACb,mBAAmB;EACnB,mBAAmB;EACnB,sCAAsC;AACxC;;AAEA;EACE,iBAAiB;AACnB;;AAEA;EACE,+BAA+B;EAC/B,kCAAkC;EAClC,oBAAoB;EACpB,gBAAgB;EAChB,qBAAqB;EACrB,iBAAiB;EACjB,gBAAgB;EAChB,kBAAkB;AACpB;;AAEA;EACE,aAAa;EACb,kCAAkC;EAClC,sCAAsC;EACtC,+BAA+B;EAC/B,kBAAkB;EAClB,sBAAsB;EACtB,gBAAgB;EAChB,gBAAgB;AAClB;;AAEA;EACE,qBAAqB;AACvB;;AAEA;EACE,aAAa;EACb,sBAAsB;EACtB,2BAA2B;EAC3B,uBAAuB;EACvB,4BAA4B;AAC9B;;AAEA;EACE,iBAAiB;AACnB","sourcesContent":[".jp-AboutNotebook .jp-Dialog-header {\n  justify-content: center;\n  padding: 0;\n}\n\n.jp-AboutNotebook-header {\n  display: flex;\n  flex-direction: row;\n  align-items: center;\n  padding: var(--jp-flat-button-padding);\n}\n\n.jp-AboutNotebook-header-text {\n  margin-left: 16px;\n}\n\n.jp-AboutNotebook-version {\n  color: var(--jp-ui-font-color1);\n  font-size: var(--jp-ui-font-size1);\n  padding-bottom: 30px;\n  font-weight: 400;\n  letter-spacing: 0.4px;\n  line-height: 1.12;\n  min-width: 360px;\n  text-align: center;\n}\n\n.jp-AboutNotebook-body {\n  display: flex;\n  font-size: var(--jp-ui-font-size2);\n  padding: var(--jp-flat-button-padding);\n  color: var(--jp-ui-font-color1);\n  text-align: center;\n  flex-direction: column;\n  min-width: 360px;\n  overflow: hidden;\n}\n\n.jp-AboutNotebook-about-body pre {\n  white-space: pre-wrap;\n}\n\n.jp-AboutNotebook-about-externalLinks {\n  display: flex;\n  flex-direction: column;\n  justify-content: flex-start;\n  align-items: flex-start;\n  color: var(--jp-warn-color0);\n}\n\n.jp-AboutNotebook-about-copyright {\n  padding-top: 25px;\n}\n"],"sourceRoot":""}]);
// Exports
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (___CSS_LOADER_EXPORT___);


/***/ }),

/***/ 71564:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(82971);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11368);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_cjs_js_variables_css__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(99686);
// Imports



var ___CSS_LOADER_EXPORT___ = _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default()((_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default()));
___CSS_LOADER_EXPORT___.i(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_cjs_js_variables_css__WEBPACK_IMPORTED_MODULE_2__/* ["default"] */ .Z);
// Module
___CSS_LOADER_EXPORT___.push([module.id, `/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/**
  Document oriented look for the notebook.
  This includes changes to the look and feel of the JupyterLab Notebook
  component like:
  - scrollbar to the right of the page
  - drop shadow on the notebook
  - smaller empty space at the bottom of the notebook
  - compact view on mobile
*/

/* Make the notebook take up the full width of the page when jp-mod-fullwidth is set */

body[data-notebook='notebooks']
  .jp-NotebookPanel.jp-mod-fullwidth
  .jp-WindowedPanel-outer {
  padding-left: unset;
  padding-right: unset !important;
  width: unset;
}

/* Keep the notebook centered on the page */

body[data-notebook='notebooks'] .jp-NotebookPanel-toolbar {
  padding-left: calc(calc(100% - var(--jp-notebook-max-width)) * 0.5);
  padding-right: calc(calc(100% - var(--jp-notebook-max-width)) * 0.5);
}

body[data-notebook='notebooks'] .jp-WindowedPanel-outer {
  width: unset !important;
  padding-top: unset;
  padding-left: calc(calc(100% - var(--jp-notebook-max-width)) * 0.5);
  padding-right: calc(
    calc(
        100% - var(--jp-notebook-max-width) - var(--jp-notebook-padding-offset)
      ) * 0.5
  ) !important;
  background: var(--jp-layout-color2);
}

body[data-notebook='notebooks'] .jp-WindowedPanel-inner {
  margin-top: var(--jp-notebook-toolbar-margin-bottom);
  /* Adjustments for the extra top and bottom notebook padding */
  margin-bottom: calc(4 * var(--jp-notebook-padding));
}

body[data-notebook='notebooks'] .jp-Notebook-cell {
  background: var(--jp-layout-color0);
  padding-left: calc(2 * var(--jp-cell-padding));
  padding-right: calc(2 * var(--jp-cell-padding));
}

/* Empty space at the bottom of the notebook (similar to classic) */
body[data-notebook='notebooks']
  .jp-Notebook.jp-mod-scrollPastEnd
  .jp-WindowedPanel-outer::after {
  min-height: 100px;
}

/* Fix background colors */

body[data-notebook='notebooks'] .jp-WindowedPanel-outer > * {
  background: var(--jp-layout-color0);
}

body[data-notebook='notebooks']
  .jp-Notebook.jp-mod-commandMode
  .jp-Cell.jp-mod-active.jp-mod-selected:not(.jp-mod-multiSelected) {
  background: var(--jp-layout-color0) !important;
}

body[data-notebook='notebooks']
  .jp-Notebook
  .jp-Notebook-cell:not(:first-child)::before {
  content: ' ';
  height: 100%;
  position: absolute;
  top: 0;
  width: 11px;
}

/* Cell toolbar adjustments */

body[data-notebook='notebooks'] .jp-cell-toolbar {
  background: unset;
  box-shadow: unset;
}

/** first code cell on mobile
    (keep the selector above the media query)
*/
body[data-notebook='notebooks']
  .jp-CodeCell[data-windowed-list-index='0']
  .jp-cell-toolbar {
  top: unset;
}

@media only screen and (max-width: 760px) {
  /* first code cell on mobile */
  body[data-notebook='notebooks']
    .jp-CodeCell[data-windowed-list-index='0']
    .jp-cell-toolbar {
    top: var(--jp-notebook-padding);
  }

  body[data-notebook='notebooks'] .jp-MarkdownCell .jp-cell-toolbar,
  body[data-notebook='notebooks'] .jp-RawCell .jp-cell-toolbar {
    top: calc(0.5 * var(--jp-notebook-padding));
  }
}

/* Tweak the notebook footer (to add a new cell) */
body[data-notebook='notebooks'] .jp-Notebook-footer {
  background: unset;
  width: 100%;
  margin-left: unset;
}

/* Mobile View */

body[data-format='mobile'] .jp-NotebookCheckpoint {
  display: none;
}

body[data-format='mobile'] .jp-WindowedPanel-outer > *:first-child {
  margin-top: 0;
}

body[data-format='mobile'] .jp-ToolbarButton .jp-DebuggerBugButton {
  display: none;
}

body[data-notebook='notebooks'] .jp-WindowedPanel-viewport {
  background: var(--jp-layout-color0);
  box-shadow: var(--jp-elevation-z4);
  padding: unset;

  /* Extra padding at the top and bottom so the notebook looks nicer */
  padding-top: calc(2 * var(--jp-notebook-padding));
  padding-bottom: calc(2 * var(--jp-notebook-padding));
}

/* Notebook box shadow */

body[data-notebook='notebooks']
  .jp-Notebook
  > *:first-child:last-child::before {
  content: '';
  position: absolute;
  top: 0;
  bottom: 0;
  left: 0;
  right: 0;
  box-shadow: 0px 0px 12px 1px var(--jp-shadow-umbra-color);
}

/* Additional customizations of the components on the notebook page */

.jp-NotebookKernelLogo {
  flex: 0 0 auto;
  display: flex;
  align-items: center;
  text-align: center;
  margin-right: 8px;
}

.jp-NotebookKernelLogo img {
  max-width: 28px;
  max-height: 28px;
  display: flex;
}

.jp-NotebookKernelStatus {
  margin: 0;
  font-weight: normal;
  font-size: var(--jp-ui-font-size1);
  color: var(--jp-ui-font-color0);
  font-family: var(--jp-ui-font-family);
  line-height: var(--jp-private-title-panel-height);
  padding-left: var(--jp-kernel-status-padding);
  padding-right: var(--jp-kernel-status-padding);
}

.jp-NotebookKernelStatus-error {
  background-color: var(--jp-error-color0);
}

.jp-NotebookKernelStatus-warn {
  background-color: var(--jp-warn-color0);
}

.jp-NotebookKernelStatus-info {
  background-color: var(--jp-info-color0);
}

.jp-NotebookKernelStatus-fade {
  animation: 0.5s fade-out forwards;
}

.jp-NotebookTrustedStatus {
  background: var(--jp-layout-color1);
  color: var(--jp-ui-font-color1);
  margin-top: 4px;
  margin-bottom: 4px;
  border: solid 1px var(--jp-border-color2);
  cursor: help;
}

.jp-NotebookTrustedStatus-not-trusted {
  cursor: pointer;
}

@keyframes fade-out {
  0% {
    opacity: 1;
  }
  100% {
    opacity: 0;
  }
}

#jp-title h1 {
  cursor: pointer;
  font-size: 18px;
  margin: 0;
  font-weight: normal;
  color: var(--jp-ui-font-color0);
  font-family: var(--jp-ui-font-family);
  line-height: calc(1.5 * var(--jp-private-title-panel-height));
  text-overflow: ellipsis;
  overflow: hidden;
  white-space: nowrap;
}

#jp-title h1:hover {
  background: var(--jp-layout-color2);
}

.jp-NotebookCheckpoint {
  font-size: 14px;
  margin-left: 5px;
  margin-right: 5px;
  font-weight: normal;
  color: var(--jp-ui-font-color0);
  font-family: var(--jp-ui-font-family);
  line-height: calc(1.5 * var(--jp-private-title-panel-height));
  text-overflow: ellipsis;
  overflow: hidden;
  white-space: nowrap;
}

.jp-skiplink {
  position: absolute;
  top: -100em;
}

.jp-skiplink:focus-within {
  position: absolute;
  z-index: 10000;
  top: 0;
  left: 46%;
  margin: 0 auto;
  padding: 1em;
  width: 15%;
  box-shadow: var(--jp-elevation-z4);
  border-radius: 4px;
  background: var(--jp-layout-color0);
  text-align: center;
}

.jp-skiplink:focus-within a {
  text-decoration: underline;
  color: var(--jp-content-link-color);
}
`, "",{"version":3,"sources":["webpack://./../packages/notebook-extension/style/base.css"],"names":[],"mappings":"AAAA;;;;8EAI8E;;AAI9E;;;;;;;;CAQC;;AAED,sFAAsF;;AAEtF;;;EAGE,mBAAmB;EACnB,+BAA+B;EAC/B,YAAY;AACd;;AAEA,2CAA2C;;AAE3C;EACE,mEAAmE;EACnE,oEAAoE;AACtE;;AAEA;EACE,uBAAuB;EACvB,kBAAkB;EAClB,mEAAmE;EACnE;;;;cAIY;EACZ,mCAAmC;AACrC;;AAEA;EACE,oDAAoD;EACpD,8DAA8D;EAC9D,mDAAmD;AACrD;;AAEA;EACE,mCAAmC;EACnC,8CAA8C;EAC9C,+CAA+C;AACjD;;AAEA,mEAAmE;AACnE;;;EAGE,iBAAiB;AACnB;;AAEA,0BAA0B;;AAE1B;EACE,mCAAmC;AACrC;;AAEA;;;EAGE,8CAA8C;AAChD;;AAEA;;;EAGE,YAAY;EACZ,YAAY;EACZ,kBAAkB;EAClB,MAAM;EACN,WAAW;AACb;;AAEA,6BAA6B;;AAE7B;EACE,iBAAiB;EACjB,iBAAiB;AACnB;;AAEA;;CAEC;AACD;;;EAGE,UAAU;AACZ;;AAEA;EACE,8BAA8B;EAC9B;;;IAGE,+BAA+B;EACjC;;EAEA;;IAEE,2CAA2C;EAC7C;AACF;;AAEA,kDAAkD;AAClD;EACE,iBAAiB;EACjB,WAAW;EACX,kBAAkB;AACpB;;AAEA,gBAAgB;;AAEhB;EACE,aAAa;AACf;;AAEA;EACE,aAAa;AACf;;AAEA;EACE,aAAa;AACf;;AAEA;EACE,mCAAmC;EACnC,kCAAkC;EAClC,cAAc;;EAEd,oEAAoE;EACpE,iDAAiD;EACjD,oDAAoD;AACtD;;AAEA,wBAAwB;;AAExB;;;EAGE,WAAW;EACX,kBAAkB;EAClB,MAAM;EACN,SAAS;EACT,OAAO;EACP,QAAQ;EACR,yDAAyD;AAC3D;;AAEA,qEAAqE;;AAErE;EACE,cAAc;EACd,aAAa;EACb,mBAAmB;EACnB,kBAAkB;EAClB,iBAAiB;AACnB;;AAEA;EACE,eAAe;EACf,gBAAgB;EAChB,aAAa;AACf;;AAEA;EACE,SAAS;EACT,mBAAmB;EACnB,kCAAkC;EAClC,+BAA+B;EAC/B,qCAAqC;EACrC,iDAAiD;EACjD,6CAA6C;EAC7C,8CAA8C;AAChD;;AAEA;EACE,wCAAwC;AAC1C;;AAEA;EACE,uCAAuC;AACzC;;AAEA;EACE,uCAAuC;AACzC;;AAEA;EACE,iCAAiC;AACnC;;AAEA;EACE,mCAAmC;EACnC,+BAA+B;EAC/B,eAAe;EACf,kBAAkB;EAClB,yCAAyC;EACzC,YAAY;AACd;;AAEA;EACE,eAAe;AACjB;;AAEA;EACE;IACE,UAAU;EACZ;EACA;IACE,UAAU;EACZ;AACF;;AAEA;EACE,eAAe;EACf,eAAe;EACf,SAAS;EACT,mBAAmB;EACnB,+BAA+B;EAC/B,qCAAqC;EACrC,6DAA6D;EAC7D,uBAAuB;EACvB,gBAAgB;EAChB,mBAAmB;AACrB;;AAEA;EACE,mCAAmC;AACrC;;AAEA;EACE,eAAe;EACf,gBAAgB;EAChB,iBAAiB;EACjB,mBAAmB;EACnB,+BAA+B;EAC/B,qCAAqC;EACrC,6DAA6D;EAC7D,uBAAuB;EACvB,gBAAgB;EAChB,mBAAmB;AACrB;;AAEA;EACE,kBAAkB;EAClB,WAAW;AACb;;AAEA;EACE,kBAAkB;EAClB,cAAc;EACd,MAAM;EACN,SAAS;EACT,cAAc;EACd,YAAY;EACZ,UAAU;EACV,kCAAkC;EAClC,kBAAkB;EAClB,mCAAmC;EACnC,kBAAkB;AACpB;;AAEA;EACE,0BAA0B;EAC1B,mCAAmC;AACrC","sourcesContent":["/*-----------------------------------------------------------------------------\n| Copyright (c) Jupyter Development Team.\n|\n| Distributed under the terms of the Modified BSD License.\n|----------------------------------------------------------------------------*/\n\n@import './variables.css';\n\n/**\n  Document oriented look for the notebook.\n  This includes changes to the look and feel of the JupyterLab Notebook\n  component like:\n  - scrollbar to the right of the page\n  - drop shadow on the notebook\n  - smaller empty space at the bottom of the notebook\n  - compact view on mobile\n*/\n\n/* Make the notebook take up the full width of the page when jp-mod-fullwidth is set */\n\nbody[data-notebook='notebooks']\n  .jp-NotebookPanel.jp-mod-fullwidth\n  .jp-WindowedPanel-outer {\n  padding-left: unset;\n  padding-right: unset !important;\n  width: unset;\n}\n\n/* Keep the notebook centered on the page */\n\nbody[data-notebook='notebooks'] .jp-NotebookPanel-toolbar {\n  padding-left: calc(calc(100% - var(--jp-notebook-max-width)) * 0.5);\n  padding-right: calc(calc(100% - var(--jp-notebook-max-width)) * 0.5);\n}\n\nbody[data-notebook='notebooks'] .jp-WindowedPanel-outer {\n  width: unset !important;\n  padding-top: unset;\n  padding-left: calc(calc(100% - var(--jp-notebook-max-width)) * 0.5);\n  padding-right: calc(\n    calc(\n        100% - var(--jp-notebook-max-width) - var(--jp-notebook-padding-offset)\n      ) * 0.5\n  ) !important;\n  background: var(--jp-layout-color2);\n}\n\nbody[data-notebook='notebooks'] .jp-WindowedPanel-inner {\n  margin-top: var(--jp-notebook-toolbar-margin-bottom);\n  /* Adjustments for the extra top and bottom notebook padding */\n  margin-bottom: calc(4 * var(--jp-notebook-padding));\n}\n\nbody[data-notebook='notebooks'] .jp-Notebook-cell {\n  background: var(--jp-layout-color0);\n  padding-left: calc(2 * var(--jp-cell-padding));\n  padding-right: calc(2 * var(--jp-cell-padding));\n}\n\n/* Empty space at the bottom of the notebook (similar to classic) */\nbody[data-notebook='notebooks']\n  .jp-Notebook.jp-mod-scrollPastEnd\n  .jp-WindowedPanel-outer::after {\n  min-height: 100px;\n}\n\n/* Fix background colors */\n\nbody[data-notebook='notebooks'] .jp-WindowedPanel-outer > * {\n  background: var(--jp-layout-color0);\n}\n\nbody[data-notebook='notebooks']\n  .jp-Notebook.jp-mod-commandMode\n  .jp-Cell.jp-mod-active.jp-mod-selected:not(.jp-mod-multiSelected) {\n  background: var(--jp-layout-color0) !important;\n}\n\nbody[data-notebook='notebooks']\n  .jp-Notebook\n  .jp-Notebook-cell:not(:first-child)::before {\n  content: ' ';\n  height: 100%;\n  position: absolute;\n  top: 0;\n  width: 11px;\n}\n\n/* Cell toolbar adjustments */\n\nbody[data-notebook='notebooks'] .jp-cell-toolbar {\n  background: unset;\n  box-shadow: unset;\n}\n\n/** first code cell on mobile\n    (keep the selector above the media query)\n*/\nbody[data-notebook='notebooks']\n  .jp-CodeCell[data-windowed-list-index='0']\n  .jp-cell-toolbar {\n  top: unset;\n}\n\n@media only screen and (max-width: 760px) {\n  /* first code cell on mobile */\n  body[data-notebook='notebooks']\n    .jp-CodeCell[data-windowed-list-index='0']\n    .jp-cell-toolbar {\n    top: var(--jp-notebook-padding);\n  }\n\n  body[data-notebook='notebooks'] .jp-MarkdownCell .jp-cell-toolbar,\n  body[data-notebook='notebooks'] .jp-RawCell .jp-cell-toolbar {\n    top: calc(0.5 * var(--jp-notebook-padding));\n  }\n}\n\n/* Tweak the notebook footer (to add a new cell) */\nbody[data-notebook='notebooks'] .jp-Notebook-footer {\n  background: unset;\n  width: 100%;\n  margin-left: unset;\n}\n\n/* Mobile View */\n\nbody[data-format='mobile'] .jp-NotebookCheckpoint {\n  display: none;\n}\n\nbody[data-format='mobile'] .jp-WindowedPanel-outer > *:first-child {\n  margin-top: 0;\n}\n\nbody[data-format='mobile'] .jp-ToolbarButton .jp-DebuggerBugButton {\n  display: none;\n}\n\nbody[data-notebook='notebooks'] .jp-WindowedPanel-viewport {\n  background: var(--jp-layout-color0);\n  box-shadow: var(--jp-elevation-z4);\n  padding: unset;\n\n  /* Extra padding at the top and bottom so the notebook looks nicer */\n  padding-top: calc(2 * var(--jp-notebook-padding));\n  padding-bottom: calc(2 * var(--jp-notebook-padding));\n}\n\n/* Notebook box shadow */\n\nbody[data-notebook='notebooks']\n  .jp-Notebook\n  > *:first-child:last-child::before {\n  content: '';\n  position: absolute;\n  top: 0;\n  bottom: 0;\n  left: 0;\n  right: 0;\n  box-shadow: 0px 0px 12px 1px var(--jp-shadow-umbra-color);\n}\n\n/* Additional customizations of the components on the notebook page */\n\n.jp-NotebookKernelLogo {\n  flex: 0 0 auto;\n  display: flex;\n  align-items: center;\n  text-align: center;\n  margin-right: 8px;\n}\n\n.jp-NotebookKernelLogo img {\n  max-width: 28px;\n  max-height: 28px;\n  display: flex;\n}\n\n.jp-NotebookKernelStatus {\n  margin: 0;\n  font-weight: normal;\n  font-size: var(--jp-ui-font-size1);\n  color: var(--jp-ui-font-color0);\n  font-family: var(--jp-ui-font-family);\n  line-height: var(--jp-private-title-panel-height);\n  padding-left: var(--jp-kernel-status-padding);\n  padding-right: var(--jp-kernel-status-padding);\n}\n\n.jp-NotebookKernelStatus-error {\n  background-color: var(--jp-error-color0);\n}\n\n.jp-NotebookKernelStatus-warn {\n  background-color: var(--jp-warn-color0);\n}\n\n.jp-NotebookKernelStatus-info {\n  background-color: var(--jp-info-color0);\n}\n\n.jp-NotebookKernelStatus-fade {\n  animation: 0.5s fade-out forwards;\n}\n\n.jp-NotebookTrustedStatus {\n  background: var(--jp-layout-color1);\n  color: var(--jp-ui-font-color1);\n  margin-top: 4px;\n  margin-bottom: 4px;\n  border: solid 1px var(--jp-border-color2);\n  cursor: help;\n}\n\n.jp-NotebookTrustedStatus-not-trusted {\n  cursor: pointer;\n}\n\n@keyframes fade-out {\n  0% {\n    opacity: 1;\n  }\n  100% {\n    opacity: 0;\n  }\n}\n\n#jp-title h1 {\n  cursor: pointer;\n  font-size: 18px;\n  margin: 0;\n  font-weight: normal;\n  color: var(--jp-ui-font-color0);\n  font-family: var(--jp-ui-font-family);\n  line-height: calc(1.5 * var(--jp-private-title-panel-height));\n  text-overflow: ellipsis;\n  overflow: hidden;\n  white-space: nowrap;\n}\n\n#jp-title h1:hover {\n  background: var(--jp-layout-color2);\n}\n\n.jp-NotebookCheckpoint {\n  font-size: 14px;\n  margin-left: 5px;\n  margin-right: 5px;\n  font-weight: normal;\n  color: var(--jp-ui-font-color0);\n  font-family: var(--jp-ui-font-family);\n  line-height: calc(1.5 * var(--jp-private-title-panel-height));\n  text-overflow: ellipsis;\n  overflow: hidden;\n  white-space: nowrap;\n}\n\n.jp-skiplink {\n  position: absolute;\n  top: -100em;\n}\n\n.jp-skiplink:focus-within {\n  position: absolute;\n  z-index: 10000;\n  top: 0;\n  left: 46%;\n  margin: 0 auto;\n  padding: 1em;\n  width: 15%;\n  box-shadow: var(--jp-elevation-z4);\n  border-radius: 4px;\n  background: var(--jp-layout-color0);\n  text-align: center;\n}\n\n.jp-skiplink:focus-within a {\n  text-decoration: underline;\n  color: var(--jp-content-link-color);\n}\n"],"sourceRoot":""}]);
// Exports
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (___CSS_LOADER_EXPORT___);


/***/ }),

/***/ 99686:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(82971);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11368);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__);
// Imports


var ___CSS_LOADER_EXPORT___ = _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default()((_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default()));
// Module
___CSS_LOADER_EXPORT___.push([module.id, `:root {
  --jp-notebook-toolbar-margin-bottom: 20px;
  --jp-notebook-padding-offset: 20px;

  --jp-kernel-status-padding: 5px;
}
`, "",{"version":3,"sources":["webpack://./../packages/notebook-extension/style/variables.css"],"names":[],"mappings":"AAAA;EACE,yCAAyC;EACzC,kCAAkC;;EAElC,+BAA+B;AACjC","sourcesContent":[":root {\n  --jp-notebook-toolbar-margin-bottom: 20px;\n  --jp-notebook-padding-offset: 20px;\n\n  --jp-kernel-status-padding: 5px;\n}\n"],"sourceRoot":""}]);
// Exports
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (___CSS_LOADER_EXPORT___);


/***/ }),

/***/ 8173:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(82971);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11368);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__);
// Imports


var ___CSS_LOADER_EXPORT___ = _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default()((_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default()));
// Module
___CSS_LOADER_EXPORT___.push([module.id, ``, "",{"version":3,"sources":[],"names":[],"mappings":"","sourceRoot":""}]);
// Exports
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (___CSS_LOADER_EXPORT___);


/***/ }),

/***/ 94623:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(82971);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11368);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__);
// Imports


var ___CSS_LOADER_EXPORT___ = _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default()((_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default()));
// Module
___CSS_LOADER_EXPORT___.push([module.id, `/*-----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
|
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

.jp-FileBrowser-toolbar .jp-Toolbar-item.jp-DropdownMenu,
.jp-FileBrowser-toolbar .jp-Toolbar-item.jp-ToolbarButton,
.jp-FileBrowser-toolbar .jp-Toolbar-item.jp-CommandToolbarButton {
  border: solid 1px var(--jp-border-color2);
  margin: 1px;
  padding: 0px;
}

.jp-FileBrowser-toolbar > .jp-Toolbar-item.jp-ToolbarButton:hover,
.jp-FileBrowser-toolbar > .jp-Toolbar-item.jp-CommandToolbarButton:hover,
.jp-FileBrowser-toolbar > .jp-Toolbar-item.jp-DropdownMenu:hover {
  background: var(--neutral-fill-stealth-hover);
}

.jp-FileBrowser-toolbar .lm-MenuBar-item {
  height: var(--jp-private-toolbar-height);
  display: inline-flex;
  align-items: center;
}

.jp-FileBrowser-toolbar .jp-ToolbarButtonComponent {
  height: var(--jp-flat-button-height);
}

.jp-FileBrowser-toolbar jp-button.jp-ToolbarButtonComponent:hover {
  background: inherit;
}

.jp-DirListing-content .jp-DirListing-checkboxWrapper {
  visibility: visible;
}

/* Action buttons */

.jp-FileBrowser-toolbar > .jp-FileAction > .jp-ToolbarButtonComponent > svg {
  display: none;
}

.jp-FileBrowser-toolbar > #fileAction-delete {
  background-color: var(--jp-error-color1);
}

.jp-FileBrowser-toolbar
  .jp-ToolbarButtonComponent[data-command='filebrowser:delete']
  .jp-ToolbarButtonComponent-label {
  color: var(--jp-ui-inverse-font-color1);
}

.jp-FileBrowser-toolbar .jp-FileAction {
  border: solid 1px var(--jp-border-color2);
  margin: 1px;
  min-height: var(--jp-private-toolbar-height);
}

body[data-format='mobile'] #fileAction-placeholder {
  display: none;
}
`, "",{"version":3,"sources":["webpack://./../packages/tree-extension/style/base.css"],"names":[],"mappings":"AAAA;;;;8EAI8E;;AAE9E;;;EAGE,yCAAyC;EACzC,WAAW;EACX,YAAY;AACd;;AAEA;;;EAGE,6CAA6C;AAC/C;;AAEA;EACE,wCAAwC;EACxC,oBAAoB;EACpB,mBAAmB;AACrB;;AAEA;EACE,oCAAoC;AACtC;;AAEA;EACE,mBAAmB;AACrB;;AAEA;EACE,mBAAmB;AACrB;;AAEA,mBAAmB;;AAEnB;EACE,aAAa;AACf;;AAEA;EACE,wCAAwC;AAC1C;;AAEA;;;EAGE,uCAAuC;AACzC;;AAEA;EACE,yCAAyC;EACzC,WAAW;EACX,4CAA4C;AAC9C;;AAEA;EACE,aAAa;AACf","sourcesContent":["/*-----------------------------------------------------------------------------\n| Copyright (c) Jupyter Development Team.\n|\n| Distributed under the terms of the Modified BSD License.\n|----------------------------------------------------------------------------*/\n\n.jp-FileBrowser-toolbar .jp-Toolbar-item.jp-DropdownMenu,\n.jp-FileBrowser-toolbar .jp-Toolbar-item.jp-ToolbarButton,\n.jp-FileBrowser-toolbar .jp-Toolbar-item.jp-CommandToolbarButton {\n  border: solid 1px var(--jp-border-color2);\n  margin: 1px;\n  padding: 0px;\n}\n\n.jp-FileBrowser-toolbar > .jp-Toolbar-item.jp-ToolbarButton:hover,\n.jp-FileBrowser-toolbar > .jp-Toolbar-item.jp-CommandToolbarButton:hover,\n.jp-FileBrowser-toolbar > .jp-Toolbar-item.jp-DropdownMenu:hover {\n  background: var(--neutral-fill-stealth-hover);\n}\n\n.jp-FileBrowser-toolbar .lm-MenuBar-item {\n  height: var(--jp-private-toolbar-height);\n  display: inline-flex;\n  align-items: center;\n}\n\n.jp-FileBrowser-toolbar .jp-ToolbarButtonComponent {\n  height: var(--jp-flat-button-height);\n}\n\n.jp-FileBrowser-toolbar jp-button.jp-ToolbarButtonComponent:hover {\n  background: inherit;\n}\n\n.jp-DirListing-content .jp-DirListing-checkboxWrapper {\n  visibility: visible;\n}\n\n/* Action buttons */\n\n.jp-FileBrowser-toolbar > .jp-FileAction > .jp-ToolbarButtonComponent > svg {\n  display: none;\n}\n\n.jp-FileBrowser-toolbar > #fileAction-delete {\n  background-color: var(--jp-error-color1);\n}\n\n.jp-FileBrowser-toolbar\n  .jp-ToolbarButtonComponent[data-command='filebrowser:delete']\n  .jp-ToolbarButtonComponent-label {\n  color: var(--jp-ui-inverse-font-color1);\n}\n\n.jp-FileBrowser-toolbar .jp-FileAction {\n  border: solid 1px var(--jp-border-color2);\n  margin: 1px;\n  min-height: var(--jp-private-toolbar-height);\n}\n\nbody[data-format='mobile'] #fileAction-placeholder {\n  display: none;\n}\n"],"sourceRoot":""}]);
// Exports
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (___CSS_LOADER_EXPORT___);


/***/ }),

/***/ 30546:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(82971);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11368);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__);
// Imports


var ___CSS_LOADER_EXPORT___ = _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default()((_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default()));
// Module
___CSS_LOADER_EXPORT___.push([module.id, `.jp-FileBrowser {
  height: 100%;
}

.lm-TabPanel {
  height: 100%;
}

.jp-TreePanel .lm-TabPanel-tabBar {
  overflow: visible;
  min-height: 32px;
  border-bottom: unset;
  height: var(--jp-private-toolbar-height);
}

.jp-TreePanel .lm-TabBar-content {
  height: 100%;
}

.jp-TreePanel .lm-TabBar-tab {
  flex: 0 1 auto;
  color: var(--jp-ui-font-color0);
  font-size: var(--jp-ui-font-size1);
  height: 100%;
}

.jp-TreePanel .lm-TabBar-tabLabel {
  padding-left: 5px;
  padding-right: 5px;
}

.jp-FileBrowser-toolbar.jp-Toolbar .jp-ToolbarButtonComponent {
  width: unset;
}

.jp-FileBrowser-toolbar > .jp-Toolbar-item {
  flex-direction: column;
  justify-content: center;
}

.jp-DropdownMenu .lm-MenuBar-itemIcon svg {
  vertical-align: sub;
}

jp-button[data-command='filebrowser:refresh'] .jp-ToolbarButtonComponent-label {
  display: none;
}

.jp-TreePanel .lm-TabBar-tabIcon svg {
  vertical-align: sub;
}
`, "",{"version":3,"sources":["webpack://./../packages/tree/style/base.css"],"names":[],"mappings":"AAAA;EACE,YAAY;AACd;;AAEA;EACE,YAAY;AACd;;AAEA;EACE,iBAAiB;EACjB,gBAAgB;EAChB,oBAAoB;EACpB,wCAAwC;AAC1C;;AAEA;EACE,YAAY;AACd;;AAEA;EACE,cAAc;EACd,+BAA+B;EAC/B,kCAAkC;EAClC,YAAY;AACd;;AAEA;EACE,iBAAiB;EACjB,kBAAkB;AACpB;;AAEA;EACE,YAAY;AACd;;AAEA;EACE,sBAAsB;EACtB,uBAAuB;AACzB;;AAEA;EACE,mBAAmB;AACrB;;AAEA;EACE,aAAa;AACf;;AAEA;EACE,mBAAmB;AACrB","sourcesContent":[".jp-FileBrowser {\n  height: 100%;\n}\n\n.lm-TabPanel {\n  height: 100%;\n}\n\n.jp-TreePanel .lm-TabPanel-tabBar {\n  overflow: visible;\n  min-height: 32px;\n  border-bottom: unset;\n  height: var(--jp-private-toolbar-height);\n}\n\n.jp-TreePanel .lm-TabBar-content {\n  height: 100%;\n}\n\n.jp-TreePanel .lm-TabBar-tab {\n  flex: 0 1 auto;\n  color: var(--jp-ui-font-color0);\n  font-size: var(--jp-ui-font-size1);\n  height: 100%;\n}\n\n.jp-TreePanel .lm-TabBar-tabLabel {\n  padding-left: 5px;\n  padding-right: 5px;\n}\n\n.jp-FileBrowser-toolbar.jp-Toolbar .jp-ToolbarButtonComponent {\n  width: unset;\n}\n\n.jp-FileBrowser-toolbar > .jp-Toolbar-item {\n  flex-direction: column;\n  justify-content: center;\n}\n\n.jp-DropdownMenu .lm-MenuBar-itemIcon svg {\n  vertical-align: sub;\n}\n\njp-button[data-command='filebrowser:refresh'] .jp-ToolbarButtonComponent-label {\n  display: none;\n}\n\n.jp-TreePanel .lm-TabBar-tabIcon svg {\n  vertical-align: sub;\n}\n"],"sourceRoot":""}]);
// Exports
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (___CSS_LOADER_EXPORT___);


/***/ }),

/***/ 99402:
/***/ ((module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(82971);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11368);
/* harmony import */ var _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1__);
// Imports


var ___CSS_LOADER_EXPORT___ = _node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_api_js__WEBPACK_IMPORTED_MODULE_1___default()((_node_modules_jupyterlab_builder_node_modules_css_loader_dist_runtime_sourceMaps_js__WEBPACK_IMPORTED_MODULE_0___default()));
// Module
___CSS_LOADER_EXPORT___.push([module.id, ``, "",{"version":3,"sources":[],"names":[],"mappings":"","sourceRoot":""}]);
// Exports
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (___CSS_LOADER_EXPORT___);


/***/ }),

/***/ 7413:
/***/ ((module) => {

"use strict";
module.exports = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAgAAAAFCAYAAAB4ka1VAAAAsElEQVQIHQGlAFr/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA7+r3zKmT0/+pk9P/7+r3zAAAAAAAAAAABAAAAAAAAAAA6OPzM+/q9wAAAAAA6OPzMwAAAAAAAAAAAgAAAAAAAAAAGR8NiRQaCgAZIA0AGR8NiQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQyoYJ/SY80UAAAAASUVORK5CYII=";

/***/ })

}]);
//# sourceMappingURL=8781.b7c620b7ec3a9833df1c.js.map?v=b7c620b7ec3a9833df1c
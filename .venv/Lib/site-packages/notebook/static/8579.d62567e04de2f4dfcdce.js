"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[8579],{

/***/ 88579:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

__webpack_require__.r(__webpack_exports__);
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   "default": () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(6646);
/* harmony import */ var _jupyterlab_application__WEBPACK_IMPORTED_MODULE_0___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__);
/* harmony import */ var _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(11922);
/* harmony import */ var _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_1___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_1__);
/* harmony import */ var _jupyterlab_console__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(64190);
/* harmony import */ var _jupyterlab_console__WEBPACK_IMPORTED_MODULE_2___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_console__WEBPACK_IMPORTED_MODULE_2__);
/* harmony import */ var _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(64275);
/* harmony import */ var _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__);
/* harmony import */ var _jupyterlab_docmanager__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(23925);
/* harmony import */ var _jupyterlab_docmanager__WEBPACK_IMPORTED_MODULE_4___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_docmanager__WEBPACK_IMPORTED_MODULE_4__);
/* harmony import */ var _jupyterlab_docregistry__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(81728);
/* harmony import */ var _jupyterlab_docregistry__WEBPACK_IMPORTED_MODULE_5___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_docregistry__WEBPACK_IMPORTED_MODULE_5__);
/* harmony import */ var _jupyterlab_mainmenu__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(81369);
/* harmony import */ var _jupyterlab_mainmenu__WEBPACK_IMPORTED_MODULE_6___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_mainmenu__WEBPACK_IMPORTED_MODULE_6__);
/* harmony import */ var _jupyterlab_rendermime__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(41338);
/* harmony import */ var _jupyterlab_rendermime__WEBPACK_IMPORTED_MODULE_7___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_rendermime__WEBPACK_IMPORTED_MODULE_7__);
/* harmony import */ var _jupyterlab_settingregistry__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(49531);
/* harmony import */ var _jupyterlab_settingregistry__WEBPACK_IMPORTED_MODULE_8___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_settingregistry__WEBPACK_IMPORTED_MODULE_8__);
/* harmony import */ var _jupyterlab_translation__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(25830);
/* harmony import */ var _jupyterlab_translation__WEBPACK_IMPORTED_MODULE_9___default = /*#__PURE__*/__webpack_require__.n(_jupyterlab_translation__WEBPACK_IMPORTED_MODULE_9__);
/* harmony import */ var _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(43223);
/* harmony import */ var _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10___default = /*#__PURE__*/__webpack_require__.n(_jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__);
/* harmony import */ var _jupyter_notebook_ui_components__WEBPACK_IMPORTED_MODULE_11__ = __webpack_require__(61646);
/* harmony import */ var _jupyter_notebook_ui_components__WEBPACK_IMPORTED_MODULE_11___default = /*#__PURE__*/__webpack_require__.n(_jupyter_notebook_ui_components__WEBPACK_IMPORTED_MODULE_11__);
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_12__ = __webpack_require__(5406);
/* harmony import */ var _lumino_coreutils__WEBPACK_IMPORTED_MODULE_12___default = /*#__PURE__*/__webpack_require__.n(_lumino_coreutils__WEBPACK_IMPORTED_MODULE_12__);
/* harmony import */ var _lumino_disposable__WEBPACK_IMPORTED_MODULE_13__ = __webpack_require__(38302);
/* harmony import */ var _lumino_disposable__WEBPACK_IMPORTED_MODULE_13___default = /*#__PURE__*/__webpack_require__.n(_lumino_disposable__WEBPACK_IMPORTED_MODULE_13__);
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_14__ = __webpack_require__(60920);
/* harmony import */ var _lumino_widgets__WEBPACK_IMPORTED_MODULE_14___default = /*#__PURE__*/__webpack_require__.n(_lumino_widgets__WEBPACK_IMPORTED_MODULE_14__);
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.















/**
 * A regular expression to match path to notebooks and documents
 */
const TREE_PATTERN = new RegExp('/(notebooks|edit)/(.*)');
/**
 * A regular expression to suppress the file extension from display for .ipynb files.
 */
const STRIP_IPYNB = /\.ipynb$/;
/**
 * The JupyterLab document manager plugin id.
 */
const JUPYTERLAB_DOCMANAGER_PLUGIN_ID = '@jupyterlab/docmanager-extension:plugin';
/**
 * The command IDs used by the application plugin.
 */
var CommandIDs;
(function (CommandIDs) {
    /**
     * Duplicate the current document and open the new document
     */
    CommandIDs.duplicate = 'application:duplicate';
    /**
     * Handle local links
     */
    CommandIDs.handleLink = 'application:handle-local-link';
    /**
     * Toggle Top Bar visibility
     */
    CommandIDs.toggleTop = 'application:toggle-top';
    /**
     * Toggle side panel visibility
     */
    CommandIDs.togglePanel = 'application:toggle-panel';
    /**
     * Toggle the Zen mode
     */
    CommandIDs.toggleZen = 'application:toggle-zen';
    /**
     * Open JupyterLab
     */
    CommandIDs.openLab = 'application:open-lab';
    /**
     * Open the tree page.
     */
    CommandIDs.openTree = 'application:open-tree';
    /**
     * Rename the current document
     */
    CommandIDs.rename = 'application:rename';
    /**
     * Resolve tree path
     */
    CommandIDs.resolveTree = 'application:resolve-tree';
})(CommandIDs || (CommandIDs = {}));
/**
 * Check if the application is dirty before closing the browser tab.
 */
const dirty = {
    id: '@jupyter-notebook/application-extension:dirty',
    description: 'Check if the application is dirty before closing the browser tab.',
    autoStart: true,
    requires: [_jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__.ILabStatus, _jupyterlab_translation__WEBPACK_IMPORTED_MODULE_9__.ITranslator],
    activate: (app, status, translator) => {
        if (!(app instanceof _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.NotebookApp)) {
            throw new Error(`${dirty.id} must be activated in Jupyter Notebook.`);
        }
        const trans = translator.load('notebook');
        const message = trans.__('Are you sure you want to exit Jupyter Notebook?\n\nAny unsaved changes will be lost.');
        window.addEventListener('beforeunload', (event) => {
            if (app.status.isDirty) {
                return (event.returnValue = message);
            }
        });
    },
};
/**
 * The application info.
 */
const info = {
    id: '@jupyter-notebook/application-extension:info',
    autoStart: true,
    provides: _jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__.JupyterLab.IInfo,
    activate: (app) => {
        if (!(app instanceof _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.NotebookApp)) {
            throw new Error(`${info.id} must be activated in Jupyter Notebook.`);
        }
        return app.info;
    },
};
/**
 * The logo plugin.
 */
const logo = {
    id: '@jupyter-notebook/application-extension:logo',
    description: 'The logo plugin.',
    autoStart: true,
    activate: (app) => {
        const baseUrl = _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.PageConfig.getBaseUrl();
        const node = document.createElement('a');
        node.href = `${baseUrl}tree`;
        node.target = '_blank';
        node.rel = 'noopener noreferrer';
        const logo = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_14__.Widget({ node });
        _jupyter_notebook_ui_components__WEBPACK_IMPORTED_MODULE_11__.jupyterIcon.element({
            container: node,
            elementPosition: 'center',
            padding: '2px 2px 2px 8px',
            height: '28px',
            width: 'auto',
            cursor: 'pointer',
            margin: 'auto',
        });
        logo.id = 'jp-NotebookLogo';
        app.shell.add(logo, 'top', { rank: 0 });
    },
};
/**
 * A plugin to open documents in the main area.
 */
const opener = {
    id: '@jupyter-notebook/application-extension:opener',
    description: 'A plugin to open documents in the main area.',
    autoStart: true,
    requires: [_jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__.IRouter, _jupyterlab_docmanager__WEBPACK_IMPORTED_MODULE_4__.IDocumentManager],
    optional: [_jupyterlab_settingregistry__WEBPACK_IMPORTED_MODULE_8__.ISettingRegistry],
    activate: (app, router, docManager, settingRegistry) => {
        const { commands, docRegistry } = app;
        const command = 'router:tree';
        commands.addCommand(command, {
            execute: (args) => {
                var _a;
                const parsed = args;
                const matches = (_a = parsed.path.match(TREE_PATTERN)) !== null && _a !== void 0 ? _a : [];
                const [, , path] = matches;
                if (!path) {
                    return;
                }
                app.started.then(async () => {
                    var _a;
                    const file = decodeURIComponent(path);
                    const urlParams = new URLSearchParams(parsed.search);
                    let defaultFactory = docRegistry.defaultWidgetFactory(path).name;
                    // Explicitly get the default viewers from the settings because
                    // JupyterLab might not have had the time to load the settings yet (race condition)
                    // Relevant code: https://github.com/jupyterlab/jupyterlab/blob/d56ff811f39b3c10c6d8b6eb27a94624b753eb53/packages/docmanager-extension/src/index.tsx#L265-L293
                    if (settingRegistry) {
                        const settings = await settingRegistry.load(JUPYTERLAB_DOCMANAGER_PLUGIN_ID);
                        const defaultViewers = settings.get('defaultViewers').composite;
                        // get the file types for the path
                        const types = docRegistry.getFileTypesForPath(path);
                        // for each file type, check if there is a default viewer and if it
                        // is available in the docRegistry. If it is the case, use it as the
                        // default factory
                        types.forEach((ft) => {
                            if (defaultViewers[ft.name] !== undefined &&
                                docRegistry.getWidgetFactory(defaultViewers[ft.name])) {
                                defaultFactory = defaultViewers[ft.name];
                            }
                        });
                    }
                    const factory = (_a = urlParams.get('factory')) !== null && _a !== void 0 ? _a : defaultFactory;
                    docManager.open(file, factory, undefined, {
                        ref: '_noref',
                    });
                });
            },
        });
        router.register({ command, pattern: TREE_PATTERN });
    },
};
/**
 * A plugin to customize menus
 */
const menus = {
    id: '@jupyter-notebook/application-extension:menus',
    description: 'A plugin to customize menus.',
    requires: [_jupyterlab_mainmenu__WEBPACK_IMPORTED_MODULE_6__.IMainMenu],
    autoStart: true,
    activate: (app, menu) => {
        // always disable the Tabs menu
        menu.tabsMenu.dispose();
        const page = _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.PageConfig.getOption('notebookPage');
        switch (page) {
            case 'consoles':
            case 'terminals':
            case 'tree':
                menu.editMenu.dispose();
                menu.kernelMenu.dispose();
                menu.runMenu.dispose();
                break;
            case 'edit':
                menu.kernelMenu.dispose();
                menu.runMenu.dispose();
                break;
            default:
                break;
        }
    },
};
/**
 * A plugin to provide a spacer at rank 900 in the menu area
 */
const menuSpacer = {
    id: '@jupyter-notebook/application-extension:menu-spacer',
    description: 'A plugin to provide a spacer at rank 900 in the menu area.',
    autoStart: true,
    activate: (app) => {
        const menu = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_14__.Widget();
        menu.id = _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_1__.DOMUtils.createDomID();
        menu.addClass('jp-NotebookSpacer');
        app.shell.add(menu, 'menu', { rank: 900 });
    },
};
/**
 * Add commands to open the tree and running pages.
 */
const pages = {
    id: '@jupyter-notebook/application-extension:pages',
    description: 'Add commands to open the tree and running pages.',
    autoStart: true,
    requires: [_jupyterlab_translation__WEBPACK_IMPORTED_MODULE_9__.ITranslator],
    optional: [_jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_1__.ICommandPalette],
    activate: (app, translator, palette) => {
        const trans = translator.load('notebook');
        const baseUrl = _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.PageConfig.getBaseUrl();
        app.commands.addCommand(CommandIDs.openLab, {
            label: trans.__('Open JupyterLab'),
            execute: () => {
                window.open(_jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.URLExt.join(baseUrl, 'lab'));
            },
        });
        const page = _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.PageConfig.getOption('notebookPage');
        app.commands.addCommand(CommandIDs.openTree, {
            label: trans.__('File Browser'),
            execute: () => {
                if (page === 'tree') {
                    app.commands.execute('filebrowser:activate');
                }
                else {
                    window.open(_jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.URLExt.join(baseUrl, 'tree'));
                }
            },
        });
        if (palette) {
            palette.addItem({ command: CommandIDs.openLab, category: 'View' });
            palette.addItem({ command: CommandIDs.openTree, category: 'View' });
        }
    },
};
/**
 * A plugin to open paths in new browser tabs.
 */
const pathOpener = {
    id: '@jupyter-notebook/application-extension:path-opener',
    description: 'A plugin to open paths in new browser tabs.',
    autoStart: true,
    provides: _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.INotebookPathOpener,
    activate: (app) => {
        return _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.defaultNotebookPathOpener;
    },
};
/**
 * The default paths for a Jupyter Notebook app.
 */
const paths = {
    id: '@jupyter-notebook/application-extension:paths',
    description: 'The default paths for a Jupyter Notebook app.',
    autoStart: true,
    provides: _jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__.JupyterFrontEnd.IPaths,
    activate: (app) => {
        if (!(app instanceof _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.NotebookApp)) {
            throw new Error(`${paths.id} must be activated in Jupyter Notebook.`);
        }
        return app.paths;
    },
};
/**
 * A plugin providing a rendermime registry.
 */
const rendermime = {
    id: '@jupyter-notebook/application-extension:rendermime',
    description: 'A plugin providing a rendermime registry.',
    autoStart: true,
    provides: _jupyterlab_rendermime__WEBPACK_IMPORTED_MODULE_7__.IRenderMimeRegistry,
    optional: [
        _jupyterlab_docmanager__WEBPACK_IMPORTED_MODULE_4__.IDocumentManager,
        _jupyterlab_rendermime__WEBPACK_IMPORTED_MODULE_7__.ILatexTypesetter,
        _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_1__.ISanitizer,
        _jupyterlab_rendermime__WEBPACK_IMPORTED_MODULE_7__.IMarkdownParser,
        _jupyterlab_translation__WEBPACK_IMPORTED_MODULE_9__.ITranslator,
        _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.INotebookPathOpener,
    ],
    activate: (app, docManager, latexTypesetter, sanitizer, markdownParser, translator, notebookPathOpener) => {
        const trans = (translator !== null && translator !== void 0 ? translator : _jupyterlab_translation__WEBPACK_IMPORTED_MODULE_9__.nullTranslator).load('jupyterlab');
        const opener = notebookPathOpener !== null && notebookPathOpener !== void 0 ? notebookPathOpener : _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.defaultNotebookPathOpener;
        if (docManager) {
            app.commands.addCommand(CommandIDs.handleLink, {
                label: trans.__('Handle Local Link'),
                execute: (args) => {
                    const path = args['path'];
                    if (path === undefined || path === null) {
                        return;
                    }
                    return docManager.services.contents
                        .get(path, { content: false })
                        .then((model) => {
                        const baseUrl = _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.PageConfig.getBaseUrl();
                        opener.open({
                            prefix: _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.URLExt.join(baseUrl, 'tree'),
                            path: model.path,
                            target: '_blank',
                        });
                    });
                },
            });
        }
        return new _jupyterlab_rendermime__WEBPACK_IMPORTED_MODULE_7__.RenderMimeRegistry({
            initialFactories: _jupyterlab_rendermime__WEBPACK_IMPORTED_MODULE_7__.standardRendererFactories,
            linkHandler: !docManager
                ? undefined
                : {
                    handleLink: (node, path, id) => {
                        // If node has the download attribute explicitly set, use the
                        // default browser downloading behavior.
                        if (node.tagName === 'A' && node.hasAttribute('download')) {
                            return;
                        }
                        app.commandLinker.connectNode(node, CommandIDs.handleLink, {
                            path,
                            id,
                        });
                    },
                },
            latexTypesetter: latexTypesetter !== null && latexTypesetter !== void 0 ? latexTypesetter : undefined,
            markdownParser: markdownParser !== null && markdownParser !== void 0 ? markdownParser : undefined,
            translator: translator !== null && translator !== void 0 ? translator : undefined,
            sanitizer: sanitizer !== null && sanitizer !== void 0 ? sanitizer : undefined,
        });
    },
};
/**
 * The default Jupyter Notebook application shell.
 */
const shell = {
    id: '@jupyter-notebook/application-extension:shell',
    description: 'The default Jupyter Notebook application shell.',
    autoStart: true,
    provides: _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.INotebookShell,
    optional: [_jupyterlab_settingregistry__WEBPACK_IMPORTED_MODULE_8__.ISettingRegistry],
    activate: (app, settingRegistry) => {
        if (!(app.shell instanceof _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.NotebookShell)) {
            throw new Error(`${shell.id} did not find a NotebookShell instance.`);
        }
        const notebookShell = app.shell;
        if (settingRegistry) {
            settingRegistry
                .load(shell.id)
                .then((settings) => {
                // Add a layer of customization to support app shell mode
                const customLayout = settings.composite['layout'];
                // Restore the layout.
                void notebookShell.restoreLayout(customLayout);
            })
                .catch((reason) => {
                console.error('Fail to load settings for the layout restorer.');
                console.error(reason);
            });
        }
        return notebookShell;
    },
};
/**
 * The default splash screen provider.
 */
const splash = {
    id: '@jupyter-notebook/application-extension:splash',
    description: 'Provides an empty splash screen.',
    autoStart: true,
    provides: _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_1__.ISplashScreen,
    activate: (app) => {
        const { restored } = app;
        const splash = document.createElement('div');
        splash.style.position = 'absolute';
        splash.style.width = '100%';
        splash.style.height = '100%';
        splash.style.zIndex = '10';
        return {
            show: (light = true) => {
                splash.style.backgroundColor = light ? 'white' : '#111111';
                document.body.appendChild(splash);
                return new _lumino_disposable__WEBPACK_IMPORTED_MODULE_13__.DisposableDelegate(async () => {
                    await restored;
                    document.body.removeChild(splash);
                });
            },
        };
    },
};
/**
 * The default JupyterLab application status provider.
 */
const status = {
    id: '@jupyter-notebook/application-extension:status',
    description: 'The default JupyterLab application status provider.',
    autoStart: true,
    provides: _jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__.ILabStatus,
    activate: (app) => {
        if (!(app instanceof _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.NotebookApp)) {
            throw new Error(`${status.id} must be activated in Jupyter Notebook.`);
        }
        return app.status;
    },
};
/**
 * A plugin to display the document title in the browser tab title
 */
const tabTitle = {
    id: '@jupyter-notebook/application-extension:tab-title',
    description: 'A plugin to display the document title in the browser tab title.',
    autoStart: true,
    requires: [_jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.INotebookShell],
    activate: (app, shell) => {
        const setTabTitle = () => {
            const current = shell.currentWidget;
            if (current instanceof _jupyterlab_console__WEBPACK_IMPORTED_MODULE_2__.ConsolePanel) {
                const update = () => {
                    const title = current.sessionContext.path || current.sessionContext.name;
                    const basename = _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.PathExt.basename(title);
                    // Strip the ".ipynb" suffix from filenames for display in tab titles.
                    document.title = basename.replace(STRIP_IPYNB, '');
                };
                current.sessionContext.sessionChanged.connect(update);
                update();
                return;
            }
            else if (current instanceof _jupyterlab_docregistry__WEBPACK_IMPORTED_MODULE_5__.DocumentWidget) {
                const update = () => {
                    const basename = _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.PathExt.basename(current.context.path);
                    document.title = basename.replace(STRIP_IPYNB, '');
                };
                current.context.pathChanged.connect(update);
                update();
            }
        };
        shell.currentChanged.connect(setTabTitle);
        setTabTitle();
    },
};
/**
 * A plugin to display and rename the title of a file
 */
const title = {
    id: '@jupyter-notebook/application-extension:title',
    description: 'A plugin to display and rename the title of a file.',
    autoStart: true,
    requires: [_jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.INotebookShell, _jupyterlab_translation__WEBPACK_IMPORTED_MODULE_9__.ITranslator],
    optional: [_jupyterlab_docmanager__WEBPACK_IMPORTED_MODULE_4__.IDocumentManager, _jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__.IRouter, _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_1__.IToolbarWidgetRegistry],
    activate: (app, shell, translator, docManager, router, toolbarRegistry) => {
        const { commands } = app;
        const trans = translator.load('notebook');
        const node = document.createElement('div');
        if (toolbarRegistry) {
            toolbarRegistry.addFactory('TopBar', 'widgetTitle', (toolbar) => {
                const widget = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_14__.Widget({ node });
                widget.id = 'jp-title';
                return widget;
            });
        }
        const addTitle = async () => {
            const current = shell.currentWidget;
            if (!current || !(current instanceof _jupyterlab_docregistry__WEBPACK_IMPORTED_MODULE_5__.DocumentWidget)) {
                return;
            }
            if (node.children.length > 0) {
                return;
            }
            const h = document.createElement('h1');
            h.textContent = current.title.label.replace(STRIP_IPYNB, '');
            node.appendChild(h);
            node.style.marginLeft = '10px';
            if (!docManager) {
                return;
            }
            const isEnabled = () => {
                const { currentWidget } = shell;
                return !!(currentWidget && docManager.contextForWidget(currentWidget));
            };
            commands.addCommand(CommandIDs.duplicate, {
                label: () => trans.__('Duplicate'),
                isEnabled,
                execute: async () => {
                    if (!isEnabled()) {
                        return;
                    }
                    // Duplicate the file, and open the new file.
                    const result = await docManager.duplicate(current.context.path);
                    await commands.execute('docmanager:open', { path: result.path });
                },
            });
            commands.addCommand(CommandIDs.rename, {
                label: () => trans.__('Rename…'),
                isEnabled,
                execute: async () => {
                    var _a;
                    if (!isEnabled()) {
                        return;
                    }
                    const result = await (0,_jupyterlab_docmanager__WEBPACK_IMPORTED_MODULE_4__.renameDialog)(docManager, current.context);
                    // activate the current widget to bring the focus
                    if (current) {
                        current.activate();
                    }
                    if (result === null) {
                        return;
                    }
                    const newPath = current.context.path;
                    const basename = _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.PathExt.basename(newPath);
                    h.textContent = basename.replace(STRIP_IPYNB, '');
                    if (!router) {
                        return;
                    }
                    const matches = (_a = router.current.path.match(TREE_PATTERN)) !== null && _a !== void 0 ? _a : [];
                    const [, route, path] = matches;
                    if (!route || !path) {
                        return;
                    }
                    const encoded = encodeURIComponent(newPath);
                    router.navigate(`/${route}/${encoded}`, {
                        skipRouting: true,
                    });
                },
            });
            node.onclick = async () => {
                void commands.execute(CommandIDs.rename);
            };
        };
        shell.currentChanged.connect(addTitle);
        void addTitle();
    },
};
/**
 * Plugin to toggle the top header visibility.
 */
const topVisibility = {
    id: '@jupyter-notebook/application-extension:top',
    description: 'Plugin to toggle the top header visibility.',
    requires: [_jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.INotebookShell, _jupyterlab_translation__WEBPACK_IMPORTED_MODULE_9__.ITranslator],
    optional: [_jupyterlab_settingregistry__WEBPACK_IMPORTED_MODULE_8__.ISettingRegistry, _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_1__.ICommandPalette],
    activate: (app, notebookShell, translator, settingRegistry, palette) => {
        const trans = translator.load('notebook');
        const top = notebookShell.top;
        const pluginId = topVisibility.id;
        app.commands.addCommand(CommandIDs.toggleTop, {
            label: trans.__('Show Header'),
            execute: () => {
                top.setHidden(top.isVisible);
                if (settingRegistry) {
                    void settingRegistry.set(pluginId, 'visible', top.isVisible ? 'yes' : 'no');
                }
            },
            isToggled: () => top.isVisible,
        });
        let adjustToScreen = false;
        if (settingRegistry) {
            const loadSettings = settingRegistry.load(pluginId);
            const updateSettings = (settings) => {
                // 'visible' property from user preferences or default settings
                let visible = settings.get('visible').composite;
                if (settings.user.visible !== undefined) {
                    visible = settings.user.visible;
                }
                top.setHidden(visible === 'no');
                // adjust to screen from user preferences or default settings
                adjustToScreen = visible === 'automatic';
            };
            Promise.all([loadSettings, app.restored])
                .then(([settings]) => {
                updateSettings(settings);
                settings.changed.connect((settings) => {
                    updateSettings(settings);
                });
            })
                .catch((reason) => {
                console.error(reason.message);
            });
        }
        if (palette) {
            palette.addItem({ command: CommandIDs.toggleTop, category: 'View' });
        }
        const onChanged = () => {
            if (!adjustToScreen) {
                return;
            }
            if (app.format === 'desktop') {
                notebookShell.expandTop();
            }
            else {
                notebookShell.collapseTop();
            }
        };
        // listen on format change (mobile and desktop) to make the view more compact
        app.formatChanged.connect(onChanged);
    },
    autoStart: true,
};
/**
 * Plugin to toggle the left or right side panel's visibility.
 */
const sidePanelVisibility = {
    id: '@jupyter-notebook/application-extension:sidepanel',
    description: 'Plugin to toggle the visibility of left or right side panel.',
    requires: [_jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.INotebookShell, _jupyterlab_translation__WEBPACK_IMPORTED_MODULE_9__.ITranslator],
    optional: [_jupyterlab_mainmenu__WEBPACK_IMPORTED_MODULE_6__.IMainMenu, _jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_1__.ICommandPalette],
    autoStart: true,
    activate: (app, notebookShell, translator, menu, palette) => {
        const trans = translator.load('notebook');
        /* Arguments for togglePanel command:
         * side, left or right area
         * title, widget title to show in the menu
         * id, widget ID to activate in the side panel
         */
        app.commands.addCommand(CommandIDs.togglePanel, {
            label: (args) => args['title'],
            caption: (args) => {
                // We do not substitute the parameter into the string because the parameter is not
                // localized (e.g., it is always 'left') even though the string is localized.
                if (args['side'] === 'left') {
                    return trans.__('Show %1 in the left sidebar', args['title']);
                }
                else if (args['side'] === 'right') {
                    return trans.__('Show %1 in the right sidebar', args['title']);
                }
                return trans.__('Show %1 in the sidebar', args['title']);
            },
            execute: (args) => {
                var _a, _b;
                switch (args['side']) {
                    case 'left':
                        if (notebookShell.leftCollapsed) {
                            notebookShell.expandLeft(args.id);
                        }
                        else if (((_a = notebookShell.leftHandler.currentWidget) === null || _a === void 0 ? void 0 : _a.id) !== args.id) {
                            notebookShell.expandLeft(args.id);
                        }
                        else {
                            notebookShell.collapseLeft();
                            if (notebookShell.currentWidget) {
                                notebookShell.activateById(notebookShell.currentWidget.id);
                            }
                        }
                        break;
                    case 'right':
                        if (notebookShell.rightCollapsed) {
                            notebookShell.expandRight(args.id);
                        }
                        else if (((_b = notebookShell.rightHandler.currentWidget) === null || _b === void 0 ? void 0 : _b.id) !== args.id) {
                            notebookShell.expandRight(args.id);
                        }
                        else {
                            notebookShell.collapseRight();
                            if (notebookShell.currentWidget) {
                                notebookShell.activateById(notebookShell.currentWidget.id);
                            }
                        }
                        break;
                }
            },
            isToggled: (args) => {
                switch (args['side']) {
                    case 'left': {
                        if (notebookShell.leftCollapsed) {
                            return false;
                        }
                        const currentWidget = notebookShell.leftHandler.currentWidget;
                        if (!currentWidget) {
                            return false;
                        }
                        return currentWidget.id === args['id'];
                    }
                    case 'right': {
                        if (notebookShell.rightCollapsed) {
                            return false;
                        }
                        const currentWidget = notebookShell.rightHandler.currentWidget;
                        if (!currentWidget) {
                            return false;
                        }
                        return currentWidget.id === args['id'];
                    }
                }
                return false;
            },
        });
        const sidePanelMenu = {
            left: null,
            right: null,
        };
        /**
         * The function which adds entries to the View menu for each widget of a side panel.
         *
         * @param area - 'left' or 'right', the area of the side panel.
         * @param entryLabel - the name of the main entry in the View menu for that side panel.
         * @returns - The disposable menu added to the View menu or null.
         */
        const updateMenu = (area, entryLabel) => {
            var _a;
            if (menu === null) {
                return null;
            }
            // Remove the previous menu entry for this side panel.
            (_a = sidePanelMenu[area]) === null || _a === void 0 ? void 0 : _a.dispose();
            // Creates a new menu entry and populates it with side panel widgets.
            const newMenu = new _lumino_widgets__WEBPACK_IMPORTED_MODULE_14__.Menu({ commands: app.commands });
            newMenu.title.label = entryLabel;
            const widgets = notebookShell.widgets(area);
            let menuToAdd = false;
            for (const widget of widgets) {
                newMenu.addItem({
                    command: CommandIDs.togglePanel,
                    args: {
                        side: area,
                        title: `Show ${widget.title.caption}`,
                        id: widget.id,
                    },
                });
                menuToAdd = true;
            }
            // If there are widgets, add the menu to the main menu entry.
            if (menuToAdd) {
                sidePanelMenu[area] = menu.viewMenu.addItem({
                    type: 'submenu',
                    submenu: newMenu,
                });
            }
        };
        app.restored.then(() => {
            // Create menu entries for the left and right panel.
            if (menu) {
                const getSidePanelLabel = (area) => {
                    if (area === 'left') {
                        return trans.__('Left Sidebar');
                    }
                    else {
                        return trans.__('Right Sidebar');
                    }
                };
                const leftArea = notebookShell.leftHandler.area;
                const leftLabel = getSidePanelLabel(leftArea);
                updateMenu(leftArea, leftLabel);
                const rightArea = notebookShell.rightHandler.area;
                const rightLabel = getSidePanelLabel(rightArea);
                updateMenu(rightArea, rightLabel);
                const handleSidePanelChange = (sidePanel, widget) => {
                    const label = getSidePanelLabel(sidePanel.area);
                    updateMenu(sidePanel.area, label);
                };
                notebookShell.leftHandler.widgetAdded.connect(handleSidePanelChange);
                notebookShell.leftHandler.widgetRemoved.connect(handleSidePanelChange);
                notebookShell.rightHandler.widgetAdded.connect(handleSidePanelChange);
                notebookShell.rightHandler.widgetRemoved.connect(handleSidePanelChange);
            }
            // Add palette entries for side panels.
            if (palette) {
                const sidePanelPalette = new _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.SidePanelPalette({
                    commandPalette: palette,
                    command: CommandIDs.togglePanel,
                });
                notebookShell.leftHandler.widgets.forEach((widget) => {
                    sidePanelPalette.addItem(widget, notebookShell.leftHandler.area);
                });
                notebookShell.rightHandler.widgets.forEach((widget) => {
                    sidePanelPalette.addItem(widget, notebookShell.rightHandler.area);
                });
                // Update menu and palette when widgets are added or removed from side panels.
                notebookShell.leftHandler.widgetAdded.connect((sidePanel, widget) => {
                    sidePanelPalette.addItem(widget, sidePanel.area);
                });
                notebookShell.leftHandler.widgetRemoved.connect((sidePanel, widget) => {
                    sidePanelPalette.removeItem(widget, sidePanel.area);
                });
                notebookShell.rightHandler.widgetAdded.connect((sidePanel, widget) => {
                    sidePanelPalette.addItem(widget, sidePanel.area);
                });
                notebookShell.rightHandler.widgetRemoved.connect((sidePanel, widget) => {
                    sidePanelPalette.removeItem(widget, sidePanel.area);
                });
            }
        });
    },
};
/**
 * The default tree route resolver plugin.
 */
const tree = {
    id: '@jupyter-notebook/application-extension:tree-resolver',
    description: 'The default tree route resolver plugin.',
    autoStart: true,
    requires: [_jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__.IRouter],
    provides: _jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__.JupyterFrontEnd.ITreeResolver,
    activate: (app, router) => {
        const { commands } = app;
        const set = new _lumino_disposable__WEBPACK_IMPORTED_MODULE_13__.DisposableSet();
        const delegate = new _lumino_coreutils__WEBPACK_IMPORTED_MODULE_12__.PromiseDelegate();
        const treePattern = new RegExp('/(/tree/.*)?');
        set.add(commands.addCommand(CommandIDs.resolveTree, {
            execute: (async (args) => {
                var _a;
                if (set.isDisposed) {
                    return;
                }
                const query = _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.URLExt.queryStringToObject((_a = args.search) !== null && _a !== void 0 ? _a : '');
                const browser = query['file-browser-path'] || '';
                // Remove the file browser path from the query string.
                delete query['file-browser-path'];
                // Clean up artifacts immediately upon routing.
                set.dispose();
                delegate.resolve({ browser, file: _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.PageConfig.getOption('treePath') });
            }),
        }));
        set.add(router.register({ command: CommandIDs.resolveTree, pattern: treePattern }));
        // If a route is handled by the router without the tree command being
        // invoked, resolve to `null` and clean up artifacts.
        const listener = () => {
            if (set.isDisposed) {
                return;
            }
            set.dispose();
            delegate.resolve(null);
        };
        router.routed.connect(listener);
        set.add(new _lumino_disposable__WEBPACK_IMPORTED_MODULE_13__.DisposableDelegate(() => {
            router.routed.disconnect(listener);
        }));
        return { paths: delegate.promise };
    },
};
/**
 * Plugin to update tree path.
 */
const treePathUpdater = {
    id: '@jupyter-notebook/application-extension:tree-updater',
    description: 'Plugin to update tree path.',
    requires: [_jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__.IRouter],
    provides: _jupyterlab_application__WEBPACK_IMPORTED_MODULE_0__.ITreePathUpdater,
    activate: (app, router) => {
        function updateTreePath(treePath) {
            if (treePath !== _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.PageConfig.getOption('treePath')) {
                const path = _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.URLExt.join(_jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.PageConfig.getOption('baseUrl') || '/', 'tree', _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.URLExt.encodeParts(treePath));
                router.navigate(path, { skipRouting: true });
                // Persist the new tree path to PageConfig as it is used elsewhere at runtime.
                _jupyterlab_coreutils__WEBPACK_IMPORTED_MODULE_3__.PageConfig.setOption('treePath', treePath);
            }
        }
        return updateTreePath;
    },
    autoStart: true,
};
/**
 * Translator plugin
 */
const translator = {
    id: '@jupyter-notebook/application-extension:translator',
    description: 'Translator plugin',
    requires: [_jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.INotebookShell, _jupyterlab_translation__WEBPACK_IMPORTED_MODULE_9__.ITranslator],
    autoStart: true,
    activate: (app, notebookShell, translator) => {
        notebookShell.translator = translator;
    },
};
/**
 * Zen mode plugin
 */
const zen = {
    id: '@jupyter-notebook/application-extension:zen',
    description: 'Zen mode plugin.',
    autoStart: true,
    requires: [_jupyterlab_translation__WEBPACK_IMPORTED_MODULE_9__.ITranslator],
    optional: [_jupyterlab_apputils__WEBPACK_IMPORTED_MODULE_1__.ICommandPalette, _jupyter_notebook_application__WEBPACK_IMPORTED_MODULE_10__.INotebookShell],
    activate: (app, translator, palette, notebookShell) => {
        const { commands } = app;
        const elem = document.documentElement;
        const trans = translator.load('notebook');
        const toggleOn = () => {
            notebookShell === null || notebookShell === void 0 ? void 0 : notebookShell.collapseTop();
            notebookShell === null || notebookShell === void 0 ? void 0 : notebookShell.menu.setHidden(true);
            zenModeEnabled = true;
        };
        const toggleOff = () => {
            notebookShell === null || notebookShell === void 0 ? void 0 : notebookShell.expandTop();
            notebookShell === null || notebookShell === void 0 ? void 0 : notebookShell.menu.setHidden(false);
            zenModeEnabled = false;
        };
        let zenModeEnabled = false;
        commands.addCommand(CommandIDs.toggleZen, {
            label: trans.__('Toggle Zen Mode'),
            execute: () => {
                if (!zenModeEnabled) {
                    elem.requestFullscreen();
                    toggleOn();
                }
                else {
                    document.exitFullscreen();
                    toggleOff();
                }
            },
        });
        document.addEventListener('fullscreenchange', () => {
            if (!document.fullscreenElement) {
                toggleOff();
            }
        });
        if (palette) {
            palette.addItem({ command: CommandIDs.toggleZen, category: 'Mode' });
        }
    },
};
/**
 * Export the plugins as default.
 */
const plugins = [
    dirty,
    info,
    logo,
    menus,
    menuSpacer,
    opener,
    pages,
    pathOpener,
    paths,
    rendermime,
    shell,
    sidePanelVisibility,
    splash,
    status,
    tabTitle,
    title,
    topVisibility,
    tree,
    treePathUpdater,
    translator,
    zen,
];
/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (plugins);


/***/ })

}]);
//# sourceMappingURL=8579.d62567e04de2f4dfcdce.js.map?v=d62567e04de2f4dfcdce
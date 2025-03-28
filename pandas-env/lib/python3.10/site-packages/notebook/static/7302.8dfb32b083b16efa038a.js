"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[7302,3768],{

/***/ 83768:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

// ESM COMPAT FLAG
__webpack_require__.r(__webpack_exports__);

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  "default": () => (/* binding */ lib)
});

// EXTERNAL MODULE: consume shared module (default) @jupyterlab/apputils@~4.4.4 (singleton) (fallback: ../node_modules/@jupyterlab/apputils/lib/index.js)
var index_js_ = __webpack_require__(11363);
// EXTERNAL MODULE: consume shared module (default) @jupyterlab/coreutils@~6.3.4 (singleton) (fallback: ../node_modules/@jupyterlab/coreutils/lib/index.js)
var lib_index_js_ = __webpack_require__(76107);
// EXTERNAL MODULE: consume shared module (default) @jupyterlab/filebrowser@~4.3.4 (singleton) (fallback: ../node_modules/@jupyterlab/filebrowser/lib/index.js)
var filebrowser_lib_index_js_ = __webpack_require__(23066);
// EXTERNAL MODULE: consume shared module (default) @jupyterlab/settingregistry@~4.3.4 (singleton) (fallback: ../node_modules/@jupyterlab/settingregistry/lib/index.js)
var settingregistry_lib_index_js_ = __webpack_require__(57534);
// EXTERNAL MODULE: consume shared module (default) @jupyterlab/running@^4.3.4 (strict) (fallback: ../node_modules/@jupyterlab/running/lib/index.js)
var running_lib_index_js_ = __webpack_require__(80407);
// EXTERNAL MODULE: consume shared module (default) @jupyterlab/settingeditor@~4.3.4 (singleton) (fallback: ../node_modules/@jupyterlab/settingeditor/lib/index.js)
var settingeditor_lib_index_js_ = __webpack_require__(2881);
// EXTERNAL MODULE: consume shared module (default) @jupyterlab/translation@~4.3.4 (singleton) (fallback: ../node_modules/@jupyterlab/translation/lib/index.js)
var translation_lib_index_js_ = __webpack_require__(7801);
// EXTERNAL MODULE: consume shared module (default) @jupyterlab/ui-components@~4.3.4 (singleton) (fallback: ../node_modules/@jupyterlab/ui-components/lib/index.js)
var ui_components_lib_index_js_ = __webpack_require__(29072);
// EXTERNAL MODULE: consume shared module (default) @lumino/signaling@~2.1.3 (singleton) (fallback: ../node_modules/@lumino/signaling/dist/index.es6.js)
var index_es6_js_ = __webpack_require__(2159);
// EXTERNAL MODULE: consume shared module (default) @lumino/widgets@~2.5.0 (singleton) (fallback: ../node_modules/@lumino/widgets/dist/index.es6.js)
var dist_index_es6_js_ = __webpack_require__(2260);
// EXTERNAL MODULE: consume shared module (default) @jupyter-notebook/tree@~7.3.2 (singleton) (fallback: ../packages/tree/lib/index.js)
var tree_lib_index_js_ = __webpack_require__(8412);
// EXTERNAL MODULE: consume shared module (default) react@~18.2.0 (singleton) (fallback: ../node_modules/react/index.js)
var react_index_js_ = __webpack_require__(78156);
var react_index_js_default = /*#__PURE__*/__webpack_require__.n(react_index_js_);
;// CONCATENATED MODULE: ../packages/tree-extension/lib/fileactions.js
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.


class FilesActionButtons {
    /**
     * The constructor of FilesActionButtons.
     * @param options
     */
    constructor(options) {
        /**
         * Triggered when the selection change in file browser.
         */
        this._onSelectionChanged = () => {
            var _a, _b, _c, _d, _e, _f;
            const selectedItems = Array.from(this._browser.selectedItems());
            const selection = selectedItems.length > 0;
            const oneFolder = selectedItems.some((item) => item.type === 'directory');
            (_a = this._widgets.get('placeholder')) === null || _a === void 0 ? void 0 : _a.setHidden(selection);
            (_b = this._widgets.get('delete')) === null || _b === void 0 ? void 0 : _b.setHidden(!selection);
            (_c = this._widgets.get('duplicate')) === null || _c === void 0 ? void 0 : _c.setHidden(!selection || oneFolder);
            (_d = this._widgets.get('download')) === null || _d === void 0 ? void 0 : _d.setHidden(!selection || oneFolder);
            (_e = this._widgets.get('open')) === null || _e === void 0 ? void 0 : _e.setHidden(!selection || oneFolder);
            (_f = this._widgets.get('rename')) === null || _f === void 0 ? void 0 : _f.setHidden(selectedItems.length !== 1);
        };
        this._widgets = new Map();
        this._browser = options.browser;
        const { commands, selectionChanged, translator } = options;
        const trans = translator.load('notebook');
        // Placeholder, when no file is selected.
        const placeholder = index_js_.ReactWidget.create(react_index_js_default().createElement("div", { key: 'placeholder' }, trans.__('Select items to perform actions on them.')));
        placeholder.id = 'fileAction-placeholder';
        this._widgets.set('placeholder', placeholder);
        // The action buttons.
        const actions = ['open', 'download', 'rename', 'duplicate', 'delete'];
        actions.forEach((action) => {
            const widget = index_js_.ReactWidget.create(react_index_js_default().createElement(index_js_.CommandToolbarButtonComponent, { key: action, commands: commands, id: `filebrowser:${action}`, args: { toolbar: true }, icon: undefined }));
            widget.id = `fileAction-${action}`;
            widget.addClass('jp-ToolbarButton');
            widget.addClass('jp-FileAction');
            this._widgets.set(action, widget);
        });
        selectionChanged.connect(this._onSelectionChanged, this);
        this._onSelectionChanged();
    }
    /**
     * Return an iterator with all the action widgets.
     */
    get widgets() {
        return this._widgets.values();
    }
}

;// CONCATENATED MODULE: ../packages/tree-extension/lib/index.js
// Copyright (c) Jupyter Development Team.
// Distributed under the terms of the Modified BSD License.












/**
 * The file browser factory.
 */
const FILE_BROWSER_FACTORY = 'FileBrowser';
/**
 * The file browser plugin id.
 */
const FILE_BROWSER_PLUGIN_ID = '@jupyterlab/filebrowser-extension:browser';
/**
 * The namespace for command IDs.
 */
var CommandIDs;
(function (CommandIDs) {
    // The command to activate the filebrowser widget in tree view.
    CommandIDs.activate = 'filebrowser:activate';
    // Activate the file filter in the file browser
    CommandIDs.toggleFileFilter = 'filebrowser:toggle-file-filter';
})(CommandIDs || (CommandIDs = {}));
/**
 * Plugin to add extra commands to the file browser to create
 * new notebooks, files, consoles and terminals
 */
const createNew = {
    id: '@jupyter-notebook/tree-extension:new',
    description: 'Plugin to add extra commands to the file browser to create new notebooks, files, consoles and terminals.',
    requires: [translation_lib_index_js_.ITranslator],
    optional: [index_js_.IToolbarWidgetRegistry],
    autoStart: true,
    activate: (app, translator, toolbarRegistry) => {
        var _a;
        const { commands, serviceManager } = app;
        const trans = translator.load('notebook');
        const overflowOptions = {
            overflowMenuOptions: { isVisible: false },
        };
        const menubar = new dist_index_es6_js_.MenuBar(overflowOptions);
        const newMenu = new dist_index_es6_js_.Menu({ commands });
        newMenu.title.label = trans.__('New');
        newMenu.title.icon = ui_components_lib_index_js_.caretDownIcon;
        menubar.addMenu(newMenu);
        const populateNewMenu = () => {
            var _a, _b;
            // create an entry per kernel spec for creating a new notebook
            const specs = (_b = (_a = serviceManager.kernelspecs) === null || _a === void 0 ? void 0 : _a.specs) === null || _b === void 0 ? void 0 : _b.kernelspecs;
            for (const name in specs) {
                newMenu.addItem({
                    args: { kernelName: name, isLauncher: true },
                    command: 'notebook:create-new',
                });
            }
            const baseCommands = [
                'terminal:create-new',
                'console:create',
                'filebrowser:create-new-file',
                'filebrowser:create-new-directory',
            ];
            baseCommands.forEach((command) => {
                newMenu.addItem({ command });
            });
        };
        (_a = serviceManager.kernelspecs) === null || _a === void 0 ? void 0 : _a.specsChanged.connect(() => {
            newMenu.clearItems();
            populateNewMenu();
        });
        populateNewMenu();
        if (toolbarRegistry) {
            toolbarRegistry.addFactory(FILE_BROWSER_FACTORY, 'new-dropdown', (browser) => {
                const menubar = new dist_index_es6_js_.MenuBar(overflowOptions);
                menubar.addMenu(newMenu);
                menubar.addClass('jp-DropdownMenu');
                return menubar;
            });
        }
    },
};
/**
 * A plugin to add file browser actions to the file browser toolbar.
 */
const fileActions = {
    id: '@jupyter-notebook/tree-extension:file-actions',
    description: 'A plugin to add file browser actions to the file browser toolbar.',
    autoStart: true,
    requires: [filebrowser_lib_index_js_.IDefaultFileBrowser, index_js_.IToolbarWidgetRegistry, translation_lib_index_js_.ITranslator],
    activate: (app, browser, toolbarRegistry, translator) => {
        // TODO: use upstream signal when available to detect selection changes
        // https://github.com/jupyterlab/jupyterlab/issues/14598
        const selectionChanged = new index_es6_js_.Signal(browser);
        const methods = [
            '_selectItem',
            '_handleMultiSelect',
            'handleFileSelect',
        ];
        methods.forEach((method) => {
            const original = browser['listing'][method];
            browser['listing'][method] = (...args) => {
                original.call(browser['listing'], ...args);
                selectionChanged.emit(void 0);
            };
        });
        browser.model.pathChanged.connect(() => {
            selectionChanged.emit(void 0);
        });
        // Create a toolbar item that adds buttons to the file browser toolbar
        // to perform actions on the files
        const { commands } = app;
        const fileActions = new FilesActionButtons({
            commands,
            browser,
            selectionChanged,
            translator,
        });
        for (const widget of fileActions.widgets) {
            toolbarRegistry.addFactory(FILE_BROWSER_FACTORY, widget.id, () => widget);
        }
    },
};
/**
 * A plugin to set the default file browser settings.
 */
const fileBrowserSettings = {
    id: '@jupyter-notebook/tree-extension:settings',
    description: 'Set up the default file browser settings',
    requires: [filebrowser_lib_index_js_.IDefaultFileBrowser],
    optional: [settingregistry_lib_index_js_.ISettingRegistry],
    autoStart: true,
    activate: (app, browser, settingRegistry) => {
        // Default config for notebook.
        // This is a different set of defaults than JupyterLab.
        const defaultFileBrowserConfig = {
            navigateToCurrentDirectory: false,
            singleClickNavigation: true,
            showLastModifiedColumn: true,
            showFileSizeColumn: true,
            showHiddenFiles: false,
            showFileCheckboxes: true,
            sortNotebooksFirst: true,
            showFullPath: false,
        };
        // Apply defaults on plugin activation
        let key;
        for (key in defaultFileBrowserConfig) {
            browser[key] = defaultFileBrowserConfig[key];
        }
        if (settingRegistry) {
            void settingRegistry.load(FILE_BROWSER_PLUGIN_ID).then((settings) => {
                function onSettingsChanged(settings) {
                    let key;
                    for (key in defaultFileBrowserConfig) {
                        const value = settings.get(key).user;
                        // only set the setting if it is defined by the user
                        if (value !== undefined) {
                            browser[key] = value;
                        }
                    }
                }
                settings.changed.connect(onSettingsChanged);
                onSettingsChanged(settings);
            });
        }
    },
};
/**
 * A plugin to add the file filter toggle command to the palette
 */
const fileFilterCommand = {
    id: '@jupyter-notebook/tree-extension:file-filter-command',
    description: 'A plugin to add file filter command to the palette.',
    autoStart: true,
    optional: [index_js_.ICommandPalette],
    activate: (app, palette) => {
        if (palette) {
            palette.addItem({
                command: CommandIDs.toggleFileFilter,
                category: 'File Browser',
            });
        }
    },
};
/**
 * Plugin to load the default plugins that are loaded on all the Notebook pages
 * (tree, edit, view, etc.) so they are visible in the settings editor.
 */
const loadPlugins = {
    id: '@jupyter-notebook/tree-extension:load-plugins',
    description: 'Plugin to load the default plugins that are loaded on all the Notebook pages (tree, edit, view, etc.) so they are visible in the settings editor.',
    autoStart: true,
    requires: [settingregistry_lib_index_js_.ISettingRegistry],
    activate(app, settingRegistry) {
        const { isDisabled } = lib_index_js_.PageConfig.Extension;
        const connector = settingRegistry.connector;
        const allPluginsOption = lib_index_js_.PageConfig.getOption('allPlugins');
        if (!allPluginsOption) {
            return;
        }
        // build the list of plugins shipped by default on the all the notebook pages
        // this avoid explicitly loading `'all'` plugins such as the ones used
        // in JupyterLab only
        const allPlugins = JSON.parse(allPluginsOption);
        const pluginsSet = new Set();
        Object.keys(allPlugins).forEach((key) => {
            const extensionsAndPlugins = allPlugins[key];
            Object.keys(extensionsAndPlugins).forEach((plugin) => {
                const value = extensionsAndPlugins[plugin];
                if (typeof value === 'boolean' && value) {
                    pluginsSet.add(plugin);
                }
                else if (Array.isArray(value)) {
                    value.forEach((v) => {
                        pluginsSet.add(v);
                    });
                }
            });
        });
        app.restored.then(async () => {
            const plugins = await connector.list('all');
            plugins.ids.forEach(async (id) => {
                const [extension] = id.split(':');
                // load the plugin if it is built-in the notebook application explicitly
                // either included as an extension or as a plugin directly
                const hasPlugin = pluginsSet.has(extension) || pluginsSet.has(id);
                if (!hasPlugin || isDisabled(id) || id in settingRegistry.plugins) {
                    return;
                }
                try {
                    await settingRegistry.load(id);
                }
                catch (error) {
                    console.warn(`Settings failed to load for (${id})`, error);
                }
            });
        });
    },
};
/**
 * A plugin to add file browser commands for the tree view.
 */
const openFileBrowser = {
    id: '@jupyter-notebook/tree-extension:open-file-browser',
    description: 'A plugin to add file browser commands for the tree view.',
    requires: [tree_lib_index_js_.INotebookTree, filebrowser_lib_index_js_.IDefaultFileBrowser],
    autoStart: true,
    activate: (app, notebookTree, browser) => {
        const { commands } = app;
        commands.addCommand(CommandIDs.activate, {
            execute: () => {
                notebookTree.currentWidget = browser;
            },
        });
    },
};
/**
 * A plugin to add the file browser widget to an INotebookShell
 */
const notebookTreeWidget = {
    id: '@jupyter-notebook/tree-extension:widget',
    description: 'A plugin to add the file browser widget to an INotebookShell.',
    requires: [
        filebrowser_lib_index_js_.IDefaultFileBrowser,
        translation_lib_index_js_.ITranslator,
        settingregistry_lib_index_js_.ISettingRegistry,
        index_js_.IToolbarWidgetRegistry,
        filebrowser_lib_index_js_.IFileBrowserFactory,
    ],
    optional: [
        running_lib_index_js_.IRunningSessionManagers,
        settingeditor_lib_index_js_.ISettingEditorTracker,
        settingeditor_lib_index_js_.IJSONSettingEditorTracker,
    ],
    autoStart: true,
    provides: tree_lib_index_js_.INotebookTree,
    activate: (app, browser, translator, settingRegistry, toolbarRegistry, factory, manager, settingEditorTracker, jsonSettingEditorTracker) => {
        const nbTreeWidget = new tree_lib_index_js_.NotebookTreeWidget();
        const trans = translator.load('notebook');
        browser.title.label = trans.__('Files');
        browser.node.setAttribute('role', 'region');
        browser.node.setAttribute('aria-label', trans.__('File Browser Section'));
        browser.title.icon = ui_components_lib_index_js_.folderIcon;
        nbTreeWidget.addWidget(browser);
        nbTreeWidget.tabBar.addTab(browser.title);
        nbTreeWidget.tabsMovable = false;
        toolbarRegistry.addFactory(FILE_BROWSER_FACTORY, 'uploader', (browser) => new filebrowser_lib_index_js_.Uploader({
            model: browser.model,
            translator,
            label: trans.__('Upload'),
        }));
        (0,index_js_.setToolbar)(browser, (0,index_js_.createToolbarFactory)(toolbarRegistry, settingRegistry, FILE_BROWSER_FACTORY, notebookTreeWidget.id, translator));
        if (manager) {
            const running = new running_lib_index_js_.RunningSessions(manager, translator);
            running.id = 'jp-running-sessions-tree';
            running.title.label = trans.__('Running');
            running.title.icon = ui_components_lib_index_js_.runningIcon;
            nbTreeWidget.addWidget(running);
            nbTreeWidget.tabBar.addTab(running.title);
        }
        app.shell.add(nbTreeWidget, 'main', { rank: 100 });
        // add a separate tab for each setting editor
        [settingEditorTracker, jsonSettingEditorTracker].forEach((editorTracker) => {
            if (editorTracker) {
                editorTracker.widgetAdded.connect((_, editor) => {
                    nbTreeWidget.addWidget(editor);
                    nbTreeWidget.tabBar.addTab(editor.title);
                    nbTreeWidget.currentWidget = editor;
                });
            }
        });
        const { tracker } = factory;
        // TODO: remove
        // Workaround to force the focus on the default file browser
        // See https://github.com/jupyterlab/jupyterlab/issues/15629 for more info
        const setCurrentToDefaultBrower = () => {
            tracker['_pool'].current = browser;
        };
        tracker.widgetAdded.connect((sender, widget) => {
            setCurrentToDefaultBrower();
        });
        setCurrentToDefaultBrower();
        return nbTreeWidget;
    },
};
/**
 * Export the plugins as default.
 */
const plugins = [
    createNew,
    fileActions,
    fileBrowserSettings,
    fileFilterCommand,
    loadPlugins,
    openFileBrowser,
    notebookTreeWidget,
];
/* harmony default export */ const lib = (plugins);


/***/ })

}]);
//# sourceMappingURL=7302.8dfb32b083b16efa038a.js.map?v=8dfb32b083b16efa038a
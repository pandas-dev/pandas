"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[35],{

/***/ 50035:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

// ESM COMPAT FLAG
__webpack_require__.r(__webpack_exports__);

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  YBaseCell: () => (/* reexport */ YBaseCell),
  YCodeCell: () => (/* reexport */ YCodeCell),
  YDocument: () => (/* reexport */ YDocument),
  YFile: () => (/* reexport */ YFile),
  YMarkdownCell: () => (/* reexport */ YMarkdownCell),
  YNotebook: () => (/* reexport */ YNotebook),
  YRawCell: () => (/* reexport */ YRawCell),
  convertYMapEventToMapChange: () => (/* reexport */ convertYMapEventToMapChange),
  createMutex: () => (/* reexport */ createMutex),
  createStandaloneCell: () => (/* reexport */ createStandaloneCell)
});

;// CONCATENATED MODULE: ../node_modules/@jupyter/ydoc/lib/utils.js
/* -----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/
function convertYMapEventToMapChange(event) {
    let changes = new Map();
    event.changes.keys.forEach((event, key) => {
        changes.set(key, {
            action: event.action,
            oldValue: event.oldValue,
            newValue: this.ymeta.get(key)
        });
    });
    return changes;
}
/**
 * Creates a mutual exclude function with the following property:
 *
 * ```js
 * const mutex = createMutex()
 * mutex(() => {
 *   // This function is immediately executed
 *   mutex(() => {
 *     // This function is not executed, as the mutex is already active.
 *   })
 * })
 * ```
 */
const createMutex = () => {
    let token = true;
    return (f) => {
        if (token) {
            token = false;
            try {
                f();
            }
            finally {
                token = true;
            }
        }
    };
};
//# sourceMappingURL=utils.js.map
// EXTERNAL MODULE: consume shared module (default) @lumino/coreutils@~2.2.2 (singleton) (fallback: ../node_modules/@lumino/coreutils/dist/index.js)
var index_js_ = __webpack_require__(72215);
// EXTERNAL MODULE: consume shared module (default) @lumino/signaling@~2.1.5 (singleton) (fallback: ../node_modules/@lumino/signaling/dist/index.es6.js)
var index_es6_js_ = __webpack_require__(46257);
// EXTERNAL MODULE: ../node_modules/lib0/time.js
var lib0_time = __webpack_require__(2431);
// EXTERNAL MODULE: ../node_modules/lib0/math.js
var math = __webpack_require__(11182);
// EXTERNAL MODULE: ../node_modules/lib0/observable.js
var observable = __webpack_require__(12330);
// EXTERNAL MODULE: ../node_modules/lib0/function.js
var lib0_function = __webpack_require__(38828);
// EXTERNAL MODULE: consume shared module (default) yjs@~13.6.8 (singleton) (fallback: ../node_modules/yjs/dist/yjs.mjs)
var yjs_mjs_ = __webpack_require__(17843);
;// CONCATENATED MODULE: ../node_modules/y-protocols/awareness.js
/**
 * @module awareness-protocol
 */







 // eslint-disable-line

const outdatedTimeout = 30000

/**
 * @typedef {Object} MetaClientState
 * @property {number} MetaClientState.clock
 * @property {number} MetaClientState.lastUpdated unix timestamp
 */

/**
 * The Awareness class implements a simple shared state protocol that can be used for non-persistent data like awareness information
 * (cursor, username, status, ..). Each client can update its own local state and listen to state changes of
 * remote clients. Every client may set a state of a remote peer to `null` to mark the client as offline.
 *
 * Each client is identified by a unique client id (something we borrow from `doc.clientID`). A client can override
 * its own state by propagating a message with an increasing timestamp (`clock`). If such a message is received, it is
 * applied if the known state of that client is older than the new state (`clock < newClock`). If a client thinks that
 * a remote client is offline, it may propagate a message with
 * `{ clock: currentClientClock, state: null, client: remoteClient }`. If such a
 * message is received, and the known clock of that client equals the received clock, it will override the state with `null`.
 *
 * Before a client disconnects, it should propagate a `null` state with an updated clock.
 *
 * Awareness states must be updated every 30 seconds. Otherwise the Awareness instance will delete the client state.
 *
 * @extends {Observable<string>}
 */
class Awareness extends observable/* Observable */.y {
  /**
   * @param {Y.Doc} doc
   */
  constructor (doc) {
    super()
    this.doc = doc
    /**
     * @type {number}
     */
    this.clientID = doc.clientID
    /**
     * Maps from client id to client state
     * @type {Map<number, Object<string, any>>}
     */
    this.states = new Map()
    /**
     * @type {Map<number, MetaClientState>}
     */
    this.meta = new Map()
    this._checkInterval = /** @type {any} */ (setInterval(() => {
      const now = lib0_time/* getUnixTime */.ZG()
      if (this.getLocalState() !== null && (outdatedTimeout / 2 <= now - /** @type {{lastUpdated:number}} */ (this.meta.get(this.clientID)).lastUpdated)) {
        // renew local clock
        this.setLocalState(this.getLocalState())
      }
      /**
       * @type {Array<number>}
       */
      const remove = []
      this.meta.forEach((meta, clientid) => {
        if (clientid !== this.clientID && outdatedTimeout <= now - meta.lastUpdated && this.states.has(clientid)) {
          remove.push(clientid)
        }
      })
      if (remove.length > 0) {
        removeAwarenessStates(this, remove, 'timeout')
      }
    }, math/* floor */.GW(outdatedTimeout / 10)))
    doc.on('destroy', () => {
      this.destroy()
    })
    this.setLocalState({})
  }

  destroy () {
    this.emit('destroy', [this])
    this.setLocalState(null)
    super.destroy()
    clearInterval(this._checkInterval)
  }

  /**
   * @return {Object<string,any>|null}
   */
  getLocalState () {
    return this.states.get(this.clientID) || null
  }

  /**
   * @param {Object<string,any>|null} state
   */
  setLocalState (state) {
    const clientID = this.clientID
    const currLocalMeta = this.meta.get(clientID)
    const clock = currLocalMeta === undefined ? 0 : currLocalMeta.clock + 1
    const prevState = this.states.get(clientID)
    if (state === null) {
      this.states.delete(clientID)
    } else {
      this.states.set(clientID, state)
    }
    this.meta.set(clientID, {
      clock,
      lastUpdated: lib0_time/* getUnixTime */.ZG()
    })
    const added = []
    const updated = []
    const filteredUpdated = []
    const removed = []
    if (state === null) {
      removed.push(clientID)
    } else if (prevState == null) {
      if (state != null) {
        added.push(clientID)
      }
    } else {
      updated.push(clientID)
      if (!lib0_function/* equalityDeep */.Hi(prevState, state)) {
        filteredUpdated.push(clientID)
      }
    }
    if (added.length > 0 || filteredUpdated.length > 0 || removed.length > 0) {
      this.emit('change', [{ added, updated: filteredUpdated, removed }, 'local'])
    }
    this.emit('update', [{ added, updated, removed }, 'local'])
  }

  /**
   * @param {string} field
   * @param {any} value
   */
  setLocalStateField (field, value) {
    const state = this.getLocalState()
    if (state !== null) {
      this.setLocalState({
        ...state,
        [field]: value
      })
    }
  }

  /**
   * @return {Map<number,Object<string,any>>}
   */
  getStates () {
    return this.states
  }
}

/**
 * Mark (remote) clients as inactive and remove them from the list of active peers.
 * This change will be propagated to remote clients.
 *
 * @param {Awareness} awareness
 * @param {Array<number>} clients
 * @param {any} origin
 */
const removeAwarenessStates = (awareness, clients, origin) => {
  const removed = []
  for (let i = 0; i < clients.length; i++) {
    const clientID = clients[i]
    if (awareness.states.has(clientID)) {
      awareness.states.delete(clientID)
      if (clientID === awareness.clientID) {
        const curMeta = /** @type {MetaClientState} */ (awareness.meta.get(clientID))
        awareness.meta.set(clientID, {
          clock: curMeta.clock + 1,
          lastUpdated: lib0_time/* getUnixTime */.ZG()
        })
      }
      removed.push(clientID)
    }
  }
  if (removed.length > 0) {
    awareness.emit('change', [{ added: [], updated: [], removed }, origin])
    awareness.emit('update', [{ added: [], updated: [], removed }, origin])
  }
}

/**
 * @param {Awareness} awareness
 * @param {Array<number>} clients
 * @return {Uint8Array}
 */
const encodeAwarenessUpdate = (awareness, clients, states = awareness.states) => {
  const len = clients.length
  const encoder = encoding.createEncoder()
  encoding.writeVarUint(encoder, len)
  for (let i = 0; i < len; i++) {
    const clientID = clients[i]
    const state = states.get(clientID) || null
    const clock = /** @type {MetaClientState} */ (awareness.meta.get(clientID)).clock
    encoding.writeVarUint(encoder, clientID)
    encoding.writeVarUint(encoder, clock)
    encoding.writeVarString(encoder, JSON.stringify(state))
  }
  return encoding.toUint8Array(encoder)
}

/**
 * Modify the content of an awareness update before re-encoding it to an awareness update.
 *
 * This might be useful when you have a central server that wants to ensure that clients
 * cant hijack somebody elses identity.
 *
 * @param {Uint8Array} update
 * @param {function(any):any} modify
 * @return {Uint8Array}
 */
const modifyAwarenessUpdate = (update, modify) => {
  const decoder = decoding.createDecoder(update)
  const encoder = encoding.createEncoder()
  const len = decoding.readVarUint(decoder)
  encoding.writeVarUint(encoder, len)
  for (let i = 0; i < len; i++) {
    const clientID = decoding.readVarUint(decoder)
    const clock = decoding.readVarUint(decoder)
    const state = JSON.parse(decoding.readVarString(decoder))
    const modifiedState = modify(state)
    encoding.writeVarUint(encoder, clientID)
    encoding.writeVarUint(encoder, clock)
    encoding.writeVarString(encoder, JSON.stringify(modifiedState))
  }
  return encoding.toUint8Array(encoder)
}

/**
 * @param {Awareness} awareness
 * @param {Uint8Array} update
 * @param {any} origin This will be added to the emitted change event
 */
const applyAwarenessUpdate = (awareness, update, origin) => {
  const decoder = decoding.createDecoder(update)
  const timestamp = time.getUnixTime()
  const added = []
  const updated = []
  const filteredUpdated = []
  const removed = []
  const len = decoding.readVarUint(decoder)
  for (let i = 0; i < len; i++) {
    const clientID = decoding.readVarUint(decoder)
    let clock = decoding.readVarUint(decoder)
    const state = JSON.parse(decoding.readVarString(decoder))
    const clientMeta = awareness.meta.get(clientID)
    const prevState = awareness.states.get(clientID)
    const currClock = clientMeta === undefined ? 0 : clientMeta.clock
    if (currClock < clock || (currClock === clock && state === null && awareness.states.has(clientID))) {
      if (state === null) {
        // never let a remote client remove this local state
        if (clientID === awareness.clientID && awareness.getLocalState() != null) {
          // remote client removed the local state. Do not remote state. Broadcast a message indicating
          // that this client still exists by increasing the clock
          clock++
        } else {
          awareness.states.delete(clientID)
        }
      } else {
        awareness.states.set(clientID, state)
      }
      awareness.meta.set(clientID, {
        clock,
        lastUpdated: timestamp
      })
      if (clientMeta === undefined && state !== null) {
        added.push(clientID)
      } else if (clientMeta !== undefined && state === null) {
        removed.push(clientID)
      } else if (state !== null) {
        if (!f.equalityDeep(state, prevState)) {
          filteredUpdated.push(clientID)
        }
        updated.push(clientID)
      }
    }
  }
  if (added.length > 0 || filteredUpdated.length > 0 || removed.length > 0) {
    awareness.emit('change', [{
      added, updated: filteredUpdated, removed
    }, origin])
  }
  if (added.length > 0 || updated.length > 0 || removed.length > 0) {
    awareness.emit('update', [{
      added, updated, removed
    }, origin])
  }
}

;// CONCATENATED MODULE: ../node_modules/@jupyter/ydoc/lib/ydocument.js
/* -----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/




/**
 * Generic shareable document.
 */
class YDocument {
    constructor(options) {
        var _a;
        /**
         * Handle a change to the ystate.
         */
        this.onStateChanged = (event) => {
            const stateChange = new Array();
            event.keysChanged.forEach(key => {
                const change = event.changes.keys.get(key);
                if (change) {
                    stateChange.push({
                        name: key,
                        oldValue: change.oldValue,
                        newValue: this.ystate.get(key)
                    });
                }
            });
            this._changed.emit({ stateChange });
        };
        this._changed = new index_es6_js_.Signal(this);
        this._isDisposed = false;
        this._disposed = new index_es6_js_.Signal(this);
        this._ydoc = (_a = options === null || options === void 0 ? void 0 : options.ydoc) !== null && _a !== void 0 ? _a : new yjs_mjs_.Doc();
        this._ystate = this._ydoc.getMap('state');
        this._undoManager = new yjs_mjs_.UndoManager([], {
            trackedOrigins: new Set([this]),
            doc: this._ydoc
        });
        this._awareness = new Awareness(this._ydoc);
        this._ystate.observe(this.onStateChanged);
    }
    /**
     * YJS document.
     */
    get ydoc() {
        return this._ydoc;
    }
    /**
     * Shared state
     */
    get ystate() {
        return this._ystate;
    }
    /**
     * YJS document undo manager
     */
    get undoManager() {
        return this._undoManager;
    }
    /**
     * Shared awareness
     */
    get awareness() {
        return this._awareness;
    }
    /**
     * The changed signal.
     */
    get changed() {
        return this._changed;
    }
    /**
     * A signal emitted when the document is disposed.
     */
    get disposed() {
        return this._disposed;
    }
    /**
     * Whether the document is disposed or not.
     */
    get isDisposed() {
        return this._isDisposed;
    }
    /**
     * Document state
     */
    get state() {
        return index_js_.JSONExt.deepCopy(this.ystate.toJSON());
    }
    /**
     * Whether the object can undo changes.
     */
    canUndo() {
        return this.undoManager.undoStack.length > 0;
    }
    /**
     * Whether the object can redo changes.
     */
    canRedo() {
        return this.undoManager.redoStack.length > 0;
    }
    /**
     * Dispose of the resources.
     */
    dispose() {
        if (this._isDisposed) {
            return;
        }
        this._isDisposed = true;
        this.ystate.unobserve(this.onStateChanged);
        this.awareness.destroy();
        this.undoManager.destroy();
        this.ydoc.destroy();
        this._disposed.emit();
        index_es6_js_.Signal.clearData(this);
    }
    /**
     * Get the value for a state attribute
     *
     * @param key Key to get
     */
    getState(key) {
        const value = this.ystate.get(key);
        return typeof value === 'undefined'
            ? value
            : index_js_.JSONExt.deepCopy(value);
    }
    /**
     * Set the value of a state attribute
     *
     * @param key Key to set
     * @param value New attribute value
     */
    setState(key, value) {
        if (!index_js_.JSONExt.deepEqual(this.ystate.get(key), value)) {
            this.ystate.set(key, value);
        }
    }
    /**
     * Get the document source
     *
     * @returns The source
     */
    get source() {
        return this.getSource();
    }
    /**
     * Set the document source
     *
     * @param value The source to set
     */
    set source(value) {
        this.setSource(value);
    }
    /**
     * Undo an operation.
     */
    undo() {
        this.undoManager.undo();
    }
    /**
     * Redo an operation.
     */
    redo() {
        this.undoManager.redo();
    }
    /**
     * Clear the change stack.
     */
    clearUndoHistory() {
        this.undoManager.clear();
    }
    /**
     * Perform a transaction. While the function f is called, all changes to the shared
     * document are bundled into a single event.
     */
    transact(f, undoable = true, origin = null) {
        this.ydoc.transact(f, undoable ? this : origin);
    }
}
//# sourceMappingURL=ydocument.js.map
;// CONCATENATED MODULE: ../node_modules/@jupyter/ydoc/lib/yfile.js
/* -----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/

/**
 * Shareable text file.
 */
class YFile extends YDocument {
    /**
     * Create a new file
     *
     * #### Notes
     * The document is empty and must be populated
     */
    constructor() {
        super();
        /**
         * Document version
         */
        this.version = '1.0.0';
        /**
         * YJS file text.
         */
        this.ysource = this.ydoc.getText('source');
        /**
         * Handle a change to the ymodel.
         */
        this._modelObserver = (event) => {
            this._changed.emit({ sourceChange: event.changes.delta });
        };
        this.undoManager.addToScope(this.ysource);
        this.ysource.observe(this._modelObserver);
    }
    /**
     * Creates a standalone YFile
     */
    static create() {
        return new YFile();
    }
    /**
     * File text
     */
    get source() {
        return this.getSource();
    }
    set source(v) {
        this.setSource(v);
    }
    /**
     * Dispose of the resources.
     */
    dispose() {
        if (this.isDisposed) {
            return;
        }
        this.ysource.unobserve(this._modelObserver);
        super.dispose();
    }
    /**
     * Get the file text.
     *
     * @returns File text.
     */
    getSource() {
        return this.ysource.toString();
    }
    /**
     * Set the file text.
     *
     * @param value New text
     */
    setSource(value) {
        this.transact(() => {
            const ytext = this.ysource;
            ytext.delete(0, ytext.length);
            ytext.insert(0, value);
        });
    }
    /**
     * Replace content from `start' to `end` with `value`.
     *
     * @param start: The start index of the range to replace (inclusive).
     * @param end: The end index of the range to replace (exclusive).
     * @param value: New source (optional).
     */
    updateSource(start, end, value = '') {
        this.transact(() => {
            const ysource = this.ysource;
            // insert and then delete.
            // This ensures that the cursor position is adjusted after the replaced content.
            ysource.insert(start, value);
            ysource.delete(start + value.length, end - start);
        });
    }
}
//# sourceMappingURL=yfile.js.map
;// CONCATENATED MODULE: ../node_modules/@jupyter/ydoc/lib/ycell.js
/* -----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/




/**
 * Create a new shared cell model given the YJS shared type.
 */
const createCellModelFromSharedType = (type, options = {}) => {
    switch (type.get('cell_type')) {
        case 'code':
            return new YCodeCell(type, type.get('source'), type.get('outputs'), options);
        case 'markdown':
            return new YMarkdownCell(type, type.get('source'), options);
        case 'raw':
            return new YRawCell(type, type.get('source'), options);
        default:
            throw new Error('Found unknown cell type');
    }
};
/**
 * Create a new cell that can be inserted in an existing shared model.
 *
 * If no notebook is specified the cell will be standalone.
 *
 * @param cell Cell JSON representation
 * @param notebook Notebook to which the cell will be added
 */
const createCell = (cell, notebook) => {
    var _a, _b;
    const ymodel = new yjs_mjs_.Map();
    const ysource = new yjs_mjs_.Text();
    const ymetadata = new yjs_mjs_.Map();
    ymodel.set('source', ysource);
    ymodel.set('metadata', ymetadata);
    ymodel.set('cell_type', cell.cell_type);
    ymodel.set('id', (_a = cell.id) !== null && _a !== void 0 ? _a : index_js_.UUID.uuid4());
    let ycell;
    switch (cell.cell_type) {
        case 'markdown': {
            ycell = new YMarkdownCell(ymodel, ysource, { notebook }, ymetadata);
            if (cell.attachments != null) {
                ycell.setAttachments(cell.attachments);
            }
            break;
        }
        case 'code': {
            const youtputs = new yjs_mjs_.Array();
            ymodel.set('outputs', youtputs);
            ycell = new YCodeCell(ymodel, ysource, youtputs, {
                notebook
            }, ymetadata);
            const cCell = cell;
            ycell.execution_count = (_b = cCell.execution_count) !== null && _b !== void 0 ? _b : null;
            if (cCell.outputs) {
                ycell.setOutputs(cCell.outputs);
            }
            break;
        }
        default: {
            // raw
            ycell = new YRawCell(ymodel, ysource, { notebook }, ymetadata);
            if (cell.attachments) {
                ycell.setAttachments(cell.attachments);
            }
            break;
        }
    }
    if (cell.metadata != null) {
        ycell.setMetadata(cell.metadata);
    }
    if (cell.source != null) {
        ycell.setSource(typeof cell.source === 'string' ? cell.source : cell.source.join(''));
    }
    return ycell;
};
/**
 * Create a new cell that cannot be inserted in an existing shared model.
 *
 * @param cell Cell JSON representation
 */
const createStandaloneCell = (cell) => createCell(cell);
class YBaseCell {
    /**
     * Create a new YCell that works standalone. It cannot be
     * inserted into a YNotebook because the Yjs model is already
     * attached to an anonymous Y.Doc instance.
     */
    static create(id) {
        return createCell({ id, cell_type: this.prototype.cell_type });
    }
    /**
     * Base cell constructor
     *
     * ### Notes
     * Don't use the constructor directly - prefer using ``YNotebook.insertCell``
     *
     * The ``ysource`` is needed because ``ymodel.get('source')`` will
     * not return the real source if the model is not yet attached to
     * a document. Requesting it explicitly allows to introspect a non-empty
     * source before the cell is attached to the document.
     *
     * @param ymodel Cell map
     * @param ysource Cell source
     * @param options \{ notebook?: The notebook the cell is attached to \}
     * @param ymetadata Cell metadata
     */
    constructor(ymodel, ysource, options = {}, ymetadata) {
        /**
         * Handle a change to the ymodel.
         */
        this._modelObserver = (events, transaction) => {
            if (transaction.origin !== 'silent-change') {
                this._changed.emit(this.getChanges(events));
            }
        };
        this._metadataChanged = new index_es6_js_.Signal(this);
        /**
         * The notebook that this cell belongs to.
         */
        this._notebook = null;
        this._changed = new index_es6_js_.Signal(this);
        this._disposed = new index_es6_js_.Signal(this);
        this._isDisposed = false;
        this._undoManager = null;
        this.ymodel = ymodel;
        this._ysource = ysource;
        this._ymetadata = ymetadata !== null && ymetadata !== void 0 ? ymetadata : this.ymodel.get('metadata');
        this._prevSourceLength = ysource ? ysource.length : 0;
        this._notebook = null;
        this._awareness = null;
        this._undoManager = null;
        if (options.notebook) {
            this._notebook = options.notebook;
            if (this._notebook.disableDocumentWideUndoRedo) {
                this._undoManager = new yjs_mjs_.UndoManager([this.ymodel], {
                    trackedOrigins: new Set([this]),
                    doc: this._notebook.ydoc
                });
            }
        }
        else {
            // Standalone cell
            const doc = new yjs_mjs_.Doc();
            doc.getArray().insert(0, [this.ymodel]);
            this._awareness = new Awareness(doc);
            this._undoManager = new yjs_mjs_.UndoManager([this.ymodel], {
                trackedOrigins: new Set([this])
            });
        }
        this.ymodel.observeDeep(this._modelObserver);
    }
    /**
     * Cell notebook awareness or null.
     */
    get awareness() {
        var _a, _b, _c;
        return (_c = (_a = this._awareness) !== null && _a !== void 0 ? _a : (_b = this.notebook) === null || _b === void 0 ? void 0 : _b.awareness) !== null && _c !== void 0 ? _c : null;
    }
    /**
     * The type of the cell.
     */
    get cell_type() {
        throw new Error('A YBaseCell must not be constructed');
    }
    /**
     * The changed signal.
     */
    get changed() {
        return this._changed;
    }
    /**
     * Signal emitted when the cell is disposed.
     */
    get disposed() {
        return this._disposed;
    }
    /**
     * Cell id
     */
    get id() {
        return this.getId();
    }
    /**
     * Whether the model has been disposed or not.
     */
    get isDisposed() {
        return this._isDisposed;
    }
    /**
     * Whether the cell is standalone or not.
     *
     * If the cell is standalone. It cannot be
     * inserted into a YNotebook because the Yjs model is already
     * attached to an anonymous Y.Doc instance.
     */
    get isStandalone() {
        return this._notebook !== null;
    }
    /**
     * Cell metadata.
     *
     * #### Notes
     * You should prefer to access and modify the specific key of interest.
     */
    get metadata() {
        return this.getMetadata();
    }
    set metadata(v) {
        this.setMetadata(v);
    }
    /**
     * Signal triggered when the cell metadata changes.
     */
    get metadataChanged() {
        return this._metadataChanged;
    }
    /**
     * The notebook that this cell belongs to.
     */
    get notebook() {
        return this._notebook;
    }
    /**
     * Cell input content.
     */
    get source() {
        return this.getSource();
    }
    set source(v) {
        this.setSource(v);
    }
    /**
     * The cell undo manager.
     */
    get undoManager() {
        var _a;
        if (!this.notebook) {
            return this._undoManager;
        }
        return ((_a = this.notebook) === null || _a === void 0 ? void 0 : _a.disableDocumentWideUndoRedo)
            ? this._undoManager
            : this.notebook.undoManager;
    }
    get ysource() {
        return this._ysource;
    }
    /**
     * Whether the object can undo changes.
     */
    canUndo() {
        return !!this.undoManager && this.undoManager.undoStack.length > 0;
    }
    /**
     * Whether the object can redo changes.
     */
    canRedo() {
        return !!this.undoManager && this.undoManager.redoStack.length > 0;
    }
    /**
     * Clear the change stack.
     */
    clearUndoHistory() {
        var _a;
        (_a = this.undoManager) === null || _a === void 0 ? void 0 : _a.clear();
    }
    /**
     * Undo an operation.
     */
    undo() {
        var _a;
        (_a = this.undoManager) === null || _a === void 0 ? void 0 : _a.undo();
    }
    /**
     * Redo an operation.
     */
    redo() {
        var _a;
        (_a = this.undoManager) === null || _a === void 0 ? void 0 : _a.redo();
    }
    /**
     * Dispose of the resources.
     */
    dispose() {
        var _a;
        if (this._isDisposed)
            return;
        this._isDisposed = true;
        this.ymodel.unobserveDeep(this._modelObserver);
        if (this._awareness) {
            // A new document is created for standalone cell.
            const doc = this._awareness.doc;
            this._awareness.destroy();
            doc.destroy();
        }
        if (this._undoManager) {
            // Be sure to not destroy the document undo manager.
            if (this._undoManager === ((_a = this.notebook) === null || _a === void 0 ? void 0 : _a.undoManager)) {
                this._undoManager = null;
            }
            else {
                this._undoManager.destroy();
            }
        }
        this._disposed.emit();
        index_es6_js_.Signal.clearData(this);
    }
    /**
     * Get cell id.
     *
     * @returns Cell id
     */
    getId() {
        return this.ymodel.get('id');
    }
    /**
     * Gets cell's source.
     *
     * @returns Cell's source.
     */
    getSource() {
        return this.ysource.toString();
    }
    /**
     * Sets cell's source.
     *
     * @param value: New source.
     */
    setSource(value) {
        this.transact(() => {
            this.ysource.delete(0, this.ysource.length);
            this.ysource.insert(0, value);
        });
        // @todo Do we need proper replace semantic? This leads to issues in editor bindings because they don't switch source.
        // this.ymodel.set('source', new Y.Text(value));
    }
    /**
     * Replace content from `start' to `end` with `value`.
     *
     * @param start: The start index of the range to replace (inclusive).
     *
     * @param end: The end index of the range to replace (exclusive).
     *
     * @param value: New source (optional).
     */
    updateSource(start, end, value = '') {
        this.transact(() => {
            const ysource = this.ysource;
            // insert and then delete.
            // This ensures that the cursor position is adjusted after the replaced content.
            ysource.insert(start, value);
            ysource.delete(start + value.length, end - start);
        });
    }
    /**
     * Delete a metadata cell.
     *
     * @param key The key to delete
     */
    deleteMetadata(key) {
        if (typeof this.getMetadata(key) === 'undefined') {
            return;
        }
        this.transact(() => {
            this._ymetadata.delete(key);
            const jupyter = this.getMetadata('jupyter');
            if (key === 'collapsed' && jupyter) {
                // eslint-disable-next-line @typescript-eslint/no-unused-vars
                const { outputs_hidden, ...others } = jupyter;
                if (Object.keys(others).length === 0) {
                    this._ymetadata.delete('jupyter');
                }
                else {
                    this._ymetadata.set('jupyter', others);
                }
            }
            else if (key === 'jupyter') {
                this._ymetadata.delete('collapsed');
            }
        }, false);
    }
    getMetadata(key) {
        const metadata = this._ymetadata;
        // Transiently the metadata can be missing - like during destruction
        if (metadata === undefined) {
            return undefined;
        }
        if (typeof key === 'string') {
            const value = metadata.get(key);
            return typeof value === 'undefined'
                ? undefined // undefined is converted to `{}` by `JSONExt.deepCopy`
                : index_js_.JSONExt.deepCopy(metadata.get(key));
        }
        else {
            return index_js_.JSONExt.deepCopy(metadata.toJSON());
        }
    }
    setMetadata(metadata, value) {
        var _a, _b;
        if (typeof metadata === 'string') {
            if (typeof value === 'undefined') {
                throw new TypeError(`Metadata value for ${metadata} cannot be 'undefined'; use deleteMetadata.`);
            }
            const key = metadata;
            // Only set metadata if we change something to avoid infinite
            // loop of signal changes.
            if (index_js_.JSONExt.deepEqual((_a = this.getMetadata(key)) !== null && _a !== void 0 ? _a : null, value)) {
                return;
            }
            this.transact(() => {
                var _a;
                this._ymetadata.set(key, value);
                if (key === 'collapsed') {
                    const jupyter = ((_a = this.getMetadata('jupyter')) !== null && _a !== void 0 ? _a : {});
                    if (jupyter.outputs_hidden !== value) {
                        this.setMetadata('jupyter', {
                            ...jupyter,
                            outputs_hidden: value
                        });
                    }
                }
                else if (key === 'jupyter') {
                    const isHidden = value['outputs_hidden'];
                    if (typeof isHidden !== 'undefined') {
                        if (this.getMetadata('collapsed') !== isHidden) {
                            this.setMetadata('collapsed', isHidden);
                        }
                    }
                    else {
                        this.deleteMetadata('collapsed');
                    }
                }
            }, false);
        }
        else {
            const clone = index_js_.JSONExt.deepCopy(metadata);
            if (clone.collapsed != null) {
                clone.jupyter = clone.jupyter || {};
                clone.jupyter.outputs_hidden = clone.collapsed;
            }
            else if (((_b = clone === null || clone === void 0 ? void 0 : clone.jupyter) === null || _b === void 0 ? void 0 : _b.outputs_hidden) != null) {
                clone.collapsed = clone.jupyter.outputs_hidden;
            }
            if (!index_js_.JSONExt.deepEqual(clone, this.getMetadata())) {
                this.transact(() => {
                    for (const [key, value] of Object.entries(clone)) {
                        this._ymetadata.set(key, value);
                    }
                }, false);
            }
        }
    }
    /**
     * Serialize the model to JSON.
     */
    toJSON() {
        return {
            id: this.getId(),
            cell_type: this.cell_type,
            source: this.getSource(),
            metadata: this.getMetadata()
        };
    }
    /**
     * Perform a transaction. While the function f is called, all changes to the shared
     * document are bundled into a single event.
     *
     * @param f Transaction to execute
     * @param undoable Whether to track the change in the action history or not (default `true`)
     */
    transact(f, undoable = true, origin = null) {
        !this.notebook || this.notebook.disableDocumentWideUndoRedo
            ? this.ymodel.doc == null
                ? f()
                : this.ymodel.doc.transact(f, undoable ? this : origin)
            : this.notebook.transact(f, undoable);
    }
    /**
     * Extract changes from YJS events
     *
     * @param events YJS events
     * @returns Cell changes
     */
    getChanges(events) {
        const changes = {};
        const sourceEvent = events.find(event => event.target === this.ymodel.get('source'));
        if (sourceEvent) {
            changes.sourceChange = sourceEvent.changes.delta;
        }
        const metadataEvents = events.find(event => event.target === this._ymetadata);
        if (metadataEvents) {
            changes.metadataChange = metadataEvents.changes.keys;
            metadataEvents.changes.keys.forEach((change, key) => {
                switch (change.action) {
                    case 'add':
                        this._metadataChanged.emit({
                            key,
                            newValue: this._ymetadata.get(key),
                            type: 'add'
                        });
                        break;
                    case 'delete':
                        this._metadataChanged.emit({
                            key,
                            oldValue: change.oldValue,
                            type: 'remove'
                        });
                        break;
                    case 'update':
                        {
                            const newValue = this._ymetadata.get(key);
                            const oldValue = change.oldValue;
                            let equal = true;
                            if (typeof oldValue == 'object' && typeof newValue == 'object') {
                                equal = index_js_.JSONExt.deepEqual(oldValue, newValue);
                            }
                            else {
                                equal = oldValue === newValue;
                            }
                            if (!equal) {
                                this._metadataChanged.emit({
                                    key,
                                    type: 'change',
                                    oldValue,
                                    newValue
                                });
                            }
                        }
                        break;
                }
            });
        }
        const modelEvent = events.find(event => event.target === this.ymodel);
        // The model allows us to replace the complete source with a new string. We express this in the Delta format
        // as a replace of the complete string.
        const ysource = this.ymodel.get('source');
        if (modelEvent && modelEvent.keysChanged.has('source')) {
            changes.sourceChange = [
                { delete: this._prevSourceLength },
                { insert: ysource.toString() }
            ];
        }
        this._prevSourceLength = ysource.length;
        return changes;
    }
}
/**
 * Shareable code cell.
 */
class YCodeCell extends YBaseCell {
    /**
     * Create a new YCodeCell that works standalone. It cannot be
     * inserted into a YNotebook because the Yjs model is already
     * attached to an anonymous Y.Doc instance.
     */
    static create(id) {
        return super.create(id);
    }
    /**
     * Code cell constructor
     *
     * ### Notes
     * Don't use the constructor directly - prefer using ``YNotebook.insertCell``
     *
     * The ``ysource`` is needed because ``ymodel.get('source')`` will
     * not return the real source if the model is not yet attached to
     * a document. Requesting it explicitly allows to introspect a non-empty
     * source before the cell is attached to the document.
     *
     * @param ymodel Cell map
     * @param ysource Cell source
     * @param youtputs Code cell outputs
     * @param options \{ notebook?: The notebook the cell is attached to \}
     * @param ymetadata Cell metadata
     */
    constructor(ymodel, ysource, youtputs, options = {}, ymetadata) {
        super(ymodel, ysource, options, ymetadata);
        this._youtputs = youtputs;
    }
    /**
     * The type of the cell.
     */
    get cell_type() {
        return 'code';
    }
    /**
     * The code cell's prompt number. Will be null if the cell has not been run.
     */
    get execution_count() {
        return this.ymodel.get('execution_count') || null;
    }
    set execution_count(count) {
        // Do not use `this.execution_count`. When initializing the
        // cell, we need to set execution_count to `null` if we compare
        // using `this.execution_count` it will return `null` and we will
        // never initialize it
        if (this.ymodel.get('execution_count') !== count) {
            this.transact(() => {
                this.ymodel.set('execution_count', count);
            }, false);
        }
    }
    /**
     * The code cell's execution state.
     */
    get executionState() {
        var _a;
        return (_a = this.ymodel.get('execution_state')) !== null && _a !== void 0 ? _a : 'idle';
    }
    set executionState(state) {
        if (this.ymodel.get('execution_state') !== state) {
            this.transact(() => {
                this.ymodel.set('execution_state', state);
            }, false);
        }
    }
    /**
     * Cell outputs.
     */
    get outputs() {
        return this.getOutputs();
    }
    set outputs(v) {
        this.setOutputs(v);
    }
    get youtputs() {
        return this._youtputs;
    }
    /**
     * Execution, display, or stream outputs.
     */
    getOutputs() {
        return index_js_.JSONExt.deepCopy(this._youtputs.toJSON());
    }
    createOutputs(outputs) {
        const newOutputs = [];
        for (const output of index_js_.JSONExt.deepCopy(outputs)) {
            let _newOutput1;
            if (output.output_type === 'stream') {
                // Set the text field as a Y.Text
                const { text, ...outputWithoutText } = output;
                _newOutput1 = outputWithoutText;
                const newText = new yjs_mjs_.Text();
                let _text = text instanceof Array ? text.join() : text;
                newText.insert(0, _text);
                _newOutput1['text'] = newText;
            }
            else {
                _newOutput1 = output;
            }
            const _newOutput2 = [];
            for (const [key, value] of Object.entries(_newOutput1)) {
                _newOutput2.push([key, value]);
            }
            const newOutput = new yjs_mjs_.Map(_newOutput2);
            newOutputs.push(newOutput);
        }
        return newOutputs;
    }
    /**
     * Replace all outputs.
     */
    setOutputs(outputs) {
        this.transact(() => {
            this._youtputs.delete(0, this._youtputs.length);
            const newOutputs = this.createOutputs(outputs);
            this._youtputs.insert(0, newOutputs);
        }, false);
    }
    /**
     * Remove text from a stream output.
     */
    removeStreamOutput(index, start, origin = null) {
        this.transact(() => {
            const output = this._youtputs.get(index);
            const prevText = output.get('text');
            const length = prevText.length - start;
            prevText.delete(start, length);
        }, false, origin);
    }
    /**
     * Append text to a stream output.
     */
    appendStreamOutput(index, text, origin = null) {
        this.transact(() => {
            const output = this._youtputs.get(index);
            const prevText = output.get('text');
            prevText.insert(prevText.length, text);
        }, false, origin);
    }
    /**
     * Replace content from `start' to `end` with `outputs`.
     *
     * @param start: The start index of the range to replace (inclusive).
     *
     * @param end: The end index of the range to replace (exclusive).
     *
     * @param outputs: New outputs (optional).
     */
    updateOutputs(start, end, outputs = [], origin = null) {
        const fin = end < this._youtputs.length ? end - start : this._youtputs.length - start;
        this.transact(() => {
            this._youtputs.delete(start, fin);
            const newOutputs = this.createOutputs(outputs);
            this._youtputs.insert(start, newOutputs);
        }, false, origin);
    }
    /**
     * Clear all outputs from the cell.
     */
    clearOutputs(origin = null) {
        this.transact(() => {
            this._youtputs.delete(0, this._youtputs.length);
        }, false, origin);
    }
    /**
     * Serialize the model to JSON.
     */
    toJSON() {
        return {
            ...super.toJSON(),
            outputs: this.getOutputs(),
            execution_count: this.execution_count
        };
    }
    /**
     * Extract changes from YJS events
     *
     * @param events YJS events
     * @returns Cell changes
     */
    getChanges(events) {
        const changes = super.getChanges(events);
        const streamOutputEvent = events.find(
        // Changes to the 'text' of a cell's stream output can be accessed like so:
        // ycell['outputs'][output_idx]['text']
        // This translates to an event path of: ['outputs', output_idx, 'text]
        event => event.path.length === 3 &&
            event.path[0] === 'outputs' &&
            event.path[2] === 'text');
        if (streamOutputEvent) {
            changes.streamOutputChange = streamOutputEvent.changes.delta;
        }
        const outputEvent = events.find(event => event.target === this.ymodel.get('outputs'));
        if (outputEvent) {
            changes.outputsChange = outputEvent.changes.delta;
        }
        const modelEvent = events.find(event => event.target === this.ymodel);
        if (modelEvent && modelEvent.keysChanged.has('execution_count')) {
            const change = modelEvent.changes.keys.get('execution_count');
            changes.executionCountChange = {
                oldValue: change.oldValue,
                newValue: this.ymodel.get('execution_count')
            };
        }
        if (modelEvent && modelEvent.keysChanged.has('execution_state')) {
            const change = modelEvent.changes.keys.get('execution_state');
            changes.executionStateChange = {
                oldValue: change.oldValue,
                newValue: this.ymodel.get('execution_state')
            };
        }
        return changes;
    }
}
class YAttachmentCell extends YBaseCell {
    /**
     * Cell attachments
     */
    get attachments() {
        return this.getAttachments();
    }
    set attachments(v) {
        this.setAttachments(v);
    }
    /**
     * Gets the cell attachments.
     *
     * @returns The cell attachments.
     */
    getAttachments() {
        return this.ymodel.get('attachments');
    }
    /**
     * Sets the cell attachments
     *
     * @param attachments: The cell attachments.
     */
    setAttachments(attachments) {
        this.transact(() => {
            if (attachments == null) {
                this.ymodel.delete('attachments');
            }
            else {
                this.ymodel.set('attachments', attachments);
            }
        }, false);
    }
    /**
     * Extract changes from YJS events
     *
     * @param events YJS events
     * @returns Cell changes
     */
    getChanges(events) {
        const changes = super.getChanges(events);
        const modelEvent = events.find(event => event.target === this.ymodel);
        if (modelEvent && modelEvent.keysChanged.has('attachments')) {
            const change = modelEvent.changes.keys.get('attachments');
            changes.attachmentsChange = {
                oldValue: change.oldValue,
                newValue: this.ymodel.get('attachments')
            };
        }
        return changes;
    }
}
/**
 * Shareable raw cell.
 */
class YRawCell extends YAttachmentCell {
    /**
     * Create a new YRawCell that works standalone. It cannot be
     * inserted into a YNotebook because the Yjs model is already
     * attached to an anonymous Y.Doc instance.
     */
    static create(id) {
        return super.create(id);
    }
    /**
     * String identifying the type of cell.
     */
    get cell_type() {
        return 'raw';
    }
    /**
     * Serialize the model to JSON.
     */
    toJSON() {
        return {
            id: this.getId(),
            cell_type: 'raw',
            source: this.getSource(),
            metadata: this.getMetadata(),
            attachments: this.getAttachments()
        };
    }
}
/**
 * Shareable markdown cell.
 */
class YMarkdownCell extends YAttachmentCell {
    /**
     * Create a new YMarkdownCell that works standalone. It cannot be
     * inserted into a YNotebook because the Yjs model is already
     * attached to an anonymous Y.Doc instance.
     */
    static create(id) {
        return super.create(id);
    }
    /**
     * String identifying the type of cell.
     */
    get cell_type() {
        return 'markdown';
    }
    /**
     * Serialize the model to JSON.
     */
    toJSON() {
        return {
            id: this.getId(),
            cell_type: 'markdown',
            source: this.getSource(),
            metadata: this.getMetadata(),
            attachments: this.getAttachments()
        };
    }
}
//# sourceMappingURL=ycell.js.map
;// CONCATENATED MODULE: ../node_modules/@jupyter/ydoc/lib/ynotebook.js
/* -----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/





/**
 * Shared implementation of the Shared Document types.
 *
 * Shared cells can be inserted into a SharedNotebook.
 * Shared cells only start emitting events when they are connected to a SharedNotebook.
 *
 * "Standalone" cells must not be inserted into a (Shared)Notebook.
 * Standalone cells emit events immediately after they have been created, but they must not
 * be included into a (Shared)Notebook.
 */
class YNotebook extends YDocument {
    /**
     * Create a new notebook
     *
     * #### Notes
     * The document is empty and must be populated
     *
     * @param options
     */
    constructor(options = {}) {
        var _a;
        super();
        /**
         * Document version
         */
        this.version = '2.0.0';
        /**
         * YJS map for the notebook metadata
         */
        this.ymeta = this.ydoc.getMap('meta');
        /**
         * Handle a change to the ystate.
         */
        this._onMetaChanged = (events) => {
            const metadataEvents = events.find(event => event.target === this.ymeta.get('metadata'));
            if (metadataEvents) {
                const metadataChange = metadataEvents.changes.keys;
                const ymetadata = this.ymeta.get('metadata');
                metadataEvents.changes.keys.forEach((change, key) => {
                    switch (change.action) {
                        case 'add':
                            this._metadataChanged.emit({
                                key,
                                type: 'add',
                                newValue: ymetadata.get(key)
                            });
                            break;
                        case 'delete':
                            this._metadataChanged.emit({
                                key,
                                type: 'remove',
                                oldValue: change.oldValue
                            });
                            break;
                        case 'update':
                            {
                                const newValue = ymetadata.get(key);
                                const oldValue = change.oldValue;
                                let equal = true;
                                if (typeof oldValue == 'object' && typeof newValue == 'object') {
                                    equal = index_js_.JSONExt.deepEqual(oldValue, newValue);
                                }
                                else {
                                    equal = oldValue === newValue;
                                }
                                if (!equal) {
                                    this._metadataChanged.emit({
                                        key,
                                        type: 'change',
                                        oldValue,
                                        newValue
                                    });
                                }
                            }
                            break;
                    }
                });
                this._changed.emit({ metadataChange });
            }
            const metaEvent = events.find(event => event.target === this.ymeta);
            if (!metaEvent) {
                return;
            }
            if (metaEvent.keysChanged.has('metadata')) {
                // Handle metadata change when adding/removing the YMap
                const change = metaEvent.changes.keys.get('metadata');
                if ((change === null || change === void 0 ? void 0 : change.action) === 'add' && !change.oldValue) {
                    const metadataChange = new Map();
                    for (const key of Object.keys(this.metadata)) {
                        metadataChange.set(key, {
                            action: 'add',
                            oldValue: undefined
                        });
                        this._metadataChanged.emit({
                            key,
                            type: 'add',
                            newValue: this.getMetadata(key)
                        });
                    }
                    this._changed.emit({ metadataChange });
                }
            }
            if (metaEvent.keysChanged.has('nbformat')) {
                const change = metaEvent.changes.keys.get('nbformat');
                const nbformatChanged = {
                    key: 'nbformat',
                    oldValue: (change === null || change === void 0 ? void 0 : change.oldValue) ? change.oldValue : undefined,
                    newValue: this.nbformat
                };
                this._changed.emit({ nbformatChanged });
            }
            if (metaEvent.keysChanged.has('nbformat_minor')) {
                const change = metaEvent.changes.keys.get('nbformat_minor');
                const nbformatChanged = {
                    key: 'nbformat_minor',
                    oldValue: (change === null || change === void 0 ? void 0 : change.oldValue) ? change.oldValue : undefined,
                    newValue: this.nbformat_minor
                };
                this._changed.emit({ nbformatChanged });
            }
        };
        /**
         * Handle a change to the list of cells.
         */
        this._onYCellsChanged = (event) => {
            // update the type cell mapping by iterating through the added/removed types
            event.changes.added.forEach(item => {
                const type = item.content.type;
                if (!this._ycellMapping.has(type)) {
                    const c = createCellModelFromSharedType(type, { notebook: this });
                    this._ycellMapping.set(type, c);
                }
            });
            event.changes.deleted.forEach(item => {
                const type = item.content.type;
                const model = this._ycellMapping.get(type);
                if (model) {
                    model.dispose();
                    this._ycellMapping.delete(type);
                }
            });
            let index = 0;
            // this reflects the event.changes.delta, but replaces the content of delta.insert with ycells
            const cellsChange = [];
            event.changes.delta.forEach((d) => {
                if (d.insert != null) {
                    const insertedCells = d.insert.map((ycell) => this._ycellMapping.get(ycell));
                    cellsChange.push({ insert: insertedCells });
                    this.cells.splice(index, 0, ...insertedCells);
                    index += d.insert.length;
                }
                else if (d.delete != null) {
                    cellsChange.push(d);
                    this.cells.splice(index, d.delete);
                }
                else if (d.retain != null) {
                    cellsChange.push(d);
                    index += d.retain;
                }
            });
            this._changed.emit({
                cellsChange: cellsChange
            });
        };
        this._metadataChanged = new index_es6_js_.Signal(this);
        /**
         * Internal Yjs cells list
         */
        this._ycells = this.ydoc.getArray('cells');
        this._ycellMapping = new WeakMap();
        this._disableDocumentWideUndoRedo =
            (_a = options.disableDocumentWideUndoRedo) !== null && _a !== void 0 ? _a : false;
        this.cells = this._ycells.toArray().map(ycell => {
            if (!this._ycellMapping.has(ycell)) {
                this._ycellMapping.set(ycell, createCellModelFromSharedType(ycell, { notebook: this }));
            }
            return this._ycellMapping.get(ycell);
        });
        this.undoManager.addToScope(this._ycells);
        this._ycells.observe(this._onYCellsChanged);
        this.ymeta.observeDeep(this._onMetaChanged);
    }
    /**
     * Creates a standalone YNotebook
     *
     * Note: This method is useful when we need to initialize
     * the YNotebook from the JavaScript side.
     */
    static create(options = {}) {
        var _a, _b, _c, _d, _e, _f, _g, _h, _j;
        const ynotebook = new YNotebook({
            disableDocumentWideUndoRedo: (_a = options.disableDocumentWideUndoRedo) !== null && _a !== void 0 ? _a : false
        });
        const data = {
            cells: (_c = (_b = options.data) === null || _b === void 0 ? void 0 : _b.cells) !== null && _c !== void 0 ? _c : [],
            nbformat: (_e = (_d = options.data) === null || _d === void 0 ? void 0 : _d.nbformat) !== null && _e !== void 0 ? _e : 4,
            nbformat_minor: (_g = (_f = options.data) === null || _f === void 0 ? void 0 : _f.nbformat_minor) !== null && _g !== void 0 ? _g : 5,
            metadata: (_j = (_h = options.data) === null || _h === void 0 ? void 0 : _h.metadata) !== null && _j !== void 0 ? _j : {}
        };
        ynotebook.fromJSON(data);
        return ynotebook;
    }
    /**
     * Wether the undo/redo logic should be
     * considered on the full document across all cells.
     *
     * Default: false
     */
    get disableDocumentWideUndoRedo() {
        return this._disableDocumentWideUndoRedo;
    }
    /**
     * Notebook metadata
     */
    get metadata() {
        return this.getMetadata();
    }
    set metadata(v) {
        this.setMetadata(v);
    }
    /**
     * Signal triggered when a metadata changes.
     */
    get metadataChanged() {
        return this._metadataChanged;
    }
    /**
     * nbformat major version
     */
    get nbformat() {
        return this.ymeta.get('nbformat');
    }
    set nbformat(value) {
        this.transact(() => {
            this.ymeta.set('nbformat', value);
        }, false);
    }
    /**
     * nbformat minor version
     */
    get nbformat_minor() {
        return this.ymeta.get('nbformat_minor');
    }
    set nbformat_minor(value) {
        this.transact(() => {
            this.ymeta.set('nbformat_minor', value);
        }, false);
    }
    /**
     * Dispose of the resources.
     */
    dispose() {
        if (this.isDisposed) {
            return;
        }
        this._ycells.unobserve(this._onYCellsChanged);
        this.ymeta.unobserveDeep(this._onMetaChanged);
        super.dispose();
    }
    /**
     * Get a shared cell by index.
     *
     * @param index: Cell's position.
     *
     * @returns The requested shared cell.
     */
    getCell(index) {
        return this.cells[index];
    }
    /**
     * Add a shared cell at the notebook bottom.
     *
     * @param cell Cell to add.
     *
     * @returns The added cell.
     */
    addCell(cell) {
        return this.insertCell(this._ycells.length, cell);
    }
    /**
     * Insert a shared cell into a specific position.
     *
     * @param index: Cell's position.
     * @param cell: Cell to insert.
     *
     * @returns The inserted cell.
     */
    insertCell(index, cell) {
        return this.insertCells(index, [cell])[0];
    }
    /**
     * Insert a list of shared cells into a specific position.
     *
     * @param index: Position to insert the cells.
     * @param cells: Array of shared cells to insert.
     *
     * @returns The inserted cells.
     */
    insertCells(index, cells) {
        const yCells = cells.map(c => {
            const cell = createCell(c, this);
            this._ycellMapping.set(cell.ymodel, cell);
            return cell;
        });
        this.transact(() => {
            this._ycells.insert(index, yCells.map(cell => cell.ymodel));
        });
        return yCells;
    }
    /**
     * Move a cell.
     *
     * @param fromIndex: Index of the cell to move.
     * @param toIndex: New position of the cell.
     */
    moveCell(fromIndex, toIndex) {
        this.moveCells(fromIndex, toIndex);
    }
    /**
     * Move cells.
     *
     * @param fromIndex: Index of the first cells to move.
     * @param toIndex: New position of the first cell (in the current array).
     * @param n: Number of cells to move (default 1)
     */
    moveCells(fromIndex, toIndex, n = 1) {
        // FIXME we need to use yjs move feature to preserve undo history
        const clones = new Array(n)
            .fill(true)
            .map((_, idx) => this.getCell(fromIndex + idx).toJSON());
        this.transact(() => {
            this._ycells.delete(fromIndex, n);
            this._ycells.insert(fromIndex > toIndex ? toIndex : toIndex - n + 1, clones.map(clone => createCell(clone, this).ymodel));
        });
    }
    /**
     * Remove a cell.
     *
     * @param index: Index of the cell to remove.
     */
    deleteCell(index) {
        this.deleteCellRange(index, index + 1);
    }
    /**
     * Remove a range of cells.
     *
     * @param from: The start index of the range to remove (inclusive).
     * @param to: The end index of the range to remove (exclusive).
     */
    deleteCellRange(from, to) {
        // Cells will be removed from the mapping in the model event listener.
        this.transact(() => {
            this._ycells.delete(from, to - from);
        });
    }
    /**
     * Delete a metadata notebook.
     *
     * @param key The key to delete
     */
    deleteMetadata(key) {
        if (typeof this.getMetadata(key) === 'undefined') {
            return;
        }
        const allMetadata = this.metadata;
        delete allMetadata[key];
        this.setMetadata(allMetadata);
    }
    getMetadata(key) {
        const ymetadata = this.ymeta.get('metadata');
        // Transiently the metadata can be missing - like during destruction
        if (ymetadata === undefined) {
            return undefined;
        }
        if (typeof key === 'string') {
            const value = ymetadata.get(key);
            return typeof value === 'undefined'
                ? undefined // undefined is converted to `{}` by `JSONExt.deepCopy`
                : index_js_.JSONExt.deepCopy(value);
        }
        else {
            return index_js_.JSONExt.deepCopy(ymetadata.toJSON());
        }
    }
    setMetadata(metadata, value) {
        var _a;
        if (typeof metadata === 'string') {
            if (typeof value === 'undefined') {
                throw new TypeError(`Metadata value for ${metadata} cannot be 'undefined'; use deleteMetadata.`);
            }
            if (index_js_.JSONExt.deepEqual((_a = this.getMetadata(metadata)) !== null && _a !== void 0 ? _a : null, value)) {
                return;
            }
            const update = {};
            update[metadata] = value;
            this.updateMetadata(update);
        }
        else {
            if (!this.metadata || !index_js_.JSONExt.deepEqual(this.metadata, metadata)) {
                const clone = index_js_.JSONExt.deepCopy(metadata);
                const ymetadata = this.ymeta.get('metadata');
                // Transiently the metadata can be missing - like during destruction
                if (ymetadata === undefined) {
                    return undefined;
                }
                this.transact(() => {
                    ymetadata.clear();
                    for (const [key, value] of Object.entries(clone)) {
                        ymetadata.set(key, value);
                    }
                });
            }
        }
    }
    /**
     * Updates the metadata associated with the notebook.
     *
     * @param value: Metadata's attribute to update.
     */
    updateMetadata(value) {
        // TODO: Maybe modify only attributes instead of replacing the whole metadata?
        const clone = index_js_.JSONExt.deepCopy(value);
        const ymetadata = this.ymeta.get('metadata');
        // Transiently the metadata can be missing - like during destruction
        if (ymetadata === undefined) {
            return undefined;
        }
        this.transact(() => {
            for (const [key, value] of Object.entries(clone)) {
                ymetadata.set(key, value);
            }
        });
    }
    /**
     * Get the notebook source
     *
     * @returns The notebook
     */
    getSource() {
        return this.toJSON();
    }
    /**
     * Set the notebook source
     *
     * @param value The notebook
     */
    setSource(value) {
        this.fromJSON(value);
    }
    /**
     * Override the notebook with a JSON-serialized document.
     *
     * @param value The notebook
     */
    fromJSON(value) {
        this.transact(() => {
            this.nbformat = value.nbformat;
            this.nbformat_minor = value.nbformat_minor;
            const metadata = value.metadata;
            if (metadata['orig_nbformat'] !== undefined) {
                delete metadata['orig_nbformat'];
            }
            if (!this.metadata) {
                const ymetadata = new yjs_mjs_.Map();
                for (const [key, value] of Object.entries(metadata)) {
                    ymetadata.set(key, value);
                }
                this.ymeta.set('metadata', ymetadata);
            }
            else {
                this.metadata = metadata;
            }
            const useId = value.nbformat === 4 && value.nbformat_minor >= 5;
            const ycells = value.cells.map(cell => {
                if (!useId) {
                    delete cell.id;
                }
                return cell;
            });
            this.insertCells(this.cells.length, ycells);
            this.deleteCellRange(0, this.cells.length);
        });
    }
    /**
     * Serialize the model to JSON.
     */
    toJSON() {
        // strip cell ids if we have notebook format 4.0-4.4
        const pruneCellId = this.nbformat === 4 && this.nbformat_minor <= 4;
        return {
            metadata: this.metadata,
            nbformat_minor: this.nbformat_minor,
            nbformat: this.nbformat,
            cells: this.cells.map(c => {
                const raw = c.toJSON();
                if (pruneCellId) {
                    delete raw.id;
                }
                return raw;
            })
        };
    }
}
//# sourceMappingURL=ynotebook.js.map
;// CONCATENATED MODULE: ../node_modules/@jupyter/ydoc/lib/index.js
/* -----------------------------------------------------------------------------
| Copyright (c) Jupyter Development Team.
| Distributed under the terms of the Modified BSD License.
|----------------------------------------------------------------------------*/
/**
 * @packageDocumentation
 * @module ydoc
 */








//# sourceMappingURL=index.js.map

/***/ }),

/***/ 79504:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Dp: () => (/* binding */ from),
/* harmony export */   G: () => (/* binding */ some),
/* harmony export */   JJ: () => (/* binding */ unfold),
/* harmony export */   Z$: () => (/* binding */ last),
/* harmony export */   kJ: () => (/* binding */ isArray),
/* harmony export */   s7: () => (/* binding */ appendTo)
/* harmony export */ });
/* unused harmony exports create, copy, every, equalFlat, flatten, fold, unique, uniqueBy, map */
/**
 * Utility module to work with Arrays.
 *
 * @module array
 */



/**
 * Return the last element of an array. The element must exist
 *
 * @template L
 * @param {ArrayLike<L>} arr
 * @return {L}
 */
const last = arr => arr[arr.length - 1]

/**
 * @template C
 * @return {Array<C>}
 */
const create = () => /** @type {Array<C>} */ ([])

/**
 * @template D
 * @param {Array<D>} a
 * @return {Array<D>}
 */
const copy = a => /** @type {Array<D>} */ (a.slice())

/**
 * Append elements from src to dest
 *
 * @template M
 * @param {Array<M>} dest
 * @param {Array<M>} src
 */
const appendTo = (dest, src) => {
  for (let i = 0; i < src.length; i++) {
    dest.push(src[i])
  }
}

/**
 * Transforms something array-like to an actual Array.
 *
 * @function
 * @template T
 * @param {ArrayLike<T>|Iterable<T>} arraylike
 * @return {T}
 */
const from = Array.from

/**
 * True iff condition holds on every element in the Array.
 *
 * @function
 * @template ITEM
 * @template {ArrayLike<ITEM>} ARR
 *
 * @param {ARR} arr
 * @param {function(ITEM, number, ARR):boolean} f
 * @return {boolean}
 */
const every = (arr, f) => {
  for (let i = 0; i < arr.length; i++) {
    if (!f(arr[i], i, arr)) {
      return false
    }
  }
  return true
}

/**
 * True iff condition holds on some element in the Array.
 *
 * @function
 * @template S
 * @template {ArrayLike<S>} ARR
 * @param {ARR} arr
 * @param {function(S, number, ARR):boolean} f
 * @return {boolean}
 */
const some = (arr, f) => {
  for (let i = 0; i < arr.length; i++) {
    if (f(arr[i], i, arr)) {
      return true
    }
  }
  return false
}

/**
 * @template ELEM
 *
 * @param {ArrayLike<ELEM>} a
 * @param {ArrayLike<ELEM>} b
 * @return {boolean}
 */
const equalFlat = (a, b) => a.length === b.length && every(a, (item, index) => item === b[index])

/**
 * @template ELEM
 * @param {Array<Array<ELEM>>} arr
 * @return {Array<ELEM>}
 */
const flatten = arr => fold(arr, /** @type {Array<ELEM>} */ ([]), (acc, val) => acc.concat(val))

/**
 * @template T
 * @param {number} len
 * @param {function(number, Array<T>):T} f
 * @return {Array<T>}
 */
const unfold = (len, f) => {
  const array = new Array(len)
  for (let i = 0; i < len; i++) {
    array[i] = f(i, array)
  }
  return array
}

/**
 * @template T
 * @template RESULT
 * @param {Array<T>} arr
 * @param {RESULT} seed
 * @param {function(RESULT, T, number):RESULT} folder
 */
const fold = (arr, seed, folder) => arr.reduce(folder, seed)

const isArray = Array.isArray

/**
 * @template T
 * @param {Array<T>} arr
 * @return {Array<T>}
 */
const unique = arr => from(set.from(arr))

/**
 * @template T
 * @template M
 * @param {ArrayLike<T>} arr
 * @param {function(T):M} mapper
 * @return {Array<T>}
 */
const uniqueBy = (arr, mapper) => {
  /**
   * @type {Set<M>}
   */
  const happened = set.create()
  /**
   * @type {Array<T>}
   */
  const result = []
  for (let i = 0; i < arr.length; i++) {
    const el = arr[i]
    const mapped = mapper(el)
    if (!happened.has(mapped)) {
      happened.add(mapped)
      result.push(el)
    }
  }
  return result
}

/**
 * @template {ArrayLike<any>} ARR
 * @template {function(ARR extends ArrayLike<infer T> ? T : never, number, ARR):any} MAPPER
 * @param {ARR} arr
 * @param {MAPPER} mapper
 * @return {Array<MAPPER extends function(...any): infer M ? M : never>}
 */
const map = (arr, mapper) => {
  /**
   * @type {Array<any>}
   */
  const res = Array(arr.length)
  for (let i = 0; i < arr.length; i++) {
    res[i] = mapper(/** @type {any} */ (arr[i]), i, /** @type {any} */ (arr))
  }
  return /** @type {any} */ (res)
}


/***/ }),

/***/ 38828:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Hi: () => (/* binding */ equalityDeep),
/* harmony export */   PP: () => (/* binding */ callAll),
/* harmony export */   gB: () => (/* binding */ isOneOf),
/* harmony export */   id: () => (/* binding */ id)
/* harmony export */ });
/* unused harmony exports nop, apply, equalityStrict, equalityFlat, isArray, isString, isNumber, is, isTemplate */
/* harmony import */ var _array_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(79504);
/* harmony import */ var _object_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(36498);
/**
 * Common functions and function call helpers.
 *
 * @module function
 */




/**
 * Calls all functions in `fs` with args. Only throws after all functions were called.
 *
 * @param {Array<function>} fs
 * @param {Array<any>} args
 */
const callAll = (fs, args, i = 0) => {
  try {
    for (; i < fs.length; i++) {
      fs[i](...args)
    }
  } finally {
    if (i < fs.length) {
      callAll(fs, args, i + 1)
    }
  }
}

const nop = () => {}

/**
 * @template T
 * @param {function():T} f
 * @return {T}
 */
const apply = f => f()

/**
 * @template A
 *
 * @param {A} a
 * @return {A}
 */
const id = a => a

/**
 * @template T
 *
 * @param {T} a
 * @param {T} b
 * @return {boolean}
 */
const equalityStrict = (a, b) => a === b

/**
 * @template T
 *
 * @param {Array<T>|object} a
 * @param {Array<T>|object} b
 * @return {boolean}
 */
const equalityFlat = (a, b) => a === b || (a != null && b != null && a.constructor === b.constructor && ((array.isArray(a) && array.equalFlat(a, /** @type {Array<T>} */ (b))) || (typeof a === 'object' && object.equalFlat(a, b))))

/* c8 ignore start */

/**
 * @param {any} a
 * @param {any} b
 * @return {boolean}
 */
const equalityDeep = (a, b) => {
  if (a == null || b == null) {
    return equalityStrict(a, b)
  }
  if (a.constructor !== b.constructor) {
    return false
  }
  if (a === b) {
    return true
  }
  switch (a.constructor) {
    case ArrayBuffer:
      a = new Uint8Array(a)
      b = new Uint8Array(b)
    // eslint-disable-next-line no-fallthrough
    case Uint8Array: {
      if (a.byteLength !== b.byteLength) {
        return false
      }
      for (let i = 0; i < a.length; i++) {
        if (a[i] !== b[i]) {
          return false
        }
      }
      break
    }
    case Set: {
      if (a.size !== b.size) {
        return false
      }
      for (const value of a) {
        if (!b.has(value)) {
          return false
        }
      }
      break
    }
    case Map: {
      if (a.size !== b.size) {
        return false
      }
      for (const key of a.keys()) {
        if (!b.has(key) || !equalityDeep(a.get(key), b.get(key))) {
          return false
        }
      }
      break
    }
    case Object:
      if (_object_js__WEBPACK_IMPORTED_MODULE_0__/* .length */ .kE(a) !== _object_js__WEBPACK_IMPORTED_MODULE_0__/* .length */ .kE(b)) {
        return false
      }
      for (const key in a) {
        if (!_object_js__WEBPACK_IMPORTED_MODULE_0__/* .hasProperty */ .l$(a, key) || !equalityDeep(a[key], b[key])) {
          return false
        }
      }
      break
    case Array:
      if (a.length !== b.length) {
        return false
      }
      for (let i = 0; i < a.length; i++) {
        if (!equalityDeep(a[i], b[i])) {
          return false
        }
      }
      break
    default:
      return false
  }
  return true
}

/**
 * @template V
 * @template {V} OPTS
 *
 * @param {V} value
 * @param {Array<OPTS>} options
 */
// @ts-ignore
const isOneOf = (value, options) => options.includes(value)
/* c8 ignore stop */

const isArray = _array_js__WEBPACK_IMPORTED_MODULE_1__/* .isArray */ .kJ

/**
 * @param {any} s
 * @return {s is String}
 */
const isString = (s) => s && s.constructor === String

/**
 * @param {any} n
 * @return {n is Number}
 */
const isNumber = n => n != null && n.constructor === Number

/**
 * @template {abstract new (...args: any) => any} TYPE
 * @param {any} n
 * @param {TYPE} T
 * @return {n is InstanceType<TYPE>}
 */
const is = (n, T) => n && n.constructor === T

/**
 * @template {abstract new (...args: any) => any} TYPE
 * @param {TYPE} T
 */
const isTemplate = (T) =>
  /**
   * @param {any} n
   * @return {n is InstanceType<TYPE>}
   **/
  n => n && n.constructor === T


/***/ }),

/***/ 22592:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   JG: () => (/* binding */ copy),
/* harmony export */   UI: () => (/* binding */ map),
/* harmony export */   Ue: () => (/* binding */ create),
/* harmony export */   Yj: () => (/* binding */ any),
/* harmony export */   Yu: () => (/* binding */ setIfUndefined)
/* harmony export */ });
/* unused harmony export all */
/**
 * Utility module to work with key-value stores.
 *
 * @module map
 */

/**
 * Creates a new Map instance.
 *
 * @function
 * @return {Map<any, any>}
 *
 * @function
 */
const create = () => new Map()

/**
 * Copy a Map object into a fresh Map object.
 *
 * @function
 * @template X,Y
 * @param {Map<X,Y>} m
 * @return {Map<X,Y>}
 */
const copy = m => {
  const r = create()
  m.forEach((v, k) => { r.set(k, v) })
  return r
}

/**
 * Get map property. Create T if property is undefined and set T on map.
 *
 * ```js
 * const listeners = map.setIfUndefined(events, 'eventName', set.create)
 * listeners.add(listener)
 * ```
 *
 * @function
 * @template V,K
 * @template {Map<K,V>} MAP
 * @param {MAP} map
 * @param {K} key
 * @param {function():V} createT
 * @return {V}
 */
const setIfUndefined = (map, key, createT) => {
  let set = map.get(key)
  if (set === undefined) {
    map.set(key, set = createT())
  }
  return set
}

/**
 * Creates an Array and populates it with the content of all key-value pairs using the `f(value, key)` function.
 *
 * @function
 * @template K
 * @template V
 * @template R
 * @param {Map<K,V>} m
 * @param {function(V,K):R} f
 * @return {Array<R>}
 */
const map = (m, f) => {
  const res = []
  for (const [key, value] of m) {
    res.push(f(value, key))
  }
  return res
}

/**
 * Tests whether any key-value pairs pass the test implemented by `f(value, key)`.
 *
 * @todo should rename to some - similarly to Array.some
 *
 * @function
 * @template K
 * @template V
 * @param {Map<K,V>} m
 * @param {function(V,K):boolean} f
 * @return {boolean}
 */
const any = (m, f) => {
  for (const [key, value] of m) {
    if (f(value, key)) {
      return true
    }
  }
  return false
}

/**
 * Tests whether all key-value pairs pass the test implemented by `f(value, key)`.
 *
 * @function
 * @template K
 * @template V
 * @param {Map<K,V>} m
 * @param {function(V,K):boolean} f
 * @return {boolean}
 */
const all = (m, f) => {
  for (const [key, value] of m) {
    if (!f(value, key)) {
      return false
    }
  }
  return true
}


/***/ }),

/***/ 11182:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Fp: () => (/* binding */ max),
/* harmony export */   GR: () => (/* binding */ isNegativeZero),
/* harmony export */   GW: () => (/* binding */ floor),
/* harmony export */   VV: () => (/* binding */ min),
/* harmony export */   Wn: () => (/* binding */ abs)
/* harmony export */ });
/* unused harmony exports ceil, imul, round, log10, log2, log, sqrt, add, isNaN, pow, exp10, sign */
/**
 * Common Math expressions.
 *
 * @module math
 */

const floor = Math.floor
const ceil = Math.ceil
const abs = Math.abs
const imul = Math.imul
const round = Math.round
const log10 = Math.log10
const log2 = Math.log2
const log = Math.log
const sqrt = Math.sqrt

/**
 * @function
 * @param {number} a
 * @param {number} b
 * @return {number} The sum of a and b
 */
const add = (a, b) => a + b

/**
 * @function
 * @param {number} a
 * @param {number} b
 * @return {number} The smaller element of a and b
 */
const min = (a, b) => a < b ? a : b

/**
 * @function
 * @param {number} a
 * @param {number} b
 * @return {number} The bigger element of a and b
 */
const max = (a, b) => a > b ? a : b

const isNaN = Number.isNaN

const pow = Math.pow
/**
 * Base 10 exponential function. Returns the value of 10 raised to the power of pow.
 *
 * @param {number} exp
 * @return {number}
 */
const exp10 = exp => Math.pow(10, exp)

const sign = Math.sign

/**
 * @param {number} n
 * @return {boolean} Wether n is negative. This function also differentiates between -0 and +0
 */
const isNegativeZero = n => n !== 0 ? n < 0 : 1 / n < 0


/***/ }),

/***/ 36498:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   $m: () => (/* binding */ equalFlat),
/* harmony export */   Ed: () => (/* binding */ forEach),
/* harmony export */   f0: () => (/* binding */ assign),
/* harmony export */   kE: () => (/* binding */ length),
/* harmony export */   l$: () => (/* binding */ hasProperty),
/* harmony export */   xb: () => (/* binding */ isEmpty)
/* harmony export */ });
/* unused harmony exports create, keys, map, some, every */
/**
 * Utility functions for working with EcmaScript objects.
 *
 * @module object
 */

/**
 * @return {Object<string,any>} obj
 */
const create = () => Object.create(null)

/**
 * Object.assign
 */
const assign = Object.assign

/**
 * @param {Object<string,any>} obj
 */
const keys = Object.keys

/**
 * @template V
 * @param {{[k:string]:V}} obj
 * @param {function(V,string):any} f
 */
const forEach = (obj, f) => {
  for (const key in obj) {
    f(obj[key], key)
  }
}

/**
 * @todo implement mapToArray & map
 *
 * @template R
 * @param {Object<string,any>} obj
 * @param {function(any,string):R} f
 * @return {Array<R>}
 */
const map = (obj, f) => {
  const results = []
  for (const key in obj) {
    results.push(f(obj[key], key))
  }
  return results
}

/**
 * @param {Object<string,any>} obj
 * @return {number}
 */
const length = obj => keys(obj).length

/**
 * @param {Object<string,any>} obj
 * @param {function(any,string):boolean} f
 * @return {boolean}
 */
const some = (obj, f) => {
  for (const key in obj) {
    if (f(obj[key], key)) {
      return true
    }
  }
  return false
}

/**
 * @param {Object|undefined} obj
 */
const isEmpty = obj => {
  // eslint-disable-next-line
  for (const _k in obj) {
    return false
  }
  return true
}

/**
 * @param {Object<string,any>} obj
 * @param {function(any,string):boolean} f
 * @return {boolean}
 */
const every = (obj, f) => {
  for (const key in obj) {
    if (!f(obj[key], key)) {
      return false
    }
  }
  return true
}

/**
 * Calls `Object.prototype.hasOwnProperty`.
 *
 * @param {any} obj
 * @param {string|symbol} key
 * @return {boolean}
 */
const hasProperty = (obj, key) => Object.prototype.hasOwnProperty.call(obj, key)

/**
 * @param {Object<string,any>} a
 * @param {Object<string,any>} b
 * @return {boolean}
 */
const equalFlat = (a, b) => a === b || (length(a) === length(b) && every(a, (val, key) => (val !== undefined || hasProperty(b, key)) && b[key] === val))


/***/ }),

/***/ 12330:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   y: () => (/* binding */ Observable)
/* harmony export */ });
/* unused harmony export ObservableV2 */
/* harmony import */ var _map_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(22592);
/* harmony import */ var _set_js__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(79287);
/* harmony import */ var _array_js__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(79504);
/**
 * Observable class prototype.
 *
 * @module observable
 */





/**
 * Handles named events.
 * @experimental
 *
 * This is basically a (better typed) duplicate of Observable, which will replace Observable in the
 * next release.
 *
 * @template {{[key: string]: function(...any):void}} EVENTS
 */
class ObservableV2 {
  constructor () {
    /**
     * Some desc.
     * @type {Map<string, Set<any>>}
     */
    this._observers = map.create()
  }

  /**
   * @template {string} NAME
   * @param {NAME} name
   * @param {EVENTS[NAME]} f
   */
  on (name, f) {
    map.setIfUndefined(this._observers, /** @type {string} */ (name), set.create).add(f)
    return f
  }

  /**
   * @template {string} NAME
   * @param {NAME} name
   * @param {EVENTS[NAME]} f
   */
  once (name, f) {
    /**
     * @param  {...any} args
     */
    const _f = (...args) => {
      this.off(name, /** @type {any} */ (_f))
      f(...args)
    }
    this.on(name, /** @type {any} */ (_f))
  }

  /**
   * @template {string} NAME
   * @param {NAME} name
   * @param {EVENTS[NAME]} f
   */
  off (name, f) {
    const observers = this._observers.get(name)
    if (observers !== undefined) {
      observers.delete(f)
      if (observers.size === 0) {
        this._observers.delete(name)
      }
    }
  }

  /**
   * Emit a named event. All registered event listeners that listen to the
   * specified name will receive the event.
   *
   * @todo This should catch exceptions
   *
   * @template {string} NAME
   * @param {NAME} name The event name.
   * @param {Parameters<EVENTS[NAME]>} args The arguments that are applied to the event listener.
   */
  emit (name, args) {
    // copy all listeners to an array first to make sure that no event is emitted to listeners that are subscribed while the event handler is called.
    return array.from((this._observers.get(name) || map.create()).values()).forEach(f => f(...args))
  }

  destroy () {
    this._observers = map.create()
  }
}

/* c8 ignore start */
/**
 * Handles named events.
 *
 * @deprecated
 * @template N
 */
class Observable {
  constructor () {
    /**
     * Some desc.
     * @type {Map<N, any>}
     */
    this._observers = _map_js__WEBPACK_IMPORTED_MODULE_0__/* .create */ .Ue()
  }

  /**
   * @param {N} name
   * @param {function} f
   */
  on (name, f) {
    _map_js__WEBPACK_IMPORTED_MODULE_0__/* .setIfUndefined */ .Yu(this._observers, name, _set_js__WEBPACK_IMPORTED_MODULE_1__/* .create */ .Ue).add(f)
  }

  /**
   * @param {N} name
   * @param {function} f
   */
  once (name, f) {
    /**
     * @param  {...any} args
     */
    const _f = (...args) => {
      this.off(name, _f)
      f(...args)
    }
    this.on(name, _f)
  }

  /**
   * @param {N} name
   * @param {function} f
   */
  off (name, f) {
    const observers = this._observers.get(name)
    if (observers !== undefined) {
      observers.delete(f)
      if (observers.size === 0) {
        this._observers.delete(name)
      }
    }
  }

  /**
   * Emit a named event. All registered event listeners that listen to the
   * specified name will receive the event.
   *
   * @todo This should catch exceptions
   *
   * @param {N} name The event name.
   * @param {Array<any>} args The arguments that are applied to the event listener.
   */
  emit (name, args) {
    // copy all listeners to an array first to make sure that no event is emitted to listeners that are subscribed while the event handler is called.
    return _array_js__WEBPACK_IMPORTED_MODULE_2__/* .from */ .Dp((this._observers.get(name) || _map_js__WEBPACK_IMPORTED_MODULE_0__/* .create */ .Ue()).values()).forEach(f => f(...args))
  }

  destroy () {
    this._observers = _map_js__WEBPACK_IMPORTED_MODULE_0__/* .create */ .Ue()
  }
}
/* c8 ignore end */


/***/ }),

/***/ 79287:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Ue: () => (/* binding */ create)
/* harmony export */ });
/* unused harmony exports toArray, first, from */
/**
 * Utility module to work with sets.
 *
 * @module set
 */

const create = () => new Set()

/**
 * @template T
 * @param {Set<T>} set
 * @return {Array<T>}
 */
const toArray = set => Array.from(set)

/**
 * @template T
 * @param {Set<T>} set
 * @return {T}
 */
const first = set =>
  set.values().next().value || undefined

/**
 * @template T
 * @param {Iterable<T>} entries
 * @return {Set<T>}
 */
const from = entries => new Set(entries)


/***/ }),

/***/ 2431:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   ZG: () => (/* binding */ getUnixTime)
/* harmony export */ });
/* unused harmony exports getDate, humanizeDuration */
/**
 * Utility module to work with time.
 *
 * @module time
 */




/**
 * Return current time.
 *
 * @return {Date}
 */
const getDate = () => new Date()

/**
 * Return current unix time.
 *
 * @return {number}
 */
const getUnixTime = Date.now

/**
 * Transform time (in ms) to a human readable format. E.g. 1100 => 1.1s. 60s => 1min. .001 => 10μs.
 *
 * @param {number} d duration in milliseconds
 * @return {string} humanized approximation of time
 */
const humanizeDuration = d => {
  if (d < 60000) {
    const p = metric.prefix(d, -1)
    return math.round(p.n * 100) / 100 + p.prefix + 's'
  }
  d = math.floor(d / 1000)
  const seconds = d % 60
  const minutes = math.floor(d / 60) % 60
  const hours = math.floor(d / 3600) % 24
  const days = math.floor(d / 86400)
  if (days > 0) {
    return days + 'd' + ((hours > 0 || minutes > 30) ? ' ' + (minutes > 30 ? hours + 1 : hours) + 'h' : '')
  }
  if (hours > 0) {
    /* c8 ignore next */
    return hours + 'h' + ((minutes > 0 || seconds > 30) ? ' ' + (seconds > 30 ? minutes + 1 : minutes) + 'min' : '')
  }
  return minutes + 'min' + (seconds > 0 ? ' ' + seconds + 's' : '')
}


/***/ })

}]);
//# sourceMappingURL=35.59a288da566759795f5b.js.map?v=59a288da566759795f5b
"use strict";
(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[7561],{

/***/ 20190:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   Z: () => (__WEBPACK_DEFAULT_EXPORT__)
/* harmony export */ });
/* harmony import */ var _baseClone_js__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(52390);


/** Used to compose bitmasks for cloning. */
var CLONE_SYMBOLS_FLAG = 4;

/**
 * Creates a shallow clone of `value`.
 *
 * **Note:** This method is loosely based on the
 * [structured clone algorithm](https://mdn.io/Structured_clone_algorithm)
 * and supports cloning arrays, array buffers, booleans, date objects, maps,
 * numbers, `Object` objects, regexes, sets, strings, symbols, and typed
 * arrays. The own enumerable properties of `arguments` objects are cloned
 * as plain objects. An empty object is returned for uncloneable values such
 * as error objects, functions, DOM nodes, and WeakMaps.
 *
 * @static
 * @memberOf _
 * @since 0.1.0
 * @category Lang
 * @param {*} value The value to clone.
 * @returns {*} Returns the cloned value.
 * @see _.cloneDeep
 * @example
 *
 * var objects = [{ 'a': 1 }, { 'b': 2 }];
 *
 * var shallow = _.clone(objects);
 * console.log(shallow[0] === objects[0]);
 * // => true
 */
function clone(value) {
  return (0,_baseClone_js__WEBPACK_IMPORTED_MODULE_0__/* ["default"] */ .Z)(value, CLONE_SYMBOLS_FLAG);
}

/* harmony default export */ const __WEBPACK_DEFAULT_EXPORT__ = (clone);


/***/ }),

/***/ 7561:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

// ESM COMPAT FLAG
__webpack_require__.r(__webpack_exports__);

// EXPORTS
__webpack_require__.d(__webpack_exports__, {
  render: () => (/* binding */ render)
});

// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-IIMUDSI4.mjs
var chunk_IIMUDSI4 = __webpack_require__(52017);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-VV3M67IP.mjs
var chunk_VV3M67IP = __webpack_require__(83133);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-HRU6DDCH.mjs
var chunk_HRU6DDCH = __webpack_require__(41921);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-K557N5IZ.mjs
var chunk_K557N5IZ = __webpack_require__(81861);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-H2D2JQ3I.mjs + 12 modules
var chunk_H2D2JQ3I = __webpack_require__(62072);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-C3MQ5ANM.mjs
var chunk_C3MQ5ANM = __webpack_require__(39769);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-O4NI6UNU.mjs
var chunk_O4NI6UNU = __webpack_require__(42626);
// EXTERNAL MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/chunk-YTJNT7DU.mjs + 3 modules
var chunk_YTJNT7DU = __webpack_require__(86906);
// EXTERNAL MODULE: ../node_modules/dagre-d3-es/src/dagre/index.js + 62 modules
var dagre = __webpack_require__(24276);
// EXTERNAL MODULE: ../node_modules/lodash-es/isUndefined.js
var isUndefined = __webpack_require__(52307);
// EXTERNAL MODULE: ../node_modules/lodash-es/clone.js
var clone = __webpack_require__(20190);
// EXTERNAL MODULE: ../node_modules/lodash-es/map.js
var map = __webpack_require__(12930);
// EXTERNAL MODULE: ../node_modules/dagre-d3-es/src/graphlib/graph.js + 1 modules
var graph = __webpack_require__(11109);
;// CONCATENATED MODULE: ../node_modules/dagre-d3-es/src/graphlib/json.js





function write(g) {
  var json = {
    options: {
      directed: g.isDirected(),
      multigraph: g.isMultigraph(),
      compound: g.isCompound(),
    },
    nodes: writeNodes(g),
    edges: writeEdges(g),
  };
  if (!isUndefined/* default */.Z(g.graph())) {
    json.value = clone/* default */.Z(g.graph());
  }
  return json;
}

function writeNodes(g) {
  return map/* default */.Z(g.nodes(), function (v) {
    var nodeValue = g.node(v);
    var parent = g.parent(v);
    var node = { v: v };
    if (!isUndefined/* default */.Z(nodeValue)) {
      node.value = nodeValue;
    }
    if (!isUndefined/* default */.Z(parent)) {
      node.parent = parent;
    }
    return node;
  });
}

function writeEdges(g) {
  return map/* default */.Z(g.edges(), function (e) {
    var edgeValue = g.edge(e);
    var edge = { v: e.v, w: e.w };
    if (!isUndefined/* default */.Z(e.name)) {
      edge.name = e.name;
    }
    if (!isUndefined/* default */.Z(edgeValue)) {
      edge.value = edgeValue;
    }
    return edge;
  });
}

function read(json) {
  var g = new Graph(json.options).setGraph(json.value);
  _.each(json.nodes, function (entry) {
    g.setNode(entry.v, entry.value);
    if (entry.parent) {
      g.setParent(entry.v, entry.parent);
    }
  });
  _.each(json.edges, function (entry) {
    g.setEdge({ v: entry.v, w: entry.w, name: entry.name }, entry.value);
  });
  return g;
}

// EXTERNAL MODULE: ../node_modules/dagre-d3-es/src/graphlib/index.js
var graphlib = __webpack_require__(67406);
;// CONCATENATED MODULE: ../node_modules/mermaid/dist/chunks/mermaid.core/dagre-OKDRZEBW.mjs









// src/rendering-util/layout-algorithms/dagre/index.js




// src/rendering-util/layout-algorithms/dagre/mermaid-graphlib.js


var clusterDb = /* @__PURE__ */ new Map();
var descendants = /* @__PURE__ */ new Map();
var parents = /* @__PURE__ */ new Map();
var clear4 = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)(() => {
  descendants.clear();
  parents.clear();
  clusterDb.clear();
}, "clear");
var isDescendant = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)((id, ancestorId) => {
  const ancestorDescendants = descendants.get(ancestorId) || [];
  chunk_YTJNT7DU/* log */.cM.trace("In isDescendant", ancestorId, " ", id, " = ", ancestorDescendants.includes(id));
  return ancestorDescendants.includes(id);
}, "isDescendant");
var edgeInCluster = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)((edge, clusterId) => {
  const clusterDescendants = descendants.get(clusterId) || [];
  chunk_YTJNT7DU/* log */.cM.info("Descendants of ", clusterId, " is ", clusterDescendants);
  chunk_YTJNT7DU/* log */.cM.info("Edge is ", edge);
  if (edge.v === clusterId || edge.w === clusterId) {
    return false;
  }
  if (!clusterDescendants) {
    chunk_YTJNT7DU/* log */.cM.debug("Tilt, ", clusterId, ",not in descendants");
    return false;
  }
  return clusterDescendants.includes(edge.v) || isDescendant(edge.v, clusterId) || isDescendant(edge.w, clusterId) || clusterDescendants.includes(edge.w);
}, "edgeInCluster");
var copy = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)((clusterId, graph, newGraph, rootId) => {
  chunk_YTJNT7DU/* log */.cM.warn(
    "Copying children of ",
    clusterId,
    "root",
    rootId,
    "data",
    graph.node(clusterId),
    rootId
  );
  const nodes = graph.children(clusterId) || [];
  if (clusterId !== rootId) {
    nodes.push(clusterId);
  }
  chunk_YTJNT7DU/* log */.cM.warn("Copying (nodes) clusterId", clusterId, "nodes", nodes);
  nodes.forEach((node) => {
    if (graph.children(node).length > 0) {
      copy(node, graph, newGraph, rootId);
    } else {
      const data = graph.node(node);
      chunk_YTJNT7DU/* log */.cM.info("cp ", node, " to ", rootId, " with parent ", clusterId);
      newGraph.setNode(node, data);
      if (rootId !== graph.parent(node)) {
        chunk_YTJNT7DU/* log */.cM.warn("Setting parent", node, graph.parent(node));
        newGraph.setParent(node, graph.parent(node));
      }
      if (clusterId !== rootId && node !== clusterId) {
        chunk_YTJNT7DU/* log */.cM.debug("Setting parent", node, clusterId);
        newGraph.setParent(node, clusterId);
      } else {
        chunk_YTJNT7DU/* log */.cM.info("In copy ", clusterId, "root", rootId, "data", graph.node(clusterId), rootId);
        chunk_YTJNT7DU/* log */.cM.debug(
          "Not Setting parent for node=",
          node,
          "cluster!==rootId",
          clusterId !== rootId,
          "node!==clusterId",
          node !== clusterId
        );
      }
      const edges = graph.edges(node);
      chunk_YTJNT7DU/* log */.cM.debug("Copying Edges", edges);
      edges.forEach((edge) => {
        chunk_YTJNT7DU/* log */.cM.info("Edge", edge);
        const data2 = graph.edge(edge.v, edge.w, edge.name);
        chunk_YTJNT7DU/* log */.cM.info("Edge data", data2, rootId);
        try {
          if (edgeInCluster(edge, rootId)) {
            chunk_YTJNT7DU/* log */.cM.info("Copying as ", edge.v, edge.w, data2, edge.name);
            newGraph.setEdge(edge.v, edge.w, data2, edge.name);
            chunk_YTJNT7DU/* log */.cM.info("newGraph edges ", newGraph.edges(), newGraph.edge(newGraph.edges()[0]));
          } else {
            chunk_YTJNT7DU/* log */.cM.info(
              "Skipping copy of edge ",
              edge.v,
              "-->",
              edge.w,
              " rootId: ",
              rootId,
              " clusterId:",
              clusterId
            );
          }
        } catch (e) {
          chunk_YTJNT7DU/* log */.cM.error(e);
        }
      });
    }
    chunk_YTJNT7DU/* log */.cM.debug("Removing node", node);
    graph.removeNode(node);
  });
}, "copy");
var extractDescendants = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)((id, graph) => {
  const children = graph.children(id);
  let res = [...children];
  for (const child of children) {
    parents.set(child, id);
    res = [...res, ...extractDescendants(child, graph)];
  }
  return res;
}, "extractDescendants");
var findCommonEdges = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)((graph, id1, id2) => {
  const edges1 = graph.edges().filter((edge) => edge.v === id1 || edge.w === id1);
  const edges2 = graph.edges().filter((edge) => edge.v === id2 || edge.w === id2);
  const edges1Prim = edges1.map((edge) => {
    return { v: edge.v === id1 ? id2 : edge.v, w: edge.w === id1 ? id1 : edge.w };
  });
  const edges2Prim = edges2.map((edge) => {
    return { v: edge.v, w: edge.w };
  });
  const result = edges1Prim.filter((edgeIn1) => {
    return edges2Prim.some((edge) => edgeIn1.v === edge.v && edgeIn1.w === edge.w);
  });
  return result;
}, "findCommonEdges");
var findNonClusterChild = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)((id, graph, clusterId) => {
  const children = graph.children(id);
  chunk_YTJNT7DU/* log */.cM.trace("Searching children of id ", id, children);
  if (children.length < 1) {
    return id;
  }
  let reserve;
  for (const child of children) {
    const _id = findNonClusterChild(child, graph, clusterId);
    const commonEdges = findCommonEdges(graph, clusterId, _id);
    if (_id) {
      if (commonEdges.length > 0) {
        reserve = _id;
      } else {
        return _id;
      }
    }
  }
  return reserve;
}, "findNonClusterChild");
var getAnchorId = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)((id) => {
  if (!clusterDb.has(id)) {
    return id;
  }
  if (!clusterDb.get(id).externalConnections) {
    return id;
  }
  if (clusterDb.has(id)) {
    return clusterDb.get(id).id;
  }
  return id;
}, "getAnchorId");
var adjustClustersAndEdges = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)((graph, depth) => {
  if (!graph || depth > 10) {
    chunk_YTJNT7DU/* log */.cM.debug("Opting out, no graph ");
    return;
  } else {
    chunk_YTJNT7DU/* log */.cM.debug("Opting in, graph ");
  }
  graph.nodes().forEach(function(id) {
    const children = graph.children(id);
    if (children.length > 0) {
      chunk_YTJNT7DU/* log */.cM.warn(
        "Cluster identified",
        id,
        " Replacement id in edges: ",
        findNonClusterChild(id, graph, id)
      );
      descendants.set(id, extractDescendants(id, graph));
      clusterDb.set(id, { id: findNonClusterChild(id, graph, id), clusterData: graph.node(id) });
    }
  });
  graph.nodes().forEach(function(id) {
    const children = graph.children(id);
    const edges = graph.edges();
    if (children.length > 0) {
      chunk_YTJNT7DU/* log */.cM.debug("Cluster identified", id, descendants);
      edges.forEach((edge) => {
        const d1 = isDescendant(edge.v, id);
        const d2 = isDescendant(edge.w, id);
        if (d1 ^ d2) {
          chunk_YTJNT7DU/* log */.cM.warn("Edge: ", edge, " leaves cluster ", id);
          chunk_YTJNT7DU/* log */.cM.warn("Descendants of XXX ", id, ": ", descendants.get(id));
          clusterDb.get(id).externalConnections = true;
        }
      });
    } else {
      chunk_YTJNT7DU/* log */.cM.debug("Not a cluster ", id, descendants);
    }
  });
  for (let id of clusterDb.keys()) {
    const nonClusterChild = clusterDb.get(id).id;
    const parent = graph.parent(nonClusterChild);
    if (parent !== id && clusterDb.has(parent) && !clusterDb.get(parent).externalConnections) {
      clusterDb.get(id).id = parent;
    }
  }
  graph.edges().forEach(function(e) {
    const edge = graph.edge(e);
    chunk_YTJNT7DU/* log */.cM.warn("Edge " + e.v + " -> " + e.w + ": " + JSON.stringify(e));
    chunk_YTJNT7DU/* log */.cM.warn("Edge " + e.v + " -> " + e.w + ": " + JSON.stringify(graph.edge(e)));
    let v = e.v;
    let w = e.w;
    chunk_YTJNT7DU/* log */.cM.warn(
      "Fix XXX",
      clusterDb,
      "ids:",
      e.v,
      e.w,
      "Translating: ",
      clusterDb.get(e.v),
      " --- ",
      clusterDb.get(e.w)
    );
    if (clusterDb.get(e.v) || clusterDb.get(e.w)) {
      chunk_YTJNT7DU/* log */.cM.warn("Fixing and trying - removing XXX", e.v, e.w, e.name);
      v = getAnchorId(e.v);
      w = getAnchorId(e.w);
      graph.removeEdge(e.v, e.w, e.name);
      if (v !== e.v) {
        const parent = graph.parent(v);
        clusterDb.get(parent).externalConnections = true;
        edge.fromCluster = e.v;
      }
      if (w !== e.w) {
        const parent = graph.parent(w);
        clusterDb.get(parent).externalConnections = true;
        edge.toCluster = e.w;
      }
      chunk_YTJNT7DU/* log */.cM.warn("Fix Replacing with XXX", v, w, e.name);
      graph.setEdge(v, w, edge, e.name);
    }
  });
  chunk_YTJNT7DU/* log */.cM.warn("Adjusted Graph", write(graph));
  extractor(graph, 0);
  chunk_YTJNT7DU/* log */.cM.trace(clusterDb);
}, "adjustClustersAndEdges");
var extractor = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)((graph, depth) => {
  chunk_YTJNT7DU/* log */.cM.warn("extractor - ", depth, write(graph), graph.children("D"));
  if (depth > 10) {
    chunk_YTJNT7DU/* log */.cM.error("Bailing out");
    return;
  }
  let nodes = graph.nodes();
  let hasChildren = false;
  for (const node of nodes) {
    const children = graph.children(node);
    hasChildren = hasChildren || children.length > 0;
  }
  if (!hasChildren) {
    chunk_YTJNT7DU/* log */.cM.debug("Done, no node has children", graph.nodes());
    return;
  }
  chunk_YTJNT7DU/* log */.cM.debug("Nodes = ", nodes, depth);
  for (const node of nodes) {
    chunk_YTJNT7DU/* log */.cM.debug(
      "Extracting node",
      node,
      clusterDb,
      clusterDb.has(node) && !clusterDb.get(node).externalConnections,
      !graph.parent(node),
      graph.node(node),
      graph.children("D"),
      " Depth ",
      depth
    );
    if (!clusterDb.has(node)) {
      chunk_YTJNT7DU/* log */.cM.debug("Not a cluster", node, depth);
    } else if (!clusterDb.get(node).externalConnections && graph.children(node) && graph.children(node).length > 0) {
      chunk_YTJNT7DU/* log */.cM.warn(
        "Cluster without external connections, without a parent and with children",
        node,
        depth
      );
      const graphSettings = graph.graph();
      let dir = graphSettings.rankdir === "TB" ? "LR" : "TB";
      if (clusterDb.get(node)?.clusterData?.dir) {
        dir = clusterDb.get(node).clusterData.dir;
        chunk_YTJNT7DU/* log */.cM.warn("Fixing dir", clusterDb.get(node).clusterData.dir, dir);
      }
      const clusterGraph = new graphlib/* Graph */.k({
        multigraph: true,
        compound: true
      }).setGraph({
        rankdir: dir,
        nodesep: 50,
        ranksep: 50,
        marginx: 8,
        marginy: 8
      }).setDefaultEdgeLabel(function() {
        return {};
      });
      chunk_YTJNT7DU/* log */.cM.warn("Old graph before copy", write(graph));
      copy(node, graph, clusterGraph, node);
      graph.setNode(node, {
        clusterNode: true,
        id: node,
        clusterData: clusterDb.get(node).clusterData,
        label: clusterDb.get(node).label,
        graph: clusterGraph
      });
      chunk_YTJNT7DU/* log */.cM.warn("New graph after copy node: (", node, ")", write(clusterGraph));
      chunk_YTJNT7DU/* log */.cM.debug("Old graph after copy", write(graph));
    } else {
      chunk_YTJNT7DU/* log */.cM.warn(
        "Cluster ** ",
        node,
        " **not meeting the criteria !externalConnections:",
        !clusterDb.get(node).externalConnections,
        " no parent: ",
        !graph.parent(node),
        " children ",
        graph.children(node) && graph.children(node).length > 0,
        graph.children("D"),
        depth
      );
      chunk_YTJNT7DU/* log */.cM.debug(clusterDb);
    }
  }
  nodes = graph.nodes();
  chunk_YTJNT7DU/* log */.cM.warn("New list of nodes", nodes);
  for (const node of nodes) {
    const data = graph.node(node);
    chunk_YTJNT7DU/* log */.cM.warn(" Now next level", node, data);
    if (data?.clusterNode) {
      extractor(data.graph, depth + 1);
    }
  }
}, "extractor");
var sorter = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)((graph, nodes) => {
  if (nodes.length === 0) {
    return [];
  }
  let result = Object.assign([], nodes);
  nodes.forEach((node) => {
    const children = graph.children(node);
    const sorted = sorter(graph, children);
    result = [...result, ...sorted];
  });
  return result;
}, "sorter");
var sortNodesByHierarchy = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)((graph) => sorter(graph, graph.children()), "sortNodesByHierarchy");

// src/rendering-util/layout-algorithms/dagre/index.js
var recursiveRender = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)(async (_elem, graph, diagramType, id, parentCluster, siteConfig) => {
  chunk_YTJNT7DU/* log */.cM.warn("Graph in recursive render:XAX", write(graph), parentCluster);
  const dir = graph.graph().rankdir;
  chunk_YTJNT7DU/* log */.cM.trace("Dir in recursive render - dir:", dir);
  const elem = _elem.insert("g").attr("class", "root");
  if (!graph.nodes()) {
    chunk_YTJNT7DU/* log */.cM.info("No nodes found for", graph);
  } else {
    chunk_YTJNT7DU/* log */.cM.info("Recursive render XXX", graph.nodes());
  }
  if (graph.edges().length > 0) {
    chunk_YTJNT7DU/* log */.cM.info("Recursive edges", graph.edge(graph.edges()[0]));
  }
  const clusters = elem.insert("g").attr("class", "clusters");
  const edgePaths = elem.insert("g").attr("class", "edgePaths");
  const edgeLabels = elem.insert("g").attr("class", "edgeLabels");
  const nodes = elem.insert("g").attr("class", "nodes");
  await Promise.all(
    graph.nodes().map(async function(v) {
      const node = graph.node(v);
      if (parentCluster !== void 0) {
        const data = JSON.parse(JSON.stringify(parentCluster.clusterData));
        chunk_YTJNT7DU/* log */.cM.trace(
          "Setting data for parent cluster XXX\n Node.id = ",
          v,
          "\n data=",
          data.height,
          "\nParent cluster",
          parentCluster.height
        );
        graph.setNode(parentCluster.id, data);
        if (!graph.parent(v)) {
          chunk_YTJNT7DU/* log */.cM.trace("Setting parent", v, parentCluster.id);
          graph.setParent(v, parentCluster.id, data);
        }
      }
      chunk_YTJNT7DU/* log */.cM.info("(Insert) Node XXX" + v + ": " + JSON.stringify(graph.node(v)));
      if (node?.clusterNode) {
        chunk_YTJNT7DU/* log */.cM.info("Cluster identified XBX", v, node.width, graph.node(v));
        const { ranksep, nodesep } = graph.graph();
        node.graph.setGraph({
          ...node.graph.graph(),
          ranksep: ranksep + 25,
          nodesep
        });
        const o = await recursiveRender(
          nodes,
          node.graph,
          diagramType,
          id,
          graph.node(v),
          siteConfig
        );
        const newEl = o.elem;
        (0,chunk_HRU6DDCH/* updateNodeBounds */.jr)(node, newEl);
        node.diff = o.diff || 0;
        chunk_YTJNT7DU/* log */.cM.info(
          "New compound node after recursive render XAX",
          v,
          "width",
          // node,
          node.width,
          "height",
          node.height
          // node.x,
          // node.y
        );
        (0,chunk_HRU6DDCH/* setNodeElem */.Yn)(newEl, node);
      } else {
        if (graph.children(v).length > 0) {
          chunk_YTJNT7DU/* log */.cM.trace(
            "Cluster - the non recursive path XBX",
            v,
            node.id,
            node,
            node.width,
            "Graph:",
            graph
          );
          chunk_YTJNT7DU/* log */.cM.trace(findNonClusterChild(node.id, graph));
          clusterDb.set(node.id, { id: findNonClusterChild(node.id, graph), node });
        } else {
          chunk_YTJNT7DU/* log */.cM.trace("Node - the non recursive path XAX", v, nodes, graph.node(v), dir);
          await (0,chunk_HRU6DDCH/* insertNode */.Lf)(nodes, graph.node(v), { config: siteConfig, dir });
        }
      }
    })
  );
  const processEdges = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)(async () => {
    const edgePromises = graph.edges().map(async function(e) {
      const edge = graph.edge(e.v, e.w, e.name);
      chunk_YTJNT7DU/* log */.cM.info("Edge " + e.v + " -> " + e.w + ": " + JSON.stringify(e));
      chunk_YTJNT7DU/* log */.cM.info("Edge " + e.v + " -> " + e.w + ": ", e, " ", JSON.stringify(graph.edge(e)));
      chunk_YTJNT7DU/* log */.cM.info(
        "Fix",
        clusterDb,
        "ids:",
        e.v,
        e.w,
        "Translating: ",
        clusterDb.get(e.v),
        clusterDb.get(e.w)
      );
      await (0,chunk_IIMUDSI4/* insertEdgeLabel */.I_)(edgeLabels, edge);
    });
    await Promise.all(edgePromises);
  }, "processEdges");
  await processEdges();
  chunk_YTJNT7DU/* log */.cM.info("Graph before layout:", JSON.stringify(write(graph)));
  chunk_YTJNT7DU/* log */.cM.info("############################################# XXX");
  chunk_YTJNT7DU/* log */.cM.info("###                Layout                 ### XXX");
  chunk_YTJNT7DU/* log */.cM.info("############################################# XXX");
  (0,dagre/* layout */.bK)(graph);
  chunk_YTJNT7DU/* log */.cM.info("Graph after layout:", JSON.stringify(write(graph)));
  let diff = 0;
  let { subGraphTitleTotalMargin } = (0,chunk_K557N5IZ/* getSubGraphTitleMargins */.L)(siteConfig);
  await Promise.all(
    sortNodesByHierarchy(graph).map(async function(v) {
      const node = graph.node(v);
      chunk_YTJNT7DU/* log */.cM.info(
        "Position XBX => " + v + ": (" + node.x,
        "," + node.y,
        ") width: ",
        node.width,
        " height: ",
        node.height
      );
      if (node?.clusterNode) {
        node.y += subGraphTitleTotalMargin;
        chunk_YTJNT7DU/* log */.cM.info(
          "A tainted cluster node XBX1",
          v,
          node.id,
          node.width,
          node.height,
          node.x,
          node.y,
          graph.parent(v)
        );
        clusterDb.get(node.id).node = node;
        (0,chunk_HRU6DDCH/* positionNode */.aH)(node);
      } else {
        if (graph.children(v).length > 0) {
          chunk_YTJNT7DU/* log */.cM.info(
            "A pure cluster node XBX1",
            v,
            node.id,
            node.x,
            node.y,
            node.width,
            node.height,
            graph.parent(v)
          );
          node.height += subGraphTitleTotalMargin;
          graph.node(node.parentId);
          const halfPadding = node?.padding / 2 || 0;
          const labelHeight = node?.labelBBox?.height || 0;
          const offsetY = labelHeight - halfPadding || 0;
          chunk_YTJNT7DU/* log */.cM.debug("OffsetY", offsetY, "labelHeight", labelHeight, "halfPadding", halfPadding);
          await (0,chunk_HRU6DDCH/* insertCluster */.us)(clusters, node);
          clusterDb.get(node.id).node = node;
        } else {
          const parent = graph.node(node.parentId);
          node.y += subGraphTitleTotalMargin / 2;
          chunk_YTJNT7DU/* log */.cM.info(
            "A regular node XBX1 - using the padding",
            node.id,
            "parent",
            node.parentId,
            node.width,
            node.height,
            node.x,
            node.y,
            "offsetY",
            node.offsetY,
            "parent",
            parent,
            parent?.offsetY,
            node
          );
          (0,chunk_HRU6DDCH/* positionNode */.aH)(node);
        }
      }
    })
  );
  graph.edges().forEach(function(e) {
    const edge = graph.edge(e);
    chunk_YTJNT7DU/* log */.cM.info("Edge " + e.v + " -> " + e.w + ": " + JSON.stringify(edge), edge);
    edge.points.forEach((point) => point.y += subGraphTitleTotalMargin / 2);
    const startNode = graph.node(e.v);
    var endNode = graph.node(e.w);
    const paths = (0,chunk_IIMUDSI4/* insertEdge */.QP)(edgePaths, edge, clusterDb, diagramType, startNode, endNode, id);
    (0,chunk_IIMUDSI4/* positionEdgeLabel */.Jj)(edge, paths);
  });
  graph.nodes().forEach(function(v) {
    const n = graph.node(v);
    chunk_YTJNT7DU/* log */.cM.info(v, n.type, n.diff);
    if (n.isGroup) {
      diff = n.diff;
    }
  });
  chunk_YTJNT7DU/* log */.cM.warn("Returning from recursive render XAX", elem, diff);
  return { elem, diff };
}, "recursiveRender");
var render = /* @__PURE__ */ (0,chunk_YTJNT7DU/* __name */.eW)(async (data4Layout, svg) => {
  const graph = new graphlib/* Graph */.k({
    multigraph: true,
    compound: true
  }).setGraph({
    rankdir: data4Layout.direction,
    nodesep: data4Layout.config?.nodeSpacing || data4Layout.config?.flowchart?.nodeSpacing || data4Layout.nodeSpacing,
    ranksep: data4Layout.config?.rankSpacing || data4Layout.config?.flowchart?.rankSpacing || data4Layout.rankSpacing,
    marginx: 8,
    marginy: 8
  }).setDefaultEdgeLabel(function() {
    return {};
  });
  const element = svg.select("g");
  (0,chunk_IIMUDSI4/* markers_default */.DQ)(element, data4Layout.markers, data4Layout.type, data4Layout.diagramId);
  (0,chunk_HRU6DDCH/* clear2 */.gU)();
  (0,chunk_IIMUDSI4/* clear */.ZH)();
  (0,chunk_HRU6DDCH/* clear */.ZH)();
  clear4();
  data4Layout.nodes.forEach((node) => {
    graph.setNode(node.id, { ...node });
    if (node.parentId) {
      graph.setParent(node.id, node.parentId);
    }
  });
  chunk_YTJNT7DU/* log */.cM.debug("Edges:", data4Layout.edges);
  data4Layout.edges.forEach((edge) => {
    if (edge.start === edge.end) {
      const nodeId = edge.start;
      const specialId1 = nodeId + "---" + nodeId + "---1";
      const specialId2 = nodeId + "---" + nodeId + "---2";
      const node = graph.node(nodeId);
      graph.setNode(specialId1, {
        domId: specialId1,
        id: specialId1,
        parentId: node.parentId,
        labelStyle: "",
        label: "",
        padding: 0,
        shape: "labelRect",
        // shape: 'rect',
        style: "",
        width: 10,
        height: 10
      });
      graph.setParent(specialId1, node.parentId);
      graph.setNode(specialId2, {
        domId: specialId2,
        id: specialId2,
        parentId: node.parentId,
        labelStyle: "",
        padding: 0,
        // shape: 'rect',
        shape: "labelRect",
        label: "",
        style: "",
        width: 10,
        height: 10
      });
      graph.setParent(specialId2, node.parentId);
      const edge1 = structuredClone(edge);
      const edgeMid = structuredClone(edge);
      const edge2 = structuredClone(edge);
      edge1.label = "";
      edge1.arrowTypeEnd = "none";
      edge1.id = nodeId + "-cyclic-special-1";
      edgeMid.arrowTypeStart = "none";
      edgeMid.arrowTypeEnd = "none";
      edgeMid.id = nodeId + "-cyclic-special-mid";
      edge2.label = "";
      if (node.isGroup) {
        edge1.fromCluster = nodeId;
        edge2.toCluster = nodeId;
      }
      edge2.id = nodeId + "-cyclic-special-2";
      edge2.arrowTypeStart = "none";
      graph.setEdge(nodeId, specialId1, edge1, nodeId + "-cyclic-special-0");
      graph.setEdge(specialId1, specialId2, edgeMid, nodeId + "-cyclic-special-1");
      graph.setEdge(specialId2, nodeId, edge2, nodeId + "-cyc<lic-special-2");
    } else {
      graph.setEdge(edge.start, edge.end, { ...edge }, edge.id);
    }
  });
  chunk_YTJNT7DU/* log */.cM.warn("Graph at first:", JSON.stringify(write(graph)));
  adjustClustersAndEdges(graph);
  chunk_YTJNT7DU/* log */.cM.warn("Graph after XAX:", JSON.stringify(write(graph)));
  const siteConfig = (0,chunk_YTJNT7DU/* getConfig2 */.nV)();
  await recursiveRender(
    element,
    graph,
    data4Layout.type,
    data4Layout.diagramId,
    void 0,
    siteConfig
  );
}, "render");



/***/ })

}]);
//# sourceMappingURL=7561.ab0621a9e054b91897f7.js.map?v=ab0621a9e054b91897f7
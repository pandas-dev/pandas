(self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] = self["webpackChunk_JUPYTERLAB_CORE_OUTPUT"] || []).push([[2425],{

/***/ 29287:
/***/ (function(module, __unused_webpack_exports, __webpack_require__) {

(function webpackUniversalModuleDefinition(root, factory) {
	if(true)
		module.exports = factory(__webpack_require__(75925));
	else {}
})(this, function(__WEBPACK_EXTERNAL_MODULE__140__) {
return /******/ (() => { // webpackBootstrap
/******/ 	"use strict";
/******/ 	var __webpack_modules__ = ({

/***/ 658:
/***/ ((module) => {



// Simple, internal Object.assign() polyfill for options objects etc.

module.exports = Object.assign != null ? Object.assign.bind(Object) : function (tgt) {
  for (var _len = arguments.length, srcs = Array(_len > 1 ? _len - 1 : 0), _key = 1; _key < _len; _key++) {
    srcs[_key - 1] = arguments[_key];
  }

  srcs.forEach(function (src) {
    Object.keys(src).forEach(function (k) {
      return tgt[k] = src[k];
    });
  });

  return tgt;
};

/***/ }),

/***/ 548:
/***/ ((module, __unused_webpack_exports, __nested_webpack_require_1129__) => {



var _slicedToArray = function () { function sliceIterator(arr, i) { var _arr = []; var _n = true; var _d = false; var _e = undefined; try { for (var _i = arr[Symbol.iterator](), _s; !(_n = (_s = _i.next()).done); _n = true) { _arr.push(_s.value); if (i && _arr.length === i) break; } } catch (err) { _d = true; _e = err; } finally { try { if (!_n && _i["return"]) _i["return"](); } finally { if (_d) throw _e; } } return _arr; } return function (arr, i) { if (Array.isArray(arr)) { return arr; } else if (Symbol.iterator in Object(arr)) { return sliceIterator(arr, i); } else { throw new TypeError("Invalid attempt to destructure non-iterable instance"); } }; }();

/*
 * Auxiliary functions
 */

var LinkedList = __nested_webpack_require_1129__(140).layoutBase.LinkedList;

var auxiliary = {};

// get the top most nodes
auxiliary.getTopMostNodes = function (nodes) {
  var nodesMap = {};
  for (var i = 0; i < nodes.length; i++) {
    nodesMap[nodes[i].id()] = true;
  }
  var roots = nodes.filter(function (ele, i) {
    if (typeof ele === "number") {
      ele = i;
    }
    var parent = ele.parent()[0];
    while (parent != null) {
      if (nodesMap[parent.id()]) {
        return false;
      }
      parent = parent.parent()[0];
    }
    return true;
  });

  return roots;
};

// find disconnected components and create dummy nodes that connect them
auxiliary.connectComponents = function (cy, eles, topMostNodes, dummyNodes) {
  var queue = new LinkedList();
  var visited = new Set();
  var visitedTopMostNodes = [];
  var currentNeighbor = void 0;
  var minDegreeNode = void 0;
  var minDegree = void 0;

  var isConnected = false;
  var count = 1;
  var nodesConnectedToDummy = [];
  var components = [];

  var _loop = function _loop() {
    var cmpt = cy.collection();
    components.push(cmpt);

    var currentNode = topMostNodes[0];
    var childrenOfCurrentNode = cy.collection();
    childrenOfCurrentNode.merge(currentNode).merge(currentNode.descendants().intersection(eles));
    visitedTopMostNodes.push(currentNode);

    childrenOfCurrentNode.forEach(function (node) {
      queue.push(node);
      visited.add(node);
      cmpt.merge(node);
    });

    var _loop2 = function _loop2() {
      currentNode = queue.shift();

      // Traverse all neighbors of this node
      var neighborNodes = cy.collection();
      currentNode.neighborhood().nodes().forEach(function (node) {
        if (eles.intersection(currentNode.edgesWith(node)).length > 0) {
          neighborNodes.merge(node);
        }
      });

      for (var i = 0; i < neighborNodes.length; i++) {
        var neighborNode = neighborNodes[i];
        currentNeighbor = topMostNodes.intersection(neighborNode.union(neighborNode.ancestors()));
        if (currentNeighbor != null && !visited.has(currentNeighbor[0])) {
          var childrenOfNeighbor = currentNeighbor.union(currentNeighbor.descendants());

          childrenOfNeighbor.forEach(function (node) {
            queue.push(node);
            visited.add(node);
            cmpt.merge(node);
            if (topMostNodes.has(node)) {
              visitedTopMostNodes.push(node);
            }
          });
        }
      }
    };

    while (queue.length != 0) {
      _loop2();
    }

    cmpt.forEach(function (node) {
      eles.intersection(node.connectedEdges()).forEach(function (e) {
        // connectedEdges() usually cached
        if (cmpt.has(e.source()) && cmpt.has(e.target())) {
          // has() is cheap
          cmpt.merge(e);
        }
      });
    });

    if (visitedTopMostNodes.length == topMostNodes.length) {
      isConnected = true;
    }

    if (!isConnected || isConnected && count > 1) {
      minDegreeNode = visitedTopMostNodes[0];
      minDegree = minDegreeNode.connectedEdges().length;
      visitedTopMostNodes.forEach(function (node) {
        if (node.connectedEdges().length < minDegree) {
          minDegree = node.connectedEdges().length;
          minDegreeNode = node;
        }
      });
      nodesConnectedToDummy.push(minDegreeNode.id());
      // TO DO: Check efficiency of this part
      var temp = cy.collection();
      temp.merge(visitedTopMostNodes[0]);
      visitedTopMostNodes.forEach(function (node) {
        temp.merge(node);
      });
      visitedTopMostNodes = [];
      topMostNodes = topMostNodes.difference(temp);
      count++;
    }
  };

  do {
    _loop();
  } while (!isConnected);

  if (dummyNodes) {
    if (nodesConnectedToDummy.length > 0) {
      dummyNodes.set('dummy' + (dummyNodes.size + 1), nodesConnectedToDummy);
    }
  }
  return components;
};

// relocates componentResult to originalCenter if there is no fixedNodeConstraint
auxiliary.relocateComponent = function (originalCenter, componentResult, options) {
  if (!options.fixedNodeConstraint) {
    var minXCoord = Number.POSITIVE_INFINITY;
    var maxXCoord = Number.NEGATIVE_INFINITY;
    var minYCoord = Number.POSITIVE_INFINITY;
    var maxYCoord = Number.NEGATIVE_INFINITY;
    if (options.quality == "draft") {
      // calculate current bounding box
      var _iteratorNormalCompletion = true;
      var _didIteratorError = false;
      var _iteratorError = undefined;

      try {
        for (var _iterator = componentResult.nodeIndexes[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
          var _ref = _step.value;

          var _ref2 = _slicedToArray(_ref, 2);

          var key = _ref2[0];
          var value = _ref2[1];

          var cyNode = options.cy.getElementById(key);
          if (cyNode) {
            var nodeBB = cyNode.boundingBox();
            var leftX = componentResult.xCoords[value] - nodeBB.w / 2;
            var rightX = componentResult.xCoords[value] + nodeBB.w / 2;
            var topY = componentResult.yCoords[value] - nodeBB.h / 2;
            var bottomY = componentResult.yCoords[value] + nodeBB.h / 2;

            if (leftX < minXCoord) minXCoord = leftX;
            if (rightX > maxXCoord) maxXCoord = rightX;
            if (topY < minYCoord) minYCoord = topY;
            if (bottomY > maxYCoord) maxYCoord = bottomY;
          }
        }
        // find difference between current and original center
      } catch (err) {
        _didIteratorError = true;
        _iteratorError = err;
      } finally {
        try {
          if (!_iteratorNormalCompletion && _iterator.return) {
            _iterator.return();
          }
        } finally {
          if (_didIteratorError) {
            throw _iteratorError;
          }
        }
      }

      var diffOnX = originalCenter.x - (maxXCoord + minXCoord) / 2;
      var diffOnY = originalCenter.y - (maxYCoord + minYCoord) / 2;
      // move component to original center
      componentResult.xCoords = componentResult.xCoords.map(function (x) {
        return x + diffOnX;
      });
      componentResult.yCoords = componentResult.yCoords.map(function (y) {
        return y + diffOnY;
      });
    } else {
      // calculate current bounding box
      Object.keys(componentResult).forEach(function (item) {
        var node = componentResult[item];
        var leftX = node.getRect().x;
        var rightX = node.getRect().x + node.getRect().width;
        var topY = node.getRect().y;
        var bottomY = node.getRect().y + node.getRect().height;

        if (leftX < minXCoord) minXCoord = leftX;
        if (rightX > maxXCoord) maxXCoord = rightX;
        if (topY < minYCoord) minYCoord = topY;
        if (bottomY > maxYCoord) maxYCoord = bottomY;
      });
      // find difference between current and original center
      var _diffOnX = originalCenter.x - (maxXCoord + minXCoord) / 2;
      var _diffOnY = originalCenter.y - (maxYCoord + minYCoord) / 2;
      // move component to original center
      Object.keys(componentResult).forEach(function (item) {
        var node = componentResult[item];
        node.setCenter(node.getCenterX() + _diffOnX, node.getCenterY() + _diffOnY);
      });
    }
  }
};

auxiliary.calcBoundingBox = function (parentNode, xCoords, yCoords, nodeIndexes) {
  // calculate bounds
  var left = Number.MAX_SAFE_INTEGER;
  var right = Number.MIN_SAFE_INTEGER;
  var top = Number.MAX_SAFE_INTEGER;
  var bottom = Number.MIN_SAFE_INTEGER;
  var nodeLeft = void 0;
  var nodeRight = void 0;
  var nodeTop = void 0;
  var nodeBottom = void 0;

  var nodes = parentNode.descendants().not(":parent");
  var s = nodes.length;
  for (var i = 0; i < s; i++) {
    var node = nodes[i];

    nodeLeft = xCoords[nodeIndexes.get(node.id())] - node.width() / 2;
    nodeRight = xCoords[nodeIndexes.get(node.id())] + node.width() / 2;
    nodeTop = yCoords[nodeIndexes.get(node.id())] - node.height() / 2;
    nodeBottom = yCoords[nodeIndexes.get(node.id())] + node.height() / 2;

    if (left > nodeLeft) {
      left = nodeLeft;
    }

    if (right < nodeRight) {
      right = nodeRight;
    }

    if (top > nodeTop) {
      top = nodeTop;
    }

    if (bottom < nodeBottom) {
      bottom = nodeBottom;
    }
  }

  var boundingBox = {};
  boundingBox.topLeftX = left;
  boundingBox.topLeftY = top;
  boundingBox.width = right - left;
  boundingBox.height = bottom - top;
  return boundingBox;
};

// This function finds and returns parent nodes whose all children are hidden
auxiliary.calcParentsWithoutChildren = function (cy, eles) {
  var parentsWithoutChildren = cy.collection();
  eles.nodes(':parent').forEach(function (parent) {
    var check = false;
    parent.children().forEach(function (child) {
      if (child.css('display') != 'none') {
        check = true;
      }
    });
    if (!check) {
      parentsWithoutChildren.merge(parent);
    }
  });

  return parentsWithoutChildren;
};

module.exports = auxiliary;

/***/ }),

/***/ 816:
/***/ ((module, __unused_webpack_exports, __nested_webpack_require_10973__) => {



/**
  The implementation of the postprocessing part that applies CoSE layout over the spectral layout
*/

var aux = __nested_webpack_require_10973__(548);
var CoSELayout = __nested_webpack_require_10973__(140).CoSELayout;
var CoSENode = __nested_webpack_require_10973__(140).CoSENode;
var PointD = __nested_webpack_require_10973__(140).layoutBase.PointD;
var DimensionD = __nested_webpack_require_10973__(140).layoutBase.DimensionD;
var LayoutConstants = __nested_webpack_require_10973__(140).layoutBase.LayoutConstants;
var FDLayoutConstants = __nested_webpack_require_10973__(140).layoutBase.FDLayoutConstants;
var CoSEConstants = __nested_webpack_require_10973__(140).CoSEConstants;

// main function that cose layout is processed
var coseLayout = function coseLayout(options, spectralResult) {

  var cy = options.cy;
  var eles = options.eles;
  var nodes = eles.nodes();
  var edges = eles.edges();

  var nodeIndexes = void 0;
  var xCoords = void 0;
  var yCoords = void 0;
  var idToLNode = {};

  if (options.randomize) {
    nodeIndexes = spectralResult["nodeIndexes"];
    xCoords = spectralResult["xCoords"];
    yCoords = spectralResult["yCoords"];
  }

  var isFn = function isFn(fn) {
    return typeof fn === 'function';
  };

  var optFn = function optFn(opt, ele) {
    if (isFn(opt)) {
      return opt(ele);
    } else {
      return opt;
    }
  };

  /**** Postprocessing functions ****/

  var parentsWithoutChildren = aux.calcParentsWithoutChildren(cy, eles);

  // transfer cytoscape nodes to cose nodes
  var processChildrenList = function processChildrenList(parent, children, layout, options) {
    var size = children.length;
    for (var i = 0; i < size; i++) {
      var theChild = children[i];
      var children_of_children = null;
      if (theChild.intersection(parentsWithoutChildren).length == 0) {
        children_of_children = theChild.children();
      }
      var theNode = void 0;

      var dimensions = theChild.layoutDimensions({
        nodeDimensionsIncludeLabels: options.nodeDimensionsIncludeLabels
      });

      if (theChild.outerWidth() != null && theChild.outerHeight() != null) {
        if (options.randomize) {
          if (!theChild.isParent()) {
            theNode = parent.add(new CoSENode(layout.graphManager, new PointD(xCoords[nodeIndexes.get(theChild.id())] - dimensions.w / 2, yCoords[nodeIndexes.get(theChild.id())] - dimensions.h / 2), new DimensionD(parseFloat(dimensions.w), parseFloat(dimensions.h))));
          } else {
            var parentInfo = aux.calcBoundingBox(theChild, xCoords, yCoords, nodeIndexes);
            if (theChild.intersection(parentsWithoutChildren).length == 0) {
              theNode = parent.add(new CoSENode(layout.graphManager, new PointD(parentInfo.topLeftX, parentInfo.topLeftY), new DimensionD(parentInfo.width, parentInfo.height)));
            } else {
              // for the parentsWithoutChildren
              theNode = parent.add(new CoSENode(layout.graphManager, new PointD(parentInfo.topLeftX, parentInfo.topLeftY), new DimensionD(parseFloat(dimensions.w), parseFloat(dimensions.h))));
            }
          }
        } else {
          theNode = parent.add(new CoSENode(layout.graphManager, new PointD(theChild.position('x') - dimensions.w / 2, theChild.position('y') - dimensions.h / 2), new DimensionD(parseFloat(dimensions.w), parseFloat(dimensions.h))));
        }
      } else {
        theNode = parent.add(new CoSENode(this.graphManager));
      }
      // Attach id to the layout node and repulsion value
      theNode.id = theChild.data("id");
      theNode.nodeRepulsion = optFn(options.nodeRepulsion, theChild);
      // Attach the paddings of cy node to layout node
      theNode.paddingLeft = parseInt(theChild.css('padding'));
      theNode.paddingTop = parseInt(theChild.css('padding'));
      theNode.paddingRight = parseInt(theChild.css('padding'));
      theNode.paddingBottom = parseInt(theChild.css('padding'));

      //Attach the label properties to both compound and simple nodes if labels will be included in node dimensions
      //These properties will be used while updating bounds of compounds during iterations or tiling
      //and will be used for simple nodes while transferring final positions to cytoscape
      if (options.nodeDimensionsIncludeLabels) {
        theNode.labelWidth = theChild.boundingBox({ includeLabels: true, includeNodes: false, includeOverlays: false }).w;
        theNode.labelHeight = theChild.boundingBox({ includeLabels: true, includeNodes: false, includeOverlays: false }).h;
        theNode.labelPosVertical = theChild.css("text-valign");
        theNode.labelPosHorizontal = theChild.css("text-halign");
      }

      // Map the layout node
      idToLNode[theChild.data("id")] = theNode;

      if (isNaN(theNode.rect.x)) {
        theNode.rect.x = 0;
      }

      if (isNaN(theNode.rect.y)) {
        theNode.rect.y = 0;
      }

      if (children_of_children != null && children_of_children.length > 0) {
        var theNewGraph = void 0;
        theNewGraph = layout.getGraphManager().add(layout.newGraph(), theNode);
        processChildrenList(theNewGraph, children_of_children, layout, options);
      }
    }
  };

  // transfer cytoscape edges to cose edges
  var processEdges = function processEdges(layout, gm, edges) {
    var idealLengthTotal = 0;
    var edgeCount = 0;
    for (var i = 0; i < edges.length; i++) {
      var edge = edges[i];
      var sourceNode = idToLNode[edge.data("source")];
      var targetNode = idToLNode[edge.data("target")];
      if (sourceNode && targetNode && sourceNode !== targetNode && sourceNode.getEdgesBetween(targetNode).length == 0) {
        var e1 = gm.add(layout.newEdge(), sourceNode, targetNode);
        e1.id = edge.id();
        e1.idealLength = optFn(options.idealEdgeLength, edge);
        e1.edgeElasticity = optFn(options.edgeElasticity, edge);
        idealLengthTotal += e1.idealLength;
        edgeCount++;
      }
    }
    // we need to update the ideal edge length constant with the avg. ideal length value after processing edges
    // in case there is no edge, use other options
    if (options.idealEdgeLength != null) {
      if (edgeCount > 0) CoSEConstants.DEFAULT_EDGE_LENGTH = FDLayoutConstants.DEFAULT_EDGE_LENGTH = idealLengthTotal / edgeCount;else if (!isFn(options.idealEdgeLength)) // in case there is no edge, but option gives a value to use
        CoSEConstants.DEFAULT_EDGE_LENGTH = FDLayoutConstants.DEFAULT_EDGE_LENGTH = options.idealEdgeLength;else // in case there is no edge and we cannot get a value from option (because it's a function)
        CoSEConstants.DEFAULT_EDGE_LENGTH = FDLayoutConstants.DEFAULT_EDGE_LENGTH = 50;
      // we need to update these constant values based on the ideal edge length constant
      CoSEConstants.MIN_REPULSION_DIST = FDLayoutConstants.MIN_REPULSION_DIST = FDLayoutConstants.DEFAULT_EDGE_LENGTH / 10.0;
      CoSEConstants.DEFAULT_RADIAL_SEPARATION = FDLayoutConstants.DEFAULT_EDGE_LENGTH;
    }
  };

  // transfer cytoscape constraints to cose layout
  var processConstraints = function processConstraints(layout, options) {
    // get nodes to be fixed
    if (options.fixedNodeConstraint) {
      layout.constraints["fixedNodeConstraint"] = options.fixedNodeConstraint;
    }
    // get nodes to be aligned
    if (options.alignmentConstraint) {
      layout.constraints["alignmentConstraint"] = options.alignmentConstraint;
    }
    // get nodes to be relatively placed
    if (options.relativePlacementConstraint) {
      layout.constraints["relativePlacementConstraint"] = options.relativePlacementConstraint;
    }
  };

  /**** Apply postprocessing ****/
  if (options.nestingFactor != null) CoSEConstants.PER_LEVEL_IDEAL_EDGE_LENGTH_FACTOR = FDLayoutConstants.PER_LEVEL_IDEAL_EDGE_LENGTH_FACTOR = options.nestingFactor;
  if (options.gravity != null) CoSEConstants.DEFAULT_GRAVITY_STRENGTH = FDLayoutConstants.DEFAULT_GRAVITY_STRENGTH = options.gravity;
  if (options.numIter != null) CoSEConstants.MAX_ITERATIONS = FDLayoutConstants.MAX_ITERATIONS = options.numIter;
  if (options.gravityRange != null) CoSEConstants.DEFAULT_GRAVITY_RANGE_FACTOR = FDLayoutConstants.DEFAULT_GRAVITY_RANGE_FACTOR = options.gravityRange;
  if (options.gravityCompound != null) CoSEConstants.DEFAULT_COMPOUND_GRAVITY_STRENGTH = FDLayoutConstants.DEFAULT_COMPOUND_GRAVITY_STRENGTH = options.gravityCompound;
  if (options.gravityRangeCompound != null) CoSEConstants.DEFAULT_COMPOUND_GRAVITY_RANGE_FACTOR = FDLayoutConstants.DEFAULT_COMPOUND_GRAVITY_RANGE_FACTOR = options.gravityRangeCompound;
  if (options.initialEnergyOnIncremental != null) CoSEConstants.DEFAULT_COOLING_FACTOR_INCREMENTAL = FDLayoutConstants.DEFAULT_COOLING_FACTOR_INCREMENTAL = options.initialEnergyOnIncremental;

  if (options.tilingCompareBy != null) CoSEConstants.TILING_COMPARE_BY = options.tilingCompareBy;

  if (options.quality == 'proof') LayoutConstants.QUALITY = 2;else LayoutConstants.QUALITY = 0;

  CoSEConstants.NODE_DIMENSIONS_INCLUDE_LABELS = FDLayoutConstants.NODE_DIMENSIONS_INCLUDE_LABELS = LayoutConstants.NODE_DIMENSIONS_INCLUDE_LABELS = options.nodeDimensionsIncludeLabels;
  CoSEConstants.DEFAULT_INCREMENTAL = FDLayoutConstants.DEFAULT_INCREMENTAL = LayoutConstants.DEFAULT_INCREMENTAL = !options.randomize;
  CoSEConstants.ANIMATE = FDLayoutConstants.ANIMATE = LayoutConstants.ANIMATE = options.animate;
  CoSEConstants.TILE = options.tile;
  CoSEConstants.TILING_PADDING_VERTICAL = typeof options.tilingPaddingVertical === 'function' ? options.tilingPaddingVertical.call() : options.tilingPaddingVertical;
  CoSEConstants.TILING_PADDING_HORIZONTAL = typeof options.tilingPaddingHorizontal === 'function' ? options.tilingPaddingHorizontal.call() : options.tilingPaddingHorizontal;

  CoSEConstants.DEFAULT_INCREMENTAL = FDLayoutConstants.DEFAULT_INCREMENTAL = LayoutConstants.DEFAULT_INCREMENTAL = true;
  CoSEConstants.PURE_INCREMENTAL = !options.randomize;
  LayoutConstants.DEFAULT_UNIFORM_LEAF_NODE_SIZES = options.uniformNodeDimensions;

  // This part is for debug/demo purpose
  if (options.step == "transformed") {
    CoSEConstants.TRANSFORM_ON_CONSTRAINT_HANDLING = true;
    CoSEConstants.ENFORCE_CONSTRAINTS = false;
    CoSEConstants.APPLY_LAYOUT = false;
  }
  if (options.step == "enforced") {
    CoSEConstants.TRANSFORM_ON_CONSTRAINT_HANDLING = false;
    CoSEConstants.ENFORCE_CONSTRAINTS = true;
    CoSEConstants.APPLY_LAYOUT = false;
  }
  if (options.step == "cose") {
    CoSEConstants.TRANSFORM_ON_CONSTRAINT_HANDLING = false;
    CoSEConstants.ENFORCE_CONSTRAINTS = false;
    CoSEConstants.APPLY_LAYOUT = true;
  }
  if (options.step == "all") {
    if (options.randomize) CoSEConstants.TRANSFORM_ON_CONSTRAINT_HANDLING = true;else CoSEConstants.TRANSFORM_ON_CONSTRAINT_HANDLING = false;
    CoSEConstants.ENFORCE_CONSTRAINTS = true;
    CoSEConstants.APPLY_LAYOUT = true;
  }

  if (options.fixedNodeConstraint || options.alignmentConstraint || options.relativePlacementConstraint) {
    CoSEConstants.TREE_REDUCTION_ON_INCREMENTAL = false;
  } else {
    CoSEConstants.TREE_REDUCTION_ON_INCREMENTAL = true;
  }

  var coseLayout = new CoSELayout();
  var gm = coseLayout.newGraphManager();

  processChildrenList(gm.addRoot(), aux.getTopMostNodes(nodes), coseLayout, options);
  processEdges(coseLayout, gm, edges);
  processConstraints(coseLayout, options);

  coseLayout.runLayout();

  return idToLNode;
};

module.exports = { coseLayout: coseLayout };

/***/ }),

/***/ 212:
/***/ ((module, __unused_webpack_exports, __nested_webpack_require_22455__) => {



var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

/**
  The implementation of the fcose layout algorithm
*/

var assign = __nested_webpack_require_22455__(658);
var aux = __nested_webpack_require_22455__(548);

var _require = __nested_webpack_require_22455__(657),
    spectralLayout = _require.spectralLayout;

var _require2 = __nested_webpack_require_22455__(816),
    coseLayout = _require2.coseLayout;

var defaults = Object.freeze({

  // 'draft', 'default' or 'proof' 
  // - 'draft' only applies spectral layout 
  // - 'default' improves the quality with subsequent CoSE layout (fast cooling rate)
  // - 'proof' improves the quality with subsequent CoSE layout (slow cooling rate) 
  quality: "default",
  // Use random node positions at beginning of layout
  // if this is set to false, then quality option must be "proof"
  randomize: true,
  // Whether or not to animate the layout
  animate: true,
  // Duration of animation in ms, if enabled
  animationDuration: 1000,
  // Easing of animation, if enabled
  animationEasing: undefined,
  // Fit the viewport to the repositioned nodes
  fit: true,
  // Padding around layout
  padding: 30,
  // Whether to include labels in node dimensions. Valid in "proof" quality
  nodeDimensionsIncludeLabels: false,
  // Whether or not simple nodes (non-compound nodes) are of uniform dimensions
  uniformNodeDimensions: false,
  // Whether to pack disconnected components - valid only if randomize: true
  packComponents: true,
  // Layout step - all, transformed, enforced, cose - for debug purpose only
  step: "all",

  /* spectral layout options */

  // False for random, true for greedy
  samplingType: true,
  // Sample size to construct distance matrix
  sampleSize: 25,
  // Separation amount between nodes
  nodeSeparation: 75,
  // Power iteration tolerance
  piTol: 0.0000001,

  /* CoSE layout options */

  // Node repulsion (non overlapping) multiplier
  nodeRepulsion: function nodeRepulsion(node) {
    return 4500;
  },
  // Ideal edge (non nested) length
  idealEdgeLength: function idealEdgeLength(edge) {
    return 50;
  },
  // Divisor to compute edge forces
  edgeElasticity: function edgeElasticity(edge) {
    return 0.45;
  },
  // Nesting factor (multiplier) to compute ideal edge length for nested edges
  nestingFactor: 0.1,
  // Gravity force (constant)
  gravity: 0.25,
  // Maximum number of iterations to perform
  numIter: 2500,
  // For enabling tiling
  tile: true,
  // The function that specifies the criteria for comparing nodes while sorting them during tiling operation.
  // Takes the node id as a parameter and the default tiling operation is perfomed when this option is not set.
  tilingCompareBy: undefined,
  // Represents the amount of the vertical space to put between the zero degree members during the tiling operation(can also be a function)
  tilingPaddingVertical: 10,
  // Represents the amount of the horizontal space to put between the zero degree members during the tiling operation(can also be a function)
  tilingPaddingHorizontal: 10,
  // Gravity range (constant) for compounds
  gravityRangeCompound: 1.5,
  // Gravity force (constant) for compounds
  gravityCompound: 1.0,
  // Gravity range (constant)
  gravityRange: 3.8,
  // Initial cooling factor for incremental layout  
  initialEnergyOnIncremental: 0.3,

  /* constraint options */

  // Fix required nodes to predefined positions
  // [{nodeId: 'n1', position: {x: 100, y: 200}, {...}]
  fixedNodeConstraint: undefined,
  // Align required nodes in vertical/horizontal direction
  // {vertical: [['n1', 'n2')], ['n3', 'n4']], horizontal: ['n2', 'n4']}
  alignmentConstraint: undefined,
  // Place two nodes relatively in vertical/horizontal direction 
  // [{top: 'n1', bottom: 'n2', gap: 100}, {left: 'n3', right: 'n4', gap: 75}]
  relativePlacementConstraint: undefined,

  /* layout event callbacks */
  ready: function ready() {}, // on layoutready
  stop: function stop() {} // on layoutstop
});

var Layout = function () {
  function Layout(options) {
    _classCallCheck(this, Layout);

    this.options = assign({}, defaults, options);
  }

  _createClass(Layout, [{
    key: 'run',
    value: function run() {
      var layout = this;
      var options = this.options;
      var cy = options.cy;
      var eles = options.eles;

      var spectralResult = [];
      var xCoords = void 0;
      var yCoords = void 0;
      var coseResult = [];
      var components = void 0;
      var componentCenters = [];

      // basic validity check for constraint inputs 
      if (options.fixedNodeConstraint && (!Array.isArray(options.fixedNodeConstraint) || options.fixedNodeConstraint.length == 0)) {
        options.fixedNodeConstraint = undefined;
      }

      if (options.alignmentConstraint) {
        if (options.alignmentConstraint.vertical && (!Array.isArray(options.alignmentConstraint.vertical) || options.alignmentConstraint.vertical.length == 0)) {
          options.alignmentConstraint.vertical = undefined;
        }
        if (options.alignmentConstraint.horizontal && (!Array.isArray(options.alignmentConstraint.horizontal) || options.alignmentConstraint.horizontal.length == 0)) {
          options.alignmentConstraint.horizontal = undefined;
        }
      }

      if (options.relativePlacementConstraint && (!Array.isArray(options.relativePlacementConstraint) || options.relativePlacementConstraint.length == 0)) {
        options.relativePlacementConstraint = undefined;
      }

      // if any constraint exists, set some options
      var constraintExist = options.fixedNodeConstraint || options.alignmentConstraint || options.relativePlacementConstraint;
      if (constraintExist) {
        // constraints work with these options
        options.tile = false;
        options.packComponents = false;
      }

      // decide component packing is enabled or not
      var layUtil = void 0;
      var packingEnabled = false;
      if (cy.layoutUtilities && options.packComponents) {
        layUtil = cy.layoutUtilities("get");
        if (!layUtil) layUtil = cy.layoutUtilities();
        packingEnabled = true;
      }

      if (eles.nodes().length > 0) {
        // if packing is not enabled, perform layout on the whole graph
        if (!packingEnabled) {
          // store component center
          var boundingBox = options.eles.boundingBox();
          componentCenters.push({ x: boundingBox.x1 + boundingBox.w / 2, y: boundingBox.y1 + boundingBox.h / 2 });
          // apply spectral layout
          if (options.randomize) {
            var result = spectralLayout(options);
            spectralResult.push(result);
          }
          // apply cose layout as postprocessing
          if (options.quality == "default" || options.quality == "proof") {
            coseResult.push(coseLayout(options, spectralResult[0]));
            aux.relocateComponent(componentCenters[0], coseResult[0], options); // relocate center to original position
          } else {
            aux.relocateComponent(componentCenters[0], spectralResult[0], options); // relocate center to original position
          }
        } else {
          // packing is enabled
          var topMostNodes = aux.getTopMostNodes(options.eles.nodes());
          components = aux.connectComponents(cy, options.eles, topMostNodes);
          // store component centers
          components.forEach(function (component) {
            var boundingBox = component.boundingBox();
            componentCenters.push({ x: boundingBox.x1 + boundingBox.w / 2, y: boundingBox.y1 + boundingBox.h / 2 });
          });

          //send each component to spectral layout if randomized
          if (options.randomize) {
            components.forEach(function (component) {
              options.eles = component;
              spectralResult.push(spectralLayout(options));
            });
          }

          if (options.quality == "default" || options.quality == "proof") {
            var toBeTiledNodes = cy.collection();
            if (options.tile) {
              // behave nodes to be tiled as one component
              var nodeIndexes = new Map();
              var _xCoords = [];
              var _yCoords = [];
              var count = 0;
              var tempSpectralResult = { nodeIndexes: nodeIndexes, xCoords: _xCoords, yCoords: _yCoords };
              var indexesToBeDeleted = [];
              components.forEach(function (component, index) {
                if (component.edges().length == 0) {
                  component.nodes().forEach(function (node, i) {
                    toBeTiledNodes.merge(component.nodes()[i]);
                    if (!node.isParent()) {
                      tempSpectralResult.nodeIndexes.set(component.nodes()[i].id(), count++);
                      tempSpectralResult.xCoords.push(component.nodes()[0].position().x);
                      tempSpectralResult.yCoords.push(component.nodes()[0].position().y);
                    }
                  });
                  indexesToBeDeleted.push(index);
                }
              });
              if (toBeTiledNodes.length > 1) {
                var _boundingBox = toBeTiledNodes.boundingBox();
                componentCenters.push({ x: _boundingBox.x1 + _boundingBox.w / 2, y: _boundingBox.y1 + _boundingBox.h / 2 });
                components.push(toBeTiledNodes);
                spectralResult.push(tempSpectralResult);
                for (var i = indexesToBeDeleted.length - 1; i >= 0; i--) {
                  components.splice(indexesToBeDeleted[i], 1);
                  spectralResult.splice(indexesToBeDeleted[i], 1);
                  componentCenters.splice(indexesToBeDeleted[i], 1);
                };
              }
            }
            components.forEach(function (component, index) {
              // send each component to cose layout
              options.eles = component;
              coseResult.push(coseLayout(options, spectralResult[index]));
              aux.relocateComponent(componentCenters[index], coseResult[index], options); // relocate center to original position
            });
          } else {
            components.forEach(function (component, index) {
              aux.relocateComponent(componentCenters[index], spectralResult[index], options); // relocate center to original position
            });
          }

          // packing
          var componentsEvaluated = new Set();
          if (components.length > 1) {
            var subgraphs = [];
            var hiddenEles = eles.filter(function (ele) {
              return ele.css('display') == 'none';
            });
            components.forEach(function (component, index) {
              var nodeIndexes = void 0;
              if (options.quality == "draft") {
                nodeIndexes = spectralResult[index].nodeIndexes;
              }

              if (component.nodes().not(hiddenEles).length > 0) {
                var subgraph = {};
                subgraph.edges = [];
                subgraph.nodes = [];
                var nodeIndex = void 0;
                component.nodes().not(hiddenEles).forEach(function (node) {
                  if (options.quality == "draft") {
                    if (!node.isParent()) {
                      nodeIndex = nodeIndexes.get(node.id());
                      subgraph.nodes.push({ x: spectralResult[index].xCoords[nodeIndex] - node.boundingbox().w / 2, y: spectralResult[index].yCoords[nodeIndex] - node.boundingbox().h / 2, width: node.boundingbox().w, height: node.boundingbox().h });
                    } else {
                      var parentInfo = aux.calcBoundingBox(node, spectralResult[index].xCoords, spectralResult[index].yCoords, nodeIndexes);
                      subgraph.nodes.push({ x: parentInfo.topLeftX, y: parentInfo.topLeftY, width: parentInfo.width, height: parentInfo.height });
                    }
                  } else {
                    if (coseResult[index][node.id()]) {
                      subgraph.nodes.push({ x: coseResult[index][node.id()].getLeft(), y: coseResult[index][node.id()].getTop(), width: coseResult[index][node.id()].getWidth(), height: coseResult[index][node.id()].getHeight() });
                    }
                  }
                });
                component.edges().forEach(function (edge) {
                  var source = edge.source();
                  var target = edge.target();
                  if (source.css("display") != "none" && target.css("display") != "none") {
                    if (options.quality == "draft") {
                      var sourceNodeIndex = nodeIndexes.get(source.id());
                      var targetNodeIndex = nodeIndexes.get(target.id());
                      var sourceCenter = [];
                      var targetCenter = [];
                      if (source.isParent()) {
                        var parentInfo = aux.calcBoundingBox(source, spectralResult[index].xCoords, spectralResult[index].yCoords, nodeIndexes);
                        sourceCenter.push(parentInfo.topLeftX + parentInfo.width / 2);
                        sourceCenter.push(parentInfo.topLeftY + parentInfo.height / 2);
                      } else {
                        sourceCenter.push(spectralResult[index].xCoords[sourceNodeIndex]);
                        sourceCenter.push(spectralResult[index].yCoords[sourceNodeIndex]);
                      }
                      if (target.isParent()) {
                        var _parentInfo = aux.calcBoundingBox(target, spectralResult[index].xCoords, spectralResult[index].yCoords, nodeIndexes);
                        targetCenter.push(_parentInfo.topLeftX + _parentInfo.width / 2);
                        targetCenter.push(_parentInfo.topLeftY + _parentInfo.height / 2);
                      } else {
                        targetCenter.push(spectralResult[index].xCoords[targetNodeIndex]);
                        targetCenter.push(spectralResult[index].yCoords[targetNodeIndex]);
                      }
                      subgraph.edges.push({ startX: sourceCenter[0], startY: sourceCenter[1], endX: targetCenter[0], endY: targetCenter[1] });
                    } else {
                      if (coseResult[index][source.id()] && coseResult[index][target.id()]) {
                        subgraph.edges.push({ startX: coseResult[index][source.id()].getCenterX(), startY: coseResult[index][source.id()].getCenterY(), endX: coseResult[index][target.id()].getCenterX(), endY: coseResult[index][target.id()].getCenterY() });
                      }
                    }
                  }
                });
                if (subgraph.nodes.length > 0) {
                  subgraphs.push(subgraph);
                  componentsEvaluated.add(index);
                }
              }
            });
            var shiftResult = layUtil.packComponents(subgraphs, options.randomize).shifts;
            if (options.quality == "draft") {
              spectralResult.forEach(function (result, index) {
                var newXCoords = result.xCoords.map(function (x) {
                  return x + shiftResult[index].dx;
                });
                var newYCoords = result.yCoords.map(function (y) {
                  return y + shiftResult[index].dy;
                });
                result.xCoords = newXCoords;
                result.yCoords = newYCoords;
              });
            } else {
              var _count = 0;
              componentsEvaluated.forEach(function (index) {
                Object.keys(coseResult[index]).forEach(function (item) {
                  var nodeRectangle = coseResult[index][item];
                  nodeRectangle.setCenter(nodeRectangle.getCenterX() + shiftResult[_count].dx, nodeRectangle.getCenterY() + shiftResult[_count].dy);
                });
                _count++;
              });
            }
          }
        }
      }

      // get each element's calculated position
      var getPositions = function getPositions(ele, i) {
        if (options.quality == "default" || options.quality == "proof") {
          if (typeof ele === "number") {
            ele = i;
          }
          var pos = void 0;
          var node = void 0;
          var theId = ele.data('id');
          coseResult.forEach(function (result) {
            if (theId in result) {
              pos = { x: result[theId].getRect().getCenterX(), y: result[theId].getRect().getCenterY() };
              node = result[theId];
            }
          });
          if (options.nodeDimensionsIncludeLabels) {
            if (node.labelWidth) {
              if (node.labelPosHorizontal == "left") {
                pos.x += node.labelWidth / 2;
              } else if (node.labelPosHorizontal == "right") {
                pos.x -= node.labelWidth / 2;
              }
            }
            if (node.labelHeight) {
              if (node.labelPosVertical == "top") {
                pos.y += node.labelHeight / 2;
              } else if (node.labelPosVertical == "bottom") {
                pos.y -= node.labelHeight / 2;
              }
            }
          }
          if (pos == undefined) pos = { x: ele.position("x"), y: ele.position("y") };
          return {
            x: pos.x,
            y: pos.y
          };
        } else {
          var _pos = void 0;
          spectralResult.forEach(function (result) {
            var index = result.nodeIndexes.get(ele.id());
            if (index != undefined) {
              _pos = { x: result.xCoords[index], y: result.yCoords[index] };
            }
          });
          if (_pos == undefined) _pos = { x: ele.position("x"), y: ele.position("y") };
          return {
            x: _pos.x,
            y: _pos.y
          };
        }
      };

      // quality = "draft" and randomize = false are contradictive so in that case positions don't change
      if (options.quality == "default" || options.quality == "proof" || options.randomize) {
        // transfer calculated positions to nodes (positions of only simple nodes are evaluated, compounds are positioned automatically)
        var parentsWithoutChildren = aux.calcParentsWithoutChildren(cy, eles);
        var _hiddenEles = eles.filter(function (ele) {
          return ele.css('display') == 'none';
        });
        options.eles = eles.not(_hiddenEles);

        eles.nodes().not(":parent").not(_hiddenEles).layoutPositions(layout, options, getPositions);

        if (parentsWithoutChildren.length > 0) {
          parentsWithoutChildren.forEach(function (ele) {
            ele.position(getPositions(ele));
          });
        }
      } else {
        console.log("If randomize option is set to false, then quality option must be 'default' or 'proof'.");
      }
    }
  }]);

  return Layout;
}();

module.exports = Layout;

/***/ }),

/***/ 657:
/***/ ((module, __unused_webpack_exports, __nested_webpack_require_41986__) => {



/**
  The implementation of the spectral layout that is the first part of the fcose layout algorithm
*/

var aux = __nested_webpack_require_41986__(548);
var Matrix = __nested_webpack_require_41986__(140).layoutBase.Matrix;
var SVD = __nested_webpack_require_41986__(140).layoutBase.SVD;

// main function that spectral layout is processed
var spectralLayout = function spectralLayout(options) {

  var cy = options.cy;
  var eles = options.eles;
  var nodes = eles.nodes();
  var parentNodes = eles.nodes(":parent");

  var dummyNodes = new Map(); // map to keep dummy nodes and their neighbors
  var nodeIndexes = new Map(); // map to keep indexes to nodes
  var parentChildMap = new Map(); // mapping btw. compound and its representative node 
  var allNodesNeighborhood = []; // array to keep neighborhood of all nodes
  var xCoords = [];
  var yCoords = [];

  var samplesColumn = []; // sampled vertices
  var minDistancesColumn = [];
  var C = []; // column sampling matrix
  var PHI = []; // intersection of column and row sampling matrices 
  var INV = []; // inverse of PHI 

  var firstSample = void 0; // the first sampled node
  var nodeSize = void 0;

  var infinity = 100000000;
  var small = 0.000000001;

  var piTol = options.piTol;
  var samplingType = options.samplingType; // false for random, true for greedy
  var nodeSeparation = options.nodeSeparation;
  var sampleSize = void 0;

  /**** Spectral-preprocessing functions ****/

  /**** Spectral layout functions ****/

  // determine which columns to be sampled
  var randomSampleCR = function randomSampleCR() {
    var sample = 0;
    var count = 0;
    var flag = false;

    while (count < sampleSize) {
      sample = Math.floor(Math.random() * nodeSize);

      flag = false;
      for (var i = 0; i < count; i++) {
        if (samplesColumn[i] == sample) {
          flag = true;
          break;
        }
      }

      if (!flag) {
        samplesColumn[count] = sample;
        count++;
      } else {
        continue;
      }
    }
  };

  // takes the index of the node(pivot) to initiate BFS as a parameter
  var BFS = function BFS(pivot, index, samplingMethod) {
    var path = []; // the front of the path
    var front = 0; // the back of the path
    var back = 0;
    var current = 0;
    var temp = void 0;
    var distance = [];

    var max_dist = 0; // the furthest node to be returned
    var max_ind = 1;

    for (var i = 0; i < nodeSize; i++) {
      distance[i] = infinity;
    }

    path[back] = pivot;
    distance[pivot] = 0;

    while (back >= front) {
      current = path[front++];
      var neighbors = allNodesNeighborhood[current];
      for (var _i = 0; _i < neighbors.length; _i++) {
        temp = nodeIndexes.get(neighbors[_i]);
        if (distance[temp] == infinity) {
          distance[temp] = distance[current] + 1;
          path[++back] = temp;
        }
      }
      C[current][index] = distance[current] * nodeSeparation;
    }

    if (samplingMethod) {
      for (var _i2 = 0; _i2 < nodeSize; _i2++) {
        if (C[_i2][index] < minDistancesColumn[_i2]) minDistancesColumn[_i2] = C[_i2][index];
      }

      for (var _i3 = 0; _i3 < nodeSize; _i3++) {
        if (minDistancesColumn[_i3] > max_dist) {
          max_dist = minDistancesColumn[_i3];
          max_ind = _i3;
        }
      }
    }
    return max_ind;
  };

  // apply BFS to all nodes or selected samples
  var allBFS = function allBFS(samplingMethod) {

    var sample = void 0;

    if (!samplingMethod) {
      randomSampleCR();

      // call BFS
      for (var i = 0; i < sampleSize; i++) {
        BFS(samplesColumn[i], i, samplingMethod, false);
      }
    } else {
      sample = Math.floor(Math.random() * nodeSize);
      firstSample = sample;

      for (var _i4 = 0; _i4 < nodeSize; _i4++) {
        minDistancesColumn[_i4] = infinity;
      }

      for (var _i5 = 0; _i5 < sampleSize; _i5++) {
        samplesColumn[_i5] = sample;
        sample = BFS(sample, _i5, samplingMethod);
      }
    }

    // form the squared distances for C
    for (var _i6 = 0; _i6 < nodeSize; _i6++) {
      for (var j = 0; j < sampleSize; j++) {
        C[_i6][j] *= C[_i6][j];
      }
    }

    // form PHI
    for (var _i7 = 0; _i7 < sampleSize; _i7++) {
      PHI[_i7] = [];
    }

    for (var _i8 = 0; _i8 < sampleSize; _i8++) {
      for (var _j = 0; _j < sampleSize; _j++) {
        PHI[_i8][_j] = C[samplesColumn[_j]][_i8];
      }
    }
  };

  // perform the SVD algorithm and apply a regularization step
  var sample = function sample() {

    var SVDResult = SVD.svd(PHI);

    var a_q = SVDResult.S;
    var a_u = SVDResult.U;
    var a_v = SVDResult.V;

    var max_s = a_q[0] * a_q[0] * a_q[0];

    var a_Sig = [];

    //  regularization
    for (var i = 0; i < sampleSize; i++) {
      a_Sig[i] = [];
      for (var j = 0; j < sampleSize; j++) {
        a_Sig[i][j] = 0;
        if (i == j) {
          a_Sig[i][j] = a_q[i] / (a_q[i] * a_q[i] + max_s / (a_q[i] * a_q[i]));
        }
      }
    }

    INV = Matrix.multMat(Matrix.multMat(a_v, a_Sig), Matrix.transpose(a_u));
  };

  // calculate final coordinates 
  var powerIteration = function powerIteration() {
    // two largest eigenvalues
    var theta1 = void 0;
    var theta2 = void 0;

    // initial guesses for eigenvectors
    var Y1 = [];
    var Y2 = [];

    var V1 = [];
    var V2 = [];

    for (var i = 0; i < nodeSize; i++) {
      Y1[i] = Math.random();
      Y2[i] = Math.random();
    }

    Y1 = Matrix.normalize(Y1);
    Y2 = Matrix.normalize(Y2);

    var count = 0;
    // to keep track of the improvement ratio in power iteration
    var current = small;
    var previous = small;

    var temp = void 0;

    while (true) {
      count++;

      for (var _i9 = 0; _i9 < nodeSize; _i9++) {
        V1[_i9] = Y1[_i9];
      }

      Y1 = Matrix.multGamma(Matrix.multL(Matrix.multGamma(V1), C, INV));
      theta1 = Matrix.dotProduct(V1, Y1);
      Y1 = Matrix.normalize(Y1);

      current = Matrix.dotProduct(V1, Y1);

      temp = Math.abs(current / previous);

      if (temp <= 1 + piTol && temp >= 1) {
        break;
      }

      previous = current;
    }

    for (var _i10 = 0; _i10 < nodeSize; _i10++) {
      V1[_i10] = Y1[_i10];
    }

    count = 0;
    previous = small;
    while (true) {
      count++;

      for (var _i11 = 0; _i11 < nodeSize; _i11++) {
        V2[_i11] = Y2[_i11];
      }

      V2 = Matrix.minusOp(V2, Matrix.multCons(V1, Matrix.dotProduct(V1, V2)));
      Y2 = Matrix.multGamma(Matrix.multL(Matrix.multGamma(V2), C, INV));
      theta2 = Matrix.dotProduct(V2, Y2);
      Y2 = Matrix.normalize(Y2);

      current = Matrix.dotProduct(V2, Y2);

      temp = Math.abs(current / previous);

      if (temp <= 1 + piTol && temp >= 1) {
        break;
      }

      previous = current;
    }

    for (var _i12 = 0; _i12 < nodeSize; _i12++) {
      V2[_i12] = Y2[_i12];
    }

    // theta1 now contains dominant eigenvalue
    // theta2 now contains the second-largest eigenvalue
    // V1 now contains theta1's eigenvector
    // V2 now contains theta2's eigenvector

    //populate the two vectors
    xCoords = Matrix.multCons(V1, Math.sqrt(Math.abs(theta1)));
    yCoords = Matrix.multCons(V2, Math.sqrt(Math.abs(theta2)));
  };

  /**** Preparation for spectral layout (Preprocessing) ****/

  // connect disconnected components (first top level, then inside of each compound node)
  aux.connectComponents(cy, eles, aux.getTopMostNodes(nodes), dummyNodes);

  parentNodes.forEach(function (ele) {
    aux.connectComponents(cy, eles, aux.getTopMostNodes(ele.descendants().intersection(eles)), dummyNodes);
  });

  // assign indexes to nodes (first real, then dummy nodes)
  var index = 0;
  for (var i = 0; i < nodes.length; i++) {
    if (!nodes[i].isParent()) {
      nodeIndexes.set(nodes[i].id(), index++);
    }
  }

  var _iteratorNormalCompletion = true;
  var _didIteratorError = false;
  var _iteratorError = undefined;

  try {
    for (var _iterator = dummyNodes.keys()[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
      var key = _step.value;

      nodeIndexes.set(key, index++);
    }

    // instantiate the neighborhood matrix
  } catch (err) {
    _didIteratorError = true;
    _iteratorError = err;
  } finally {
    try {
      if (!_iteratorNormalCompletion && _iterator.return) {
        _iterator.return();
      }
    } finally {
      if (_didIteratorError) {
        throw _iteratorError;
      }
    }
  }

  for (var _i13 = 0; _i13 < nodeIndexes.size; _i13++) {
    allNodesNeighborhood[_i13] = [];
  }

  // form a parent-child map to keep representative node of each compound node  
  parentNodes.forEach(function (ele) {
    var children = ele.children().intersection(eles);

    //      let random = 0;
    while (children.nodes(":childless").length == 0) {
      //        random = Math.floor(Math.random() * children.nodes().length); // if all children are compound then proceed randomly
      children = children.nodes()[0].children().intersection(eles);
    }
    //  select the representative node - we can apply different methods here
    //      random = Math.floor(Math.random() * children.nodes(":childless").length);
    var index = 0;
    var min = children.nodes(":childless")[0].connectedEdges().length;
    children.nodes(":childless").forEach(function (ele2, i) {
      if (ele2.connectedEdges().length < min) {
        min = ele2.connectedEdges().length;
        index = i;
      }
    });
    parentChildMap.set(ele.id(), children.nodes(":childless")[index].id());
  });

  // add neighborhood relations (first real, then dummy nodes)
  nodes.forEach(function (ele) {
    var eleIndex = void 0;

    if (ele.isParent()) eleIndex = nodeIndexes.get(parentChildMap.get(ele.id()));else eleIndex = nodeIndexes.get(ele.id());

    ele.neighborhood().nodes().forEach(function (node) {
      if (eles.intersection(ele.edgesWith(node)).length > 0) {
        if (node.isParent()) allNodesNeighborhood[eleIndex].push(parentChildMap.get(node.id()));else allNodesNeighborhood[eleIndex].push(node.id());
      }
    });
  });

  var _loop = function _loop(_key) {
    var eleIndex = nodeIndexes.get(_key);
    var disconnectedId = void 0;
    dummyNodes.get(_key).forEach(function (id) {
      if (cy.getElementById(id).isParent()) disconnectedId = parentChildMap.get(id);else disconnectedId = id;

      allNodesNeighborhood[eleIndex].push(disconnectedId);
      allNodesNeighborhood[nodeIndexes.get(disconnectedId)].push(_key);
    });
  };

  var _iteratorNormalCompletion2 = true;
  var _didIteratorError2 = false;
  var _iteratorError2 = undefined;

  try {
    for (var _iterator2 = dummyNodes.keys()[Symbol.iterator](), _step2; !(_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done); _iteratorNormalCompletion2 = true) {
      var _key = _step2.value;

      _loop(_key);
    }

    // nodeSize now only considers the size of transformed graph
  } catch (err) {
    _didIteratorError2 = true;
    _iteratorError2 = err;
  } finally {
    try {
      if (!_iteratorNormalCompletion2 && _iterator2.return) {
        _iterator2.return();
      }
    } finally {
      if (_didIteratorError2) {
        throw _iteratorError2;
      }
    }
  }

  nodeSize = nodeIndexes.size;

  var spectralResult = void 0;

  // If number of nodes in transformed graph is 1 or 2, either SVD or powerIteration causes problem
  // So skip spectral and layout the graph with cose
  if (nodeSize > 2) {
    // if # of nodes in transformed graph is smaller than sample size,
    // then use # of nodes as sample size
    sampleSize = nodeSize < options.sampleSize ? nodeSize : options.sampleSize;

    // instantiates the partial matrices that will be used in spectral layout
    for (var _i14 = 0; _i14 < nodeSize; _i14++) {
      C[_i14] = [];
    }
    for (var _i15 = 0; _i15 < sampleSize; _i15++) {
      INV[_i15] = [];
    }

    /**** Apply spectral layout ****/

    if (options.quality == "draft" || options.step == "all") {
      allBFS(samplingType);
      sample();
      powerIteration();

      spectralResult = { nodeIndexes: nodeIndexes, xCoords: xCoords, yCoords: yCoords };
    } else {
      nodeIndexes.forEach(function (value, key) {
        xCoords.push(cy.getElementById(key).position("x"));
        yCoords.push(cy.getElementById(key).position("y"));
      });
      spectralResult = { nodeIndexes: nodeIndexes, xCoords: xCoords, yCoords: yCoords };
    }
    return spectralResult;
  } else {
    var iterator = nodeIndexes.keys();
    var firstNode = cy.getElementById(iterator.next().value);
    var firstNodePos = firstNode.position();
    var firstNodeWidth = firstNode.outerWidth();
    xCoords.push(firstNodePos.x);
    yCoords.push(firstNodePos.y);
    if (nodeSize == 2) {
      var secondNode = cy.getElementById(iterator.next().value);
      var secondNodeWidth = secondNode.outerWidth();
      xCoords.push(firstNodePos.x + firstNodeWidth / 2 + secondNodeWidth / 2 + options.idealEdgeLength);
      yCoords.push(firstNodePos.y);
    }

    spectralResult = { nodeIndexes: nodeIndexes, xCoords: xCoords, yCoords: yCoords };
    return spectralResult;
  }
};

module.exports = { spectralLayout: spectralLayout };

/***/ }),

/***/ 579:
/***/ ((module, __unused_webpack_exports, __nested_webpack_require_55337__) => {



var impl = __nested_webpack_require_55337__(212);

// registers the extension on a cytoscape lib ref
var register = function register(cytoscape) {
  if (!cytoscape) {
    return;
  } // can't register if cytoscape unspecified

  cytoscape('layout', 'fcose', impl); // register with cytoscape.js
};

if (typeof cytoscape !== 'undefined') {
  // expose to global cytoscape (i.e. window.cytoscape)
  register(cytoscape);
}

module.exports = register;

/***/ }),

/***/ 140:
/***/ ((module) => {

module.exports = __WEBPACK_EXTERNAL_MODULE__140__;

/***/ })

/******/ 	});
/************************************************************************/
/******/ 	// The module cache
/******/ 	var __webpack_module_cache__ = {};
/******/ 	
/******/ 	// The require function
/******/ 	function __nested_webpack_require_56127__(moduleId) {
/******/ 		// Check if module is in cache
/******/ 		var cachedModule = __webpack_module_cache__[moduleId];
/******/ 		if (cachedModule !== undefined) {
/******/ 			return cachedModule.exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = __webpack_module_cache__[moduleId] = {
/******/ 			// no module.id needed
/******/ 			// no module.loaded needed
/******/ 			exports: {}
/******/ 		};
/******/ 	
/******/ 		// Execute the module function
/******/ 		__webpack_modules__[moduleId](module, module.exports, __nested_webpack_require_56127__);
/******/ 	
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/ 	
/************************************************************************/
/******/ 	
/******/ 	// startup
/******/ 	// Load entry module and return exports
/******/ 	// This entry module is referenced by other modules so it can't be inlined
/******/ 	var __nested_webpack_exports__ = __nested_webpack_require_56127__(579);
/******/ 	
/******/ 	return __nested_webpack_exports__;
/******/ })()
;
});

/***/ }),

/***/ 75925:
/***/ (function(module, __unused_webpack_exports, __webpack_require__) {

(function webpackUniversalModuleDefinition(root, factory) {
	if(true)
		module.exports = factory(__webpack_require__(12591));
	else {}
})(this, function(__WEBPACK_EXTERNAL_MODULE__551__) {
return /******/ (() => { // webpackBootstrap
/******/ 	"use strict";
/******/ 	var __webpack_modules__ = ({

/***/ 45:
/***/ ((module, __unused_webpack_exports, __nested_webpack_require_630__) => {



var coseBase = {};

coseBase.layoutBase = __nested_webpack_require_630__(551);
coseBase.CoSEConstants = __nested_webpack_require_630__(806);
coseBase.CoSEEdge = __nested_webpack_require_630__(767);
coseBase.CoSEGraph = __nested_webpack_require_630__(880);
coseBase.CoSEGraphManager = __nested_webpack_require_630__(578);
coseBase.CoSELayout = __nested_webpack_require_630__(765);
coseBase.CoSENode = __nested_webpack_require_630__(991);
coseBase.ConstraintHandler = __nested_webpack_require_630__(902);

module.exports = coseBase;

/***/ }),

/***/ 806:
/***/ ((module, __unused_webpack_exports, __nested_webpack_require_1167__) => {



var FDLayoutConstants = __nested_webpack_require_1167__(551).FDLayoutConstants;

function CoSEConstants() {}

//CoSEConstants inherits static props in FDLayoutConstants
for (var prop in FDLayoutConstants) {
  CoSEConstants[prop] = FDLayoutConstants[prop];
}

CoSEConstants.DEFAULT_USE_MULTI_LEVEL_SCALING = false;
CoSEConstants.DEFAULT_RADIAL_SEPARATION = FDLayoutConstants.DEFAULT_EDGE_LENGTH;
CoSEConstants.DEFAULT_COMPONENT_SEPERATION = 60;
CoSEConstants.TILE = true;
CoSEConstants.TILING_PADDING_VERTICAL = 10;
CoSEConstants.TILING_PADDING_HORIZONTAL = 10;
CoSEConstants.TRANSFORM_ON_CONSTRAINT_HANDLING = true;
CoSEConstants.ENFORCE_CONSTRAINTS = true;
CoSEConstants.APPLY_LAYOUT = true;
CoSEConstants.RELAX_MOVEMENT_ON_CONSTRAINTS = true;
CoSEConstants.TREE_REDUCTION_ON_INCREMENTAL = true; // this should be set to false if there will be a constraint
// This constant is for differentiating whether actual layout algorithm that uses cose-base wants to apply only incremental layout or 
// an incremental layout on top of a randomized layout. If it is only incremental layout, then this constant should be true.
CoSEConstants.PURE_INCREMENTAL = CoSEConstants.DEFAULT_INCREMENTAL;

module.exports = CoSEConstants;

/***/ }),

/***/ 767:
/***/ ((module, __unused_webpack_exports, __nested_webpack_require_2468__) => {



var FDLayoutEdge = __nested_webpack_require_2468__(551).FDLayoutEdge;

function CoSEEdge(source, target, vEdge) {
  FDLayoutEdge.call(this, source, target, vEdge);
}

CoSEEdge.prototype = Object.create(FDLayoutEdge.prototype);
for (var prop in FDLayoutEdge) {
  CoSEEdge[prop] = FDLayoutEdge[prop];
}

module.exports = CoSEEdge;

/***/ }),

/***/ 880:
/***/ ((module, __unused_webpack_exports, __nested_webpack_require_2879__) => {



var LGraph = __nested_webpack_require_2879__(551).LGraph;

function CoSEGraph(parent, graphMgr, vGraph) {
  LGraph.call(this, parent, graphMgr, vGraph);
}

CoSEGraph.prototype = Object.create(LGraph.prototype);
for (var prop in LGraph) {
  CoSEGraph[prop] = LGraph[prop];
}

module.exports = CoSEGraph;

/***/ }),

/***/ 578:
/***/ ((module, __unused_webpack_exports, __nested_webpack_require_3264__) => {



var LGraphManager = __nested_webpack_require_3264__(551).LGraphManager;

function CoSEGraphManager(layout) {
  LGraphManager.call(this, layout);
}

CoSEGraphManager.prototype = Object.create(LGraphManager.prototype);
for (var prop in LGraphManager) {
  CoSEGraphManager[prop] = LGraphManager[prop];
}

module.exports = CoSEGraphManager;

/***/ }),

/***/ 765:
/***/ ((module, __unused_webpack_exports, __nested_webpack_require_3683__) => {



var FDLayout = __nested_webpack_require_3683__(551).FDLayout;
var CoSEGraphManager = __nested_webpack_require_3683__(578);
var CoSEGraph = __nested_webpack_require_3683__(880);
var CoSENode = __nested_webpack_require_3683__(991);
var CoSEEdge = __nested_webpack_require_3683__(767);
var CoSEConstants = __nested_webpack_require_3683__(806);
var ConstraintHandler = __nested_webpack_require_3683__(902);
var FDLayoutConstants = __nested_webpack_require_3683__(551).FDLayoutConstants;
var LayoutConstants = __nested_webpack_require_3683__(551).LayoutConstants;
var Point = __nested_webpack_require_3683__(551).Point;
var PointD = __nested_webpack_require_3683__(551).PointD;
var DimensionD = __nested_webpack_require_3683__(551).DimensionD;
var Layout = __nested_webpack_require_3683__(551).Layout;
var Integer = __nested_webpack_require_3683__(551).Integer;
var IGeometry = __nested_webpack_require_3683__(551).IGeometry;
var LGraph = __nested_webpack_require_3683__(551).LGraph;
var Transform = __nested_webpack_require_3683__(551).Transform;
var LinkedList = __nested_webpack_require_3683__(551).LinkedList;

function CoSELayout() {
  FDLayout.call(this);

  this.toBeTiled = {}; // Memorize if a node is to be tiled or is tiled
  this.constraints = {}; // keep layout constraints
}

CoSELayout.prototype = Object.create(FDLayout.prototype);

for (var prop in FDLayout) {
  CoSELayout[prop] = FDLayout[prop];
}

CoSELayout.prototype.newGraphManager = function () {
  var gm = new CoSEGraphManager(this);
  this.graphManager = gm;
  return gm;
};

CoSELayout.prototype.newGraph = function (vGraph) {
  return new CoSEGraph(null, this.graphManager, vGraph);
};

CoSELayout.prototype.newNode = function (vNode) {
  return new CoSENode(this.graphManager, vNode);
};

CoSELayout.prototype.newEdge = function (vEdge) {
  return new CoSEEdge(null, null, vEdge);
};

CoSELayout.prototype.initParameters = function () {
  FDLayout.prototype.initParameters.call(this, arguments);
  if (!this.isSubLayout) {
    if (CoSEConstants.DEFAULT_EDGE_LENGTH < 10) {
      this.idealEdgeLength = 10;
    } else {
      this.idealEdgeLength = CoSEConstants.DEFAULT_EDGE_LENGTH;
    }

    this.useSmartIdealEdgeLengthCalculation = CoSEConstants.DEFAULT_USE_SMART_IDEAL_EDGE_LENGTH_CALCULATION;
    this.gravityConstant = FDLayoutConstants.DEFAULT_GRAVITY_STRENGTH;
    this.compoundGravityConstant = FDLayoutConstants.DEFAULT_COMPOUND_GRAVITY_STRENGTH;
    this.gravityRangeFactor = FDLayoutConstants.DEFAULT_GRAVITY_RANGE_FACTOR;
    this.compoundGravityRangeFactor = FDLayoutConstants.DEFAULT_COMPOUND_GRAVITY_RANGE_FACTOR;

    // variables for tree reduction support
    this.prunedNodesAll = [];
    this.growTreeIterations = 0;
    this.afterGrowthIterations = 0;
    this.isTreeGrowing = false;
    this.isGrowthFinished = false;
  }
};

// This method is used to set CoSE related parameters used by spring embedder.
CoSELayout.prototype.initSpringEmbedder = function () {
  FDLayout.prototype.initSpringEmbedder.call(this);

  // variables for cooling
  this.coolingCycle = 0;
  this.maxCoolingCycle = this.maxIterations / FDLayoutConstants.CONVERGENCE_CHECK_PERIOD;
  this.finalTemperature = 0.04;
  this.coolingAdjuster = 1;
};

CoSELayout.prototype.layout = function () {
  var createBendsAsNeeded = LayoutConstants.DEFAULT_CREATE_BENDS_AS_NEEDED;
  if (createBendsAsNeeded) {
    this.createBendpoints();
    this.graphManager.resetAllEdges();
  }

  this.level = 0;
  return this.classicLayout();
};

CoSELayout.prototype.classicLayout = function () {
  this.nodesWithGravity = this.calculateNodesToApplyGravitationTo();
  this.graphManager.setAllNodesToApplyGravitation(this.nodesWithGravity);
  this.calcNoOfChildrenForAllNodes();
  this.graphManager.calcLowestCommonAncestors();
  this.graphManager.calcInclusionTreeDepths();
  this.graphManager.getRoot().calcEstimatedSize();
  this.calcIdealEdgeLengths();

  if (!this.incremental) {
    var forest = this.getFlatForest();

    // The graph associated with this layout is flat and a forest
    if (forest.length > 0) {
      this.positionNodesRadially(forest);
    }
    // The graph associated with this layout is not flat or a forest
    else {
        // Reduce the trees when incremental mode is not enabled and graph is not a forest 
        this.reduceTrees();
        // Update nodes that gravity will be applied
        this.graphManager.resetAllNodesToApplyGravitation();
        var allNodes = new Set(this.getAllNodes());
        var intersection = this.nodesWithGravity.filter(function (x) {
          return allNodes.has(x);
        });
        this.graphManager.setAllNodesToApplyGravitation(intersection);

        this.positionNodesRandomly();
      }
  } else {
    if (CoSEConstants.TREE_REDUCTION_ON_INCREMENTAL) {
      // Reduce the trees in incremental mode if only this constant is set to true 
      this.reduceTrees();
      // Update nodes that gravity will be applied
      this.graphManager.resetAllNodesToApplyGravitation();
      var allNodes = new Set(this.getAllNodes());
      var intersection = this.nodesWithGravity.filter(function (x) {
        return allNodes.has(x);
      });
      this.graphManager.setAllNodesToApplyGravitation(intersection);
    }
  }

  if (Object.keys(this.constraints).length > 0) {
    ConstraintHandler.handleConstraints(this);
    this.initConstraintVariables();
  }

  this.initSpringEmbedder();
  if (CoSEConstants.APPLY_LAYOUT) {
    this.runSpringEmbedder();
  }

  return true;
};

CoSELayout.prototype.tick = function () {
  this.totalIterations++;

  if (this.totalIterations === this.maxIterations && !this.isTreeGrowing && !this.isGrowthFinished) {
    if (this.prunedNodesAll.length > 0) {
      this.isTreeGrowing = true;
    } else {
      return true;
    }
  }

  if (this.totalIterations % FDLayoutConstants.CONVERGENCE_CHECK_PERIOD == 0 && !this.isTreeGrowing && !this.isGrowthFinished) {
    if (this.isConverged()) {
      if (this.prunedNodesAll.length > 0) {
        this.isTreeGrowing = true;
      } else {
        return true;
      }
    }

    this.coolingCycle++;

    if (this.layoutQuality == 0) {
      // quality - "draft"
      this.coolingAdjuster = this.coolingCycle;
    } else if (this.layoutQuality == 1) {
      // quality - "default"
      this.coolingAdjuster = this.coolingCycle / 3;
    }

    // cooling schedule is based on http://www.btluke.com/simanf1.html -> cooling schedule 3
    this.coolingFactor = Math.max(this.initialCoolingFactor - Math.pow(this.coolingCycle, Math.log(100 * (this.initialCoolingFactor - this.finalTemperature)) / Math.log(this.maxCoolingCycle)) / 100 * this.coolingAdjuster, this.finalTemperature);
    this.animationPeriod = Math.ceil(this.initialAnimationPeriod * Math.sqrt(this.coolingFactor));
  }
  // Operations while tree is growing again 
  if (this.isTreeGrowing) {
    if (this.growTreeIterations % 10 == 0) {
      if (this.prunedNodesAll.length > 0) {
        this.graphManager.updateBounds();
        this.updateGrid();
        this.growTree(this.prunedNodesAll);
        // Update nodes that gravity will be applied
        this.graphManager.resetAllNodesToApplyGravitation();
        var allNodes = new Set(this.getAllNodes());
        var intersection = this.nodesWithGravity.filter(function (x) {
          return allNodes.has(x);
        });
        this.graphManager.setAllNodesToApplyGravitation(intersection);

        this.graphManager.updateBounds();
        this.updateGrid();
        if (CoSEConstants.PURE_INCREMENTAL) this.coolingFactor = FDLayoutConstants.DEFAULT_COOLING_FACTOR_INCREMENTAL / 2;else this.coolingFactor = FDLayoutConstants.DEFAULT_COOLING_FACTOR_INCREMENTAL;
      } else {
        this.isTreeGrowing = false;
        this.isGrowthFinished = true;
      }
    }
    this.growTreeIterations++;
  }
  // Operations after growth is finished
  if (this.isGrowthFinished) {
    if (this.isConverged()) {
      return true;
    }
    if (this.afterGrowthIterations % 10 == 0) {
      this.graphManager.updateBounds();
      this.updateGrid();
    }
    if (CoSEConstants.PURE_INCREMENTAL) this.coolingFactor = FDLayoutConstants.DEFAULT_COOLING_FACTOR_INCREMENTAL / 2 * ((100 - this.afterGrowthIterations) / 100);else this.coolingFactor = FDLayoutConstants.DEFAULT_COOLING_FACTOR_INCREMENTAL * ((100 - this.afterGrowthIterations) / 100);
    this.afterGrowthIterations++;
  }

  var gridUpdateAllowed = !this.isTreeGrowing && !this.isGrowthFinished;
  var forceToNodeSurroundingUpdate = this.growTreeIterations % 10 == 1 && this.isTreeGrowing || this.afterGrowthIterations % 10 == 1 && this.isGrowthFinished;

  this.totalDisplacement = 0;
  this.graphManager.updateBounds();
  this.calcSpringForces();
  this.calcRepulsionForces(gridUpdateAllowed, forceToNodeSurroundingUpdate);
  this.calcGravitationalForces();
  this.moveNodes();
  this.animate();

  return false; // Layout is not ended yet return false
};

CoSELayout.prototype.getPositionsData = function () {
  var allNodes = this.graphManager.getAllNodes();
  var pData = {};
  for (var i = 0; i < allNodes.length; i++) {
    var rect = allNodes[i].rect;
    var id = allNodes[i].id;
    pData[id] = {
      id: id,
      x: rect.getCenterX(),
      y: rect.getCenterY(),
      w: rect.width,
      h: rect.height
    };
  }

  return pData;
};

CoSELayout.prototype.runSpringEmbedder = function () {
  this.initialAnimationPeriod = 25;
  this.animationPeriod = this.initialAnimationPeriod;
  var layoutEnded = false;

  // If aminate option is 'during' signal that layout is supposed to start iterating
  if (FDLayoutConstants.ANIMATE === 'during') {
    this.emit('layoutstarted');
  } else {
    // If aminate option is 'during' tick() function will be called on index.js
    while (!layoutEnded) {
      layoutEnded = this.tick();
    }

    this.graphManager.updateBounds();
  }
};

// overrides moveNodes method in FDLayout
CoSELayout.prototype.moveNodes = function () {
  var lNodes = this.getAllNodes();
  var node;

  // calculate displacement for each node 
  for (var i = 0; i < lNodes.length; i++) {
    node = lNodes[i];
    node.calculateDisplacement();
  }

  if (Object.keys(this.constraints).length > 0) {
    this.updateDisplacements();
  }

  // move each node
  for (var i = 0; i < lNodes.length; i++) {
    node = lNodes[i];
    node.move();
  }
};

// constraint related methods: initConstraintVariables and updateDisplacements

// initialize constraint related variables
CoSELayout.prototype.initConstraintVariables = function () {
  var self = this;
  this.idToNodeMap = new Map();
  this.fixedNodeSet = new Set();

  var allNodes = this.graphManager.getAllNodes();

  // fill idToNodeMap
  for (var i = 0; i < allNodes.length; i++) {
    var node = allNodes[i];
    this.idToNodeMap.set(node.id, node);
  }

  // calculate fixed node weight for given compound node
  var calculateCompoundWeight = function calculateCompoundWeight(compoundNode) {
    var nodes = compoundNode.getChild().getNodes();
    var node;
    var fixedNodeWeight = 0;
    for (var i = 0; i < nodes.length; i++) {
      node = nodes[i];
      if (node.getChild() == null) {
        if (self.fixedNodeSet.has(node.id)) {
          fixedNodeWeight += 100;
        }
      } else {
        fixedNodeWeight += calculateCompoundWeight(node);
      }
    }
    return fixedNodeWeight;
  };

  if (this.constraints.fixedNodeConstraint) {
    // fill fixedNodeSet
    this.constraints.fixedNodeConstraint.forEach(function (nodeData) {
      self.fixedNodeSet.add(nodeData.nodeId);
    });

    // assign fixed node weights to compounds if they contain fixed nodes
    var allNodes = this.graphManager.getAllNodes();
    var node;

    for (var i = 0; i < allNodes.length; i++) {
      node = allNodes[i];
      if (node.getChild() != null) {
        var fixedNodeWeight = calculateCompoundWeight(node);
        if (fixedNodeWeight > 0) {
          node.fixedNodeWeight = fixedNodeWeight;
        }
      }
    }
  }

  if (this.constraints.relativePlacementConstraint) {
    var nodeToDummyForVerticalAlignment = new Map();
    var nodeToDummyForHorizontalAlignment = new Map();
    this.dummyToNodeForVerticalAlignment = new Map();
    this.dummyToNodeForHorizontalAlignment = new Map();
    this.fixedNodesOnHorizontal = new Set();
    this.fixedNodesOnVertical = new Set();

    // fill maps and sets
    this.fixedNodeSet.forEach(function (nodeId) {
      self.fixedNodesOnHorizontal.add(nodeId);
      self.fixedNodesOnVertical.add(nodeId);
    });

    if (this.constraints.alignmentConstraint) {
      if (this.constraints.alignmentConstraint.vertical) {
        var verticalAlignment = this.constraints.alignmentConstraint.vertical;
        for (var i = 0; i < verticalAlignment.length; i++) {
          this.dummyToNodeForVerticalAlignment.set("dummy" + i, []);
          verticalAlignment[i].forEach(function (nodeId) {
            nodeToDummyForVerticalAlignment.set(nodeId, "dummy" + i);
            self.dummyToNodeForVerticalAlignment.get("dummy" + i).push(nodeId);
            if (self.fixedNodeSet.has(nodeId)) {
              self.fixedNodesOnHorizontal.add("dummy" + i);
            }
          });
        }
      }
      if (this.constraints.alignmentConstraint.horizontal) {
        var horizontalAlignment = this.constraints.alignmentConstraint.horizontal;
        for (var i = 0; i < horizontalAlignment.length; i++) {
          this.dummyToNodeForHorizontalAlignment.set("dummy" + i, []);
          horizontalAlignment[i].forEach(function (nodeId) {
            nodeToDummyForHorizontalAlignment.set(nodeId, "dummy" + i);
            self.dummyToNodeForHorizontalAlignment.get("dummy" + i).push(nodeId);
            if (self.fixedNodeSet.has(nodeId)) {
              self.fixedNodesOnVertical.add("dummy" + i);
            }
          });
        }
      }
    }

    if (CoSEConstants.RELAX_MOVEMENT_ON_CONSTRAINTS) {

      this.shuffle = function (array) {
        var j, x, i;
        for (i = array.length - 1; i >= 2 * array.length / 3; i--) {
          j = Math.floor(Math.random() * (i + 1));
          x = array[i];
          array[i] = array[j];
          array[j] = x;
        }
        return array;
      };

      this.nodesInRelativeHorizontal = [];
      this.nodesInRelativeVertical = [];
      this.nodeToRelativeConstraintMapHorizontal = new Map();
      this.nodeToRelativeConstraintMapVertical = new Map();
      this.nodeToTempPositionMapHorizontal = new Map();
      this.nodeToTempPositionMapVertical = new Map();

      // fill arrays and maps
      this.constraints.relativePlacementConstraint.forEach(function (constraint) {
        if (constraint.left) {
          var nodeIdLeft = nodeToDummyForVerticalAlignment.has(constraint.left) ? nodeToDummyForVerticalAlignment.get(constraint.left) : constraint.left;
          var nodeIdRight = nodeToDummyForVerticalAlignment.has(constraint.right) ? nodeToDummyForVerticalAlignment.get(constraint.right) : constraint.right;

          if (!self.nodesInRelativeHorizontal.includes(nodeIdLeft)) {
            self.nodesInRelativeHorizontal.push(nodeIdLeft);
            self.nodeToRelativeConstraintMapHorizontal.set(nodeIdLeft, []);
            if (self.dummyToNodeForVerticalAlignment.has(nodeIdLeft)) {
              self.nodeToTempPositionMapHorizontal.set(nodeIdLeft, self.idToNodeMap.get(self.dummyToNodeForVerticalAlignment.get(nodeIdLeft)[0]).getCenterX());
            } else {
              self.nodeToTempPositionMapHorizontal.set(nodeIdLeft, self.idToNodeMap.get(nodeIdLeft).getCenterX());
            }
          }
          if (!self.nodesInRelativeHorizontal.includes(nodeIdRight)) {
            self.nodesInRelativeHorizontal.push(nodeIdRight);
            self.nodeToRelativeConstraintMapHorizontal.set(nodeIdRight, []);
            if (self.dummyToNodeForVerticalAlignment.has(nodeIdRight)) {
              self.nodeToTempPositionMapHorizontal.set(nodeIdRight, self.idToNodeMap.get(self.dummyToNodeForVerticalAlignment.get(nodeIdRight)[0]).getCenterX());
            } else {
              self.nodeToTempPositionMapHorizontal.set(nodeIdRight, self.idToNodeMap.get(nodeIdRight).getCenterX());
            }
          }

          self.nodeToRelativeConstraintMapHorizontal.get(nodeIdLeft).push({ right: nodeIdRight, gap: constraint.gap });
          self.nodeToRelativeConstraintMapHorizontal.get(nodeIdRight).push({ left: nodeIdLeft, gap: constraint.gap });
        } else {
          var nodeIdTop = nodeToDummyForHorizontalAlignment.has(constraint.top) ? nodeToDummyForHorizontalAlignment.get(constraint.top) : constraint.top;
          var nodeIdBottom = nodeToDummyForHorizontalAlignment.has(constraint.bottom) ? nodeToDummyForHorizontalAlignment.get(constraint.bottom) : constraint.bottom;

          if (!self.nodesInRelativeVertical.includes(nodeIdTop)) {
            self.nodesInRelativeVertical.push(nodeIdTop);
            self.nodeToRelativeConstraintMapVertical.set(nodeIdTop, []);
            if (self.dummyToNodeForHorizontalAlignment.has(nodeIdTop)) {
              self.nodeToTempPositionMapVertical.set(nodeIdTop, self.idToNodeMap.get(self.dummyToNodeForHorizontalAlignment.get(nodeIdTop)[0]).getCenterY());
            } else {
              self.nodeToTempPositionMapVertical.set(nodeIdTop, self.idToNodeMap.get(nodeIdTop).getCenterY());
            }
          }
          if (!self.nodesInRelativeVertical.includes(nodeIdBottom)) {
            self.nodesInRelativeVertical.push(nodeIdBottom);
            self.nodeToRelativeConstraintMapVertical.set(nodeIdBottom, []);
            if (self.dummyToNodeForHorizontalAlignment.has(nodeIdBottom)) {
              self.nodeToTempPositionMapVertical.set(nodeIdBottom, self.idToNodeMap.get(self.dummyToNodeForHorizontalAlignment.get(nodeIdBottom)[0]).getCenterY());
            } else {
              self.nodeToTempPositionMapVertical.set(nodeIdBottom, self.idToNodeMap.get(nodeIdBottom).getCenterY());
            }
          }
          self.nodeToRelativeConstraintMapVertical.get(nodeIdTop).push({ bottom: nodeIdBottom, gap: constraint.gap });
          self.nodeToRelativeConstraintMapVertical.get(nodeIdBottom).push({ top: nodeIdTop, gap: constraint.gap });
        }
      });
    } else {
      var subGraphOnHorizontal = new Map(); // subgraph from vertical RP constraints
      var subGraphOnVertical = new Map(); // subgraph from vertical RP constraints

      // construct subgraphs from relative placement constraints 
      this.constraints.relativePlacementConstraint.forEach(function (constraint) {
        if (constraint.left) {
          var left = nodeToDummyForVerticalAlignment.has(constraint.left) ? nodeToDummyForVerticalAlignment.get(constraint.left) : constraint.left;
          var right = nodeToDummyForVerticalAlignment.has(constraint.right) ? nodeToDummyForVerticalAlignment.get(constraint.right) : constraint.right;
          if (subGraphOnHorizontal.has(left)) {
            subGraphOnHorizontal.get(left).push(right);
          } else {
            subGraphOnHorizontal.set(left, [right]);
          }
          if (subGraphOnHorizontal.has(right)) {
            subGraphOnHorizontal.get(right).push(left);
          } else {
            subGraphOnHorizontal.set(right, [left]);
          }
        } else {
          var top = nodeToDummyForHorizontalAlignment.has(constraint.top) ? nodeToDummyForHorizontalAlignment.get(constraint.top) : constraint.top;
          var bottom = nodeToDummyForHorizontalAlignment.has(constraint.bottom) ? nodeToDummyForHorizontalAlignment.get(constraint.bottom) : constraint.bottom;
          if (subGraphOnVertical.has(top)) {
            subGraphOnVertical.get(top).push(bottom);
          } else {
            subGraphOnVertical.set(top, [bottom]);
          }
          if (subGraphOnVertical.has(bottom)) {
            subGraphOnVertical.get(bottom).push(top);
          } else {
            subGraphOnVertical.set(bottom, [top]);
          }
        }
      });

      // function to construct components from a given graph 
      // also returns an array that keeps whether each component contains fixed node
      var constructComponents = function constructComponents(graph, fixedNodes) {
        var components = [];
        var isFixed = [];
        var queue = new LinkedList();
        var visited = new Set();
        var count = 0;

        graph.forEach(function (value, key) {
          if (!visited.has(key)) {
            components[count] = [];
            isFixed[count] = false;
            var currentNode = key;
            queue.push(currentNode);
            visited.add(currentNode);
            components[count].push(currentNode);

            while (queue.length != 0) {
              currentNode = queue.shift();
              if (fixedNodes.has(currentNode)) {
                isFixed[count] = true;
              }
              var neighbors = graph.get(currentNode);
              neighbors.forEach(function (neighbor) {
                if (!visited.has(neighbor)) {
                  queue.push(neighbor);
                  visited.add(neighbor);
                  components[count].push(neighbor);
                }
              });
            }
            count++;
          }
        });

        return { components: components, isFixed: isFixed };
      };

      var resultOnHorizontal = constructComponents(subGraphOnHorizontal, self.fixedNodesOnHorizontal);
      this.componentsOnHorizontal = resultOnHorizontal.components;
      this.fixedComponentsOnHorizontal = resultOnHorizontal.isFixed;
      var resultOnVertical = constructComponents(subGraphOnVertical, self.fixedNodesOnVertical);
      this.componentsOnVertical = resultOnVertical.components;
      this.fixedComponentsOnVertical = resultOnVertical.isFixed;
    }
  }
};

// updates node displacements based on constraints
CoSELayout.prototype.updateDisplacements = function () {
  var self = this;
  if (this.constraints.fixedNodeConstraint) {
    this.constraints.fixedNodeConstraint.forEach(function (nodeData) {
      var fixedNode = self.idToNodeMap.get(nodeData.nodeId);
      fixedNode.displacementX = 0;
      fixedNode.displacementY = 0;
    });
  }

  if (this.constraints.alignmentConstraint) {
    if (this.constraints.alignmentConstraint.vertical) {
      var allVerticalAlignments = this.constraints.alignmentConstraint.vertical;
      for (var i = 0; i < allVerticalAlignments.length; i++) {
        var totalDisplacementX = 0;
        for (var j = 0; j < allVerticalAlignments[i].length; j++) {
          if (this.fixedNodeSet.has(allVerticalAlignments[i][j])) {
            totalDisplacementX = 0;
            break;
          }
          totalDisplacementX += this.idToNodeMap.get(allVerticalAlignments[i][j]).displacementX;
        }
        var averageDisplacementX = totalDisplacementX / allVerticalAlignments[i].length;
        for (var j = 0; j < allVerticalAlignments[i].length; j++) {
          this.idToNodeMap.get(allVerticalAlignments[i][j]).displacementX = averageDisplacementX;
        }
      }
    }
    if (this.constraints.alignmentConstraint.horizontal) {
      var allHorizontalAlignments = this.constraints.alignmentConstraint.horizontal;
      for (var i = 0; i < allHorizontalAlignments.length; i++) {
        var totalDisplacementY = 0;
        for (var j = 0; j < allHorizontalAlignments[i].length; j++) {
          if (this.fixedNodeSet.has(allHorizontalAlignments[i][j])) {
            totalDisplacementY = 0;
            break;
          }
          totalDisplacementY += this.idToNodeMap.get(allHorizontalAlignments[i][j]).displacementY;
        }
        var averageDisplacementY = totalDisplacementY / allHorizontalAlignments[i].length;
        for (var j = 0; j < allHorizontalAlignments[i].length; j++) {
          this.idToNodeMap.get(allHorizontalAlignments[i][j]).displacementY = averageDisplacementY;
        }
      }
    }
  }

  if (this.constraints.relativePlacementConstraint) {

    if (CoSEConstants.RELAX_MOVEMENT_ON_CONSTRAINTS) {
      // shuffle array to randomize node processing order
      if (this.totalIterations % 10 == 0) {
        this.shuffle(this.nodesInRelativeHorizontal);
        this.shuffle(this.nodesInRelativeVertical);
      }

      this.nodesInRelativeHorizontal.forEach(function (nodeId) {
        if (!self.fixedNodesOnHorizontal.has(nodeId)) {
          var displacement = 0;
          if (self.dummyToNodeForVerticalAlignment.has(nodeId)) {
            displacement = self.idToNodeMap.get(self.dummyToNodeForVerticalAlignment.get(nodeId)[0]).displacementX;
          } else {
            displacement = self.idToNodeMap.get(nodeId).displacementX;
          }
          self.nodeToRelativeConstraintMapHorizontal.get(nodeId).forEach(function (constraint) {
            if (constraint.right) {
              var diff = self.nodeToTempPositionMapHorizontal.get(constraint.right) - self.nodeToTempPositionMapHorizontal.get(nodeId) - displacement;
              if (diff < constraint.gap) {
                displacement -= constraint.gap - diff;
              }
            } else {
              var diff = self.nodeToTempPositionMapHorizontal.get(nodeId) - self.nodeToTempPositionMapHorizontal.get(constraint.left) + displacement;
              if (diff < constraint.gap) {
                displacement += constraint.gap - diff;
              }
            }
          });
          self.nodeToTempPositionMapHorizontal.set(nodeId, self.nodeToTempPositionMapHorizontal.get(nodeId) + displacement);
          if (self.dummyToNodeForVerticalAlignment.has(nodeId)) {
            self.dummyToNodeForVerticalAlignment.get(nodeId).forEach(function (nodeId) {
              self.idToNodeMap.get(nodeId).displacementX = displacement;
            });
          } else {
            self.idToNodeMap.get(nodeId).displacementX = displacement;
          }
        }
      });

      this.nodesInRelativeVertical.forEach(function (nodeId) {
        if (!self.fixedNodesOnHorizontal.has(nodeId)) {
          var displacement = 0;
          if (self.dummyToNodeForHorizontalAlignment.has(nodeId)) {
            displacement = self.idToNodeMap.get(self.dummyToNodeForHorizontalAlignment.get(nodeId)[0]).displacementY;
          } else {
            displacement = self.idToNodeMap.get(nodeId).displacementY;
          }
          self.nodeToRelativeConstraintMapVertical.get(nodeId).forEach(function (constraint) {
            if (constraint.bottom) {
              var diff = self.nodeToTempPositionMapVertical.get(constraint.bottom) - self.nodeToTempPositionMapVertical.get(nodeId) - displacement;
              if (diff < constraint.gap) {
                displacement -= constraint.gap - diff;
              }
            } else {
              var diff = self.nodeToTempPositionMapVertical.get(nodeId) - self.nodeToTempPositionMapVertical.get(constraint.top) + displacement;
              if (diff < constraint.gap) {
                displacement += constraint.gap - diff;
              }
            }
          });
          self.nodeToTempPositionMapVertical.set(nodeId, self.nodeToTempPositionMapVertical.get(nodeId) + displacement);
          if (self.dummyToNodeForHorizontalAlignment.has(nodeId)) {
            self.dummyToNodeForHorizontalAlignment.get(nodeId).forEach(function (nodeId) {
              self.idToNodeMap.get(nodeId).displacementY = displacement;
            });
          } else {
            self.idToNodeMap.get(nodeId).displacementY = displacement;
          }
        }
      });
    } else {
      for (var i = 0; i < this.componentsOnHorizontal.length; i++) {
        var component = this.componentsOnHorizontal[i];
        if (this.fixedComponentsOnHorizontal[i]) {
          for (var j = 0; j < component.length; j++) {
            if (this.dummyToNodeForVerticalAlignment.has(component[j])) {
              this.dummyToNodeForVerticalAlignment.get(component[j]).forEach(function (nodeId) {
                self.idToNodeMap.get(nodeId).displacementX = 0;
              });
            } else {
              this.idToNodeMap.get(component[j]).displacementX = 0;
            }
          }
        } else {
          var sum = 0;
          var count = 0;
          for (var j = 0; j < component.length; j++) {
            if (this.dummyToNodeForVerticalAlignment.has(component[j])) {
              var actualNodes = this.dummyToNodeForVerticalAlignment.get(component[j]);
              sum += actualNodes.length * this.idToNodeMap.get(actualNodes[0]).displacementX;
              count += actualNodes.length;
            } else {
              sum += this.idToNodeMap.get(component[j]).displacementX;
              count++;
            }
          }
          var averageDisplacement = sum / count;
          for (var j = 0; j < component.length; j++) {
            if (this.dummyToNodeForVerticalAlignment.has(component[j])) {
              this.dummyToNodeForVerticalAlignment.get(component[j]).forEach(function (nodeId) {
                self.idToNodeMap.get(nodeId).displacementX = averageDisplacement;
              });
            } else {
              this.idToNodeMap.get(component[j]).displacementX = averageDisplacement;
            }
          }
        }
      }

      for (var i = 0; i < this.componentsOnVertical.length; i++) {
        var component = this.componentsOnVertical[i];
        if (this.fixedComponentsOnVertical[i]) {
          for (var j = 0; j < component.length; j++) {
            if (this.dummyToNodeForHorizontalAlignment.has(component[j])) {
              this.dummyToNodeForHorizontalAlignment.get(component[j]).forEach(function (nodeId) {
                self.idToNodeMap.get(nodeId).displacementY = 0;
              });
            } else {
              this.idToNodeMap.get(component[j]).displacementY = 0;
            }
          }
        } else {
          var sum = 0;
          var count = 0;
          for (var j = 0; j < component.length; j++) {
            if (this.dummyToNodeForHorizontalAlignment.has(component[j])) {
              var actualNodes = this.dummyToNodeForHorizontalAlignment.get(component[j]);
              sum += actualNodes.length * this.idToNodeMap.get(actualNodes[0]).displacementY;
              count += actualNodes.length;
            } else {
              sum += this.idToNodeMap.get(component[j]).displacementY;
              count++;
            }
          }
          var averageDisplacement = sum / count;
          for (var j = 0; j < component.length; j++) {
            if (this.dummyToNodeForHorizontalAlignment.has(component[j])) {
              this.dummyToNodeForHorizontalAlignment.get(component[j]).forEach(function (nodeId) {
                self.idToNodeMap.get(nodeId).displacementY = averageDisplacement;
              });
            } else {
              this.idToNodeMap.get(component[j]).displacementY = averageDisplacement;
            }
          }
        }
      }
    }
  }
};

CoSELayout.prototype.calculateNodesToApplyGravitationTo = function () {
  var nodeList = [];
  var graph;

  var graphs = this.graphManager.getGraphs();
  var size = graphs.length;
  var i;
  for (i = 0; i < size; i++) {
    graph = graphs[i];

    graph.updateConnected();

    if (!graph.isConnected) {
      nodeList = nodeList.concat(graph.getNodes());
    }
  }

  return nodeList;
};

CoSELayout.prototype.createBendpoints = function () {
  var edges = [];
  edges = edges.concat(this.graphManager.getAllEdges());
  var visited = new Set();
  var i;
  for (i = 0; i < edges.length; i++) {
    var edge = edges[i];

    if (!visited.has(edge)) {
      var source = edge.getSource();
      var target = edge.getTarget();

      if (source == target) {
        edge.getBendpoints().push(new PointD());
        edge.getBendpoints().push(new PointD());
        this.createDummyNodesForBendpoints(edge);
        visited.add(edge);
      } else {
        var edgeList = [];

        edgeList = edgeList.concat(source.getEdgeListToNode(target));
        edgeList = edgeList.concat(target.getEdgeListToNode(source));

        if (!visited.has(edgeList[0])) {
          if (edgeList.length > 1) {
            var k;
            for (k = 0; k < edgeList.length; k++) {
              var multiEdge = edgeList[k];
              multiEdge.getBendpoints().push(new PointD());
              this.createDummyNodesForBendpoints(multiEdge);
            }
          }
          edgeList.forEach(function (edge) {
            visited.add(edge);
          });
        }
      }
    }

    if (visited.size == edges.length) {
      break;
    }
  }
};

CoSELayout.prototype.positionNodesRadially = function (forest) {
  // We tile the trees to a grid row by row; first tree starts at (0,0)
  var currentStartingPoint = new Point(0, 0);
  var numberOfColumns = Math.ceil(Math.sqrt(forest.length));
  var height = 0;
  var currentY = 0;
  var currentX = 0;
  var point = new PointD(0, 0);

  for (var i = 0; i < forest.length; i++) {
    if (i % numberOfColumns == 0) {
      // Start of a new row, make the x coordinate 0, increment the
      // y coordinate with the max height of the previous row
      currentX = 0;
      currentY = height;

      if (i != 0) {
        currentY += CoSEConstants.DEFAULT_COMPONENT_SEPERATION;
      }

      height = 0;
    }

    var tree = forest[i];

    // Find the center of the tree
    var centerNode = Layout.findCenterOfTree(tree);

    // Set the staring point of the next tree
    currentStartingPoint.x = currentX;
    currentStartingPoint.y = currentY;

    // Do a radial layout starting with the center
    point = CoSELayout.radialLayout(tree, centerNode, currentStartingPoint);

    if (point.y > height) {
      height = Math.floor(point.y);
    }

    currentX = Math.floor(point.x + CoSEConstants.DEFAULT_COMPONENT_SEPERATION);
  }

  this.transform(new PointD(LayoutConstants.WORLD_CENTER_X - point.x / 2, LayoutConstants.WORLD_CENTER_Y - point.y / 2));
};

CoSELayout.radialLayout = function (tree, centerNode, startingPoint) {
  var radialSep = Math.max(this.maxDiagonalInTree(tree), CoSEConstants.DEFAULT_RADIAL_SEPARATION);
  CoSELayout.branchRadialLayout(centerNode, null, 0, 359, 0, radialSep);
  var bounds = LGraph.calculateBounds(tree);

  var transform = new Transform();
  transform.setDeviceOrgX(bounds.getMinX());
  transform.setDeviceOrgY(bounds.getMinY());
  transform.setWorldOrgX(startingPoint.x);
  transform.setWorldOrgY(startingPoint.y);

  for (var i = 0; i < tree.length; i++) {
    var node = tree[i];
    node.transform(transform);
  }

  var bottomRight = new PointD(bounds.getMaxX(), bounds.getMaxY());

  return transform.inverseTransformPoint(bottomRight);
};

CoSELayout.branchRadialLayout = function (node, parentOfNode, startAngle, endAngle, distance, radialSeparation) {
  // First, position this node by finding its angle.
  var halfInterval = (endAngle - startAngle + 1) / 2;

  if (halfInterval < 0) {
    halfInterval += 180;
  }

  var nodeAngle = (halfInterval + startAngle) % 360;
  var teta = nodeAngle * IGeometry.TWO_PI / 360;

  // Make polar to java cordinate conversion.
  var cos_teta = Math.cos(teta);
  var x_ = distance * Math.cos(teta);
  var y_ = distance * Math.sin(teta);

  node.setCenter(x_, y_);

  // Traverse all neighbors of this node and recursively call this
  // function.
  var neighborEdges = [];
  neighborEdges = neighborEdges.concat(node.getEdges());
  var childCount = neighborEdges.length;

  if (parentOfNode != null) {
    childCount--;
  }

  var branchCount = 0;

  var incEdgesCount = neighborEdges.length;
  var startIndex;

  var edges = node.getEdgesBetween(parentOfNode);

  // If there are multiple edges, prune them until there remains only one
  // edge.
  while (edges.length > 1) {
    //neighborEdges.remove(edges.remove(0));
    var temp = edges[0];
    edges.splice(0, 1);
    var index = neighborEdges.indexOf(temp);
    if (index >= 0) {
      neighborEdges.splice(index, 1);
    }
    incEdgesCount--;
    childCount--;
  }

  if (parentOfNode != null) {
    //assert edges.length == 1;
    startIndex = (neighborEdges.indexOf(edges[0]) + 1) % incEdgesCount;
  } else {
    startIndex = 0;
  }

  var stepAngle = Math.abs(endAngle - startAngle) / childCount;

  for (var i = startIndex; branchCount != childCount; i = ++i % incEdgesCount) {
    var currentNeighbor = neighborEdges[i].getOtherEnd(node);

    // Don't back traverse to root node in current tree.
    if (currentNeighbor == parentOfNode) {
      continue;
    }

    var childStartAngle = (startAngle + branchCount * stepAngle) % 360;
    var childEndAngle = (childStartAngle + stepAngle) % 360;

    CoSELayout.branchRadialLayout(currentNeighbor, node, childStartAngle, childEndAngle, distance + radialSeparation, radialSeparation);

    branchCount++;
  }
};

CoSELayout.maxDiagonalInTree = function (tree) {
  var maxDiagonal = Integer.MIN_VALUE;

  for (var i = 0; i < tree.length; i++) {
    var node = tree[i];
    var diagonal = node.getDiagonal();

    if (diagonal > maxDiagonal) {
      maxDiagonal = diagonal;
    }
  }

  return maxDiagonal;
};

CoSELayout.prototype.calcRepulsionRange = function () {
  // formula is 2 x (level + 1) x idealEdgeLength
  return 2 * (this.level + 1) * this.idealEdgeLength;
};

// Tiling methods

// Group zero degree members whose parents are not to be tiled, create dummy parents where needed and fill memberGroups by their dummp parent id's
CoSELayout.prototype.groupZeroDegreeMembers = function () {
  var self = this;
  // array of [parent_id x oneDegreeNode_id]
  var tempMemberGroups = {}; // A temporary map of parent node and its zero degree members
  this.memberGroups = {}; // A map of dummy parent node and its zero degree members whose parents are not to be tiled
  this.idToDummyNode = {}; // A map of id to dummy node 

  var zeroDegree = []; // List of zero degree nodes whose parents are not to be tiled
  var allNodes = this.graphManager.getAllNodes();

  // Fill zero degree list
  for (var i = 0; i < allNodes.length; i++) {
    var node = allNodes[i];
    var parent = node.getParent();
    // If a node has zero degree and its parent is not to be tiled if exists add that node to zeroDegres list
    if (this.getNodeDegreeWithChildren(node) === 0 && (parent.id == undefined || !this.getToBeTiled(parent))) {
      zeroDegree.push(node);
    }
  }

  // Create a map of parent node and its zero degree members
  for (var i = 0; i < zeroDegree.length; i++) {
    var node = zeroDegree[i]; // Zero degree node itself
    var p_id = node.getParent().id; // Parent id

    if (typeof tempMemberGroups[p_id] === "undefined") tempMemberGroups[p_id] = [];

    tempMemberGroups[p_id] = tempMemberGroups[p_id].concat(node); // Push node to the list belongs to its parent in tempMemberGroups
  }

  // If there are at least two nodes at a level, create a dummy compound for them
  Object.keys(tempMemberGroups).forEach(function (p_id) {
    if (tempMemberGroups[p_id].length > 1) {
      var dummyCompoundId = "DummyCompound_" + p_id; // The id of dummy compound which will be created soon
      self.memberGroups[dummyCompoundId] = tempMemberGroups[p_id]; // Add dummy compound to memberGroups

      var parent = tempMemberGroups[p_id][0].getParent(); // The parent of zero degree nodes will be the parent of new dummy compound

      // Create a dummy compound with calculated id
      var dummyCompound = new CoSENode(self.graphManager);
      dummyCompound.id = dummyCompoundId;
      dummyCompound.paddingLeft = parent.paddingLeft || 0;
      dummyCompound.paddingRight = parent.paddingRight || 0;
      dummyCompound.paddingBottom = parent.paddingBottom || 0;
      dummyCompound.paddingTop = parent.paddingTop || 0;

      self.idToDummyNode[dummyCompoundId] = dummyCompound;

      var dummyParentGraph = self.getGraphManager().add(self.newGraph(), dummyCompound);
      var parentGraph = parent.getChild();

      // Add dummy compound to parent the graph
      parentGraph.add(dummyCompound);

      // For each zero degree node in this level remove it from its parent graph and add it to the graph of dummy parent
      for (var i = 0; i < tempMemberGroups[p_id].length; i++) {
        var node = tempMemberGroups[p_id][i];

        parentGraph.remove(node);
        dummyParentGraph.add(node);
      }
    }
  });
};

CoSELayout.prototype.clearCompounds = function () {
  var childGraphMap = {};
  var idToNode = {};

  // Get compound ordering by finding the inner one first
  this.performDFSOnCompounds();

  for (var i = 0; i < this.compoundOrder.length; i++) {

    idToNode[this.compoundOrder[i].id] = this.compoundOrder[i];
    childGraphMap[this.compoundOrder[i].id] = [].concat(this.compoundOrder[i].getChild().getNodes());

    // Remove children of compounds
    this.graphManager.remove(this.compoundOrder[i].getChild());
    this.compoundOrder[i].child = null;
  }

  this.graphManager.resetAllNodes();

  // Tile the removed children
  this.tileCompoundMembers(childGraphMap, idToNode);
};

CoSELayout.prototype.clearZeroDegreeMembers = function () {
  var self = this;
  var tiledZeroDegreePack = this.tiledZeroDegreePack = [];

  Object.keys(this.memberGroups).forEach(function (id) {
    var compoundNode = self.idToDummyNode[id]; // Get the dummy compound

    tiledZeroDegreePack[id] = self.tileNodes(self.memberGroups[id], compoundNode.paddingLeft + compoundNode.paddingRight);

    // Set the width and height of the dummy compound as calculated
    compoundNode.rect.width = tiledZeroDegreePack[id].width;
    compoundNode.rect.height = tiledZeroDegreePack[id].height;
    compoundNode.setCenter(tiledZeroDegreePack[id].centerX, tiledZeroDegreePack[id].centerY);

    // compound left and top margings for labels
    // when node labels are included, these values may be set to different values below and are used in tilingPostLayout,
    // otherwise they stay as zero
    compoundNode.labelMarginLeft = 0;
    compoundNode.labelMarginTop = 0;

    // Update compound bounds considering its label properties and set label margins for left and top
    if (CoSEConstants.NODE_DIMENSIONS_INCLUDE_LABELS) {

      var width = compoundNode.rect.width;
      var height = compoundNode.rect.height;

      if (compoundNode.labelWidth) {
        if (compoundNode.labelPosHorizontal == "left") {
          compoundNode.rect.x -= compoundNode.labelWidth;
          compoundNode.setWidth(width + compoundNode.labelWidth);
          compoundNode.labelMarginLeft = compoundNode.labelWidth;
        } else if (compoundNode.labelPosHorizontal == "center" && compoundNode.labelWidth > width) {
          compoundNode.rect.x -= (compoundNode.labelWidth - width) / 2;
          compoundNode.setWidth(compoundNode.labelWidth);
          compoundNode.labelMarginLeft = (compoundNode.labelWidth - width) / 2;
        } else if (compoundNode.labelPosHorizontal == "right") {
          compoundNode.setWidth(width + compoundNode.labelWidth);
        }
      }

      if (compoundNode.labelHeight) {
        if (compoundNode.labelPosVertical == "top") {
          compoundNode.rect.y -= compoundNode.labelHeight;
          compoundNode.setHeight(height + compoundNode.labelHeight);
          compoundNode.labelMarginTop = compoundNode.labelHeight;
        } else if (compoundNode.labelPosVertical == "center" && compoundNode.labelHeight > height) {
          compoundNode.rect.y -= (compoundNode.labelHeight - height) / 2;
          compoundNode.setHeight(compoundNode.labelHeight);
          compoundNode.labelMarginTop = (compoundNode.labelHeight - height) / 2;
        } else if (compoundNode.labelPosVertical == "bottom") {
          compoundNode.setHeight(height + compoundNode.labelHeight);
        }
      }
    }
  });
};

CoSELayout.prototype.repopulateCompounds = function () {
  for (var i = this.compoundOrder.length - 1; i >= 0; i--) {
    var lCompoundNode = this.compoundOrder[i];
    var id = lCompoundNode.id;
    var horizontalMargin = lCompoundNode.paddingLeft;
    var verticalMargin = lCompoundNode.paddingTop;
    var labelMarginLeft = lCompoundNode.labelMarginLeft;
    var labelMarginTop = lCompoundNode.labelMarginTop;

    this.adjustLocations(this.tiledMemberPack[id], lCompoundNode.rect.x, lCompoundNode.rect.y, horizontalMargin, verticalMargin, labelMarginLeft, labelMarginTop);
  }
};

CoSELayout.prototype.repopulateZeroDegreeMembers = function () {
  var self = this;
  var tiledPack = this.tiledZeroDegreePack;

  Object.keys(tiledPack).forEach(function (id) {
    var compoundNode = self.idToDummyNode[id]; // Get the dummy compound by its id
    var horizontalMargin = compoundNode.paddingLeft;
    var verticalMargin = compoundNode.paddingTop;
    var labelMarginLeft = compoundNode.labelMarginLeft;
    var labelMarginTop = compoundNode.labelMarginTop;

    // Adjust the positions of nodes wrt its compound
    self.adjustLocations(tiledPack[id], compoundNode.rect.x, compoundNode.rect.y, horizontalMargin, verticalMargin, labelMarginLeft, labelMarginTop);
  });
};

CoSELayout.prototype.getToBeTiled = function (node) {
  var id = node.id;
  //firstly check the previous results
  if (this.toBeTiled[id] != null) {
    return this.toBeTiled[id];
  }

  //only compound nodes are to be tiled
  var childGraph = node.getChild();
  if (childGraph == null) {
    this.toBeTiled[id] = false;
    return false;
  }

  var children = childGraph.getNodes(); // Get the children nodes

  //a compound node is not to be tiled if all of its compound children are not to be tiled
  for (var i = 0; i < children.length; i++) {
    var theChild = children[i];

    if (this.getNodeDegree(theChild) > 0) {
      this.toBeTiled[id] = false;
      return false;
    }

    //pass the children not having the compound structure
    if (theChild.getChild() == null) {
      this.toBeTiled[theChild.id] = false;
      continue;
    }

    if (!this.getToBeTiled(theChild)) {
      this.toBeTiled[id] = false;
      return false;
    }
  }
  this.toBeTiled[id] = true;
  return true;
};

// Get degree of a node depending of its edges and independent of its children
CoSELayout.prototype.getNodeDegree = function (node) {
  var id = node.id;
  var edges = node.getEdges();
  var degree = 0;

  // For the edges connected
  for (var i = 0; i < edges.length; i++) {
    var edge = edges[i];
    if (edge.getSource().id !== edge.getTarget().id) {
      degree = degree + 1;
    }
  }
  return degree;
};

// Get degree of a node with its children
CoSELayout.prototype.getNodeDegreeWithChildren = function (node) {
  var degree = this.getNodeDegree(node);
  if (node.getChild() == null) {
    return degree;
  }
  var children = node.getChild().getNodes();
  for (var i = 0; i < children.length; i++) {
    var child = children[i];
    degree += this.getNodeDegreeWithChildren(child);
  }
  return degree;
};

CoSELayout.prototype.performDFSOnCompounds = function () {
  this.compoundOrder = [];
  this.fillCompexOrderByDFS(this.graphManager.getRoot().getNodes());
};

CoSELayout.prototype.fillCompexOrderByDFS = function (children) {
  for (var i = 0; i < children.length; i++) {
    var child = children[i];
    if (child.getChild() != null) {
      this.fillCompexOrderByDFS(child.getChild().getNodes());
    }
    if (this.getToBeTiled(child)) {
      this.compoundOrder.push(child);
    }
  }
};

/**
* This method places each zero degree member wrt given (x,y) coordinates (top left).
*/
CoSELayout.prototype.adjustLocations = function (organization, x, y, compoundHorizontalMargin, compoundVerticalMargin, compoundLabelMarginLeft, compoundLabelMarginTop) {
  x += compoundHorizontalMargin + compoundLabelMarginLeft;
  y += compoundVerticalMargin + compoundLabelMarginTop;

  var left = x;

  for (var i = 0; i < organization.rows.length; i++) {
    var row = organization.rows[i];
    x = left;
    var maxHeight = 0;

    for (var j = 0; j < row.length; j++) {
      var lnode = row[j];

      lnode.rect.x = x; // + lnode.rect.width / 2;
      lnode.rect.y = y; // + lnode.rect.height / 2;

      x += lnode.rect.width + organization.horizontalPadding;

      if (lnode.rect.height > maxHeight) maxHeight = lnode.rect.height;
    }

    y += maxHeight + organization.verticalPadding;
  }
};

CoSELayout.prototype.tileCompoundMembers = function (childGraphMap, idToNode) {
  var self = this;
  this.tiledMemberPack = [];

  Object.keys(childGraphMap).forEach(function (id) {
    // Get the compound node
    var compoundNode = idToNode[id];

    self.tiledMemberPack[id] = self.tileNodes(childGraphMap[id], compoundNode.paddingLeft + compoundNode.paddingRight);

    compoundNode.rect.width = self.tiledMemberPack[id].width;
    compoundNode.rect.height = self.tiledMemberPack[id].height;
    compoundNode.setCenter(self.tiledMemberPack[id].centerX, self.tiledMemberPack[id].centerY);

    // compound left and top margings for labels
    // when node labels are included, these values may be set to different values below and are used in tilingPostLayout,
    // otherwise they stay as zero
    compoundNode.labelMarginLeft = 0;
    compoundNode.labelMarginTop = 0;

    // Update compound bounds considering its label properties and set label margins for left and top
    if (CoSEConstants.NODE_DIMENSIONS_INCLUDE_LABELS) {

      var width = compoundNode.rect.width;
      var height = compoundNode.rect.height;

      if (compoundNode.labelWidth) {
        if (compoundNode.labelPosHorizontal == "left") {
          compoundNode.rect.x -= compoundNode.labelWidth;
          compoundNode.setWidth(width + compoundNode.labelWidth);
          compoundNode.labelMarginLeft = compoundNode.labelWidth;
        } else if (compoundNode.labelPosHorizontal == "center" && compoundNode.labelWidth > width) {
          compoundNode.rect.x -= (compoundNode.labelWidth - width) / 2;
          compoundNode.setWidth(compoundNode.labelWidth);
          compoundNode.labelMarginLeft = (compoundNode.labelWidth - width) / 2;
        } else if (compoundNode.labelPosHorizontal == "right") {
          compoundNode.setWidth(width + compoundNode.labelWidth);
        }
      }

      if (compoundNode.labelHeight) {
        if (compoundNode.labelPosVertical == "top") {
          compoundNode.rect.y -= compoundNode.labelHeight;
          compoundNode.setHeight(height + compoundNode.labelHeight);
          compoundNode.labelMarginTop = compoundNode.labelHeight;
        } else if (compoundNode.labelPosVertical == "center" && compoundNode.labelHeight > height) {
          compoundNode.rect.y -= (compoundNode.labelHeight - height) / 2;
          compoundNode.setHeight(compoundNode.labelHeight);
          compoundNode.labelMarginTop = (compoundNode.labelHeight - height) / 2;
        } else if (compoundNode.labelPosVertical == "bottom") {
          compoundNode.setHeight(height + compoundNode.labelHeight);
        }
      }
    }
  });
};

CoSELayout.prototype.tileNodes = function (nodes, minWidth) {
  var horizontalOrg = this.tileNodesByFavoringDim(nodes, minWidth, true);
  var verticalOrg = this.tileNodesByFavoringDim(nodes, minWidth, false);

  var horizontalRatio = this.getOrgRatio(horizontalOrg);
  var verticalRatio = this.getOrgRatio(verticalOrg);
  var bestOrg;

  // the best ratio is the one that is closer to 1 since the ratios are already normalized
  // and the best organization is the one that has the best ratio
  if (verticalRatio < horizontalRatio) {
    bestOrg = verticalOrg;
  } else {
    bestOrg = horizontalOrg;
  }

  return bestOrg;
};

// get the width/height ratio of the organization that is normalized so that it will not be less than 1
CoSELayout.prototype.getOrgRatio = function (organization) {
  // get dimensions and calculate the initial ratio
  var width = organization.width;
  var height = organization.height;
  var ratio = width / height;

  // if the initial ratio is less then 1 then inverse it
  if (ratio < 1) {
    ratio = 1 / ratio;
  }

  // return the normalized ratio
  return ratio;
};

/*
 * Calculates the ideal width for the rows. This method assumes that
 * each node has the same sizes and calculates the ideal row width that
 * approximates a square shaped complex accordingly. However, since nodes would
 * have different sizes some rows would have different sizes and the resulting
 * shape would not be an exact square.
 */
CoSELayout.prototype.calcIdealRowWidth = function (members, favorHorizontalDim) {
  // To approximate a square shaped complex we need to make complex width equal to complex height.
  // To achieve this we need to solve the following equation system for hc:
  // (x + bx) * hc - bx = (y + by) * vc - by, hc * vc = n
  // where x is the avarage width of the nodes, y is the avarage height of nodes
  // bx and by are the buffer sizes in horizontal and vertical dimensions accordingly,
  // hc and vc are the number of rows in horizontal and vertical dimensions
  // n is number of members.

  var verticalPadding = CoSEConstants.TILING_PADDING_VERTICAL;
  var horizontalPadding = CoSEConstants.TILING_PADDING_HORIZONTAL;

  // number of members
  var membersSize = members.length;

  // sum of the width of all members
  var totalWidth = 0;

  // sum of the height of all members
  var totalHeight = 0;

  var maxWidth = 0;

  // traverse all members to calculate total width and total height and get the maximum members width
  members.forEach(function (node) {
    totalWidth += node.getWidth();
    totalHeight += node.getHeight();

    if (node.getWidth() > maxWidth) {
      maxWidth = node.getWidth();
    }
  });

  // average width of the members
  var averageWidth = totalWidth / membersSize;

  // average height of the members
  var averageHeight = totalHeight / membersSize;

  // solving the initial equation system for the hc yields the following second degree equation:
  // hc^2 * (x+bx) + hc * (by - bx) - n * (y + by) = 0

  // the delta value to solve the equation above for hc
  var delta = Math.pow(verticalPadding - horizontalPadding, 2) + 4 * (averageWidth + horizontalPadding) * (averageHeight + verticalPadding) * membersSize;

  // solve the equation using delta value to calculate the horizontal count
  // that represents the number of nodes in an ideal row
  var horizontalCountDouble = (horizontalPadding - verticalPadding + Math.sqrt(delta)) / (2 * (averageWidth + horizontalPadding));
  // round the calculated horizontal count up or down according to the favored dimension
  var horizontalCount;

  if (favorHorizontalDim) {
    horizontalCount = Math.ceil(horizontalCountDouble);
    // if horizontalCount count is not a float value then both of rounding to floor and ceil
    // will yield the same values. Instead of repeating the same calculation try going up
    // while favoring horizontal dimension in such cases
    if (horizontalCount == horizontalCountDouble) {
      horizontalCount++;
    }
  } else {
    horizontalCount = Math.floor(horizontalCountDouble);
  }

  // ideal width to be calculated
  var idealWidth = horizontalCount * (averageWidth + horizontalPadding) - horizontalPadding;

  // if max width is bigger than calculated ideal width reset ideal width to it
  if (maxWidth > idealWidth) {
    idealWidth = maxWidth;
  }

  // add the left-right margins to the ideal row width
  idealWidth += horizontalPadding * 2;

  // return the ideal row width1
  return idealWidth;
};

CoSELayout.prototype.tileNodesByFavoringDim = function (nodes, minWidth, favorHorizontalDim) {
  var verticalPadding = CoSEConstants.TILING_PADDING_VERTICAL;
  var horizontalPadding = CoSEConstants.TILING_PADDING_HORIZONTAL;
  var tilingCompareBy = CoSEConstants.TILING_COMPARE_BY;
  var organization = {
    rows: [],
    rowWidth: [],
    rowHeight: [],
    width: 0,
    height: minWidth, // assume minHeight equals to minWidth
    verticalPadding: verticalPadding,
    horizontalPadding: horizontalPadding,
    centerX: 0,
    centerY: 0
  };

  if (tilingCompareBy) {
    organization.idealRowWidth = this.calcIdealRowWidth(nodes, favorHorizontalDim);
  }

  var getNodeArea = function getNodeArea(n) {
    return n.rect.width * n.rect.height;
  };

  var areaCompareFcn = function areaCompareFcn(n1, n2) {
    return getNodeArea(n2) - getNodeArea(n1);
  };

  // Sort the nodes in descending order of their areas
  nodes.sort(function (n1, n2) {
    var cmpBy = areaCompareFcn;
    if (organization.idealRowWidth) {
      cmpBy = tilingCompareBy;
      return cmpBy(n1.id, n2.id);
    }
    return cmpBy(n1, n2);
  });

  // Create the organization -> calculate compound center
  var sumCenterX = 0;
  var sumCenterY = 0;
  for (var i = 0; i < nodes.length; i++) {
    var lNode = nodes[i];

    sumCenterX += lNode.getCenterX();
    sumCenterY += lNode.getCenterY();
  }

  organization.centerX = sumCenterX / nodes.length;
  organization.centerY = sumCenterY / nodes.length;

  // Create the organization -> tile members
  for (var i = 0; i < nodes.length; i++) {
    var lNode = nodes[i];

    if (organization.rows.length == 0) {
      this.insertNodeToRow(organization, lNode, 0, minWidth);
    } else if (this.canAddHorizontal(organization, lNode.rect.width, lNode.rect.height)) {
      var rowIndex = organization.rows.length - 1;
      if (!organization.idealRowWidth) {
        rowIndex = this.getShortestRowIndex(organization);
      }
      this.insertNodeToRow(organization, lNode, rowIndex, minWidth);
    } else {
      this.insertNodeToRow(organization, lNode, organization.rows.length, minWidth);
    }

    this.shiftToLastRow(organization);
  }

  return organization;
};

CoSELayout.prototype.insertNodeToRow = function (organization, node, rowIndex, minWidth) {
  var minCompoundSize = minWidth;

  // Add new row if needed
  if (rowIndex == organization.rows.length) {
    var secondDimension = [];

    organization.rows.push(secondDimension);
    organization.rowWidth.push(minCompoundSize);
    organization.rowHeight.push(0);
  }

  // Update row width
  var w = organization.rowWidth[rowIndex] + node.rect.width;

  if (organization.rows[rowIndex].length > 0) {
    w += organization.horizontalPadding;
  }

  organization.rowWidth[rowIndex] = w;
  // Update compound width
  if (organization.width < w) {
    organization.width = w;
  }

  // Update height
  var h = node.rect.height;
  if (rowIndex > 0) h += organization.verticalPadding;

  var extraHeight = 0;
  if (h > organization.rowHeight[rowIndex]) {
    extraHeight = organization.rowHeight[rowIndex];
    organization.rowHeight[rowIndex] = h;
    extraHeight = organization.rowHeight[rowIndex] - extraHeight;
  }

  organization.height += extraHeight;

  // Insert node
  organization.rows[rowIndex].push(node);
};

//Scans the rows of an organization and returns the one with the min width
CoSELayout.prototype.getShortestRowIndex = function (organization) {
  var r = -1;
  var min = Number.MAX_VALUE;

  for (var i = 0; i < organization.rows.length; i++) {
    if (organization.rowWidth[i] < min) {
      r = i;
      min = organization.rowWidth[i];
    }
  }
  return r;
};

//Scans the rows of an organization and returns the one with the max width
CoSELayout.prototype.getLongestRowIndex = function (organization) {
  var r = -1;
  var max = Number.MIN_VALUE;

  for (var i = 0; i < organization.rows.length; i++) {

    if (organization.rowWidth[i] > max) {
      r = i;
      max = organization.rowWidth[i];
    }
  }

  return r;
};

/**
* This method checks whether adding extra width to the organization violates
* the aspect ratio(1) or not.
*/
CoSELayout.prototype.canAddHorizontal = function (organization, extraWidth, extraHeight) {

  // if there is an ideal row width specified use it instead of checking the aspect ratio
  if (organization.idealRowWidth) {
    var lastRowIndex = organization.rows.length - 1;
    var lastRowWidth = organization.rowWidth[lastRowIndex];

    // check and return if ideal row width will be exceed if the node is added to the row
    return lastRowWidth + extraWidth + organization.horizontalPadding <= organization.idealRowWidth;
  }

  var sri = this.getShortestRowIndex(organization);

  if (sri < 0) {
    return true;
  }

  var min = organization.rowWidth[sri];

  if (min + organization.horizontalPadding + extraWidth <= organization.width) return true;

  var hDiff = 0;

  // Adding to an existing row
  if (organization.rowHeight[sri] < extraHeight) {
    if (sri > 0) hDiff = extraHeight + organization.verticalPadding - organization.rowHeight[sri];
  }

  var add_to_row_ratio;
  if (organization.width - min >= extraWidth + organization.horizontalPadding) {
    add_to_row_ratio = (organization.height + hDiff) / (min + extraWidth + organization.horizontalPadding);
  } else {
    add_to_row_ratio = (organization.height + hDiff) / organization.width;
  }

  // Adding a new row for this node
  hDiff = extraHeight + organization.verticalPadding;
  var add_new_row_ratio;
  if (organization.width < extraWidth) {
    add_new_row_ratio = (organization.height + hDiff) / extraWidth;
  } else {
    add_new_row_ratio = (organization.height + hDiff) / organization.width;
  }

  if (add_new_row_ratio < 1) add_new_row_ratio = 1 / add_new_row_ratio;

  if (add_to_row_ratio < 1) add_to_row_ratio = 1 / add_to_row_ratio;

  return add_to_row_ratio < add_new_row_ratio;
};

//If moving the last node from the longest row and adding it to the last
//row makes the bounding box smaller, do it.
CoSELayout.prototype.shiftToLastRow = function (organization) {
  var longest = this.getLongestRowIndex(organization);
  var last = organization.rowWidth.length - 1;
  var row = organization.rows[longest];
  var node = row[row.length - 1];

  var diff = node.width + organization.horizontalPadding;

  // Check if there is enough space on the last row
  if (organization.width - organization.rowWidth[last] > diff && longest != last) {
    // Remove the last element of the longest row
    row.splice(-1, 1);

    // Push it to the last row
    organization.rows[last].push(node);

    organization.rowWidth[longest] = organization.rowWidth[longest] - diff;
    organization.rowWidth[last] = organization.rowWidth[last] + diff;
    organization.width = organization.rowWidth[instance.getLongestRowIndex(organization)];

    // Update heights of the organization
    var maxHeight = Number.MIN_VALUE;
    for (var i = 0; i < row.length; i++) {
      if (row[i].height > maxHeight) maxHeight = row[i].height;
    }
    if (longest > 0) maxHeight += organization.verticalPadding;

    var prevTotal = organization.rowHeight[longest] + organization.rowHeight[last];

    organization.rowHeight[longest] = maxHeight;
    if (organization.rowHeight[last] < node.height + organization.verticalPadding) organization.rowHeight[last] = node.height + organization.verticalPadding;

    var finalTotal = organization.rowHeight[longest] + organization.rowHeight[last];
    organization.height += finalTotal - prevTotal;

    this.shiftToLastRow(organization);
  }
};

CoSELayout.prototype.tilingPreLayout = function () {
  if (CoSEConstants.TILE) {
    // Find zero degree nodes and create a compound for each level
    this.groupZeroDegreeMembers();
    // Tile and clear children of each compound
    this.clearCompounds();
    // Separately tile and clear zero degree nodes for each level
    this.clearZeroDegreeMembers();
  }
};

CoSELayout.prototype.tilingPostLayout = function () {
  if (CoSEConstants.TILE) {
    this.repopulateZeroDegreeMembers();
    this.repopulateCompounds();
  }
};

// -----------------------------------------------------------------------------
// Section: Tree Reduction methods
// -----------------------------------------------------------------------------
// Reduce trees 
CoSELayout.prototype.reduceTrees = function () {
  var prunedNodesAll = [];
  var containsLeaf = true;
  var node;

  while (containsLeaf) {
    var allNodes = this.graphManager.getAllNodes();
    var prunedNodesInStepTemp = [];
    containsLeaf = false;

    for (var i = 0; i < allNodes.length; i++) {
      node = allNodes[i];
      if (node.getEdges().length == 1 && !node.getEdges()[0].isInterGraph && node.getChild() == null) {
        if (CoSEConstants.PURE_INCREMENTAL) {
          var otherEnd = node.getEdges()[0].getOtherEnd(node);
          var relativePosition = new DimensionD(node.getCenterX() - otherEnd.getCenterX(), node.getCenterY() - otherEnd.getCenterY());
          prunedNodesInStepTemp.push([node, node.getEdges()[0], node.getOwner(), relativePosition]);
        } else {
          prunedNodesInStepTemp.push([node, node.getEdges()[0], node.getOwner()]);
        }
        containsLeaf = true;
      }
    }
    if (containsLeaf == true) {
      var prunedNodesInStep = [];
      for (var j = 0; j < prunedNodesInStepTemp.length; j++) {
        if (prunedNodesInStepTemp[j][0].getEdges().length == 1) {
          prunedNodesInStep.push(prunedNodesInStepTemp[j]);
          prunedNodesInStepTemp[j][0].getOwner().remove(prunedNodesInStepTemp[j][0]);
        }
      }
      prunedNodesAll.push(prunedNodesInStep);
      this.graphManager.resetAllNodes();
      this.graphManager.resetAllEdges();
    }
  }
  this.prunedNodesAll = prunedNodesAll;
};

// Grow tree one step 
CoSELayout.prototype.growTree = function (prunedNodesAll) {
  var lengthOfPrunedNodesInStep = prunedNodesAll.length;
  var prunedNodesInStep = prunedNodesAll[lengthOfPrunedNodesInStep - 1];

  var nodeData;
  for (var i = 0; i < prunedNodesInStep.length; i++) {
    nodeData = prunedNodesInStep[i];

    this.findPlaceforPrunedNode(nodeData);

    nodeData[2].add(nodeData[0]);
    nodeData[2].add(nodeData[1], nodeData[1].source, nodeData[1].target);
  }

  prunedNodesAll.splice(prunedNodesAll.length - 1, 1);
  this.graphManager.resetAllNodes();
  this.graphManager.resetAllEdges();
};

// Find an appropriate position to replace pruned node, this method can be improved
CoSELayout.prototype.findPlaceforPrunedNode = function (nodeData) {

  var gridForPrunedNode;
  var nodeToConnect;
  var prunedNode = nodeData[0];
  if (prunedNode == nodeData[1].source) {
    nodeToConnect = nodeData[1].target;
  } else {
    nodeToConnect = nodeData[1].source;
  }

  if (CoSEConstants.PURE_INCREMENTAL) {
    prunedNode.setCenter(nodeToConnect.getCenterX() + nodeData[3].getWidth(), nodeToConnect.getCenterY() + nodeData[3].getHeight());
  } else {
    var startGridX = nodeToConnect.startX;
    var finishGridX = nodeToConnect.finishX;
    var startGridY = nodeToConnect.startY;
    var finishGridY = nodeToConnect.finishY;

    var upNodeCount = 0;
    var downNodeCount = 0;
    var rightNodeCount = 0;
    var leftNodeCount = 0;
    var controlRegions = [upNodeCount, rightNodeCount, downNodeCount, leftNodeCount];

    if (startGridY > 0) {
      for (var i = startGridX; i <= finishGridX; i++) {
        controlRegions[0] += this.grid[i][startGridY - 1].length + this.grid[i][startGridY].length - 1;
      }
    }
    if (finishGridX < this.grid.length - 1) {
      for (var i = startGridY; i <= finishGridY; i++) {
        controlRegions[1] += this.grid[finishGridX + 1][i].length + this.grid[finishGridX][i].length - 1;
      }
    }
    if (finishGridY < this.grid[0].length - 1) {
      for (var i = startGridX; i <= finishGridX; i++) {
        controlRegions[2] += this.grid[i][finishGridY + 1].length + this.grid[i][finishGridY].length - 1;
      }
    }
    if (startGridX > 0) {
      for (var i = startGridY; i <= finishGridY; i++) {
        controlRegions[3] += this.grid[startGridX - 1][i].length + this.grid[startGridX][i].length - 1;
      }
    }
    var min = Integer.MAX_VALUE;
    var minCount;
    var minIndex;
    for (var j = 0; j < controlRegions.length; j++) {
      if (controlRegions[j] < min) {
        min = controlRegions[j];
        minCount = 1;
        minIndex = j;
      } else if (controlRegions[j] == min) {
        minCount++;
      }
    }

    if (minCount == 3 && min == 0) {
      if (controlRegions[0] == 0 && controlRegions[1] == 0 && controlRegions[2] == 0) {
        gridForPrunedNode = 1;
      } else if (controlRegions[0] == 0 && controlRegions[1] == 0 && controlRegions[3] == 0) {
        gridForPrunedNode = 0;
      } else if (controlRegions[0] == 0 && controlRegions[2] == 0 && controlRegions[3] == 0) {
        gridForPrunedNode = 3;
      } else if (controlRegions[1] == 0 && controlRegions[2] == 0 && controlRegions[3] == 0) {
        gridForPrunedNode = 2;
      }
    } else if (minCount == 2 && min == 0) {
      var random = Math.floor(Math.random() * 2);
      if (controlRegions[0] == 0 && controlRegions[1] == 0) {
        ;
        if (random == 0) {
          gridForPrunedNode = 0;
        } else {
          gridForPrunedNode = 1;
        }
      } else if (controlRegions[0] == 0 && controlRegions[2] == 0) {
        if (random == 0) {
          gridForPrunedNode = 0;
        } else {
          gridForPrunedNode = 2;
        }
      } else if (controlRegions[0] == 0 && controlRegions[3] == 0) {
        if (random == 0) {
          gridForPrunedNode = 0;
        } else {
          gridForPrunedNode = 3;
        }
      } else if (controlRegions[1] == 0 && controlRegions[2] == 0) {
        if (random == 0) {
          gridForPrunedNode = 1;
        } else {
          gridForPrunedNode = 2;
        }
      } else if (controlRegions[1] == 0 && controlRegions[3] == 0) {
        if (random == 0) {
          gridForPrunedNode = 1;
        } else {
          gridForPrunedNode = 3;
        }
      } else {
        if (random == 0) {
          gridForPrunedNode = 2;
        } else {
          gridForPrunedNode = 3;
        }
      }
    } else if (minCount == 4 && min == 0) {
      var random = Math.floor(Math.random() * 4);
      gridForPrunedNode = random;
    } else {
      gridForPrunedNode = minIndex;
    }

    if (gridForPrunedNode == 0) {
      prunedNode.setCenter(nodeToConnect.getCenterX(), nodeToConnect.getCenterY() - nodeToConnect.getHeight() / 2 - FDLayoutConstants.DEFAULT_EDGE_LENGTH - prunedNode.getHeight() / 2);
    } else if (gridForPrunedNode == 1) {
      prunedNode.setCenter(nodeToConnect.getCenterX() + nodeToConnect.getWidth() / 2 + FDLayoutConstants.DEFAULT_EDGE_LENGTH + prunedNode.getWidth() / 2, nodeToConnect.getCenterY());
    } else if (gridForPrunedNode == 2) {
      prunedNode.setCenter(nodeToConnect.getCenterX(), nodeToConnect.getCenterY() + nodeToConnect.getHeight() / 2 + FDLayoutConstants.DEFAULT_EDGE_LENGTH + prunedNode.getHeight() / 2);
    } else {
      prunedNode.setCenter(nodeToConnect.getCenterX() - nodeToConnect.getWidth() / 2 - FDLayoutConstants.DEFAULT_EDGE_LENGTH - prunedNode.getWidth() / 2, nodeToConnect.getCenterY());
    }
  }
};

module.exports = CoSELayout;

/***/ }),

/***/ 991:
/***/ ((module, __unused_webpack_exports, __nested_webpack_require_74190__) => {



var FDLayoutNode = __nested_webpack_require_74190__(551).FDLayoutNode;
var IMath = __nested_webpack_require_74190__(551).IMath;

function CoSENode(gm, loc, size, vNode) {
  FDLayoutNode.call(this, gm, loc, size, vNode);
}

CoSENode.prototype = Object.create(FDLayoutNode.prototype);
for (var prop in FDLayoutNode) {
  CoSENode[prop] = FDLayoutNode[prop];
}

CoSENode.prototype.calculateDisplacement = function () {
  var layout = this.graphManager.getLayout();
  // this check is for compound nodes that contain fixed nodes
  if (this.getChild() != null && this.fixedNodeWeight) {
    this.displacementX += layout.coolingFactor * (this.springForceX + this.repulsionForceX + this.gravitationForceX) / this.fixedNodeWeight;
    this.displacementY += layout.coolingFactor * (this.springForceY + this.repulsionForceY + this.gravitationForceY) / this.fixedNodeWeight;
  } else {
    this.displacementX += layout.coolingFactor * (this.springForceX + this.repulsionForceX + this.gravitationForceX) / this.noOfChildren;
    this.displacementY += layout.coolingFactor * (this.springForceY + this.repulsionForceY + this.gravitationForceY) / this.noOfChildren;
  }

  if (Math.abs(this.displacementX) > layout.coolingFactor * layout.maxNodeDisplacement) {
    this.displacementX = layout.coolingFactor * layout.maxNodeDisplacement * IMath.sign(this.displacementX);
  }

  if (Math.abs(this.displacementY) > layout.coolingFactor * layout.maxNodeDisplacement) {
    this.displacementY = layout.coolingFactor * layout.maxNodeDisplacement * IMath.sign(this.displacementY);
  }

  // non-empty compound node, propogate movement to children as well
  if (this.child && this.child.getNodes().length > 0) {
    this.propogateDisplacementToChildren(this.displacementX, this.displacementY);
  }
};

CoSENode.prototype.propogateDisplacementToChildren = function (dX, dY) {
  var nodes = this.getChild().getNodes();
  var node;
  for (var i = 0; i < nodes.length; i++) {
    node = nodes[i];
    if (node.getChild() == null) {
      node.displacementX += dX;
      node.displacementY += dY;
    } else {
      node.propogateDisplacementToChildren(dX, dY);
    }
  }
};

CoSENode.prototype.move = function () {
  var layout = this.graphManager.getLayout();

  // a simple node or an empty compound node, move it
  if (this.child == null || this.child.getNodes().length == 0) {
    this.moveBy(this.displacementX, this.displacementY);

    layout.totalDisplacement += Math.abs(this.displacementX) + Math.abs(this.displacementY);
  }

  this.springForceX = 0;
  this.springForceY = 0;
  this.repulsionForceX = 0;
  this.repulsionForceY = 0;
  this.gravitationForceX = 0;
  this.gravitationForceY = 0;
  this.displacementX = 0;
  this.displacementY = 0;
};

CoSENode.prototype.setPred1 = function (pred1) {
  this.pred1 = pred1;
};

CoSENode.prototype.getPred1 = function () {
  return pred1;
};

CoSENode.prototype.getPred2 = function () {
  return pred2;
};

CoSENode.prototype.setNext = function (next) {
  this.next = next;
};

CoSENode.prototype.getNext = function () {
  return next;
};

CoSENode.prototype.setProcessed = function (processed) {
  this.processed = processed;
};

CoSENode.prototype.isProcessed = function () {
  return processed;
};

module.exports = CoSENode;

/***/ }),

/***/ 902:
/***/ ((module, __unused_webpack_exports, __nested_webpack_require_77514__) => {



function _toConsumableArray(arr) { if (Array.isArray(arr)) { for (var i = 0, arr2 = Array(arr.length); i < arr.length; i++) { arr2[i] = arr[i]; } return arr2; } else { return Array.from(arr); } }

var CoSEConstants = __nested_webpack_require_77514__(806);
var LinkedList = __nested_webpack_require_77514__(551).LinkedList;
var Matrix = __nested_webpack_require_77514__(551).Matrix;
var SVD = __nested_webpack_require_77514__(551).SVD;

function ConstraintHandler() {}

ConstraintHandler.handleConstraints = function (layout) {
  //  let layout = this.graphManager.getLayout();

  // get constraints from layout
  var constraints = {};
  constraints.fixedNodeConstraint = layout.constraints.fixedNodeConstraint;
  constraints.alignmentConstraint = layout.constraints.alignmentConstraint;
  constraints.relativePlacementConstraint = layout.constraints.relativePlacementConstraint;

  var idToNodeMap = new Map();
  var nodeIndexes = new Map();
  var xCoords = [];
  var yCoords = [];

  var allNodes = layout.getAllNodes();
  var index = 0;
  // fill index map and coordinates
  for (var i = 0; i < allNodes.length; i++) {
    var node = allNodes[i];
    if (node.getChild() == null) {
      nodeIndexes.set(node.id, index++);
      xCoords.push(node.getCenterX());
      yCoords.push(node.getCenterY());
      idToNodeMap.set(node.id, node);
    }
  }

  // if there exists relative placement constraint without gap value, set it to default 
  if (constraints.relativePlacementConstraint) {
    constraints.relativePlacementConstraint.forEach(function (constraint) {
      if (!constraint.gap && constraint.gap != 0) {
        if (constraint.left) {
          constraint.gap = CoSEConstants.DEFAULT_EDGE_LENGTH + idToNodeMap.get(constraint.left).getWidth() / 2 + idToNodeMap.get(constraint.right).getWidth() / 2;
        } else {
          constraint.gap = CoSEConstants.DEFAULT_EDGE_LENGTH + idToNodeMap.get(constraint.top).getHeight() / 2 + idToNodeMap.get(constraint.bottom).getHeight() / 2;
        }
      }
    });
  }

  /* auxiliary functions */

  // calculate difference between two position objects
  var calculatePositionDiff = function calculatePositionDiff(pos1, pos2) {
    return { x: pos1.x - pos2.x, y: pos1.y - pos2.y };
  };

  // calculate average position of the nodes
  var calculateAvgPosition = function calculateAvgPosition(nodeIdSet) {
    var xPosSum = 0;
    var yPosSum = 0;
    nodeIdSet.forEach(function (nodeId) {
      xPosSum += xCoords[nodeIndexes.get(nodeId)];
      yPosSum += yCoords[nodeIndexes.get(nodeId)];
    });

    return { x: xPosSum / nodeIdSet.size, y: yPosSum / nodeIdSet.size };
  };

  // find an appropriate positioning for the nodes in a given graph according to relative placement constraints
  // this function also takes the fixed nodes and alignment constraints into account
  // graph: dag to be evaluated, direction: "horizontal" or "vertical", 
  // fixedNodes: set of fixed nodes to consider during evaluation, dummyPositions: appropriate coordinates of the dummy nodes  
  var findAppropriatePositionForRelativePlacement = function findAppropriatePositionForRelativePlacement(graph, direction, fixedNodes, dummyPositions, componentSources) {

    // find union of two sets
    function setUnion(setA, setB) {
      var union = new Set(setA);
      var _iteratorNormalCompletion = true;
      var _didIteratorError = false;
      var _iteratorError = undefined;

      try {
        for (var _iterator = setB[Symbol.iterator](), _step; !(_iteratorNormalCompletion = (_step = _iterator.next()).done); _iteratorNormalCompletion = true) {
          var elem = _step.value;

          union.add(elem);
        }
      } catch (err) {
        _didIteratorError = true;
        _iteratorError = err;
      } finally {
        try {
          if (!_iteratorNormalCompletion && _iterator.return) {
            _iterator.return();
          }
        } finally {
          if (_didIteratorError) {
            throw _iteratorError;
          }
        }
      }

      return union;
    }

    // find indegree count for each node
    var inDegrees = new Map();

    graph.forEach(function (value, key) {
      inDegrees.set(key, 0);
    });
    graph.forEach(function (value, key) {
      value.forEach(function (adjacent) {
        inDegrees.set(adjacent.id, inDegrees.get(adjacent.id) + 1);
      });
    });

    var positionMap = new Map(); // keeps the position for each node
    var pastMap = new Map(); // keeps the predecessors(past) of a node
    var queue = new LinkedList();
    inDegrees.forEach(function (value, key) {
      if (value == 0) {
        queue.push(key);
        if (!fixedNodes) {
          if (direction == "horizontal") {
            positionMap.set(key, nodeIndexes.has(key) ? xCoords[nodeIndexes.get(key)] : dummyPositions.get(key));
          } else {
            positionMap.set(key, nodeIndexes.has(key) ? yCoords[nodeIndexes.get(key)] : dummyPositions.get(key));
          }
        }
      } else {
        positionMap.set(key, Number.NEGATIVE_INFINITY);
      }
      if (fixedNodes) {
        pastMap.set(key, new Set([key]));
      }
    });

    // align sources of each component in enforcement phase
    if (fixedNodes) {
      componentSources.forEach(function (component) {
        var fixedIds = [];
        component.forEach(function (nodeId) {
          if (fixedNodes.has(nodeId)) {
            fixedIds.push(nodeId);
          }
        });
        if (fixedIds.length > 0) {
          var position = 0;
          fixedIds.forEach(function (fixedId) {
            if (direction == "horizontal") {
              positionMap.set(fixedId, nodeIndexes.has(fixedId) ? xCoords[nodeIndexes.get(fixedId)] : dummyPositions.get(fixedId));
              position += positionMap.get(fixedId);
            } else {
              positionMap.set(fixedId, nodeIndexes.has(fixedId) ? yCoords[nodeIndexes.get(fixedId)] : dummyPositions.get(fixedId));
              position += positionMap.get(fixedId);
            }
          });
          position = position / fixedIds.length;
          component.forEach(function (nodeId) {
            if (!fixedNodes.has(nodeId)) {
              positionMap.set(nodeId, position);
            }
          });
        } else {
          var _position = 0;
          component.forEach(function (nodeId) {
            if (direction == "horizontal") {
              _position += nodeIndexes.has(nodeId) ? xCoords[nodeIndexes.get(nodeId)] : dummyPositions.get(nodeId);
            } else {
              _position += nodeIndexes.has(nodeId) ? yCoords[nodeIndexes.get(nodeId)] : dummyPositions.get(nodeId);
            }
          });
          _position = _position / component.length;
          component.forEach(function (nodeId) {
            positionMap.set(nodeId, _position);
          });
        }
      });
    }

    // calculate positions of the nodes

    var _loop = function _loop() {
      var currentNode = queue.shift();
      var neighbors = graph.get(currentNode);
      neighbors.forEach(function (neighbor) {
        if (positionMap.get(neighbor.id) < positionMap.get(currentNode) + neighbor.gap) {
          if (fixedNodes && fixedNodes.has(neighbor.id)) {
            var fixedPosition = void 0;
            if (direction == "horizontal") {
              fixedPosition = nodeIndexes.has(neighbor.id) ? xCoords[nodeIndexes.get(neighbor.id)] : dummyPositions.get(neighbor.id);
            } else {
              fixedPosition = nodeIndexes.has(neighbor.id) ? yCoords[nodeIndexes.get(neighbor.id)] : dummyPositions.get(neighbor.id);
            }
            positionMap.set(neighbor.id, fixedPosition); // TODO: may do unnecessary work
            if (fixedPosition < positionMap.get(currentNode) + neighbor.gap) {
              var diff = positionMap.get(currentNode) + neighbor.gap - fixedPosition;
              pastMap.get(currentNode).forEach(function (nodeId) {
                positionMap.set(nodeId, positionMap.get(nodeId) - diff);
              });
            }
          } else {
            positionMap.set(neighbor.id, positionMap.get(currentNode) + neighbor.gap);
          }
        }
        inDegrees.set(neighbor.id, inDegrees.get(neighbor.id) - 1);
        if (inDegrees.get(neighbor.id) == 0) {
          queue.push(neighbor.id);
        }
        if (fixedNodes) {
          pastMap.set(neighbor.id, setUnion(pastMap.get(currentNode), pastMap.get(neighbor.id)));
        }
      });
    };

    while (queue.length != 0) {
      _loop();
    }

    // readjust position of the nodes after enforcement
    if (fixedNodes) {
      // find indegree count for each node
      var sinkNodes = new Set();

      graph.forEach(function (value, key) {
        if (value.length == 0) {
          sinkNodes.add(key);
        }
      });

      var _components = [];
      pastMap.forEach(function (value, key) {
        if (sinkNodes.has(key)) {
          var isFixedComponent = false;
          var _iteratorNormalCompletion2 = true;
          var _didIteratorError2 = false;
          var _iteratorError2 = undefined;

          try {
            for (var _iterator2 = value[Symbol.iterator](), _step2; !(_iteratorNormalCompletion2 = (_step2 = _iterator2.next()).done); _iteratorNormalCompletion2 = true) {
              var nodeId = _step2.value;

              if (fixedNodes.has(nodeId)) {
                isFixedComponent = true;
              }
            }
          } catch (err) {
            _didIteratorError2 = true;
            _iteratorError2 = err;
          } finally {
            try {
              if (!_iteratorNormalCompletion2 && _iterator2.return) {
                _iterator2.return();
              }
            } finally {
              if (_didIteratorError2) {
                throw _iteratorError2;
              }
            }
          }

          if (!isFixedComponent) {
            var isExist = false;
            var existAt = void 0;
            _components.forEach(function (component, index) {
              if (component.has([].concat(_toConsumableArray(value))[0])) {
                isExist = true;
                existAt = index;
              }
            });
            if (!isExist) {
              _components.push(new Set(value));
            } else {
              value.forEach(function (ele) {
                _components[existAt].add(ele);
              });
            }
          }
        }
      });

      _components.forEach(function (component, index) {
        var minBefore = Number.POSITIVE_INFINITY;
        var minAfter = Number.POSITIVE_INFINITY;
        var maxBefore = Number.NEGATIVE_INFINITY;
        var maxAfter = Number.NEGATIVE_INFINITY;

        var _iteratorNormalCompletion3 = true;
        var _didIteratorError3 = false;
        var _iteratorError3 = undefined;

        try {
          for (var _iterator3 = component[Symbol.iterator](), _step3; !(_iteratorNormalCompletion3 = (_step3 = _iterator3.next()).done); _iteratorNormalCompletion3 = true) {
            var nodeId = _step3.value;

            var posBefore = void 0;
            if (direction == "horizontal") {
              posBefore = nodeIndexes.has(nodeId) ? xCoords[nodeIndexes.get(nodeId)] : dummyPositions.get(nodeId);
            } else {
              posBefore = nodeIndexes.has(nodeId) ? yCoords[nodeIndexes.get(nodeId)] : dummyPositions.get(nodeId);
            }
            var posAfter = positionMap.get(nodeId);
            if (posBefore < minBefore) {
              minBefore = posBefore;
            }
            if (posBefore > maxBefore) {
              maxBefore = posBefore;
            }
            if (posAfter < minAfter) {
              minAfter = posAfter;
            }
            if (posAfter > maxAfter) {
              maxAfter = posAfter;
            }
          }
        } catch (err) {
          _didIteratorError3 = true;
          _iteratorError3 = err;
        } finally {
          try {
            if (!_iteratorNormalCompletion3 && _iterator3.return) {
              _iterator3.return();
            }
          } finally {
            if (_didIteratorError3) {
              throw _iteratorError3;
            }
          }
        }

        var diff = (minBefore + maxBefore) / 2 - (minAfter + maxAfter) / 2;

        var _iteratorNormalCompletion4 = true;
        var _didIteratorError4 = false;
        var _iteratorError4 = undefined;

        try {
          for (var _iterator4 = component[Symbol.iterator](), _step4; !(_iteratorNormalCompletion4 = (_step4 = _iterator4.next()).done); _iteratorNormalCompletion4 = true) {
            var _nodeId = _step4.value;

            positionMap.set(_nodeId, positionMap.get(_nodeId) + diff);
          }
        } catch (err) {
          _didIteratorError4 = true;
          _iteratorError4 = err;
        } finally {
          try {
            if (!_iteratorNormalCompletion4 && _iterator4.return) {
              _iterator4.return();
            }
          } finally {
            if (_didIteratorError4) {
              throw _iteratorError4;
            }
          }
        }
      });
    }

    return positionMap;
  };

  // find transformation based on rel. placement constraints if there are both alignment and rel. placement constraints
  // or if there are only rel. placement contraints where the largest component isn't sufficiently large
  var applyReflectionForRelativePlacement = function applyReflectionForRelativePlacement(relativePlacementConstraints) {
    // variables to count votes
    var reflectOnY = 0,
        notReflectOnY = 0;
    var reflectOnX = 0,
        notReflectOnX = 0;

    relativePlacementConstraints.forEach(function (constraint) {
      if (constraint.left) {
        xCoords[nodeIndexes.get(constraint.left)] - xCoords[nodeIndexes.get(constraint.right)] >= 0 ? reflectOnY++ : notReflectOnY++;
      } else {
        yCoords[nodeIndexes.get(constraint.top)] - yCoords[nodeIndexes.get(constraint.bottom)] >= 0 ? reflectOnX++ : notReflectOnX++;
      }
    });

    if (reflectOnY > notReflectOnY && reflectOnX > notReflectOnX) {
      for (var _i = 0; _i < nodeIndexes.size; _i++) {
        xCoords[_i] = -1 * xCoords[_i];
        yCoords[_i] = -1 * yCoords[_i];
      }
    } else if (reflectOnY > notReflectOnY) {
      for (var _i2 = 0; _i2 < nodeIndexes.size; _i2++) {
        xCoords[_i2] = -1 * xCoords[_i2];
      }
    } else if (reflectOnX > notReflectOnX) {
      for (var _i3 = 0; _i3 < nodeIndexes.size; _i3++) {
        yCoords[_i3] = -1 * yCoords[_i3];
      }
    }
  };

  // find weakly connected components in undirected graph
  var findComponents = function findComponents(graph) {
    // find weakly connected components in dag
    var components = [];
    var queue = new LinkedList();
    var visited = new Set();
    var count = 0;

    graph.forEach(function (value, key) {
      if (!visited.has(key)) {
        components[count] = [];
        var _currentNode = key;
        queue.push(_currentNode);
        visited.add(_currentNode);
        components[count].push(_currentNode);

        while (queue.length != 0) {
          _currentNode = queue.shift();
          var neighbors = graph.get(_currentNode);
          neighbors.forEach(function (neighbor) {
            if (!visited.has(neighbor.id)) {
              queue.push(neighbor.id);
              visited.add(neighbor.id);
              components[count].push(neighbor.id);
            }
          });
        }
        count++;
      }
    });
    return components;
  };

  // return undirected version of given dag
  var dagToUndirected = function dagToUndirected(dag) {
    var undirected = new Map();

    dag.forEach(function (value, key) {
      undirected.set(key, []);
    });

    dag.forEach(function (value, key) {
      value.forEach(function (adjacent) {
        undirected.get(key).push(adjacent);
        undirected.get(adjacent.id).push({ id: key, gap: adjacent.gap, direction: adjacent.direction });
      });
    });

    return undirected;
  };

  // return reversed (directions inverted) version of given dag
  var dagToReversed = function dagToReversed(dag) {
    var reversed = new Map();

    dag.forEach(function (value, key) {
      reversed.set(key, []);
    });

    dag.forEach(function (value, key) {
      value.forEach(function (adjacent) {
        reversed.get(adjacent.id).push({ id: key, gap: adjacent.gap, direction: adjacent.direction });
      });
    });

    return reversed;
  };

  /****  apply transformation to the initial draft layout to better align with constrained nodes ****/
  // solve the Orthogonal Procrustean Problem to rotate and/or reflect initial draft layout
  // here we follow the solution in Chapter 20.2 of Borg, I. & Groenen, P. (2005) Modern Multidimensional Scaling: Theory and Applications 

  /* construct source and target configurations */

  var targetMatrix = []; // A - target configuration
  var sourceMatrix = []; // B - source configuration 
  var standardTransformation = false; // false for no transformation, true for standart (Procrustes) transformation (rotation and/or reflection)
  var reflectionType = false; // false/true for reflection check, 'reflectOnX', 'reflectOnY' or 'reflectOnBoth' for reflection type if necessary
  var fixedNodes = new Set();
  var dag = new Map(); // adjacency list to keep directed acyclic graph (dag) that consists of relative placement constraints
  var dagUndirected = new Map(); // undirected version of the dag
  var components = []; // weakly connected components

  // fill fixedNodes collection to use later
  if (constraints.fixedNodeConstraint) {
    constraints.fixedNodeConstraint.forEach(function (nodeData) {
      fixedNodes.add(nodeData.nodeId);
    });
  }

  // construct dag from relative placement constraints 
  if (constraints.relativePlacementConstraint) {
    // construct both directed and undirected version of the dag
    constraints.relativePlacementConstraint.forEach(function (constraint) {
      if (constraint.left) {
        if (dag.has(constraint.left)) {
          dag.get(constraint.left).push({ id: constraint.right, gap: constraint.gap, direction: "horizontal" });
        } else {
          dag.set(constraint.left, [{ id: constraint.right, gap: constraint.gap, direction: "horizontal" }]);
        }
        if (!dag.has(constraint.right)) {
          dag.set(constraint.right, []);
        }
      } else {
        if (dag.has(constraint.top)) {
          dag.get(constraint.top).push({ id: constraint.bottom, gap: constraint.gap, direction: "vertical" });
        } else {
          dag.set(constraint.top, [{ id: constraint.bottom, gap: constraint.gap, direction: "vertical" }]);
        }
        if (!dag.has(constraint.bottom)) {
          dag.set(constraint.bottom, []);
        }
      }
    });

    dagUndirected = dagToUndirected(dag);
    components = findComponents(dagUndirected);
  }

  if (CoSEConstants.TRANSFORM_ON_CONSTRAINT_HANDLING) {
    // first check fixed node constraint
    if (constraints.fixedNodeConstraint && constraints.fixedNodeConstraint.length > 1) {
      constraints.fixedNodeConstraint.forEach(function (nodeData, i) {
        targetMatrix[i] = [nodeData.position.x, nodeData.position.y];
        sourceMatrix[i] = [xCoords[nodeIndexes.get(nodeData.nodeId)], yCoords[nodeIndexes.get(nodeData.nodeId)]];
      });
      standardTransformation = true;
    } else if (constraints.alignmentConstraint) {
      (function () {
        // then check alignment constraint
        var count = 0;
        if (constraints.alignmentConstraint.vertical) {
          var verticalAlign = constraints.alignmentConstraint.vertical;

          var _loop2 = function _loop2(_i4) {
            var alignmentSet = new Set();
            verticalAlign[_i4].forEach(function (nodeId) {
              alignmentSet.add(nodeId);
            });
            var intersection = new Set([].concat(_toConsumableArray(alignmentSet)).filter(function (x) {
              return fixedNodes.has(x);
            }));
            var xPos = void 0;
            if (intersection.size > 0) xPos = xCoords[nodeIndexes.get(intersection.values().next().value)];else xPos = calculateAvgPosition(alignmentSet).x;

            verticalAlign[_i4].forEach(function (nodeId) {
              targetMatrix[count] = [xPos, yCoords[nodeIndexes.get(nodeId)]];
              sourceMatrix[count] = [xCoords[nodeIndexes.get(nodeId)], yCoords[nodeIndexes.get(nodeId)]];
              count++;
            });
          };

          for (var _i4 = 0; _i4 < verticalAlign.length; _i4++) {
            _loop2(_i4);
          }
          standardTransformation = true;
        }
        if (constraints.alignmentConstraint.horizontal) {
          var horizontalAlign = constraints.alignmentConstraint.horizontal;

          var _loop3 = function _loop3(_i5) {
            var alignmentSet = new Set();
            horizontalAlign[_i5].forEach(function (nodeId) {
              alignmentSet.add(nodeId);
            });
            var intersection = new Set([].concat(_toConsumableArray(alignmentSet)).filter(function (x) {
              return fixedNodes.has(x);
            }));
            var yPos = void 0;
            if (intersection.size > 0) yPos = xCoords[nodeIndexes.get(intersection.values().next().value)];else yPos = calculateAvgPosition(alignmentSet).y;

            horizontalAlign[_i5].forEach(function (nodeId) {
              targetMatrix[count] = [xCoords[nodeIndexes.get(nodeId)], yPos];
              sourceMatrix[count] = [xCoords[nodeIndexes.get(nodeId)], yCoords[nodeIndexes.get(nodeId)]];
              count++;
            });
          };

          for (var _i5 = 0; _i5 < horizontalAlign.length; _i5++) {
            _loop3(_i5);
          }
          standardTransformation = true;
        }
        if (constraints.relativePlacementConstraint) {
          reflectionType = true;
        }
      })();
    } else if (constraints.relativePlacementConstraint) {
      // finally check relative placement constraint
      // find largest component in dag
      var largestComponentSize = 0;
      var largestComponentIndex = 0;
      for (var _i6 = 0; _i6 < components.length; _i6++) {
        if (components[_i6].length > largestComponentSize) {
          largestComponentSize = components[_i6].length;
          largestComponentIndex = _i6;
        }
      }
      // if largest component isn't dominant, then take the votes for reflection
      if (largestComponentSize < dagUndirected.size / 2) {
        applyReflectionForRelativePlacement(constraints.relativePlacementConstraint);
        standardTransformation = false;
        reflectionType = false;
      } else {
        // use largest component for transformation
        // construct horizontal and vertical subgraphs in the largest component
        var subGraphOnHorizontal = new Map();
        var subGraphOnVertical = new Map();
        var constraintsInlargestComponent = [];

        components[largestComponentIndex].forEach(function (nodeId) {
          dag.get(nodeId).forEach(function (adjacent) {
            if (adjacent.direction == "horizontal") {
              if (subGraphOnHorizontal.has(nodeId)) {
                subGraphOnHorizontal.get(nodeId).push(adjacent);
              } else {
                subGraphOnHorizontal.set(nodeId, [adjacent]);
              }
              if (!subGraphOnHorizontal.has(adjacent.id)) {
                subGraphOnHorizontal.set(adjacent.id, []);
              }
              constraintsInlargestComponent.push({ left: nodeId, right: adjacent.id });
            } else {
              if (subGraphOnVertical.has(nodeId)) {
                subGraphOnVertical.get(nodeId).push(adjacent);
              } else {
                subGraphOnVertical.set(nodeId, [adjacent]);
              }
              if (!subGraphOnVertical.has(adjacent.id)) {
                subGraphOnVertical.set(adjacent.id, []);
              }
              constraintsInlargestComponent.push({ top: nodeId, bottom: adjacent.id });
            }
          });
        });

        applyReflectionForRelativePlacement(constraintsInlargestComponent);
        reflectionType = false;

        // calculate appropriate positioning for subgraphs
        var positionMapHorizontal = findAppropriatePositionForRelativePlacement(subGraphOnHorizontal, "horizontal");
        var positionMapVertical = findAppropriatePositionForRelativePlacement(subGraphOnVertical, "vertical");

        // construct source and target configuration
        components[largestComponentIndex].forEach(function (nodeId, i) {
          sourceMatrix[i] = [xCoords[nodeIndexes.get(nodeId)], yCoords[nodeIndexes.get(nodeId)]];
          targetMatrix[i] = [];
          if (positionMapHorizontal.has(nodeId)) {
            targetMatrix[i][0] = positionMapHorizontal.get(nodeId);
          } else {
            targetMatrix[i][0] = xCoords[nodeIndexes.get(nodeId)];
          }
          if (positionMapVertical.has(nodeId)) {
            targetMatrix[i][1] = positionMapVertical.get(nodeId);
          } else {
            targetMatrix[i][1] = yCoords[nodeIndexes.get(nodeId)];
          }
        });

        standardTransformation = true;
      }
    }

    // if transformation is required, then calculate and apply transformation matrix
    if (standardTransformation) {
      /* calculate transformation matrix */
      var transformationMatrix = void 0;
      var targetMatrixTranspose = Matrix.transpose(targetMatrix); // A'
      var sourceMatrixTranspose = Matrix.transpose(sourceMatrix); // B'

      // centralize transpose matrices
      for (var _i7 = 0; _i7 < targetMatrixTranspose.length; _i7++) {
        targetMatrixTranspose[_i7] = Matrix.multGamma(targetMatrixTranspose[_i7]);
        sourceMatrixTranspose[_i7] = Matrix.multGamma(sourceMatrixTranspose[_i7]);
      }

      // do actual calculation for transformation matrix
      var tempMatrix = Matrix.multMat(targetMatrixTranspose, Matrix.transpose(sourceMatrixTranspose)); // tempMatrix = A'B
      var SVDResult = SVD.svd(tempMatrix); // SVD(A'B) = USV', svd function returns U, S and V 
      transformationMatrix = Matrix.multMat(SVDResult.V, Matrix.transpose(SVDResult.U)); // transformationMatrix = T = VU'

      /* apply found transformation matrix to obtain final draft layout */
      for (var _i8 = 0; _i8 < nodeIndexes.size; _i8++) {
        var temp1 = [xCoords[_i8], yCoords[_i8]];
        var temp2 = [transformationMatrix[0][0], transformationMatrix[1][0]];
        var temp3 = [transformationMatrix[0][1], transformationMatrix[1][1]];
        xCoords[_i8] = Matrix.dotProduct(temp1, temp2);
        yCoords[_i8] = Matrix.dotProduct(temp1, temp3);
      }

      // applied only both alignment and rel. placement constraints exist
      if (reflectionType) {
        applyReflectionForRelativePlacement(constraints.relativePlacementConstraint);
      }
    }
  }

  if (CoSEConstants.ENFORCE_CONSTRAINTS) {
    /****  enforce constraints on the transformed draft layout ****/

    /* first enforce fixed node constraint */

    if (constraints.fixedNodeConstraint && constraints.fixedNodeConstraint.length > 0) {
      var translationAmount = { x: 0, y: 0 };
      constraints.fixedNodeConstraint.forEach(function (nodeData, i) {
        var posInTheory = { x: xCoords[nodeIndexes.get(nodeData.nodeId)], y: yCoords[nodeIndexes.get(nodeData.nodeId)] };
        var posDesired = nodeData.position;
        var posDiff = calculatePositionDiff(posDesired, posInTheory);
        translationAmount.x += posDiff.x;
        translationAmount.y += posDiff.y;
      });
      translationAmount.x /= constraints.fixedNodeConstraint.length;
      translationAmount.y /= constraints.fixedNodeConstraint.length;

      xCoords.forEach(function (value, i) {
        xCoords[i] += translationAmount.x;
      });

      yCoords.forEach(function (value, i) {
        yCoords[i] += translationAmount.y;
      });

      constraints.fixedNodeConstraint.forEach(function (nodeData) {
        xCoords[nodeIndexes.get(nodeData.nodeId)] = nodeData.position.x;
        yCoords[nodeIndexes.get(nodeData.nodeId)] = nodeData.position.y;
      });
    }

    /* then enforce alignment constraint */

    if (constraints.alignmentConstraint) {
      if (constraints.alignmentConstraint.vertical) {
        var xAlign = constraints.alignmentConstraint.vertical;

        var _loop4 = function _loop4(_i9) {
          var alignmentSet = new Set();
          xAlign[_i9].forEach(function (nodeId) {
            alignmentSet.add(nodeId);
          });
          var intersection = new Set([].concat(_toConsumableArray(alignmentSet)).filter(function (x) {
            return fixedNodes.has(x);
          }));
          var xPos = void 0;
          if (intersection.size > 0) xPos = xCoords[nodeIndexes.get(intersection.values().next().value)];else xPos = calculateAvgPosition(alignmentSet).x;

          alignmentSet.forEach(function (nodeId) {
            if (!fixedNodes.has(nodeId)) xCoords[nodeIndexes.get(nodeId)] = xPos;
          });
        };

        for (var _i9 = 0; _i9 < xAlign.length; _i9++) {
          _loop4(_i9);
        }
      }
      if (constraints.alignmentConstraint.horizontal) {
        var yAlign = constraints.alignmentConstraint.horizontal;

        var _loop5 = function _loop5(_i10) {
          var alignmentSet = new Set();
          yAlign[_i10].forEach(function (nodeId) {
            alignmentSet.add(nodeId);
          });
          var intersection = new Set([].concat(_toConsumableArray(alignmentSet)).filter(function (x) {
            return fixedNodes.has(x);
          }));
          var yPos = void 0;
          if (intersection.size > 0) yPos = yCoords[nodeIndexes.get(intersection.values().next().value)];else yPos = calculateAvgPosition(alignmentSet).y;

          alignmentSet.forEach(function (nodeId) {
            if (!fixedNodes.has(nodeId)) yCoords[nodeIndexes.get(nodeId)] = yPos;
          });
        };

        for (var _i10 = 0; _i10 < yAlign.length; _i10++) {
          _loop5(_i10);
        }
      }
    }

    /* finally enforce relative placement constraint */

    if (constraints.relativePlacementConstraint) {
      (function () {
        var nodeToDummyForVerticalAlignment = new Map();
        var nodeToDummyForHorizontalAlignment = new Map();
        var dummyToNodeForVerticalAlignment = new Map();
        var dummyToNodeForHorizontalAlignment = new Map();
        var dummyPositionsForVerticalAlignment = new Map();
        var dummyPositionsForHorizontalAlignment = new Map();
        var fixedNodesOnHorizontal = new Set();
        var fixedNodesOnVertical = new Set();

        // fill maps and sets      
        fixedNodes.forEach(function (nodeId) {
          fixedNodesOnHorizontal.add(nodeId);
          fixedNodesOnVertical.add(nodeId);
        });

        if (constraints.alignmentConstraint) {
          if (constraints.alignmentConstraint.vertical) {
            var verticalAlignment = constraints.alignmentConstraint.vertical;

            var _loop6 = function _loop6(_i11) {
              dummyToNodeForVerticalAlignment.set("dummy" + _i11, []);
              verticalAlignment[_i11].forEach(function (nodeId) {
                nodeToDummyForVerticalAlignment.set(nodeId, "dummy" + _i11);
                dummyToNodeForVerticalAlignment.get("dummy" + _i11).push(nodeId);
                if (fixedNodes.has(nodeId)) {
                  fixedNodesOnHorizontal.add("dummy" + _i11);
                }
              });
              dummyPositionsForVerticalAlignment.set("dummy" + _i11, xCoords[nodeIndexes.get(verticalAlignment[_i11][0])]);
            };

            for (var _i11 = 0; _i11 < verticalAlignment.length; _i11++) {
              _loop6(_i11);
            }
          }
          if (constraints.alignmentConstraint.horizontal) {
            var horizontalAlignment = constraints.alignmentConstraint.horizontal;

            var _loop7 = function _loop7(_i12) {
              dummyToNodeForHorizontalAlignment.set("dummy" + _i12, []);
              horizontalAlignment[_i12].forEach(function (nodeId) {
                nodeToDummyForHorizontalAlignment.set(nodeId, "dummy" + _i12);
                dummyToNodeForHorizontalAlignment.get("dummy" + _i12).push(nodeId);
                if (fixedNodes.has(nodeId)) {
                  fixedNodesOnVertical.add("dummy" + _i12);
                }
              });
              dummyPositionsForHorizontalAlignment.set("dummy" + _i12, yCoords[nodeIndexes.get(horizontalAlignment[_i12][0])]);
            };

            for (var _i12 = 0; _i12 < horizontalAlignment.length; _i12++) {
              _loop7(_i12);
            }
          }
        }

        // construct horizontal and vertical dags (subgraphs) from overall dag
        var dagOnHorizontal = new Map();
        var dagOnVertical = new Map();

        var _loop8 = function _loop8(nodeId) {
          dag.get(nodeId).forEach(function (adjacent) {
            var sourceId = void 0;
            var targetNode = void 0;
            if (adjacent["direction"] == "horizontal") {
              sourceId = nodeToDummyForVerticalAlignment.get(nodeId) ? nodeToDummyForVerticalAlignment.get(nodeId) : nodeId;
              if (nodeToDummyForVerticalAlignment.get(adjacent.id)) {
                targetNode = { id: nodeToDummyForVerticalAlignment.get(adjacent.id), gap: adjacent.gap, direction: adjacent.direction };
              } else {
                targetNode = adjacent;
              }
              if (dagOnHorizontal.has(sourceId)) {
                dagOnHorizontal.get(sourceId).push(targetNode);
              } else {
                dagOnHorizontal.set(sourceId, [targetNode]);
              }
              if (!dagOnHorizontal.has(targetNode.id)) {
                dagOnHorizontal.set(targetNode.id, []);
              }
            } else {
              sourceId = nodeToDummyForHorizontalAlignment.get(nodeId) ? nodeToDummyForHorizontalAlignment.get(nodeId) : nodeId;
              if (nodeToDummyForHorizontalAlignment.get(adjacent.id)) {
                targetNode = { id: nodeToDummyForHorizontalAlignment.get(adjacent.id), gap: adjacent.gap, direction: adjacent.direction };
              } else {
                targetNode = adjacent;
              }
              if (dagOnVertical.has(sourceId)) {
                dagOnVertical.get(sourceId).push(targetNode);
              } else {
                dagOnVertical.set(sourceId, [targetNode]);
              }
              if (!dagOnVertical.has(targetNode.id)) {
                dagOnVertical.set(targetNode.id, []);
              }
            }
          });
        };

        var _iteratorNormalCompletion5 = true;
        var _didIteratorError5 = false;
        var _iteratorError5 = undefined;

        try {
          for (var _iterator5 = dag.keys()[Symbol.iterator](), _step5; !(_iteratorNormalCompletion5 = (_step5 = _iterator5.next()).done); _iteratorNormalCompletion5 = true) {
            var nodeId = _step5.value;

            _loop8(nodeId);
          }

          // find source nodes of each component in horizontal and vertical dags
        } catch (err) {
          _didIteratorError5 = true;
          _iteratorError5 = err;
        } finally {
          try {
            if (!_iteratorNormalCompletion5 && _iterator5.return) {
              _iterator5.return();
            }
          } finally {
            if (_didIteratorError5) {
              throw _iteratorError5;
            }
          }
        }

        var undirectedOnHorizontal = dagToUndirected(dagOnHorizontal);
        var undirectedOnVertical = dagToUndirected(dagOnVertical);
        var componentsOnHorizontal = findComponents(undirectedOnHorizontal);
        var componentsOnVertical = findComponents(undirectedOnVertical);
        var reversedDagOnHorizontal = dagToReversed(dagOnHorizontal);
        var reversedDagOnVertical = dagToReversed(dagOnVertical);
        var componentSourcesOnHorizontal = [];
        var componentSourcesOnVertical = [];

        componentsOnHorizontal.forEach(function (component, index) {
          componentSourcesOnHorizontal[index] = [];
          component.forEach(function (nodeId) {
            if (reversedDagOnHorizontal.get(nodeId).length == 0) {
              componentSourcesOnHorizontal[index].push(nodeId);
            }
          });
        });

        componentsOnVertical.forEach(function (component, index) {
          componentSourcesOnVertical[index] = [];
          component.forEach(function (nodeId) {
            if (reversedDagOnVertical.get(nodeId).length == 0) {
              componentSourcesOnVertical[index].push(nodeId);
            }
          });
        });

        // calculate appropriate positioning for subgraphs
        var positionMapHorizontal = findAppropriatePositionForRelativePlacement(dagOnHorizontal, "horizontal", fixedNodesOnHorizontal, dummyPositionsForVerticalAlignment, componentSourcesOnHorizontal);
        var positionMapVertical = findAppropriatePositionForRelativePlacement(dagOnVertical, "vertical", fixedNodesOnVertical, dummyPositionsForHorizontalAlignment, componentSourcesOnVertical);

        // update positions of the nodes based on relative placement constraints

        var _loop9 = function _loop9(key) {
          if (dummyToNodeForVerticalAlignment.get(key)) {
            dummyToNodeForVerticalAlignment.get(key).forEach(function (nodeId) {
              xCoords[nodeIndexes.get(nodeId)] = positionMapHorizontal.get(key);
            });
          } else {
            xCoords[nodeIndexes.get(key)] = positionMapHorizontal.get(key);
          }
        };

        var _iteratorNormalCompletion6 = true;
        var _didIteratorError6 = false;
        var _iteratorError6 = undefined;

        try {
          for (var _iterator6 = positionMapHorizontal.keys()[Symbol.iterator](), _step6; !(_iteratorNormalCompletion6 = (_step6 = _iterator6.next()).done); _iteratorNormalCompletion6 = true) {
            var key = _step6.value;

            _loop9(key);
          }
        } catch (err) {
          _didIteratorError6 = true;
          _iteratorError6 = err;
        } finally {
          try {
            if (!_iteratorNormalCompletion6 && _iterator6.return) {
              _iterator6.return();
            }
          } finally {
            if (_didIteratorError6) {
              throw _iteratorError6;
            }
          }
        }

        var _loop10 = function _loop10(key) {
          if (dummyToNodeForHorizontalAlignment.get(key)) {
            dummyToNodeForHorizontalAlignment.get(key).forEach(function (nodeId) {
              yCoords[nodeIndexes.get(nodeId)] = positionMapVertical.get(key);
            });
          } else {
            yCoords[nodeIndexes.get(key)] = positionMapVertical.get(key);
          }
        };

        var _iteratorNormalCompletion7 = true;
        var _didIteratorError7 = false;
        var _iteratorError7 = undefined;

        try {
          for (var _iterator7 = positionMapVertical.keys()[Symbol.iterator](), _step7; !(_iteratorNormalCompletion7 = (_step7 = _iterator7.next()).done); _iteratorNormalCompletion7 = true) {
            var key = _step7.value;

            _loop10(key);
          }
        } catch (err) {
          _didIteratorError7 = true;
          _iteratorError7 = err;
        } finally {
          try {
            if (!_iteratorNormalCompletion7 && _iterator7.return) {
              _iterator7.return();
            }
          } finally {
            if (_didIteratorError7) {
              throw _iteratorError7;
            }
          }
        }
      })();
    }
  }

  // assign new coordinates to nodes after constraint handling
  for (var _i13 = 0; _i13 < allNodes.length; _i13++) {
    var _node = allNodes[_i13];
    if (_node.getChild() == null) {
      _node.setCenter(xCoords[nodeIndexes.get(_node.id)], yCoords[nodeIndexes.get(_node.id)]);
    }
  }
};

module.exports = ConstraintHandler;

/***/ }),

/***/ 551:
/***/ ((module) => {

module.exports = __WEBPACK_EXTERNAL_MODULE__551__;

/***/ })

/******/ 	});
/************************************************************************/
/******/ 	// The module cache
/******/ 	var __webpack_module_cache__ = {};
/******/ 	
/******/ 	// The require function
/******/ 	function __nested_webpack_require_117795__(moduleId) {
/******/ 		// Check if module is in cache
/******/ 		var cachedModule = __webpack_module_cache__[moduleId];
/******/ 		if (cachedModule !== undefined) {
/******/ 			return cachedModule.exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = __webpack_module_cache__[moduleId] = {
/******/ 			// no module.id needed
/******/ 			// no module.loaded needed
/******/ 			exports: {}
/******/ 		};
/******/ 	
/******/ 		// Execute the module function
/******/ 		__webpack_modules__[moduleId](module, module.exports, __nested_webpack_require_117795__);
/******/ 	
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/ 	
/************************************************************************/
/******/ 	
/******/ 	// startup
/******/ 	// Load entry module and return exports
/******/ 	// This entry module is referenced by other modules so it can't be inlined
/******/ 	var __nested_webpack_exports__ = __nested_webpack_require_117795__(45);
/******/ 	
/******/ 	return __nested_webpack_exports__;
/******/ })()
;
});

/***/ }),

/***/ 12591:
/***/ (function(module) {

(function webpackUniversalModuleDefinition(root, factory) {
	if(true)
		module.exports = factory();
	else {}
})(this, function() {
return /******/ (function(modules) { // webpackBootstrap
/******/ 	// The module cache
/******/ 	var installedModules = {};
/******/
/******/ 	// The require function
/******/ 	function __nested_webpack_require_543__(moduleId) {
/******/
/******/ 		// Check if module is in cache
/******/ 		if(installedModules[moduleId]) {
/******/ 			return installedModules[moduleId].exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = installedModules[moduleId] = {
/******/ 			i: moduleId,
/******/ 			l: false,
/******/ 			exports: {}
/******/ 		};
/******/
/******/ 		// Execute the module function
/******/ 		modules[moduleId].call(module.exports, module, module.exports, __nested_webpack_require_543__);
/******/
/******/ 		// Flag the module as loaded
/******/ 		module.l = true;
/******/
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/
/******/
/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__nested_webpack_require_543__.m = modules;
/******/
/******/ 	// expose the module cache
/******/ 	__nested_webpack_require_543__.c = installedModules;
/******/
/******/ 	// identity function for calling harmony imports with the correct context
/******/ 	__nested_webpack_require_543__.i = function(value) { return value; };
/******/
/******/ 	// define getter function for harmony exports
/******/ 	__nested_webpack_require_543__.d = function(exports, name, getter) {
/******/ 		if(!__nested_webpack_require_543__.o(exports, name)) {
/******/ 			Object.defineProperty(exports, name, {
/******/ 				configurable: false,
/******/ 				enumerable: true,
/******/ 				get: getter
/******/ 			});
/******/ 		}
/******/ 	};
/******/
/******/ 	// getDefaultExport function for compatibility with non-harmony modules
/******/ 	__nested_webpack_require_543__.n = function(module) {
/******/ 		var getter = module && module.__esModule ?
/******/ 			function getDefault() { return module['default']; } :
/******/ 			function getModuleExports() { return module; };
/******/ 		__nested_webpack_require_543__.d(getter, 'a', getter);
/******/ 		return getter;
/******/ 	};
/******/
/******/ 	// Object.prototype.hasOwnProperty.call
/******/ 	__nested_webpack_require_543__.o = function(object, property) { return Object.prototype.hasOwnProperty.call(object, property); };
/******/
/******/ 	// __webpack_public_path__
/******/ 	__nested_webpack_require_543__.p = "";
/******/
/******/ 	// Load entry module and return exports
/******/ 	return __nested_webpack_require_543__(__nested_webpack_require_543__.s = 28);
/******/ })
/************************************************************************/
/******/ ([
/* 0 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


function LayoutConstants() {}

/**
 * Layout Quality: 0:draft, 1:default, 2:proof
 */
LayoutConstants.QUALITY = 1;

/**
 * Default parameters
 */
LayoutConstants.DEFAULT_CREATE_BENDS_AS_NEEDED = false;
LayoutConstants.DEFAULT_INCREMENTAL = false;
LayoutConstants.DEFAULT_ANIMATION_ON_LAYOUT = true;
LayoutConstants.DEFAULT_ANIMATION_DURING_LAYOUT = false;
LayoutConstants.DEFAULT_ANIMATION_PERIOD = 50;
LayoutConstants.DEFAULT_UNIFORM_LEAF_NODE_SIZES = false;

// -----------------------------------------------------------------------------
// Section: General other constants
// -----------------------------------------------------------------------------
/*
 * Margins of a graph to be applied on bouding rectangle of its contents. We
 * assume margins on all four sides to be uniform.
 */
LayoutConstants.DEFAULT_GRAPH_MARGIN = 15;

/*
 * Whether to consider labels in node dimensions or not
 */
LayoutConstants.NODE_DIMENSIONS_INCLUDE_LABELS = false;

/*
 * Default dimension of a non-compound node.
 */
LayoutConstants.SIMPLE_NODE_SIZE = 40;

/*
 * Default dimension of a non-compound node.
 */
LayoutConstants.SIMPLE_NODE_HALF_SIZE = LayoutConstants.SIMPLE_NODE_SIZE / 2;

/*
 * Empty compound node size. When a compound node is empty, its both
 * dimensions should be of this value.
 */
LayoutConstants.EMPTY_COMPOUND_NODE_SIZE = 40;

/*
 * Minimum length that an edge should take during layout
 */
LayoutConstants.MIN_EDGE_LENGTH = 1;

/*
 * World boundaries that layout operates on
 */
LayoutConstants.WORLD_BOUNDARY = 1000000;

/*
 * World boundaries that random positioning can be performed with
 */
LayoutConstants.INITIAL_WORLD_BOUNDARY = LayoutConstants.WORLD_BOUNDARY / 1000;

/*
 * Coordinates of the world center
 */
LayoutConstants.WORLD_CENTER_X = 1200;
LayoutConstants.WORLD_CENTER_Y = 900;

module.exports = LayoutConstants;

/***/ }),
/* 1 */
/***/ (function(module, exports, __nested_webpack_require_4947__) {

"use strict";


var LGraphObject = __nested_webpack_require_4947__(2);
var IGeometry = __nested_webpack_require_4947__(8);
var IMath = __nested_webpack_require_4947__(9);

function LEdge(source, target, vEdge) {
  LGraphObject.call(this, vEdge);

  this.isOverlapingSourceAndTarget = false;
  this.vGraphObject = vEdge;
  this.bendpoints = [];
  this.source = source;
  this.target = target;
}

LEdge.prototype = Object.create(LGraphObject.prototype);

for (var prop in LGraphObject) {
  LEdge[prop] = LGraphObject[prop];
}

LEdge.prototype.getSource = function () {
  return this.source;
};

LEdge.prototype.getTarget = function () {
  return this.target;
};

LEdge.prototype.isInterGraph = function () {
  return this.isInterGraph;
};

LEdge.prototype.getLength = function () {
  return this.length;
};

LEdge.prototype.isOverlapingSourceAndTarget = function () {
  return this.isOverlapingSourceAndTarget;
};

LEdge.prototype.getBendpoints = function () {
  return this.bendpoints;
};

LEdge.prototype.getLca = function () {
  return this.lca;
};

LEdge.prototype.getSourceInLca = function () {
  return this.sourceInLca;
};

LEdge.prototype.getTargetInLca = function () {
  return this.targetInLca;
};

LEdge.prototype.getOtherEnd = function (node) {
  if (this.source === node) {
    return this.target;
  } else if (this.target === node) {
    return this.source;
  } else {
    throw "Node is not incident with this edge";
  }
};

LEdge.prototype.getOtherEndInGraph = function (node, graph) {
  var otherEnd = this.getOtherEnd(node);
  var root = graph.getGraphManager().getRoot();

  while (true) {
    if (otherEnd.getOwner() == graph) {
      return otherEnd;
    }

    if (otherEnd.getOwner() == root) {
      break;
    }

    otherEnd = otherEnd.getOwner().getParent();
  }

  return null;
};

LEdge.prototype.updateLength = function () {
  var clipPointCoordinates = new Array(4);

  this.isOverlapingSourceAndTarget = IGeometry.getIntersection(this.target.getRect(), this.source.getRect(), clipPointCoordinates);

  if (!this.isOverlapingSourceAndTarget) {
    this.lengthX = clipPointCoordinates[0] - clipPointCoordinates[2];
    this.lengthY = clipPointCoordinates[1] - clipPointCoordinates[3];

    if (Math.abs(this.lengthX) < 1.0) {
      this.lengthX = IMath.sign(this.lengthX);
    }

    if (Math.abs(this.lengthY) < 1.0) {
      this.lengthY = IMath.sign(this.lengthY);
    }

    this.length = Math.sqrt(this.lengthX * this.lengthX + this.lengthY * this.lengthY);
  }
};

LEdge.prototype.updateLengthSimple = function () {
  this.lengthX = this.target.getCenterX() - this.source.getCenterX();
  this.lengthY = this.target.getCenterY() - this.source.getCenterY();

  if (Math.abs(this.lengthX) < 1.0) {
    this.lengthX = IMath.sign(this.lengthX);
  }

  if (Math.abs(this.lengthY) < 1.0) {
    this.lengthY = IMath.sign(this.lengthY);
  }

  this.length = Math.sqrt(this.lengthX * this.lengthX + this.lengthY * this.lengthY);
};

module.exports = LEdge;

/***/ }),
/* 2 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


function LGraphObject(vGraphObject) {
  this.vGraphObject = vGraphObject;
}

module.exports = LGraphObject;

/***/ }),
/* 3 */
/***/ (function(module, exports, __nested_webpack_require_8167__) {

"use strict";


var LGraphObject = __nested_webpack_require_8167__(2);
var Integer = __nested_webpack_require_8167__(10);
var RectangleD = __nested_webpack_require_8167__(13);
var LayoutConstants = __nested_webpack_require_8167__(0);
var RandomSeed = __nested_webpack_require_8167__(16);
var PointD = __nested_webpack_require_8167__(5);

function LNode(gm, loc, size, vNode) {
  //Alternative constructor 1 : LNode(LGraphManager gm, Point loc, Dimension size, Object vNode)
  if (size == null && vNode == null) {
    vNode = loc;
  }

  LGraphObject.call(this, vNode);

  //Alternative constructor 2 : LNode(Layout layout, Object vNode)
  if (gm.graphManager != null) gm = gm.graphManager;

  this.estimatedSize = Integer.MIN_VALUE;
  this.inclusionTreeDepth = Integer.MAX_VALUE;
  this.vGraphObject = vNode;
  this.edges = [];
  this.graphManager = gm;

  if (size != null && loc != null) this.rect = new RectangleD(loc.x, loc.y, size.width, size.height);else this.rect = new RectangleD();
}

LNode.prototype = Object.create(LGraphObject.prototype);
for (var prop in LGraphObject) {
  LNode[prop] = LGraphObject[prop];
}

LNode.prototype.getEdges = function () {
  return this.edges;
};

LNode.prototype.getChild = function () {
  return this.child;
};

LNode.prototype.getOwner = function () {
  //  if (this.owner != null) {
  //    if (!(this.owner == null || this.owner.getNodes().indexOf(this) > -1)) {
  //      throw "assert failed";
  //    }
  //  }

  return this.owner;
};

LNode.prototype.getWidth = function () {
  return this.rect.width;
};

LNode.prototype.setWidth = function (width) {
  this.rect.width = width;
};

LNode.prototype.getHeight = function () {
  return this.rect.height;
};

LNode.prototype.setHeight = function (height) {
  this.rect.height = height;
};

LNode.prototype.getCenterX = function () {
  return this.rect.x + this.rect.width / 2;
};

LNode.prototype.getCenterY = function () {
  return this.rect.y + this.rect.height / 2;
};

LNode.prototype.getCenter = function () {
  return new PointD(this.rect.x + this.rect.width / 2, this.rect.y + this.rect.height / 2);
};

LNode.prototype.getLocation = function () {
  return new PointD(this.rect.x, this.rect.y);
};

LNode.prototype.getRect = function () {
  return this.rect;
};

LNode.prototype.getDiagonal = function () {
  return Math.sqrt(this.rect.width * this.rect.width + this.rect.height * this.rect.height);
};

/**
 * This method returns half the diagonal length of this node.
 */
LNode.prototype.getHalfTheDiagonal = function () {
  return Math.sqrt(this.rect.height * this.rect.height + this.rect.width * this.rect.width) / 2;
};

LNode.prototype.setRect = function (upperLeft, dimension) {
  this.rect.x = upperLeft.x;
  this.rect.y = upperLeft.y;
  this.rect.width = dimension.width;
  this.rect.height = dimension.height;
};

LNode.prototype.setCenter = function (cx, cy) {
  this.rect.x = cx - this.rect.width / 2;
  this.rect.y = cy - this.rect.height / 2;
};

LNode.prototype.setLocation = function (x, y) {
  this.rect.x = x;
  this.rect.y = y;
};

LNode.prototype.moveBy = function (dx, dy) {
  this.rect.x += dx;
  this.rect.y += dy;
};

LNode.prototype.getEdgeListToNode = function (to) {
  var edgeList = [];
  var edge;
  var self = this;

  self.edges.forEach(function (edge) {

    if (edge.target == to) {
      if (edge.source != self) throw "Incorrect edge source!";

      edgeList.push(edge);
    }
  });

  return edgeList;
};

LNode.prototype.getEdgesBetween = function (other) {
  var edgeList = [];
  var edge;

  var self = this;
  self.edges.forEach(function (edge) {

    if (!(edge.source == self || edge.target == self)) throw "Incorrect edge source and/or target";

    if (edge.target == other || edge.source == other) {
      edgeList.push(edge);
    }
  });

  return edgeList;
};

LNode.prototype.getNeighborsList = function () {
  var neighbors = new Set();

  var self = this;
  self.edges.forEach(function (edge) {

    if (edge.source == self) {
      neighbors.add(edge.target);
    } else {
      if (edge.target != self) {
        throw "Incorrect incidency!";
      }

      neighbors.add(edge.source);
    }
  });

  return neighbors;
};

LNode.prototype.withChildren = function () {
  var withNeighborsList = new Set();
  var childNode;
  var children;

  withNeighborsList.add(this);

  if (this.child != null) {
    var nodes = this.child.getNodes();
    for (var i = 0; i < nodes.length; i++) {
      childNode = nodes[i];
      children = childNode.withChildren();
      children.forEach(function (node) {
        withNeighborsList.add(node);
      });
    }
  }

  return withNeighborsList;
};

LNode.prototype.getNoOfChildren = function () {
  var noOfChildren = 0;
  var childNode;

  if (this.child == null) {
    noOfChildren = 1;
  } else {
    var nodes = this.child.getNodes();
    for (var i = 0; i < nodes.length; i++) {
      childNode = nodes[i];

      noOfChildren += childNode.getNoOfChildren();
    }
  }

  if (noOfChildren == 0) {
    noOfChildren = 1;
  }
  return noOfChildren;
};

LNode.prototype.getEstimatedSize = function () {
  if (this.estimatedSize == Integer.MIN_VALUE) {
    throw "assert failed";
  }
  return this.estimatedSize;
};

LNode.prototype.calcEstimatedSize = function () {
  if (this.child == null) {
    return this.estimatedSize = (this.rect.width + this.rect.height) / 2;
  } else {
    this.estimatedSize = this.child.calcEstimatedSize();
    this.rect.width = this.estimatedSize;
    this.rect.height = this.estimatedSize;

    return this.estimatedSize;
  }
};

LNode.prototype.scatter = function () {
  var randomCenterX;
  var randomCenterY;

  var minX = -LayoutConstants.INITIAL_WORLD_BOUNDARY;
  var maxX = LayoutConstants.INITIAL_WORLD_BOUNDARY;
  randomCenterX = LayoutConstants.WORLD_CENTER_X + RandomSeed.nextDouble() * (maxX - minX) + minX;

  var minY = -LayoutConstants.INITIAL_WORLD_BOUNDARY;
  var maxY = LayoutConstants.INITIAL_WORLD_BOUNDARY;
  randomCenterY = LayoutConstants.WORLD_CENTER_Y + RandomSeed.nextDouble() * (maxY - minY) + minY;

  this.rect.x = randomCenterX;
  this.rect.y = randomCenterY;
};

LNode.prototype.updateBounds = function () {
  if (this.getChild() == null) {
    throw "assert failed";
  }
  if (this.getChild().getNodes().length != 0) {
    // wrap the children nodes by re-arranging the boundaries
    var childGraph = this.getChild();
    childGraph.updateBounds(true);

    this.rect.x = childGraph.getLeft();
    this.rect.y = childGraph.getTop();

    this.setWidth(childGraph.getRight() - childGraph.getLeft());
    this.setHeight(childGraph.getBottom() - childGraph.getTop());

    // Update compound bounds considering its label properties    
    if (LayoutConstants.NODE_DIMENSIONS_INCLUDE_LABELS) {

      var width = childGraph.getRight() - childGraph.getLeft();
      var height = childGraph.getBottom() - childGraph.getTop();

      if (this.labelWidth) {
        if (this.labelPosHorizontal == "left") {
          this.rect.x -= this.labelWidth;
          this.setWidth(width + this.labelWidth);
        } else if (this.labelPosHorizontal == "center" && this.labelWidth > width) {
          this.rect.x -= (this.labelWidth - width) / 2;
          this.setWidth(this.labelWidth);
        } else if (this.labelPosHorizontal == "right") {
          this.setWidth(width + this.labelWidth);
        }
      }

      if (this.labelHeight) {
        if (this.labelPosVertical == "top") {
          this.rect.y -= this.labelHeight;
          this.setHeight(height + this.labelHeight);
        } else if (this.labelPosVertical == "center" && this.labelHeight > height) {
          this.rect.y -= (this.labelHeight - height) / 2;
          this.setHeight(this.labelHeight);
        } else if (this.labelPosVertical == "bottom") {
          this.setHeight(height + this.labelHeight);
        }
      }
    }
  }
};

LNode.prototype.getInclusionTreeDepth = function () {
  if (this.inclusionTreeDepth == Integer.MAX_VALUE) {
    throw "assert failed";
  }
  return this.inclusionTreeDepth;
};

LNode.prototype.transform = function (trans) {
  var left = this.rect.x;

  if (left > LayoutConstants.WORLD_BOUNDARY) {
    left = LayoutConstants.WORLD_BOUNDARY;
  } else if (left < -LayoutConstants.WORLD_BOUNDARY) {
    left = -LayoutConstants.WORLD_BOUNDARY;
  }

  var top = this.rect.y;

  if (top > LayoutConstants.WORLD_BOUNDARY) {
    top = LayoutConstants.WORLD_BOUNDARY;
  } else if (top < -LayoutConstants.WORLD_BOUNDARY) {
    top = -LayoutConstants.WORLD_BOUNDARY;
  }

  var leftTop = new PointD(left, top);
  var vLeftTop = trans.inverseTransformPoint(leftTop);

  this.setLocation(vLeftTop.x, vLeftTop.y);
};

LNode.prototype.getLeft = function () {
  return this.rect.x;
};

LNode.prototype.getRight = function () {
  return this.rect.x + this.rect.width;
};

LNode.prototype.getTop = function () {
  return this.rect.y;
};

LNode.prototype.getBottom = function () {
  return this.rect.y + this.rect.height;
};

LNode.prototype.getParent = function () {
  if (this.owner == null) {
    return null;
  }

  return this.owner.getParent();
};

module.exports = LNode;

/***/ }),
/* 4 */
/***/ (function(module, exports, __nested_webpack_require_17290__) {

"use strict";


var LayoutConstants = __nested_webpack_require_17290__(0);

function FDLayoutConstants() {}

//FDLayoutConstants inherits static props in LayoutConstants
for (var prop in LayoutConstants) {
  FDLayoutConstants[prop] = LayoutConstants[prop];
}

FDLayoutConstants.MAX_ITERATIONS = 2500;

FDLayoutConstants.DEFAULT_EDGE_LENGTH = 50;
FDLayoutConstants.DEFAULT_SPRING_STRENGTH = 0.45;
FDLayoutConstants.DEFAULT_REPULSION_STRENGTH = 4500.0;
FDLayoutConstants.DEFAULT_GRAVITY_STRENGTH = 0.4;
FDLayoutConstants.DEFAULT_COMPOUND_GRAVITY_STRENGTH = 1.0;
FDLayoutConstants.DEFAULT_GRAVITY_RANGE_FACTOR = 3.8;
FDLayoutConstants.DEFAULT_COMPOUND_GRAVITY_RANGE_FACTOR = 1.5;
FDLayoutConstants.DEFAULT_USE_SMART_IDEAL_EDGE_LENGTH_CALCULATION = true;
FDLayoutConstants.DEFAULT_USE_SMART_REPULSION_RANGE_CALCULATION = true;
FDLayoutConstants.DEFAULT_COOLING_FACTOR_INCREMENTAL = 0.3;
FDLayoutConstants.COOLING_ADAPTATION_FACTOR = 0.33;
FDLayoutConstants.ADAPTATION_LOWER_NODE_LIMIT = 1000;
FDLayoutConstants.ADAPTATION_UPPER_NODE_LIMIT = 5000;
FDLayoutConstants.MAX_NODE_DISPLACEMENT_INCREMENTAL = 100.0;
FDLayoutConstants.MAX_NODE_DISPLACEMENT = FDLayoutConstants.MAX_NODE_DISPLACEMENT_INCREMENTAL * 3;
FDLayoutConstants.MIN_REPULSION_DIST = FDLayoutConstants.DEFAULT_EDGE_LENGTH / 10.0;
FDLayoutConstants.CONVERGENCE_CHECK_PERIOD = 100;
FDLayoutConstants.PER_LEVEL_IDEAL_EDGE_LENGTH_FACTOR = 0.1;
FDLayoutConstants.MIN_EDGE_LENGTH = 1;
FDLayoutConstants.GRID_CALCULATION_CHECK_PERIOD = 10;

module.exports = FDLayoutConstants;

/***/ }),
/* 5 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


function PointD(x, y) {
  if (x == null && y == null) {
    this.x = 0;
    this.y = 0;
  } else {
    this.x = x;
    this.y = y;
  }
}

PointD.prototype.getX = function () {
  return this.x;
};

PointD.prototype.getY = function () {
  return this.y;
};

PointD.prototype.setX = function (x) {
  this.x = x;
};

PointD.prototype.setY = function (y) {
  this.y = y;
};

PointD.prototype.getDifference = function (pt) {
  return new DimensionD(this.x - pt.x, this.y - pt.y);
};

PointD.prototype.getCopy = function () {
  return new PointD(this.x, this.y);
};

PointD.prototype.translate = function (dim) {
  this.x += dim.width;
  this.y += dim.height;
  return this;
};

module.exports = PointD;

/***/ }),
/* 6 */
/***/ (function(module, exports, __nested_webpack_require_19670__) {

"use strict";


var LGraphObject = __nested_webpack_require_19670__(2);
var Integer = __nested_webpack_require_19670__(10);
var LayoutConstants = __nested_webpack_require_19670__(0);
var LGraphManager = __nested_webpack_require_19670__(7);
var LNode = __nested_webpack_require_19670__(3);
var LEdge = __nested_webpack_require_19670__(1);
var RectangleD = __nested_webpack_require_19670__(13);
var Point = __nested_webpack_require_19670__(12);
var LinkedList = __nested_webpack_require_19670__(11);

function LGraph(parent, obj2, vGraph) {
  LGraphObject.call(this, vGraph);
  this.estimatedSize = Integer.MIN_VALUE;
  this.margin = LayoutConstants.DEFAULT_GRAPH_MARGIN;
  this.edges = [];
  this.nodes = [];
  this.isConnected = false;
  this.parent = parent;

  if (obj2 != null && obj2 instanceof LGraphManager) {
    this.graphManager = obj2;
  } else if (obj2 != null && obj2 instanceof Layout) {
    this.graphManager = obj2.graphManager;
  }
}

LGraph.prototype = Object.create(LGraphObject.prototype);
for (var prop in LGraphObject) {
  LGraph[prop] = LGraphObject[prop];
}

LGraph.prototype.getNodes = function () {
  return this.nodes;
};

LGraph.prototype.getEdges = function () {
  return this.edges;
};

LGraph.prototype.getGraphManager = function () {
  return this.graphManager;
};

LGraph.prototype.getParent = function () {
  return this.parent;
};

LGraph.prototype.getLeft = function () {
  return this.left;
};

LGraph.prototype.getRight = function () {
  return this.right;
};

LGraph.prototype.getTop = function () {
  return this.top;
};

LGraph.prototype.getBottom = function () {
  return this.bottom;
};

LGraph.prototype.isConnected = function () {
  return this.isConnected;
};

LGraph.prototype.add = function (obj1, sourceNode, targetNode) {
  if (sourceNode == null && targetNode == null) {
    var newNode = obj1;
    if (this.graphManager == null) {
      throw "Graph has no graph mgr!";
    }
    if (this.getNodes().indexOf(newNode) > -1) {
      throw "Node already in graph!";
    }
    newNode.owner = this;
    this.getNodes().push(newNode);

    return newNode;
  } else {
    var newEdge = obj1;
    if (!(this.getNodes().indexOf(sourceNode) > -1 && this.getNodes().indexOf(targetNode) > -1)) {
      throw "Source or target not in graph!";
    }

    if (!(sourceNode.owner == targetNode.owner && sourceNode.owner == this)) {
      throw "Both owners must be this graph!";
    }

    if (sourceNode.owner != targetNode.owner) {
      return null;
    }

    // set source and target
    newEdge.source = sourceNode;
    newEdge.target = targetNode;

    // set as intra-graph edge
    newEdge.isInterGraph = false;

    // add to graph edge list
    this.getEdges().push(newEdge);

    // add to incidency lists
    sourceNode.edges.push(newEdge);

    if (targetNode != sourceNode) {
      targetNode.edges.push(newEdge);
    }

    return newEdge;
  }
};

LGraph.prototype.remove = function (obj) {
  var node = obj;
  if (obj instanceof LNode) {
    if (node == null) {
      throw "Node is null!";
    }
    if (!(node.owner != null && node.owner == this)) {
      throw "Owner graph is invalid!";
    }
    if (this.graphManager == null) {
      throw "Owner graph manager is invalid!";
    }
    // remove incident edges first (make a copy to do it safely)
    var edgesToBeRemoved = node.edges.slice();
    var edge;
    var s = edgesToBeRemoved.length;
    for (var i = 0; i < s; i++) {
      edge = edgesToBeRemoved[i];

      if (edge.isInterGraph) {
        this.graphManager.remove(edge);
      } else {
        edge.source.owner.remove(edge);
      }
    }

    // now the node itself
    var index = this.nodes.indexOf(node);
    if (index == -1) {
      throw "Node not in owner node list!";
    }

    this.nodes.splice(index, 1);
  } else if (obj instanceof LEdge) {
    var edge = obj;
    if (edge == null) {
      throw "Edge is null!";
    }
    if (!(edge.source != null && edge.target != null)) {
      throw "Source and/or target is null!";
    }
    if (!(edge.source.owner != null && edge.target.owner != null && edge.source.owner == this && edge.target.owner == this)) {
      throw "Source and/or target owner is invalid!";
    }

    var sourceIndex = edge.source.edges.indexOf(edge);
    var targetIndex = edge.target.edges.indexOf(edge);
    if (!(sourceIndex > -1 && targetIndex > -1)) {
      throw "Source and/or target doesn't know this edge!";
    }

    edge.source.edges.splice(sourceIndex, 1);

    if (edge.target != edge.source) {
      edge.target.edges.splice(targetIndex, 1);
    }

    var index = edge.source.owner.getEdges().indexOf(edge);
    if (index == -1) {
      throw "Not in owner's edge list!";
    }

    edge.source.owner.getEdges().splice(index, 1);
  }
};

LGraph.prototype.updateLeftTop = function () {
  var top = Integer.MAX_VALUE;
  var left = Integer.MAX_VALUE;
  var nodeTop;
  var nodeLeft;
  var margin;

  var nodes = this.getNodes();
  var s = nodes.length;

  for (var i = 0; i < s; i++) {
    var lNode = nodes[i];
    nodeTop = lNode.getTop();
    nodeLeft = lNode.getLeft();

    if (top > nodeTop) {
      top = nodeTop;
    }

    if (left > nodeLeft) {
      left = nodeLeft;
    }
  }

  // Do we have any nodes in this graph?
  if (top == Integer.MAX_VALUE) {
    return null;
  }

  if (nodes[0].getParent().paddingLeft != undefined) {
    margin = nodes[0].getParent().paddingLeft;
  } else {
    margin = this.margin;
  }

  this.left = left - margin;
  this.top = top - margin;

  // Apply the margins and return the result
  return new Point(this.left, this.top);
};

LGraph.prototype.updateBounds = function (recursive) {
  // calculate bounds
  var left = Integer.MAX_VALUE;
  var right = -Integer.MAX_VALUE;
  var top = Integer.MAX_VALUE;
  var bottom = -Integer.MAX_VALUE;
  var nodeLeft;
  var nodeRight;
  var nodeTop;
  var nodeBottom;
  var margin;

  var nodes = this.nodes;
  var s = nodes.length;
  for (var i = 0; i < s; i++) {
    var lNode = nodes[i];

    if (recursive && lNode.child != null) {
      lNode.updateBounds();
    }
    nodeLeft = lNode.getLeft();
    nodeRight = lNode.getRight();
    nodeTop = lNode.getTop();
    nodeBottom = lNode.getBottom();

    if (left > nodeLeft) {
      left = nodeLeft;
    }

    if (right < nodeRight) {
      right = nodeRight;
    }

    if (top > nodeTop) {
      top = nodeTop;
    }

    if (bottom < nodeBottom) {
      bottom = nodeBottom;
    }
  }

  var boundingRect = new RectangleD(left, top, right - left, bottom - top);
  if (left == Integer.MAX_VALUE) {
    this.left = this.parent.getLeft();
    this.right = this.parent.getRight();
    this.top = this.parent.getTop();
    this.bottom = this.parent.getBottom();
  }

  if (nodes[0].getParent().paddingLeft != undefined) {
    margin = nodes[0].getParent().paddingLeft;
  } else {
    margin = this.margin;
  }

  this.left = boundingRect.x - margin;
  this.right = boundingRect.x + boundingRect.width + margin;
  this.top = boundingRect.y - margin;
  this.bottom = boundingRect.y + boundingRect.height + margin;
};

LGraph.calculateBounds = function (nodes) {
  var left = Integer.MAX_VALUE;
  var right = -Integer.MAX_VALUE;
  var top = Integer.MAX_VALUE;
  var bottom = -Integer.MAX_VALUE;
  var nodeLeft;
  var nodeRight;
  var nodeTop;
  var nodeBottom;

  var s = nodes.length;

  for (var i = 0; i < s; i++) {
    var lNode = nodes[i];
    nodeLeft = lNode.getLeft();
    nodeRight = lNode.getRight();
    nodeTop = lNode.getTop();
    nodeBottom = lNode.getBottom();

    if (left > nodeLeft) {
      left = nodeLeft;
    }

    if (right < nodeRight) {
      right = nodeRight;
    }

    if (top > nodeTop) {
      top = nodeTop;
    }

    if (bottom < nodeBottom) {
      bottom = nodeBottom;
    }
  }

  var boundingRect = new RectangleD(left, top, right - left, bottom - top);

  return boundingRect;
};

LGraph.prototype.getInclusionTreeDepth = function () {
  if (this == this.graphManager.getRoot()) {
    return 1;
  } else {
    return this.parent.getInclusionTreeDepth();
  }
};

LGraph.prototype.getEstimatedSize = function () {
  if (this.estimatedSize == Integer.MIN_VALUE) {
    throw "assert failed";
  }
  return this.estimatedSize;
};

LGraph.prototype.calcEstimatedSize = function () {
  var size = 0;
  var nodes = this.nodes;
  var s = nodes.length;

  for (var i = 0; i < s; i++) {
    var lNode = nodes[i];
    size += lNode.calcEstimatedSize();
  }

  if (size == 0) {
    this.estimatedSize = LayoutConstants.EMPTY_COMPOUND_NODE_SIZE;
  } else {
    this.estimatedSize = size / Math.sqrt(this.nodes.length);
  }

  return this.estimatedSize;
};

LGraph.prototype.updateConnected = function () {
  var self = this;
  if (this.nodes.length == 0) {
    this.isConnected = true;
    return;
  }

  var queue = new LinkedList();
  var visited = new Set();
  var currentNode = this.nodes[0];
  var neighborEdges;
  var currentNeighbor;
  var childrenOfNode = currentNode.withChildren();
  childrenOfNode.forEach(function (node) {
    queue.push(node);
    visited.add(node);
  });

  while (queue.length !== 0) {
    currentNode = queue.shift();

    // Traverse all neighbors of this node
    neighborEdges = currentNode.getEdges();
    var size = neighborEdges.length;
    for (var i = 0; i < size; i++) {
      var neighborEdge = neighborEdges[i];
      currentNeighbor = neighborEdge.getOtherEndInGraph(currentNode, this);

      // Add unvisited neighbors to the list to visit
      if (currentNeighbor != null && !visited.has(currentNeighbor)) {
        var childrenOfNeighbor = currentNeighbor.withChildren();

        childrenOfNeighbor.forEach(function (node) {
          queue.push(node);
          visited.add(node);
        });
      }
    }
  }

  this.isConnected = false;

  if (visited.size >= this.nodes.length) {
    var noOfVisitedInThisGraph = 0;

    visited.forEach(function (visitedNode) {
      if (visitedNode.owner == self) {
        noOfVisitedInThisGraph++;
      }
    });

    if (noOfVisitedInThisGraph == this.nodes.length) {
      this.isConnected = true;
    }
  }
};

module.exports = LGraph;

/***/ }),
/* 7 */
/***/ (function(module, exports, __nested_webpack_require_29738__) {

"use strict";


var LGraph;
var LEdge = __nested_webpack_require_29738__(1);

function LGraphManager(layout) {
  LGraph = __nested_webpack_require_29738__(6); // It may be better to initilize this out of this function but it gives an error (Right-hand side of 'instanceof' is not callable) now.
  this.layout = layout;

  this.graphs = [];
  this.edges = [];
}

LGraphManager.prototype.addRoot = function () {
  var ngraph = this.layout.newGraph();
  var nnode = this.layout.newNode(null);
  var root = this.add(ngraph, nnode);
  this.setRootGraph(root);
  return this.rootGraph;
};

LGraphManager.prototype.add = function (newGraph, parentNode, newEdge, sourceNode, targetNode) {
  //there are just 2 parameters are passed then it adds an LGraph else it adds an LEdge
  if (newEdge == null && sourceNode == null && targetNode == null) {
    if (newGraph == null) {
      throw "Graph is null!";
    }
    if (parentNode == null) {
      throw "Parent node is null!";
    }
    if (this.graphs.indexOf(newGraph) > -1) {
      throw "Graph already in this graph mgr!";
    }

    this.graphs.push(newGraph);

    if (newGraph.parent != null) {
      throw "Already has a parent!";
    }
    if (parentNode.child != null) {
      throw "Already has a child!";
    }

    newGraph.parent = parentNode;
    parentNode.child = newGraph;

    return newGraph;
  } else {
    //change the order of the parameters
    targetNode = newEdge;
    sourceNode = parentNode;
    newEdge = newGraph;
    var sourceGraph = sourceNode.getOwner();
    var targetGraph = targetNode.getOwner();

    if (!(sourceGraph != null && sourceGraph.getGraphManager() == this)) {
      throw "Source not in this graph mgr!";
    }
    if (!(targetGraph != null && targetGraph.getGraphManager() == this)) {
      throw "Target not in this graph mgr!";
    }

    if (sourceGraph == targetGraph) {
      newEdge.isInterGraph = false;
      return sourceGraph.add(newEdge, sourceNode, targetNode);
    } else {
      newEdge.isInterGraph = true;

      // set source and target
      newEdge.source = sourceNode;
      newEdge.target = targetNode;

      // add edge to inter-graph edge list
      if (this.edges.indexOf(newEdge) > -1) {
        throw "Edge already in inter-graph edge list!";
      }

      this.edges.push(newEdge);

      // add edge to source and target incidency lists
      if (!(newEdge.source != null && newEdge.target != null)) {
        throw "Edge source and/or target is null!";
      }

      if (!(newEdge.source.edges.indexOf(newEdge) == -1 && newEdge.target.edges.indexOf(newEdge) == -1)) {
        throw "Edge already in source and/or target incidency list!";
      }

      newEdge.source.edges.push(newEdge);
      newEdge.target.edges.push(newEdge);

      return newEdge;
    }
  }
};

LGraphManager.prototype.remove = function (lObj) {
  if (lObj instanceof LGraph) {
    var graph = lObj;
    if (graph.getGraphManager() != this) {
      throw "Graph not in this graph mgr";
    }
    if (!(graph == this.rootGraph || graph.parent != null && graph.parent.graphManager == this)) {
      throw "Invalid parent node!";
    }

    // first the edges (make a copy to do it safely)
    var edgesToBeRemoved = [];

    edgesToBeRemoved = edgesToBeRemoved.concat(graph.getEdges());

    var edge;
    var s = edgesToBeRemoved.length;
    for (var i = 0; i < s; i++) {
      edge = edgesToBeRemoved[i];
      graph.remove(edge);
    }

    // then the nodes (make a copy to do it safely)
    var nodesToBeRemoved = [];

    nodesToBeRemoved = nodesToBeRemoved.concat(graph.getNodes());

    var node;
    s = nodesToBeRemoved.length;
    for (var i = 0; i < s; i++) {
      node = nodesToBeRemoved[i];
      graph.remove(node);
    }

    // check if graph is the root
    if (graph == this.rootGraph) {
      this.setRootGraph(null);
    }

    // now remove the graph itself
    var index = this.graphs.indexOf(graph);
    this.graphs.splice(index, 1);

    // also reset the parent of the graph
    graph.parent = null;
  } else if (lObj instanceof LEdge) {
    edge = lObj;
    if (edge == null) {
      throw "Edge is null!";
    }
    if (!edge.isInterGraph) {
      throw "Not an inter-graph edge!";
    }
    if (!(edge.source != null && edge.target != null)) {
      throw "Source and/or target is null!";
    }

    // remove edge from source and target nodes' incidency lists

    if (!(edge.source.edges.indexOf(edge) != -1 && edge.target.edges.indexOf(edge) != -1)) {
      throw "Source and/or target doesn't know this edge!";
    }

    var index = edge.source.edges.indexOf(edge);
    edge.source.edges.splice(index, 1);
    index = edge.target.edges.indexOf(edge);
    edge.target.edges.splice(index, 1);

    // remove edge from owner graph manager's inter-graph edge list

    if (!(edge.source.owner != null && edge.source.owner.getGraphManager() != null)) {
      throw "Edge owner graph or owner graph manager is null!";
    }
    if (edge.source.owner.getGraphManager().edges.indexOf(edge) == -1) {
      throw "Not in owner graph manager's edge list!";
    }

    var index = edge.source.owner.getGraphManager().edges.indexOf(edge);
    edge.source.owner.getGraphManager().edges.splice(index, 1);
  }
};

LGraphManager.prototype.updateBounds = function () {
  this.rootGraph.updateBounds(true);
};

LGraphManager.prototype.getGraphs = function () {
  return this.graphs;
};

LGraphManager.prototype.getAllNodes = function () {
  if (this.allNodes == null) {
    var nodeList = [];
    var graphs = this.getGraphs();
    var s = graphs.length;
    for (var i = 0; i < s; i++) {
      nodeList = nodeList.concat(graphs[i].getNodes());
    }
    this.allNodes = nodeList;
  }
  return this.allNodes;
};

LGraphManager.prototype.resetAllNodes = function () {
  this.allNodes = null;
};

LGraphManager.prototype.resetAllEdges = function () {
  this.allEdges = null;
};

LGraphManager.prototype.resetAllNodesToApplyGravitation = function () {
  this.allNodesToApplyGravitation = null;
};

LGraphManager.prototype.getAllEdges = function () {
  if (this.allEdges == null) {
    var edgeList = [];
    var graphs = this.getGraphs();
    var s = graphs.length;
    for (var i = 0; i < graphs.length; i++) {
      edgeList = edgeList.concat(graphs[i].getEdges());
    }

    edgeList = edgeList.concat(this.edges);

    this.allEdges = edgeList;
  }
  return this.allEdges;
};

LGraphManager.prototype.getAllNodesToApplyGravitation = function () {
  return this.allNodesToApplyGravitation;
};

LGraphManager.prototype.setAllNodesToApplyGravitation = function (nodeList) {
  if (this.allNodesToApplyGravitation != null) {
    throw "assert failed";
  }

  this.allNodesToApplyGravitation = nodeList;
};

LGraphManager.prototype.getRoot = function () {
  return this.rootGraph;
};

LGraphManager.prototype.setRootGraph = function (graph) {
  if (graph.getGraphManager() != this) {
    throw "Root not in this graph mgr!";
  }

  this.rootGraph = graph;
  // root graph must have a root node associated with it for convenience
  if (graph.parent == null) {
    graph.parent = this.layout.newNode("Root node");
  }
};

LGraphManager.prototype.getLayout = function () {
  return this.layout;
};

LGraphManager.prototype.isOneAncestorOfOther = function (firstNode, secondNode) {
  if (!(firstNode != null && secondNode != null)) {
    throw "assert failed";
  }

  if (firstNode == secondNode) {
    return true;
  }
  // Is second node an ancestor of the first one?
  var ownerGraph = firstNode.getOwner();
  var parentNode;

  do {
    parentNode = ownerGraph.getParent();

    if (parentNode == null) {
      break;
    }

    if (parentNode == secondNode) {
      return true;
    }

    ownerGraph = parentNode.getOwner();
    if (ownerGraph == null) {
      break;
    }
  } while (true);
  // Is first node an ancestor of the second one?
  ownerGraph = secondNode.getOwner();

  do {
    parentNode = ownerGraph.getParent();

    if (parentNode == null) {
      break;
    }

    if (parentNode == firstNode) {
      return true;
    }

    ownerGraph = parentNode.getOwner();
    if (ownerGraph == null) {
      break;
    }
  } while (true);

  return false;
};

LGraphManager.prototype.calcLowestCommonAncestors = function () {
  var edge;
  var sourceNode;
  var targetNode;
  var sourceAncestorGraph;
  var targetAncestorGraph;

  var edges = this.getAllEdges();
  var s = edges.length;
  for (var i = 0; i < s; i++) {
    edge = edges[i];

    sourceNode = edge.source;
    targetNode = edge.target;
    edge.lca = null;
    edge.sourceInLca = sourceNode;
    edge.targetInLca = targetNode;

    if (sourceNode == targetNode) {
      edge.lca = sourceNode.getOwner();
      continue;
    }

    sourceAncestorGraph = sourceNode.getOwner();

    while (edge.lca == null) {
      edge.targetInLca = targetNode;
      targetAncestorGraph = targetNode.getOwner();

      while (edge.lca == null) {
        if (targetAncestorGraph == sourceAncestorGraph) {
          edge.lca = targetAncestorGraph;
          break;
        }

        if (targetAncestorGraph == this.rootGraph) {
          break;
        }

        if (edge.lca != null) {
          throw "assert failed";
        }
        edge.targetInLca = targetAncestorGraph.getParent();
        targetAncestorGraph = edge.targetInLca.getOwner();
      }

      if (sourceAncestorGraph == this.rootGraph) {
        break;
      }

      if (edge.lca == null) {
        edge.sourceInLca = sourceAncestorGraph.getParent();
        sourceAncestorGraph = edge.sourceInLca.getOwner();
      }
    }

    if (edge.lca == null) {
      throw "assert failed";
    }
  }
};

LGraphManager.prototype.calcLowestCommonAncestor = function (firstNode, secondNode) {
  if (firstNode == secondNode) {
    return firstNode.getOwner();
  }
  var firstOwnerGraph = firstNode.getOwner();

  do {
    if (firstOwnerGraph == null) {
      break;
    }
    var secondOwnerGraph = secondNode.getOwner();

    do {
      if (secondOwnerGraph == null) {
        break;
      }

      if (secondOwnerGraph == firstOwnerGraph) {
        return secondOwnerGraph;
      }
      secondOwnerGraph = secondOwnerGraph.getParent().getOwner();
    } while (true);

    firstOwnerGraph = firstOwnerGraph.getParent().getOwner();
  } while (true);

  return firstOwnerGraph;
};

LGraphManager.prototype.calcInclusionTreeDepths = function (graph, depth) {
  if (graph == null && depth == null) {
    graph = this.rootGraph;
    depth = 1;
  }
  var node;

  var nodes = graph.getNodes();
  var s = nodes.length;
  for (var i = 0; i < s; i++) {
    node = nodes[i];
    node.inclusionTreeDepth = depth;

    if (node.child != null) {
      this.calcInclusionTreeDepths(node.child, depth + 1);
    }
  }
};

LGraphManager.prototype.includesInvalidEdge = function () {
  var edge;
  var edgesToRemove = [];

  var s = this.edges.length;
  for (var i = 0; i < s; i++) {
    edge = this.edges[i];

    if (this.isOneAncestorOfOther(edge.source, edge.target)) {
      edgesToRemove.push(edge);
    }
  }

  // Remove invalid edges from graph manager
  for (var i = 0; i < edgesToRemove.length; i++) {
    this.remove(edgesToRemove[i]);
  }

  // Invalid edges are cleared, so return false
  return false;
};

module.exports = LGraphManager;

/***/ }),
/* 8 */
/***/ (function(module, exports, __nested_webpack_require_41052__) {

"use strict";


/**
 * This class maintains a list of static geometry related utility methods.
 *
 *
 * Copyright: i-Vis Research Group, Bilkent University, 2007 - present
 */

var Point = __nested_webpack_require_41052__(12);

function IGeometry() {}

/**
 * This method calculates *half* the amount in x and y directions of the two
 * input rectangles needed to separate them keeping their respective
 * positioning, and returns the result in the input array. An input
 * separation buffer added to the amount in both directions. We assume that
 * the two rectangles do intersect.
 */
IGeometry.calcSeparationAmount = function (rectA, rectB, overlapAmount, separationBuffer) {
  if (!rectA.intersects(rectB)) {
    throw "assert failed";
  }

  var directions = new Array(2);

  this.decideDirectionsForOverlappingNodes(rectA, rectB, directions);

  overlapAmount[0] = Math.min(rectA.getRight(), rectB.getRight()) - Math.max(rectA.x, rectB.x);
  overlapAmount[1] = Math.min(rectA.getBottom(), rectB.getBottom()) - Math.max(rectA.y, rectB.y);

  // update the overlapping amounts for the following cases:
  if (rectA.getX() <= rectB.getX() && rectA.getRight() >= rectB.getRight()) {
    /* Case x.1:
    *
    * rectA
    * 	|                       |
    * 	|        _________      |
    * 	|        |       |      |
    * 	|________|_______|______|
    * 			 |       |
    *           |       |
    *        rectB
    */
    overlapAmount[0] += Math.min(rectB.getX() - rectA.getX(), rectA.getRight() - rectB.getRight());
  } else if (rectB.getX() <= rectA.getX() && rectB.getRight() >= rectA.getRight()) {
    /* Case x.2:
    *
    * rectB
    * 	|                       |
    * 	|        _________      |
    * 	|        |       |      |
    * 	|________|_______|______|
    * 			 |       |
    *           |       |
    *        rectA
    */
    overlapAmount[0] += Math.min(rectA.getX() - rectB.getX(), rectB.getRight() - rectA.getRight());
  }
  if (rectA.getY() <= rectB.getY() && rectA.getBottom() >= rectB.getBottom()) {
    /* Case y.1:
     *          ________ rectA
     *         |
     *         |
     *   ______|____  rectB
     *         |    |
     *         |    |
     *   ______|____|
     *         |
     *         |
     *         |________
     *
     */
    overlapAmount[1] += Math.min(rectB.getY() - rectA.getY(), rectA.getBottom() - rectB.getBottom());
  } else if (rectB.getY() <= rectA.getY() && rectB.getBottom() >= rectA.getBottom()) {
    /* Case y.2:
    *          ________ rectB
    *         |
    *         |
    *   ______|____  rectA
    *         |    |
    *         |    |
    *   ______|____|
    *         |
    *         |
    *         |________
    *
    */
    overlapAmount[1] += Math.min(rectA.getY() - rectB.getY(), rectB.getBottom() - rectA.getBottom());
  }

  // find slope of the line passes two centers
  var slope = Math.abs((rectB.getCenterY() - rectA.getCenterY()) / (rectB.getCenterX() - rectA.getCenterX()));
  // if centers are overlapped
  if (rectB.getCenterY() === rectA.getCenterY() && rectB.getCenterX() === rectA.getCenterX()) {
    // assume the slope is 1 (45 degree)
    slope = 1.0;
  }

  var moveByY = slope * overlapAmount[0];
  var moveByX = overlapAmount[1] / slope;
  if (overlapAmount[0] < moveByX) {
    moveByX = overlapAmount[0];
  } else {
    moveByY = overlapAmount[1];
  }
  // return half the amount so that if each rectangle is moved by these
  // amounts in opposite directions, overlap will be resolved
  overlapAmount[0] = -1 * directions[0] * (moveByX / 2 + separationBuffer);
  overlapAmount[1] = -1 * directions[1] * (moveByY / 2 + separationBuffer);
};

/**
 * This method decides the separation direction of overlapping nodes
 *
 * if directions[0] = -1, then rectA goes left
 * if directions[0] = 1,  then rectA goes right
 * if directions[1] = -1, then rectA goes up
 * if directions[1] = 1,  then rectA goes down
 */
IGeometry.decideDirectionsForOverlappingNodes = function (rectA, rectB, directions) {
  if (rectA.getCenterX() < rectB.getCenterX()) {
    directions[0] = -1;
  } else {
    directions[0] = 1;
  }

  if (rectA.getCenterY() < rectB.getCenterY()) {
    directions[1] = -1;
  } else {
    directions[1] = 1;
  }
};

/**
 * This method calculates the intersection (clipping) points of the two
 * input rectangles with line segment defined by the centers of these two
 * rectangles. The clipping points are saved in the input double array and
 * whether or not the two rectangles overlap is returned.
 */
IGeometry.getIntersection2 = function (rectA, rectB, result) {
  //result[0-1] will contain clipPoint of rectA, result[2-3] will contain clipPoint of rectB
  var p1x = rectA.getCenterX();
  var p1y = rectA.getCenterY();
  var p2x = rectB.getCenterX();
  var p2y = rectB.getCenterY();

  //if two rectangles intersect, then clipping points are centers
  if (rectA.intersects(rectB)) {
    result[0] = p1x;
    result[1] = p1y;
    result[2] = p2x;
    result[3] = p2y;
    return true;
  }
  //variables for rectA
  var topLeftAx = rectA.getX();
  var topLeftAy = rectA.getY();
  var topRightAx = rectA.getRight();
  var bottomLeftAx = rectA.getX();
  var bottomLeftAy = rectA.getBottom();
  var bottomRightAx = rectA.getRight();
  var halfWidthA = rectA.getWidthHalf();
  var halfHeightA = rectA.getHeightHalf();
  //variables for rectB
  var topLeftBx = rectB.getX();
  var topLeftBy = rectB.getY();
  var topRightBx = rectB.getRight();
  var bottomLeftBx = rectB.getX();
  var bottomLeftBy = rectB.getBottom();
  var bottomRightBx = rectB.getRight();
  var halfWidthB = rectB.getWidthHalf();
  var halfHeightB = rectB.getHeightHalf();

  //flag whether clipping points are found
  var clipPointAFound = false;
  var clipPointBFound = false;

  // line is vertical
  if (p1x === p2x) {
    if (p1y > p2y) {
      result[0] = p1x;
      result[1] = topLeftAy;
      result[2] = p2x;
      result[3] = bottomLeftBy;
      return false;
    } else if (p1y < p2y) {
      result[0] = p1x;
      result[1] = bottomLeftAy;
      result[2] = p2x;
      result[3] = topLeftBy;
      return false;
    } else {
      //not line, return null;
    }
  }
  // line is horizontal
  else if (p1y === p2y) {
      if (p1x > p2x) {
        result[0] = topLeftAx;
        result[1] = p1y;
        result[2] = topRightBx;
        result[3] = p2y;
        return false;
      } else if (p1x < p2x) {
        result[0] = topRightAx;
        result[1] = p1y;
        result[2] = topLeftBx;
        result[3] = p2y;
        return false;
      } else {
        //not valid line, return null;
      }
    } else {
      //slopes of rectA's and rectB's diagonals
      var slopeA = rectA.height / rectA.width;
      var slopeB = rectB.height / rectB.width;

      //slope of line between center of rectA and center of rectB
      var slopePrime = (p2y - p1y) / (p2x - p1x);
      var cardinalDirectionA = void 0;
      var cardinalDirectionB = void 0;
      var tempPointAx = void 0;
      var tempPointAy = void 0;
      var tempPointBx = void 0;
      var tempPointBy = void 0;

      //determine whether clipping point is the corner of nodeA
      if (-slopeA === slopePrime) {
        if (p1x > p2x) {
          result[0] = bottomLeftAx;
          result[1] = bottomLeftAy;
          clipPointAFound = true;
        } else {
          result[0] = topRightAx;
          result[1] = topLeftAy;
          clipPointAFound = true;
        }
      } else if (slopeA === slopePrime) {
        if (p1x > p2x) {
          result[0] = topLeftAx;
          result[1] = topLeftAy;
          clipPointAFound = true;
        } else {
          result[0] = bottomRightAx;
          result[1] = bottomLeftAy;
          clipPointAFound = true;
        }
      }

      //determine whether clipping point is the corner of nodeB
      if (-slopeB === slopePrime) {
        if (p2x > p1x) {
          result[2] = bottomLeftBx;
          result[3] = bottomLeftBy;
          clipPointBFound = true;
        } else {
          result[2] = topRightBx;
          result[3] = topLeftBy;
          clipPointBFound = true;
        }
      } else if (slopeB === slopePrime) {
        if (p2x > p1x) {
          result[2] = topLeftBx;
          result[3] = topLeftBy;
          clipPointBFound = true;
        } else {
          result[2] = bottomRightBx;
          result[3] = bottomLeftBy;
          clipPointBFound = true;
        }
      }

      //if both clipping points are corners
      if (clipPointAFound && clipPointBFound) {
        return false;
      }

      //determine Cardinal Direction of rectangles
      if (p1x > p2x) {
        if (p1y > p2y) {
          cardinalDirectionA = this.getCardinalDirection(slopeA, slopePrime, 4);
          cardinalDirectionB = this.getCardinalDirection(slopeB, slopePrime, 2);
        } else {
          cardinalDirectionA = this.getCardinalDirection(-slopeA, slopePrime, 3);
          cardinalDirectionB = this.getCardinalDirection(-slopeB, slopePrime, 1);
        }
      } else {
        if (p1y > p2y) {
          cardinalDirectionA = this.getCardinalDirection(-slopeA, slopePrime, 1);
          cardinalDirectionB = this.getCardinalDirection(-slopeB, slopePrime, 3);
        } else {
          cardinalDirectionA = this.getCardinalDirection(slopeA, slopePrime, 2);
          cardinalDirectionB = this.getCardinalDirection(slopeB, slopePrime, 4);
        }
      }
      //calculate clipping Point if it is not found before
      if (!clipPointAFound) {
        switch (cardinalDirectionA) {
          case 1:
            tempPointAy = topLeftAy;
            tempPointAx = p1x + -halfHeightA / slopePrime;
            result[0] = tempPointAx;
            result[1] = tempPointAy;
            break;
          case 2:
            tempPointAx = bottomRightAx;
            tempPointAy = p1y + halfWidthA * slopePrime;
            result[0] = tempPointAx;
            result[1] = tempPointAy;
            break;
          case 3:
            tempPointAy = bottomLeftAy;
            tempPointAx = p1x + halfHeightA / slopePrime;
            result[0] = tempPointAx;
            result[1] = tempPointAy;
            break;
          case 4:
            tempPointAx = bottomLeftAx;
            tempPointAy = p1y + -halfWidthA * slopePrime;
            result[0] = tempPointAx;
            result[1] = tempPointAy;
            break;
        }
      }
      if (!clipPointBFound) {
        switch (cardinalDirectionB) {
          case 1:
            tempPointBy = topLeftBy;
            tempPointBx = p2x + -halfHeightB / slopePrime;
            result[2] = tempPointBx;
            result[3] = tempPointBy;
            break;
          case 2:
            tempPointBx = bottomRightBx;
            tempPointBy = p2y + halfWidthB * slopePrime;
            result[2] = tempPointBx;
            result[3] = tempPointBy;
            break;
          case 3:
            tempPointBy = bottomLeftBy;
            tempPointBx = p2x + halfHeightB / slopePrime;
            result[2] = tempPointBx;
            result[3] = tempPointBy;
            break;
          case 4:
            tempPointBx = bottomLeftBx;
            tempPointBy = p2y + -halfWidthB * slopePrime;
            result[2] = tempPointBx;
            result[3] = tempPointBy;
            break;
        }
      }
    }
  return false;
};

/**
 * This method returns in which cardinal direction does input point stays
 * 1: North
 * 2: East
 * 3: South
 * 4: West
 */
IGeometry.getCardinalDirection = function (slope, slopePrime, line) {
  if (slope > slopePrime) {
    return line;
  } else {
    return 1 + line % 4;
  }
};

/**
 * This method calculates the intersection of the two lines defined by
 * point pairs (s1,s2) and (f1,f2).
 */
IGeometry.getIntersection = function (s1, s2, f1, f2) {
  if (f2 == null) {
    return this.getIntersection2(s1, s2, f1);
  }

  var x1 = s1.x;
  var y1 = s1.y;
  var x2 = s2.x;
  var y2 = s2.y;
  var x3 = f1.x;
  var y3 = f1.y;
  var x4 = f2.x;
  var y4 = f2.y;
  var x = void 0,
      y = void 0; // intersection point
  var a1 = void 0,
      a2 = void 0,
      b1 = void 0,
      b2 = void 0,
      c1 = void 0,
      c2 = void 0; // coefficients of line eqns.
  var denom = void 0;

  a1 = y2 - y1;
  b1 = x1 - x2;
  c1 = x2 * y1 - x1 * y2; // { a1*x + b1*y + c1 = 0 is line 1 }

  a2 = y4 - y3;
  b2 = x3 - x4;
  c2 = x4 * y3 - x3 * y4; // { a2*x + b2*y + c2 = 0 is line 2 }

  denom = a1 * b2 - a2 * b1;

  if (denom === 0) {
    return null;
  }

  x = (b1 * c2 - b2 * c1) / denom;
  y = (a2 * c1 - a1 * c2) / denom;

  return new Point(x, y);
};

/**
 * This method finds and returns the angle of the vector from the + x-axis
 * in clockwise direction (compatible w/ Java coordinate system!).
 */
IGeometry.angleOfVector = function (Cx, Cy, Nx, Ny) {
  var C_angle = void 0;

  if (Cx !== Nx) {
    C_angle = Math.atan((Ny - Cy) / (Nx - Cx));

    if (Nx < Cx) {
      C_angle += Math.PI;
    } else if (Ny < Cy) {
      C_angle += this.TWO_PI;
    }
  } else if (Ny < Cy) {
    C_angle = this.ONE_AND_HALF_PI; // 270 degrees
  } else {
    C_angle = this.HALF_PI; // 90 degrees
  }

  return C_angle;
};

/**
 * This method checks whether the given two line segments (one with point
 * p1 and p2, the other with point p3 and p4) intersect at a point other
 * than these points.
 */
IGeometry.doIntersect = function (p1, p2, p3, p4) {
  var a = p1.x;
  var b = p1.y;
  var c = p2.x;
  var d = p2.y;
  var p = p3.x;
  var q = p3.y;
  var r = p4.x;
  var s = p4.y;
  var det = (c - a) * (s - q) - (r - p) * (d - b);

  if (det === 0) {
    return false;
  } else {
    var lambda = ((s - q) * (r - a) + (p - r) * (s - b)) / det;
    var gamma = ((b - d) * (r - a) + (c - a) * (s - b)) / det;
    return 0 < lambda && lambda < 1 && 0 < gamma && gamma < 1;
  }
};

/**
 * This method checks and calculates the intersection of 
 * a line segment and a circle.
 */
IGeometry.findCircleLineIntersections = function (Ex, Ey, Lx, Ly, Cx, Cy, r) {

  // E is the starting point of the ray,
  // L is the end point of the ray,
  // C is the center of sphere you're testing against
  // r is the radius of that sphere

  // Compute:
  // d = L - E ( Direction vector of ray, from start to end )
  // f = E - C ( Vector from center sphere to ray start )

  // Then the intersection is found by..
  // P = E + t * d
  // This is a parametric equation:
  // Px = Ex + tdx
  // Py = Ey + tdy

  // get a, b, c values
  var a = (Lx - Ex) * (Lx - Ex) + (Ly - Ey) * (Ly - Ey);
  var b = 2 * ((Ex - Cx) * (Lx - Ex) + (Ey - Cy) * (Ly - Ey));
  var c = (Ex - Cx) * (Ex - Cx) + (Ey - Cy) * (Ey - Cy) - r * r;

  // get discriminant
  var disc = b * b - 4 * a * c;
  if (disc >= 0) {
    // insert into quadratic formula
    var t1 = (-b + Math.sqrt(b * b - 4 * a * c)) / (2 * a);
    var t2 = (-b - Math.sqrt(b * b - 4 * a * c)) / (2 * a);
    var intersections = null;
    if (t1 >= 0 && t1 <= 1) {
      // t1 is the intersection, and it's closer than t2
      // (since t1 uses -b - discriminant)
      // Impale, Poke
      return [t1];
    }

    // here t1 didn't intersect so we are either started
    // inside the sphere or completely past it
    if (t2 >= 0 && t2 <= 1) {
      // ExitWound
      return [t2];
    }

    return intersections;
  } else return null;
};

// -----------------------------------------------------------------------------
// Section: Class Constants
// -----------------------------------------------------------------------------
/**
 * Some useful pre-calculated constants
 */
IGeometry.HALF_PI = 0.5 * Math.PI;
IGeometry.ONE_AND_HALF_PI = 1.5 * Math.PI;
IGeometry.TWO_PI = 2.0 * Math.PI;
IGeometry.THREE_PI = 3.0 * Math.PI;

module.exports = IGeometry;

/***/ }),
/* 9 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


function IMath() {}

/**
 * This method returns the sign of the input value.
 */
IMath.sign = function (value) {
  if (value > 0) {
    return 1;
  } else if (value < 0) {
    return -1;
  } else {
    return 0;
  }
};

IMath.floor = function (value) {
  return value < 0 ? Math.ceil(value) : Math.floor(value);
};

IMath.ceil = function (value) {
  return value < 0 ? Math.floor(value) : Math.ceil(value);
};

module.exports = IMath;

/***/ }),
/* 10 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


function Integer() {}

Integer.MAX_VALUE = 2147483647;
Integer.MIN_VALUE = -2147483648;

module.exports = Integer;

/***/ }),
/* 11 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

var nodeFrom = function nodeFrom(value) {
  return { value: value, next: null, prev: null };
};

var add = function add(prev, node, next, list) {
  if (prev !== null) {
    prev.next = node;
  } else {
    list.head = node;
  }

  if (next !== null) {
    next.prev = node;
  } else {
    list.tail = node;
  }

  node.prev = prev;
  node.next = next;

  list.length++;

  return node;
};

var _remove = function _remove(node, list) {
  var prev = node.prev,
      next = node.next;


  if (prev !== null) {
    prev.next = next;
  } else {
    list.head = next;
  }

  if (next !== null) {
    next.prev = prev;
  } else {
    list.tail = prev;
  }

  node.prev = node.next = null;

  list.length--;

  return node;
};

var LinkedList = function () {
  function LinkedList(vals) {
    var _this = this;

    _classCallCheck(this, LinkedList);

    this.length = 0;
    this.head = null;
    this.tail = null;

    if (vals != null) {
      vals.forEach(function (v) {
        return _this.push(v);
      });
    }
  }

  _createClass(LinkedList, [{
    key: "size",
    value: function size() {
      return this.length;
    }
  }, {
    key: "insertBefore",
    value: function insertBefore(val, otherNode) {
      return add(otherNode.prev, nodeFrom(val), otherNode, this);
    }
  }, {
    key: "insertAfter",
    value: function insertAfter(val, otherNode) {
      return add(otherNode, nodeFrom(val), otherNode.next, this);
    }
  }, {
    key: "insertNodeBefore",
    value: function insertNodeBefore(newNode, otherNode) {
      return add(otherNode.prev, newNode, otherNode, this);
    }
  }, {
    key: "insertNodeAfter",
    value: function insertNodeAfter(newNode, otherNode) {
      return add(otherNode, newNode, otherNode.next, this);
    }
  }, {
    key: "push",
    value: function push(val) {
      return add(this.tail, nodeFrom(val), null, this);
    }
  }, {
    key: "unshift",
    value: function unshift(val) {
      return add(null, nodeFrom(val), this.head, this);
    }
  }, {
    key: "remove",
    value: function remove(node) {
      return _remove(node, this);
    }
  }, {
    key: "pop",
    value: function pop() {
      return _remove(this.tail, this).value;
    }
  }, {
    key: "popNode",
    value: function popNode() {
      return _remove(this.tail, this);
    }
  }, {
    key: "shift",
    value: function shift() {
      return _remove(this.head, this).value;
    }
  }, {
    key: "shiftNode",
    value: function shiftNode() {
      return _remove(this.head, this);
    }
  }, {
    key: "get_object_at",
    value: function get_object_at(index) {
      if (index <= this.length()) {
        var i = 1;
        var current = this.head;
        while (i < index) {
          current = current.next;
          i++;
        }
        return current.value;
      }
    }
  }, {
    key: "set_object_at",
    value: function set_object_at(index, value) {
      if (index <= this.length()) {
        var i = 1;
        var current = this.head;
        while (i < index) {
          current = current.next;
          i++;
        }
        current.value = value;
      }
    }
  }]);

  return LinkedList;
}();

module.exports = LinkedList;

/***/ }),
/* 12 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


/*
 *This class is the javascript implementation of the Point.java class in jdk
 */
function Point(x, y, p) {
  this.x = null;
  this.y = null;
  if (x == null && y == null && p == null) {
    this.x = 0;
    this.y = 0;
  } else if (typeof x == 'number' && typeof y == 'number' && p == null) {
    this.x = x;
    this.y = y;
  } else if (x.constructor.name == 'Point' && y == null && p == null) {
    p = x;
    this.x = p.x;
    this.y = p.y;
  }
}

Point.prototype.getX = function () {
  return this.x;
};

Point.prototype.getY = function () {
  return this.y;
};

Point.prototype.getLocation = function () {
  return new Point(this.x, this.y);
};

Point.prototype.setLocation = function (x, y, p) {
  if (x.constructor.name == 'Point' && y == null && p == null) {
    p = x;
    this.setLocation(p.x, p.y);
  } else if (typeof x == 'number' && typeof y == 'number' && p == null) {
    //if both parameters are integer just move (x,y) location
    if (parseInt(x) == x && parseInt(y) == y) {
      this.move(x, y);
    } else {
      this.x = Math.floor(x + 0.5);
      this.y = Math.floor(y + 0.5);
    }
  }
};

Point.prototype.move = function (x, y) {
  this.x = x;
  this.y = y;
};

Point.prototype.translate = function (dx, dy) {
  this.x += dx;
  this.y += dy;
};

Point.prototype.equals = function (obj) {
  if (obj.constructor.name == "Point") {
    var pt = obj;
    return this.x == pt.x && this.y == pt.y;
  }
  return this == obj;
};

Point.prototype.toString = function () {
  return new Point().constructor.name + "[x=" + this.x + ",y=" + this.y + "]";
};

module.exports = Point;

/***/ }),
/* 13 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


function RectangleD(x, y, width, height) {
  this.x = 0;
  this.y = 0;
  this.width = 0;
  this.height = 0;

  if (x != null && y != null && width != null && height != null) {
    this.x = x;
    this.y = y;
    this.width = width;
    this.height = height;
  }
}

RectangleD.prototype.getX = function () {
  return this.x;
};

RectangleD.prototype.setX = function (x) {
  this.x = x;
};

RectangleD.prototype.getY = function () {
  return this.y;
};

RectangleD.prototype.setY = function (y) {
  this.y = y;
};

RectangleD.prototype.getWidth = function () {
  return this.width;
};

RectangleD.prototype.setWidth = function (width) {
  this.width = width;
};

RectangleD.prototype.getHeight = function () {
  return this.height;
};

RectangleD.prototype.setHeight = function (height) {
  this.height = height;
};

RectangleD.prototype.getRight = function () {
  return this.x + this.width;
};

RectangleD.prototype.getBottom = function () {
  return this.y + this.height;
};

RectangleD.prototype.intersects = function (a) {
  if (this.getRight() < a.x) {
    return false;
  }

  if (this.getBottom() < a.y) {
    return false;
  }

  if (a.getRight() < this.x) {
    return false;
  }

  if (a.getBottom() < this.y) {
    return false;
  }

  return true;
};

RectangleD.prototype.getCenterX = function () {
  return this.x + this.width / 2;
};

RectangleD.prototype.getMinX = function () {
  return this.getX();
};

RectangleD.prototype.getMaxX = function () {
  return this.getX() + this.width;
};

RectangleD.prototype.getCenterY = function () {
  return this.y + this.height / 2;
};

RectangleD.prototype.getMinY = function () {
  return this.getY();
};

RectangleD.prototype.getMaxY = function () {
  return this.getY() + this.height;
};

RectangleD.prototype.getWidthHalf = function () {
  return this.width / 2;
};

RectangleD.prototype.getHeightHalf = function () {
  return this.height / 2;
};

module.exports = RectangleD;

/***/ }),
/* 14 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var _typeof = typeof Symbol === "function" && typeof Symbol.iterator === "symbol" ? function (obj) { return typeof obj; } : function (obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; };

function UniqueIDGeneretor() {}

UniqueIDGeneretor.lastID = 0;

UniqueIDGeneretor.createID = function (obj) {
  if (UniqueIDGeneretor.isPrimitive(obj)) {
    return obj;
  }
  if (obj.uniqueID != null) {
    return obj.uniqueID;
  }
  obj.uniqueID = UniqueIDGeneretor.getString();
  UniqueIDGeneretor.lastID++;
  return obj.uniqueID;
};

UniqueIDGeneretor.getString = function (id) {
  if (id == null) id = UniqueIDGeneretor.lastID;
  return "Object#" + id + "";
};

UniqueIDGeneretor.isPrimitive = function (arg) {
  var type = typeof arg === "undefined" ? "undefined" : _typeof(arg);
  return arg == null || type != "object" && type != "function";
};

module.exports = UniqueIDGeneretor;

/***/ }),
/* 15 */
/***/ (function(module, exports, __nested_webpack_require_66321__) {

"use strict";


function _toConsumableArray(arr) { if (Array.isArray(arr)) { for (var i = 0, arr2 = Array(arr.length); i < arr.length; i++) { arr2[i] = arr[i]; } return arr2; } else { return Array.from(arr); } }

var LayoutConstants = __nested_webpack_require_66321__(0);
var LGraphManager = __nested_webpack_require_66321__(7);
var LNode = __nested_webpack_require_66321__(3);
var LEdge = __nested_webpack_require_66321__(1);
var LGraph = __nested_webpack_require_66321__(6);
var PointD = __nested_webpack_require_66321__(5);
var Transform = __nested_webpack_require_66321__(17);
var Emitter = __nested_webpack_require_66321__(29);

function Layout(isRemoteUse) {
  Emitter.call(this);

  //Layout Quality: 0:draft, 1:default, 2:proof
  this.layoutQuality = LayoutConstants.QUALITY;
  //Whether layout should create bendpoints as needed or not
  this.createBendsAsNeeded = LayoutConstants.DEFAULT_CREATE_BENDS_AS_NEEDED;
  //Whether layout should be incremental or not
  this.incremental = LayoutConstants.DEFAULT_INCREMENTAL;
  //Whether we animate from before to after layout node positions
  this.animationOnLayout = LayoutConstants.DEFAULT_ANIMATION_ON_LAYOUT;
  //Whether we animate the layout process or not
  this.animationDuringLayout = LayoutConstants.DEFAULT_ANIMATION_DURING_LAYOUT;
  //Number iterations that should be done between two successive animations
  this.animationPeriod = LayoutConstants.DEFAULT_ANIMATION_PERIOD;
  /**
   * Whether or not leaf nodes (non-compound nodes) are of uniform sizes. When
   * they are, both spring and repulsion forces between two leaf nodes can be
   * calculated without the expensive clipping point calculations, resulting
   * in major speed-up.
   */
  this.uniformLeafNodeSizes = LayoutConstants.DEFAULT_UNIFORM_LEAF_NODE_SIZES;
  /**
   * This is used for creation of bendpoints by using dummy nodes and edges.
   * Maps an LEdge to its dummy bendpoint path.
   */
  this.edgeToDummyNodes = new Map();
  this.graphManager = new LGraphManager(this);
  this.isLayoutFinished = false;
  this.isSubLayout = false;
  this.isRemoteUse = false;

  if (isRemoteUse != null) {
    this.isRemoteUse = isRemoteUse;
  }
}

Layout.RANDOM_SEED = 1;

Layout.prototype = Object.create(Emitter.prototype);

Layout.prototype.getGraphManager = function () {
  return this.graphManager;
};

Layout.prototype.getAllNodes = function () {
  return this.graphManager.getAllNodes();
};

Layout.prototype.getAllEdges = function () {
  return this.graphManager.getAllEdges();
};

Layout.prototype.getAllNodesToApplyGravitation = function () {
  return this.graphManager.getAllNodesToApplyGravitation();
};

Layout.prototype.newGraphManager = function () {
  var gm = new LGraphManager(this);
  this.graphManager = gm;
  return gm;
};

Layout.prototype.newGraph = function (vGraph) {
  return new LGraph(null, this.graphManager, vGraph);
};

Layout.prototype.newNode = function (vNode) {
  return new LNode(this.graphManager, vNode);
};

Layout.prototype.newEdge = function (vEdge) {
  return new LEdge(null, null, vEdge);
};

Layout.prototype.checkLayoutSuccess = function () {
  return this.graphManager.getRoot() == null || this.graphManager.getRoot().getNodes().length == 0 || this.graphManager.includesInvalidEdge();
};

Layout.prototype.runLayout = function () {
  this.isLayoutFinished = false;

  if (this.tilingPreLayout) {
    this.tilingPreLayout();
  }

  this.initParameters();
  var isLayoutSuccessfull;

  if (this.checkLayoutSuccess()) {
    isLayoutSuccessfull = false;
  } else {
    isLayoutSuccessfull = this.layout();
  }

  if (LayoutConstants.ANIMATE === 'during') {
    // If this is a 'during' layout animation. Layout is not finished yet. 
    // We need to perform these in index.js when layout is really finished.
    return false;
  }

  if (isLayoutSuccessfull) {
    if (!this.isSubLayout) {
      this.doPostLayout();
    }
  }

  if (this.tilingPostLayout) {
    this.tilingPostLayout();
  }

  this.isLayoutFinished = true;

  return isLayoutSuccessfull;
};

/**
 * This method performs the operations required after layout.
 */
Layout.prototype.doPostLayout = function () {
  //assert !isSubLayout : "Should not be called on sub-layout!";
  // Propagate geometric changes to v-level objects
  if (!this.incremental) {
    this.transform();
  }
  this.update();
};

/**
 * This method updates the geometry of the target graph according to
 * calculated layout.
 */
Layout.prototype.update2 = function () {
  // update bend points
  if (this.createBendsAsNeeded) {
    this.createBendpointsFromDummyNodes();

    // reset all edges, since the topology has changed
    this.graphManager.resetAllEdges();
  }

  // perform edge, node and root updates if layout is not called
  // remotely
  if (!this.isRemoteUse) {
    // update all edges
    var edge;
    var allEdges = this.graphManager.getAllEdges();
    for (var i = 0; i < allEdges.length; i++) {
      edge = allEdges[i];
      //      this.update(edge);
    }

    // recursively update nodes
    var node;
    var nodes = this.graphManager.getRoot().getNodes();
    for (var i = 0; i < nodes.length; i++) {
      node = nodes[i];
      //      this.update(node);
    }

    // update root graph
    this.update(this.graphManager.getRoot());
  }
};

Layout.prototype.update = function (obj) {
  if (obj == null) {
    this.update2();
  } else if (obj instanceof LNode) {
    var node = obj;
    if (node.getChild() != null) {
      // since node is compound, recursively update child nodes
      var nodes = node.getChild().getNodes();
      for (var i = 0; i < nodes.length; i++) {
        update(nodes[i]);
      }
    }

    // if the l-level node is associated with a v-level graph object,
    // then it is assumed that the v-level node implements the
    // interface Updatable.
    if (node.vGraphObject != null) {
      // cast to Updatable without any type check
      var vNode = node.vGraphObject;

      // call the update method of the interface
      vNode.update(node);
    }
  } else if (obj instanceof LEdge) {
    var edge = obj;
    // if the l-level edge is associated with a v-level graph object,
    // then it is assumed that the v-level edge implements the
    // interface Updatable.

    if (edge.vGraphObject != null) {
      // cast to Updatable without any type check
      var vEdge = edge.vGraphObject;

      // call the update method of the interface
      vEdge.update(edge);
    }
  } else if (obj instanceof LGraph) {
    var graph = obj;
    // if the l-level graph is associated with a v-level graph object,
    // then it is assumed that the v-level object implements the
    // interface Updatable.

    if (graph.vGraphObject != null) {
      // cast to Updatable without any type check
      var vGraph = graph.vGraphObject;

      // call the update method of the interface
      vGraph.update(graph);
    }
  }
};

/**
 * This method is used to set all layout parameters to default values
 * determined at compile time.
 */
Layout.prototype.initParameters = function () {
  if (!this.isSubLayout) {
    this.layoutQuality = LayoutConstants.QUALITY;
    this.animationDuringLayout = LayoutConstants.DEFAULT_ANIMATION_DURING_LAYOUT;
    this.animationPeriod = LayoutConstants.DEFAULT_ANIMATION_PERIOD;
    this.animationOnLayout = LayoutConstants.DEFAULT_ANIMATION_ON_LAYOUT;
    this.incremental = LayoutConstants.DEFAULT_INCREMENTAL;
    this.createBendsAsNeeded = LayoutConstants.DEFAULT_CREATE_BENDS_AS_NEEDED;
    this.uniformLeafNodeSizes = LayoutConstants.DEFAULT_UNIFORM_LEAF_NODE_SIZES;
  }

  if (this.animationDuringLayout) {
    this.animationOnLayout = false;
  }
};

Layout.prototype.transform = function (newLeftTop) {
  if (newLeftTop == undefined) {
    this.transform(new PointD(0, 0));
  } else {
    // create a transformation object (from Eclipse to layout). When an
    // inverse transform is applied, we get upper-left coordinate of the
    // drawing or the root graph at given input coordinate (some margins
    // already included in calculation of left-top).

    var trans = new Transform();
    var leftTop = this.graphManager.getRoot().updateLeftTop();

    if (leftTop != null) {
      trans.setWorldOrgX(newLeftTop.x);
      trans.setWorldOrgY(newLeftTop.y);

      trans.setDeviceOrgX(leftTop.x);
      trans.setDeviceOrgY(leftTop.y);

      var nodes = this.getAllNodes();
      var node;

      for (var i = 0; i < nodes.length; i++) {
        node = nodes[i];
        node.transform(trans);
      }
    }
  }
};

Layout.prototype.positionNodesRandomly = function (graph) {

  if (graph == undefined) {
    //assert !this.incremental;
    this.positionNodesRandomly(this.getGraphManager().getRoot());
    this.getGraphManager().getRoot().updateBounds(true);
  } else {
    var lNode;
    var childGraph;

    var nodes = graph.getNodes();
    for (var i = 0; i < nodes.length; i++) {
      lNode = nodes[i];
      childGraph = lNode.getChild();

      if (childGraph == null) {
        lNode.scatter();
      } else if (childGraph.getNodes().length == 0) {
        lNode.scatter();
      } else {
        this.positionNodesRandomly(childGraph);
        lNode.updateBounds();
      }
    }
  }
};

/**
 * This method returns a list of trees where each tree is represented as a
 * list of l-nodes. The method returns a list of size 0 when:
 * - The graph is not flat or
 * - One of the component(s) of the graph is not a tree.
 */
Layout.prototype.getFlatForest = function () {
  var flatForest = [];
  var isForest = true;

  // Quick reference for all nodes in the graph manager associated with
  // this layout. The list should not be changed.
  var allNodes = this.graphManager.getRoot().getNodes();

  // First be sure that the graph is flat
  var isFlat = true;

  for (var i = 0; i < allNodes.length; i++) {
    if (allNodes[i].getChild() != null) {
      isFlat = false;
    }
  }

  // Return empty forest if the graph is not flat.
  if (!isFlat) {
    return flatForest;
  }

  // Run BFS for each component of the graph.

  var visited = new Set();
  var toBeVisited = [];
  var parents = new Map();
  var unProcessedNodes = [];

  unProcessedNodes = unProcessedNodes.concat(allNodes);

  // Each iteration of this loop finds a component of the graph and
  // decides whether it is a tree or not. If it is a tree, adds it to the
  // forest and continued with the next component.

  while (unProcessedNodes.length > 0 && isForest) {
    toBeVisited.push(unProcessedNodes[0]);

    // Start the BFS. Each iteration of this loop visits a node in a
    // BFS manner.
    while (toBeVisited.length > 0 && isForest) {
      //pool operation
      var currentNode = toBeVisited[0];
      toBeVisited.splice(0, 1);
      visited.add(currentNode);

      // Traverse all neighbors of this node
      var neighborEdges = currentNode.getEdges();

      for (var i = 0; i < neighborEdges.length; i++) {
        var currentNeighbor = neighborEdges[i].getOtherEnd(currentNode);

        // If BFS is not growing from this neighbor.
        if (parents.get(currentNode) != currentNeighbor) {
          // We haven't previously visited this neighbor.
          if (!visited.has(currentNeighbor)) {
            toBeVisited.push(currentNeighbor);
            parents.set(currentNeighbor, currentNode);
          }
          // Since we have previously visited this neighbor and
          // this neighbor is not parent of currentNode, given
          // graph contains a component that is not tree, hence
          // it is not a forest.
          else {
              isForest = false;
              break;
            }
        }
      }
    }

    // The graph contains a component that is not a tree. Empty
    // previously found trees. The method will end.
    if (!isForest) {
      flatForest = [];
    }
    // Save currently visited nodes as a tree in our forest. Reset
    // visited and parents lists. Continue with the next component of
    // the graph, if any.
    else {
        var temp = [].concat(_toConsumableArray(visited));
        flatForest.push(temp);
        //flatForest = flatForest.concat(temp);
        //unProcessedNodes.removeAll(visited);
        for (var i = 0; i < temp.length; i++) {
          var value = temp[i];
          var index = unProcessedNodes.indexOf(value);
          if (index > -1) {
            unProcessedNodes.splice(index, 1);
          }
        }
        visited = new Set();
        parents = new Map();
      }
  }

  return flatForest;
};

/**
 * This method creates dummy nodes (an l-level node with minimal dimensions)
 * for the given edge (one per bendpoint). The existing l-level structure
 * is updated accordingly.
 */
Layout.prototype.createDummyNodesForBendpoints = function (edge) {
  var dummyNodes = [];
  var prev = edge.source;

  var graph = this.graphManager.calcLowestCommonAncestor(edge.source, edge.target);

  for (var i = 0; i < edge.bendpoints.length; i++) {
    // create new dummy node
    var dummyNode = this.newNode(null);
    dummyNode.setRect(new Point(0, 0), new Dimension(1, 1));

    graph.add(dummyNode);

    // create new dummy edge between prev and dummy node
    var dummyEdge = this.newEdge(null);
    this.graphManager.add(dummyEdge, prev, dummyNode);

    dummyNodes.add(dummyNode);
    prev = dummyNode;
  }

  var dummyEdge = this.newEdge(null);
  this.graphManager.add(dummyEdge, prev, edge.target);

  this.edgeToDummyNodes.set(edge, dummyNodes);

  // remove real edge from graph manager if it is inter-graph
  if (edge.isInterGraph()) {
    this.graphManager.remove(edge);
  }
  // else, remove the edge from the current graph
  else {
      graph.remove(edge);
    }

  return dummyNodes;
};

/**
 * This method creates bendpoints for edges from the dummy nodes
 * at l-level.
 */
Layout.prototype.createBendpointsFromDummyNodes = function () {
  var edges = [];
  edges = edges.concat(this.graphManager.getAllEdges());
  edges = [].concat(_toConsumableArray(this.edgeToDummyNodes.keys())).concat(edges);

  for (var k = 0; k < edges.length; k++) {
    var lEdge = edges[k];

    if (lEdge.bendpoints.length > 0) {
      var path = this.edgeToDummyNodes.get(lEdge);

      for (var i = 0; i < path.length; i++) {
        var dummyNode = path[i];
        var p = new PointD(dummyNode.getCenterX(), dummyNode.getCenterY());

        // update bendpoint's location according to dummy node
        var ebp = lEdge.bendpoints.get(i);
        ebp.x = p.x;
        ebp.y = p.y;

        // remove the dummy node, dummy edges incident with this
        // dummy node is also removed (within the remove method)
        dummyNode.getOwner().remove(dummyNode);
      }

      // add the real edge to graph
      this.graphManager.add(lEdge, lEdge.source, lEdge.target);
    }
  }
};

Layout.transform = function (sliderValue, defaultValue, minDiv, maxMul) {
  if (minDiv != undefined && maxMul != undefined) {
    var value = defaultValue;

    if (sliderValue <= 50) {
      var minValue = defaultValue / minDiv;
      value -= (defaultValue - minValue) / 50 * (50 - sliderValue);
    } else {
      var maxValue = defaultValue * maxMul;
      value += (maxValue - defaultValue) / 50 * (sliderValue - 50);
    }

    return value;
  } else {
    var a, b;

    if (sliderValue <= 50) {
      a = 9.0 * defaultValue / 500.0;
      b = defaultValue / 10.0;
    } else {
      a = 9.0 * defaultValue / 50.0;
      b = -8 * defaultValue;
    }

    return a * sliderValue + b;
  }
};

/**
 * This method finds and returns the center of the given nodes, assuming
 * that the given nodes form a tree in themselves.
 */
Layout.findCenterOfTree = function (nodes) {
  var list = [];
  list = list.concat(nodes);

  var removedNodes = [];
  var remainingDegrees = new Map();
  var foundCenter = false;
  var centerNode = null;

  if (list.length == 1 || list.length == 2) {
    foundCenter = true;
    centerNode = list[0];
  }

  for (var i = 0; i < list.length; i++) {
    var node = list[i];
    var degree = node.getNeighborsList().size;
    remainingDegrees.set(node, node.getNeighborsList().size);

    if (degree == 1) {
      removedNodes.push(node);
    }
  }

  var tempList = [];
  tempList = tempList.concat(removedNodes);

  while (!foundCenter) {
    var tempList2 = [];
    tempList2 = tempList2.concat(tempList);
    tempList = [];

    for (var i = 0; i < list.length; i++) {
      var node = list[i];

      var index = list.indexOf(node);
      if (index >= 0) {
        list.splice(index, 1);
      }

      var neighbours = node.getNeighborsList();

      neighbours.forEach(function (neighbour) {
        if (removedNodes.indexOf(neighbour) < 0) {
          var otherDegree = remainingDegrees.get(neighbour);
          var newDegree = otherDegree - 1;

          if (newDegree == 1) {
            tempList.push(neighbour);
          }

          remainingDegrees.set(neighbour, newDegree);
        }
      });
    }

    removedNodes = removedNodes.concat(tempList);

    if (list.length == 1 || list.length == 2) {
      foundCenter = true;
      centerNode = list[0];
    }
  }

  return centerNode;
};

/**
 * During the coarsening process, this layout may be referenced by two graph managers
 * this setter function grants access to change the currently being used graph manager
 */
Layout.prototype.setGraphManager = function (gm) {
  this.graphManager = gm;
};

module.exports = Layout;

/***/ }),
/* 16 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


function RandomSeed() {}
// adapted from: https://stackoverflow.com/a/19303725
RandomSeed.seed = 1;
RandomSeed.x = 0;

RandomSeed.nextDouble = function () {
  RandomSeed.x = Math.sin(RandomSeed.seed++) * 10000;
  return RandomSeed.x - Math.floor(RandomSeed.x);
};

module.exports = RandomSeed;

/***/ }),
/* 17 */
/***/ (function(module, exports, __nested_webpack_require_84109__) {

"use strict";


var PointD = __nested_webpack_require_84109__(5);

function Transform(x, y) {
  this.lworldOrgX = 0.0;
  this.lworldOrgY = 0.0;
  this.ldeviceOrgX = 0.0;
  this.ldeviceOrgY = 0.0;
  this.lworldExtX = 1.0;
  this.lworldExtY = 1.0;
  this.ldeviceExtX = 1.0;
  this.ldeviceExtY = 1.0;
}

Transform.prototype.getWorldOrgX = function () {
  return this.lworldOrgX;
};

Transform.prototype.setWorldOrgX = function (wox) {
  this.lworldOrgX = wox;
};

Transform.prototype.getWorldOrgY = function () {
  return this.lworldOrgY;
};

Transform.prototype.setWorldOrgY = function (woy) {
  this.lworldOrgY = woy;
};

Transform.prototype.getWorldExtX = function () {
  return this.lworldExtX;
};

Transform.prototype.setWorldExtX = function (wex) {
  this.lworldExtX = wex;
};

Transform.prototype.getWorldExtY = function () {
  return this.lworldExtY;
};

Transform.prototype.setWorldExtY = function (wey) {
  this.lworldExtY = wey;
};

/* Device related */

Transform.prototype.getDeviceOrgX = function () {
  return this.ldeviceOrgX;
};

Transform.prototype.setDeviceOrgX = function (dox) {
  this.ldeviceOrgX = dox;
};

Transform.prototype.getDeviceOrgY = function () {
  return this.ldeviceOrgY;
};

Transform.prototype.setDeviceOrgY = function (doy) {
  this.ldeviceOrgY = doy;
};

Transform.prototype.getDeviceExtX = function () {
  return this.ldeviceExtX;
};

Transform.prototype.setDeviceExtX = function (dex) {
  this.ldeviceExtX = dex;
};

Transform.prototype.getDeviceExtY = function () {
  return this.ldeviceExtY;
};

Transform.prototype.setDeviceExtY = function (dey) {
  this.ldeviceExtY = dey;
};

Transform.prototype.transformX = function (x) {
  var xDevice = 0.0;
  var worldExtX = this.lworldExtX;
  if (worldExtX != 0.0) {
    xDevice = this.ldeviceOrgX + (x - this.lworldOrgX) * this.ldeviceExtX / worldExtX;
  }

  return xDevice;
};

Transform.prototype.transformY = function (y) {
  var yDevice = 0.0;
  var worldExtY = this.lworldExtY;
  if (worldExtY != 0.0) {
    yDevice = this.ldeviceOrgY + (y - this.lworldOrgY) * this.ldeviceExtY / worldExtY;
  }

  return yDevice;
};

Transform.prototype.inverseTransformX = function (x) {
  var xWorld = 0.0;
  var deviceExtX = this.ldeviceExtX;
  if (deviceExtX != 0.0) {
    xWorld = this.lworldOrgX + (x - this.ldeviceOrgX) * this.lworldExtX / deviceExtX;
  }

  return xWorld;
};

Transform.prototype.inverseTransformY = function (y) {
  var yWorld = 0.0;
  var deviceExtY = this.ldeviceExtY;
  if (deviceExtY != 0.0) {
    yWorld = this.lworldOrgY + (y - this.ldeviceOrgY) * this.lworldExtY / deviceExtY;
  }
  return yWorld;
};

Transform.prototype.inverseTransformPoint = function (inPoint) {
  var outPoint = new PointD(this.inverseTransformX(inPoint.x), this.inverseTransformY(inPoint.y));
  return outPoint;
};

module.exports = Transform;

/***/ }),
/* 18 */
/***/ (function(module, exports, __nested_webpack_require_86996__) {

"use strict";


function _toConsumableArray(arr) { if (Array.isArray(arr)) { for (var i = 0, arr2 = Array(arr.length); i < arr.length; i++) { arr2[i] = arr[i]; } return arr2; } else { return Array.from(arr); } }

var Layout = __nested_webpack_require_86996__(15);
var FDLayoutConstants = __nested_webpack_require_86996__(4);
var LayoutConstants = __nested_webpack_require_86996__(0);
var IGeometry = __nested_webpack_require_86996__(8);
var IMath = __nested_webpack_require_86996__(9);

function FDLayout() {
  Layout.call(this);

  this.useSmartIdealEdgeLengthCalculation = FDLayoutConstants.DEFAULT_USE_SMART_IDEAL_EDGE_LENGTH_CALCULATION;
  this.gravityConstant = FDLayoutConstants.DEFAULT_GRAVITY_STRENGTH;
  this.compoundGravityConstant = FDLayoutConstants.DEFAULT_COMPOUND_GRAVITY_STRENGTH;
  this.gravityRangeFactor = FDLayoutConstants.DEFAULT_GRAVITY_RANGE_FACTOR;
  this.compoundGravityRangeFactor = FDLayoutConstants.DEFAULT_COMPOUND_GRAVITY_RANGE_FACTOR;
  this.displacementThresholdPerNode = 3.0 * FDLayoutConstants.DEFAULT_EDGE_LENGTH / 100;
  this.coolingFactor = FDLayoutConstants.DEFAULT_COOLING_FACTOR_INCREMENTAL;
  this.initialCoolingFactor = FDLayoutConstants.DEFAULT_COOLING_FACTOR_INCREMENTAL;
  this.totalDisplacement = 0.0;
  this.oldTotalDisplacement = 0.0;
  this.maxIterations = FDLayoutConstants.MAX_ITERATIONS;
}

FDLayout.prototype = Object.create(Layout.prototype);

for (var prop in Layout) {
  FDLayout[prop] = Layout[prop];
}

FDLayout.prototype.initParameters = function () {
  Layout.prototype.initParameters.call(this, arguments);

  this.totalIterations = 0;
  this.notAnimatedIterations = 0;

  this.useFRGridVariant = FDLayoutConstants.DEFAULT_USE_SMART_REPULSION_RANGE_CALCULATION;

  this.grid = [];
};

FDLayout.prototype.calcIdealEdgeLengths = function () {
  var edge;
  var originalIdealLength;
  var lcaDepth;
  var source;
  var target;
  var sizeOfSourceInLca;
  var sizeOfTargetInLca;

  var allEdges = this.getGraphManager().getAllEdges();
  for (var i = 0; i < allEdges.length; i++) {
    edge = allEdges[i];

    originalIdealLength = edge.idealLength;

    if (edge.isInterGraph) {
      source = edge.getSource();
      target = edge.getTarget();

      sizeOfSourceInLca = edge.getSourceInLca().getEstimatedSize();
      sizeOfTargetInLca = edge.getTargetInLca().getEstimatedSize();

      if (this.useSmartIdealEdgeLengthCalculation) {
        edge.idealLength += sizeOfSourceInLca + sizeOfTargetInLca - 2 * LayoutConstants.SIMPLE_NODE_SIZE;
      }

      lcaDepth = edge.getLca().getInclusionTreeDepth();

      edge.idealLength += originalIdealLength * FDLayoutConstants.PER_LEVEL_IDEAL_EDGE_LENGTH_FACTOR * (source.getInclusionTreeDepth() + target.getInclusionTreeDepth() - 2 * lcaDepth);
    }
  }
};

FDLayout.prototype.initSpringEmbedder = function () {

  var s = this.getAllNodes().length;
  if (this.incremental) {
    if (s > FDLayoutConstants.ADAPTATION_LOWER_NODE_LIMIT) {
      this.coolingFactor = Math.max(this.coolingFactor * FDLayoutConstants.COOLING_ADAPTATION_FACTOR, this.coolingFactor - (s - FDLayoutConstants.ADAPTATION_LOWER_NODE_LIMIT) / (FDLayoutConstants.ADAPTATION_UPPER_NODE_LIMIT - FDLayoutConstants.ADAPTATION_LOWER_NODE_LIMIT) * this.coolingFactor * (1 - FDLayoutConstants.COOLING_ADAPTATION_FACTOR));
    }
    this.maxNodeDisplacement = FDLayoutConstants.MAX_NODE_DISPLACEMENT_INCREMENTAL;
  } else {
    if (s > FDLayoutConstants.ADAPTATION_LOWER_NODE_LIMIT) {
      this.coolingFactor = Math.max(FDLayoutConstants.COOLING_ADAPTATION_FACTOR, 1.0 - (s - FDLayoutConstants.ADAPTATION_LOWER_NODE_LIMIT) / (FDLayoutConstants.ADAPTATION_UPPER_NODE_LIMIT - FDLayoutConstants.ADAPTATION_LOWER_NODE_LIMIT) * (1 - FDLayoutConstants.COOLING_ADAPTATION_FACTOR));
    } else {
      this.coolingFactor = 1.0;
    }
    this.initialCoolingFactor = this.coolingFactor;
    this.maxNodeDisplacement = FDLayoutConstants.MAX_NODE_DISPLACEMENT;
  }

  this.maxIterations = Math.max(this.getAllNodes().length * 5, this.maxIterations);

  // Reassign this attribute by using new constant value
  this.displacementThresholdPerNode = 3.0 * FDLayoutConstants.DEFAULT_EDGE_LENGTH / 100;
  this.totalDisplacementThreshold = this.displacementThresholdPerNode * this.getAllNodes().length;

  this.repulsionRange = this.calcRepulsionRange();
};

FDLayout.prototype.calcSpringForces = function () {
  var lEdges = this.getAllEdges();
  var edge;

  for (var i = 0; i < lEdges.length; i++) {
    edge = lEdges[i];

    this.calcSpringForce(edge, edge.idealLength);
  }
};

FDLayout.prototype.calcRepulsionForces = function () {
  var gridUpdateAllowed = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : true;
  var forceToNodeSurroundingUpdate = arguments.length > 1 && arguments[1] !== undefined ? arguments[1] : false;

  var i, j;
  var nodeA, nodeB;
  var lNodes = this.getAllNodes();
  var processedNodeSet;

  if (this.useFRGridVariant) {
    if (this.totalIterations % FDLayoutConstants.GRID_CALCULATION_CHECK_PERIOD == 1 && gridUpdateAllowed) {
      this.updateGrid();
    }

    processedNodeSet = new Set();

    // calculate repulsion forces between each nodes and its surrounding
    for (i = 0; i < lNodes.length; i++) {
      nodeA = lNodes[i];
      this.calculateRepulsionForceOfANode(nodeA, processedNodeSet, gridUpdateAllowed, forceToNodeSurroundingUpdate);
      processedNodeSet.add(nodeA);
    }
  } else {
    for (i = 0; i < lNodes.length; i++) {
      nodeA = lNodes[i];

      for (j = i + 1; j < lNodes.length; j++) {
        nodeB = lNodes[j];

        // If both nodes are not members of the same graph, skip.
        if (nodeA.getOwner() != nodeB.getOwner()) {
          continue;
        }

        this.calcRepulsionForce(nodeA, nodeB);
      }
    }
  }
};

FDLayout.prototype.calcGravitationalForces = function () {
  var node;
  var lNodes = this.getAllNodesToApplyGravitation();

  for (var i = 0; i < lNodes.length; i++) {
    node = lNodes[i];
    this.calcGravitationalForce(node);
  }
};

FDLayout.prototype.moveNodes = function () {
  var lNodes = this.getAllNodes();
  var node;

  for (var i = 0; i < lNodes.length; i++) {
    node = lNodes[i];
    node.move();
  }
};

FDLayout.prototype.calcSpringForce = function (edge, idealLength) {
  var sourceNode = edge.getSource();
  var targetNode = edge.getTarget();

  var length;
  var springForce;
  var springForceX;
  var springForceY;

  // Update edge length
  if (this.uniformLeafNodeSizes && sourceNode.getChild() == null && targetNode.getChild() == null) {
    edge.updateLengthSimple();
  } else {
    edge.updateLength();

    if (edge.isOverlapingSourceAndTarget) {
      return;
    }
  }

  length = edge.getLength();

  if (length == 0) return;

  // Calculate spring forces
  springForce = edge.edgeElasticity * (length - idealLength);

  // Project force onto x and y axes
  springForceX = springForce * (edge.lengthX / length);
  springForceY = springForce * (edge.lengthY / length);

  // Apply forces on the end nodes
  sourceNode.springForceX += springForceX;
  sourceNode.springForceY += springForceY;
  targetNode.springForceX -= springForceX;
  targetNode.springForceY -= springForceY;
};

FDLayout.prototype.calcRepulsionForce = function (nodeA, nodeB) {
  var rectA = nodeA.getRect();
  var rectB = nodeB.getRect();
  var overlapAmount = new Array(2);
  var clipPoints = new Array(4);
  var distanceX;
  var distanceY;
  var distanceSquared;
  var distance;
  var repulsionForce;
  var repulsionForceX;
  var repulsionForceY;

  if (rectA.intersects(rectB)) // two nodes overlap
    {
      // calculate separation amount in x and y directions
      IGeometry.calcSeparationAmount(rectA, rectB, overlapAmount, FDLayoutConstants.DEFAULT_EDGE_LENGTH / 2.0);

      repulsionForceX = 2 * overlapAmount[0];
      repulsionForceY = 2 * overlapAmount[1];

      var childrenConstant = nodeA.noOfChildren * nodeB.noOfChildren / (nodeA.noOfChildren + nodeB.noOfChildren);

      // Apply forces on the two nodes
      nodeA.repulsionForceX -= childrenConstant * repulsionForceX;
      nodeA.repulsionForceY -= childrenConstant * repulsionForceY;
      nodeB.repulsionForceX += childrenConstant * repulsionForceX;
      nodeB.repulsionForceY += childrenConstant * repulsionForceY;
    } else // no overlap
    {
      // calculate distance

      if (this.uniformLeafNodeSizes && nodeA.getChild() == null && nodeB.getChild() == null) // simply base repulsion on distance of node centers
        {
          distanceX = rectB.getCenterX() - rectA.getCenterX();
          distanceY = rectB.getCenterY() - rectA.getCenterY();
        } else // use clipping points
        {
          IGeometry.getIntersection(rectA, rectB, clipPoints);

          distanceX = clipPoints[2] - clipPoints[0];
          distanceY = clipPoints[3] - clipPoints[1];
        }

      // No repulsion range. FR grid variant should take care of this.
      if (Math.abs(distanceX) < FDLayoutConstants.MIN_REPULSION_DIST) {
        distanceX = IMath.sign(distanceX) * FDLayoutConstants.MIN_REPULSION_DIST;
      }

      if (Math.abs(distanceY) < FDLayoutConstants.MIN_REPULSION_DIST) {
        distanceY = IMath.sign(distanceY) * FDLayoutConstants.MIN_REPULSION_DIST;
      }

      distanceSquared = distanceX * distanceX + distanceY * distanceY;
      distance = Math.sqrt(distanceSquared);

      // Here we use half of the nodes' repulsion values for backward compatibility
      repulsionForce = (nodeA.nodeRepulsion / 2 + nodeB.nodeRepulsion / 2) * nodeA.noOfChildren * nodeB.noOfChildren / distanceSquared;

      // Project force onto x and y axes
      repulsionForceX = repulsionForce * distanceX / distance;
      repulsionForceY = repulsionForce * distanceY / distance;

      // Apply forces on the two nodes    
      nodeA.repulsionForceX -= repulsionForceX;
      nodeA.repulsionForceY -= repulsionForceY;
      nodeB.repulsionForceX += repulsionForceX;
      nodeB.repulsionForceY += repulsionForceY;
    }
};

FDLayout.prototype.calcGravitationalForce = function (node) {
  var ownerGraph;
  var ownerCenterX;
  var ownerCenterY;
  var distanceX;
  var distanceY;
  var absDistanceX;
  var absDistanceY;
  var estimatedSize;
  ownerGraph = node.getOwner();

  ownerCenterX = (ownerGraph.getRight() + ownerGraph.getLeft()) / 2;
  ownerCenterY = (ownerGraph.getTop() + ownerGraph.getBottom()) / 2;
  distanceX = node.getCenterX() - ownerCenterX;
  distanceY = node.getCenterY() - ownerCenterY;
  absDistanceX = Math.abs(distanceX) + node.getWidth() / 2;
  absDistanceY = Math.abs(distanceY) + node.getHeight() / 2;

  if (node.getOwner() == this.graphManager.getRoot()) // in the root graph
    {
      estimatedSize = ownerGraph.getEstimatedSize() * this.gravityRangeFactor;

      if (absDistanceX > estimatedSize || absDistanceY > estimatedSize) {
        node.gravitationForceX = -this.gravityConstant * distanceX;
        node.gravitationForceY = -this.gravityConstant * distanceY;
      }
    } else // inside a compound
    {
      estimatedSize = ownerGraph.getEstimatedSize() * this.compoundGravityRangeFactor;

      if (absDistanceX > estimatedSize || absDistanceY > estimatedSize) {
        node.gravitationForceX = -this.gravityConstant * distanceX * this.compoundGravityConstant;
        node.gravitationForceY = -this.gravityConstant * distanceY * this.compoundGravityConstant;
      }
    }
};

FDLayout.prototype.isConverged = function () {
  var converged;
  var oscilating = false;

  if (this.totalIterations > this.maxIterations / 3) {
    oscilating = Math.abs(this.totalDisplacement - this.oldTotalDisplacement) < 2;
  }

  converged = this.totalDisplacement < this.totalDisplacementThreshold;

  this.oldTotalDisplacement = this.totalDisplacement;

  return converged || oscilating;
};

FDLayout.prototype.animate = function () {
  if (this.animationDuringLayout && !this.isSubLayout) {
    if (this.notAnimatedIterations == this.animationPeriod) {
      this.update();
      this.notAnimatedIterations = 0;
    } else {
      this.notAnimatedIterations++;
    }
  }
};

//This method calculates the number of children (weight) for all nodes
FDLayout.prototype.calcNoOfChildrenForAllNodes = function () {
  var node;
  var allNodes = this.graphManager.getAllNodes();

  for (var i = 0; i < allNodes.length; i++) {
    node = allNodes[i];
    node.noOfChildren = node.getNoOfChildren();
  }
};

// -----------------------------------------------------------------------------
// Section: FR-Grid Variant Repulsion Force Calculation
// -----------------------------------------------------------------------------

FDLayout.prototype.calcGrid = function (graph) {

  var sizeX = 0;
  var sizeY = 0;

  sizeX = parseInt(Math.ceil((graph.getRight() - graph.getLeft()) / this.repulsionRange));
  sizeY = parseInt(Math.ceil((graph.getBottom() - graph.getTop()) / this.repulsionRange));

  var grid = new Array(sizeX);

  for (var i = 0; i < sizeX; i++) {
    grid[i] = new Array(sizeY);
  }

  for (var i = 0; i < sizeX; i++) {
    for (var j = 0; j < sizeY; j++) {
      grid[i][j] = new Array();
    }
  }

  return grid;
};

FDLayout.prototype.addNodeToGrid = function (v, left, top) {

  var startX = 0;
  var finishX = 0;
  var startY = 0;
  var finishY = 0;

  startX = parseInt(Math.floor((v.getRect().x - left) / this.repulsionRange));
  finishX = parseInt(Math.floor((v.getRect().width + v.getRect().x - left) / this.repulsionRange));
  startY = parseInt(Math.floor((v.getRect().y - top) / this.repulsionRange));
  finishY = parseInt(Math.floor((v.getRect().height + v.getRect().y - top) / this.repulsionRange));

  for (var i = startX; i <= finishX; i++) {
    for (var j = startY; j <= finishY; j++) {
      this.grid[i][j].push(v);
      v.setGridCoordinates(startX, finishX, startY, finishY);
    }
  }
};

FDLayout.prototype.updateGrid = function () {
  var i;
  var nodeA;
  var lNodes = this.getAllNodes();

  this.grid = this.calcGrid(this.graphManager.getRoot());

  // put all nodes to proper grid cells
  for (i = 0; i < lNodes.length; i++) {
    nodeA = lNodes[i];
    this.addNodeToGrid(nodeA, this.graphManager.getRoot().getLeft(), this.graphManager.getRoot().getTop());
  }
};

FDLayout.prototype.calculateRepulsionForceOfANode = function (nodeA, processedNodeSet, gridUpdateAllowed, forceToNodeSurroundingUpdate) {

  if (this.totalIterations % FDLayoutConstants.GRID_CALCULATION_CHECK_PERIOD == 1 && gridUpdateAllowed || forceToNodeSurroundingUpdate) {
    var surrounding = new Set();
    nodeA.surrounding = new Array();
    var nodeB;
    var grid = this.grid;

    for (var i = nodeA.startX - 1; i < nodeA.finishX + 2; i++) {
      for (var j = nodeA.startY - 1; j < nodeA.finishY + 2; j++) {
        if (!(i < 0 || j < 0 || i >= grid.length || j >= grid[0].length)) {
          for (var k = 0; k < grid[i][j].length; k++) {
            nodeB = grid[i][j][k];

            // If both nodes are not members of the same graph, 
            // or both nodes are the same, skip.
            if (nodeA.getOwner() != nodeB.getOwner() || nodeA == nodeB) {
              continue;
            }

            // check if the repulsion force between
            // nodeA and nodeB has already been calculated
            if (!processedNodeSet.has(nodeB) && !surrounding.has(nodeB)) {
              var distanceX = Math.abs(nodeA.getCenterX() - nodeB.getCenterX()) - (nodeA.getWidth() / 2 + nodeB.getWidth() / 2);
              var distanceY = Math.abs(nodeA.getCenterY() - nodeB.getCenterY()) - (nodeA.getHeight() / 2 + nodeB.getHeight() / 2);

              // if the distance between nodeA and nodeB 
              // is less then calculation range
              if (distanceX <= this.repulsionRange && distanceY <= this.repulsionRange) {
                //then add nodeB to surrounding of nodeA
                surrounding.add(nodeB);
              }
            }
          }
        }
      }
    }

    nodeA.surrounding = [].concat(_toConsumableArray(surrounding));
  }
  for (i = 0; i < nodeA.surrounding.length; i++) {
    this.calcRepulsionForce(nodeA, nodeA.surrounding[i]);
  }
};

FDLayout.prototype.calcRepulsionRange = function () {
  return 0.0;
};

module.exports = FDLayout;

/***/ }),
/* 19 */
/***/ (function(module, exports, __nested_webpack_require_103214__) {

"use strict";


var LEdge = __nested_webpack_require_103214__(1);
var FDLayoutConstants = __nested_webpack_require_103214__(4);

function FDLayoutEdge(source, target, vEdge) {
  LEdge.call(this, source, target, vEdge);

  // Ideal length and elasticity value for this edge
  this.idealLength = FDLayoutConstants.DEFAULT_EDGE_LENGTH;
  this.edgeElasticity = FDLayoutConstants.DEFAULT_SPRING_STRENGTH;
}

FDLayoutEdge.prototype = Object.create(LEdge.prototype);

for (var prop in LEdge) {
  FDLayoutEdge[prop] = LEdge[prop];
}

module.exports = FDLayoutEdge;

/***/ }),
/* 20 */
/***/ (function(module, exports, __nested_webpack_require_103820__) {

"use strict";


var LNode = __nested_webpack_require_103820__(3);
var FDLayoutConstants = __nested_webpack_require_103820__(4);

function FDLayoutNode(gm, loc, size, vNode) {
  // alternative constructor is handled inside LNode
  LNode.call(this, gm, loc, size, vNode);

  // Repulsion value of this node
  this.nodeRepulsion = FDLayoutConstants.DEFAULT_REPULSION_STRENGTH;

  //Spring, repulsion and gravitational forces acting on this node
  this.springForceX = 0;
  this.springForceY = 0;
  this.repulsionForceX = 0;
  this.repulsionForceY = 0;
  this.gravitationForceX = 0;
  this.gravitationForceY = 0;
  //Amount by which this node is to be moved in this iteration
  this.displacementX = 0;
  this.displacementY = 0;

  //Start and finish grid coordinates that this node is fallen into
  this.startX = 0;
  this.finishX = 0;
  this.startY = 0;
  this.finishY = 0;

  //Geometric neighbors of this node
  this.surrounding = [];
}

FDLayoutNode.prototype = Object.create(LNode.prototype);

for (var prop in LNode) {
  FDLayoutNode[prop] = LNode[prop];
}

FDLayoutNode.prototype.setGridCoordinates = function (_startX, _finishX, _startY, _finishY) {
  this.startX = _startX;
  this.finishX = _finishX;
  this.startY = _startY;
  this.finishY = _finishY;
};

module.exports = FDLayoutNode;

/***/ }),
/* 21 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


function DimensionD(width, height) {
  this.width = 0;
  this.height = 0;
  if (width !== null && height !== null) {
    this.height = height;
    this.width = width;
  }
}

DimensionD.prototype.getWidth = function () {
  return this.width;
};

DimensionD.prototype.setWidth = function (width) {
  this.width = width;
};

DimensionD.prototype.getHeight = function () {
  return this.height;
};

DimensionD.prototype.setHeight = function (height) {
  this.height = height;
};

module.exports = DimensionD;

/***/ }),
/* 22 */
/***/ (function(module, exports, __nested_webpack_require_105759__) {

"use strict";


var UniqueIDGeneretor = __nested_webpack_require_105759__(14);

function HashMap() {
  this.map = {};
  this.keys = [];
}

HashMap.prototype.put = function (key, value) {
  var theId = UniqueIDGeneretor.createID(key);
  if (!this.contains(theId)) {
    this.map[theId] = value;
    this.keys.push(key);
  }
};

HashMap.prototype.contains = function (key) {
  var theId = UniqueIDGeneretor.createID(key);
  return this.map[key] != null;
};

HashMap.prototype.get = function (key) {
  var theId = UniqueIDGeneretor.createID(key);
  return this.map[theId];
};

HashMap.prototype.keySet = function () {
  return this.keys;
};

module.exports = HashMap;

/***/ }),
/* 23 */
/***/ (function(module, exports, __nested_webpack_require_106487__) {

"use strict";


var UniqueIDGeneretor = __nested_webpack_require_106487__(14);

function HashSet() {
  this.set = {};
}
;

HashSet.prototype.add = function (obj) {
  var theId = UniqueIDGeneretor.createID(obj);
  if (!this.contains(theId)) this.set[theId] = obj;
};

HashSet.prototype.remove = function (obj) {
  delete this.set[UniqueIDGeneretor.createID(obj)];
};

HashSet.prototype.clear = function () {
  this.set = {};
};

HashSet.prototype.contains = function (obj) {
  return this.set[UniqueIDGeneretor.createID(obj)] == obj;
};

HashSet.prototype.isEmpty = function () {
  return this.size() === 0;
};

HashSet.prototype.size = function () {
  return Object.keys(this.set).length;
};

//concats this.set to the given list
HashSet.prototype.addAllTo = function (list) {
  var keys = Object.keys(this.set);
  var length = keys.length;
  for (var i = 0; i < length; i++) {
    list.push(this.set[keys[i]]);
  }
};

HashSet.prototype.size = function () {
  return Object.keys(this.set).length;
};

HashSet.prototype.addAll = function (list) {
  var s = list.length;
  for (var i = 0; i < s; i++) {
    var v = list[i];
    this.add(v);
  }
};

module.exports = HashSet;

/***/ }),
/* 24 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


// Some matrix (1d and 2d array) operations
function Matrix() {}

/**
 * matrix multiplication
 * array1, array2 and result are 2d arrays
 */
Matrix.multMat = function (array1, array2) {
  var result = [];

  for (var i = 0; i < array1.length; i++) {
    result[i] = [];
    for (var j = 0; j < array2[0].length; j++) {
      result[i][j] = 0;
      for (var k = 0; k < array1[0].length; k++) {
        result[i][j] += array1[i][k] * array2[k][j];
      }
    }
  }
  return result;
};

/**
 * matrix transpose
 * array and result are 2d arrays
 */
Matrix.transpose = function (array) {
  var result = [];

  for (var i = 0; i < array[0].length; i++) {
    result[i] = [];
    for (var j = 0; j < array.length; j++) {
      result[i][j] = array[j][i];
    }
  }

  return result;
};

/**
 * multiply array with constant
 * array and result are 1d arrays
 */
Matrix.multCons = function (array, constant) {
  var result = [];

  for (var i = 0; i < array.length; i++) {
    result[i] = array[i] * constant;
  }

  return result;
};

/**
 * substract two arrays
 * array1, array2 and result are 1d arrays
 */
Matrix.minusOp = function (array1, array2) {
  var result = [];

  for (var i = 0; i < array1.length; i++) {
    result[i] = array1[i] - array2[i];
  }

  return result;
};

/**
 * dot product of two arrays with same size
 * array1 and array2 are 1d arrays
 */
Matrix.dotProduct = function (array1, array2) {
  var product = 0;

  for (var i = 0; i < array1.length; i++) {
    product += array1[i] * array2[i];
  }

  return product;
};

/**
 * magnitude of an array
 * array is 1d array
 */
Matrix.mag = function (array) {
  return Math.sqrt(this.dotProduct(array, array));
};

/**
 * normalization of an array
 * array and result are 1d array
 */
Matrix.normalize = function (array) {
  var result = [];
  var magnitude = this.mag(array);

  for (var i = 0; i < array.length; i++) {
    result[i] = array[i] / magnitude;
  }

  return result;
};

/**
 * multiply an array with centering matrix
 * array and result are 1d array
 */
Matrix.multGamma = function (array) {
  var result = [];
  var sum = 0;

  for (var i = 0; i < array.length; i++) {
    sum += array[i];
  }

  sum *= -1 / array.length;

  for (var _i = 0; _i < array.length; _i++) {
    result[_i] = sum + array[_i];
  }
  return result;
};

/**
 * a special matrix multiplication
 * result = 0.5 * C * INV * C^T * array
 * array and result are 1d, C and INV are 2d arrays
 */
Matrix.multL = function (array, C, INV) {
  var result = [];
  var temp1 = [];
  var temp2 = [];

  // multiply by C^T
  for (var i = 0; i < C[0].length; i++) {
    var sum = 0;
    for (var j = 0; j < C.length; j++) {
      sum += -0.5 * C[j][i] * array[j];
    }
    temp1[i] = sum;
  }
  // multiply the result by INV
  for (var _i2 = 0; _i2 < INV.length; _i2++) {
    var _sum = 0;
    for (var _j = 0; _j < INV.length; _j++) {
      _sum += INV[_i2][_j] * temp1[_j];
    }
    temp2[_i2] = _sum;
  }
  // multiply the result by C
  for (var _i3 = 0; _i3 < C.length; _i3++) {
    var _sum2 = 0;
    for (var _j2 = 0; _j2 < C[0].length; _j2++) {
      _sum2 += C[_i3][_j2] * temp2[_j2];
    }
    result[_i3] = _sum2;
  }

  return result;
};

module.exports = Matrix;

/***/ }),
/* 25 */
/***/ (function(module, exports, __nested_webpack_require_111039__) {

"use strict";


var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

/**
 * A classic Quicksort algorithm with Hoare's partition
 * - Works also on LinkedList objects
 *
 * Copyright: i-Vis Research Group, Bilkent University, 2007 - present
 */

var LinkedList = __nested_webpack_require_111039__(11);

var Quicksort = function () {
    function Quicksort(A, compareFunction) {
        _classCallCheck(this, Quicksort);

        if (compareFunction !== null || compareFunction !== undefined) this.compareFunction = this._defaultCompareFunction;

        var length = void 0;
        if (A instanceof LinkedList) length = A.size();else length = A.length;

        this._quicksort(A, 0, length - 1);
    }

    _createClass(Quicksort, [{
        key: '_quicksort',
        value: function _quicksort(A, p, r) {
            if (p < r) {
                var q = this._partition(A, p, r);
                this._quicksort(A, p, q);
                this._quicksort(A, q + 1, r);
            }
        }
    }, {
        key: '_partition',
        value: function _partition(A, p, r) {
            var x = this._get(A, p);
            var i = p;
            var j = r;
            while (true) {
                while (this.compareFunction(x, this._get(A, j))) {
                    j--;
                }while (this.compareFunction(this._get(A, i), x)) {
                    i++;
                }if (i < j) {
                    this._swap(A, i, j);
                    i++;
                    j--;
                } else return j;
            }
        }
    }, {
        key: '_get',
        value: function _get(object, index) {
            if (object instanceof LinkedList) return object.get_object_at(index);else return object[index];
        }
    }, {
        key: '_set',
        value: function _set(object, index, value) {
            if (object instanceof LinkedList) object.set_object_at(index, value);else object[index] = value;
        }
    }, {
        key: '_swap',
        value: function _swap(A, i, j) {
            var temp = this._get(A, i);
            this._set(A, i, this._get(A, j));
            this._set(A, j, temp);
        }
    }, {
        key: '_defaultCompareFunction',
        value: function _defaultCompareFunction(a, b) {
            return b > a;
        }
    }]);

    return Quicksort;
}();

module.exports = Quicksort;

/***/ }),
/* 26 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


// Singular Value Decomposition implementation
function SVD() {};

/* Below singular value decomposition (svd) code including hypot function is adopted from https://github.com/dragonfly-ai/JamaJS
   Some changes are applied to make the code compatible with the fcose code and to make it independent from Jama.
   Input matrix is changed to a 2D array instead of Jama matrix. Matrix dimensions are taken according to 2D array instead of using Jama functions.
   An object that includes singular value components is created for return. 
   The types of input parameters of the hypot function are removed. 
   let is used instead of var for the variable initialization.
*/
/*
                               Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent
      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for reasonable and customary use in describing the
      origin of the Work and reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Additional Liability. While redistributing
      the Work or Derivative Works thereof, You may choose to offer,
      and charge a fee for, acceptance of support, warranty, indemnity,
      or other liability obligations and/or rights consistent with this
      License. However, in accepting such obligations, You may act only
      on Your own behalf and on Your sole responsibility, not on behalf
      of any other Contributor, and only if You agree to indemnify,
      defend, and hold each Contributor harmless for any liability
      incurred by, or claims asserted against, such Contributor by reason
      of your accepting any such warranty or additional liability.

   END OF TERMS AND CONDITIONS

   APPENDIX: How to apply the Apache License to your work.

      To apply the Apache License to your work, attach the following
      boilerplate notice, with the fields enclosed by brackets "{}"
      replaced with your own identifying information. (Don't include
      the brackets!)  The text should be enclosed in the appropriate
      comment syntax for the file format. We also recommend that a
      file or class name and description of purpose be included on the
      same "printed page" as the copyright notice for easier
      identification within third-party archives.

   Copyright {yyyy} {name of copyright owner}

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

SVD.svd = function (A) {
  this.U = null;
  this.V = null;
  this.s = null;
  this.m = 0;
  this.n = 0;
  this.m = A.length;
  this.n = A[0].length;
  var nu = Math.min(this.m, this.n);
  this.s = function (s) {
    var a = [];
    while (s-- > 0) {
      a.push(0);
    }return a;
  }(Math.min(this.m + 1, this.n));
  this.U = function (dims) {
    var allocate = function allocate(dims) {
      if (dims.length == 0) {
        return 0;
      } else {
        var array = [];
        for (var i = 0; i < dims[0]; i++) {
          array.push(allocate(dims.slice(1)));
        }
        return array;
      }
    };
    return allocate(dims);
  }([this.m, nu]);
  this.V = function (dims) {
    var allocate = function allocate(dims) {
      if (dims.length == 0) {
        return 0;
      } else {
        var array = [];
        for (var i = 0; i < dims[0]; i++) {
          array.push(allocate(dims.slice(1)));
        }
        return array;
      }
    };
    return allocate(dims);
  }([this.n, this.n]);
  var e = function (s) {
    var a = [];
    while (s-- > 0) {
      a.push(0);
    }return a;
  }(this.n);
  var work = function (s) {
    var a = [];
    while (s-- > 0) {
      a.push(0);
    }return a;
  }(this.m);
  var wantu = true;
  var wantv = true;
  var nct = Math.min(this.m - 1, this.n);
  var nrt = Math.max(0, Math.min(this.n - 2, this.m));
  for (var k = 0; k < Math.max(nct, nrt); k++) {
    if (k < nct) {
      this.s[k] = 0;
      for (var i = k; i < this.m; i++) {
        this.s[k] = SVD.hypot(this.s[k], A[i][k]);
      }
      ;
      if (this.s[k] !== 0.0) {
        if (A[k][k] < 0.0) {
          this.s[k] = -this.s[k];
        }
        for (var _i = k; _i < this.m; _i++) {
          A[_i][k] /= this.s[k];
        }
        ;
        A[k][k] += 1.0;
      }
      this.s[k] = -this.s[k];
    }
    for (var j = k + 1; j < this.n; j++) {
      if (function (lhs, rhs) {
        return lhs && rhs;
      }(k < nct, this.s[k] !== 0.0)) {
        var t = 0;
        for (var _i2 = k; _i2 < this.m; _i2++) {
          t += A[_i2][k] * A[_i2][j];
        }
        ;
        t = -t / A[k][k];
        for (var _i3 = k; _i3 < this.m; _i3++) {
          A[_i3][j] += t * A[_i3][k];
        }
        ;
      }
      e[j] = A[k][j];
    }
    ;
    if (function (lhs, rhs) {
      return lhs && rhs;
    }(wantu, k < nct)) {
      for (var _i4 = k; _i4 < this.m; _i4++) {
        this.U[_i4][k] = A[_i4][k];
      }
      ;
    }
    if (k < nrt) {
      e[k] = 0;
      for (var _i5 = k + 1; _i5 < this.n; _i5++) {
        e[k] = SVD.hypot(e[k], e[_i5]);
      }
      ;
      if (e[k] !== 0.0) {
        if (e[k + 1] < 0.0) {
          e[k] = -e[k];
        }
        for (var _i6 = k + 1; _i6 < this.n; _i6++) {
          e[_i6] /= e[k];
        }
        ;
        e[k + 1] += 1.0;
      }
      e[k] = -e[k];
      if (function (lhs, rhs) {
        return lhs && rhs;
      }(k + 1 < this.m, e[k] !== 0.0)) {
        for (var _i7 = k + 1; _i7 < this.m; _i7++) {
          work[_i7] = 0.0;
        }
        ;
        for (var _j = k + 1; _j < this.n; _j++) {
          for (var _i8 = k + 1; _i8 < this.m; _i8++) {
            work[_i8] += e[_j] * A[_i8][_j];
          }
          ;
        }
        ;
        for (var _j2 = k + 1; _j2 < this.n; _j2++) {
          var _t = -e[_j2] / e[k + 1];
          for (var _i9 = k + 1; _i9 < this.m; _i9++) {
            A[_i9][_j2] += _t * work[_i9];
          }
          ;
        }
        ;
      }
      if (wantv) {
        for (var _i10 = k + 1; _i10 < this.n; _i10++) {
          this.V[_i10][k] = e[_i10];
        };
      }
    }
  };
  var p = Math.min(this.n, this.m + 1);
  if (nct < this.n) {
    this.s[nct] = A[nct][nct];
  }
  if (this.m < p) {
    this.s[p - 1] = 0.0;
  }
  if (nrt + 1 < p) {
    e[nrt] = A[nrt][p - 1];
  }
  e[p - 1] = 0.0;
  if (wantu) {
    for (var _j3 = nct; _j3 < nu; _j3++) {
      for (var _i11 = 0; _i11 < this.m; _i11++) {
        this.U[_i11][_j3] = 0.0;
      }
      ;
      this.U[_j3][_j3] = 1.0;
    };
    for (var _k = nct - 1; _k >= 0; _k--) {
      if (this.s[_k] !== 0.0) {
        for (var _j4 = _k + 1; _j4 < nu; _j4++) {
          var _t2 = 0;
          for (var _i12 = _k; _i12 < this.m; _i12++) {
            _t2 += this.U[_i12][_k] * this.U[_i12][_j4];
          };
          _t2 = -_t2 / this.U[_k][_k];
          for (var _i13 = _k; _i13 < this.m; _i13++) {
            this.U[_i13][_j4] += _t2 * this.U[_i13][_k];
          };
        };
        for (var _i14 = _k; _i14 < this.m; _i14++) {
          this.U[_i14][_k] = -this.U[_i14][_k];
        };
        this.U[_k][_k] = 1.0 + this.U[_k][_k];
        for (var _i15 = 0; _i15 < _k - 1; _i15++) {
          this.U[_i15][_k] = 0.0;
        };
      } else {
        for (var _i16 = 0; _i16 < this.m; _i16++) {
          this.U[_i16][_k] = 0.0;
        };
        this.U[_k][_k] = 1.0;
      }
    };
  }
  if (wantv) {
    for (var _k2 = this.n - 1; _k2 >= 0; _k2--) {
      if (function (lhs, rhs) {
        return lhs && rhs;
      }(_k2 < nrt, e[_k2] !== 0.0)) {
        for (var _j5 = _k2 + 1; _j5 < nu; _j5++) {
          var _t3 = 0;
          for (var _i17 = _k2 + 1; _i17 < this.n; _i17++) {
            _t3 += this.V[_i17][_k2] * this.V[_i17][_j5];
          };
          _t3 = -_t3 / this.V[_k2 + 1][_k2];
          for (var _i18 = _k2 + 1; _i18 < this.n; _i18++) {
            this.V[_i18][_j5] += _t3 * this.V[_i18][_k2];
          };
        };
      }
      for (var _i19 = 0; _i19 < this.n; _i19++) {
        this.V[_i19][_k2] = 0.0;
      };
      this.V[_k2][_k2] = 1.0;
    };
  }
  var pp = p - 1;
  var iter = 0;
  var eps = Math.pow(2.0, -52.0);
  var tiny = Math.pow(2.0, -966.0);
  while (p > 0) {
    var _k3 = void 0;
    var kase = void 0;
    for (_k3 = p - 2; _k3 >= -1; _k3--) {
      if (_k3 === -1) {
        break;
      }
      if (Math.abs(e[_k3]) <= tiny + eps * (Math.abs(this.s[_k3]) + Math.abs(this.s[_k3 + 1]))) {
        e[_k3] = 0.0;
        break;
      }
    };
    if (_k3 === p - 2) {
      kase = 4;
    } else {
      var ks = void 0;
      for (ks = p - 1; ks >= _k3; ks--) {
        if (ks === _k3) {
          break;
        }
        var _t4 = (ks !== p ? Math.abs(e[ks]) : 0.0) + (ks !== _k3 + 1 ? Math.abs(e[ks - 1]) : 0.0);
        if (Math.abs(this.s[ks]) <= tiny + eps * _t4) {
          this.s[ks] = 0.0;
          break;
        }
      };
      if (ks === _k3) {
        kase = 3;
      } else if (ks === p - 1) {
        kase = 1;
      } else {
        kase = 2;
        _k3 = ks;
      }
    }
    _k3++;
    switch (kase) {
      case 1:
        {
          var f = e[p - 2];
          e[p - 2] = 0.0;
          for (var _j6 = p - 2; _j6 >= _k3; _j6--) {
            var _t5 = SVD.hypot(this.s[_j6], f);
            var cs = this.s[_j6] / _t5;
            var sn = f / _t5;
            this.s[_j6] = _t5;
            if (_j6 !== _k3) {
              f = -sn * e[_j6 - 1];
              e[_j6 - 1] = cs * e[_j6 - 1];
            }
            if (wantv) {
              for (var _i20 = 0; _i20 < this.n; _i20++) {
                _t5 = cs * this.V[_i20][_j6] + sn * this.V[_i20][p - 1];
                this.V[_i20][p - 1] = -sn * this.V[_i20][_j6] + cs * this.V[_i20][p - 1];
                this.V[_i20][_j6] = _t5;
              };
            }
          };
        };
        break;
      case 2:
        {
          var _f = e[_k3 - 1];
          e[_k3 - 1] = 0.0;
          for (var _j7 = _k3; _j7 < p; _j7++) {
            var _t6 = SVD.hypot(this.s[_j7], _f);
            var _cs = this.s[_j7] / _t6;
            var _sn = _f / _t6;
            this.s[_j7] = _t6;
            _f = -_sn * e[_j7];
            e[_j7] = _cs * e[_j7];
            if (wantu) {
              for (var _i21 = 0; _i21 < this.m; _i21++) {
                _t6 = _cs * this.U[_i21][_j7] + _sn * this.U[_i21][_k3 - 1];
                this.U[_i21][_k3 - 1] = -_sn * this.U[_i21][_j7] + _cs * this.U[_i21][_k3 - 1];
                this.U[_i21][_j7] = _t6;
              };
            }
          };
        };
        break;
      case 3:
        {
          var scale = Math.max(Math.max(Math.max(Math.max(Math.abs(this.s[p - 1]), Math.abs(this.s[p - 2])), Math.abs(e[p - 2])), Math.abs(this.s[_k3])), Math.abs(e[_k3]));
          var sp = this.s[p - 1] / scale;
          var spm1 = this.s[p - 2] / scale;
          var epm1 = e[p - 2] / scale;
          var sk = this.s[_k3] / scale;
          var ek = e[_k3] / scale;
          var b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
          var c = sp * epm1 * (sp * epm1);
          var shift = 0.0;
          if (function (lhs, rhs) {
            return lhs || rhs;
          }(b !== 0.0, c !== 0.0)) {
            shift = Math.sqrt(b * b + c);
            if (b < 0.0) {
              shift = -shift;
            }
            shift = c / (b + shift);
          }
          var _f2 = (sk + sp) * (sk - sp) + shift;
          var g = sk * ek;
          for (var _j8 = _k3; _j8 < p - 1; _j8++) {
            var _t7 = SVD.hypot(_f2, g);
            var _cs2 = _f2 / _t7;
            var _sn2 = g / _t7;
            if (_j8 !== _k3) {
              e[_j8 - 1] = _t7;
            }
            _f2 = _cs2 * this.s[_j8] + _sn2 * e[_j8];
            e[_j8] = _cs2 * e[_j8] - _sn2 * this.s[_j8];
            g = _sn2 * this.s[_j8 + 1];
            this.s[_j8 + 1] = _cs2 * this.s[_j8 + 1];
            if (wantv) {
              for (var _i22 = 0; _i22 < this.n; _i22++) {
                _t7 = _cs2 * this.V[_i22][_j8] + _sn2 * this.V[_i22][_j8 + 1];
                this.V[_i22][_j8 + 1] = -_sn2 * this.V[_i22][_j8] + _cs2 * this.V[_i22][_j8 + 1];
                this.V[_i22][_j8] = _t7;
              };
            }
            _t7 = SVD.hypot(_f2, g);
            _cs2 = _f2 / _t7;
            _sn2 = g / _t7;
            this.s[_j8] = _t7;
            _f2 = _cs2 * e[_j8] + _sn2 * this.s[_j8 + 1];
            this.s[_j8 + 1] = -_sn2 * e[_j8] + _cs2 * this.s[_j8 + 1];
            g = _sn2 * e[_j8 + 1];
            e[_j8 + 1] = _cs2 * e[_j8 + 1];
            if (wantu && _j8 < this.m - 1) {
              for (var _i23 = 0; _i23 < this.m; _i23++) {
                _t7 = _cs2 * this.U[_i23][_j8] + _sn2 * this.U[_i23][_j8 + 1];
                this.U[_i23][_j8 + 1] = -_sn2 * this.U[_i23][_j8] + _cs2 * this.U[_i23][_j8 + 1];
                this.U[_i23][_j8] = _t7;
              };
            }
          };
          e[p - 2] = _f2;
          iter = iter + 1;
        };
        break;
      case 4:
        {
          if (this.s[_k3] <= 0.0) {
            this.s[_k3] = this.s[_k3] < 0.0 ? -this.s[_k3] : 0.0;
            if (wantv) {
              for (var _i24 = 0; _i24 <= pp; _i24++) {
                this.V[_i24][_k3] = -this.V[_i24][_k3];
              };
            }
          }
          while (_k3 < pp) {
            if (this.s[_k3] >= this.s[_k3 + 1]) {
              break;
            }
            var _t8 = this.s[_k3];
            this.s[_k3] = this.s[_k3 + 1];
            this.s[_k3 + 1] = _t8;
            if (wantv && _k3 < this.n - 1) {
              for (var _i25 = 0; _i25 < this.n; _i25++) {
                _t8 = this.V[_i25][_k3 + 1];
                this.V[_i25][_k3 + 1] = this.V[_i25][_k3];
                this.V[_i25][_k3] = _t8;
              };
            }
            if (wantu && _k3 < this.m - 1) {
              for (var _i26 = 0; _i26 < this.m; _i26++) {
                _t8 = this.U[_i26][_k3 + 1];
                this.U[_i26][_k3 + 1] = this.U[_i26][_k3];
                this.U[_i26][_k3] = _t8;
              };
            }
            _k3++;
          };
          iter = 0;
          p--;
        };
        break;
    }
  };
  var result = { U: this.U, V: this.V, S: this.s };
  return result;
};

// sqrt(a^2 + b^2) without under/overflow.
SVD.hypot = function (a, b) {
  var r = void 0;
  if (Math.abs(a) > Math.abs(b)) {
    r = b / a;
    r = Math.abs(a) * Math.sqrt(1 + r * r);
  } else if (b != 0) {
    r = a / b;
    r = Math.abs(b) * Math.sqrt(1 + r * r);
  } else {
    r = 0.0;
  }
  return r;
};

module.exports = SVD;

/***/ }),
/* 27 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

/**
 *   Needleman-Wunsch algorithm is an procedure to compute the optimal global alignment of two string
 *   sequences by S.B.Needleman and C.D.Wunsch (1970).
 *
 *   Aside from the inputs, you can assign the scores for,
 *   - Match: The two characters at the current index are same.
 *   - Mismatch: The two characters at the current index are different.
 *   - Insertion/Deletion(gaps): The best alignment involves one letter aligning to a gap in the other string.
 */

var NeedlemanWunsch = function () {
    function NeedlemanWunsch(sequence1, sequence2) {
        var match_score = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : 1;
        var mismatch_penalty = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : -1;
        var gap_penalty = arguments.length > 4 && arguments[4] !== undefined ? arguments[4] : -1;

        _classCallCheck(this, NeedlemanWunsch);

        this.sequence1 = sequence1;
        this.sequence2 = sequence2;
        this.match_score = match_score;
        this.mismatch_penalty = mismatch_penalty;
        this.gap_penalty = gap_penalty;

        // Just the remove redundancy
        this.iMax = sequence1.length + 1;
        this.jMax = sequence2.length + 1;

        // Grid matrix of scores
        this.grid = new Array(this.iMax);
        for (var i = 0; i < this.iMax; i++) {
            this.grid[i] = new Array(this.jMax);

            for (var j = 0; j < this.jMax; j++) {
                this.grid[i][j] = 0;
            }
        }

        // Traceback matrix (2D array, each cell is an array of boolean values for [`Diag`, `Up`, `Left`] positions)
        this.tracebackGrid = new Array(this.iMax);
        for (var _i = 0; _i < this.iMax; _i++) {
            this.tracebackGrid[_i] = new Array(this.jMax);

            for (var _j = 0; _j < this.jMax; _j++) {
                this.tracebackGrid[_i][_j] = [null, null, null];
            }
        }

        // The aligned sequences (return multiple possibilities)
        this.alignments = [];

        // Final alignment score
        this.score = -1;

        // Calculate scores and tracebacks
        this.computeGrids();
    }

    _createClass(NeedlemanWunsch, [{
        key: "getScore",
        value: function getScore() {
            return this.score;
        }
    }, {
        key: "getAlignments",
        value: function getAlignments() {
            return this.alignments;
        }

        // Main dynamic programming procedure

    }, {
        key: "computeGrids",
        value: function computeGrids() {
            // Fill in the first row
            for (var j = 1; j < this.jMax; j++) {
                this.grid[0][j] = this.grid[0][j - 1] + this.gap_penalty;
                this.tracebackGrid[0][j] = [false, false, true];
            }

            // Fill in the first column
            for (var i = 1; i < this.iMax; i++) {
                this.grid[i][0] = this.grid[i - 1][0] + this.gap_penalty;
                this.tracebackGrid[i][0] = [false, true, false];
            }

            // Fill the rest of the grid
            for (var _i2 = 1; _i2 < this.iMax; _i2++) {
                for (var _j2 = 1; _j2 < this.jMax; _j2++) {
                    // Find the max score(s) among [`Diag`, `Up`, `Left`]
                    var diag = void 0;
                    if (this.sequence1[_i2 - 1] === this.sequence2[_j2 - 1]) diag = this.grid[_i2 - 1][_j2 - 1] + this.match_score;else diag = this.grid[_i2 - 1][_j2 - 1] + this.mismatch_penalty;

                    var up = this.grid[_i2 - 1][_j2] + this.gap_penalty;
                    var left = this.grid[_i2][_j2 - 1] + this.gap_penalty;

                    // If there exists multiple max values, capture them for multiple paths
                    var maxOf = [diag, up, left];
                    var indices = this.arrayAllMaxIndexes(maxOf);

                    // Update Grids
                    this.grid[_i2][_j2] = maxOf[indices[0]];
                    this.tracebackGrid[_i2][_j2] = [indices.includes(0), indices.includes(1), indices.includes(2)];
                }
            }

            // Update alignment score
            this.score = this.grid[this.iMax - 1][this.jMax - 1];
        }

        // Gets all possible valid sequence combinations

    }, {
        key: "alignmentTraceback",
        value: function alignmentTraceback() {
            var inProcessAlignments = [];

            inProcessAlignments.push({ pos: [this.sequence1.length, this.sequence2.length],
                seq1: "",
                seq2: ""
            });

            while (inProcessAlignments[0]) {
                var current = inProcessAlignments[0];
                var directions = this.tracebackGrid[current.pos[0]][current.pos[1]];

                if (directions[0]) {
                    inProcessAlignments.push({ pos: [current.pos[0] - 1, current.pos[1] - 1],
                        seq1: this.sequence1[current.pos[0] - 1] + current.seq1,
                        seq2: this.sequence2[current.pos[1] - 1] + current.seq2
                    });
                }
                if (directions[1]) {
                    inProcessAlignments.push({ pos: [current.pos[0] - 1, current.pos[1]],
                        seq1: this.sequence1[current.pos[0] - 1] + current.seq1,
                        seq2: '-' + current.seq2
                    });
                }
                if (directions[2]) {
                    inProcessAlignments.push({ pos: [current.pos[0], current.pos[1] - 1],
                        seq1: '-' + current.seq1,
                        seq2: this.sequence2[current.pos[1] - 1] + current.seq2
                    });
                }

                if (current.pos[0] === 0 && current.pos[1] === 0) this.alignments.push({ sequence1: current.seq1,
                    sequence2: current.seq2
                });

                inProcessAlignments.shift();
            }

            return this.alignments;
        }

        // Helper Functions

    }, {
        key: "getAllIndexes",
        value: function getAllIndexes(arr, val) {
            var indexes = [],
                i = -1;
            while ((i = arr.indexOf(val, i + 1)) !== -1) {
                indexes.push(i);
            }
            return indexes;
        }
    }, {
        key: "arrayAllMaxIndexes",
        value: function arrayAllMaxIndexes(array) {
            return this.getAllIndexes(array, Math.max.apply(null, array));
        }
    }]);

    return NeedlemanWunsch;
}();

module.exports = NeedlemanWunsch;

/***/ }),
/* 28 */
/***/ (function(module, exports, __nested_webpack_require_145732__) {

"use strict";


var layoutBase = function layoutBase() {
  return;
};

layoutBase.FDLayout = __nested_webpack_require_145732__(18);
layoutBase.FDLayoutConstants = __nested_webpack_require_145732__(4);
layoutBase.FDLayoutEdge = __nested_webpack_require_145732__(19);
layoutBase.FDLayoutNode = __nested_webpack_require_145732__(20);
layoutBase.DimensionD = __nested_webpack_require_145732__(21);
layoutBase.HashMap = __nested_webpack_require_145732__(22);
layoutBase.HashSet = __nested_webpack_require_145732__(23);
layoutBase.IGeometry = __nested_webpack_require_145732__(8);
layoutBase.IMath = __nested_webpack_require_145732__(9);
layoutBase.Integer = __nested_webpack_require_145732__(10);
layoutBase.Point = __nested_webpack_require_145732__(12);
layoutBase.PointD = __nested_webpack_require_145732__(5);
layoutBase.RandomSeed = __nested_webpack_require_145732__(16);
layoutBase.RectangleD = __nested_webpack_require_145732__(13);
layoutBase.Transform = __nested_webpack_require_145732__(17);
layoutBase.UniqueIDGeneretor = __nested_webpack_require_145732__(14);
layoutBase.Quicksort = __nested_webpack_require_145732__(25);
layoutBase.LinkedList = __nested_webpack_require_145732__(11);
layoutBase.LGraphObject = __nested_webpack_require_145732__(2);
layoutBase.LGraph = __nested_webpack_require_145732__(6);
layoutBase.LEdge = __nested_webpack_require_145732__(1);
layoutBase.LGraphManager = __nested_webpack_require_145732__(7);
layoutBase.LNode = __nested_webpack_require_145732__(3);
layoutBase.Layout = __nested_webpack_require_145732__(15);
layoutBase.LayoutConstants = __nested_webpack_require_145732__(0);
layoutBase.NeedlemanWunsch = __nested_webpack_require_145732__(27);
layoutBase.Matrix = __nested_webpack_require_145732__(24);
layoutBase.SVD = __nested_webpack_require_145732__(26);

module.exports = layoutBase;

/***/ }),
/* 29 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


function Emitter() {
  this.listeners = [];
}

var p = Emitter.prototype;

p.addListener = function (event, callback) {
  this.listeners.push({
    event: event,
    callback: callback
  });
};

p.removeListener = function (event, callback) {
  for (var i = this.listeners.length; i >= 0; i--) {
    var l = this.listeners[i];

    if (l.event === event && l.callback === callback) {
      this.listeners.splice(i, 1);
    }
  }
};

p.emit = function (event, data) {
  for (var i = 0; i < this.listeners.length; i++) {
    var l = this.listeners[i];

    if (event === l.event) {
      l.callback(data);
    }
  }
};

module.exports = Emitter;

/***/ })
/******/ ]);
});

/***/ }),

/***/ 72425:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   diagram: () => (/* binding */ diagram)
/* harmony export */ });
/* harmony import */ var _chunk_H2D2JQ3I_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(62072);
/* harmony import */ var _chunk_C3MQ5ANM_mjs__WEBPACK_IMPORTED_MODULE_1__ = __webpack_require__(39769);
/* harmony import */ var _chunk_4BMEZGHF_mjs__WEBPACK_IMPORTED_MODULE_2__ = __webpack_require__(80910);
/* harmony import */ var _chunk_XZIHB7SX_mjs__WEBPACK_IMPORTED_MODULE_3__ = __webpack_require__(82002);
/* harmony import */ var _chunk_O4NI6UNU_mjs__WEBPACK_IMPORTED_MODULE_4__ = __webpack_require__(42626);
/* harmony import */ var _chunk_7B677QYD_mjs__WEBPACK_IMPORTED_MODULE_5__ = __webpack_require__(30103);
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__ = __webpack_require__(86906);
/* harmony import */ var _mermaid_js_parser__WEBPACK_IMPORTED_MODULE_7__ = __webpack_require__(13197);
/* harmony import */ var cytoscape__WEBPACK_IMPORTED_MODULE_8__ = __webpack_require__(23207);
/* harmony import */ var cytoscape_fcose__WEBPACK_IMPORTED_MODULE_9__ = __webpack_require__(29287);
/* harmony import */ var cytoscape_fcose__WEBPACK_IMPORTED_MODULE_9___default = /*#__PURE__*/__webpack_require__.n(cytoscape_fcose__WEBPACK_IMPORTED_MODULE_9__);
/* harmony import */ var d3__WEBPACK_IMPORTED_MODULE_10__ = __webpack_require__(83619);








// src/diagrams/architecture/architectureParser.ts


// src/diagrams/architecture/architectureTypes.ts
var ArchitectureDirectionName = {
  L: "left",
  R: "right",
  T: "top",
  B: "bottom"
};
var ArchitectureDirectionArrow = {
  L: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((scale) => `${scale},${scale / 2} 0,${scale} 0,0`, "L"),
  R: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((scale) => `0,${scale / 2} ${scale},0 ${scale},${scale}`, "R"),
  T: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((scale) => `0,0 ${scale},0 ${scale / 2},${scale}`, "T"),
  B: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((scale) => `${scale / 2},0 ${scale},${scale} 0,${scale}`, "B")
};
var ArchitectureDirectionArrowShift = {
  L: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((orig, arrowSize) => orig - arrowSize + 2, "L"),
  R: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((orig, _arrowSize) => orig - 2, "R"),
  T: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((orig, arrowSize) => orig - arrowSize + 2, "T"),
  B: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((orig, _arrowSize) => orig - 2, "B")
};
var getOppositeArchitectureDirection = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function(x) {
  if (isArchitectureDirectionX(x)) {
    return x === "L" ? "R" : "L";
  } else {
    return x === "T" ? "B" : "T";
  }
}, "getOppositeArchitectureDirection");
var isArchitectureDirection = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function(x) {
  const temp = x;
  return temp === "L" || temp === "R" || temp === "T" || temp === "B";
}, "isArchitectureDirection");
var isArchitectureDirectionX = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function(x) {
  const temp = x;
  return temp === "L" || temp === "R";
}, "isArchitectureDirectionX");
var isArchitectureDirectionY = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function(x) {
  const temp = x;
  return temp === "T" || temp === "B";
}, "isArchitectureDirectionY");
var isArchitectureDirectionXY = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function(a, b) {
  const aX_bY = isArchitectureDirectionX(a) && isArchitectureDirectionY(b);
  const aY_bX = isArchitectureDirectionY(a) && isArchitectureDirectionX(b);
  return aX_bY || aY_bX;
}, "isArchitectureDirectionXY");
var isArchitecturePairXY = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function(pair) {
  const lhs = pair[0];
  const rhs = pair[1];
  const aX_bY = isArchitectureDirectionX(lhs) && isArchitectureDirectionY(rhs);
  const aY_bX = isArchitectureDirectionY(lhs) && isArchitectureDirectionX(rhs);
  return aX_bY || aY_bX;
}, "isArchitecturePairXY");
var isValidArchitectureDirectionPair = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function(x) {
  return x !== "LL" && x !== "RR" && x !== "TT" && x !== "BB";
}, "isValidArchitectureDirectionPair");
var getArchitectureDirectionPair = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function(sourceDir, targetDir) {
  const pair = `${sourceDir}${targetDir}`;
  return isValidArchitectureDirectionPair(pair) ? pair : void 0;
}, "getArchitectureDirectionPair");
var shiftPositionByArchitectureDirectionPair = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function([x, y], pair) {
  const lhs = pair[0];
  const rhs = pair[1];
  if (isArchitectureDirectionX(lhs)) {
    if (isArchitectureDirectionY(rhs)) {
      return [x + (lhs === "L" ? -1 : 1), y + (rhs === "T" ? 1 : -1)];
    } else {
      return [x + (lhs === "L" ? -1 : 1), y];
    }
  } else {
    if (isArchitectureDirectionX(rhs)) {
      return [x + (rhs === "L" ? 1 : -1), y + (lhs === "T" ? 1 : -1)];
    } else {
      return [x, y + (lhs === "T" ? 1 : -1)];
    }
  }
}, "shiftPositionByArchitectureDirectionPair");
var getArchitectureDirectionXYFactors = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function(pair) {
  if (pair === "LT" || pair === "TL") {
    return [1, 1];
  } else if (pair === "BL" || pair === "LB") {
    return [1, -1];
  } else if (pair === "BR" || pair === "RB") {
    return [-1, -1];
  } else {
    return [-1, 1];
  }
}, "getArchitectureDirectionXYFactors");
var getArchitectureDirectionAlignment = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function(a, b) {
  if (isArchitectureDirectionXY(a, b)) {
    return "bend";
  } else if (isArchitectureDirectionX(a)) {
    return "horizontal";
  }
  return "vertical";
}, "getArchitectureDirectionAlignment");
var isArchitectureService = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function(x) {
  const temp = x;
  return temp.type === "service";
}, "isArchitectureService");
var isArchitectureJunction = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function(x) {
  const temp = x;
  return temp.type === "junction";
}, "isArchitectureJunction");
var edgeData = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((edge) => {
  return edge.data();
}, "edgeData");
var nodeData = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((node) => {
  return node.data();
}, "nodeData");

// src/diagrams/architecture/architectureDb.ts
var DEFAULT_ARCHITECTURE_CONFIG = _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .defaultConfig_default */ .vZ.architecture;
var state = new _chunk_XZIHB7SX_mjs__WEBPACK_IMPORTED_MODULE_3__/* .ImperativeState */ .A(() => ({
  nodes: {},
  groups: {},
  edges: [],
  registeredIds: {},
  config: DEFAULT_ARCHITECTURE_CONFIG,
  dataStructures: void 0,
  elements: {}
}));
var clear2 = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(() => {
  state.reset();
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .clear */ .ZH)();
}, "clear");
var addService = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function({
  id,
  icon,
  in: parent,
  title,
  iconText
}) {
  if (state.records.registeredIds[id] !== void 0) {
    throw new Error(
      `The service id [${id}] is already in use by another ${state.records.registeredIds[id]}`
    );
  }
  if (parent !== void 0) {
    if (id === parent) {
      throw new Error(`The service [${id}] cannot be placed within itself`);
    }
    if (state.records.registeredIds[parent] === void 0) {
      throw new Error(
        `The service [${id}]'s parent does not exist. Please make sure the parent is created before this service`
      );
    }
    if (state.records.registeredIds[parent] === "node") {
      throw new Error(`The service [${id}]'s parent is not a group`);
    }
  }
  state.records.registeredIds[id] = "node";
  state.records.nodes[id] = {
    id,
    type: "service",
    icon,
    iconText,
    title,
    edges: [],
    in: parent
  };
}, "addService");
var getServices = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(() => Object.values(state.records.nodes).filter(isArchitectureService), "getServices");
var addJunction = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function({ id, in: parent }) {
  state.records.registeredIds[id] = "node";
  state.records.nodes[id] = {
    id,
    type: "junction",
    edges: [],
    in: parent
  };
}, "addJunction");
var getJunctions = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(() => Object.values(state.records.nodes).filter(isArchitectureJunction), "getJunctions");
var getNodes = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(() => Object.values(state.records.nodes), "getNodes");
var getNode = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((id) => state.records.nodes[id], "getNode");
var addGroup = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function({ id, icon, in: parent, title }) {
  if (state.records.registeredIds[id] !== void 0) {
    throw new Error(
      `The group id [${id}] is already in use by another ${state.records.registeredIds[id]}`
    );
  }
  if (parent !== void 0) {
    if (id === parent) {
      throw new Error(`The group [${id}] cannot be placed within itself`);
    }
    if (state.records.registeredIds[parent] === void 0) {
      throw new Error(
        `The group [${id}]'s parent does not exist. Please make sure the parent is created before this group`
      );
    }
    if (state.records.registeredIds[parent] === "node") {
      throw new Error(`The group [${id}]'s parent is not a group`);
    }
  }
  state.records.registeredIds[id] = "group";
  state.records.groups[id] = {
    id,
    icon,
    title,
    in: parent
  };
}, "addGroup");
var getGroups = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(() => {
  return Object.values(state.records.groups);
}, "getGroups");
var addEdge = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function({
  lhsId,
  rhsId,
  lhsDir,
  rhsDir,
  lhsInto,
  rhsInto,
  lhsGroup,
  rhsGroup,
  title
}) {
  if (!isArchitectureDirection(lhsDir)) {
    throw new Error(
      `Invalid direction given for left hand side of edge ${lhsId}--${rhsId}. Expected (L,R,T,B) got ${lhsDir}`
    );
  }
  if (!isArchitectureDirection(rhsDir)) {
    throw new Error(
      `Invalid direction given for right hand side of edge ${lhsId}--${rhsId}. Expected (L,R,T,B) got ${rhsDir}`
    );
  }
  if (state.records.nodes[lhsId] === void 0 && state.records.groups[lhsId] === void 0) {
    throw new Error(
      `The left-hand id [${lhsId}] does not yet exist. Please create the service/group before declaring an edge to it.`
    );
  }
  if (state.records.nodes[rhsId] === void 0 && state.records.groups[lhsId] === void 0) {
    throw new Error(
      `The right-hand id [${rhsId}] does not yet exist. Please create the service/group before declaring an edge to it.`
    );
  }
  const lhsGroupId = state.records.nodes[lhsId].in;
  const rhsGroupId = state.records.nodes[rhsId].in;
  if (lhsGroup && lhsGroupId && rhsGroupId && lhsGroupId == rhsGroupId) {
    throw new Error(
      `The left-hand id [${lhsId}] is modified to traverse the group boundary, but the edge does not pass through two groups.`
    );
  }
  if (rhsGroup && lhsGroupId && rhsGroupId && lhsGroupId == rhsGroupId) {
    throw new Error(
      `The right-hand id [${rhsId}] is modified to traverse the group boundary, but the edge does not pass through two groups.`
    );
  }
  const edge = {
    lhsId,
    lhsDir,
    lhsInto,
    lhsGroup,
    rhsId,
    rhsDir,
    rhsInto,
    rhsGroup,
    title
  };
  state.records.edges.push(edge);
  if (state.records.nodes[lhsId] && state.records.nodes[rhsId]) {
    state.records.nodes[lhsId].edges.push(state.records.edges[state.records.edges.length - 1]);
    state.records.nodes[rhsId].edges.push(state.records.edges[state.records.edges.length - 1]);
  }
}, "addEdge");
var getEdges = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(() => state.records.edges, "getEdges");
var getDataStructures = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(() => {
  if (state.records.dataStructures === void 0) {
    const groupAlignments = {};
    const adjList = Object.entries(state.records.nodes).reduce((prevOuter, [id, service]) => {
      prevOuter[id] = service.edges.reduce((prevInner, edge) => {
        const lhsGroupId = getNode(edge.lhsId)?.in;
        const rhsGroupId = getNode(edge.rhsId)?.in;
        if (lhsGroupId && rhsGroupId && lhsGroupId !== rhsGroupId) {
          const alignment = getArchitectureDirectionAlignment(edge.lhsDir, edge.rhsDir);
          if (alignment !== "bend") {
            groupAlignments[lhsGroupId] ??= {};
            groupAlignments[lhsGroupId][rhsGroupId] = alignment;
            groupAlignments[rhsGroupId] ??= {};
            groupAlignments[rhsGroupId][lhsGroupId] = alignment;
          }
        }
        if (edge.lhsId === id) {
          const pair = getArchitectureDirectionPair(edge.lhsDir, edge.rhsDir);
          if (pair) {
            prevInner[pair] = edge.rhsId;
          }
        } else {
          const pair = getArchitectureDirectionPair(edge.rhsDir, edge.lhsDir);
          if (pair) {
            prevInner[pair] = edge.lhsId;
          }
        }
        return prevInner;
      }, {});
      return prevOuter;
    }, {});
    const firstId = Object.keys(adjList)[0];
    const visited = { [firstId]: 1 };
    const notVisited = Object.keys(adjList).reduce(
      (prev, id) => id === firstId ? prev : { ...prev, [id]: 1 },
      {}
    );
    const BFS = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((startingId) => {
      const spatialMap = { [startingId]: [0, 0] };
      const queue = [startingId];
      while (queue.length > 0) {
        const id = queue.shift();
        if (id) {
          visited[id] = 1;
          delete notVisited[id];
          const adj = adjList[id];
          const [posX, posY] = spatialMap[id];
          Object.entries(adj).forEach(([dir, rhsId]) => {
            if (!visited[rhsId]) {
              spatialMap[rhsId] = shiftPositionByArchitectureDirectionPair(
                [posX, posY],
                dir
              );
              queue.push(rhsId);
            }
          });
        }
      }
      return spatialMap;
    }, "BFS");
    const spatialMaps = [BFS(firstId)];
    while (Object.keys(notVisited).length > 0) {
      spatialMaps.push(BFS(Object.keys(notVisited)[0]));
    }
    state.records.dataStructures = {
      adjList,
      spatialMaps,
      groupAlignments
    };
  }
  return state.records.dataStructures;
}, "getDataStructures");
var setElementForId = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((id, element) => {
  state.records.elements[id] = element;
}, "setElementForId");
var getElementById = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((id) => state.records.elements[id], "getElementById");
var db = {
  clear: clear2,
  setDiagramTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .setDiagramTitle */ .g2,
  getDiagramTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .getDiagramTitle */ .Kr,
  setAccTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .setAccTitle */ .GN,
  getAccTitle: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .getAccTitle */ .eu,
  setAccDescription: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .setAccDescription */ .U$,
  getAccDescription: _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .getAccDescription */ .Mx,
  addService,
  getServices,
  addJunction,
  getJunctions,
  getNodes,
  getNode,
  addGroup,
  getGroups,
  addEdge,
  getEdges,
  setElementForId,
  getElementById,
  getDataStructures
};
function getConfigField(field) {
  const arch = (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .getConfig2 */ .nV)().architecture;
  if (arch?.[field]) {
    return arch[field];
  }
  return DEFAULT_ARCHITECTURE_CONFIG[field];
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(getConfigField, "getConfigField");

// src/diagrams/architecture/architectureParser.ts
var populateDb = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((ast, db2) => {
  (0,_chunk_4BMEZGHF_mjs__WEBPACK_IMPORTED_MODULE_2__/* .populateCommonDb */ .A)(ast, db2);
  ast.groups.map(db2.addGroup);
  ast.services.map((service) => db2.addService({ ...service, type: "service" }));
  ast.junctions.map((service) => db2.addJunction({ ...service, type: "junction" }));
  ast.edges.map(db2.addEdge);
}, "populateDb");
var parser = {
  parse: /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(async (input) => {
    const ast = await (0,_mermaid_js_parser__WEBPACK_IMPORTED_MODULE_7__/* .parse */ .Qc)("architecture", input);
    _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .log */ .cM.debug(ast);
    populateDb(ast, db);
  }, "parse")
};

// src/diagrams/architecture/architectureStyles.ts
var getStyles = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((options) => `
  .edge {
    stroke-width: ${options.archEdgeWidth};
    stroke: ${options.archEdgeColor};
    fill: none;
  }

  .arrow {
    fill: ${options.archEdgeArrowColor};
  }

  .node-bkg {
    fill: none;
    stroke: ${options.archGroupBorderColor};
    stroke-width: ${options.archGroupBorderWidth};
    stroke-dasharray: 8;
  }
  .node-icon-text {
    display: flex; 
    align-items: center;
  }
  
  .node-icon-text > div {
    color: #fff;
    margin: 1px;
    height: fit-content;
    text-align: center;
    overflow: hidden;
    display: -webkit-box;
    -webkit-box-orient: vertical;
  }
`, "getStyles");
var architectureStyles_default = getStyles;

// src/diagrams/architecture/architectureRenderer.ts




// src/diagrams/architecture/architectureIcons.ts
var wrapIcon = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((icon) => {
  return `<g><rect width="80" height="80" style="fill: #087ebf; stroke-width: 0px;"/>${icon}</g>`;
}, "wrapIcon");
var architectureIcons = {
  prefix: "mermaid-architecture",
  height: 80,
  width: 80,
  icons: {
    database: {
      body: wrapIcon(
        '<path id="b" data-name="4" d="m20,57.86c0,3.94,8.95,7.14,20,7.14s20-3.2,20-7.14" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><path id="c" data-name="3" d="m20,45.95c0,3.94,8.95,7.14,20,7.14s20-3.2,20-7.14" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><path id="d" data-name="2" d="m20,34.05c0,3.94,8.95,7.14,20,7.14s20-3.2,20-7.14" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><ellipse id="e" data-name="1" cx="40" cy="22.14" rx="20" ry="7.14" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><line x1="20" y1="57.86" x2="20" y2="22.14" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><line x1="60" y1="57.86" x2="60" y2="22.14" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/>'
      )
    },
    server: {
      body: wrapIcon(
        '<rect x="17.5" y="17.5" width="45" height="45" rx="2" ry="2" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><line x1="17.5" y1="32.5" x2="62.5" y2="32.5" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><line x1="17.5" y1="47.5" x2="62.5" y2="47.5" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><g><path d="m56.25,25c0,.27-.45.5-1,.5h-10.5c-.55,0-1-.23-1-.5s.45-.5,1-.5h10.5c.55,0,1,.23,1,.5Z" style="fill: #fff; stroke-width: 0px;"/><path d="m56.25,25c0,.27-.45.5-1,.5h-10.5c-.55,0-1-.23-1-.5s.45-.5,1-.5h10.5c.55,0,1,.23,1,.5Z" style="fill: none; stroke: #fff; stroke-miterlimit: 10;"/></g><g><path d="m56.25,40c0,.27-.45.5-1,.5h-10.5c-.55,0-1-.23-1-.5s.45-.5,1-.5h10.5c.55,0,1,.23,1,.5Z" style="fill: #fff; stroke-width: 0px;"/><path d="m56.25,40c0,.27-.45.5-1,.5h-10.5c-.55,0-1-.23-1-.5s.45-.5,1-.5h10.5c.55,0,1,.23,1,.5Z" style="fill: none; stroke: #fff; stroke-miterlimit: 10;"/></g><g><path d="m56.25,55c0,.27-.45.5-1,.5h-10.5c-.55,0-1-.23-1-.5s.45-.5,1-.5h10.5c.55,0,1,.23,1,.5Z" style="fill: #fff; stroke-width: 0px;"/><path d="m56.25,55c0,.27-.45.5-1,.5h-10.5c-.55,0-1-.23-1-.5s.45-.5,1-.5h10.5c.55,0,1,.23,1,.5Z" style="fill: none; stroke: #fff; stroke-miterlimit: 10;"/></g><g><circle cx="32.5" cy="25" r=".75" style="fill: #fff; stroke: #fff; stroke-miterlimit: 10;"/><circle cx="27.5" cy="25" r=".75" style="fill: #fff; stroke: #fff; stroke-miterlimit: 10;"/><circle cx="22.5" cy="25" r=".75" style="fill: #fff; stroke: #fff; stroke-miterlimit: 10;"/></g><g><circle cx="32.5" cy="40" r=".75" style="fill: #fff; stroke: #fff; stroke-miterlimit: 10;"/><circle cx="27.5" cy="40" r=".75" style="fill: #fff; stroke: #fff; stroke-miterlimit: 10;"/><circle cx="22.5" cy="40" r=".75" style="fill: #fff; stroke: #fff; stroke-miterlimit: 10;"/></g><g><circle cx="32.5" cy="55" r=".75" style="fill: #fff; stroke: #fff; stroke-miterlimit: 10;"/><circle cx="27.5" cy="55" r=".75" style="fill: #fff; stroke: #fff; stroke-miterlimit: 10;"/><circle cx="22.5" cy="55" r=".75" style="fill: #fff; stroke: #fff; stroke-miterlimit: 10;"/></g>'
      )
    },
    disk: {
      body: wrapIcon(
        '<rect x="20" y="15" width="40" height="50" rx="1" ry="1" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><ellipse cx="24" cy="19.17" rx=".8" ry=".83" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><ellipse cx="56" cy="19.17" rx=".8" ry=".83" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><ellipse cx="24" cy="60.83" rx=".8" ry=".83" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><ellipse cx="56" cy="60.83" rx=".8" ry=".83" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><ellipse cx="40" cy="33.75" rx="14" ry="14.58" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><ellipse cx="40" cy="33.75" rx="4" ry="4.17" style="fill: #fff; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><path d="m37.51,42.52l-4.83,13.22c-.26.71-1.1,1.02-1.76.64l-4.18-2.42c-.66-.38-.81-1.26-.33-1.84l9.01-10.8c.88-1.05,2.56-.08,2.09,1.2Z" style="fill: #fff; stroke-width: 0px;"/>'
      )
    },
    internet: {
      body: wrapIcon(
        '<circle cx="40" cy="40" r="22.5" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><line x1="40" y1="17.5" x2="40" y2="62.5" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><line x1="17.5" y1="40" x2="62.5" y2="40" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><path d="m39.99,17.51c-15.28,11.1-15.28,33.88,0,44.98" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><path d="m40.01,17.51c15.28,11.1,15.28,33.88,0,44.98" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><line x1="19.75" y1="30.1" x2="60.25" y2="30.1" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/><line x1="19.75" y1="49.9" x2="60.25" y2="49.9" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/>'
      )
    },
    cloud: {
      body: wrapIcon(
        '<path d="m65,47.5c0,2.76-2.24,5-5,5H20c-2.76,0-5-2.24-5-5,0-1.87,1.03-3.51,2.56-4.36-.04-.21-.06-.42-.06-.64,0-2.6,2.48-4.74,5.65-4.97,1.65-4.51,6.34-7.76,11.85-7.76.86,0,1.69.08,2.5.23,2.09-1.57,4.69-2.5,7.5-2.5,6.1,0,11.19,4.38,12.28,10.17,2.14.56,3.72,2.51,3.72,4.83,0,.03,0,.07-.01.1,2.29.46,4.01,2.48,4.01,4.9Z" style="fill: none; stroke: #fff; stroke-miterlimit: 10; stroke-width: 2px;"/>'
      )
    },
    unknown: _chunk_H2D2JQ3I_mjs__WEBPACK_IMPORTED_MODULE_0__/* .unknownIcon */ .cN,
    blank: {
      body: wrapIcon("")
    }
  }
};

// src/diagrams/architecture/svgDraw.ts
var drawEdges = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(async function(edgesEl, cy) {
  const padding = getConfigField("padding");
  const iconSize = getConfigField("iconSize");
  const halfIconSize = iconSize / 2;
  const arrowSize = iconSize / 6;
  const halfArrowSize = arrowSize / 2;
  await Promise.all(
    cy.edges().map(async (edge) => {
      const {
        source,
        sourceDir,
        sourceArrow,
        sourceGroup,
        target,
        targetDir,
        targetArrow,
        targetGroup,
        label
      } = edgeData(edge);
      let { x: startX, y: startY } = edge[0].sourceEndpoint();
      const { x: midX, y: midY } = edge[0].midpoint();
      let { x: endX, y: endY } = edge[0].targetEndpoint();
      const groupEdgeShift = padding + 4;
      if (sourceGroup) {
        if (isArchitectureDirectionX(sourceDir)) {
          startX += sourceDir === "L" ? -groupEdgeShift : groupEdgeShift;
        } else {
          startY += sourceDir === "T" ? -groupEdgeShift : groupEdgeShift + 18;
        }
      }
      if (targetGroup) {
        if (isArchitectureDirectionX(targetDir)) {
          endX += targetDir === "L" ? -groupEdgeShift : groupEdgeShift;
        } else {
          endY += targetDir === "T" ? -groupEdgeShift : groupEdgeShift + 18;
        }
      }
      if (!sourceGroup && db.getNode(source)?.type === "junction") {
        if (isArchitectureDirectionX(sourceDir)) {
          startX += sourceDir === "L" ? halfIconSize : -halfIconSize;
        } else {
          startY += sourceDir === "T" ? halfIconSize : -halfIconSize;
        }
      }
      if (!targetGroup && db.getNode(target)?.type === "junction") {
        if (isArchitectureDirectionX(targetDir)) {
          endX += targetDir === "L" ? halfIconSize : -halfIconSize;
        } else {
          endY += targetDir === "T" ? halfIconSize : -halfIconSize;
        }
      }
      if (edge[0]._private.rscratch) {
        const g = edgesEl.insert("g");
        g.insert("path").attr("d", `M ${startX},${startY} L ${midX},${midY} L${endX},${endY} `).attr("class", "edge");
        if (sourceArrow) {
          const xShift = isArchitectureDirectionX(sourceDir) ? ArchitectureDirectionArrowShift[sourceDir](startX, arrowSize) : startX - halfArrowSize;
          const yShift = isArchitectureDirectionY(sourceDir) ? ArchitectureDirectionArrowShift[sourceDir](startY, arrowSize) : startY - halfArrowSize;
          g.insert("polygon").attr("points", ArchitectureDirectionArrow[sourceDir](arrowSize)).attr("transform", `translate(${xShift},${yShift})`).attr("class", "arrow");
        }
        if (targetArrow) {
          const xShift = isArchitectureDirectionX(targetDir) ? ArchitectureDirectionArrowShift[targetDir](endX, arrowSize) : endX - halfArrowSize;
          const yShift = isArchitectureDirectionY(targetDir) ? ArchitectureDirectionArrowShift[targetDir](endY, arrowSize) : endY - halfArrowSize;
          g.insert("polygon").attr("points", ArchitectureDirectionArrow[targetDir](arrowSize)).attr("transform", `translate(${xShift},${yShift})`).attr("class", "arrow");
        }
        if (label) {
          const axis = !isArchitectureDirectionXY(sourceDir, targetDir) ? isArchitectureDirectionX(sourceDir) ? "X" : "Y" : "XY";
          let width = 0;
          if (axis === "X") {
            width = Math.abs(startX - endX);
          } else if (axis === "Y") {
            width = Math.abs(startY - endY) / 1.5;
          } else {
            width = Math.abs(startX - endX) / 2;
          }
          const textElem = g.append("g");
          await (0,_chunk_C3MQ5ANM_mjs__WEBPACK_IMPORTED_MODULE_1__/* .createText */ .rw)(
            textElem,
            label,
            {
              useHtmlLabels: false,
              width,
              classes: "architecture-service-label"
            },
            (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .getConfig2 */ .nV)()
          );
          textElem.attr("dy", "1em").attr("alignment-baseline", "middle").attr("dominant-baseline", "middle").attr("text-anchor", "middle");
          if (axis === "X") {
            textElem.attr("transform", "translate(" + midX + ", " + midY + ")");
          } else if (axis === "Y") {
            textElem.attr("transform", "translate(" + midX + ", " + midY + ") rotate(-90)");
          } else if (axis === "XY") {
            const pair = getArchitectureDirectionPair(sourceDir, targetDir);
            if (pair && isArchitecturePairXY(pair)) {
              const bboxOrig = textElem.node().getBoundingClientRect();
              const [x, y] = getArchitectureDirectionXYFactors(pair);
              textElem.attr("dominant-baseline", "auto").attr("transform", `rotate(${-1 * x * y * 45})`);
              const bboxNew = textElem.node().getBoundingClientRect();
              textElem.attr(
                "transform",
                `
                translate(${midX}, ${midY - bboxOrig.height / 2})
                translate(${x * bboxNew.width / 2}, ${y * bboxNew.height / 2})
                rotate(${-1 * x * y * 45}, 0, ${bboxOrig.height / 2})
              `
              );
            }
          }
        }
      }
    })
  );
}, "drawEdges");
var drawGroups = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(async function(groupsEl, cy) {
  const padding = getConfigField("padding");
  const groupIconSize = padding * 0.75;
  const fontSize = getConfigField("fontSize");
  const iconSize = getConfigField("iconSize");
  const halfIconSize = iconSize / 2;
  await Promise.all(
    cy.nodes().map(async (node) => {
      const data = nodeData(node);
      if (data.type === "group") {
        const { h, w, x1, y1 } = node.boundingBox();
        groupsEl.append("rect").attr("x", x1 + halfIconSize).attr("y", y1 + halfIconSize).attr("width", w).attr("height", h).attr("class", "node-bkg");
        const groupLabelContainer = groupsEl.append("g");
        let shiftedX1 = x1;
        let shiftedY1 = y1;
        if (data.icon) {
          const bkgElem = groupLabelContainer.append("g");
          bkgElem.html(
            `<g>${await (0,_chunk_H2D2JQ3I_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getIconSVG */ .s4)(data.icon, { height: groupIconSize, width: groupIconSize, fallbackPrefix: architectureIcons.prefix })}</g>`
          );
          bkgElem.attr(
            "transform",
            "translate(" + (shiftedX1 + halfIconSize + 1) + ", " + (shiftedY1 + halfIconSize + 1) + ")"
          );
          shiftedX1 += groupIconSize;
          shiftedY1 += fontSize / 2 - 1 - 2;
        }
        if (data.label) {
          const textElem = groupLabelContainer.append("g");
          await (0,_chunk_C3MQ5ANM_mjs__WEBPACK_IMPORTED_MODULE_1__/* .createText */ .rw)(
            textElem,
            data.label,
            {
              useHtmlLabels: false,
              width: w,
              classes: "architecture-service-label"
            },
            (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .getConfig2 */ .nV)()
          );
          textElem.attr("dy", "1em").attr("alignment-baseline", "middle").attr("dominant-baseline", "start").attr("text-anchor", "start");
          textElem.attr(
            "transform",
            "translate(" + (shiftedX1 + halfIconSize + 4) + ", " + (shiftedY1 + halfIconSize + 2) + ")"
          );
        }
      }
    })
  );
}, "drawGroups");
var drawServices = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(async function(db2, elem, services) {
  for (const service of services) {
    const serviceElem = elem.append("g");
    const iconSize = getConfigField("iconSize");
    if (service.title) {
      const textElem = serviceElem.append("g");
      await (0,_chunk_C3MQ5ANM_mjs__WEBPACK_IMPORTED_MODULE_1__/* .createText */ .rw)(
        textElem,
        service.title,
        {
          useHtmlLabels: false,
          width: iconSize * 1.5,
          classes: "architecture-service-label"
        },
        (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .getConfig2 */ .nV)()
      );
      textElem.attr("dy", "1em").attr("alignment-baseline", "middle").attr("dominant-baseline", "middle").attr("text-anchor", "middle");
      textElem.attr("transform", "translate(" + iconSize / 2 + ", " + iconSize + ")");
    }
    const bkgElem = serviceElem.append("g");
    if (service.icon) {
      bkgElem.html(
        `<g>${await (0,_chunk_H2D2JQ3I_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getIconSVG */ .s4)(service.icon, { height: iconSize, width: iconSize, fallbackPrefix: architectureIcons.prefix })}</g>`
      );
    } else if (service.iconText) {
      bkgElem.html(
        `<g>${await (0,_chunk_H2D2JQ3I_mjs__WEBPACK_IMPORTED_MODULE_0__/* .getIconSVG */ .s4)("blank", { height: iconSize, width: iconSize, fallbackPrefix: architectureIcons.prefix })}</g>`
      );
      const textElemContainer = bkgElem.append("g");
      const fo = textElemContainer.append("foreignObject").attr("width", iconSize).attr("height", iconSize);
      const divElem = fo.append("div").attr("class", "node-icon-text").attr("style", `height: ${iconSize}px;`).append("div").html(service.iconText);
      const fontSize = parseInt(
        window.getComputedStyle(divElem.node(), null).getPropertyValue("font-size").replace(/\D/g, "")
      ) ?? 16;
      divElem.attr("style", `-webkit-line-clamp: ${Math.floor((iconSize - 2) / fontSize)};`);
    } else {
      bkgElem.append("path").attr("class", "node-bkg").attr("id", "node-" + service.id).attr(
        "d",
        `M0 ${iconSize} v${-iconSize} q0,-5 5,-5 h${iconSize} q5,0 5,5 v${iconSize} H0 Z`
      );
    }
    serviceElem.attr("class", "architecture-service");
    const { width, height } = serviceElem._groups[0][0].getBBox();
    service.width = width;
    service.height = height;
    db2.setElementForId(service.id, serviceElem);
  }
  return 0;
}, "drawServices");
var drawJunctions = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(function(db2, elem, junctions) {
  junctions.forEach((junction) => {
    const junctionElem = elem.append("g");
    const iconSize = getConfigField("iconSize");
    const bkgElem = junctionElem.append("g");
    bkgElem.append("rect").attr("id", "node-" + junction.id).attr("fill-opacity", "0").attr("width", iconSize).attr("height", iconSize);
    junctionElem.attr("class", "architecture-junction");
    const { width, height } = junctionElem._groups[0][0].getBBox();
    junctionElem.width = width;
    junctionElem.height = height;
    db2.setElementForId(junction.id, junctionElem);
  });
}, "drawJunctions");

// src/diagrams/architecture/architectureRenderer.ts
(0,_chunk_H2D2JQ3I_mjs__WEBPACK_IMPORTED_MODULE_0__/* .registerIconPacks */ .ef)([
  {
    name: architectureIcons.prefix,
    icons: architectureIcons
  }
]);
cytoscape__WEBPACK_IMPORTED_MODULE_8__/* ["default"] */ .Z.use((cytoscape_fcose__WEBPACK_IMPORTED_MODULE_9___default()));
function addServices(services, cy) {
  services.forEach((service) => {
    cy.add({
      group: "nodes",
      data: {
        type: "service",
        id: service.id,
        icon: service.icon,
        label: service.title,
        parent: service.in,
        width: getConfigField("iconSize"),
        height: getConfigField("iconSize")
      },
      classes: "node-service"
    });
  });
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(addServices, "addServices");
function addJunctions(junctions, cy) {
  junctions.forEach((junction) => {
    cy.add({
      group: "nodes",
      data: {
        type: "junction",
        id: junction.id,
        parent: junction.in,
        width: getConfigField("iconSize"),
        height: getConfigField("iconSize")
      },
      classes: "node-junction"
    });
  });
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(addJunctions, "addJunctions");
function positionNodes(db2, cy) {
  cy.nodes().map((node) => {
    const data = nodeData(node);
    if (data.type === "group") {
      return;
    }
    data.x = node.position().x;
    data.y = node.position().y;
    const nodeElem = db2.getElementById(data.id);
    nodeElem.attr("transform", "translate(" + (data.x || 0) + "," + (data.y || 0) + ")");
  });
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(positionNodes, "positionNodes");
function addGroups(groups, cy) {
  groups.forEach((group) => {
    cy.add({
      group: "nodes",
      data: {
        type: "group",
        id: group.id,
        icon: group.icon,
        label: group.title,
        parent: group.in
      },
      classes: "node-group"
    });
  });
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(addGroups, "addGroups");
function addEdges(edges, cy) {
  edges.forEach((parsedEdge) => {
    const { lhsId, rhsId, lhsInto, lhsGroup, rhsInto, lhsDir, rhsDir, rhsGroup, title } = parsedEdge;
    const edgeType = isArchitectureDirectionXY(parsedEdge.lhsDir, parsedEdge.rhsDir) ? "segments" : "straight";
    const edge = {
      id: `${lhsId}-${rhsId}`,
      label: title,
      source: lhsId,
      sourceDir: lhsDir,
      sourceArrow: lhsInto,
      sourceGroup: lhsGroup,
      sourceEndpoint: lhsDir === "L" ? "0 50%" : lhsDir === "R" ? "100% 50%" : lhsDir === "T" ? "50% 0" : "50% 100%",
      target: rhsId,
      targetDir: rhsDir,
      targetArrow: rhsInto,
      targetGroup: rhsGroup,
      targetEndpoint: rhsDir === "L" ? "0 50%" : rhsDir === "R" ? "100% 50%" : rhsDir === "T" ? "50% 0" : "50% 100%"
    };
    cy.add({
      group: "edges",
      data: edge,
      classes: edgeType
    });
  });
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(addEdges, "addEdges");
function getAlignments(db2, spatialMaps, groupAlignments) {
  const flattenAlignments = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((alignmentObj, alignmentDir) => {
    return Object.entries(alignmentObj).reduce(
      (prev, [dir, alignments2]) => {
        let cnt = 0;
        const arr = Object.entries(alignments2);
        if (arr.length === 1) {
          prev[dir] = arr[0][1];
          return prev;
        }
        for (let i = 0; i < arr.length - 1; i++) {
          for (let j = i + 1; j < arr.length; j++) {
            const [aGroupId, aNodeIds] = arr[i];
            const [bGroupId, bNodeIds] = arr[j];
            const alignment = groupAlignments[aGroupId]?.[bGroupId];
            if (alignment === alignmentDir) {
              prev[dir] ??= [];
              prev[dir] = [...prev[dir], ...aNodeIds, ...bNodeIds];
            } else if (aGroupId === "default" || bGroupId === "default") {
              prev[dir] ??= [];
              prev[dir] = [...prev[dir], ...aNodeIds, ...bNodeIds];
            } else {
              const keyA = `${dir}-${cnt++}`;
              prev[keyA] = aNodeIds;
              const keyB = `${dir}-${cnt++}`;
              prev[keyB] = bNodeIds;
            }
          }
        }
        return prev;
      },
      {}
    );
  }, "flattenAlignments");
  const alignments = spatialMaps.map((spatialMap) => {
    const horizontalAlignments = {};
    const verticalAlignments = {};
    Object.entries(spatialMap).forEach(([id, [x, y]]) => {
      const nodeGroup = db2.getNode(id)?.in ?? "default";
      horizontalAlignments[y] ??= {};
      horizontalAlignments[y][nodeGroup] ??= [];
      horizontalAlignments[y][nodeGroup].push(id);
      verticalAlignments[x] ??= {};
      verticalAlignments[x][nodeGroup] ??= [];
      verticalAlignments[x][nodeGroup].push(id);
    });
    return {
      horiz: Object.values(flattenAlignments(horizontalAlignments, "horizontal")).filter(
        (arr) => arr.length > 1
      ),
      vert: Object.values(flattenAlignments(verticalAlignments, "vertical")).filter(
        (arr) => arr.length > 1
      )
    };
  });
  const [horizontal, vertical] = alignments.reduce(
    ([prevHoriz, prevVert], { horiz, vert }) => {
      return [
        [...prevHoriz, ...horiz],
        [...prevVert, ...vert]
      ];
    },
    [[], []]
  );
  return {
    horizontal,
    vertical
  };
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(getAlignments, "getAlignments");
function getRelativeConstraints(spatialMaps) {
  const relativeConstraints = [];
  const posToStr = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((pos) => `${pos[0]},${pos[1]}`, "posToStr");
  const strToPos = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)((pos) => pos.split(",").map((p) => parseInt(p)), "strToPos");
  spatialMaps.forEach((spatialMap) => {
    const invSpatialMap = Object.fromEntries(
      Object.entries(spatialMap).map(([id, pos]) => [posToStr(pos), id])
    );
    const queue = [posToStr([0, 0])];
    const visited = {};
    const directions = {
      L: [-1, 0],
      R: [1, 0],
      T: [0, 1],
      B: [0, -1]
    };
    while (queue.length > 0) {
      const curr = queue.shift();
      if (curr) {
        visited[curr] = 1;
        const currId = invSpatialMap[curr];
        if (currId) {
          const currPos = strToPos(curr);
          Object.entries(directions).forEach(([dir, shift]) => {
            const newPos = posToStr([currPos[0] + shift[0], currPos[1] + shift[1]]);
            const newId = invSpatialMap[newPos];
            if (newId && !visited[newPos]) {
              queue.push(newPos);
              relativeConstraints.push({
                [ArchitectureDirectionName[dir]]: newId,
                [ArchitectureDirectionName[getOppositeArchitectureDirection(dir)]]: currId,
                gap: 1.5 * getConfigField("iconSize")
              });
            }
          });
        }
      }
    }
  });
  return relativeConstraints;
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(getRelativeConstraints, "getRelativeConstraints");
function layoutArchitecture(services, junctions, groups, edges, db2, { spatialMaps, groupAlignments }) {
  return new Promise((resolve) => {
    const renderEl = (0,d3__WEBPACK_IMPORTED_MODULE_10__/* .select */ .Ys)("body").append("div").attr("id", "cy").attr("style", "display:none");
    const cy = (0,cytoscape__WEBPACK_IMPORTED_MODULE_8__/* ["default"] */ .Z)({
      container: document.getElementById("cy"),
      style: [
        {
          selector: "edge",
          style: {
            "curve-style": "straight",
            label: "data(label)",
            "source-endpoint": "data(sourceEndpoint)",
            "target-endpoint": "data(targetEndpoint)"
          }
        },
        {
          selector: "edge.segments",
          style: {
            "curve-style": "segments",
            "segment-weights": "0",
            "segment-distances": [0.5],
            // @ts-ignore Incorrect library types
            "edge-distances": "endpoints",
            "source-endpoint": "data(sourceEndpoint)",
            "target-endpoint": "data(targetEndpoint)"
          }
        },
        {
          selector: "node",
          style: {
            // @ts-ignore Incorrect library types
            "compound-sizing-wrt-labels": "include"
          }
        },
        {
          selector: "node[label]",
          style: {
            "text-valign": "bottom",
            "text-halign": "center",
            "font-size": `${getConfigField("fontSize")}px`
          }
        },
        {
          selector: ".node-service",
          style: {
            label: "data(label)",
            width: "data(width)",
            height: "data(height)"
          }
        },
        {
          selector: ".node-junction",
          style: {
            width: "data(width)",
            height: "data(height)"
          }
        },
        {
          selector: ".node-group",
          style: {
            // @ts-ignore Incorrect library types
            padding: `${getConfigField("padding")}px`
          }
        }
      ]
    });
    renderEl.remove();
    addGroups(groups, cy);
    addServices(services, cy);
    addJunctions(junctions, cy);
    addEdges(edges, cy);
    const alignmentConstraint = getAlignments(db2, spatialMaps, groupAlignments);
    const relativePlacementConstraint = getRelativeConstraints(spatialMaps);
    const layout = cy.layout({
      name: "fcose",
      quality: "proof",
      styleEnabled: false,
      animate: false,
      nodeDimensionsIncludeLabels: false,
      // Adjust the edge parameters if it passes through the border of a group
      // Hacky fix for: https://github.com/iVis-at-Bilkent/cytoscape.js-fcose/issues/67
      idealEdgeLength(edge) {
        const [nodeA, nodeB] = edge.connectedNodes();
        const { parent: parentA } = nodeData(nodeA);
        const { parent: parentB } = nodeData(nodeB);
        const elasticity = parentA === parentB ? 1.5 * getConfigField("iconSize") : 0.5 * getConfigField("iconSize");
        return elasticity;
      },
      edgeElasticity(edge) {
        const [nodeA, nodeB] = edge.connectedNodes();
        const { parent: parentA } = nodeData(nodeA);
        const { parent: parentB } = nodeData(nodeB);
        const elasticity = parentA === parentB ? 0.45 : 1e-3;
        return elasticity;
      },
      alignmentConstraint,
      relativePlacementConstraint
    });
    layout.one("layoutstop", () => {
      function getSegmentWeights(source, target, pointX, pointY) {
        let W, D;
        const { x: sX, y: sY } = source;
        const { x: tX, y: tY } = target;
        D = (pointY - sY + (sX - pointX) * (sY - tY) / (sX - tX)) / Math.sqrt(1 + Math.pow((sY - tY) / (sX - tX), 2));
        W = Math.sqrt(Math.pow(pointY - sY, 2) + Math.pow(pointX - sX, 2) - Math.pow(D, 2));
        const distAB = Math.sqrt(Math.pow(tX - sX, 2) + Math.pow(tY - sY, 2));
        W = W / distAB;
        let delta1 = (tX - sX) * (pointY - sY) - (tY - sY) * (pointX - sX);
        switch (true) {
          case delta1 >= 0:
            delta1 = 1;
            break;
          case delta1 < 0:
            delta1 = -1;
            break;
        }
        let delta2 = (tX - sX) * (pointX - sX) + (tY - sY) * (pointY - sY);
        switch (true) {
          case delta2 >= 0:
            delta2 = 1;
            break;
          case delta2 < 0:
            delta2 = -1;
            break;
        }
        D = Math.abs(D) * delta1;
        W = W * delta2;
        return {
          distances: D,
          weights: W
        };
      }
      (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(getSegmentWeights, "getSegmentWeights");
      cy.startBatch();
      for (const edge of Object.values(cy.edges())) {
        if (edge.data?.()) {
          const { x: sX, y: sY } = edge.source().position();
          const { x: tX, y: tY } = edge.target().position();
          if (sX !== tX && sY !== tY) {
            const sEP = edge.sourceEndpoint();
            const tEP = edge.targetEndpoint();
            const { sourceDir } = edgeData(edge);
            const [pointX, pointY] = isArchitectureDirectionY(sourceDir) ? [sEP.x, tEP.y] : [tEP.x, sEP.y];
            const { weights, distances } = getSegmentWeights(sEP, tEP, pointX, pointY);
            edge.style("segment-distances", distances);
            edge.style("segment-weights", weights);
          }
        }
      }
      cy.endBatch();
      layout.run();
    });
    layout.run();
    cy.ready((e) => {
      _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .log */ .cM.info("Ready", e);
      resolve(cy);
    });
  });
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(layoutArchitecture, "layoutArchitecture");
var draw = /* @__PURE__ */ (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .__name */ .eW)(async (text, id, _version, diagObj) => {
  const db2 = diagObj.db;
  const services = db2.getServices();
  const junctions = db2.getJunctions();
  const groups = db2.getGroups();
  const edges = db2.getEdges();
  const ds = db2.getDataStructures();
  const svg = (0,_chunk_7B677QYD_mjs__WEBPACK_IMPORTED_MODULE_5__/* .selectSvgElement */ .P)(id);
  const edgesElem = svg.append("g");
  edgesElem.attr("class", "architecture-edges");
  const servicesElem = svg.append("g");
  servicesElem.attr("class", "architecture-services");
  const groupElem = svg.append("g");
  groupElem.attr("class", "architecture-groups");
  await drawServices(db2, servicesElem, services);
  drawJunctions(db2, servicesElem, junctions);
  const cy = await layoutArchitecture(services, junctions, groups, edges, db2, ds);
  await drawEdges(edgesElem, cy);
  await drawGroups(groupElem, cy);
  positionNodes(db2, cy);
  (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_6__/* .setupGraphViewbox */ .j7)(void 0, svg, getConfigField("padding"), getConfigField("useMaxWidth"));
}, "draw");
var renderer = { draw };

// src/diagrams/architecture/architectureDiagram.ts
var diagram = {
  parser,
  db,
  renderer,
  styles: architectureStyles_default
};



/***/ }),

/***/ 80910:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   A: () => (/* binding */ populateCommonDb)
/* harmony export */ });
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(86906);


// src/diagrams/common/populateCommonDb.ts
function populateCommonDb(ast, db) {
  if (ast.accDescr) {
    db.setAccDescription?.(ast.accDescr);
  }
  if (ast.accTitle) {
    db.setAccTitle?.(ast.accTitle);
  }
  if (ast.title) {
    db.setDiagramTitle?.(ast.title);
  }
}
(0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(populateCommonDb, "populateCommonDb");




/***/ }),

/***/ 82002:
/***/ ((__unused_webpack_module, __webpack_exports__, __webpack_require__) => {

"use strict";
/* harmony export */ __webpack_require__.d(__webpack_exports__, {
/* harmony export */   A: () => (/* binding */ ImperativeState)
/* harmony export */ });
/* harmony import */ var _chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__ = __webpack_require__(86906);


// src/utils/imperativeState.ts
var ImperativeState = class {
  /**
   * @param init - Function that creates the default state.
   */
  constructor(init) {
    this.init = init;
    this.records = this.init();
  }
  static {
    (0,_chunk_YTJNT7DU_mjs__WEBPACK_IMPORTED_MODULE_0__/* .__name */ .eW)(this, "ImperativeState");
  }
  reset() {
    this.records = this.init();
  }
};




/***/ })

}]);
//# sourceMappingURL=2425.6472013fa9798ece0dc9.js.map?v=6472013fa9798ece0dc9
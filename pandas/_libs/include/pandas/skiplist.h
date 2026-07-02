/*
Copyright (c) 2016, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.

Flexibly-sized, index-able skiplist data structure for maintaining a sorted
list of values

Port of Wes McKinney's Cython version of Raymond Hettinger's original pure
Python recipe (https://rhettinger.wordpress.com/2010/02/06/lost-knowledge/)
*/

#pragma once

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct node_t node_t;

struct node_t {
  node_t **next;
  int *width;
  double value;
  int is_nil;
  int levels;
};

// Arena chunk header; node allocations are bump-allocated from the space
// following the header.
typedef struct chunk_t chunk_t;

struct chunk_t {
  chunk_t *next;
  size_t used;
  size_t capacity;
};

// First arena chunk for large windows. Windows whose full node budget fits in
// one chunk of this size are allocated exactly (identical to a single-chunk
// arena); larger windows start here and grow geometrically, so a window that
// never fills (sparse/NaN data) only commits what it touches instead of the
// whole window up front.
#ifndef SKIPLIST_INIT_CHUNK
#  define SKIPLIST_INIT_CHUNK (128 * 1024)
#endif
#define SKIPLIST_MAX_CHUNK (4 << 20)

typedef struct {
  node_t *head;
  node_t *nil;
  node_t **tmp_chain;
  int *tmp_steps;
  int size;
  int maxlevels;
  // Per-level free lists of removed nodes, chained through next[0]:
  // freelist[k] holds nodes with levels == k, so a recycled node's
  // next/width arrays always have exactly the capacity needed.
  node_t **freelist;
  chunk_t *chunks;
  size_t chunk_bytes;
  uint64_t rng_state;
} skiplist_t;

// Draw a level count with P(levels == k) = 2**-k, capped at maxlevels.
// Uses a per-skiplist xorshift64 generator: cheaper than libc rand()+log2
// and does not contend on shared libc state across threads.
static inline int skiplist_random_levels(skiplist_t *skp) {
  uint64_t x = skp->rng_state;
  int levels = 1;

  x ^= x << 13;
  x ^= x >> 7;
  x ^= x << 17;
  skp->rng_state = x;

  while ((x & 1) != 0 && levels < skp->maxlevels) {
    ++levels;
    x >>= 1;
  }
  return levels;
}

// A node and its next/width arrays live in a single allocation.
static inline size_t node_alloc_size(int levels) {
  size_t size =
      sizeof(node_t) + (size_t)levels * (sizeof(node_t *) + sizeof(int));
  // keep the arena bump pointer 8-aligned
  return (size + 7) & ~(size_t)7;
}

static inline node_t *node_init_arrays(node_t *node, int levels) {
  node->next = (node_t **)((char *)node + sizeof(node_t));
  node->width = (int *)(node->next + levels);
  node->levels = levels;
  return node;
}

static inline node_t *skiplist_alloc_node(skiplist_t *skp, int levels) {
  node_t *node;
  chunk_t *chunk;
  size_t nbytes;

  node = skp->freelist[levels];
  if (node != NULL) {
    // recycled node: next/width already point at arrays of this level count
    skp->freelist[levels] = node->next[0];
    return node;
  }

  nbytes = node_alloc_size(levels);
  chunk = skp->chunks;
  if (chunk == NULL || chunk->capacity - chunk->used < nbytes) {
    size_t capacity = skp->chunk_bytes;
    if (capacity < nbytes) {
      capacity = nbytes;
    }
    chunk = (chunk_t *)malloc(sizeof(chunk_t) + capacity);
    if (chunk == NULL) {
      return NULL;
    }
    chunk->used = 0;
    chunk->capacity = capacity;
    chunk->next = skp->chunks;
    skp->chunks = chunk;
    // next chunk doubles (capped), so filling a full window costs O(log)
    // allocations while an underfilled one stops early
    if (skp->chunk_bytes < SKIPLIST_MAX_CHUNK) {
      size_t grown = skp->chunk_bytes * 2;
      skp->chunk_bytes =
          grown > SKIPLIST_MAX_CHUNK ? SKIPLIST_MAX_CHUNK : grown;
    }
  }
  node = (node_t *)((char *)chunk + sizeof(chunk_t) + chunk->used);
  chunk->used += nbytes;
  return node_init_arrays(node, levels);
}

static inline void skiplist_free_node(skiplist_t *skp, node_t *node) {
  node->next[0] = skp->freelist[node->levels];
  skp->freelist[node->levels] = node;
}

static inline void skiplist_destroy(skiplist_t *skp) {
  if (skp) {
    chunk_t *chunk = skp->chunks;
    while (chunk != NULL) {
      chunk_t *next = chunk->next;
      free(chunk);
      chunk = next;
    }
    // head is the start of the combined head/NIL/tmp/freelist allocation
    free(skp->head);
    free(skp);
  }
}

static inline skiplist_t *skiplist_init(int expected_size) {
  skiplist_t *result;
  node_t *NIL, *head;
  int maxlevels, i;
  size_t head_bytes, nil_bytes, chain_bytes, freelist_bytes;
  uint64_t chunk_target;
  char *buf;

  maxlevels = 1 + (int)log2((double)expected_size);
  result = (skiplist_t *)malloc(sizeof(skiplist_t));
  if (result == NULL) {
    return NULL;
  }

  // One allocation for the head node, NIL node, and the fixed-size
  // scratch/freelist arrays; value nodes come from the arena chunks.
  head_bytes = node_alloc_size(maxlevels);
  nil_bytes = node_alloc_size(0);
  chain_bytes = (size_t)maxlevels * sizeof(node_t *);
  freelist_bytes = (size_t)(maxlevels + 1) * sizeof(node_t *);
  buf = (char *)malloc(head_bytes + nil_bytes + chain_bytes + freelist_bytes +
                       (size_t)maxlevels * sizeof(int));
  if (buf == NULL) {
    free(result);
    return NULL;
  }

  head = node_init_arrays((node_t *)buf, maxlevels);
  NIL = node_init_arrays((node_t *)(buf + head_bytes), 0);
  result->tmp_chain = (node_t **)(buf + head_bytes + nil_bytes);
  result->freelist = (node_t **)(buf + head_bytes + nil_bytes + chain_bytes);
  result->tmp_steps =
      (int *)(buf + head_bytes + nil_bytes + chain_bytes + freelist_bytes);

  head->value = NAN;
  head->is_nil = 0;
  NIL->value = 0.0;
  NIL->is_nil = 1;

  memset(result->freelist, 0, freelist_bytes);
  for (i = 0; i < maxlevels; ++i) {
    head->next[i] = NIL;
    head->width[i] = 1;
  }

  result->head = head;
  result->nil = NIL;
  result->maxlevels = maxlevels;
  result->size = 0;
  result->chunks = NULL;
  result->rng_state = 0x9E3779B97F4A7C15ULL ^ (uint64_t)expected_size;
  // First chunk covers a full window's node budget (~96 B/node incl. free-list
  // slack) but is capped at SKIPLIST_INIT_CHUNK; bigger windows grow from there
  // geometrically. Computed in 64 bits so expected_size cannot overflow size_t.
  chunk_target = (uint64_t)expected_size * 96 + 64;
  if (chunk_target > SKIPLIST_INIT_CHUNK) {
    chunk_target = SKIPLIST_INIT_CHUNK;
  }
  result->chunk_bytes = (size_t)chunk_target;

  return result;
}

// Empty the skiplist, retaining the already-allocated arena chunks.
static inline void skiplist_reset(skiplist_t *skp) {
  int i;
  chunk_t *chunk;

  for (chunk = skp->chunks; chunk != NULL; chunk = chunk->next) {
    chunk->used = 0;
  }
  memset(skp->freelist, 0, (size_t)(skp->maxlevels + 1) * sizeof(node_t *));
  for (i = 0; i < skp->maxlevels; ++i) {
    skp->head->next[i] = skp->nil;
    skp->head->width[i] = 1;
  }
  skp->size = 0;
}

// 1 if left < right, 0 if left == right, -1 if left > right
static inline int _node_cmp(node_t *node, double value) {
  if (node->is_nil || node->value > value) {
    return -1;
  } else if (node->value < value) {
    return 1;
  } else {
    return 0;
  }
}

static inline double skiplist_get(skiplist_t *skp, int i, int *ret) {
  node_t *node;
  int level;

  if (i < 0 || i >= skp->size) {
    *ret = 0;
    return 0;
  }

  node = skp->head;
  ++i;
  for (level = skp->maxlevels - 1; level >= 0; --level) {
    while (node->width[level] <= i) {
      i -= node->width[level];
      node = node->next[level];
    }
  }

  *ret = 1;
  return node->value;
}

// Fetch the values at indices i and i + 1 in a single traversal; equivalent
// to two skiplist_get calls. Returns 0 if either index is out of bounds.
static inline int skiplist_get_pair(skiplist_t *skp, int i, double *lower,
                                    double *higher) {
  node_t *node;
  int level;

  if (i < 0 || i + 1 >= skp->size) {
    *lower = 0;
    *higher = 0;
    return 0;
  }

  node = skp->head;
  ++i;
  for (level = skp->maxlevels - 1; level >= 0; --level) {
    while (node->width[level] <= i) {
      i -= node->width[level];
      node = node->next[level];
    }
  }

  *lower = node->value;
  *higher = node->next[0]->value;
  return 1;
}

// Returns the lowest rank of all elements with value `value`, as opposed to the
// highest rank returned by `skiplist_insert`.
static inline int skiplist_min_rank(skiplist_t *skp, double value) {
  node_t *node;
  int level, rank = 0;

  node = skp->head;
  for (level = skp->maxlevels - 1; level >= 0; --level) {
    while (_node_cmp(node->next[level], value) > 0) {
      rank += node->width[level];
      node = node->next[level];
    }
  }

  return rank + 1;
}

// Returns the rank of the inserted element. When there are duplicates,
// `rank` is the highest of the group, i.e. the 'max' method of
// https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.rank.html
static inline int skiplist_insert(skiplist_t *skp, double value) {
  node_t *node, *prevnode, *newnode, *next_at_level;
  int *steps_at_level;
  int size, steps, level, rank = 0;
  node_t **chain;

  chain = skp->tmp_chain;

  steps_at_level = skp->tmp_steps;
  memset(steps_at_level, 0, skp->maxlevels * sizeof(int));

  node = skp->head;

  for (level = skp->maxlevels - 1; level >= 0; --level) {
    next_at_level = node->next[level];
    while (_node_cmp(next_at_level, value) >= 0) {
      steps_at_level[level] += node->width[level];
      rank += node->width[level];
      node = next_at_level;
      next_at_level = node->next[level];
    }
    chain[level] = node;
  }

  size = skiplist_random_levels(skp);

  newnode = skiplist_alloc_node(skp, size);
  if (newnode == NULL) {
    return -1;
  }
  newnode->value = value;
  newnode->is_nil = 0;
  steps = 0;

  for (level = 0; level < size; ++level) {
    prevnode = chain[level];
    newnode->next[level] = prevnode->next[level];

    prevnode->next[level] = newnode;

    newnode->width[level] = prevnode->width[level] - steps;
    prevnode->width[level] = steps + 1;

    steps += steps_at_level[level];
  }

  for (level = size; level < skp->maxlevels; ++level) {
    chain[level]->width[level] += 1;
  }

  ++(skp->size);

  return rank + 1;
}

static inline int skiplist_remove(skiplist_t *skp, double value) {
  int level, size;
  node_t *node, *prevnode, *tmpnode, *next_at_level;
  node_t **chain;

  chain = skp->tmp_chain;
  node = skp->head;

  for (level = skp->maxlevels - 1; level >= 0; --level) {
    next_at_level = node->next[level];
    while (_node_cmp(next_at_level, value) > 0) {
      node = next_at_level;
      next_at_level = node->next[level];
    }
    chain[level] = node;
  }

  tmpnode = chain[0]->next[0];
  if (value != tmpnode->value) {
    return 0;
  }

  size = tmpnode->levels;

  for (level = 0; level < size; ++level) {
    prevnode = chain[level];

    prevnode->width[level] += tmpnode->width[level] - 1;
    prevnode->next[level] = tmpnode->next[level];
  }

  for (level = size; level < skp->maxlevels; ++level) {
    --(chain[level]->width[level]);
  }

  skiplist_free_node(skp, tmpnode);
  --(skp->size);
  return 1;
}

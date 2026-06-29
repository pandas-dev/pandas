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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline float __skiplist_nanf(void) {
  const union {
    int __i;
    float __f;
  } __bint = {0x7fc00000UL};
  return __bint.__f;
}
#define PANDAS_NAN ((double)__skiplist_nanf())

static inline double Log2(double val) { return log(val) / log(2.); }

#if defined(__GNUC__) || defined(__clang__)
#define SL_LIKELY(x) __builtin_expect(!!(x), 1)
#define SL_UNLIKELY(x) __builtin_expect(!!(x), 0)
#define SL_PREFETCH_R(addr) __builtin_prefetch((const void *)(addr), 0, 3)
#else
#define SL_LIKELY(x) (x)
#define SL_UNLIKELY(x) (x)
#define SL_PREFETCH_R(addr) ((void)0)
#endif

typedef struct node_t node_t;

struct node_t {
  double value;
  int is_nil;
  int levels;
  int ref_count;
  int _pad;
  node_t **next;
  int *width;
};

typedef struct {
  node_t *head;
  node_t **tmp_chain;
  int *tmp_steps;
  int size;
  int maxlevels;
} skiplist_t;

static inline double urand(void) {
  return ((double)rand() + 1) / ((double)RAND_MAX + 2);
}

static inline int int_min(int a, int b) { return a < b ? a : b; }

static inline node_t *node_init(double value, int levels) {
  size_t next_off = sizeof(node_t);
  size_t width_off = next_off + (size_t)levels * sizeof(node_t *);
  size_t total = width_off + (size_t)levels * sizeof(int);

  char *mem = (char *)malloc(total);
  if (SL_UNLIKELY(!mem)) {
    return NULL;
  }

  node_t *result = (node_t *)mem;
  result->value = value;
  result->levels = levels;
  result->is_nil = 0;
  result->ref_count = 0;

  if (levels > 0) {
    result->next = (node_t **)(mem + next_off);
    result->width = (int *)(mem + width_off);
  } else {
    result->next = NULL;
    result->width = NULL;
  }

  return result;
}

// do this ourselves
static inline void node_incref(node_t *node) { ++(node->ref_count); }

static inline void node_decref(node_t *node) { --(node->ref_count); }

static void node_destroy(node_t *node) {
  int i;
  if (node) {
    if (node->ref_count <= 1) {
      for (i = 0; i < node->levels; ++i) {
        node_destroy(node->next[i]);
      }
      free(node);
    } else {
      node_decref(node);
    }
  }
}

static inline void skiplist_destroy(skiplist_t *skp) {
  if (skp) {
    node_destroy(skp->head);
    free(skp->tmp_steps);
    free(skp->tmp_chain);
    free(skp);
  }
}

static inline skiplist_t *skiplist_init(int expected_size) {
  skiplist_t *result;
  node_t *NIL, *head;
  int maxlevels, i;

  maxlevels = 1 + Log2((double)expected_size);
  result = (skiplist_t *)malloc(sizeof(skiplist_t));
  if (SL_UNLIKELY(!result)) {
    return NULL;
  }
  result->tmp_chain = (node_t **)malloc(maxlevels * sizeof(node_t *));
  result->tmp_steps = (int *)malloc(maxlevels * sizeof(int));
  result->maxlevels = maxlevels;
  result->size = 0;

  head = result->head = node_init(PANDAS_NAN, maxlevels);
  NIL = node_init(0.0, 0);

  if (!(result->tmp_chain && result->tmp_steps && result->head && NIL)) {
    skiplist_destroy(result);
    node_destroy(NIL);
    return NULL;
  }

  node_incref(head);

  NIL->is_nil = 1;

  for (i = 0; i < maxlevels; ++i) {
    head->next[i] = NIL;
    head->width[i] = 1;
    node_incref(NIL);
  }

  return result;
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

  if (SL_UNLIKELY(i < 0 || i >= skp->size)) {
    *ret = 0;
    return 0;
  }

  node = skp->head;
  ++i;
  for (level = skp->maxlevels - 1; level >= 0; --level) {
    while (node->width[level] <= i) {
      i -= node->width[level];
      node = node->next[level];
      if (i > 1) {
        SL_PREFETCH_R(node->next[level]);
      }
    }
    if (SL_LIKELY(level > 0)) {
      SL_PREFETCH_R(node->next[level - 1]);
    }
  }

  *ret = 1;
  return node->value;
}

// Returns the lowest rank of all elements with value `value`, as opposed to the
// highest rank returned by `skiplist_insert`.
static inline int skiplist_min_rank(skiplist_t *skp, double value) {
  node_t *node, *next_at_level;
  int level, rank = 0;

  node = skp->head;
  for (level = skp->maxlevels - 1; level >= 0; --level) {
    next_at_level = node->next[level];
    SL_PREFETCH_R(next_at_level);
    while (_node_cmp(next_at_level, value) > 0) {
      rank += node->width[level];
      node = next_at_level;
      next_at_level = node->next[level];
      SL_PREFETCH_R(next_at_level);
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
    SL_PREFETCH_R(next_at_level);

    while (1) {
      int nil = next_at_level->is_nil;
      double nval = next_at_level->value;

      if (nil) {
        break;
      }
      if (nval < value) {
        break;
      }

      {
        int w = node->width[level];
        steps_at_level[level] += w;
        rank += w;
      }
      node = next_at_level;
      next_at_level = node->next[level];
      SL_PREFETCH_R(next_at_level);
    }
    chain[level] = node;

    if (SL_LIKELY(level > 0)) {
      SL_PREFETCH_R(node->next[level - 1]);
    }
  }

  size = int_min(skp->maxlevels, 1 - ((int)Log2(urand())));

  newnode = node_init(value, size);
  if (SL_UNLIKELY(!newnode)) {
    return -1;
  }
  steps = 0;

  for (level = 0; level < size; ++level) {
    prevnode = chain[level];

    if (level + 1 < size) {
      SL_PREFETCH_R(chain[level + 1]);
    }

    newnode->next[level] = prevnode->next[level];

    prevnode->next[level] = newnode;
    node_incref(newnode); // increment the reference count

    {
      int pw = prevnode->width[level];
      newnode->width[level] = pw - steps;
      prevnode->width[level] = steps + 1;
      steps += steps_at_level[level];
    }
  }

  {
    int lv = size;
    int ml = skp->maxlevels;
    for (; lv + 3 < ml; lv += 4) {
      chain[lv]->width[lv] += 1;
      chain[lv + 1]->width[lv + 1] += 1;
      chain[lv + 2]->width[lv + 2] += 1;
      chain[lv + 3]->width[lv + 3] += 1;
    }
    for (; lv < ml; ++lv) {
      chain[lv]->width[lv] += 1;
    }
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
    SL_PREFETCH_R(next_at_level);

    while (1) {
      int nil = next_at_level->is_nil;
      double nval = next_at_level->value;

      if (nil) {
        break;
      }
      if (nval <= value) {
        break;
      }

      node = next_at_level;
      next_at_level = node->next[level];
      SL_PREFETCH_R(next_at_level);
    }
    chain[level] = node;

    if (SL_LIKELY(level > 0)) {
      SL_PREFETCH_R(node->next[level - 1]);
    }
  }

  if (SL_UNLIKELY(value != chain[0]->next[0]->value)) {
    return 0;
  }

  size = chain[0]->next[0]->levels;

  for (level = 0; level < size; ++level) {
    prevnode = chain[level];

    tmpnode = prevnode->next[level];

    if (level + 1 < size) {
      SL_PREFETCH_R(chain[level + 1]);
    }

    {
      int tw = tmpnode->width[level];
      prevnode->width[level] += tw - 1;
    }
    prevnode->next[level] = tmpnode->next[level];

    tmpnode->next[level] = NULL;
    node_destroy(tmpnode); // decrement refcount or free
  }

  {
    int lv = size;
    int ml = skp->maxlevels;
    for (; lv + 3 < ml; lv += 4) {
      chain[lv]->width[lv] -= 1;
      chain[lv + 1]->width[lv + 1] -= 1;
      chain[lv + 2]->width[lv + 2] -= 1;
      chain[lv + 3]->width[lv + 3] -= 1;
    }
    for (; lv < ml; ++lv) {
      chain[lv]->width[lv] -= 1;
    }
  }

  --(skp->size);
  return 1;
}

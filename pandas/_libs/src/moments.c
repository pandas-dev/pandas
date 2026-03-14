/*
Copyright (c) 2026, PyData Development Team
All rights reserved.

Distributed under the terms of the BSD Simplified License.

The full license is in the LICENSE file, distributed with this software.
*/

#include "pandas/moments.h"
#include <stdint.h>

/// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics
void moments_merge(Moments *acc, const Moments *src, int max_moment) {
  if (acc->n == 0) {
    acc->n = src->n;
    acc->mean = src->mean;
    acc->m2 = src->m2;
    acc->m3 = src->m3;
    acc->m4 = src->m4;
    return;
  }
  if (src->n == 0) {
    return;
  }

  double n_a = (double)acc->n;
  double n_b = (double)src->n;
  acc->n += src->n;
  double delta = src->mean - acc->mean;
  double delta_n = delta / (double)acc->n;
  double term1 = delta * delta_n * n_a * n_b;

  if (max_moment >= 4) {
    acc->m4 +=
        src->m4 +
        delta_n *
            (4.0 * (n_a * src->m3 - n_b * acc->m3) +
             delta_n * (6.0 * (n_a * n_a * src->m2 + n_b * n_b * acc->m2) +
                        term1 * (n_a * n_a - n_a * n_b + n_b * n_b)));
  }
  if (max_moment >= 3) {
    acc->m3 += src->m3 + delta_n * (3.0 * (n_a * src->m2 - n_b * acc->m2) +
                                    term1 * (n_a - n_b));
  }
  acc->m2 += src->m2 + term1;
  acc->mean += delta_n * n_b;
}

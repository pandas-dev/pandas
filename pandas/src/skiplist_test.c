#include "skiplist.h"

double urand2() {
  return rand() / ((double) RAND_MAX + 1);
}

int main(void) {
  skiplist_t *skp;
  double *data;
  int i, n, win, ret;

  n = 10000;
  win = 1000;
  skp = skiplist_init(win);

  data = (double*) malloc(n * sizeof(double));
  for (i = 0; i < n; ++i) {
    data[i] = urand2();
  }

  for (i = 0; i < win; ++i) {
    skiplist_insert(skp, data[i]);
  }

  double middle;
  for (i = win; i < n; ++i) {
    skiplist_remove(skp, data[i - win]);
    skiplist_insert(skp, data[i]);
    middle = skiplist_get(skp, win / 2, &ret);
  }

  skiplist_destroy(skp);
  free(data);
  return 0;
}


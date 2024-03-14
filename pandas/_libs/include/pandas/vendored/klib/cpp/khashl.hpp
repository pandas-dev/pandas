#ifndef __AC_KHASHL_HPP
#define __AC_KHASHL_HPP

#include <cstdlib>    // for malloc() etc
#include <cstring>    // for memset()
#include <functional> // for std::equal_to
#include <stdint.h>   // for uint32_t

/* // ==> Code example <==
#include "khashl.hpp"
#include <cstdio>

int main(void)
{
        klib::KHashMap<uint32_t, int, std::hash<uint32_t> > h; // NB: C++98
doesn't have std::hash uint32_t k; int absent; h[43] = 1, h[53] = 2, h[63] = 3,
h[73] = 4;       // one way to insert k = h.put(53, &absent), h.value(k) = -2;
// another way to insert if (!absent) printf("already in the table\n");    //
which allows to test presence if (h.get(33) == h.end()) printf("not found!\n");
// test presence without insertion h.del(h.get(43));               // deletion
        for (k = 0; k != h.end(); ++k)  // traversal
                if (h.occupied(k))          // some buckets are not occupied;
skip them printf("%u => %d\n", h.key(k), h.value(k)); return 0;
}
*/

namespace klib {

/***********
 * HashSet *
 ***********/

template <class T, class Hash, class Eq = std::equal_to<T>,
          typename khint_t = uint32_t>
class KHashSet {
  khint_t bits, count;
  uint32_t *used;
  T *keys;
  static inline uint32_t __kh_used(const uint32_t *flag, khint_t i) {
    return flag[i >> 5] >> (i & 0x1fU) & 1U;
  };
  static inline void __kh_set_used(uint32_t *flag, khint_t i) {
    flag[i >> 5] |= 1U << (i & 0x1fU);
  };
  static inline void __kh_set_unused(uint32_t *flag, khint_t i) {
    flag[i >> 5] &= ~(1U << (i & 0x1fU));
  };
  static inline khint_t __kh_fsize(khint_t m) { return m < 32 ? 1 : m >> 5; }
  static inline khint_t __kh_h2b(uint32_t hash, khint_t bits) {
    return hash * 2654435769U >> (32 - bits);
  }
  static inline khint_t __kh_h2b(uint64_t hash, khint_t bits) {
    return hash * 11400714819323198485ULL >> (64 - bits);
  }

public:
  KHashSet() : bits(0), count(0), used(0), keys(0){};
  ~KHashSet() {
    std::free(used);
    std::free(keys);
  };
  inline khint_t n_buckets() const { return used ? khint_t(1) << bits : 0; }
  inline khint_t end() const { return n_buckets(); }
  inline khint_t size() const { return count; }
  inline T &key(khint_t x) { return keys[x]; };
  inline bool occupied(khint_t x) const { return (__kh_used(used, x) != 0); }
  void clear(void) {
    if (!used)
      return;
    memset(used, 0, __kh_fsize(n_buckets()) * sizeof(uint32_t));
    count = 0;
  }
  khint_t get(const T &key) const {
    khint_t i, last, mask, nb;
    if (keys == 0)
      return 0;
    nb = n_buckets();
    mask = nb - khint_t(1);
    i = last = __kh_h2b(Hash()(key), bits);
    while (__kh_used(used, i) && !Eq()(keys[i], key)) {
      i = (i + khint_t(1)) & mask;
      if (i == last)
        return nb;
    }
    return !__kh_used(used, i) ? nb : i;
  }
  int resize(khint_t new_nb) {
    uint32_t *new_used = 0;
    khint_t j = 0, x = new_nb, nb, new_bits, new_mask;
    while ((x >>= khint_t(1)) != 0)
      ++j;
    if (new_nb & (new_nb - 1))
      ++j;
    new_bits = j > 2 ? j : 2;
    new_nb = khint_t(1) << new_bits;
    if (count > (new_nb >> 1) + (new_nb >> 2))
      return 0; /* requested size is too small */
    new_used = (uint32_t *)std::malloc(__kh_fsize(new_nb) * sizeof(uint32_t));
    memset(new_used, 0, __kh_fsize(new_nb) * sizeof(uint32_t));
    if (!new_used)
      return -1; /* not enough memory */
    nb = n_buckets();
    if (nb < new_nb) { /* expand */
      T *new_keys = (T *)std::realloc(keys, new_nb * sizeof(T));
      if (!new_keys) {
        std::free(new_used);
        return -1;
      }
      keys = new_keys;
    } /* otherwise shrink */
    new_mask = new_nb - 1;
    for (j = 0; j != nb; ++j) {
      if (!__kh_used(used, j))
        continue;
      T key = keys[j];
      __kh_set_unused(used, j);
      while (1) { /* kick-out process; sort of like in Cuckoo hashing */
        khint_t i;
        i = __kh_h2b(Hash()(key), new_bits);
        while (__kh_used(new_used, i))
          i = (i + khint_t(1)) & new_mask;
        __kh_set_used(new_used, i);
        if (i < nb && __kh_used(used, i)) { /* kick out the existing element */
          {
            T tmp = keys[i];
            keys[i] = key;
            key = tmp;
          }
          __kh_set_unused(used,
                          i); /* mark it as deleted in the old hash table */
        } else {              /* write the element and jump out of the loop */
          keys[i] = key;
          break;
        }
      }
    }
    if (nb > new_nb) /* shrink the hash table */
      keys = (T *)std::realloc(keys, new_nb * sizeof(T));
    std::free(used); /* free the working space */
    used = new_used, bits = new_bits;
    return 0;
  }
  khint_t put(const T &key, int *absent_ = 0) {
    khint_t nb, i, last, mask;
    int absent = -1;
    nb = n_buckets();
    if (count >= (nb >> 1) + (nb >> 2)) { /* rehashing */
      if (resize(nb + khint_t(1)) < 0) {
        if (absent_)
          *absent_ = -1;
        return nb;
      }
      nb = n_buckets();
    } /* TODO: to implement automatically shrinking; resize() already support
         shrinking */
    mask = nb - 1;
    i = last = __kh_h2b(Hash()(key), bits);
    while (__kh_used(used, i) && !Eq()(keys[i], key)) {
      i = (i + 1U) & mask;
      if (i == last)
        break;
    }
    if (!__kh_used(used, i)) { /* not present at all */
      keys[i] = key;
      __kh_set_used(used, i);
      ++count, absent = 1;
    } else
      absent = 0; /* Don't touch keys[i] if present */
    if (absent_)
      *absent_ = absent;
    return i;
  }
  int del(khint_t i) {
    khint_t j = i, k, mask, nb = n_buckets();
    if (keys == 0 || i >= nb)
      return 0;
    mask = nb - khint_t(1);
    while (1) {
      j = (j + khint_t(1)) & mask;
      if (j == i || !__kh_used(used, j))
        break; /* j==i only when the table is completely full */
      k = __kh_h2b(Hash()(keys[j]), bits);
      if ((j > i && (k <= i || k > j)) || (j < i && (k <= i && k > j)))
        keys[i] = keys[j], i = j;
    }
    __kh_set_unused(used, i);
    --count;
    return 1;
  }
};

/***********
 * HashMap *
 ***********/

template <class KType, class VType> struct KHashMapBucket {
  KType key;
  VType val;
};

template <class T, class Hash, typename khint_t> struct KHashMapHash {
  khint_t operator()(const T &a) const { return Hash()(a.key); }
};

template <class T, class Eq> struct KHashMapEq {
  bool operator()(const T &a, const T &b) const { return Eq()(a.key, b.key); }
};

template <class KType, class VType, class Hash, class Eq = std::equal_to<KType>,
          typename khint_t = uint32_t>
class KHashMap
    : public KHashSet<KHashMapBucket<KType, VType>,
                      KHashMapHash<KHashMapBucket<KType, VType>, Hash, khint_t>,
                      KHashMapEq<KHashMapBucket<KType, VType>, Eq>, khint_t> {
  typedef KHashMapBucket<KType, VType> bucket_t;
  typedef KHashSet<bucket_t, KHashMapHash<bucket_t, Hash, khint_t>,
                   KHashMapEq<bucket_t, Eq>, khint_t>
      hashset_t;

public:
  khint_t get(const KType &key) const {
    bucket_t t = {key, VType()};
    return hashset_t::get(t);
  }
  khint_t put(const KType &key, int *absent) {
    bucket_t t = {key, VType()};
    return hashset_t::put(t, absent);
  }
  inline KType &key(khint_t i) { return hashset_t::key(i).key; }
  inline VType &value(khint_t i) { return hashset_t::key(i).val; }
  inline VType &operator[](const KType &key) {
    bucket_t t = {key, VType()};
    return value(hashset_t::put(t));
  }
};

/****************************
 * HashSet with cached hash *
 ****************************/

template <class KType, typename khint_t> struct KHashSetCachedBucket {
  KType key;
  khint_t hash;
};

template <class T, typename khint_t> struct KHashCachedHash {
  khint_t operator()(const T &a) const { return a.hash; }
};

template <class T, class Eq> struct KHashCachedEq {
  bool operator()(const T &a, const T &b) const {
    return a.hash == b.hash && Eq()(a.key, b.key);
  }
};

template <class KType, class Hash, class Eq = std::equal_to<KType>,
          typename khint_t = uint32_t>
class KHashSetCached
    : public KHashSet<
          KHashSetCachedBucket<KType, khint_t>,
          KHashCachedHash<KHashSetCachedBucket<KType, khint_t>, khint_t>,
          KHashCachedEq<KHashSetCachedBucket<KType, khint_t>, Eq>, khint_t> {
  typedef KHashSetCachedBucket<KType, khint_t> bucket_t;
  typedef KHashSet<bucket_t, KHashCachedHash<bucket_t, khint_t>,
                   KHashCachedEq<bucket_t, Eq>, khint_t>
      hashset_t;

public:
  khint_t get(const KType &key) const {
    bucket_t t = {key, Hash()(key)};
    return hashset_t::get(t);
  }
  khint_t put(const KType &key, int *absent) {
    bucket_t t = {key, Hash()(key)};
    return hashset_t::put(t, absent);
  }
  inline KType &key(khint_t i) { return hashset_t::key(i).key; }
};

/****************************
 * HashMap with cached hash *
 ****************************/

template <class KType, class VType, typename khint_t>
struct KHashMapCachedBucket {
  KType key;
  VType val;
  khint_t hash;
};

template <class KType, class VType, class Hash, class Eq = std::equal_to<KType>,
          typename khint_t = uint32_t>
class KHashMapCached
    : public KHashSet<
          KHashMapCachedBucket<KType, VType, khint_t>,
          KHashCachedHash<KHashMapCachedBucket<KType, VType, khint_t>, khint_t>,
          KHashCachedEq<KHashMapCachedBucket<KType, VType, khint_t>, Eq>,
          khint_t> {
  typedef KHashMapCachedBucket<KType, VType, khint_t> bucket_t;
  typedef KHashSet<bucket_t, KHashCachedHash<bucket_t, khint_t>,
                   KHashCachedEq<bucket_t, Eq>, khint_t>
      hashset_t;

public:
  khint_t get(const KType &key) const {
    bucket_t t = {key, VType(), Hash()(key)};
    return hashset_t::get(t);
  }
  khint_t put(const KType &key, int *absent) {
    bucket_t t = {key, VType(), Hash()(key)};
    return hashset_t::put(t, absent);
  }
  inline KType &key(khint_t i) { return hashset_t::key(i).key; }
  inline VType &value(khint_t i) { return hashset_t::key(i).val; }
  inline VType &operator[](const KType &key) {
    bucket_t t = {key, VType(), Hash()(key)};
    return value(hashset_t::put(t));
  }
};

} // namespace klib

#endif /* __AC_KHASHL_HPP */

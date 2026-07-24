
#define SET_MINSIZE 8

#include "setobject.h"


#if defined(_MSC_VER)
#   if _MSC_VER <= 1900  /* Visual Studio 2014 */
        typedef __int8 int8_t;
        typedef __int16 int16_t;
        typedef __int32 int32_t;
        typedef __int64 int64_t;
#   endif
    /* Use _alloca() to dynamically allocate on the stack on MSVC */
    #define STACK_ALLOC(Type, Name, Size) Type * const Name = _alloca(Size);
#else
    #define STACK_ALLOC(Type, Name, Size) Type Name[Size];
#endif


/*[clinic input]
class set "PysetObject *" "&Pyset_Type"
[clinic start generated code]*/
/*[clinic end generated code: output=da39a3ee5e6b4b0d input=f157a5a0ce9589d6]*/


/*
To ensure the lookup algorithm terminates, there must be at least one Unused
slot (NULL key) in the table.
To avoid slowing down lookups on a near-full table, we resize the table when
it's USABLE_FRACTION (currently two-thirds) full.
*/

#define PERTURB_SHIFT 5

/*
Major subtleties ahead:  Most hash schemes depend on having a "good" hash
function, in the sense of simulating randomness.  Python doesn't:  its most
important hash functions (for ints) are very regular in common
cases:

  >>>[hash(i) for i in range(4)]
  [0, 1, 2, 3]

This isn't necessarily bad!  To the contrary, in a table of size 2**i, taking
the low-order i bits as the initial table index is extremely fast, and there
are no collisions at all for sets indexed by a contiguous range of ints. So
this gives better-than-random behavior in common cases, and that's very
desirable.

OTOH, when collisions occur, the tendency to fill contiguous slices of the
hash table makes a good collision resolution strategy crucial.  Taking only
the last i bits of the hash code is also vulnerable:  for example, consider
the list [i << 16 for i in range(20000)] as a set of keys.  Since ints are
their own hash codes, and this fits in a set of size 2**15, the last 15 bits
 of every hash code are all 0:  they *all* map to the same table index.

But catering to unusual cases should not slow the usual ones, so we just take
the last i bits anyway.  It's up to collision resolution to do the rest.  If
we *usually* find the key we're looking for on the first try (and, it turns
out, we usually do -- the table load factor is kept under 2/3, so the odds
are solidly in our favor), then it makes best sense to keep the initial index
computation dirt cheap.

The first half of collision resolution is to visit table indices via this
recurrence:

    j = ((5*j) + 1) mod 2**i

For any initial j in range(2**i), repeating that 2**i times generates each
int in range(2**i) exactly once (see any text on random-number generation for
proof).  By itself, this doesn't help much:  like linear probing (setting
j += 1, or j -= 1, on each loop trip), it scans the table entries in a fixed
order.  This would be bad, except that's not the only thing we do, and it's
actually *good* in the common cases where hash keys are consecutive.  In an
example that's really too small to make this entirely clear, for a table of
size 2**3 the order of indices is:

    0 -> 1 -> 6 -> 7 -> 4 -> 5 -> 2 -> 3 -> 0 [and here it's repeating]

If two things come in at index 5, the first place we look after is index 2,
not 6, so if another comes in at index 6 the collision at 5 didn't hurt it.
Linear probing is deadly in this case because there the fixed probe order
is the *same* as the order consecutive keys are likely to arrive.  But it's
extremely unlikely hash codes will follow a 5*j+1 recurrence by accident,
and certain that consecutive hash codes do not.

The other half of the strategy is to get the other bits of the hash code
into play.  This is done by initializing a (unsigned) vrbl "perturb" to the
full hash code, and changing the recurrence to:

    perturb >>= PERTURB_SHIFT;
    j = (5*j) + 1 + perturb;
    use j % 2**i as the next table index;

Now the probe sequence depends (eventually) on every bit in the hash code,
and the pseudo-scrambling property of recurring on 5*j+1 is more valuable,
because it quickly magnifies small differences in the bits that didn't affect
the initial index.  Note that because perturb is unsigned, if the recurrence
is executed often enough perturb eventually becomes and remains 0.  At that
point (very rarely reached) the recurrence is on (just) 5*j+1 again, and
that's certain to find an empty slot eventually (since it generates every int
in range(2**i), and we make sure there's always at least one empty slot).

Selecting a good value for PERTURB_SHIFT is a balancing act.  You want it
small so that the high bits of the hash code continue to affect the probe
sequence across iterations; but you want it large so that in really bad cases
the high-order hash bits have an effect on early iterations.  5 was "the
best" in minimizing total collisions across experiments Tim Peters ran (on
both normal and pathological cases), but 4 and 6 weren't significantly worse.

Historical: Reimer Behrends contributed the idea of using a polynomial-based
approach, using repeated multiplication by x in GF(2**n) where an irreducible
polynomial for each table size was chosen such that x was a primitive root.
Christian Tismer later extended that to use division by x instead, as an
efficient way to get the high bits of the hash code into play.  This scheme
also gave excellent collision statistics, but was more expensive:  two
if-tests were required inside the loop; computing "the next" index took about
the same number of operations but without as much potential parallelism
(e.g., computing 5*j can go on at the same time as computing 1+perturb in the
above, and then shifting perturb can be done while the table index is being
masked); and the PysetObject struct required a member to hold the table's
polynomial.  In Tim's experiments the current scheme ran faster, produced
equally good collision statistics, needed less code & used less memory.

*/

#define SETK_EMPTY (-1)
#define SETK_DUMMY (-2)  /* Used internally */
#define SETK_ERROR (-3)

typedef enum {
    OK = 0,
    OK_REPLACED = 1,
    ERR_NO_MEMORY = -1,
    ERR_SET_MUTATED = -2,
    ERR_ITER_EXHAUSTED = -3,
    ERR_SET_EMPTY = -4,
    ERR_CMP_FAILED = -5,
} Status;


#ifndef NDEBUG
static
int mem_cmp_zeros(void *obj, size_t n){
    int diff = 0;
    char *mem = obj;
    char *it;
    for (it = mem; it < mem + n; ++it) {
        if (*it != 0) diff += 1;
    }
    return diff;
}
#endif

#define SET_MASK(sk) ((sk)->size-1)
#define SET_GROWTH_RATE(s) ((s)->used*3)

static int
ix_size(Py_ssize_t size) {
    if ( size < 0xff ) return 1;
    if ( size < 0xffff ) return 2;
    if ( size < 0xffffffff ) return 4;
    return sizeof(int64_t);
}

#ifndef NDEBUG
/* NOTE: This function is only used in assert()s */
/* Align pointer *ptr* to pointer size */
static void*
aligned_pointer(void *ptr) {
    return (void*)aligned_size((size_t)ptr);
}
#endif

/* lookup indices.  returns SETK_EMPTY, SETK_DUMMY, or ix >=0 */
static Py_ssize_t
get_index(NB_SetKeys *sk, Py_ssize_t i)
{
    Py_ssize_t s = sk->size;
    Py_ssize_t ix;

    if (s==0) {
        /* Empty set */
        return SETK_EMPTY;
    }
    else if(s <= 0xff) {
        int8_t *indices = (int8_t*)(sk->indices);
        assert (i < sk->size);
        ix = indices[i];
    }
    else if (s <= 0xffff) {
        int16_t *indices = (int16_t*)(sk->indices);
        ix = indices[i];
    }
#if SIZEOF_VOID_P > 4
    else if (s > 0xffffffff) {
        int64_t *indices = (int64_t*)(sk->indices);
        ix = indices[i];
    }
#endif
    else {
        int32_t *indices = (int32_t*)(sk->indices);
        ix = indices[i];
    }
    assert(ix >= SETK_DUMMY);
    return ix;
}

/* write to indices. */
static void
set_index(NB_SetKeys *sk, Py_ssize_t i, Py_ssize_t ix)
{
    Py_ssize_t s = sk->size;

    assert(ix >= SETK_DUMMY);

    if (s <= 0xff) {
        int8_t *indices = (int8_t*)(sk->indices);
        assert(ix <= 0x7f);
        indices[i] = (char)ix;
    }
    else if (s <= 0xffff) {
        int16_t *indices = (int16_t*)(sk->indices);
        assert(ix <= 0x7fff);
        indices[i] = (int16_t)ix;
    }
#if SIZEOF_VOID_P > 4
    else if (s > 0xffffffff) {
        int64_t *indices = (int64_t*)(sk->indices);
        indices[i] = ix;
    }
#endif
    else {
        int32_t *indices = (int32_t*)(sk->indices);
        assert(ix <= 0x7fffffff);
        indices[i] = (int32_t)ix;
    }
}


/* USABLE_FRACTION is the maximum set load.
 * Increasing this ratio makes sets more dense resulting in more
 * collisions.  Decreasing it improves sparseness at the expense of spreading
 * indices over more cache lines and at the cost of total memory consumed.
 *
 * USABLE_FRACTION must obey the following:
 *     (0 < USABLE_FRACTION(n) < n) for all n >= 2
 *
 * USABLE_FRACTION should be quick to calculate.
 * Fractions around 1/2 to 2/3 seem to work well in practice.
 */

#define USABLE_FRACTION(n) (((n) << 1)/3)  // ratio: 2/3

/* Alternative fraction that is otherwise close enough to 2n/3 to make
 * little difference. 8 * 2/3 == 8 * 5/8 == 5. 16 * 2/3 == 16 * 5/8 == 10.
 * 32 * 2/3 = 21, 32 * 5/8 = 20.
 * Its advantage is that it is faster to compute on machines with slow division.
 * #define USABLE_FRACTION(n) (((n) >> 1) + ((n) >> 2) - ((n) >> 3))  // ratio: 5/8
 */


/* INV_USABLE_FRACTION gives the inverse of USABLE_FRACTION.
 * Used for sizing a new set to a specified number of keys.
 * 
 * NOTE: If the denominator of the USABLE_FRACTION ratio is not a power
 * of 2, must add 1 to the result of the inverse for correct sizing.
 * 
 * For example, when USABLE_FRACTION ratio = 5/8 (8 is a power of 2):
 * #define INV_USABLE_FRACTION(n) (((n) << 3)/5)  // inv_ratio: 8/5
 * 
 * When USABLE_FRACTION ratio = 5/7 (7 is not a power of 2):
 * #define INV_USABLE_FRACTION(n) ((7*(n))/5 + 1)  // inv_ratio: 7/5
 */

#define INV_USABLE_FRACTION(n) ((n) + ((n) >> 1) + 1)  // inv_ratio: 3/2


/* GROWTH_RATE. Growth rate upon hitting maximum load.
 * Currently set to used*3.
 * This means that sets double in size when growing without deletions,
 * but have more head room when the number of deletions is on a par with the
 * number of insertions.  See also bpo-17563 and bpo-33205.
 *
 * GROWTH_RATE was set to used*4 up to version 3.2.
 * GROWTH_RATE was set to used*2 in version 3.3.0
 * GROWTH_RATE was set to used*2 + capacity/2 in 3.4.0-3.6.0.
 */
#define GROWTH_RATE(s) ((s)->used*3)


static NB_SetEntry*
get_entry(NB_SetKeys *sk, Py_ssize_t idx) {
    Py_ssize_t offset;
    char *ptr;

    assert (idx < sk->size);
    offset = idx * sk->entry_size;
    ptr = sk->indices + sk->entry_offset + offset;
    return (NB_SetEntry*)ptr;
}

static void
zero_key(NB_SetKeys *sk, char *data){
    memset(data, 0, sk->key_size);
}

static void
copy_key(NB_SetKeys *sk, char *dst, const char *src){
    memcpy(dst, src, sk->key_size);
}

/* Returns -1 for error; 0 for not equal; 1 for equal */
static int
key_equal(NB_SetKeys *sk, const char *lhs, const char *rhs) {
    if ( sk->methods.key_equal ) {
        return sk->methods.key_equal(lhs, rhs);
    } else {
        return memcmp(lhs, rhs, sk->key_size) == 0;
    }
}

static char *
entry_get_key(NB_SetKeys *sk, NB_SetEntry* entry) {
    char * out = entry->key;
    assert (out == aligned_pointer(out));
    return out;
}

static void
sk_incref_key(NB_SetKeys *sk, const char *key) {
    if ( sk->methods.key_incref ) {
        sk->methods.key_incref(key);
    }
}

static void
sk_decref_key(NB_SetKeys *sk, const char *key) {
    if ( sk->methods.key_decref ) {
        sk->methods.key_decref(key);
    }
}


void
numba_setkeys_free(NB_SetKeys *sk) {
    /* Clear all references from the entries */
    Py_ssize_t i;
    NB_SetEntry *ep;

    for (i = 0; i < sk->nentries; i++) {
        ep = get_entry(sk, i);
        if (ep->hash != SETK_EMPTY) {
            sk_decref_key(sk, entry_get_key(sk, ep));
        }
    }
    /* Deallocate */
    free(sk);
}

void
numba_set_free(NB_Set *s) {
    numba_setkeys_free(s->keys);
    free(s);
}

Py_ssize_t
numba_set_length(NB_Set *s) {
    return s->used;
}

/* Allocate new set keys

Adapted from CPython's new_keys_object().
*/
int
numba_setkeys_new(NB_SetKeys **out, Py_ssize_t size, Py_ssize_t key_size) {
    Py_ssize_t usable = USABLE_FRACTION(size);
    Py_ssize_t index_size = ix_size(size);
    Py_ssize_t entry_size = aligned_size(sizeof(NB_SetEntry) + aligned_size(key_size));
    Py_ssize_t entry_offset = aligned_size(index_size * size);
    Py_ssize_t alloc_size = sizeof(NB_SetKeys) + entry_offset + entry_size * usable;

    NB_SetKeys *sk = malloc(aligned_size(alloc_size));
    if (!sk) return ERR_NO_MEMORY;

    assert ( size >= SET_MINSIZE );

    sk->size = size;
    sk->usable = usable;
    sk->nentries = 0;
    sk->key_size = key_size;
    sk->entry_offset = entry_offset;
    sk->entry_size = entry_size;

    assert (aligned_pointer(sk->indices) == sk->indices );
    /* Ensure that the method table is all nulls */
    memset(&sk->methods, 0x00, sizeof(set_type_based_methods_table));
    /* Ensure hash is (-1) for empty entry */
    memset(sk->indices, 0xff, entry_offset + entry_size * usable);

    *out = sk;
    return OK;
}


/* Allocate new set */
int
numba_set_new(NB_Set **out, Py_ssize_t key_size, Py_ssize_t size) {
    NB_SetKeys *sk;
    NB_Set *s;
    int status = numba_setkeys_new(&sk, size, key_size);
    if (status != OK) return status;

    s = malloc(sizeof(NB_Set));
    if (!s) {
        numba_setkeys_free(sk);
        return ERR_NO_MEMORY;
    }

    s->used = 0;
    s->keys = sk;
    *out = s;
    return OK;
}


/*
Adapted from CPython lookset_index().

Search index of hash table from offset of entry table
*/
static Py_ssize_t
lookset_index(NB_SetKeys *sk, Py_hash_t hash, Py_ssize_t index)
{
    size_t mask = SET_MASK(sk);
    size_t perturb = (size_t)hash;
    size_t i = (size_t)hash & mask;

    for (;;) {
        Py_ssize_t ix = get_index(sk, i);
        if (ix == index) {
            return i;
        }
        if (ix == SETK_EMPTY) {
            return SETK_EMPTY;
        }
        perturb >>= PERTURB_SHIFT;
        i = mask & (i*5 + perturb + 1);
    }
    assert(0 && "unreachable");
}

/*

Adapted from the CPython3.7 lookset().

The basic lookup function used by all operations.
This is based on Algorithm D from Knuth Vol. 3, Sec. 6.4.
Open addressing is preferred over chaining since the link overhead for
chaining would be substantial (100% with typical malloc overhead).

The initial probe index is computed as hash mod the table size. Subsequent
probe indices are computed as explained earlier.

All arithmetic on hash should ignore overflow.

The details in this version are due to Tim Peters, building on many past
contributions by Reimer Behrends, Jyrki Alakuijala, Vladimir Marangozov and
Christian Tismer.

lookset() is general-purpose, and may return SETK_ERROR if (and only if) a
comparison raises an exception.
lookset_unicode() below is specialized to string keys, comparison of which can
never raise an exception; that function can never return SETK_ERROR when key
is string.  Otherwise, it falls back to lookset().
lookset_unicode_nodummy is further specialized for string keys that cannot be
the <dummy> value.
For both, when the key isn't found a SETK_EMPTY is returned.
*/
Py_ssize_t
numba_set_lookup(NB_Set *s, const char *key_bytes, Py_hash_t hash)
{
    NB_SetKeys *sk = s->keys;
    size_t mask = SET_MASK(sk);
    size_t perturb = hash;
    size_t i = (size_t)hash & mask;

    for (;;) {
        Py_ssize_t ix = get_index(sk, i);
        if (ix == SETK_EMPTY) {
            return ix;
        }
        if (ix >= 0) {
            NB_SetEntry *ep = get_entry(sk, ix);
            const char *startkey = NULL;
            if (ep->hash == hash) {
                int cmp;

                startkey = entry_get_key(sk, ep);
                cmp = key_equal(sk, startkey, key_bytes);
                if (cmp < 0) {
                    // error'ed in comparison
                    return SETK_ERROR;
                }
                if (cmp > 0) {
                    // key is equal; retrieve the value.
                    return ix;
                }
            }
        }
        perturb >>= PERTURB_SHIFT;
        i = (i*5 + perturb + 1) & mask;
    }
    assert(0 && "unreachable");
}

Py_ssize_t
numba_set_contains(NB_Set *setp, char *key, Py_hash_t hash)
{
    Py_ssize_t ix = numba_set_lookup(setp, key, hash);
    return ix != SETK_EMPTY && ix != SETK_DUMMY;
}


/* Internal function to find slot for an item from its hash
   when it is known that the key is not present in the set.

   The set must be combined. */
static Py_ssize_t
find_empty_slot(NB_SetKeys *sk, Py_hash_t hash){
    size_t mask;
    size_t i;
    Py_ssize_t ix;
    size_t perturb;

    assert(sk != NULL);

    mask = SET_MASK(sk);
    i = hash & mask;
    ix = get_index(sk, i);
    for (perturb = hash; ix >= 0;) {
        perturb >>= PERTURB_SHIFT;
        i = (i*5 + perturb + 1) & mask;
        ix = get_index(sk, i);
    }
    return i;
}


/*
Adapted from build_indices_set().
Internal routine used by setresize() to build a hashtable of entries.
*/
void
build_indices_set(NB_SetKeys *keys, Py_ssize_t n) {
    size_t mask = (size_t)SET_MASK(keys);
    Py_ssize_t ix;
    for (ix = 0; ix != n; ix++) {
        size_t perturb;
        Py_hash_t hash = get_entry(keys, ix)->hash;
        size_t i = hash & mask;
        for (perturb = hash; get_index(keys, i) != SETK_EMPTY;) {
            perturb >>= PERTURB_SHIFT;
            i = mask & (i*5 + perturb + 1);
        }
        set_index(keys, i, ix);
    }
}

/*

Adapted from CPython setresize().

Restructure the table by allocating a new table and reinserting all
items again.  When entries have been deleted, the new table may
actually be smaller than the old one.
If a table is split (its keys and hashes are shared, its values are not),
then the values are temporarily copied into the table, it is resized as
a combined table, then the me_value slots in the old table are NULLed out.
After resizing a table is always combined,
but can be resplit by make_keys_shared().
*/
int
numba_set_resize(NB_Set *s, Py_ssize_t minsize) {
    Py_ssize_t newsize, numentries;
    NB_SetKeys *oldkeys;
    int status;

    /* Find the smallest table size > minused. */
    for (newsize = SET_MINSIZE;
         newsize < minsize && newsize > 0;
         newsize <<= 1)
        ;
    if (newsize <= 0) {
        return ERR_NO_MEMORY;
    }
    oldkeys = s->keys;

    /* Allocate a new table. */
    status = numba_setkeys_new(
        &s->keys, newsize, oldkeys->key_size
    );
    if (status != OK) {
        s->keys = oldkeys;
        return status;
    }
    // New table must be large enough.
    assert(s->keys->usable >= s->used);
    // Copy method table
    memcpy(&s->keys->methods, &oldkeys->methods, sizeof(set_type_based_methods_table));

    numentries = s->used;

    if (oldkeys->nentries == numentries) {
        NB_SetEntry *oldentries, *newentries;

        oldentries = get_entry(oldkeys, 0);
        newentries = get_entry(s->keys, 0);
        memcpy(newentries, oldentries, numentries * oldkeys->entry_size);
        // to avoid decref
        memset(oldentries, 0xff, numentries * oldkeys->entry_size);
    }
    else {
        Py_ssize_t i;
        size_t epi = 0;
        for (i=0; i<numentries; ++i) {
            /*
                ep->hash == (-1) hash means it is empty

                Here, we skip until a non empty entry is encountered.
            */
            while( get_entry(oldkeys, epi)->hash == SETK_EMPTY ) {
                epi += 1;
            }
            memcpy(
                get_entry(s->keys, i),
                get_entry(oldkeys, epi),
                oldkeys->entry_size
            );
            get_entry(oldkeys, epi)->hash = SETK_EMPTY;  // to avoid decref
            epi += 1;

        }

    }
    numba_setkeys_free(oldkeys);

    build_indices_set(s->keys, numentries);
    s->keys->usable -= numentries;
    s->keys->nentries = numentries;
    return OK;
}

static int
insertion_resize(NB_Set *s)
{
    return numba_set_resize(s, SET_GROWTH_RATE(s));
}

int
numba_set_add(
    NB_Set    *s,
    const char *key_bytes,
    Py_hash_t   hash)
{
    NB_SetKeys *sk = s->keys;

    Py_ssize_t ix = numba_set_lookup(s, key_bytes, hash);
    if (ix == SETK_ERROR) {
        // exception in key comparison in lookup.
        return ERR_CMP_FAILED;
    }

    if (ix == SETK_EMPTY) {
        /* Insert into new slot */
        Py_ssize_t hashpos;
        NB_SetEntry *ep;

        if (sk->usable <= 0) {
            /* Need to resize */
            if (insertion_resize(s) != OK)
                return ERR_NO_MEMORY;
            else
                sk = s->keys;     // reload
        }
        hashpos = find_empty_slot(sk, hash);
        ep = get_entry(sk, sk->nentries);
        set_index(sk, hashpos, sk->nentries);
        copy_key(sk, entry_get_key(sk, ep), key_bytes);
        assert ( hash != -1 );
        ep->hash = hash;

        /* incref */
        sk_incref_key(sk, key_bytes);

        s->used += 1;
        sk->usable -= 1;
        sk->nentries += 1;
        assert (sk->usable >= 0);
        return OK;
    } else {
        return OK_REPLACED;
    }
}

int
numba_set_discard(NB_Set *s, char *key_bytes, Py_hash_t hash)
{
    Py_ssize_t hashpos;
    NB_SetEntry *ep;
    NB_SetKeys *sk = s->keys;

    Py_ssize_t ix = numba_set_lookup(s, key_bytes, hash);
    hashpos = lookset_index(sk, hash, ix);
    assert(hashpos >= 0);

    if (ix == SETK_EMPTY || ix == SETK_DUMMY) {
        /* Not found */
        return OK;
    } else {
        s->used -= 1;
    }
    ep = get_entry(sk, ix);
    set_index(sk, hashpos, SETK_DUMMY);

    /* decref */
    sk_decref_key(sk, entry_get_key(sk, ep));

    /* zero the entries */
    zero_key(sk, entry_get_key(sk, ep));
    ep->hash = SETK_EMPTY; // to mark it as empty;

    return OK;
}

void
numba_set_dump(NB_Set *s) {
    long long i, j, k;
    long long size, n;
    char *cp;
    NB_SetEntry *ep;
    NB_SetKeys *sk = s->keys;

    n = s->used;
    size = sk->nentries;

    printf("set dump\n");
    printf("   key_size = %lld\n", (long long)s->keys->key_size);

    for (i = 0, j = 0; i < size; i++) {
        ep = get_entry(sk, i);
        if (ep->hash != SETK_EMPTY) {
            long long hash = ep->hash;
            printf("  key=");
            for (cp=entry_get_key(sk, ep), k=0; k < s->keys->key_size; ++k, ++cp){
                printf("%02x ", ((int)*cp) & 0xff);
            }
            printf(" hash=%llu", hash);
            printf("\n");
            j++;
        }
    }
    printf("j = %lld; n = %lld\n", j, n);
    assert(j == n);
}

size_t
numba_set_iter_sizeof() {
    return sizeof(NB_SetIter);
}

void
numba_set_iter(NB_SetIter *it, NB_Set *s) {
    it->parent = s;
    it->parent_keys = s->keys;
    it->size = s->used;
    it->pos = 0;
}

int
numba_set_iter_next(NB_SetIter *it, const char **key_ptr) {
    /* Detect set mutation during iteration */
    NB_SetKeys *sk;
    if (it->parent->keys != it->parent_keys ||
        it->parent->used != it->size) {
        return ERR_SET_MUTATED;
    }
    sk = it->parent_keys;
    while ( it->pos < sk->nentries ) {
        NB_SetEntry *ep = get_entry(sk, it->pos++);
        if ( ep->hash != SETK_EMPTY ) {
            *key_ptr = entry_get_key(sk, ep);
            return OK;
        }
    }
    return ERR_ITER_EXHAUSTED;
}


/* Allocate a new set with enough space to hold n_keys without resizes */
int
numba_set_new_sized(NB_Set **out, Py_ssize_t key_size, Py_ssize_t n_keys) {

    /* Respect SET_MINSIZE */
    if (n_keys <= USABLE_FRACTION(SET_MINSIZE)) {
        return numba_set_new(out, key_size, SET_MINSIZE);
    }

    /*  Adjust for load factor */
    Py_ssize_t size = INV_USABLE_FRACTION(n_keys) - 1;

    /* Round up size to the nearest power of 2. */
    for (unsigned int shift = 1; shift < sizeof(Py_ssize_t) * CHAR_BIT; shift <<= 1) {
        size |= (size >> shift);
    }
    size++;

    /* Handle overflows */
    if (size <= 0) {
        return ERR_NO_MEMORY;
    }

    return numba_set_new(out, key_size, size);
}


void
numba_set_set_method_table(NB_Set *s, set_type_based_methods_table *methods)
{
    memcpy(&s->keys->methods, methods, sizeof(set_type_based_methods_table));
}


#define CHECK(CASE) {                                                   \
    if ( !(CASE) ) {                                                    \
        printf("'%s' failed file %s:%d\n", #CASE, __FILE__, __LINE__);   \
        return 1;                                                       \
    }                                                                   \
}

int
numba_test_set(void) {
    NB_Set *s;
    int status;
    Py_ssize_t ix;
    Py_ssize_t usable;
    Py_ssize_t it_count;
    const char *it_key;
    NB_SetIter iter;

    puts("test_set");

    status = numba_set_new(&s, 4, SET_MINSIZE);
    CHECK(status == OK);
    CHECK(s->keys->size == SET_MINSIZE);
    CHECK(s->keys->key_size == 4);
    CHECK(ix_size(s->keys->size) == 1);
    printf("aligned_size(index_size * size) = %d\n", (int)(aligned_size(ix_size(s->keys->size) * s->keys->size)));

    printf("s %p\n", s);
    printf("s->usable = %u\n", (int)s->keys->usable);
    usable = s->keys->usable;
    printf("s[0] %d\n", (int)((char*)get_entry(s->keys, 0) - (char*)s->keys));
    CHECK ((char*)get_entry(s->keys, 0) - (char*)s->keys->indices == s->keys->entry_offset);
    printf("s[1] %d\n", (int)((char*)get_entry(s->keys, 1) - (char*)s->keys));
    CHECK ((char*)get_entry(s->keys, 1) - (char*)s->keys->indices == s->keys->entry_offset + s->keys->entry_size);

    ix = numba_set_lookup(s, "bef", 0xbeef);
    printf("ix = %d\n", (int)ix);
    CHECK (ix == SETK_EMPTY);

    // insert 1st key
    status = numba_set_add(s, "bef", 0xbeef);
    CHECK (status == OK);
    CHECK (s->used == 1);
    CHECK (s->keys->usable == usable - s->used);

    // insert same key
    status = numba_set_add(s, "bef", 0xbeef);
    CHECK (status == OK_REPLACED);
    CHECK (s->used == 1);
    CHECK (s->keys->usable == usable - s->used);

    // insert 2nd key
    status = numba_set_add(s, "beg", 0xbeef);
    CHECK (status == OK);
    CHECK (s->used == 2);
    CHECK (s->keys->usable == usable - s->used);

    // insert 3rd key
    status = numba_set_add(s, "beh", 0xcafe);
    CHECK (status == OK);
    CHECK (s->used == 3);
    CHECK (s->keys->usable == usable - s->used);

    status = numba_set_add(s, "bef", 0xbeef);
    CHECK (status == OK_REPLACED);
    CHECK (s->used == 3);
    CHECK (s->keys->usable == usable - s->used);

    // insert 4th key
    status = numba_set_add(s, "bei", 0xcafe);
    CHECK (status == OK);
    CHECK (s->used == 4);
    CHECK (s->keys->usable == usable - s->used);

    // insert 5th key
    status = numba_set_add(s, "bej", 0xcafe);
    CHECK (status == OK);
    CHECK (s->used == 5);
    CHECK (s->keys->usable == usable - s->used);

    // insert 6th key & triggers resize
    status = numba_set_add(s, "bek", 0xcafe);
    CHECK (status == OK);
    CHECK (s->used == 6);
    CHECK (s->keys->usable == USABLE_FRACTION(s->keys->size) - s->used);

    // Dump
    numba_set_dump(s);

    // Make sure everything are still in there
    ix = numba_set_lookup(s, "bef", 0xbeef);
    CHECK (ix >= 0);

    ix = numba_set_lookup(s, "beg", 0xbeef);
    CHECK (ix >= 0);

    ix = numba_set_lookup(s, "beh", 0xcafe);
    printf("ix = %d\n", (int)ix);
    CHECK (ix >= 0);

    ix = numba_set_lookup(s, "bei", 0xcafe);
    CHECK (ix >= 0);

    ix = numba_set_lookup(s, "bej", 0xcafe);
    CHECK (ix >= 0);

    ix = numba_set_lookup(s, "bek", 0xcafe);
    CHECK (ix >= 0);

    // Test delete
    ix = numba_set_lookup(s, "beg", 0xbeef);
    CHECK (ix >= 0);

    status = numba_set_discard(s, "beg", 0xbeef);
    CHECK (status == OK);

    ix = numba_set_lookup(s, "bef", 0xbeef);
    CHECK (ix >= 0);
    ix = numba_set_lookup(s, "beh", 0xcafe);
    CHECK (ix >= 0);

    // Test iterator
    CHECK( s->used > 0 );
    numba_set_iter(&iter, s);
    it_count = 0;
    while ( (status = numba_set_iter_next(&iter, &it_key)) == OK) {
        it_count += 1;  // valid items
        CHECK(it_key != NULL);
    }

    CHECK(status == ERR_ITER_EXHAUSTED);
    CHECK(s->used == it_count);

    numba_set_free(s);

    /* numba_set_new_sized() */

    Py_ssize_t target_size;
    Py_ssize_t n_keys;

    // Test if minsize set returned with n_keys=0
    target_size = SET_MINSIZE;
    n_keys = 0;

    numba_set_new_sized(&s, 1, n_keys);
    CHECK(s->keys->size == target_size);
    CHECK(s->keys->usable == USABLE_FRACTION(target_size));
    numba_set_free(s);

    // Test sizing at power of 2 boundary
    target_size = SET_MINSIZE * 2;
    n_keys = USABLE_FRACTION(target_size);

    numba_set_new_sized(&s, 1, n_keys);
    CHECK(s->keys->size == target_size);
    CHECK(s->keys->usable == n_keys);
    numba_set_free(s);

    target_size *= 2;
    n_keys++;

    numba_set_new_sized(&s, 1, n_keys);
    CHECK(s->keys->size == target_size);
    CHECK(s->keys->usable > n_keys);
    CHECK(s->keys->usable == USABLE_FRACTION(target_size));
    numba_set_free(s);

    return 0;

}

#undef CHECK

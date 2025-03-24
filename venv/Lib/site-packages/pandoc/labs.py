"""
Design
================================================================================

TODO: first step: forget about lazyness and performance.

Fundamental object: Query. NO, see below

  - Query defined independently of the root(s) it will be applied to.

  - Query is a callable ? Once defined, applied with `query(doc)` for
    example ? Or stuff = query(doc) is the aim is some extraction ?
    Or for item in query(doc): do_stuff_with(item). Root(s?) can be
    a Pandoc item or a list or tuple of stuff, etc.
    
    Could also envision query methods directly on Pandoc Elements ...
    dunno, have to think of it. Would be convenient (shorter code), 
    but messy (Elements are not "algebraic datatypes" anymore).

  - Query construction uses chaining (fluent) interface : chaining of
    transforms, chaining of filters, etc.

  - Query definition is lazy, does not compute anything: just abstract
    description of what we intend to do.

  - Query applied to a root : by default all items (in doc order)

  - Query filters: may reduce the list of items. Chaining is an "AND" thus
    a single-filter should be likely "OR"-based. E.g.
    query.filter((Para, Plain)) means keep Para OR Plain.
    Possible filters : filter(types), filter(predicate), has_attr(...),
    id_is(id=...), etc. Need to have stuff to match the existence of a feature
    has well as its value (ex: has id vs id is "main"). Also, class matching,
    etc. Essentially, many syntaxic sugar on top of "filter(predicate)".
    May need a meta stuff: .either(pred1, pred2, ...) ? Have a look at
    shader graph to see how stuff is forked/join with such interfaces.

  - Navigation: .parent, .children, .next_sibling, get_item[i]. Extract a 
    general "navigate" function such that all that stuff is merely syntaxic
    sugar ?

  - Mutation: .remove(), .replace_with(stuff), etc. Probably complex to get it
    right here. All this stuff should be TERMINAL : if any, nothing can be
    chained after that. Also: set attributes, ids, etc.
    Again, can we come up with a single, root, general mutation operation ?
    Or is it too complex here?
    Apply mutation in reverse order by default?

  - Iteration ? Query as a container "removes" path info by default so that
    we can deal with a list of "document items", make comprehensions, etc.


New Take
--------------------------------------------------------------------------------

  - We have PREDICATES first: stuff that takes an elt and returns a bool.

  - We have a list of "fancy", parametrized predicates (match the id, match
    the type, etc.).

  - Predicates are functions (elt -> True) that have some extra support for
    |, & and ~ (and/or the or_, and_ and not_ methods?). But not that a
    "filter" method (or __call__) on Elements/Collections could/does work
    as a "and", so multiple args in () could work has a "or" and we merely
    need a not to be complete (we have a conjonctive normal form).

  - On an elt or a collection of elements, we can either query using a 
    predicate or filter using a predicate. Query is get everything + filter.
    Have a look at rethinkdb for the terminology? And/or pandas?
    Terminology: "run" a query. Rethinkdb uses a chain with "." + final "run()"
    to get the results of a query that runs.
    In db-land, there is an iteration by default?
    But not a recursive structure ... JQuery would be a better reference I 
    guess, because of its tree structure.

    The find vs filter distinction is more JQuery-ish.


  - Anyway, we have a elts + predicates -> elts function, that we should be
    able to chain.

  - Also, transforms (elts -> elts, such a .children) and apply (returns
    anything).



"""

# 🚧: Document that multiple arguments + `not_` + chaining calls allows to
#     express arbitrary boolean logic in [conjunctive normal form][CNF].
#
#     [CNF]: https://en.wikipedia.org/wiki/Conjunctive_normal_form

# Python Standard Library
pass

# Third-Party Libraries
pass

# Pandoc
import pandoc


# Logic
# ------------------------------------------------------------------------------

# 🤔 Consider support for:
#
#      - Typing.Union (supported by | in Python >=3.10) ?
#        We probably should support whatever isinstance supports.
#        That means tuples too? Yes, see the concept of "classinfo" here:
#        <https://docs.python.org/3/library/functions.html#isinstance>.
#        Note: tuples can be recursively !!!
#        Note: if we support tuples, we need to unpack them (do we? We may
#        also consider that a tuple means a classinfo and defer the error
#        management to isinstance).
#
def to_function(*predicates):
    """
    🚧 Document:

        - pass-through for function

        - type to isinstance(., type) conversion

        - vectorization by disjunction of the predicates
    """
    if not predicates:
        predicates = [lambda _: True]
    if len(predicates) == 1:
        predicate = predicates[0]
        if isinstance(predicate, type):
            return lambda elt: isinstance(elt, predicate)
        elif callable(predicate):
            return predicate
        else:
            error = "predicate should be a type or a function"
            error += f", not {predicate!r}"
            raise TypeError(error)
    else:
        return lambda elt: any(to_function(p)(elt) for p in predicates)


def not_(predicate):
    predicate = to_function(predicate)
    return lambda elt: not predicate(elt)


# Queries & Results
# ------------------------------------------------------------------------------
def query(root):
    return Query([(root, [])])


def _getitem(sequence, indices):
    if not hasattr(sequence, "__getitem__") or isinstance(sequence, str):
        raise TypeError()
    if isinstance(sequence, dict):
        sequence = list(sequence.items())
    return sequence[indices]


# 💭 "Query" probably inappropriate since the queries are the method calls,
#    not the object itself.
#    The object itself would be a "collection", "elements",
#    "forest", "results", that kind of thing. "Library", "Fragments", anything
#     like that is cool. I like "Library", as a collection of documents
#     (and doc fragments). The branding/concept would be quite nice.
#    "Corpus" would also work pretty well.


class Query:
    def __init__(self, results):
        self._elts = []
        if isinstance(results, tuple):
            results = [results]
        if isinstance(results, Query):
            self._elts.extend(results._elts)  # Mmmm ownership issues. Copy?
        else:  # "raw results": list of (elt, path)
            self._elts.extend(results)

    # ℹ️ The call `find(object)` is public and equivalent to `_iter()`.
    def _iter(self):
        results = []
        for root, root_path in self._elts:
            for elt, path in pandoc.iter(root, path=True):
                path = root_path + path
                results.append((elt, path))
        return Query(results)

    # 🤔 Think of a rename given that we now can expose this as a property
    #    (optionally restricted with a call). We can keep find and the
    #    "functionally flavor", but a more content-oriented alias would
    #    be nice (descendants? But we also return the node itself.
    #    Subtree? Tree? Contents?)
    #
    # 🤔 Replace "find" with "search"?
    def find(self, *predicates):
        return self._iter().filter(*predicates)

    def filter(self, *predicates):
        predicate = to_function(*predicates)
        results = []
        for elt, path in self._elts:
            if predicate(elt):
                results.append((elt, path))
        return Query(results)

    # ✨ This is sweet! The intended usage is `.property(test)`,
    #    which emulates the jquery API where functions can restrict the match.
    #    It also allows us call the "functions" without parentheses
    #    when no restriction is needed.
    def __call__(self, *predicates):
        return self.filter(*predicates)

    def get_children(self):
        # 🚧 TODO: use _getitem
        results = []
        for elt, path in self._elts:
            if isinstance(elt, dict):
                for i, child in enumerate(elt.items()):
                    child_path = path.copy() + [(elt, i)]
                    results.append((child, child_path))
            elif hasattr(elt, "__iter__") and not isinstance(elt, str):
                for i, child in enumerate(elt):
                    child_path = path.copy() + [(elt, i)]
                    results.append((child, child_path))
        return Query(results)

    children = property(get_children)

    def get_child(self, i):
        results = []
        for elt, elt_path in self._elts:
            children = []
            if isinstance(elt, pandoc.types.String) or not hasattr(elt, "__iter__"):
                children = []
            elif isinstance(elt, dict):
                children = list(elt.items())
            else:
                children = elt[:]

            if isinstance(i, int):
                try:
                    if i < 0:
                        i += len(elt)
                    child = children[i]
                    results.append((child, elt_path.copy() + [(elt, i)]))
                except IndexError:
                    pass
            elif isinstance(i, slice):
                indices = range(len(children))[i]
                for index in indices:
                    child = children[index]
                    results.append((child, elt_path.copy() + [(elt, index)]))
        return Query(results)

    def get_parent(self):
        results = []
        for _, path in self._elts:
            if path != []:
                results.append((path[-1][0], path[:-1]))
        return Query(results)

    parent = property(get_parent)

    def get_next(self):
        results = []
        for elt, elt_path in self._elts:
            try:  # try to go inside elt first
                first_child = _getitem(elt, 0)
                results.append((first_child, elt_path + [(elt, 0)]))
                continue
            except (IndexError, TypeError):
                pass

            q = Query([(elt, elt_path)])
            # Try the next sibling, if it doesn't exist the parent next sibling, etc.
            while True:
                next_ = q.next_sibling
                if next_:
                    results.append(next_._elts[0])
                    break
                else:
                    parent = q.parent
                    if not parent:  # we're back at the root
                        break
                    else:
                        q = parent

        return Query(results)

    next = property(get_next)

    # 🚧 TODO: previous

    def get_previous(self):
        results = []
        for elt, elt_path in self._elts:

            # TODO: try deeper in the previous sibling first
            #       or the previous sibling
            #       or of there is no previous sibling, the parent

            previous_sibling = Query([(elt, elt_path)]).previous_sibling
            if previous_sibling:
                last_child = previous_sibling
                while child := last_child.get_child(-1):
                    last_child = child
                    pass
                results.append(last_child._elts[0])
            elif parent := Query([(elt, elt_path)]).parent:
                results.append(parent._elts[0])
            else:
                pass

        return Query(results)

    previous = property(get_previous)

    def get_next_sibling(self):
        indices = [path[-1][1] for elt, path in self._elts if path != []]
        results = []
        for (parent_elt, parent_path), index in zip(self.parent._elts, indices):
            try:  # 🚧 TODO: adaptation of [] for dicts and strings (use _getitem).
                next_element = parent_elt[index + 1]
                results.append(
                    (next_element, parent_path.copy() + [(parent_elt, index + 1)])
                )
            except IndexError:
                pass
        return Query(results)

    next_sibling = property(get_next_sibling)

    def get_previous_sibling(self):
        indices = [path[-1][1] for elt, path in self._elts if path != []]
        results = []
        for (parent_elt, parent_path), index in zip(self.parent._elts, indices):
            if index > 0:
                try:  # 🚧 TODO: adaptation of [] for dicts and strings (& factor out).
                    next_element = parent_elt[index - 1]
                    results.append(
                        (next_element, parent_path.copy() + [(parent_elt, index - 1)])
                    )
                except IndexError:
                    pass
        return Query(results)

    previous_sibling = property(get_previous_sibling)

    # Query container (hide path info). Or maybe not? Be more explicit?
    # --------------------------------------------------------------------------
    def __len__(self):
        return len(self._elts)

    def __bool__(self):
        return len(self._elts) != 0

    def __getitem__(self, i):
        return Query(self._elts[i])

    def __iter__(self):  # unwrap or not? Mmmm maybe no.
        return (elt for elt, _ in self._elts)

    def __repr__(self):
        return "\n".join("- " + repr(elt) for elt, path in self._elts)

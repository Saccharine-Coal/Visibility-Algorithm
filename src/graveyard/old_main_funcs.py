
def is_nested(iterable) -> bool:
    def is_iter(element) -> bool:
        if isinstance(element, list):
            return True
        if isinstance(element, tuple):
            return True
        if isinstance(element, np.ndarray):
            return True
        return False
    return any(is_iter(val) for val in iterable)
def trim_arcs(arcs, caster):
    # trim arcs that are covered by other arcs
    orig_size = len(arcs)
    trimmed_arcs = arcs.copy()
    for i, arc in enumerate(arcs):
        for j, test_arc in enumerate(arcs):
            if i != j:
                # do not need to test equal arcs
                s1, e1 = tuple(x.as_polar(caster.xy) for x in test_arc)[:]
                s2, e2 = tuple(x.as_polar(caster.xy) for x in arc)[:]
                start, end = s2[1], e2[1]
                if polar.is_in_arc(s1[1], start, end) and polar.is_in_arc(e1[1], start, end):
                    # is in arc, need to check if the radius is larger
                    if s1[0] > s2[0]  and e1[0] > e2[0]:
                        trimmed_arcs[j] = None
    _arcs = tuple(filter(lambda x: x is not None, trimmed_arcs))
    return _arcs

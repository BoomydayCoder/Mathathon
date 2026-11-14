"""Exhaustive verification of typed, slot-preserving cycle splicing."""
from itertools import permutations, product
from math import factorial, prod


def run_case(p, ks):
    e = len(ks)
    ls = [int(i == 0)+p*ks[(i-1) % e]-ks[i] for i in range(e)]
    assert min(ls) >= 0 and ks[0] > 0
    types, internal = [], []
    for i, count in enumerate(ks):
        types.extend([i]*count)
        internal.extend([True]*count)
    for i, count in enumerate(ls):
        types.extend([i]*count)
        internal.extend([False]*count)
    vertices = range(len(types))
    roots = [v for v in vertices if internal[v] and types[v] == 0]
    r = max(roots)
    slots = [[(v, j) for v in vertices if internal[v] and (types[v]+1) % e == i
              for j in range(p)] for i in range(e)]

    def assignments(root):
        groups = [[v for v in vertices if types[v] == i and v != root]
                  for i in range(e)]
        assert all(len(a) == len(b) for a, b in zip(groups, slots))
        for perm_groups in product(*(permutations(g) for g in groups)):
            yield {v: slot for group, slot_group in zip(perm_groups, slots)
                   for v, slot in zip(group, slot_group)}

    def cycles(parent, root):
        done = {root}
        out = []
        for start in vertices:
            path, index = [], {}
            v = start
            while v not in done and v not in index:
                index[v] = len(path)
                path.append(v)
                v = parent[v][0]
            if v in index:
                out.append(path[index[v]:])
            done.update(path)
        return out

    def forward(parent):
        cs = sorted((max(v for v in cy if types[v] == 0)
                     for cy in cycles(parent, r)), reverse=True)
        out = dict(parent)
        root = r
        for c in cs:
            slot = out.pop(c)
            out[root] = slot
            root = c
        return out, root

    def inverse(parent, root):
        reverse_path = [r]
        while reverse_path[-1] != root:
            reverse_path.append(parent[reverse_path[-1]][0])
        records, greatest = [], -1
        for v in reversed(reverse_path):
            if types[v] == 0 and v > greatest:
                records.append(v)
                greatest = v
        assert records[-1] == r and records[0] == root
        out = dict(parent)
        for c, next_record in zip(records, records[1:]):
            slot = parent[next_record]
            del out[next_record]
            out[c] = slot
        return out

    def key(parent):
        return tuple(parent.get(v, (-1, -1)) for v in vertices)

    images = set()
    max_cycles = 0
    for parent in assignments(r):
        max_cycles = max(max_cycles, len(cycles(parent, r)))
        tree, root = forward(parent)
        assert not cycles(tree, root)
        assert set(tree.values()) == set(parent.values())
        assert inverse(tree, root) == parent
        image = key(tree)
        assert image not in images
        images.add(image)
    expected = prod(factorial(p*k) for k in ks)
    assert len(images) == expected
    all_trees = 0
    for root in roots:
        for parent in assignments(root):
            if not cycles(parent, root):
                all_trees += 1
                assert key(parent) in images
                recovered = inverse(parent, root)
                assert forward(recovered) == (parent, root)
    assert all_trees == expected
    print(f'PASS p={p}, k={ks}, l={ls}: {expected} assignments/trees; max {max_cycles} cycles')
    return expected


total = 0
for p, ks in [(2, (3,)), (3, (2,)), (2, (2, 1)), (2, (3, 2)),
              (2, (2, 2, 1)), (3, (2, 1)), (2, (1, 0, 0))]:
    total += run_case(p, ks)
print(f'All {total} fixed-root assignments verified bijectively, including surjectivity')

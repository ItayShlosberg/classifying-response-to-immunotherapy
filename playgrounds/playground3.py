


l = [['CD3E'], ['CD28'], ['CD27'], ['CCR7'], ['ID3'], ['TCF7'], ['HNF1A'], ['BCL2L11'], ['CD8A', 'CD8B', 'CD4']]


def rec(leftover, all_combinations=[[]]):
    def meiosis_of_combination(and_markers, or_markers):
        new_combinations = []
        for or_marker in or_markers:
            for and_marker in and_markers:
                new_combinations.append(and_marker + [or_marker])
        return new_combinations

    if len(leftover) == 0:
        return all_combinations
    all_currrent_markers = leftover[0]
    del leftover[0]
    all_combination = meiosis_of_combination(all_combinations, all_currrent_markers)
    return rec(leftover, all_combination)

r = rec(l)
print(r)
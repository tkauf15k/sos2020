_hierarchy = {
    ord('A'): {ord('R'), ord('W'), ord('M')},
    ord('C'): {ord('Y'), ord('S'), ord('M')},
    ord('G'): {ord('R'), ord('S'), ord('K')},
    ord('T'): {ord('Y'), ord('W'), ord('K')},

    ord('R'): {ord('D'), ord('V')},
    ord('Y'): {ord('B'), ord('H')},
    ord('S'): {ord('B'), ord('V')},
    ord('W'): {ord('D'), ord('H')},
    ord('M'): {ord('H'), ord('V')},
    ord('K'): {ord('B'), ord('D')},

    ord('B'): {ord('N')},
    ord('D'): {ord('N')},
    ord('H'): {ord('N')},
    ord('V'): {ord('N')},
    ord('-'): {ord('N')},
}

_levels = {
    ord('A'): 0, ord('C'): 0, ord('G'): 0, ord('T'): 0,
    ord('R'): 1, ord('Y'): 1, ord('S'): 1, ord('W'): 1, ord('M'): 1, ord('K'): 1,
    ord('B'): 2, ord('D'): 2, ord('H'): 2, ord('V'): 2, ord('-'): 2,
    ord('N'): 3
}

from collections import Counter


class Anchor:
    def __init__(self, scaffold, pos, ID):
        self.scaffold = scaffold
        if pos:
            self.pos = int(pos)
        else:
            self.pos = None
        self.ID = ID

    def __str__(self):
        return '\t'.join(map(str,[self.scaffold, self.pos, self.ID]))

empty_anchor = Anchor(None, None, None)

class Breakpoint:
    def __init__(self, left, right, specie):
        self.left = left
        self.right = right
        self.specie = specie

    def __str__(self):
        return "\t".join(map(str, [self.left, self.right]))


def parse_counted(path):
    counted = Counter()
    with open(path) as f:
        for line in f:
            line = line.strip().split()
            counted[line[3][3:]] += 1
    return set(filter(lambda x: counted[x] == 1, counted))
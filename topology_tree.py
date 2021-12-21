"""=================================================================================================
Construct a tree of stems for the in the topology

Michael Gribskov     17 December 2021
================================================================================================="""
from topology import Topology, RNAstructure, Stem


class TopologyNode:
    """=============================================================================================
    one node in the topology tree

    ============================================================================================="""

    def __init__(self, stem):
        """-----------------------------------------------------------------------------------------
        Constructor
        -----------------------------------------------------------------------------------------"""
        self.stem = stem
        self.parent = None
        self.children = []

    def contains(self, node):
        """-----------------------------------------------------------------------------------------
        returns true if the stem in node is entirely contained in the stem represented by self
        :param stem: Stem object
        :return: boolean, True if stem is contained
        -----------------------------------------------------------------------------------------"""
        stem = node.stem
        # if stem.lbegin >= self.stem.lbegin and stem.rend <= self.stem.rend:
        if stem.lbegin >= self.stem.lend and stem.rend <= self.stem.rbegin:
            return True
        return False

    def __str__(self):
        """-----------------------------------------------------------------------------------------

        :return:
        -----------------------------------------------------------------------------------------"""
        return f'{self.stem.name}: {self.stem.lbegin}, {self.stem.rend}'


####################################################################################################

class TopologyTree:
    """=============================================================================================
    Build a tree of stems such that the children are entirely contained within their parents.
    Each stem is nested as deeply as possible, i.e., is not a child of both its parent and
    grandparent.

    TODO: what happens when a stem is inside two parents?
    ============================================================================================="""

    def __init__(self, topology):
        """-----------------------------------------------------------------------------------------
        Constructor

        :param topology: Topology object (required)
        -----------------------------------------------------------------------------------------"""
        self.topology = topology

        # initialize the tree with a dummy stem that includes the entire range of the tree
        root = Stem()
        root.lbegin = root.lend = 1
        root.rbegin = root.rend = self.topology.sequence_length
        root.name = 'root'
        self.tree = TopologyNode(root)

        for stem in sorted(self.topology.stem_list,
                           key=lambda s: (s.rend - s.lbegin),
                           reverse=True):
            s = TopologyNode(stem)
            nodestack = [self.tree]
            while nodestack:
                node = nodestack.pop()
                if node.contains(s):
                    print(f'stem {str(s)} is contained in {str(node)}')
                    contained_child = []
                    for child in node.children:
                        if child.contains(s):
                            contained_child.append(child)

                    if contained_child:
                        nodestack = contained_child
                    else:
                        node.children.append(s)

    def dump(self):
        """-----------------------------------------------------------------------------------------
        Print out the entire tree
        :return: True
        -----------------------------------------------------------------------------------------"""
        nodestack = [self.tree]
        seen = []
        while nodestack:
            node = nodestack.pop()
            if node.children and not node in seen:
                print(f'\nparent={str(node)}\t{node.stem.lvienna}\t{node.stem.rvienna}')
                seen.append(node)
                for child in node.children:
                    print(f'\tchild={str(child)}')
                    nodestack.append(child)

        return True

    def merge1(self, maxgap=3):
        """-----------------------------------------------------------------------------------------
        Merge nodes that are only children of a parent.
        :return:
        -----------------------------------------------------------------------------------------"""
        nodestack = [self.tree]
        while nodestack:
            node = nodestack.pop()
            nchildren = len(node.children)
            if nchildren == 0:
                continue

            if nchildren == 1:
                child = node.children[0]
                if child.stem.lbegin - node.stem.lend - 1 < maxgap and \
                        child.stem.rend - node.stem.rbegin - 1 < maxgap:
                    node.stem.lvienna, node.stem.rvienna = TopologyTree.merge_vienna(node.stem,
                                                                                     child.stem)
                    node.stem.lend = max(node.stem.lend, child.stem.lend)
                    node.stem.rbegin = min(node.stem.rbegin, child.stem.rbegin)
                    node.children = child.children
                    nodestack.append(node)
                else:
                    nodestack.append(child)

            else:
                for child in node.children:
                    nodestack.append(child)

        return True

    def to_stemlist(self):
        """-----------------------------------------------------------------------------------------
        Convert to list of stems
        :return:
        -----------------------------------------------------------------------------------------"""
        stemlist = []
        nodestack = [self.tree]
        seen = []
        while nodestack:
            node = nodestack.pop()
            stemlist.append(node.stem)
            for child in node.children:
                nodestack.append(child)

        return stemlist

    @staticmethod
    def merge_vienna(outer, inner):
        """-----------------------------------------------------------------------------------------
        based on the existing dot bracket structure, merge the strings for the outer and inner stems

        :param outer: Stem, outer stem
        :param inner: Stem, inner stem
        :return: str, str - left and right vienna strings
        -----------------------------------------------------------------------------------------"""
        lvienna = outer.lvienna
        rvienna = inner.rvienna
        pos = outer.lend + 1
        while pos < inner.lbegin:
            lvienna += '.'
            pos += 1
        lvienna += inner.lvienna

        pos = inner.rend + 1
        while pos < outer.rbegin:
            rvienna += '.'
            pos += 1
        rvienna += outer.rvienna

        return lvienna, rvienna


# --------------------------------------------------------------------------------------------------
# testing
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    import sys

    test = RNAstructure()
    ddG = 1.75
    test.CTRead('data/mr_s129.fold.ct', ddG)
    test.XIOSwrite(sys.stdout)

    tree = TopologyTree(test)
    tree.dump()
    print('\n\n')
    tree.merge1(maxgap=2)
    tree.dump()

    test.stemlist = tree.to_stemlist()
    test.XIOSwrite(sys.stdout)

    exit(0)

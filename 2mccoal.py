#!/usr/bin/env python
# Takes a tree and population size and generates the content of  MCcoal.ctl file
#   for use with Yang and Rannala's MCcoal simulator of a contained coalescent.
# All populations have the same population size

import sys
import dendropy
from cStringIO import StringIO

def mccoalNewick(o, nd):
    mccoalNewick_r(o, nd)
    o.write(";")

def mccoalNewick_r(o, nd):
    child_nodes = nd.child_nodes()
    if child_nodes:
        o.write("(")
        for nn, n in enumerate(child_nodes):
            if nn > 0:
                o.write(", ")
            mccoalNewick_r(o, n)
        o.write(") : %f # %f " % (nd.depth, nd.pop_size))
    else:
        o.write("%s # %f" % (nd.taxon.label, nd.pop_size))

o = StringIO()
if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-s", "--seed", dest="seed", default=-1,
        type="int",
        help="The random number generator seed")
    parser.add_option("-p", "--pop-size", dest="p", default=1.0,
        type="float",
        help="The effective population size of each population")
    parser.add_option("-n", "--num-per-pop", dest="n", default=1,
        type="int",
        help="The number of individuals sampled from each population")
    (options, args) = parser.parse_args()
    if len(args) > 1:
        sys.exit("At most one argument (a newick tree string with branch lengths) can be specified")


    o = sys.stdout
    o.write("""SimulatedData.txt
%d
""" % options.seed)
    if len(args) == 1:
        newick = args[0]
    else:
        newick = sys.stdin.readline()
    tree =  dendropy.Tree.get_from_string(newick, schema="NEWICK")

    for node in tree.postorder_node_iter():
        ch = node.child_nodes()
        if len(ch) == 0:
            node.depth = 0.0
        else:
            first_child = ch[0]
            node.depth = first_child.depth + first_child.edge.length
            for nnd in ch[1:]:
                ocnd = nnd.depth + nnd.edge.length
                if abs(node.depth - ocnd) > 0.00001:
                    sys.exit("Tree is not ultrametric")
        node.pop_size = options.p
    l = tree.leaf_nodes()
    sample_n = [str(options.n)] * len(l)
    o.write("%d %s\n  %s\n" % (len(l),
                             " ".join([n.taxon.label for n in l]),
                             " ".join(sample_n)))
    mccoalNewick(o, tree.seed_node)

    o.write("""
//end of file
""")




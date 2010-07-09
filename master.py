#!/usr/bin/env python
import sys
import os
import random
from subprocess import Popen, PIPE
from cStringIO import StringIO
import dendropy

CULL_THE_YOUNGEST_CHERRY = True

#from dendropy.treecalc import pybus_harvey_gamma
def db_show_tree(t, prefix="tmp"):
    newick = str(t)
    pd = os.path.join(DENDROPY_USER_DIR, "tmp")
    if not os.path.exists(pd):
        os.makedirs(pd)
    t = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', prefix=prefix, delete=False, dir=pd)
    n = t.name
    t.write(newick)
    t.close()
    subprocess.Popen(["open", "-a", "/Applications/FigTree v1.2.1.app/", n])
    tflf = open(TMP_FILE_LIST_FILENAME, "a")
    tflf.write("%s\n" % n)
    tflf.close()



_logger_initialized = False

def get_logger(s):
    """Wrapper around logging.getLogger that make sure that the dendropy
    logging configuration file is read (or a default is applied)
    """
    global _logger_initialized
    import logging
    msg = ""
    if not _logger_initialized:
        import logging.config
        level = logging.INFO
        try:
            lfe =  os.environ.get("MASTER_PY_LOGGING_LEVEL")
            level = eval(lfe)
        except:
            pass
        logger = logging.getLogger()
        logger.setLevel(level)
        ch = logging.StreamHandler()
        ch.setLevel(level)
        default_fmt_str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        formatter = logging.Formatter(default_fmt_str)
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        _logger_initialized = True
    l = logging.getLogger(s)
    if msg:
        l.debug(msg)
    return l
_LOG = get_logger("master_py")

sp_tree_simualtors = ["bdtree"]
gene_tree_simualtors = ["mccoal"]

_RNG = random.Random()
def ext_seed(rng):
    return rng.getrandbits(32)

def exit(estring):
    sys.exit("master.py: %s" % estring)


def getProc(exe, args, stdout=None, stdin=None):
    invoc = [exe] + [str(i) for i in args]
    try:
        p = Popen(invoc, stdout=stdout, stdin=stdin)
    except:
        exit("%s was not on the PATH" % exe)
    return p


def bdtree_sim_sp_tree(options):
    args = ["--%s=%s" % (i, str(eval("options.%s" % i ))) for i in ["leaves", "sample", "birth", "death", "ultimate"]]
    full_args = ["--seed=%d" % ext_seed(_RNG)] + args
    p = getProc("bdtree.py", full_args, PIPE)
    v = p.communicate()[0]
    if p.wait() != 0:
        exit("Error in:\nbdtree %s" % " ".join(full_args))
    fn = ".bd_outgroups.txt"
    try:
        f = open(fn, "rU")
    except:
        exit("Error in bdtree.  Expected file %s was not produced" % fn)
    outgroups = [i.strip() for i in f]
    f.close()
    os.remove(fn)
    sp_tree = dendropy.Tree.get_from_string(v, schema="NEWICK")
    return sp_tree, outgroups

def rtree_outgroup_labels(tree):
    """Takes a tree (which will be treated as rooted), and returns a list of labels
    of the nodes that would serve as the "outgroup" if you were to define the
    largest clade in the tree to be the "ingroup".

    Adds "n_leaves_under" and "in_biggest" attributes to the nodes in the tree."""
    node = None
    # add an n_leaves_under attribute
    for node in tree.postorder_node_iter():
        e = node.edge
        p = getattr(e, "tail_node", None)
        if p:
            p.n_leaves_under = getattr(p, "n_leaves_under", 0) +  getattr(node, "n_leaves_under", 1)

    # find the child of the root with the largest number of descendants
    seed_node = tree.seed_node
    ch = seed_node.child_nodes()
    f = ch[0]
    f.in_biggest = False
    biggest_clade, bc_size = f, getattr(f, "n_leaves_under", 1)
    for nd in ch[1:]:
        nk = getattr(nd, "n_leaves_under", 1)
        if nd > bc_size:
            biggest_clade, bc_size = nd, nk
        nd.in_biggest = False
    # Mark the biggest clade, and accumulate out all unmarked leaf names
    biggest_clade.in_biggest = True
    outgroup_labels = []
    for node in tree.preorder_node_iter():
        par = node.parent_node
        if node == seed_node or par == seed_node:
            continue
        node.in_biggest = par.in_biggest
        if (not node.in_biggest) and (not node.child_nodes()):
            outgroup_labels.append(node.label)
    return outgroup_labels


def mccoal_sim_gene_tree(sp_tree, options):
    ctlFile = open("MCcoal.ctl", "w")
    full_args = ["--seed=%d" % ext_seed(_RNG), "--pop-size=%f" % options.pop_size, str(sp_tree)]
    p = getProc("2mccoal.py", full_args, ctlFile)
    if p.wait() != 0:
        exit("Error in:\n2mccoal.py %s" % " ".join(full_args))
    ctlFile.close()
    p = getProc("MCcoal", [], PIPE, PIPE)
    w = p.communicate("""1 1
""")[0]
    fn = "out.trees"
    if p.wait() != 0 or not os.path.exists(fn):
        exit("Error in:\n2mccoal.py %s" % " ".join(full_args))
    gene_tree_file = open(fn, "rU")
    gene_tree_newick = gene_tree_file.readline().strip()

    culled_gene_tree =  dendropy.Tree.get_from_string(gene_tree_newick, schema="NEWICK")
    non_culled_gene_tree =  dendropy.Tree.get_from_string(gene_tree_newick, schema="NEWICK")
    if CULL_THE_YOUNGEST_CHERRY:
        mrt = float('inf')
        youngest_mrca = None
        for nd in culled_gene_tree.leaf_iter():
            e = nd.edge
            if e.length < mrt:
                mrt = e.length
                youngest_mrca = nd.parent_node
        assert(youngest_mrca is not None)
        c = youngest_mrca.child_nodes()
        for child in c:
            youngest_mrca.remove_child(child)
        for nd in culled_gene_tree.leaf_iter():
            if nd is not youngest_mrca:
                e = nd.edge
                e.length = e.length - mrt

    non_culled_outgroup_labels = rtree_outgroup_labels(non_culled_gene_tree)
    outgroup_labels = rtree_outgroup_labels(culled_gene_tree)

    return non_culled_gene_tree, non_culled_outgroup_labels, culled_gene_tree, outgroup_labels


def divers_rate_stats(options, tree, tag):
    g = tree.pybus_harvey_gamma()
    sys.stdout.write("%s\tgamma\t%f\n" % (tag, g))

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-s", "--seed", dest="seed", default=0,
        type="int",
        help="The random number generator seed")
    parser.add_option("-l", "--leaves", dest="leaves", default=4,
        type="int",
        help="The number of leaves to generate")
    parser.add_option("-c", "--sample", dest="sample", default=-1,
        type="int",
        help="The number of leaves to sample (-1 to sample all leaves generated)")
    parser.add_option("-b", "--birth", dest="birth", default=1.0,
        type="float",
        help="The birth (speciation rate). default 1.")
    parser.add_option("-d", "--death", dest="death", default=0.0,
        type="float",
        help="The death (extinction rate). default 0.")
    parser.add_option("-u", "--ultimate", dest="ultimate", default="random",
        type="string",
        help='The handling of the waiting time between the last speciation and the final time point: either "random" or "next_event". If "random" is selected another event waiting time will be simulated and the final time point will be a random U(0,1) variate multiplied by this waiting time.' )
    parser.add_option("-p", "--pop-size", dest="pop_size", default=1.0,
        type="float",
        help="The effective population size of each population")
    parser.add_option("-i", "--initial-mut-rate", dest="mu_rate", default=1.0e-8,
        type="float",
        help="The initial rate of sequence evolution for the root of the tree.")
    parser.add_option("-r", "--roeotroe", dest="roeotroe", default=1.0,
        type="float",
        help="The rate of evolution of the rate of evolution. This number will be multiplied by the gene tree's length (in terms of time) to determine the standard deviation of the log of the rate of evolution. The rate for the branch will then be the mean of the values at the endpoints.")
    parser.add_option("-f", "--min-mut-rate", dest="min_rate", default=1e-300,
        type="float",
        help="The minimum rate of sequence evolution")
    parser.add_option("-x", "--max-mut-rate", dest="max_rate", default=1e300,
        type="float",
        help="The maximum rate of sequence evolution")
    parser.add_option("--sp-tree-sim", dest="sp_tree_sim", default="bdtree",
        type="string",
        help="The name of the species tree simulator. One of the following: %s" % ", ".join(sp_tree_simualtors))
    parser.add_option("--gene-tree-sim", dest="gene_tree_sim", default="mccoal",
        type="string",
        help="The name of the gene tree simulator. One of the following: %s" % ", ".join(gene_tree_simualtors))
    parser.add_option("--n-reps", dest="n_reps", default=1,
        type="int",
        help="The number of replicate simulations to conduct")
    (options, args) = parser.parse_args()
    if args:
        exit("No arguments accepted.  Use the -h flag to see the list of command line options (flags) that are accepted")
    if options.seed < 1:
        import time
        options.seed = int(time.time()*1000)
    _LOG.debug("%s -s%d -l%d -c%d -b%g -d%g -u%s -p%f --sp-tree-sim=%s --gene-tree-sim=%s --n-reps%d -i%f\n" % (
                     sys.argv[0],
                     options.seed,
                     options.leaves,
                     options.sample,
                     options.birth,
                     options.death,
                     options.ultimate,
                     options.pop_size,
                     options.sp_tree_sim,
                     options.gene_tree_sim,
                     options.n_reps,
                     options.mu_rate
                     ))
    _RNG.seed(options.seed)
    simulate_sequences = False
    for rep_n in xrange(options.n_reps):
        sp_tree_time_output = open("sp-%d.tre" % rep_n, 'w')
        gene_tree_time_output = open("timegene-%d.tre" % rep_n, 'w')
        if simulate_sequences:
            sp_outgroup_output = open("sp-%d-outgroups.txt" % rep_n, 'w')
            gene_outgroup_output = open("gene-%d-outgroups.txt" % rep_n, 'w')
            rate_mod_tree_time_output = open("modelgene-%d.tre" % rep_n, 'w')
            inferred_tree_time_output = open("inferred-%d.tre" % rep_n, 'w')

        sys.stderr.write("rep %d\n" % rep_n)
        # Generate a Species tree
        sts = options.sp_tree_sim.lower()
        if sts == "bdtree":
            sp_tree, sp_outgroup = bdtree_sim_sp_tree(options)
            for e in sp_tree.preorder_edge_iter():
                if e.length is not None:
                    e.length = e.length * options.mu_rate
        else:
            exit("Expecting the species tree simulator to be one of the following: %s" % ", ".join(sp_tree_simualtors))


        sp_tree_time_output.write("%s;\n" % sp_tree)

        # Calculate the speciation rate stats for the true species tree
        divers_rate_stats(options, sp_tree, ("sp-%d" % rep_n))

        # Generate a Gene tree for the species tree
        gts = options.gene_tree_sim.lower()
        if gts == "mccoal":
            t = mccoal_sim_gene_tree(sp_tree, options)
            gene_tree, gene_out, culled, culled_out = t
        else:
            exit("Expecting the gene tree simulator to be one of the following: %s" % ", ".join(gene_tree_simualtors))

        gene_tree_time_output.write("%s;\n" % gene_tree)

        # Calculate the speciation rate stats for the true gene tree
        divers_rate_stats(options, gene_tree, ("fulltimegene-%d" % rep_n))
        # Calculate the speciation rate stats for the true gene tree
        divers_rate_stats(options, culled, ("culledtimegene-%d" % rep_n))

        if simulate_sequences:
            sp_outgroup_output.write("%s\n" % "\n".join(sp_outgroup))
            gene_outgroup_output.write("%s\n" % "\n".join(gene_out))

            # Make the gene tree non-ultrametric

            # Generate sequence data for the gene tree

            # Estimate an ultrametric tree for the simulated sequences

            # Calculate the speciation rate stats for the inferred branch lengths
            sp_outgroup_output.close()
            gene_outgroup_output.close()
            rate_mod_tree_time_output.close()
            inferred_tree_time_output.close()

        # Close files
        sp_tree_time_output.close()
        gene_tree_time_output.close()


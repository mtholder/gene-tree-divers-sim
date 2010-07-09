#!/usr/bin/env python
# Generates tree according to a birth-death process.
#   The amount of debugging output is determined by BDTREE_LOGGING_LEVEL which
#       can be set to "logging.DEBUG" or "logging.NOTSET" to increase the
#       verbosity

import sys, os, random, math
from cStringIO import StringIO
_RNG = random.Random()
FORMAT_LIST = ["newick", "nexus"]
class FORMAT_CODE:
    NEWICK = FORMAT_LIST.index("newick")
    NEXUS = FORMAT_LIST.index("nexus")
LAST_INTERVAL_LIST = ["random", "next_event"]
class LAST_INTERVAL:
    RANDOM = LAST_INTERVAL_LIST.index("random")
    NEXT_EVENT = LAST_INTERVAL_LIST.index("next_event")
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
            lfe =  os.environ.get("BDTREE_LOGGING_LEVEL")
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
_LOG = get_logger("bdtree")


class Node(object):
    def __init__(self, parent=None):
        self.name = str(id(self))
        self.parent = parent
        self.children = []
        self.duration = 0.0

    def count_tips(self):
       return max(sum([c.count_tips() for c in self.children]), 1)

    def write_tip_names(self, out, sep='\n'):
        c = self.children
        if c:
            for child in c:
                child.write_tip_names(out)
        else:
            out.write("%s%s" % (self.name, sep))

    def speciate(self):
        self.children.append(Node(parent=self))
        self.children.append(Node(parent=self))
        return self.children

    def prune_self(self):
        if self.parent:
            p = self.parent
            p.children.remove(self)
            if not p.children:
                p.prune_self()

    def add_time(self, t):
        self.duration += t

    def add_child(self, c):
        if c.parent:
            try:
                c.parent.children.remove(c)
            except: pass
        self.children.append(c)
        c.parent = self

    def suppress_deg_two(self):
        n = 0
        ch = list(self.children) # we need to make a shallow copy in case our kids modify our children list
        ch = self.children # we need to make a shallow copy in case our kids modify our children list
        for c in ch:
            n += c.suppress_deg_two()
        if len(self.children) == 1:
            c = self.children[0]
            sd = self.duration + c.duration
            self.duration = sd
            self.children = c.children
            for gc in self.children:
                gc.parent = self
            n += 1
            if c.name:
                self.name = c.name
        return n

    def write_newick(self, o, edge_lens=True):
        if self.children:
            f = True
            o.write("(")
            for c in self.children:
                if f:
                    f = False
                else:
                    o.write(",")
                c.write_newick(o, edge_lens)
            o.write(")")
        elif self.name:
            o.write(self.name)
        if edge_lens and self.parent:
            o.write(":%g" % (self.edge_len))

    def get_edge_len(self):
        return self.duration

    def __str__(self):
        return self.newick(edge_lens=False)

    def newick(self, edge_lens=True):
        o = StringIO()
        self.write_newick(o, edge_lens=edge_lens)
        return o.getvalue()

    edge_len = property(get_edge_len)

def remove_random(l):
    n = len(l)
    ind = _RNG.randrange(n)
    el = l[ind]
    l.pop(ind)
    return el


def sim_species_tree(leaves_to_gen, last_interval, options):
    birth = options.birth
    if birth <= 0.0:
        sys.exit("Speciation (birth) rate must be positive")
    death = max(0.0, options.death)
    event_rate = birth + death
    dprob = death/event_rate

    root = Node()
    evolving = list(root.speciate())
    _LOG.debug(root.newick(edge_lens=False))
    ntips = 2
    while ntips < leaves_to_gen:
        total_rate = ntips*event_rate
        waiting_time = _RNG.expovariate(total_rate)
        for i in evolving:
            i.add_time(waiting_time)
        to_mod = remove_random(evolving)
        if _RNG.random() < dprob:
            to_mod.prune_self()
            _LOG.debug("Deleting %s after %g" % (to_mod.name, waiting_time))
            _LOG.debug(root.newick(edge_lens=False))
            ntips -= 1
            if ntips == 0:
                return None, []
        else:
            _LOG.debug("Speciating %s after %g" % (to_mod.name, waiting_time))
            evolving.extend(to_mod.speciate())
            _LOG.debug(root.newick(edge_lens=False))
            ntips += 1
    # add a duration to all terminal branches.  Duration is U[0,1] * waiting time to the next event
    total_rate = ntips*event_rate
    curr_time = _RNG.expovariate(total_rate)
    if last_interval == LAST_INTERVAL.RANDOM:
        curr_time = _RNG.random() * curr_time
    for i in evolving:
        i.add_time(curr_time)
    return root, evolving

def nexusHeader(f, n):
    f.write("#NEXUS\nbegin taxa ;\n   dimensions ntax = %d ;\n   taxlabels" % n)
    for i in range(n):
        f.write(" t%d"% (i+1))
    f.write(";\nEND;\nbegin trees;\n")

def write_tree(f, t, fmt_code, tree_ind):
    if fmt_code == FORMAT_CODE.NEXUS:
        f.write("tree sim%d = [&R] " % (tree_ind + 1))
    f.write("%s;\n" % t)


if __name__ == '__main__':
    import re
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-s", "--seed", dest="seed", default=0,
        type="int",
        help="The random number generator seed")
    parser.add_option("-n", "--num", dest="n", default=1,
        type="int",
        help="The number of trees to generate")
    parser.add_option("-l", "--leaves", dest="ntax", default=4,
        type="int",
        help="The number of leaves to generate")
    parser.add_option("-c", "--sample", dest="nchosen", default=-1,
        type="int",
        help="The number of leaves to sample (-1 to sample all leaves generated)")
    parser.add_option("-b", "--birth", dest="birth", default=1.0,
        type="float",
        help="The birth (speciation rate). default 1.")
    parser.add_option("-d", "--death", dest="death", default=0.0,
        type="float",
        help="The death (extinction rate). default 0.")
    parser.add_option("-f", "--format", dest="format", default="newick",
        type="string",
        help='The format for the tree: either "newick" or "NEXUS"')
    parser.add_option("-u", "--ultimate", dest="ultimate", default="random",
        type="string",
        help='The handling of the waiting time between the last speciation and the final time point: either "random" or "next_event". If "random" is selected another event waiting time will be simulated and the final time point will be a random U(0,1) variate multiplied by this waiting time.' )
    (options, args) = parser.parse_args()
    if args:
        sys.exit("No arguments accepted.  Use the -h flag to see the list of command line options (flags) that are accepted")
    if options.seed < 1:
        import time
        options.seed = int(time.time()*1000)
    try:
        fmt_code = FORMAT_LIST.index(options.format.lower())
    except:
        sys.exit("Expecting the format to be one of the following:\n\t%s" % "\n\t".join(FORMAT_LIST))
    try:
        last_interval = LAST_INTERVAL_LIST.index(options.ultimate.lower())
    except:
        sys.exit("Expecting the format to be one of the following:\n\t%s" % "\n\t".join(LAST_INTERVAL_LIST))

    _LOG.debug("%s -s%d -n%d -l%d -c%d -b%g -d%g -f%s -u%s\n" % (
                     sys.argv[0],
                     options.seed,
                     options.n,
                     options.ntax,
                     options.nchosen,
                     options.birth,
                     options.death,
                     FORMAT_LIST[fmt_code],
                     LAST_INTERVAL_LIST[last_interval]
                     ))
    _RNG.seed(options.seed)
    n = options.n
    leaves_to_gen = options.ntax
    leaves_to_include = options.nchosen
    if leaves_to_include < 0:
        leaves_to_include = leaves_to_gen
    if leaves_to_include > leaves_to_gen:
        sys.exit("--sample value cannot be larger than --leaves value.")

    if fmt_code == FORMAT_CODE.NEXUS:
        nexusHeader(sys.stdout, leaves_to_include)

    #   Also write a PAUP outgroup command to bd_outgroups_paup.nex
    outgroupsfile = open(".bd_outgroups.txt", "w")
    #outgroupsfile.write("#NEXUS\n")

    no_brlen_pattern = re.compile(r't\d+(,|\))')
    for i in xrange(n):

        bdseed = i + options.seed
        _RNG.seed(bdseed)
        _LOG.warn("debugging seed %d is being used!"% bdseed)

        sp_tree = None
        if leaves_to_include < 2:
            if leaves_to_include == 0:
                sp_tree = "();"
            else:
                sp_tree = "(t1);"
        else:
            while sp_tree is None:
                sp_tree, tips = sim_species_tree(leaves_to_gen, last_interval, options)
            _LOG.debug("preprune :%s\n" % str(sp_tree))
            rleaves = _RNG.sample(tips, leaves_to_include)
            if leaves_to_include < leaves_to_gen:
                for l in tips:
                    if not (l in rleaves):
                        _LOG.debug("Deleting %s" % l.name)
                        l.prune_self()
                        _LOG.debug(sp_tree.newick(edge_lens=False))
            for n, l in enumerate(rleaves):
                l.name = "t%d" % (n+1)
            sp_tree.suppress_deg_two()
        newick_str = sp_tree.newick(edge_lens=True)

        if no_brlen_pattern.search(newick_str):
            sys.exit("No branch length")
        from dendropy import dataio
        if False:
            reread_sp_tree = dataio.trees_from_string(string=newick_str, format="NEWICK")[0]
            reread_newick = str(reread_sp_tree)
            _LOG.warn(reread_newick)
            if no_brlen_pattern.search(reread_newick):
                sys.exit("No branch length in reread string:\n%s" % reread_newick)

        write_tree(sys.stdout, newick_str, fmt_code, i)

        c = sp_tree.children[0]
        if sp_tree.children[0].count_tips() > sp_tree.children[1].count_tips():
            c = sp_tree.children[1]
        c.write_tip_names(outgroupsfile)
    outgroupsfile.close()

    if fmt_code == FORMAT_CODE.NEXUS:
        sys.stdout.write("END;\n")

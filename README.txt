MCcoal must be on your PATH.

Download http://abacus.gene.ucl.ac.uk/software/MCMCcoal.html
unpack it and build it.
Put the MCcoal subdirectory on your PATH. In bash this is something like:
     export PATH=$HOME/builds/phylo/MCMCcoal/MCcoal:$PATH
     
invoc.sh is simply a call to master.py

Use:
    master.py -h
to see all of the help options.

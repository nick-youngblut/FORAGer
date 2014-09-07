FORAGer - Find Orthologous Reads And Genes
==========================================

## WARNING: this project is still under development its use is not recommended at this time.

FORAGer is a set of scripts that can be used to
try and identify unassembled genes in draft assemblies
by utilizing gene clusters from related genomes. 

The process involves mapping reads to the gene clusters,
assembling the reads that map, and determining the likelihood
that the resulting contigs contain unassembled genes orthologous
to the focal gene cluster.


# INSTALLATION

Add the bin folder to your $PATH if you would like.


# SUPPORT AND DOCUMENTATION

All scripts provide good documentation.
To access the documentation, use:

    perldoc script.pl

Or you can use:

    script.pl -h


LICENSE AND COPYRIGHT

Copyright (C) 2013 Nick Youngblut

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See L<http://dev.perl.org/licenses/> for more information.


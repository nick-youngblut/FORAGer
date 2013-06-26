#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Parallel::ForkManager;
use File::Path qw/remove_tree/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $runID, $fig, $PA_in, $screen_in);
my $overlap = 0.05;
my $evalue = "1e-30";
my $length = 0.75;
my $outdir = "passed_tblastn";
GetOptions(
		"runID=s" => \$runID,
		"fig=s" => \$fig,
		"overlap=f" => \$overlap,
		"evalue=s" => \$evalue,
		"name=s" => \$outdir,
		"PA=s" => \$PA_in,
		"screen=s" => \$screen_in,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a directory!\n" unless $ARGV[0];
die " ERROR: provide a cluster runID!\n" unless $runID;
die " ERROR: provide an organism ID!\n" unless $fig;
$ARGV[0] = File::Spec->rel2abs($ARGV[0]);
$outdir = File::Spec->rel2abs($outdir);


### MAIN
# getting passed contig file names #
my $files_r = get_file_names($ARGV[0]);

# tblastn #
make_tblastn_dir($outdir);
my %res;
foreach my $infile (keys %$files_r){
	call_tblastn_wrapper($infile, $outdir, $files_r->{$infile}, $runID, $fig, \%res);
	}

# writing whether gene exists where cluster blasted #
foreach my $cluster (keys %res){
	print join("\t", $cluster, $res{$cluster}), "\n";
	}

### Subroutines
sub make_tblastn_dir{
	my ($outdir) = @_;
	remove_tree($outdir) if -d $outdir;
	mkdir $outdir or die $!;
	}

sub call_tblastn_wrapper{
	my ($infile, $outdir, $cluster, $runID, $fig, $res_r) = @_;
	
	# opening pipe for db_TBlastN_wrapper.py #
	my $cmd = "printf \"$runID\\t$cluster\\n\" | db_getClusterGeneInformation.py | db_TBlastN_wrapper.py -o $fig |";
		#print STDERR "$cmd\n" unless $verbose;
	open PIPE, $cmd or die $!;
	
	# making tblast output file #
	(my $tblastn_out = $infile) =~ s/\.[^.]+$|$/_blastn.txt/;
	open OUT, ">$outdir/$tblastn_out" or die $!;
	
	# reading from PIPE #
	while(<PIPE>){
		print OUT;
		next if $res_r->{$cluster} && $res_r->{$cluster} eq "exists"; 	# skip if already determined to exist
		
		chomp;
		my @line = split /\t/;
		
		if($line[11] eq "SAMESTRAND" && $line[15] >= $overlap && $line[8] <= $evalue){
			$res_r->{$cluster} = "exists";
			}
		else{ $res_r->{$cluster} = "new gene" };
		}
	close PIPE;
	close OUT;
	}

sub get_file_names{
# getting cluster files names from provided directory #
	my ($dir) = @_;
	opendir IN, $dir or die $!;
	my @files = grep /^[^.]/, readdir IN;
	closedir IN;
	
	my %files;
	foreach my $file (@files){
		(my $tmp = $file) =~ s/clust(\d+).+/$1/;
		$files{$file} = $tmp;
		}
	
	return \%files;
	}


__END__

=pod

=head1 NAME

FORAGer_tblastn_filter.pl -- use tblastn to filter out contigs of existing genes

=head1 SYNOPSIS

FORAGer_tblastn_filter.pl [flags] directory

=head2 required flags

=over

=item -runID

Cluster runID.

=item -fig

ITEP organim FIG ID.

=back

=head2 optional flags

=over

=item -overlap

Fraction overlap with existing gene to call the gene pre-existing. [0.05]

=item -evalue

tblastn evalue cutoff. [1e-30]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc FORAGer_tblastn_filter.pl

=head1 DESCRIPTION

Use tblastn to determine which of the 'passed'
contigs overlap with pre-existing genes in the 
query genome. This will help screen out
false postives.

=head2 Requires:

db_TBlastN_wrapper.py

=head1 EXAMPLES

=head2 Usage: 

FORAGer_tblastn_filter.pl Mapped2Cluster_query_passed/ -r mazei_I_2.0_c_0.4_m_maxbit -org 2209.17

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/FORAGer/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


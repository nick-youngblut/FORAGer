#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Parallel::ForkManager;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $runID, $org);
my $overlap = 0.05;
my $evalue = "1e-30";
GetOptions(
		"runID=s" => \$runID,
		"organism=s" => \$org,
		"overlap=f" => \$overlap,
		"evalue=s" => \$evalue,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a directory!\n" unless $ARGV[0];
$ARGV[0] = File::Spec->rel2abs($ARGV[0]);
die " ERROR: provide a cluster runID!\n" unless $runID;
die " ERROR: provide an organism ID!\n" unless $org;

### MAIN
my $files_r = get_file_names($ARGV[0]);

# tblastn #
my %res;
foreach my $file (keys %$files_r){
	call_tblastn_wrapper($file, $files_r->{$file}, $runID, $org, \%res);
	}

foreach my $cluster (keys %res){
	print join("\t", $cluster, $res{$cluster}), "\n";
	}

### Subroutines
sub call_tblastn_wrapper{
	my ($file, $cluster, $runID, $org, $res_r) = @_;
	my $cmd = "printf \"$runID\\t$cluster\\n\" | db_getClusterGeneInformation.py | db_TBlastN_wrapper.py -o $org |";
	print STDERR "$cmd\n" unless $verbose;
	open PIPE, $cmd or die $!;
	
	while(<PIPE>){
		next if $res_r->{$cluster} && $res_r->{$cluster} eq "exists"; 
		
		chomp;
		my @line = split /\t/;
		
		if($line[11] eq "SAMESTRAND" && $line[15] >= $overlap && $line[8] <= $evalue){
			$res_r->{$cluster} = "exists";
			}
		else{ $res_r->{$cluster} = "new gene" };
		}
	close PIPE;	
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

=item -organism

ITEP organim ID (FIG).

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


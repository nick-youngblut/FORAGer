#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path qw/remove_tree/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $runID, $fig, $PA_in, $screen_in);
my $overlap = 0.05;
my $evalue = "1e-30";
my $length = 0.80;
my $outdir = "passed_tblastn";
GetOptions(
		"runID=s" => \$runID,
		"fig=s" => \$fig,
		"overlap=f" => \$overlap,
		"evalue=s" => \$evalue,
		"length=i" => \$length,
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
print STDERR "\n...tblastn output files written to $outdir\n\n" unless $verbose;

# writing whether gene exists where cluster blasted #
foreach my $cluster (keys %res){
	print join("\t", $cluster, $res{$cluster}), "\n";
	}

# updating PA and FORAGer_screen summary #
if($PA_in){
	$PA_in = File::Spec->rel2abs($PA_in);
	die " ERROR: $PA_in not found!\n" unless $PA_in;
	update_PA($PA_in, \%res);
	}
if($screen_in){
	$screen_in = File::Spec->rel2abs($screen_in);
	die " ERROR: $screen_in not found!\n" unless $screen_in;
	update_screen($screen_in, \%res);
	}


### Subroutines
sub update_screen{
# updating FORAGer_screen summary file; including columns tblastn #
	my ($screen_in, $res_r) = @_;
	open IN, $screen_in or die $!;
	(my $screen_out = $screen_in) =~ s/\.[^.]+$|$/_filter.txt/;
	open OUT, ">$screen_out" or die $!;
	
	while(<IN>){
		chomp;
		my @line = split /\t/;
		(my $clust = $line[1]) =~ s/clust|\.[^.]+$//g;
		die " ERROR: connot determine cluster number from $line[1]! Formatting incorrect!\n" 
			unless $clust =~ /^\d+$/;
			
		if (exists $res_r->{$clust} ){
			if($res_r->{$clust} eq "overlapping_gene"){ 		# if tblastn says exist, modify to say did not pass
				$line[2] = 0;		# overall failed 
				print OUT join("\t", @line, "tblast_filter_fail"), "\n";
				}
			else{
				print OUT join("\t", @line, "tblast_filter_pass"), "\n";
				}
			}
		else{
			print OUT join("\t", @line, "NA"), "\n";
			}		
		}
	close IN;
	close OUT;
	
	print STDERR "...tblastn-filtered FORAGer_screen summary file written: $screen_out\n\n" unless $verbose;	
	}

sub update_PA{
# updating PA; just writing contigs that appear to be read based on tblastn #
	my ($PA_in, $res_r) = @_;

	open IN, $PA_in or die $!;
	(my $PA_out = $PA_in) =~ s/\.[^.]+$|$/_filter.txt/;
	open OUT, ">$PA_out" or die $!;
	
	while(<IN>){
		my @line = split /\t/;
		(my $contig = $line[0]) =~ s/FORAGer__|__\d+$//g;

		print OUT unless exists $res_r->{$line[7]} &&
				$res_r->{$line[7]} eq "overlapping_gene";		# if tblastn says exist, don't include as new gene
		}
	close IN;
	close OUT;
	
	print STDERR "...tblastn-filtered Pres-Abs file written: $PA_out\n\n" unless $verbose;
	}

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
		last if exists $res_r->{$cluster} && $res_r->{$cluster} eq "overlapping_gene";  	# skip if already determined to exist
		
		chomp;
		my @line = split /\t/;
		
		# good hit? #
		if($line[8] <= $evalue &&				# tblastn hit <= e-value cutoff
			$line[6]/$line[1] >= $length){		# tblastn hit >= X% query length

			# overlapping with alread called gene? #
			if($line[11] eq "SAMESTRAND" && 			# gene already called on same strand
				$line[15] >= $overlap){					# gene overlaps w/ tblastn hit
				$res_r->{$cluster} = "overlapping_gene";		
				}	
			else{
				$res_r->{$cluster} = "new_gene";		
				}
			}
		else{ $res_r->{$cluster} = "new_gene" };
		}
	close PIPE;
	close OUT;
		#print Dumper %$res_r; exit;
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

FORAGer_tblastn_filter.pl -- use tblastn to filter out false 'passed' FORAGer contigs

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

=item -PA

Presence-absence file from FORAGer_screen.pl (will be updated with tblastn results).

=item -screen

FORAGer_screen.pl summary file (will be updated with tblastn results).

=item -overlap

Fraction overlap with existing gene to call the gene pre-existing. [0.05]

=item -evalue

tblastn evalue cutoff. [1e-30]

=item -length

minimum tblast hit length (fraction of query length). [0.80]

=item -name

Name of tblastn output file directoy. ['passed_tblastn']

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc FORAGer_tblastn_filter.pl

=head1 DESCRIPTION

Use tblastn to determine which of the 'passed'
contigs for FORAGer_screen.pl overlap with 
pre-existing genes in the query genome. 
All of the genes in the focal cluster will be
used for the tblastn query. 

=head3 Output to STDOUT: 

=over

=item 1) cluster_id

=item 2) 'overlapping_gene' OR 'new_gene'

=back

Any of the queries have good hits to a genome region
overlapping a pre-existing gene, the query contig
from FORAGer is considered pre-existing. 

If the presence-absence & screening summary
files from FORAGer_screen.pl are provided, both
tables are updated with the tblastn information:

=over

=item * 

Only tblastn-passed contigs will be written to 
the PA table.

=item * 

The screen summary table will have a final tblastn
summary & the binary pass/fail column will reflect
the tblastn results.

=back

=head2 Requires:

db_TBlastN_wrapper.py

=head1 EXAMPLES

=head2 Basic Usage:

FORAGer_tblastn_filter.pl -r mazei_I_2.0_c_0.4_m_maxbit -fig 2209.17 Mapped2Cluster_query_passed/ 

=head2 Updating PA and screen files

FORAGer_tblastn_filter.pl -r mazei_I_2.0_c_0.4_m_maxbit -fig 2209.17 -PA FORAGer_PA.txt -screen FORAGer_screen.txt Mapped2Cluster_query_passed/ 

=head2 Altering tblastn 'overlapping_gene' cutoffs

FORAGer_tblastn_filter.pl -r mazei_I_2.0_c_0.4_m_maxbit -fig 2209.17 -o 0.2 -l 0.9 Mapped2Cluster_query_passed/ 

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/FORAGer/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


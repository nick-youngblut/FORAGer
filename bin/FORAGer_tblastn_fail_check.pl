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

my ($verbose, $runID, $fig, $PA_in, $screen_bool);
my $overlap = 0.80;
my $evalue = "1e-30";
my $length = 0.8;
my $outdir = "fail_check_tblastn";
my $hit_frac = 0.5;
GetOptions(
		"runID=s" => \$runID,
		"fig=s" => \$fig,					# FIG query
		"overlap=f" => \$overlap,			# overlap
		"evalue=s" => \$evalue,				# evalue cutoff for 'good' hit
		"length=f" => \$length,				# length cutoff for 'good' hit
		"name=s" => \$outdir,				# output ditrection name
		"screen" => \$screen_bool,			# update screen? [TRUE]
		"x=f" => \$hit_frac,			# fraction of PEGs in cluster that must have an good overlapping hit
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a FORAGer_screen.pl summary file (*screen.txt)!\n" unless $ARGV[0];
die " ERROR: provide a cluster runID!\n" unless $runID;
die " ERROR: provide an organism ID!\n" unless $fig;
$ARGV[0] = File::Spec->rel2abs($ARGV[0]);
$outdir = File::Spec->rel2abs($outdir);


### MAIN
# getting passed contig file names #
my $clusters_r = get_failed_clusters($ARGV[0]);

# tblastn #
make_tblastn_dir($outdir);

my %res;
foreach my $cluster (keys %$clusters_r){
	next if $clusters_r->{$cluster} eq "PASSED";
	call_tblastn_wrapper($cluster, $outdir, $runID, $fig, \%res);
	}
print STDERR "\n...tblastn output files written to $outdir\n\n" unless $verbose;

# writing whether gene exists where cluster blasted #
foreach my $cluster (sort {$a<=>$b} keys %res){
	print join("\t", $cluster, $res{$cluster}), "\n";
	}

# making tblastn-passed table #
make_tblastn_passed_PA(\%res);

# updating screen summary table #
update_screen($ARGV[0], \%res) unless $screen_bool;



### Subroutines
sub update_screen{
# updating FORAGer_screen summary file; including columns tblastn #
	my ($screen_in, $res_r) = @_;
	
	open IN, $screen_in or die $!;
	my $screen_out = "FORAGer_screen_TblastnCheck.txt";
	open OUT, ">$screen_out" or die $!;
	
	while(<IN>){
		chomp;
		my @line = split /\t/;
		(my $clust = $line[1]) =~ s/clust|\.[^.]+$//g;
		die " ERROR: connot determine cluster number from $line[1]! Formatting incorrect!\n" 
			unless $clust =~ /^\d+$/;
			
		if ( exists $res_r->{$clust} ){
			if($res_r->{$clust} eq "gene_exists"){ 		# if tblastn says exist, modify to say did not pass
				$line[2] = 1;		# overall failed 
				print OUT join("\t", @line, "tblastn_check_pass"), "\n";
				}
			else{
				print OUT join("\t", @line, "tblast_check_fail"), "\n";
				}
			}
		else{
			print OUT join("\t", @line, "NA"), "\n";
			}		
		}
	close IN;
	close OUT;
	
	print STDERR "...tblastn-checked FORAGer_screen summary file written: $screen_out\n\n" unless $verbose;	
	}

sub make_tblastn_passed_PA{
# updating PA; just writing contigs that appear to be read based on tblastn #
	my ($res_r) = @_;

	my $PA_out = "FORAGer_PA_TblastnCheck.txt";
	open OUT, ">$PA_out" or die $!;
	
	foreach my $cluster (keys %$res_r){
		next if $cluster eq "gene_not_found";
		my $user_geneid = join("__", "FORAGer", "cluster$cluster", time());
		print OUT join("\t", $user_geneid, $fig, "CDS",
					"", "", "",
					$runID, $cluster, "",
					"", ""), "\n";		
		}
	close OUT;
	
	print STDERR "...tblastn-checked Pres-Abs file written: $PA_out\n\n" unless $verbose;
	}

sub make_tblastn_dir{
	my ($outdir) = @_;
	remove_tree($outdir) if -d $outdir;
	mkdir $outdir or die $!;
	}

sub call_tblastn_wrapper{
	my ($cluster, $outdir, $runID, $fig, $res_r) = @_;
	
	# getting number of pegs in cluster #
	my $cmd = "printf \"$runID\\t$cluster\\n\" | db_getGenesInClusters.py | wc -l |";
	open PIPE, $cmd or die $!;
	chomp(my $n_pegs = <PIPE>);
	close PIPE;
	
	# opening pipe for db_TBlastN_wrapper.py #
	$cmd = "printf \"$runID\\t$cluster\\n\" | db_getClusterGeneInformation.py | db_TBlastN_wrapper.py -o $fig |";
	open PIPE, $cmd or die $!;
	
	# making tblast output file #
	(my $tblastn_out = $cluster) =~ s/\.[^.]+$|$/_blastn.txt/;
	open OUT, ">$outdir/$tblastn_out" or die $!;
	
	# reading from PIPE #
	my %overlapping;
	while(<PIPE>){
		print OUT;
		
		chomp;
		my @line = split /\t/;
		
		# good hit? #
		if($line[8] <= $evalue &&				# tblastn hit <= e-value cutoff
			$line[6]/$line[1] >= $length){		# tblastn hit >= X% query length

			# overlapping with alread called gene? #
			if($line[11] eq "SAMESTRAND" && 			# gene already called on same strand
				$line[15] >= $overlap){					# gene overlaps w/ tblastn hit	
				$overlapping{$line[0]} = 1;				# peg has good overlapping hit
				}	
			}
		}
	close PIPE;
	close OUT;
	
	# determine if number of overlapping hits meets cutoff #
	my $overlap_frac = (scalar keys %overlapping) / $n_pegs;
	print STDERR "\nCluster: $cluster; Fraction_PEGs_with_good_overlapping_tblastn_hits: $overlap_frac\n\n"
		unless $verbose;
	if( $overlap_frac >= $hit_frac){		# if cutoff met
		$res_r->{$cluster} = "gene_exists";
		}
	else{ $res_r->{$cluster} = "gene_not_found"; }
		
		#print Dumper %$res_r; exit;
	}

sub get_failed_clusters{
# getting clusters that failed according to *screen.txt file #
	my ($screen_in) = @_;
	
	my %clusters;
	open IN, $screen_in or die $!;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;
		#next unless $line[3] == 0;		# just failed clusters
		
		(my $cluster = $line[1]) =~ s/clust|\.[^.]+$//g;
		die " ERROR: cannot determine cluster number from $line[1]\n" 
			unless $cluster =~ /^\d+$/;
		
		if(exists $clusters{$cluster}){
			if($line[3] eq "NA"){
				$clusters{$cluster} = "FAILED";
				}
			elsif($line[3] == 1 && $clusters{$cluster} eq "FAILED"){
				$clusters{$cluster} = "PASSED";
				}
			}
		else{
			if($line[3] eq "NA"){
				$clusters{$cluster} = "FAILED";
				}
			elsif($line[3] == 1){
				$clusters{$cluster} = "PASSED";
				}
			else{
				$clusters{$cluster} = "FAILED";
				}
			}
		}

		#print Dumper %clusters; exit;
	return \%clusters;
	}


__END__

=pod

=head1 NAME

FORAGer_tblastn_filter.pl -- use tblastn to filter out false 'passed' FORAGer contigs

=head1 SYNOPSIS

FORAGer_tblastn_filter.pl [flags] FORAGer_screen.txt

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

Fraction overlap with existing gene to call the gene pre-existing. [0.8]

=item -evalue

tblastn evalue cutoff. [1e-30]

=item -length

minimum tblast hit length (fraction of query length). [0.8]

=item -x

Fraction of query PEGs that must have good overlapping hits (>=). [0.5]

=item -name

Name of tblastn output file directoy. ['fail_check_tblastn']

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc FORAGer_tblastn_filter.pl

=head1 DESCRIPTION

Some 'failed' contigs from FORAGer may
be caused by multiple copies of the gene
in the genome, which breaks the assembly.

This script uses tblastn to check to see
if the gene is actually present in the genome.

=head3 Gene tblastn hit 'presence' criteria:

=over

=item * 	The e-value cutoff must be met.

=item * 	The hit length cutoff must be met.

=item * 	The hit must be on the same strand as a pre-existing gene.

=item * 	The hit must overlap to X% as the pre-existing gene.

=back 

=head3 Output to STDOUT: 

=over

=item 1) cluster_id

=item 2) 'gene_exists' OR 'gene_not_found'

=back

If '-x' fraction of query PEGs have good overlapping
hits to pre-existing genes, the query contig
from FORAGer is considered pre-existing or an 
artefact (i.e. not a real new gene). 

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

ITEP & db_TBlastN_wrapper.py

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


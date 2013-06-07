#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use List::Util qw/min/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my $verbose;
my $len_cutoff = 0;
my $count_cutoff = 1;
my $multi_cutoff = 1;
my $annotation_filter;
GetOptions(
	   "length=i" => \$len_cutoff,
	   "count=i" => \$count_cutoff,
	   "multi=i" => \$multi_cutoff,
	   "annotation=s" => \$annotation_filter,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
my @annotation_filter = split / *, */, $annotation_filter if $annotation_filter;
@annotation_filter = ("conserved protein", "hypothetical")
	if $annotation_filter && $annotation_filter eq "DEFAULT";

### MAIN
my $tbl_r = load_gene_info_tbl();
filter_by_length($tbl_r, $len_cutoff) if $len_cutoff;
filter_by_count($tbl_r, $count_cutoff) if $count_cutoff;
filter_by_multi($tbl_r, $multi_cutoff) if $multi_cutoff;
filter_by_annotation($tbl_r, \@annotation_filter) if @annotation_filter;
write_table($tbl_r);

### Subroutines
sub write_table{
	my ($tbl_r) = @_;
	foreach my $cluster (keys %$tbl_r){
		foreach my $row ( @{$tbl_r->{$cluster}{"row"}} ){ 
			print join("\t", @$row), "\n";
			}
		}
	}

sub filter_by_annotation{
	my ($tbl_r, $annot_r) = @_;
	foreach my $cluster (keys %$tbl_r){ 
		my $bad_annot_cnt = 0;
		# checking for annotations #
		foreach my $annot ( @{$tbl_r->{$cluster}{"annotation"}} ){
			foreach my $filt (@$annot_r){
				if( $annot =~ $filt ){
					$bad_annot_cnt++;
					last;
					}
				}
			}
		# removing contig if all pegs have annotation #
		 if( $bad_annot_cnt == scalar(@{$tbl_r->{$cluster}{"row"}}) ){
		 	delete $tbl_r->{$cluster};
		 	}
		}	
	}

sub filter_by_multi{
# filtering by multiple copies in taxa #
	my ($tbl_r, $multi_cutoff) = @_;
	foreach my $cluster (keys %$tbl_r){
		my $max_multi = 0;
		foreach my $taxon ( keys %{$tbl_r->{$cluster}{"copy"}} ){
			$max_multi = $tbl_r->{$cluster}{"copy"}{$taxon} if
				$tbl_r->{$cluster}{"copy"}{$taxon} > $max_multi;
			}
		delete $tbl_r->{$cluster} if $max_multi < $multi_cutoff;
		}
		#print Dumper $tbl_r; exit;
	}

sub filter_by_count{
# filtering by number in cluster #
	my ($tbl_r, $count_cutoff) = @_;
	foreach my $cluster (keys %$tbl_r){ 
		delete $tbl_r->{$cluster} if scalar(@{$tbl_r->{$cluster}{"row"}}) < $count_cutoff;
		}
	}

sub filter_by_length{
# filtering by gene length #
	my ($tbl_r, $len_cutoff) = @_;
	foreach my $cluster (keys %$tbl_r){
		delete $tbl_r->{$cluster} if min(@{$tbl_r->{$cluster}{"length"}}) < $len_cutoff;
		}
	}

sub load_gene_info_tbl{
# loading gene info table #
	my %ginfo;			# cluster => rows|lengths|annotation
	while(<>){
		chomp;
		next if /^\s*$/;
	
		my @line = split /\t/;		
		die " ERROR: the last column should have the cluster ID!\n" 
			unless $line[$#line] =~ /^\d+$/;
		
		# loading table #
		push(@{$ginfo{$line[$#line]}{"row"}}, \@line);
		push(@{$ginfo{$line[$#line]}{"length"}}, abs($line[6] - $line[5]));
		push(@{$ginfo{$line[$#line]}{"annotation"}}, $line[9]);
		$ginfo{$line[$#line]}{"copy"}{$line[1]}++;							# copy number
		}
		#print Dumper %ginfo; exit;
	return \%ginfo;
	}

__END__

=pod

=head1 NAME

db_filterGeneInfoTable.pl -- filtering clusters in ITEP gene info table

=head1 SYNOPSIS

db_filterGeneInfoTable.pl [options] < geneinfo.txt > geneinfo_edit.txt

=head2 options

=over

=item -length

Minimum length cutoff (bp) of all pegs in a cluster. [0]

=item -count

Minimum number of pegs in a cluster. [1]

-item -multi

Minimum number of gene copies in >= 1 genome (for selecting multi-copy gene clusters). [1]

=item -annotation

Filter out any clusters where all pegs have provided annotations. 
Comma delimited list. ['DEFAULT' = "conserved protein, hypothetical"]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc db_filterGeneInfoTable.pl

=head1 DESCRIPTION

Convience script for filtering out clusters in a gene info table
produced by ITEP.

=head1 EXAMPLES

=head2 Usage: filtering by length (>= 300bp)

db_filterGeneInfoTable.pl -l 300 < gene_info.txt > filtered_gene_info.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/FORAGer/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


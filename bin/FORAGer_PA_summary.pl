#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $runID, $taxa_in);
GetOptions(
	   "runID=s" => \$runID,
	   "taxa=s" => \$taxa_in,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a runID!\n" unless $runID;

### MAIN
# loading taxa of interest #
my $taxa_r = load_taxa_list($taxa_in) if $taxa_in;

# parsing PA table #
my ($PA_r, $toi_r) = parse_PA($runID, $taxa_r);

# core/variable; multi-copy #
core_var_copy_PA($PA_r, $toi_r);

# overlapping presence #



### subroutines
# quantify:
## N-taxa w/ ITEP gene
## N-taxa w/ FORAGer gene
## max number ITEP genes 
## N overlap in presence
## N conflict in presence (FORAGer: A, ITEP: P)
## N conflict in presence (FORAGer: P, ITEP: A)


sub core_var_copy_PA{
# determining which clusters are core, variable, and multi-copy #

	my ($PA_r, $toi_r) = @_;
	
		print Dumper $PA_r; exit;
	
	my %core_var;
	foreach my $cluster (keys %$PA_r){
		# core/variable #
		my $ntaxa_ITEP = scalar keys %{$PA_r->{$cluster}{"ITEP"}};
		my $max_N = 0;
		
		# number of genes per taxon #
		foreach my $taxon (keys %{$PA_r->{$cluster}{"ITEP"}} ){
			# max number of ITEP genes 
			my $N_genes = scalar @{$PA_r->{$cluster}{"ITEP"}{$taxon}};
			$max_N = $N_genes if $N_genes > $max_N;
			}
		
		# adding to hash #
		## core ##
		if($ntaxa_ITEP == scalar @$toi_r){
			$core_var{$cluster}{"core_var"} = "core";
			}
		else{
			$core_var{$cluster}{"core_var"} = "variable";
			}
		## multi-copy ##
		if($max_N == 1){
			$core_var{$cluster}{"copy"} = "single";
			}
		else{
			$core_var{$cluster}{"copy"} = "multi";
			}
		}
	
		print Dumper %core_var; exit;
	}

sub parse_PA{
# parsing PA table from ITEP #
	my ($runID, $taxa_r) = @_;

	my %PA;			# quantify as above
	my @header;
	while(<>){
		chomp;
		next if /^\s*$/;
		
		# header #
		if($.==1){
			@header = split /\t/;
			}
		# body #
		else{
			my @line = split /\t/;
			next unless $line[0] eq $runID;
			
			for my $i (3..$#line){			# each organism
				# sanity check #
				die " ERROR: line $. is not the same length as the header!\n"
					unless $header[$i];
				
				# skipping if don't care #
				if(%$taxa_r){
					next unless exists $taxa_r->{$header[$i]};
					}
					
				# laoding hash #
				my @entries = split /;/, $line[$i];
				
					#print Dumper $header[$i], @entries; 
				
				foreach my $j (@entries){
					if($j =~ /^NODE/){
						push(@{$PA{$line[1]}{"FORAGer"}{$header[$i]}}, $j);
						}
					else{
						push(@{$PA{$line[1]}{"ITEP"}{$header[$i]}}, $j) unless $j eq "NONE";
						}
					}				
				}
			}
		}
	
	# N-taxa of interest #
	my @toi;
	if(%$taxa_r){
		@toi = keys %$taxa_r;
		}
	else{ @toi = @header[3..$#header]; }
	
		#print Dumper %PA; exit;
	return \%PA, \@toi;
	}
	
	
sub load_taxa_list{
# loading list of taxa of interest #
	my ($taxa_in) = @_;
	open IN, $taxa_in or die $!;

	my %taxa;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		$taxa{$_} = 1;
		}
	close IN;

		#print Dumper @taxa; exit;
	return \%taxa;
	}

__END__

=pod

=head1 NAME

template.pl -- script template

=head1 SYNOPSIS

template.pl [options] < input > output

=head2 options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc template.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2
   3) Step 3

=head1 EXAMPLES

=head2 Usage method 1

template.pl <read1.fastq> <read2.fastq> <output directory or basename>

=head2 Usage method 2

template.pl <library file> <output directory or basename>

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


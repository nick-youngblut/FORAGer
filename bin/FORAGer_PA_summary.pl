#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use List::Util qw/sum/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $taxa_in);
GetOptions(
	   "taxa=s" => \$taxa_in,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
# loading taxa of interest #
my $taxa_r = load_taxa_list($taxa_in) if $taxa_in;

# parsing PA tables #
my ($PA_ITEP_r, $toi_ITEP_r) = parse_PA($ARGV[0], $taxa_r);
my ($PA_user_r, $toi_user_r) = parse_PA($ARGV[1], $taxa_r);

# core/variable; multi-copy #
## determining core/variable & copy ##
my $core_var_ITEP_r = core_var_copy_PA($PA_ITEP_r, $toi_ITEP_r);
my $core_var_user_r = core_var_copy_PA($PA_user_r, $toi_user_r);
## writing summary ##
write_core_var_summary($core_var_ITEP_r, "ITEP");
write_core_var_summary($core_var_user_r, "user");
## determining conflicts in core/variable ##
core_var_conflicts($core_var_ITEP_r, $core_var_user_r);


### subroutines
sub core_var_conflicts{
# determining conflicts in core/var & copy #
	my ($core_var_ITEP_r, $core_var_user_r) = @_;
	
	my %conflicts;
	foreach my $cluster (keys %$core_var_ITEP_r){
		if(exists $core_var_user_r->{$cluster}){		# intersection
			if($core_var_ITEP_r->{$cluster}{'core_var'} ne
				$core_var_user_r->{$cluster}{'core_var'} ||
				$core_var_ITEP_r->{$cluster}{'copy'} ne
				$core_var_user_r->{$cluster}{'copy'} ){
				
				print join("\t", "Conflict:", $cluster, 
					"Present", "Present", 
					$core_var_ITEP_r->{$cluster}{'core_var'}, 
					$core_var_ITEP_r->{$cluster}{'copy'},
					$core_var_user_r->{$cluster}{'core_var'}, 
					$core_var_user_r->{$cluster}{'copy'} ), "\n";	
				}
			}
		else{									# ITEP-specific
			print join("\t", "Conflict:", $cluster, "Present", "Absent",
					$core_var_ITEP_r->{$cluster}{'core_var'}, 
					$core_var_ITEP_r->{$cluster}{'copy'}, "NA", "NA"), "\n";
			}
		}
	foreach my $cluster (keys %$core_var_user_r){			# user-specific
		print join("\t", "Conflict:", $cluster, "Absent", "Present", "NA", "NA",
						$core_var_user_r->{$cluster}{'core_var'}, 
						$core_var_user_r->{$cluster}{'copy'} ), "\n";
		}
	}


sub write_core_var_summary{
# counting number of core/variable & single/multi
	my ($core_var_r, $source) = @_;
	
	my %summary;
	foreach my $cluster (keys %$core_var_r){
		die " LOGIC ERROR: $!\n" unless exists $core_var_r->{$cluster};
		$summary{ join("_", $core_var_r->{$cluster}{'core_var'}, 
						   $core_var_r->{$cluster}{'copy'}) }++;
		} 
	
	# Summary #
	print STDERR "\n" if $source eq "user";
	print STDERR "### $source ###\n";
	foreach my $cat (keys %summary){
		print STDERR join("\t", "$cat:", $summary{$cat}), "\n";
		}
	print STDERR join("\t", "total:", sum values %summary), "\n";
	}

sub core_var_copy_PA{
# determining which clusters are core, variable, and multi-copy #
	my ($PA_r, $toi_r) = @_;
	
		#print Dumper $PA_r; exit;
	
	my %core_var;
	foreach my $cluster (keys %$PA_r){
		# core/variable #
		my $ntaxa_ITEP = scalar keys %{$PA_r->{$cluster}};
		my $max_N = 0;
		
		# number of genes per taxon #
		foreach my $taxon (keys %{$PA_r->{$cluster}} ){
			# max number of ITEP genes 
			my $N_genes = scalar @{$PA_r->{$cluster}{$taxon}};
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
	
		#print Dumper %core_var; exit;
	return \%core_var;
	}

sub parse_PA{
# parsing PA table from ITEP #
	my ($infile, $taxa_r) = @_;

	open IN, $infile or die $!;

	my %PA;			# quantify as above
	my @header;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		# header #
		if($.==1){
			@header = split /\t/;
			}
		# body #
		else{
			my @line = split /\t/;
			
			for my $i (3..$#line){			# each organism
				# sanity check #
				die " ERROR: line $. is not the same length as the header!\n"
					unless $header[$i];
				
				# skipping taxa if don't care #
				if($taxa_r){
					next unless exists $taxa_r->{$header[$i]};
					}
					
				# laoding hash #
				my @entries = split /;/, $line[$i];
				
				foreach my $j (@entries){
					push(@{$PA{$line[1]}{$header[$i]}}, $j) unless $j eq "NONE";
					}			
				}
			}
		}
	close IN;
	
	# N-taxa of interest #
	my @toi;
	if($taxa_r){
		@toi = keys %$taxa_r;
		}
	else{ @toi = @header[3..$#header]; }
	
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


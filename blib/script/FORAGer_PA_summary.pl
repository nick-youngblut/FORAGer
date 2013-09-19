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

my ($verbose, $query_in, $subject_in);
GetOptions(
	   "query=s" => \$query_in,
	   "subject=s" => \$subject_in,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
# loading taxa of interest #
my $taxa_r = load_taxa_list($query_in) if $query_in;
my $subject_r = load_taxa_list($subject_in) if $subject_in;

# parsing PA tables #
my ($PA_ITEP_r, $toi_ITEP_r) = parse_PA($ARGV[0], $subject_r);
my ($PA_user_r, $toi_user_r) = parse_PA($ARGV[1], $subject_r);

# core/variable; multi-copy #
## determining core/variable & copy ##
my $core_var_ITEP_r = core_var_copy_PA($PA_ITEP_r, $toi_ITEP_r, $subject_r);
my $core_var_user_r = core_var_copy_PA($PA_user_r, $toi_user_r, $subject_r);

# sanity check #
die " ERROR: no taxon names in ITEP PA table matched taxa_list!\n" unless %$core_var_ITEP_r;
die " ERROR: no taxon names in user PA table matched taxa_list!\n" unless %$core_var_user_r;

## writing summary ##
write_core_var_summary($core_var_ITEP_r, "ITEP");
write_core_var_summary($core_var_user_r, "user");

## trim to just query taxa ##
trim_to_query($PA_ITEP_r, $taxa_r);
trim_to_query($PA_user_r, $taxa_r); 

## determining conflicts in core/variable ##
core_var_conflicts($PA_ITEP_r, $PA_user_r, $core_var_ITEP_r, $core_var_user_r);


### subroutines
sub trim_to_query{
	my ($PA_r, $taxa_r) = @_;

	foreach my $cluster (keys %$PA_r){
		foreach my $taxon (keys %{$PA_r->{$cluster}}){
			delete $PA_r->{$cluster}{$taxon} unless exists $taxa_r->{$taxon};
			unless( %{$PA_r->{$cluster}} ){
				delete $PA_r->{$cluster};
				last;
				}
			}
		}
	}

sub core_var_conflicts{
# determining conflicts in core/var & copy #
	my ($PA_ITEP_r, $PA_user_r, $core_var_ITEP_r, $core_var_user_r) = @_;

		#print Dumper $core_var_ITEP_r; exit;

	foreach my $cluster (keys %$PA_ITEP_r){
		unless( exists $PA_user_r->{$cluster} ){			# present in ITEP; not in user
			print join("\t", "Conflict:", $cluster, "Present", "Absent", "complete",
					$core_var_ITEP_r->{$cluster}{'core_var'}, 
					$core_var_ITEP_r->{$cluster}{'copy'}, "NA", "NA"), "\n";			
			}
		# present in both; checking for overlap of taxa #
		else{			
			# intersection of taxa names? # 
			my (%union, %isect);
			foreach my $e (keys %{$PA_ITEP_r->{$cluster}} , %{$PA_user_r->{$cluster}}){
				$union{$e}++ && $isect{$e}++; 
				}
			if(scalar keys %isect == 0){			# same, but not the same taxa
				print join("\t", "Cluster_overlap:", $cluster, "Present", "Present", "partial",
					$core_var_ITEP_r->{$cluster}{'core_var'}, 
					$core_var_ITEP_r->{$cluster}{'copy'},
					$core_var_user_r->{$cluster}{'core_var'}, 
					$core_var_user_r->{$cluster}{'copy'}), "\n";		
				}
			elsif(scalar keys %isect == scalar keys %{$PA_ITEP_r->{$cluster}} &&
				scalar keys %isect == scalar keys %{$PA_user_r->{$cluster}}){		# both have all of the taxa
				print join("\t", "Consistent:", $cluster, "Present", "Present", "complete",
					$core_var_ITEP_r->{$cluster}{'core_var'}, 
					$core_var_ITEP_r->{$cluster}{'copy'},
					$core_var_user_r->{$cluster}{'core_var'}, 
					$core_var_user_r->{$cluster}{'copy'}), "\n";					
				}
			else{
				print join("\t", "Taxa_overlap:", $cluster, "Present", "Present", "partial",
					$core_var_ITEP_r->{$cluster}{'core_var'}, 
					$core_var_ITEP_r->{$cluster}{'copy'},
					$core_var_user_r->{$cluster}{'core_var'}, 
					$core_var_user_r->{$cluster}{'copy'}), "\n";				
				}
			}
		}
	# present in user, but not ITEP?
	foreach my $cluster (keys %$PA_user_r){					
		unless( exists $PA_ITEP_r->{$cluster} ){
			print join("\t", "Conflict:", $cluster, "Absent", "Present", "complete",
					"NA", "NA",
					$core_var_user_r->{$cluster}{'core_var'}, 
					$core_var_user_r->{$cluster}{'copy'}), "\n";	
			
			}
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
	my ($PA_r, $toi_r, $taxa_r) = @_;
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
		my @line = split /\t/;
		
		$taxa{$line[0]} = 1;
		}
	close IN;

		#print Dumper %taxa; exit;
	return \%taxa;
	}

__END__

=pod

=head1 NAME

FORAGer_PA_summary.pl -- Compare ITEP & FORAGer presence-absence tables

=head1 SYNOPSIS

FORAGer_PA_summary.pl [options] ITEP_PA.txt FORAGer_PA.txt > summary.txt

=head2 options

=over

=item -query

A file with a list of taxon names that determine how pres-abs is called.

=item -subject

A file with a list of taxon names that determine how core-variable is called.

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc FORAGer_PA_summary.pl

=head1 DESCRIPTION

Compare ITEP and FORAGer gene cluster presence-absence tables.

=head2 Adding FORAGer pres-abs data to ITEP DB

FORAGer_screen.pl outputs a *PA.txt file.
Move *PA.txt file to $ITEP_HOME/userdata/ & rename
'user_genes'. Run ./main5.sh

ITEP table: db_getPresenceAbsenceTable.py -r mazei_I_2.0_c_0.4_m_maxbit -i

FORAGer table: db_getPresenceAbsenceTable.py -r mazei_I_2.0_c_0.4_m_maxbit -u

=head2 Output

=over

=item Conflict or consistent pres/abs between tables (also overlap of genes or just cluster_ID)?

=item Gene cluster ID

=item Present in ITEP table?

=item Present in user_genes table?

=item All genes shared in the cluster (complete|partial)?

=item ITEP: Core (in all genomes) or variable (only some genomes) gene?

=itme ITEP: Single copy or multi-copy gene?

=item User: Core (in all genomes) or variable (only some genomes) gene?

=itme User: Single copy or multi-copy gene?

=back

=head1 EXAMPLES

=head2 Basic usage:

FORAGer_PA_summary.pl <(db_getPresenceAbsenceTable.py -r mazei_I_2.0_c_0.4_m_maxbit -i ) <(db_getPresenceAbsenceTable.py -r mazei_I_2.0_c_0.4_m_maxbit -u) -q taxa_list.txt -s taxa_list.txt > summary.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/FORAGer/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


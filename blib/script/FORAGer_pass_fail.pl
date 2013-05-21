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

my ($verbose, $write_fails);
GetOptions(
	   "fail" => \$write_fails, 		# writing failed clusters
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
if($write_fails){
	write_fails();
	}
else{
	my $res_r = summing_pres_abs();
	write_pass_fail($res_r);
	write_summary($res_r);
	}

### Subroutines
sub write_summary{
	my ($res_r) = @_;
	
	my %summary;
	foreach my $query (keys %$res_r){
		foreach my $cluster (keys %{$res_r->{$query}}){
			$summary{$query}{"total"}++;
			$summary{$query}{"pass"} += $res_r->{$query}{$cluster};
			}
		}
	
	foreach my $query (keys %summary){
		print STDERR join("\t", $query, $summary{$query}{"pass"}, $summary{$query}{"total"},
				sprintf("%.2f", $summary{$query}{"pass"}/$summary{$query}{"total"} * 100)), "\n";
		}
	}

sub write_pass_fail{
# writing summary of pass fail #
	my ($res_r) = @_;

	foreach my $query (keys %$res_r){
		foreach my $cluster (keys %{$res_r->{$query}}){
			#print join("\t", $query, $cluster, $res_r->{$query}{$cluster}), "\n";
			}
		}
	}

sub summing_pres_abs{
# loading presence-absence table from FORAGer_screen #
	my %res;
	while (<>){
		chomp;
		my @line = split /\t/;
		if($line[3]){
			$res{$line[0]}{$line[1]} = 1;
			}
		else{
			$res{$line[0]}{$line[1]} = 0 unless $res{$line[0]}{$line[1]};
			}
		}
	return \%res;
	}

sub write_fails{
# writing failed clusters #
	my %res;
	my %fail_sum;
	while (<>){
		chomp;
		my @line = split /\t/;
		if($line[3]){
			$res{$line[0]}{$line[1]}{"PA"} = 1;
			}
		else{
			$res{$line[0]}{$line[1]}{"PA"} = 0 unless $res{$line[0]}{$line[1]};
			push(@{$res{$line[0]}{$line[1]}{"line"}}, $_);			# failed 
			}
		
		# summarizing what failed #
		$fail_sum{$line[0]}{"Nhit"}{"fail"}++ if $line[4] eq "FAILED";
		$fail_sum{$line[0]}{"Nhit"}{"NA"}++ if $line[4] eq "NA";
		$fail_sum{$line[0]}{"length"}{"fail"}++ if $line[7] eq "FAILED";
		$fail_sum{$line[0]}{"length"}{"NA"}++ if $line[7] eq "NA";
		$fail_sum{$line[0]}{"bitscore"}{"fail"}++ if $line[10] eq "FAILED";
		$fail_sum{$line[0]}{"bitscore"}{"NA"}++ if $line[10] eq "NA";
		}
	
	# writing out failed summary #
	foreach my $query (keys %fail_sum){
			$fail_sum{$query}{"Nhit"}{"fail"} = 0 unless $fail_sum{$query}{"Nhit"}{"fail"};
			$fail_sum{$query}{"length"}{"fail"} = 0 unless $fail_sum{$query}{"length"}{"fail"};
			$fail_sum{$query}{"bitscore"}{"fail"} = 0 unless $fail_sum{$query}{"bitscore"}{"fail"};
			$fail_sum{$query}{"Nhit"}{"NA"} = 0 unless $fail_sum{$query}{"Nhit"}{"NA"};
			$fail_sum{$query}{"length"}{"NA"} = 0 unless $fail_sum{$query}{"length"}{"NA"};
			$fail_sum{$query}{"bitscore"}{"NA"} = 0 unless $fail_sum{$query}{"bitscore"}{"NA"};
						
			
			print STDERR join("\t", $query, "FAIL", $fail_sum{$query}{"Nhit"}{"fail"},
						$fail_sum{$query}{"length"}{"fail"},
						$fail_sum{$query}{"bitscore"}{"fail"}), "\n";
			print STDERR join("\t", $query, "NA", $fail_sum{$query}{"Nhit"}{"NA"},
						$fail_sum{$query}{"length"}{"NA"},
						$fail_sum{$query}{"bitscore"}{"NA"}), "\n";
		}
	
	# writing out failed lines #
	foreach my $query (keys %res){
		foreach my $contig (keys %{$res{$query}}){
			print join("\n", @{$res{$query}{$contig}{"line"}}), "\n"
				unless $res{$query}{$contig}{"PA"};
			}
		}
	
	exit;
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


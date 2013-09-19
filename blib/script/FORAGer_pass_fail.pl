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

my ($verbose);
GetOptions(
	   "fail" => \&write_fails, 		# writing failed clusters
	   "pass" => \&write_passes, 			# writing passed clusters
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
my $res_r = summing_pres_abs();
write_pass_fail($res_r);
write_summary($res_r);


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
		if($line[3] eq "1"){
			$res{$line[0]}{$line[1]} = 1;
			}
		else{
			$res{$line[0]}{$line[1]} = 0 unless $res{$line[0]}{$line[1]};
			}
		}
		#print Dumper %res; exit;
	return \%res;
	}

sub write_fails{
# writing failed clusters #
	my %res;
	my %fail_sum;
	while (<>){
		chomp;
		my @line = split /\t/;
		if($line[3] eq "1"){
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
	
sub write_passes{
# writing passed clusters #
	my %res;
	my %pass_sum;
	while (<>){
		chomp;
		my @line = split /\t/;
		if($line[3] eq "1"){
			push(@{$res{$line[0]}{$line[1]}{"line"}}, $_);			# passed 
			}
		
		# summarizing what passed #
		$pass_sum{$line[0]}{"Nhit"}{"pass"}++ if $line[4] eq "PASSED";
		$pass_sum{$line[0]}{"Nhit"}{"NA"}++ if $line[4] eq "NA";
		$pass_sum{$line[0]}{"length"}{"pass"}++ if $line[7] eq "PASSED";
		$pass_sum{$line[0]}{"length"}{"NA"}++ if $line[7] eq "NA";
		$pass_sum{$line[0]}{"bitscore"}{"pass"}++ if $line[10] eq "PASSED";
		$pass_sum{$line[0]}{"bitscore"}{"NA"}++ if $line[10] eq "NA";
		}
	
	# writing out passed summary #
	foreach my $query (keys %pass_sum){
			$pass_sum{$query}{"Nhit"}{"pass"} = 0 unless $pass_sum{$query}{"Nhit"}{"pass"};
			$pass_sum{$query}{"length"}{"pass"} = 0 unless $pass_sum{$query}{"length"}{"pass"};
			$pass_sum{$query}{"bitscore"}{"pass"} = 0 unless $pass_sum{$query}{"bitscore"}{"pass"};
			
			$pass_sum{$query}{"Nhit"}{"NA"} = 0 unless $pass_sum{$query}{"Nhit"}{"NA"};
			$pass_sum{$query}{"length"}{"NA"} = 0 unless $pass_sum{$query}{"length"}{"NA"};
			$pass_sum{$query}{"bitscore"}{"NA"} = 0 unless $pass_sum{$query}{"bitscore"}{"NA"};
						
			
			print STDERR join("\t", $query, "PASS", $pass_sum{$query}{"Nhit"}{"pass"},
						$pass_sum{$query}{"length"}{"pass"},
						$pass_sum{$query}{"bitscore"}{"pass"}), "\n";
			print STDERR join("\t", $query, "NA", $pass_sum{$query}{"Nhit"}{"NA"},
						$pass_sum{$query}{"length"}{"NA"},
						$pass_sum{$query}{"bitscore"}{"NA"}), "\n";
		}
	
	# writing out passed lines #
	foreach my $query (keys %res){
		foreach my $contig (keys %{$res{$query}}){
			print join("\n", @{$res{$query}{$contig}{"line"}}), "\n";
			}
		}
	
	exit;
	}


__END__

=pod

=head1 NAME

FORAGer_pass_fail.pl -- investigate Passing/Failing of contigs

=head1 SYNOPSIS

FORAGer_screen.pl | FORAGer_pass_fail.pl [options]

=head2 options

=over

=item -fail

Write out all failed contig lines

=item -pass

Write out all passed contig lines

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc FORAGer_pass_fail.pl

=head1 DESCRIPTION

Investigate which contigs produced by the FORAGer pipeline
have passed or failed the screening processes. 

=head1 EXAMPLES

=head2 Usage method 1

FORAGer_pass_fail.pl -p < WWM610_pres-abs_summary.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/FORAGer/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


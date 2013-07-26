#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path;
use Parallel::ForkManager;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $qlist_in, $slist_in);
my $extra_params = "";
my $bowtie_params = "-k 10";
my $fork = 0;
GetOptions(
	   "query_list=s" => \$qlist_in,			# query directory
	   "subject_list=s" => \$slist_in,		# subject directory
	   "forks=i" => \$fork, 		# number of forked mappings
	   "params=s" => \$bowtie_params,	# params passed to bowtie2
	   "extra=s" => \$extra_params, 		# extra bowtie2 parameters (appended)
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
my $qdir = $ARGV[0];
my $sdir = $ARGV[1];
#die " ERROR: provide a query (read files) directory!\n" unless $qdir;
#die " ERROR: provide a subject (reference genomes) directory!\n" unless $sdir;
map{ $_ = File::Spec->rel2abs($_) if $_ } ($qdir, $sdir);
map{ die " ERROR: $_ not found!\n" if $_ && ! -d $_ } ($qdir, $sdir);

die " ERROR: provide a query list file!\n" unless $qlist_in;
die " ERROR: provide a subject list file!\n" unless $slist_in;
map{ die " ERROR: $_ not found!\n" unless -e $_ } ($qlist_in, $slist_in);

### MAIN
# loading list #
my $qlist_r = load_query_list($qlist_in);
my $slist_r = load_subject_list($slist_in);

# adding directory to files #
add_dir($qdir, $qlist_r) if $qdir;
add_dir($sdir, $slist_r) if $sdir;

# building indices for bowtie2 #
my $pmb = new Parallel::ForkManager($fork);
foreach my $fig (keys %$slist_r){
	my $pid = $pmb->start and next;
	die " ERROR: $slist_r->{$fig}{'ref'} not found!\n" unless -e $slist_r->{$fig}{"ref"};
	call_bowtie2_build($slist_r->{$fig}{"ref"});
	$pmb->finish;
	}
$pmb->wait_all_children;

# making output dir #
my $outdir = make_outdir();

# read mapping #
my $pm = new Parallel::ForkManager($fork);
foreach my $org (keys %$qlist_r){			# each set of read files
	foreach my $fig (keys %$slist_r){		# each subject genome
		
		#forking #
		my $pid = $pm->start and next;
	
		# bowtie2 mapping #
		my $sam_out = call_bowtie2($org, $slist_r->{$fig}{"ref"}, 
								$qlist_r, $bowtie_params, $extra_params,
								$outdir);
		
		# writing out 
		print join("\t", $sam_out, $fig, $org), "\n";
						
		# end fork #
		$pm->finish;
		}
	}
$pm->wait_all_children;


### Subroutines
sub make_outdir{
	my $outdir = File::Spec->rel2abs("./sam/");
	rmtree($outdir) if -d $outdir;
	mkdir $outdir or die $!;
	
	print STDERR "Writing files to: $outdir\n" unless $verbose;
	return $outdir;
	}

sub call_bowtie2{
	my ($org, $subject, $qlist_r, $bowtie_params, $extra_params, $outdir) = @_;
	
	my @sparts = File::Spec->splitpath($subject);
	
	my $cmd;
	if(exists $qlist_r->{$org}{"F"} && exists $qlist_r->{$org}{"R"}){	# 2 read files
		$cmd = "bowtie2 -x $subject -S $outdir/$org-$sparts[2].sam -1 $qlist_r->{$org}{'F'} -2 $qlist_r->{$org}{'R'} $bowtie_params $extra_params"; 
		}
	elsif(exists $qlist_r->{$org}{"FR"}){
		$cmd = "bowtie2 -x $subject -S $outdir/$org-$sparts[2].sam -U $qlist_r->{$org}{'FR'} $bowtie_params $extra_params"; 
		}
	print STDERR $cmd, "\n" unless $verbose;
	`$cmd`;
	
	return "$outdir/$org-$sparts[2].sam";
	}

sub call_bowtie2_build{
	my ($subject, $sdir) = @_;
	
	my $cmd = "bowtie2-build $subject $subject";
	print STDERR $cmd, "\n" unless $verbose;
	`$cmd`;
	}	

sub load_query_list{
	my ($list_in) = @_;
	
	open IN, $list_in or die $!;
	my %qlist;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;

		die " ERROR: $line[0] needs 1 or 2 reads file names in the query list file!\n"
			unless $line[1];
		if(! $line[2]){ $qlist{$line[0]}{"FR"} = $line[1]; }
		else{
			$qlist{$line[0]}{"F"} = $line[1];
			$qlist{$line[0]}{"R"} = $line[2];
			}
		}
	close IN;

		#print Dumper %qlist; exit;
	return \%qlist;
	}
	
sub load_subject_list{
# list = FIG -> fasta file name #
	my ($list_in) = @_;
	
	open IN, $list_in or die $!;
	my %slist;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;
		die " ERROR: provide a fasta file for $line[0] in subject list!\n"
			unless $line[1];
		
		$slist{$line[0]}{"ref"} = $line[1]; 		# FIG ID => fasta
		}
	close IN;
	
	return \%slist;
	}

sub add_dir{
# adding directory to file names #
	my ($dir, $list_r) = @_;

	foreach my $q (keys %$list_r){
		$list_r->{$q}{"F"} = join("/", $dir, $list_r->{$q}{"F"}) if exists $list_r->{$q}{"F"};
		$list_r->{$q}{"R"} = join("/", $dir, $list_r->{$q}{"R"}) if exists $list_r->{$q}{"R"};
		$list_r->{$q}{"FR"} = join("/", $dir, $list_r->{$q}{"FR"}) if exists $list_r->{$q}{"FR"};
		
		$list_r->{$q}{"ref"} = join("/", $dir, $list_r->{$q}{"ref"}) if exists $list_r->{$q}{"ref"};
		}
	}



__END__

=pod

=head1 NAME

FORAGer_mapper.pl -- Mapping reads from query genomes to all reference genomes

=head1 SYNOPSIS

FORAGer_mapper.pl [flags] [query_directory] [subject_directory] > sam_index.txt

=head2 Required flags

=over

=item -query

List file; tab-delimited; no header; columns:

=over

=item * query_organism_name 

=item * query_forward_read_file (or all reads file)

=item * query_reverse_read_file (or blank)

=back

=item -subject

List file; tab-delimited; no header; columns:

=over

=item * subject_FIG_ID

=item * subject_fasta

=back

=back

=head2 Optional flags

=over

=item -params

Parameters passed to bowtie2 (besides '-x, -S, -1, -2'). [-k 10]

=item -extra

Parameters appended to '-params'. []

=item -forks

Number of parallel calls of bowtie2. [1]

=item -verbose

Verbose output. [TRUE]

=item -help

This help message.

=back

=head2 For more information:

perldoc FORAGer_mapper.pl

=head1 DESCRIPTION

Just a wrapper for bowtie2 to facilitate mapping of the query reads
to all of the subject (reference) genomes.

Output SAM file naming: 'query'-'subject'.sam

By default, the top 10 hits for each read are kept (bowtie2 param: '-k 10'). 

=head2 'sam_index.txt'

3 columns: sam_file, subject_FIG_ID, query_organism_name

=head1 EXAMPLES

=head2 Basic Usage

FORAGer_mapper.pl -q ./query/ -s ./subject/ -list query_list.txt

=head2 Forking & multiple bowtie2 threads

FORAGer_mapper.pl -q ./query/ -s ./subject/ -list query_list.txt -f 5 -p "-k 10 -p 4"

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/ITEP_PopGen/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


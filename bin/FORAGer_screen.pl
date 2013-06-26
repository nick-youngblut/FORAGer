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
use List::Util qw/max min/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $nuc_clust_dir, $aa_clust_dir, @contig_dirs, $write_cluster, $header_bool);
my $fork = 0;
my $len_cutoff = 3;				# range expansion factor [3x]
my $floor = 300;				# min range (bp)
my $bit_cutoff = 0.4;			# homology cutoff
my $truncate = 0;				# truncating contig lengths to max tblastn range
my $out_prefix = "FORAGer";
my $runID = "NA";
my $query_name = "NA";
GetOptions(
	   "nuc=s{,}" => \$nuc_clust_dir,
	   "aa=s{,}" => \$aa_clust_dir,
	   "contigs=s{,}" => \@contig_dirs,
	   "length=f" => \$len_cutoff, 			# cutoff length for a contig (%)
	   "bitscore=f" => \$bit_cutoff, 		# normalized bitscore cutoff
	   "minimum=i" => \$floor,				# min range (bp)
	   "write" => \$write_cluster, 			# write cluster w/ contig
	   "truncate=s" => \$truncate, 			# truncating 'passed' contigs to max cluster sequence length [TRUE]
	   "x" => \$header_bool, 				# write header? [FALSE]
	   "runID=s" => \$runID,				# ITEP runID
	   "query=s" => \$query_name, 			# query genome name
	   "prefix=s" => \$out_prefix, 			# outfile prefix
	   "fork=i" => \$fork,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a directory containing gene cluster nucleotide fasta files!\n"
	unless $nuc_clust_dir;
die " ERROR: provide a directory containing gene cluster amino acid fasta file!\n"
	unless $aa_clust_dir;
die " ERROR: provide >=1 directory containing FORAGer contig files (1 directory per query organism)!\n"
	unless @contig_dirs;
print STDERR " WARNING: the ITEP run_ID used for obtaining the clusters was not provided! P/A table field with be NA's!\n"
	unless $runID ne "NA";
print STDERR " WARNING: the query genome FIG_ID was not provided! P/A table field with be NA's!\n"
	unless $query_name ne "NA";
map{ die " ERROR: $_ not found!\n" unless -d $_; $_=File::Spec->rel2abs($_)} ($nuc_clust_dir, $aa_clust_dir, @contig_dirs);

### MAIN
write_header() if $header_bool;

my $pm = new Parallel::ForkManager($fork);
foreach my $contig_dir (@contig_dirs){
	# loading files #
	my $files_r = get_file_names($nuc_clust_dir, $aa_clust_dir, $contig_dir);
	die " ERROR: no files found!\n" unless %$files_r;

	# outdir for fasta files of passed contigs #
	my $outdir = make_outdir($contig_dir, "_passed");

	# opening output file handles #
	my $filter_out_name = "$out_prefix\_screen.txt";
	open my $filter_fh, ">$filter_out_name" or die $!;
	my $PA_out_name = "$out_prefix\_PA.txt";
	open my $PA_fh, ">$PA_out_name" or die $!;

	# blasting, filtering, writing results #
	foreach my $clust_file (keys %$files_r){			# foreach cluster: basename of all 3 needed files
		$pm->start and next;
		
		# if no contig file; filter table gets 'NA' #
		if($files_r->{$clust_file} eq "NA"){		
			write_filter_table($clust_file, 0, $contig_dir, $filter_fh);
			$pm->finish;
			next;
			}
		
		# blasting #
		## tblastn cluster vs contig ##
		my $tblastn_r = tblastn_cluster_contig($clust_file, $files_r->{$clust_file}, $aa_clust_dir, $contig_dir);
	
		## tblastn cluster (AA) vs cluster (nuc) ##
		my $tblastn_self_r = self_tblastn($clust_file, $aa_clust_dir, $nuc_clust_dir);
	
		## normalizing bit score: bit(contig vs gene) / max-bit(self-contig | self-gene) ##
		norm_blast($tblastn_r, $tblastn_self_r);
	
		# loading fasta of contigs & clusters #
		my $contigs_r = load_fasta($files_r->{$clust_file}, $contig_dir);
		my $clusters_r = load_fasta($clust_file, $nuc_clust_dir);
	
		# filtering #
		my %summary;
		## filtering by number of tblastn hits (each gene must hit the contig) ##
		filter_by_Nhits(scalar keys %$clusters_r, $contigs_r, $tblastn_r, \%summary);
		
		## filtering by length cutoff ##	
		my $clust_range_r = get_clust_len_range($clusters_r, $len_cutoff, $floor);
		filter_by_length($clust_range_r, $contigs_r, $clusters_r, $tblastn_r, \%summary);	
	
		## filtering by score and length ##
		filter_by_bitscore($contigs_r, $tblastn_r, $bit_cutoff, \%summary);
	
		## writing out contig filter table ##
		my $passed = write_filter_table($clust_file, \%summary, $contig_dir, $filter_fh);
		
		## truncating contigs by max tblastn hit to cluster ##
		truncate_contigs($contigs_r, $tblastn_r, $truncate) if $truncate =~ /^\d+$/;		# must be integer
		
		## writing out passed contigs ##
		if($passed){			# writing out passed contig files
			my $passed_file_name = write_passed_contig_fasta(\%summary, $contigs_r, $files_r->{$clust_file}, $outdir);
			append_cluster_to_contig($nuc_clust_dir, $clust_file, $outdir, $passed_file_name) if $write_cluster;
			}

		## writing PA table (for ITEP) ##
		write_PA_table(\%summary, $clust_file, $contigs_r, $runID, $query_name, $PA_fh) if $passed;
		
		$pm->finish;
		}
	
	# closing file handles #
	close $filter_fh;
	print STDERR " Filter summary file written: $filter_out_name\n";
	close $PA_fh;
	print STDERR " Presence-absence file written: $PA_out_name\n";
	}
$pm->wait_all_children;


### Subroutines
sub write_PA_table{
	my ($summary_r, $clust_file, $contigs_r, $runID, $query_name, $PA_fh) = @_;
	
	# cluster ID #
	(my $cluster_id = $clust_file) =~ s/.*clust(\d+).+/$1/;
	
	# writing entries for PA table #
	foreach my $passed_contig (keys %$summary_r){
		if($summary_r->{$passed_contig}{"PA"}){		# passed; need to write out
			die " ERROR: $passed_contig not found in contig file!\n" 
				unless exists $contigs_r->{$passed_contig};
			
			my $user_geneid = join("__", "FORAGer", $passed_contig, time());
			print $PA_fh join("\t", $user_geneid, $query_name, "CDS",
					"", "", "",
					$runID, $cluster_id, $contigs_r->{$passed_contig},
					""), "\n";		
			}
		}	
	close OUT;
	}

sub write_header{
# writing the header to the PA table #
	my @header = qw|Query_file_name Cluster_file_name Query_contig_name Presence/Absence_(binary) Number_tblastn_hits_PASS/FAIL Number_of_tblastn_hits Number_of_genes_in_cluster Length_PASS/FAIL tblastn_hit_length_range Cluster_length_range Normalized_bit_score_PASS/FAIL Minimum_normalized_bit_score Normalized_bit_score_cutoff|;
	print join("\t", @header), "\n";
	}

sub truncate_contigs{
# truncating contigs to max tblastn length #
	my ($contigs_r, $tblastn_r, $truncate) = @_;
	
	foreach my $contig (keys %$contigs_r){ 		
		# getting max hit length #
		my $max_len = 0;
		my @start_stop;
		foreach my $query (keys %{$tblastn_r->{$contig}}){				
			my $start = ${$tblastn_r->{$contig}{$query}}[0]; 
			my $end = ${$tblastn_r->{$contig}{$query}}[1];
			my $hit_len = abs($end - $start);		# length in nuc				
			if($max_len < $hit_len){				# if longer hit
				if($start < $end){					# start-end
					@start_stop = ($start - 4 - $truncate, $end - $start + $truncate);		# start (0-indexed + 1 codon, length)
					}
				else{								# end-start
					@start_stop = ($end - 4 - $truncate, $start - $end + $truncate);
					}
				$start_stop[0] = 0 if $start_stop[0] < 0;
				$max_len = $hit_len;
				}
			}
		
		# truncating #
		$contigs_r->{$contig} = substr($contigs_r->{$contig}, $start_stop[0], $start_stop[1])
			if @start_stop;
		}
	}

sub make_outdir{
# making output directory #
	my ($outdir_name, $append) = @_;
	
	$outdir_name .= "$append";
	
	$outdir_name = File::Spec->rel2abs($outdir_name);

	rmtree($outdir_name) if -d $outdir_name;
	mkdir $outdir_name or die $!;
	
	return $outdir_name;
	}

sub append_cluster_to_contig{
# appending clusters (nuc) to contig file #
	my ($nuc_clust_dir, $clust_file, $outdir, $passed_file_name) = @_;
	
	my $cmd = "cat $nuc_clust_dir/$clust_file >> $outdir/$passed_file_name";
	print "$cmd\n" unless $verbose;
	`$cmd`;
	}

sub write_passed_contig_fasta{
# writing out a fasta of cluster & contig #
	my ($summary_r, $contigs_r, $contig_file, $outdir) = @_;

	(my $outfile = $contig_file) =~ s/\.[^\.]{1,6}$|$/_pass.fasta/; 
	open OUT, ">$outdir/$outfile" or die $!;

	foreach my $contig (keys %$summary_r){
		if($summary_r->{$contig}{"PA"}){		# passed; need to write out
			die " ERROR: $contig not found in contig file!\n" 
				unless exists $contigs_r->{$contig};
			print OUT join("\n", ">$contig", $contigs_r->{$contig}), "\n";
			}
		}
	
	close OUT;
	return $outfile;
	}
	
sub write_filter_table{
# writing out PA table to STDOUT #
	my ($clust_file, $summary_r, $contig_dir, $filter_fh) = @_;
	
	my @stats = qw/PA N_tblastn_hits_cutoff N_tblastn_hits N_cluster_genes length_cutoff hit_length_range cluster_length_range norm_bit_score min_bit_score bit_score_cutoff/;
	
	# writing body #
	my $passed;
	if($summary_r){	
		foreach my $contig (keys %$summary_r){
			print $filter_fh join("\t", $contig_dir, $clust_file, $contig);			# query, cluster file name, contig_name
			map{exists $summary_r->{$contig}{$_} ? print $filter_fh "\t$summary_r->{$contig}{$_}" : print $filter_fh "\tNA" } @stats;
			print $filter_fh "\n";
			
			$passed = 1 if $summary_r->{$contig}{"PA"};
			}
		}
	else{			# no contig file, failed assembly
		print $filter_fh join("\t", $contig_dir, $clust_file, "NO_CONTIG_FILE", ("NA") x scalar @stats), "\n";		 
		}
		
	return $passed;
	}

sub filter_by_bitscore{
# filtering by normalized bitscore #
	my ($contigs_r, $tblastn_r, $bit_cutoff, $summary_r) = @_;
	
	foreach my $subject (keys %$contigs_r){		# each contig
		if($bit_cutoff >= 0 && exists $tblastn_r->{$subject}){
			my @norm_bits;
			foreach my $query (keys %{$tblastn_r->{$subject}}){		# getting norm bit scores
				push(@norm_bits, ${$tblastn_r->{$subject}{$query}}[3]);
				}
			# pass/fail #
			my $x = min(@norm_bits);
			if($x >= $bit_cutoff){
				$summary_r->{$subject}{"norm_bit_score"} = "PASSED";
				}
			else{
				$summary_r->{$subject}{"PA"} = 0;
				$summary_r->{$subject}{"norm_bit_score"} = "FAILED";
				}
			$summary_r->{$subject}{"min_bit_score"} = sprintf("%.3f", $x);
			$summary_r->{$subject}{"bit_score_cutoff"} = sprintf("%.3f", $bit_cutoff);
			}
		else{
			$summary_r->{$subject}{"norm_bit_score"} = "NA";
			$summary_r->{$subject}{"min_bit_score"} = "NA";
			$summary_r->{$subject}{"PA"} = 0 if $bit_cutoff >= 0; 		# don't fail if not filtering by bit score
			}
		
		# PRESENT if passed all filters (no '0' for 'PA') #
		$summary_r->{$subject}{"PA"} = 1 unless exists $summary_r->{$subject}{"PA"};
		}
	
		#print Dumper %$summary_r; 
	}

sub filter_by_length{
# filtering the contigs by length relative to cluster #
	my ($clust_range_r, $contigs_r, $clusters_r, $tblastn_r, $summary_r) = @_;

	foreach my $contig (keys %$contigs_r){ 							# checking contig
		if($len_cutoff && exists $tblastn_r->{$contig}){			# if tblastn hit(s)
			
			# checking all tblastn hit lengths #
			my $N_passed = 0;		# default = fail
			my @hit_lens;
			foreach my $query (keys %{$tblastn_r->{$contig}}){				# 1 hit per gene in cluster must meet length requirement
				my $hit_len = abs(${$tblastn_r->{$contig}{$query}}[1] - ${$tblastn_r->{$contig}{$query}}[0]);		# length in nuc				
				$N_passed++ if $hit_len >= $$clust_range_r[0] &&			# expanding negative range
							   $hit_len <= $$clust_range_r[1];				# expanding positive range
				push(@hit_lens, $hit_len);
				}
				
			## hit length range ##
			$summary_r->{$contig}{"hit_length_range"} = join(":", sprintf("%.0f", min(@hit_lens)), 
												sprintf("%.0f", max(@hit_lens)));		
			
			# pass/fail length cutoff #
			if($N_passed == scalar keys %$clusters_r){		# all hits must pass
				$summary_r->{$contig}{"length_cutoff"} = "PASSED";
				}
			else{
				$summary_r->{$contig}{"PA"} = 0;
				$summary_r->{$contig}{"length_cutoff"} = "FAILED";
				}
			}
		else{								# no tblastn hits at all
			$summary_r->{$contig}{"length_cutoff"} = "NA";
			$summary_r->{$contig}{"hit_length_range"} = "NA";			
			$summary_r->{$contig}{"PA"} = 0 if $len_cutoff; 		# no tblastn hits, so not PASSING
			}
		
		## cluster length range ##
		$summary_r->{$contig}{"cluster_length_range"} = join(":", sprintf("%.0f", $$clust_range_r[0]), 
															sprintf("%.0f",$$clust_range_r[1]));	
		}
		#print Dumper $summary_r; exit;
	}

sub filter_by_Nhits{
# contig must be hit (tblastn) by all genes in cluster #
	my ($N_genes, $contigs_r, $tblastn_r, $summary_r) = @_;
	
	foreach my $contig (keys %$contigs_r){
		if(exists $tblastn_r->{$contig}){			# if tblastn hit(s)
			my $N_hits = scalar keys %{$tblastn_r->{$contig}};
			if($N_hits < $N_genes){
				$summary_r->{$contig}{"PA"} = 0;
				$summary_r->{$contig}{"N_tblastn_hits_cutoff"} = "FAILED";		
				}
			else{
				$summary_r->{$contig}{"N_tblastn_hits_cutoff"} = "PASSED";		
				}
			$summary_r->{$contig}{"N_tblastn_hits"} = $N_hits;		
			}
		else{
			$summary_r->{$contig}{"PA"} = 0;
			$summary_r->{$contig}{"N_tblastn_hits_cutoff"} = "FAILED";				
			$summary_r->{$contig}{"N_tblastn_hits"} = 0;		
			}
		$summary_r->{$contig}{"N_cluster_genes"} = $N_genes;
		}
	}

sub get_clust_len_range{
	my ($clusters_r, $len_cutoff, $floor) = @_;
	
	# getting lengths #
	my @lens;
	map{ push(@lens, length($_)) } values %$clusters_r;

	# getting stdev of lengths #
	my $N = scalar @lens;
	if($N > 1){  			# if >1 gene in cluster
		my $min = min(@lens) - 3;
		my $max = max(@lens);
		my $range = $max - $min;
		$range = $floor if $range < $floor; 							# floor of range
		my $exp = (($range * $len_cutoff) - $range)/2;					# expansion of range
		$exp += ($floor - ($max - $min))/2 if $max - $min < $floor; 	# adding to min/max of range
			#print Dumper $range;
			#print Dumper $exp; 
			#print Dumper $min - $exp, $max + $exp if ($max + $exp) - ($min - $exp) < 9;
		return [$min - $exp, $max + $exp];						# range +/- expansion
		}		# if >1 peg in cluster, max - min
	else{ return [($lens[0] -3) * 0.5, $lens[0] * 1.5]; }	# just length of the gene +/- 50%
	}

sub load_fasta{
	# version: 2.0
	# usage: load_fasta($fasta_file_name); returns a hash of sequences
	my ($fasta_in, $dir) = @_;
	open IN, "$dir/$fasta_in" or die $!;
	my (%fasta, $tmpkey);
	while(<IN>){
		chomp;
 		 s/#.+//;
 		next if  /^\s*$/;	
 		if(/>.+/){
 			s/^>//;
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{$fasta{$tmpkey} .= $_; }
		}
	close IN;
		#print Dumper(%fasta); exit;
	return \%fasta;
	}

sub norm_blast{
# normalizing all tblastn bit scores: bit(contig vs gene) / max-bit(self-contig | self-gene) #
	my ($tblastn_r, $tblastn_self_r, $blastn_res_r) = @_;
	
	foreach my $subject (keys %$tblastn_r){
		foreach my $query (keys %{$tblastn_r->{$subject}}){
			# sanity check #
			die " ERROR: cannot find self blastp hit for $query!\n"
				unless exists $tblastn_self_r->{$query};

			# normalizing #
			${$tblastn_r->{$subject}{$query}}[3] = ${$tblastn_r->{$subject}{$query}}[3] / 
				$tblastn_self_r->{$query};
			}
		}
		#print Dumper %$tblastn_r; exit;
	}

sub self_tblastn{
# tblastn against self (aa vs nuc); finding self hit w/ highest bit score #
	my ($file, $aa_dir, $nuc_dir) = @_;
	
	my $cmd = "tblastn -query $aa_dir/$file -subject $nuc_dir/$file  -soft_masking true -outfmt '6 qseqid sseqid sstart send evalue bitscore sframe'";
	print STDERR $cmd, "\n" unless $verbose;
	open PIPE, "$cmd |" or die $!;	
	
	my %blast_clust;
	while(<PIPE>){					# returns length of gene in bp
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;
		next unless $line[0] eq $line[1];
		$line[5] =~ s/\s+//g;
		if(exists $blast_clust{$line[0]}){ 		# getting max bit score if multiple self hits
			$blast_clust{$line[0]} = $line[5] if $line[5] > $blast_clust{$line[0]};
			}
		else{ 			# getting bit score if 1st self hit
			$blast_clust{$line[0]} = $line[5];
			}
		}
	close PIPE;

		#print Dumper $min_bit; exit;
	return \%blast_clust;	#query=>subject=>blast_res
	}

sub tblastn_cluster_contig{
# tblastn of contig vs cluster; parsing results #
	my ($clust_file, $contig_file, $clust_dir, $contig_dir) = @_;
	
	my $cmd = "tblastn -query $clust_dir/$clust_file -subject $contig_dir/$contig_file -soft_masking true -outfmt '6 qseqid sseqid sstart send evalue bitscore sframe'";
	
	print STDERR $cmd, "\n" unless $verbose;
	open PIPE, "$cmd |" or die $!;
	
	my %blast_res;
	while(<PIPE>){
		chomp;
			#print $_, "\n";
		my @line = split /\t/;
		if( exists $blast_res{$line[1]}{$line[0]} ){			
			$blast_res{$line[1]}{$line[0]} = [@line[2..$#line]] 
				if $line[5] > ${$blast_res{$line[1]}{$line[0]}}[3];
			}
		else{ $blast_res{$line[1]}{$line[0]} = [@line[2..$#line]]; }
		}
	close PIPE;
		#print Dumper %blast_res; exit;
	return \%blast_res;		# query=>subject=>blast_res
	}

sub get_file_names{
	my ($nuc_clust_dir, $aa_clust_dir, $contig_dir) = @_;

	# getting nucleotide cluster fasta files #
	opendir DIR, $nuc_clust_dir or die $!;
	my @nuc_clust_files = grep(/clust\d+\.(fna|fasta)$/, readdir DIR);
	closedir DIR;
	die " ERROR: no cluster fasta files (*fna|*fasta) found in $nuc_clust_dir!\n"
		unless @nuc_clust_files;
		
	# getting AA cluster fasta files #
	opendir DIR, $aa_clust_dir or die $!;
	my @aa_clust_files = grep(/clust\d+\.(fa|faa|fasta)$/, readdir DIR);
	closedir DIR;
	die " ERROR: no cluster fasta files (*fa|*faa|*fasta) found in $aa_clust_dir!\n"
		unless @aa_clust_files;	
		
	# sanity check #
	my $nuc_cnt = scalar @nuc_clust_files;
	my $aa_cnt = scalar @aa_clust_files;
	die " ERROR: the number of nucleotide and AA fasta cluster files doesn't match ($nuc_cnt vs $aa_cnt)!\n"
		unless $nuc_cnt == $aa_cnt;
	
	# getting contig files #
	opendir DIR, $contig_dir or die $!;
	my @contig_files = grep(/-contig\.(fna|fasta)$/, readdir DIR);
	closedir DIR;
	die " ERROR: no contig files found in $contig_dir!\n" 
		unless @contig_files;
	
	# making hash #
	my %files;
	foreach my $cluster (@nuc_clust_files){
		(my $q = $cluster) =~ s/\.[^.]+$//;				# stripping extension
		my @hits = grep(/^$q\_/, @contig_files);

		if(@hits){	# if >= 1 hit
			die " ERROR: multiple contig files found for $q\n" if scalar @hits > 1;
			$files{$cluster} = $hits[0];
			}
		else{		# if contig file not found 
			print STDERR " WARNING: no contig file found for $q! Failed assembly due to low coverage???\n" unless @hits || $verbose;
			$files{$cluster} = "NA";
			}
		}
		
		#print Dumper %files; exit;
	return \%files;		# cluster_file => contig_file
	}


__END__

=pod

=head1 NAME

FORAGer_screen.pl -- screening contigs to determine homology to the target gene cluster

=head1 SYNOPSIS

FORAGer_screen.pl [flags] > pres-abs_summary.txt

=head2 Required flags

=over

=item -contig

>=1 diretories of contig files produced by FORAGer_assemble.pl

=item -nuc

>=1 diretories of cluster files (nucleotide) produced by FORAGer.pl 

=item -aa

>=1 diretories of cluster files (amino acid) produced by FORAGer.pl 

=back

=head2 Semi-optional flags

=over

=item -runID

The ITEP runID used to make the clusters. (NA's written to P/A table if not provided)

=item -query 

The FIG_ID of the query genome. (NA's written to P/A table if not provided)

=back

=head2 Optional flags

=over

=item -prefix

Output file prefix. [FORAGer]

=item -bitscore

Normalized bit score cutoff (negative value to skip filtering). [0.4]

=item -length

Length range cutoff (min-max range of genes in the gene cluster * '-length'). 
'-length 0' skips filtering. [3]

=item -minimum

The minimum length range for filtering (bp). [300]

=item -fork

Number of parallel cluster comparisons to perform. [1]

=item -truncate

Truncate contigs to the longest tblastn hit ('F'=FALSE; provided value added to range)? [0] 

=item -write

Write clusters with passed contigs? [FALSE]

=item -x 	Write header for presence-absence table? [FALSE]

=item -v	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc FORAGer_screen.pl

=head1 DESCRIPTION

Determine which contigs produced by FORAGer_assemble.pl actually
should be included in the target gene cluster. Filtering is based on
sequence length & homology.

=head2 Normalized bit score cutoff

Normalized bit score = (tblastn_bit_score gene vs contig) / (blastp_bit_score gene vs gene).
The contig must hit all genes in cluster with a normalized bit score >= the cutoff.
The cutoff should probably be the same as used for the original gene clustering.

=head2 Length cutoff

The length range is defined as the min-max of gene lengths (bp) in the cluster.
'-length' is a scaling factor for how much to expand or contract this range.
By default, min-max length range values are used for the cutoff. If only one gene
in the cluster, the min-max range is 0.5*gene_length to 1.5*gene_length.

=head2 Output: *screen.txt file

Screening summary describing which filters were passed/failed. 
"NA" = not applicable.
"NO_CONTIG_FILE" = no contigs produced by the assembler. Columns:

=over

=item * 	Query file name

=item * 	Cluster file name

=item * 	Query contig name

=item * 	Presence/Absence (binary)

=item * 	Number tblastn hits PASS/FAIL

=item * 	Number of tblastn hits

=item * 	Number of genes in cluster

=item * 	Length PASS/FAIL

=item * 	tblastn hit length range (min:max)

=item * 	Cluster length range (min:max * '-length')

=item * 	Normalized bit score PASS/FAIL

=item * 	Minimum normalized bit score

=item * 	Normalized bit score cutoff

=back

=head2 Output: *PA.txt file

Presence/Absence file for updating ITEP.

=over

=item * 	user_geneid

=item * 	organismid

=item * 	genetype

=item * 	contigid

=item * 	startloc

=item * 	stoploc

=item * 	runid

=item * 	clusterid

=item * 	seq

=item * 	type

=item * 	annotation

=back

=head1 EXAMPLES

=head2 Basic Usage:

FORAGer_screen.pl -contig Mapped2Cluster_query -nuc cluster_nuc -aa cluster_aa -query 2209.25 -run mazei_I_2.0_c_0.4_m_maxbit

=head2 Append clusters to contigs:

FORAGer_screen.pl -contig Mapped2Cluster_query -nuc cluster_nuc -aa cluster_aa -w -query 2209.25 -run mazei_I_2.0_c_0.4_m_maxbit

=head2 Do not truncate contigs:

FORAGer_screen.pl -contig Mapped2Cluster_query -nuc cluster_nuc -aa cluster_aa -w -t F -query 2209.25 -run mazei_I_2.0_c_0.4_m_maxbit

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/FORAGer/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


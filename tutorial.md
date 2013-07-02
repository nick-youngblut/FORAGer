FORAGer tutorial
================

FORAGer example
---------------

### Set up
FORAGer_mapper.pl needs the following files & directories to begin:

* A 'query_list' file; tab-delimited; 3 columns: 

	* query_organism_name

	* query_forward_read_file (or all reads file)

	* query_reverse_read_file (or blank)

* A 'subject_list' file; tab-delimited; 2 columns:
           
	* subject_FIG_ID

	* subject_fasta
	
The query list associates the query organism with the reads used for
mapping and also tells FORAGer where to find the read files.

The subject list associates the subject organism FIG (needed for ITEP)
with its genome sequence file.


#### example of a query_list file

>Isolate1	Isolate1_F.fq	Isolate1_R.fq

>Isolate2	Isolate2_F.fq	Isolate2_R.fq

>Isolate3	Isolate3_F.fq	Isolate3_R.fq


#### example of a subject_list file

>2209.24	E_coli.fna

>2209.27	S_enterica.fna

>213585.6	B_subtilis.fna



### Example files for this exercise:

* query list = 'query_list.txt'

* subject list = 'subject_list.txt'

* query read files = 'query_reads/'

* subject genome files = 'subject_genomes/'


### Mapping reads






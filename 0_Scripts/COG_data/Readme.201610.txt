############################################################
0.	General remarks
############################################################

This is a December 2014 release of 2003-2014 COGs constructed by
Eugene Koonin's group at the National Center for Biotechnology
Information (NCBI), National Library of Medicine (NLM), National
Institutes of Health (NIH).

#-----------------------------------------------------------
0.1.	Citation

Galperin MY, Makarova KS, Wolf YI, Koonin EV.

Expanded microbial genome coverage and improved protein family
annotation in the COG database.

Nucleic Acids Res. 43, D261-D269, 2015
<http://www.ncbi.nlm.nih.gov/pubmed/25428365>

#-----------------------------------------------------------
0.2.	Contacts

<cogs@ncbi.nlm.nih.gov>

############################################################
1.	Notes
############################################################

#-----------------------------------------------------------
1.1.	2003-2014 COGs

This release contains 2003 COGs assigned to a representative set of
bacterial and archaeal genomes, available at February 2014. No new
COGs were constructed.

#-----------------------------------------------------------
1.2.	GIs and Refseq IDs

Sequences in COGs are identified by GenBank GI numbers. GI numbers
generally are transient. There are two ways to make a more permanent
link between the protein in COGs and the outside databases: via the
RefSeq accession codes (see 2.5) and via the protein sequences (see
2.6).

Note, however, that at the moment (April 02, 2015) RefSeq database is
in a state of transition; some of the <refseq-acc> entries are not
accessible. This accession table will be updated as soon as RefSeq is
stable.

############################################################
2.	Files
############################################################

#-----------------------------------------------------------
2.1.	genomes2003-2014.tab

Contains list of genomes (711). Tab-delimited, format:

<genome-code> <ncbi-taxid> <ncbi-ftp-name>

* Example:

Acihos	933801	Acidianus_hospitalis_W1_uid66875

* Comments:

The field <genome-name> serves as a subdirectory name on NCBI Genomes
FTP site 

ftp://ftp.ncbi.nih.gov/genomes/Bacteria/[ncbi-ftp-name]

where the genome sequence was availavle as of February 2014. The
field <ncbi-taxid> is a semicolon-delimited lineage according to the
NCBI Taxonomy database (as of February 2014).

#-----------------------------------------------------------
2.2.	cognames2003-2014.tab

Contains list of COG annotations. Tab-delimited, format:

<COG-id> <functional-class> <COG-annotation>

* Example:

COG0001	H	Glutamate-1-semialdehyde aminotransferase

* Comments:

Functional classes (categories) are described in the file
fun2003-2014.tab (see 2.3). Some COGs belong to more than one
functional class; in these cases the class listed first is considered
to be primary.

#-----------------------------------------------------------
2.3.	fun2003-2014.tab

Contains list of functional classes. Tab-delimited, format:

<class-id> <description>

* Example:

J	Translation, ribosomal structure and biogenesis

* Comments:

N/A

#-----------------------------------------------------------
2.4.	cog2003-2014.csv

Contains list of orthology domains. Comma-delimited, format:

<domain-id>, <genome-name>, <protein-id>,<protein-length>,
<domain-start>, <domain-end>, <COG-id>, <membership-class>,

* Example:

333894695,Alteromonas_SN2_uid67349,333894695,427,1,427,COG0001,0,

* Comments:

In this version the fields <domain-id> and <protein-id> are identical
and both normally refer to GenBank GIs. Thus neither <domain-id> nor
<protein-id> are necessarily unique in this file (this happens when a
protein consists of more than one orthology domains, e.g. 48478501).

The <membership-class> field indicates the nature of the match
between the sequence and the COG consensus:

0 - the domain matches the COG consensus;

1 - the domain is significantly longer than the COG consensus;

2 - the domain is significantly shorter than the COG consensus;

3 - partial match between the domain and the COG consensus.

#-----------------------------------------------------------
2.5.	prot2003-2014.tab

Contains RefSeq accession codes for all proteins with assigned COG
domains. Format:

<protein-id> <refseq-acc>

* Example:

103485499	YP_615060

* Comments:

As of September 2016 RefSeq and Genbank database has changed; many of
the <refseq-acc> entries are not accessible. This table is updated;
see prot2003-2014.gi2gbk.tab (2.6).

#-----------------------------------------------------------
2.6.	prot2003-2014.gi2gbk.tab

Contains RefSeq and GenBank accession codes for all proteins with
assigned COG domains. Format:

<protein-id> <refseq-acc> <genbank-acc>

* Example:

103485499	YP_615060	ABF51727

* Comments:

As of now (October 13, 2016) 5391 out of 1785722 proteins miss
GenBank accession codes. We will upload the updated table as soon as
possible.

#-----------------------------------------------------------
2.7.	prot2003-2014.fa.gz

* Contains:

Sequences of all proteins with assigned COG domains in FASTA format
(gzipped)

* Example:

>gi|118430838|ref|NP_146899.2| putative mercury ion binding protein
[Aeropyrum pernix K1]
MIIFKRHSQAILFSHNKQEKALLGIEGMHCEGCAIAIETALKNVKGIIDTKVNYSRGSAI
VTFDDTLVSINDILEHYIFKVPSNYRAKLVSFIS

* Comments:

The first word of the defline always starts with "gi|<protein-id>".

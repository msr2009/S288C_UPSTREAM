S288C REGIONS UPSTREAM OF CODING SEQUENCES
------------------------------------------

To the best of my knowledge, no database stores the sequences of ALL regions of
the yeast genome that are upstream of a gene. There are many transcription
factor or promoter databases (like SCPD, YEASTRACT, YeTFaSCo, etc) but these
are curated to include well-annotated promoters and factors, rather than simply
the complete set of "promoters" located in the genome.

Because we don't know the functional extent of every gene's promoter, I decided
to extract sequences in the genome that are upstream of a gene, extending to
the next closest transcribed feature or autonomously replicating sequence. To
illustrate, I've drawn a couple examples below.

               pGENE_A                   pGENE_B
[GENE_X>>>]~~~~~~~~~~~~~~~|>>>GENE_A|===============|>>>GENE_B|

This is an example of three genes, all on the same strand. The region from 
directly upstream of the start codon for GENE_A to the stop codon of GENE_X
(~~~) is extracted as pGENE_A, and the region extending from GENE_B to the stop
codon of GENE_A (===) is extracted as pGENE_B

                   pGENE_D
|GENE_C<<<|=======================|>>>GENE_D|
                   pGENE_C

In this example, we have a intergenic region driving bidirectional
transcription (e.g., the GAL1/10 promoter). This interval is counted twice,
once for each strand. As such, the FASTA records for pGENE_C and pGENE_D are
reverse complements of one another.

Because I only am interested in regions upstream of genes, the region between
two convergently transcribed genes are not reported.

All sequences reported here span the first base upstream of the start codon of
a gene to the first base before the next annotated genomic feature, whether
that is the start codon of a neighboring gene (as in the second example above),
the stop codon of a neighboring gene (like GENE_B in the first example), or
another feature (GENE_A in example 1).

The general workflow to extract these sequences was:

	1) 	Filter GFF files from SGD for only genes, ARSs, and noncoding RNAs
	2) 	Use BEDtools's "closest" function to find the closest upstream feature
		to each gene.
	3) 	Use BEDtools's getfasta to extract FASTA sequences of the intervals
		found in (2) and remove any duplicate sequences (not intervals, because
		we have many intervals that are on different strands) with fastx_tools.

DESCRIPTION OF FILES
--------------------

This folder contains five files:

	1) saccharomyces_cerevisiae.gff
	
		The S288C GFF file downloaded from SGD
		http://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff	
	
	2) saccharomyces_cerevisiae.filter.gff

		A GFF file containing only regions annotated as gene, ARS, or
		noncoding RNA.
		-------------------------------------------------------------
		for x in ARS gene ncRNA_gene rRNA_gene noncoding_exon; do 
			bioawk -v f=${x} -c gff '{if($feature == f) print $0}' saccharomyces_cerevisiae.gff >> saccharomyces_cerevisiae.filter.gff; 
		done

	3) saccharomyces_cerevisiae.gene.gff
		
		This GFF file contains only regions annotated as genes 
		(filtered to remove dubious orfs)
		------------------------------------------------------
		bioawk -c gff '{if($feature == "gene") print $0}'saccharomyces_cerevisiae.gff | grep -iv dubious > saccharomyces_cerevisiae.gene.gff

	4) saccharomyces_cerevisiae.upstream.bed

		This BED file contains upstream intergenic intervals. Columns are
		standard BED6 (as described here:
		http://bedtools.readthedocs.org/en/latest/content/general-usage.html),
		
		chrom, start, end, name, score, and strand

		chrom is in the format of chrI..chrXVI and chrmt,
		name is the "p" + the ORFID (e.g., pYBR294W),
		score is 0 (simply used as a placeholder).
		--------------------------------------------------
		Both GFF files were first sorted using BEDtools

		bedtools sort -i saccharomyces_cerevisiae.filter.gff > saccharomyces_cerevisiae.filter.sort.gff
		bedtools sort -i saccharomyces_cerevisiae.gene.gff > saccharomyces_cerevisiae.gene.sort.gff
		
		BEDtools then extracted the closest upstream feature for each gene

		bedtools closest -id -io -D a -a saccharomyces_cerevisiae.gene.sort.gff -b saccharomyces_cerevisiae.filter.sort.gff > saccharomyces_cerevisiae.closest.txt
		
		This creates a file that's two GFF's smashed together (think, R's
		merge), so we had to pull the relevant columns (and do a little
		arithmetic) to get the correct coordinates for a new BED file 
		containing all the upstream intervals. (I apologize for the
		ridiculousness that is the following bioawk one-liner)

		bioawk -t '{if($NF != -1) print $0}' saccharomyces_cerevisiae.closest.txt |
		bioawk -c gff '{d=substr($NF,2); split($9,a,";"); split(a[1],n,"="); 
			if($7=="+"){print $1, $start-d, $start-1, "p"n[2], 0, $7} else{print $1, $end, $end+d-1, "p"n[2], 0, $7} }' > saccharomyces_cerevisiae.upstream.bed

	5) 	saccharomyces_cerevisiae.upstream.fa

		This FASTA file contains the sequences of all the intervals found in
		saccharomyces_cerevisiae.upstream.bed
		--------------------------------------------------------------------
		Extract sequences using BEDtools's getfasta. 

		NB: This required a bit of editing to SGD's reference genome file. The
		FASTA headers of that file (S288C_reference_genome_R64-2-1_20150113.fsa)
		are long and detailed:

		>ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]

		In order for getfasta to work, these headers MUST be identical to the
		chrom columns of saccharomyces_cerevisiae.upstream.bed. I therefore
		changed all the FASTA headers of the genome to the "chrX" format I
		used. This is S288C.fa
		
		bedtools getfasta -s -name -fi S288C.fa -bed saccharomyces_cerevisiae.upstream.bed -fo saccharomyces_cerevisiae.upstream.fa
		
		NB^2: If you make these changes by hand (like I did), you will probably
		need to "repair" your FASTA, as there will be an unseen difference
		(probably whitespace somewhere) between your headers and chrom column. 
		I did this by reading and printing the FASTA with bioawk, which worked 
		just fine.

-------------------------------------------------------------------------------

Matt Rich (mattrich@uw.edu)
Stan Fields Lab
Department of Genome Sciences
University of Washington


go_slim_mapping.tab	This file is TAB delimited and contains the mapping of all yeast gene products (protein or RNA)
to a GO-Slim term.

Columns:	

1) ORF (mandatory) 		- Systematic name of the gene
2) Gene (optional) 		- Gene name, if one exists
3) SGDID (mandatory) 		- the SGDID, unique database identifier for the gene
4) GO_Aspect (mandatory) 	- which ontology: P=Process, F=Function, C=Component
5) GO Slim term (mandatory) 	- the name of the GO term that was selected as a GO Slim term
6) GOID (optional) 		- the unique numerical identifier of the GO term
7) Feature type (mandatory) 	- a description of the sequence feature, such as ORF or tRNA

A GO Slim is a subset of GO Terms that can be from the Biological
Process, Molecular Function, and Cellular Component ontologies.  These
may be general, high-level GO terms that represent major branches in
each ontology, as they are in go_slim_mapping.tab, or they may be more
granular terms that are used for a specific purpose (as in
the go_protein_complex.tab file below).

To determine the correct GO Slim term in the go_slim_mapping.tab file,
all GO annotations for a gene product are traced to a GO Slim term.

As of December 2007, please note that the go_slim_mapping.tab file, SGD's "GO Slim Mapper" tool
(http://www.yeastgenome.org/cgi-bin/GO/goSlimMapper.pl) and the "Genome
Snapshot" (http://www.yeastgenome.org/cache/genomeSnapshot.html)
handle parentage in the same way. Annotations are
mapped to all available GO Slim terms, regardless of parentage. For
example, if something is annotated to "meiosis" it will also be
annotated to "cell cycle". 

Each line contains the selected GO Slim term from the stated ontology
(P, F, or C) for that gene product.  Due to the structure of GO, the
GO annotations for a gene product may map to multiple GO Slim terms
for a single ontology.  Therefore, a gene product may be associated
with multiple GO Slim terms for a single ontology.

For those genes annotated to a non-GO Slim term, column 5 (GO_slim
term) displays 'Other' and in these cases the GOID column is blank.

Only annotations made by manually curated and high-throughput methods are included in this file.

This file is updated weekly.

For more information on the Gene Ontology (GO) project, see: 

  http://www.geneontology.org/
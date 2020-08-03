
FILES = CCLE_RNAseq_genes_counts_20180929.gct.gz \
	CCLE_metabolomics_20190502.csv \
	CCLE_RPPA_20181003.csv \
	CCLE_RPPA_Ab_info_20181226.csv \
	protein_quant_current_normalized.csv.gz \
	Table_S1_Sample_Information.xlsx \
	CCLE_mutations.csv

all: $(FILES)

clean: $(FILES)
	rm $(FILES)

CCLE_RNAseq_genes_counts_20180929.gct.gz: 
	wget https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_counts_20180929.gct.gz

CCLE_metabolomics_20190502.csv:
	wget https://data.broadinstitute.org/ccle/CCLE_metabolomics_20190502.csv

CCLE_RPPA_20181003.csv:
	wget https://data.broadinstitute.org/ccle/CCLE_RPPA_20181003.csv

CCLE_RPPA_Ab_info_20181226.csv:
	wget https://data.broadinstitute.org/ccle/CCLE_RPPA_Ab_info_20181226.csv

protein_quant_current_normalized.csv.gz:
	wget https://gygi.med.harvard.edu/sites/gygi.med.harvard.edu/files/documents/protein_quant_current_normalized.csv.gz

Table_S1_Sample_Information.xlsx:
	wget https://gygi.med.harvard.edu/sites/gygi.med.harvard.edu/files/documents/Table_S1_Sample_Information.xlsx

CCLE_mutations.csv:
	wget -O $@ https://ndownloader.figshare.com/files/22629110

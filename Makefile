

#--- Requirements
FILES = raw/CCLE_RNAseq_genes_counts_20180929.gct.gz \
	raw/CCLE_metabolomics_20190502.csv \
	raw/CCLE_RPPA_20181003.csv \
	raw/CCLE_RPPA_Ab_info_20181226.csv \
	raw/protein_quant_current_normalized.csv.gz \
	raw/Table_S1_Sample_Information.xlsx \
	raw/CCLE_mutations.csv \
	raw/CCLE_GCN.csv \
	raw/Achilles_GeneEssentiality.csv \
	raw/sample_info.csv


all: download MultiAssayExperiment

# shortcuts
download: $(FILES)

MultiAssayExperiment: data.md


#--- Download/Remove raw data
clean: $(FILES)
	rm $(FILES)

raw/CCLE_RNAseq_genes_counts_20180929.gct.gz: 
	wget -O $@ https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_counts_20180929.gct.gz

raw/CCLE_metabolomics_20190502.csv:
	wget -O $@ https://data.broadinstitute.org/ccle/CCLE_metabolomics_20190502.csv

raw/CCLE_RPPA_20181003.csv:
	wget -O $@ https://data.broadinstitute.org/ccle/CCLE_RPPA_20181003.csv

raw/CCLE_RPPA_Ab_info_20181226.csv:
	wget -O $@ https://data.broadinstitute.org/ccle/CCLE_RPPA_Ab_info_20181226.csv

raw/protein_quant_current_normalized.csv.gz:
	wget -O $@ https://gygi.med.harvard.edu/sites/gygi.med.harvard.edu/files/documents/protein_quant_current_normalized.csv.gz

raw/Table_S1_Sample_Information.xlsx:
	wget -O $@ https://gygi.med.harvard.edu/sites/gygi.med.harvard.edu/files/documents/Table_S1_Sample_Information.xlsx

raw/CCLE_mutations.csv:
	# DepMap Public 20Q2
	wget -O $@ https://ndownloader.figshare.com/files/22629110

raw/CCLE_GCN.csv:
	# DepMap Public 20Q2
	wget -O $@ https://ndownloader.figshare.com/files/22629107

raw/sample_info.csv:
	# DepMap Public 20Q2
	wget -O $@ https://ndownloader.figshare.com/files/22629137

raw/Achilles_GeneEssentiality.csv:
	# Depmap Public 20Q2
	wget -O $@ https://ndownloader.figshare.com/files/22629068

#--- Build MAE: processed data
data/MAE.rds: data.md

data.md: data.Rmd
	${CONDA_PREFIX}/bin/R -e "rmarkdown::render('$<', output_file='$@')"

#--- Conda environments
envs/mae:
	conda create -p $@ -f envs/mae.txt

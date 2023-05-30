1) TCGA_BRCA_gene_towork_new_withgenes_reduced.csv is the initial dataset with all the features for all the subtypes from which the results of Table I is generated.

2) In each subtype folder, the file such as TCGA_BRCA_gene_towork_new_withgenes_reduced_Basalagainstall.csv in Basal folder is the initial dataset where each subtype is considered in one-against-all fashion.

3) LASSOCV.py is run on TCGA_BRCA_gene_towork_new_withgenes_reduced_Basalagainstall.csv to find the value of alpha as well as the initial reduced set of features (TCGA_BRCA_gene_towork_new_withgenes_reduced_Basalagainstall_LASSO.csv).

4) mCGA.py is then executed on TCGA_BRCA_gene_towork_new_withgenes_reduced_Basalagainstall_LASSO.csv to find the final set of biomarkers. The result of this code is stored in TCGA_BRCA_gene_towork_new_withgenes_reduced_Basalagainstall_mCGA.

5) The final list of biomarkers are listed in Genes_SubtypeName.csv
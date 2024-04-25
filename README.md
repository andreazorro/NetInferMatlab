# NetInferMatlab
Pipeline for GRN inference from transcriptomic data in MATLAB

Tested on MATLAB R2023b with the Parallel Computing Toolbox. 

Methods 
--------------

* Anova: Implementation of the algorithm for inferring gene regulatory networks by ANOVA [1]
* Friedman: Variation of Anova applying a non-parametric Friedman Test [2] instead of a Two-way ANOVA [3].
* Statmodel [3]: Modification for GRN Inference of the method Statistical Modeling and Analysis of Experiments without ANOVA [4].
* Correlation: Pearson and Spearman correlation. 
* PCorr: Partial correlation. 

Required Input 
--------------

* One or multiple *.tsv files of the transcriptomic data in a matrix from (conditions,genes). First row must be the genes names. (See file Scoelicolor_rmabatch.tsv)
* One *.tfs file with a list of the transcription factor names, the same type of gene name as *.tsv files. (See file Scoelicolor_TF.tfs)

Output
------

Multiple files "<method>_<data>.tsv" with the inferred network interactions. 

* 1st Column: Transcription factor
* 2nd Column: Target gene
* 3rd Column: Inference confidence score   
		 
Bibliography		 
------------

[1] Küffner R, Petri T, Tavakkolkhah P, Windhager L, Zimmer R. Inferring gene regulatory networks by ANOVA. Bioinformatics. 2012;28:1376–82.

[2] Hoffman JIE. Chapter 26 - Analysis of Variance II. More Complex Forms. In: Hoffman JIE, editor. Biostatistics for Medical and Biomedical Practitioners. Academic Press; 2015. p. 421–47. doi:10.1016/B978-0-12-802387-7.00026-3.

[3] Zorro-Aranda, A., Escorcia-Rodríguez, J.M., González-Kise, J.K. et al. Curation, inference, and assessment of a globally reconstructed gene regulatory network for Streptomyces coelicolor. Sci Rep 12, 2840 (2022). https://doi.org/10.1038/s41598-022-06658-x

[4] Hernandez, H. Statistical Modeling and Analysis of Experiments without ANOVA. (2018)


License
-------

This project is licensed under the GNU General Public License. For the exact terms please see the [LICENSE file]

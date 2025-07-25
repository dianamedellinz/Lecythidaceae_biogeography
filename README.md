# Lecythidaceae_biogeography

Supplementary data for the manuscript: *A new biogeographic reconstruction for Lecythidaceae across bio and ecoregions of the American tropics supports the connectivity between the Amazon and Andean forests*

## Phylogenetic Analyses

We extracted and assembled coding sequences and introns for the nuclear markers targeted by the Angiosperms353 probe set[^1] using HybPiper v. 2.3.1[^2]. The Angiosperms353 dataset was compared with 343 nuclear loci previously identified as phylogenetically informative for Lecythidaceae[^3], identifying 31 overlapped genes. We merged the two nuclear datasets by combining their non-redundant loci, resulting in a final dataset of 665 nuclear loci. Exons and introns for each gene were separated and aligned using MAFFT v7.310[^4]. Poorly aligned regions were trimmed with trimAl[^5]. The alignments were concatenated using the _pxcat_ command[^6] to generate the supermatrix.

A Maximum likelihood (ML) phylogeny was inferred using IQ-TREE v.1.4.2[^7]. The phylogenetic inference was based on 665 nuclear markers across 202 species, representing the three subfamilies and 16 genera of Lecythidaceae[^8][^9]. To assess congruence among gene trees, we estimated unrooted gene trees in IQ-TREE, collapsing branches with bootstrap support <10%. A multi-species coalescent tree was inferred in ASTRAL[^10], and gene tree concordance and conflict across the species tree were visualized using PhyParts[^11].

## Ancestral Area Reconstruction

We time-calibrated the ML nuclear phylogeny using the penalized likelihood approach implemented in treePL[^12]. Four fossil constraints were used as calibration points: two from Asia and two from the Neotropics. The current geographic distribution of each species was compiled from different resources: the [NYBG/Lecythidaceae page](https://sweetgum.nybg.org/science/projects/lp/), [Plants of the World Online](https://powo.science.kew.org/), [JSTOR Global Plants](https://plants.jstor.org/), and [GBIF](https://www.gbif.org/). We prepared a presence/absence matrix for the species across seven ancestral areas following the bioregionalization proposed by Morrone et al. (2022)[^13].

We conducted biogeographic reconstructions for the Neotropical Lecythidaceae and the *Grias* + *Gustavia* clade in BioGeoBEARS[^14] using three models: DEC[^15], DIVALIKE (a likelihood interpretation of the DIVA model of Ronquist, 1997)[^16], and BAYAREALIKE[^17]. The Akaike Information Criterion (AIC) was used to select the best-fitting model and most reliable reconstructions.

### References

[^1]: Johnson, M. G., Pokorny, L., Dodsworth, S., Botigué, L. R., Cowan, R. S., Devault, A., … Wickett, N. J. (2018). A Universal Probe Set for Targeted Sequencing of 353 Nuclear Genes from Any Flowering Plant Designed Using k-medoids Clustering. *Systematic Biology Syst. Biol.* 68(4):594–606. https://doi.org/10.1093/sysbio/syy086
[^2]: Johnson, M.G., Gardner, E.M., Liu, Y., Medina, R., Goffinet, B., Shaw, A.J., Zerega, N.J., Wickett, N.J. (2016). HybPiper: Extracting coding sequence and introns for phylogenetics from high-throughput sequencing reads using target enrichment. *Appl Plant Sci.*, 4(7):apps.1600016. https://doi.org/10.3732/apps.1600016
[^3]: Vargas, O.M., Heuertz, M., Smith, S.A. and Dick. C.W. (2019). Target sequence capture in the Brazil nut family (Lecythidaceae): marker selection and in silico capture from genome skimming data. *Mol Phylogenet Evol*, 135:98–104. https://doi.org/10.1016/j.ympev.2019.02.020
[^4]: Katoh, K., and Standley, D. M. (2013). MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability. *Mol. Biol. Evol.*, 30(4), 772–780. https://doi.org/10.1093/molbev/mst010
[^5]: Capella-Gutiérrez, S., Silla-Martínez, J. M., Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics*, 25(15), 1972–1973. https://doi.org/10.1093/bioinformatics/btp348
[^6]: Brown, J. W., Walker, J. F., and Smith, S. A. (2017). Phyx: phylogenetic tools for unix. *Bioinformatics*, 33(12), 1886–1888. https://doi.org/10.1093/bioinformatics/btx063
[^7]: Nguyen, L.T., Schmidt, H.A., von Haeseler, A. and and Minh, B.Q. (2015). IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies. *Mol. Biol. Evol.*, 32, 268-274. https://doi.org/10.1093/molbev/msu300
[^8]: Huang, Y.Y., Mori, S.A., Kelly, L.M. (2015). Toward a phylogenetic-based generic classification of neotropical Lecythidaceae–I. Status of Bertholletia, Corythophora, Eschweilera, and Lecythis. *Phytotaxa* 203:85–121. https://doi.org/10.11646/phytotaxa.203.2.2
[^9]: Vargas, O. M., Larson, D. A., Batista, J., Cornejo, X., Garcia Luize, B., Medellín-Zabala, D., Ribeiro, M., Smith, N. P., Smith, S. A., Vincentini, A. and Dick, C. W. (2024). Reclassification of the Bertholletia clade of the Brazil nut family (Lecythidaceae) based on a phylogenetic analysis of plastome and target sequence capture data. *Harvard Papers in Botany*, 29 (1), 159–179. https://doi.org/10.3100/hpib.v29iss1.2024.n18
[^10]: Mirarab, S., Reaz, R., Bayzid, Md. S., Zimmermann, T., Swenson, M. S. and Warnow, T. (2014). ASTRAL: genome-scale coalescent-based species tree estimation. *Bioinformatics*, 30(17), i541–i548. https://doi.org/10.1093/bioinformatics/btu462
[^11]: Smith, S. A., Moore, M. J., Brown,   J. W. and Yang, Y. (2015). Analysis of phylogenomic datasets reveals conflict, concordance, and gene duplications with examples from animals and plants. *BMC Evolutionary Biology* 15: 150. https://doi.org/10.1186/s12862-015-0423-0
[^12]: Smith, S. A., and O’Meara, B. C. (2012). treePL: divergence time estimation using penalized likelihood for large phylogenies. *Bioinformatics*, 28(20), 2689–2690. https://doi.org/10.1093/bioinformatics/bts492
[^13]: Morrone, J.J., Escalante, T., Rodríguez-Tapia, G., Carmona, A., Arana, M. and Mercado-Gómez, J.D. (2022). Biogeographic regionalization of the Neotropical region: New map and shapefile. *An Acad Bras Cienc* (2022) 94(1): e20211167. https://doi.org/10.1590/0001-3765202220211167
[^14]: Matzke, N.J. (2013). BioGeoBEARS: Biogeography with Bayesian (and likelihood) evolutionary analysis in R scripts, CRAN: The Comprehensive R Archive Network, Vienna, Austria. http://cran.r-project.org/package.BioGeoBEARS
[^15]: Ree, R.H. and Smith, S.A. (2008). Maximum Likelihood Inference of Geographic Range Evolution by Dispersal, Local Extinction, and Cladogenesis. *Systematic Biology*. 57:4–14. https://doi.org/10.1080/10635150701883881
[^16]: Ronquist, F. (1997). Dispersal–vicariance analysis: a new approach to the quantification of historical biogeography. *Systematic Biology*, 46(1), 195–203. https://doi.org/10.1093/sysbio/46.1.195
[^17]: Landis, M.J., Matzke, N.J., Moore, B.R. and Huelsenbeck, J.P. (2013). Bayesian Analysis of Biogeography when the Number of Areas is Large. *Systematic Biology*. 62:789–804. https://doi.org/10.1093/sysbio/syt040

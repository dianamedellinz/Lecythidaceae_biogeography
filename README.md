# Lecythidaceae_biogeography

Supplementary data for the manuscript titled: "A new biogeographic reconstruction for Lecythidaceae across bio and ecoregions of the American tropics supports the connectivity between the Amazon and Andean forests".

## Phylogenetic Analyses

We extracted and assembled coding sequences and introns for the nuclear markers targeted by the Angiosperms353 probe set (Johnson et al., 2018) using HybPiper v. 2.3.1 (Johnson et al., 2016). The Angiosperms353 dataset was compared with 343 nuclear loci previously identified as phylogenetically informative for Lecythidaceae (Vargas et al., 2019), identifying 31 overlapped genes. We merged the two nuclear datasets by combining their non-redundant loci, resulting in a final dataset of 665 nuclear loci. Exons and introns for each gene were separated and aligned using MAFFT v7.310 (Katoh and Standley, 2013). Poorly aligned regions were trimmed with trimAl (Capella-Gutiérrez et al., 2009). The alignments were concatenated using the _pxcat_ command (Brown et al., 2017) to generate the supermatrix.

A Maximum likelihood (ML) phylogeny was inferred using IQ-TREE v.1.4.2 (Nguyen et al., 2015). The phylogenetic inference was based on 665 nuclear markers across 202 species, representing the three subfamilies (Huang et al. 2015) and the 16 genera described for Lecythidoideae in Vargas et al. (2024). To asses congruence among gene trees, we estimated unrooted gene trees in IQ-TREE, collapsing branches with bootstrap support <10%. A multi-species coalescent tree was inferred in ASTRAL (Mirarab et al., 2014), and gene tree concordance and conflict across the species tree were visualized using PhyParts (Smith et al., 2015).

## Ancestral Area Reconstruction

We time-calibrated the ML nuclear phylogeny using the penalized likelihood approach implemented in treePL (Smith and O’Meara, 2012). Four fossil constraints were used as calibration points: two from Asia and two from the Neotropics. The current geographic distribution of each species was compiled from different resources: the NYBG/Lecythidaceae page (https://sweetgum.nybg.org/science/projects/lp/), Plants of the World Online (https://powo.science.kew.org/), JSTOR Global Plants (https://plants.jstor.org/), and GBIF (https://www.gbif.org/). We prepared a presence/absence matrix for the species across seven ancestral areas following the bioregionalization proposed by Morrone et al. (2022).

We conducted biogeographic reconstructions for the Neotropical Lecythidaceae and the _Grias_ + _Gustavia_ clade in BioGeoBEARS (Matzke, 2013) using three models: DEC (Ree and Smith, 2008), DIVALIKE (a likelihood interpretation of the DIVA model of Ronquist, 1997), and BAYAREALIKE (Landis et al., 2013). The Akaike Information Criterion (AIC) was used to select the best-fitting model and most reliable reconstructions.

### References



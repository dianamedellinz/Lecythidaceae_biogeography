#######################################################
#######################################################

rm(list = ls())
dir()
getwd()

setwd("/Users/dianamedellin/Documents/Lecythidaceae/Biogeography/Biogeo_analyses_2025/BAYAREALIKE")

#######################################################
#######################################################
# From BioGeoBEARS Manual:
install.packages("Rcpp", dependencies=TRUE)
install.packages("RcppArmadillo", dependencies=TRUE)
install.packages("gdata", dependencies=TRUE)
install.packages("gtools", dependencies=TRUE)
install.packages("xtable", dependencies=TRUE)
install.packages("plotrix", dependencies=TRUE)
install.packages("vegan", dependencies=TRUE)
install.packages("FD", dependencies=TRUE)
install.packages("SparseM", dependencies=TRUE)
install.packages("ape", dependencies=TRUE)
install.packages("phylobase", dependencies=TRUE)
install.packages("rexpokit", dependencies=TRUE)
install.packages("cladoRcpp", dependencies=TRUE)
install.packages("optimx", dependencies=TRUE)
install.packages("snow", dependencies=TRUE)

# Additional ones
install.packages("GenSA", dependencies=TRUE)
install.packages("phytools", dependencies=TRUE)
install.packages("fdrtool", dependencies=TRUE)
install.packages("statmod", dependencies=TRUE)
install.packages("spam", dependencies=TRUE)
install.packages("MultinomialCI", dependencies=TRUE)

# From October 2018 onwards, install BioGeoBEARS from GitHub:
# https://github.com/nmatzke/BioGeoBEARS
install.packages("devtools")

library(devtools)
install.packages("usethis")
library(usethis)
devtools::install_github(repo="nmatzke/BioGeoBEARS", INSTALL_opts="--byte-compile", dependencies=FALSE)

#######################################################
#######################################################
library(BioGeoBEARS)

# Load other packages BioGeoBEARS depends on
library(ape)
library(cladoRcpp)

#-------------------
# Load phylogenetic tree
trfn =  "Lecy_ML_dated_NT.nwk"
moref(trfn)
tr = read.tree(trfn)
head(tr)

# Let's see how the tree looks like
plot(tr)
axisPhylo()

#-------------------
# Load and check geographic file in .txt format
geogfn = "Lecy_sample_areas.txt"

# This file should look like this: one column with the tip labels of your tree and
# a second column with the areas scored as presence (1) or absence (0).

# NOTE: Remember that the way in which you divide areas and score presences and absences
# depend on the geographic distribution of your group and which hypothesis/question you
# are interested in.

moref(geogfn)

# NOTE2: If it's the first time you're running this for your group, you can get a table of
# tip labels by running (remove hashtag to run):
# write.csv(tr$tip.label, file="tip_label_file.txt", row.names = F)

tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

#-------------------
# Check the maximum range size in your geog_data
max(rowSums(dfnums_to_numeric(tipranges@df)))
# In the examples with Myrtaceae, that is 1, meaning that no tip occurs in more than one area

# Set the maximum range size based on your geog_data
max_range_size = 7
# A value of 8 means that even though none of my modern taxa occur in more than one area,
# ancestral lineages can occupy up to four areas. That allows the model to find
# widespread ranges in the past, which can then be used to find evidence for viacariance 

# Check number of states which will influence the total analysis time (if more than 500-600 calculations will get really slow)
numstates_from_numareas(numareas=7, maxareas=7, include_null_range=TRUE)

#######################################################
#the BioGeoBEARS BayArea-like model is 
# not identical with the full Bayesian model implemented 
# in the "BayArea" program of Landis et al. (2013). 
# Instead, this is a simplified likelihood interpretation
# of the model.  Basically, in BayArea and BioGeoBEARS-BAYAREALIKE, 
# "d" and "e" work like they do in the DEC model of Lagrange 
# (and BioGeoBEARS), and then BayArea's cladogenesis assumption
# (which is that nothing in particular happens at cladogenesis) is 
# replicated by BioGeoBEARS.
######################################################

#######################################################
# Run BAYAREALIKE (preparing objects for BAYAREALIKE run)
######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    
BioGeoBEARS_run_object$include_null_range = TRUE    
BioGeoBEARS_run_object$on_NaN_error = -1e50   
BioGeoBEARS_run_object$speedup = TRUE         
BioGeoBEARS_run_object$use_optimx = "GenSA"   
BioGeoBEARS_run_object$num_cores_to_use = 4 # this might change
BioGeoBEARS_run_object$force_sparse = FALSE    
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)

# Check if run is ready -- should return "TRUE"
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#-------------------
resfn = "lecy_BAYAREALIKE_2025_2.Rdata" # Name file in which to save the run

#-------------------
# Run BAYAREALIKE
res = bears_optim_run(BioGeoBEARS_run_object)
save(res, file=resfn) # Save results
resBAYAREALIKE = res

LnL_3 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
# -546.8611

# The output of the run should include a lit of the arguments used, the estimated 
# parameters and inferred values of areas with highest probability in different parts of the
# tree. We will use these values to plot the likely ancestral ranges in the tree.

#-------------------
#Defining color by areas

# 1. Define area names
areanames <- c("P", "A", "G", "W", "E", "C", "M")

# 2. Define list of colors
mypalette <- c("red", "orange", "green", "turquoise", "blue", "#8000FF", "#FF00FF")

# Mixing colors for states
include_null_range = TRUE
states_list_0based_index = rcpp_areas_list_to_states_list(areas=areanames, maxareas=max_range_size, include_null_range=include_null_range)

statenames = areas_list_to_states_list_new(areas=areanames, maxareas=max_range_size, include_null_range=include_null_range, split_ABC=FALSE)
statenames

colors_matrix <- col2rgb(mypalette)
my_colors_list = mix_colors_for_states(colors_matrix, states_list_0based_index, plot_null_range=include_null_range)
my_colors_list

# Plot PDF of reconstruction
pdffn = "Lecythidaceae_BAYAREALIKE_2025.pdf"
pdf(pdffn, width=18, height=42)

# Plot ancestral states - DEC
analysis_titletxt ="BioGeoBEARS BAYAREALIKE on Lecythidaceae"

# Setup
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
plot_BioGeoBEARS_results(results_object, colors_list_for_states = my_colors_list, analysis_titletxt, plotwhat="text", label.offset=1.2, tipcex=1.2, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, colors_list_for_states = my_colors_list, analysis_titletxt, plotwhat="pie", label.offset=1.2, tipcex=1.2, statecex=0.5, splitcex=0.8, titlecex=0.8, plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  

# cmdstr = paste("open ", pdffn, sep="")
# system(cmdstr) 

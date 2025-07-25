#########################################################################
# CALCULATE SUMMARY STATISTICS TO COMPARE
# DEC, DIVALIKE, BAYAREALIKE
#########################################################################

setwd("/Users/dianamedellin/Documents/Lecythidaceae/Biogeography/Biogeo_analyses_2025")

# Set up empty tables to hold the statistical results
restable1 = NULL
teststable1 = NULL

restable2 = NULL
teststable2 = NULL

LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_3 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)

numparams1 = 2
numparams2 = 2
numparams3 = 2

stats1 = AICstats_2models(LnL_2, LnL_1,numparams2, numparams1)
stats1

stats2 = AICstats_2models(LnL_3, LnL_1,numparams1, numparams3)
stats2

res1 = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
res2 = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
res3 = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

#DEC is the null model. Comparing DIVALIKE with DEC
rbind(res1, res2, res3)

tmp_tests1 = conditional_format_table(stats1)

restable1 = rbind(restable1, res1, res2)
teststable1 = rbind(teststable1, tmp_tests1)

#Comparing BAYAREALIKE with DEC

tmp_tests2 = conditional_format_table(stats2)

restable2 = rbind(restable2, res1, res3)
teststable2 = rbind(teststable2, tmp_tests2)

teststable1$null = c("DEC")
teststable1$alt = c("DIVALIKE")
teststable1

teststable2$null = c("DEC")
teststable2$alt = c("BAYAREALIKE")
teststable2


# With AICs:
AICtable1 = calc_AIC_column(LnL_vals=restable1$LnL, nparam_vals=restable1$numparams)
restable1 = cbind(restable1, AICtable1)
restable_AIC_rellike1 = AkaikeWeights_on_summary_table(restable=restable1, colname_to_use="AIC")
restable_AIC_rellike1 = put_jcol_after_ecol(restable_AIC_rellike1)
restable_AIC_rellike1

AICtable2 = calc_AIC_column(LnL_vals=restable2$LnL, nparam_vals=restable2$numparams)
restable2 = cbind(restable2, AICtable2)
restable_AIC_rellike2 = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AIC")
restable_AIC_rellike2 = put_jcol_after_ecol(restable_AIC_rellike2)
restable_AIC_rellike2


######################################################
# ASSEMBLE RESULTS TABLES: DEC, DIVALIKE, BAYAREALIKE
######################################################

final <- rbind(restable_AIC_rellike1, restable_AIC_rellike2[2,])

rownames(final) <- c("DEC", "DIVALIKE", "BAYAREALIKE")
final

write.csv(final, file = "AIC_lecy.csv")
write.table(final, file = "AIC_lecy.txt", quote = F, sep = "\t")

# Also save to text files
write.table(restable_AIC_rellike1, file="restable_AIC_rellike1.txt", quote=FALSE, sep="\t")
write.table(restable_AIC_rellike2, file="restable_AIC_rellike2.txt", quote=FALSE, sep="\t")

# Save with nice conditional formatting
write.table(conditional_format_table(restable_AIC_rellike1), file="restable_AIC_rellike_formatted1.txt", quote=FALSE, sep="\t")
write.table(conditional_format_table(restable_AIC_rellike2), file="restable_AICc_rellike_formatted2.txt", quote=FALSE, sep="\t")


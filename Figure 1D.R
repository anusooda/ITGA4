
Gapmer_df = read.csv("sig_GapmervsUT.csv", header=TRUE)
Gapmer_control_df = read.csv("sig_GapmerControlvsUT.csv", header=TRUE)
MOE_df = read.csv("sig_MOEvsUT.csv", header=TRUE)
MOE_GTC_df = read.csv("sig_MOE_GTCvsUT.csv", header=TRUE)
PMO_GTC_df = read.csv("sig_PMO_GTCvsUT.csv", header=TRUE)


library(STRINGdb)
# create a new STRING_db object
string_db <- STRINGdb$new(version = "12", species = 9606, score_threshold = 200, input_directory="")

# map to STRING
Gapmer_mapped <- string_db$map(Gapmer_df, "X", removeUnmappedRows = TRUE )
Gapmer_control_mapped <- string_db$map(Gapmer_control_df, "X", removeUnmappedRows = TRUE )
MOE_mapped <- string_db$map(MOE_df, "X", removeUnmappedRows = TRUE )
MOE_GTC_mapped <- string_db$map(MOE_GTC_df, "X", removeUnmappedRows = TRUE )
PMO_GTC_mapped <- string_db$map(PMO_GTC_df, "X", removeUnmappedRows = TRUE )


# get the best 200 hits
hits1 <- Gapmer_mapped$STRING_id[1:200]
hits2 <- Gapmer_control_mapped$STRING_id[1:200]
hits3 <- MOE_mapped$STRING_id[1:200]
hits4 <- MOE_GTC_mapped$STRING_id[1:200]
hits5 <- PMO_GTC_mapped$STRING_id[1:200]

####################################### Figure 1D

# plot the STRING network png 
string_db$plot_network(hits1)
string_db$plot_network(hits3)

####################################### Supplementary Figure
string_db$plot_network(hits2)
string_db$plot_network(hits4)
string_db$plot_network(hits5)

#

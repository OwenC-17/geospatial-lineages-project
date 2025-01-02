library(ape)

sc2_tree <- read.tree("../input/cog_global_tree.newick")


"B.1.1.7" %in% sc2_tree$tip.label

cog_metadata <- read_csv("../input/cog_metadata.csv/cog_metadata.csv")
head(cog_metadata)

lineages_not_in_cog_meta <- long_ww_lin_w_sample_info %>%
  select(full_lineage_id, aliases_removed, named_variant_id) %>%
  distinct() %>%
  filter(!(full_lineage_id %in% cog_metadata$lineage))


str(cog_metadata)

sum(lineages_not_in_cog_meta$full_lineage_id %in% cog_metadata$scorpio_call)

str_detect(cog_metadata$scorpio_call, lineages_not_in_cog_meta$full_lineage_id)

big_string <- str_c(lineages_not_in_cog_meta$full_lineage_id, collapse = "|")

sum(str_detect(cog_metadata$scorpio_call, big_string), na.rm = TRUE)

str_extract(cog_metadata$lineage, "B\\.1\\.1\\.7")

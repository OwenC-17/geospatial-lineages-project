library(tidyverse)
library(ape)

sc2_tree <- read.tree("../input/cog_global_tree.newick")


"B.1.1.7" %in% sc2_tree$tip.label

cog_metadata <- read_csv("../input/cog_metadata.csv/cog_metadata.csv")
head(cog_metadata)

full_lineages_not_in_cog_meta <- idph_ww_lin_data %>%
  select(full_lineage_id, aliases_removed, named_variant_id) %>%
  distinct() %>%
  filter(!(full_lineage_id %in% cog_metadata$lineage | aliases_removed %in% cog_metadata$lineage))

t3_lineages_not_in_cog_meta <- idph_ww_lin_data %>%
  select(full_lineage_id, tier3id, tier2id, top_lin_id, aliases_removed, named_variant_id) %>%
  distinct() %>%
  filter(!(tier3id %in% cog_metadata$lineage | aliases_removed %in% cog_metadata$lineage | full_lineage_id %in% cog_metadata$lineage))

t2_lineages_not_in_cog_meta <- idph_ww_lin_data %>%
  select(full_lineage_id, tier3id, tier2id, top_lin_id, aliases_removed, named_variant_id) %>%
  distinct() %>%
  filter(!(tier2id %in% cog_metadata$lineage | aliases_removed %in% cog_metadata$lineage | full_lineage_id %in% cog_metadata$lineage))

any_level_lineages_not_in_cog_meta <- idph_ww_lin_data %>%
  select(full_lineage_id, tier4id, tier3id, tier2id, top_lin_id, aliases_removed, named_variant_id) %>%
  distinct() %>%
  filter(!(top_lin_id %in% cog_metadata$lineage | aliases_removed %in% cog_metadata$lineage | full_lineage_id %in% cog_metadata$lineage
           | tier2id %in% cog_metadata$lineage | tier3id %in% cog_metadata$lineage | tier4id %in% cog_metadata$lineage))

str(cog_metadata)


big_string <- str_c(full_lineages_not_in_cog_meta$full_lineage_id, collapse = "|")

sum(str_detect(cog_metadata$lineage, big_string), na.rm = TRUE)

b117s <- str_subset(cog_metadata$lineage, "B\\.1\\.1\\.7")
b117s <- data.frame(b117s) %>%
  filter(!is.na(b117s))
unique(b117s$b117s)


pango_designation_lineages <- read_csv("../input/lineages.csv")

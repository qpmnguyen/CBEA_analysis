library(HMP16SData)
library(curatedMetagenomicsData)

stool <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Stool" & VISITNO == 1) %>% as_phyloseq()
stool <- subset_samples(stool,!duplicated(RSID))



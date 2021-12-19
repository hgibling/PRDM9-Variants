## Step 2. Compile znf sequences, get list of unique znfs, and give them standardized names
Append author & year to each znf name and save in new file:
```
for FILE in intermediate-files/*znf-sequences*
do
NAME=$(basename ${FILE%-znf*})
awk -v NAME="$NAME" '{print NAME "_" $1 "\t" $2}' $FILE >> intermediate-files/publication-znf-sequences.tsv
done
```

Check which sequences are identical/unique:
```
# done in R

library(tidyr)
library(dplyr)

pub.znf <- read.table("publication-znf-sequences.tsv", header=F,
                  col.names=c("PubZnfName", "Sequence"))

# All unique sequences
znf.sequences.sperm <- pub.znf %>%
  # separate author from znf name
  separate(PubZnfName, sep="_", into=c("Publication", "ZnfName"), extra="merge") %>%
  mutate(Publication=sub("-", ".", Publication)) %>%
  # shift Oliver 2009 sequences to match the rest
  mutate(Sequence=ifelse(Publication=="oliver.2009", sub("TGCAGGGAG", "", Sequence), Sequence)) %>%
  mutate(Sequence=ifelse(Publication=="oliver.2009", sub("$", "TGCAGGGAG", Sequence), Sequence)) %>%
  # pivot wider to obtain unique znf sequences as rows
  pivot_wider(names_from=Publication, values_from=ZnfName) %>%
  # arrange in publication order (not including Oliver)
  arrange(berg.2010, berg.2011, borel.2012, jeffreys.2013, hussin.2013, alleva.2021)

# get list of znfs only observed in somatic blood/sperm
# a-z, symbols in jeffreys.2013
remove.germline <-znf.sequences.sperm %>%
  filter(!(jeffreys.2013 %in% LETTERS) & !is.na(jeffreys.2013)) %>%
  # arrange in jeffreys.2013 figure S2 order
  arrange(factor(jeffreys.2013, levels=c(letters, 1:9, 
                                c("!", "@", "£", "$", "%", "&", "§", "*", ":", "±")))) %>%
  mutate(StandardName=paste0("zn", str_pad(row_number(), 3, pad="0")), .before=1)

# All unique sequences observed in populations
znf.sequences <- znf.sequences.sperm %>%
  anti_join(remove.germline) %>%
  mutate(StandardName=paste0("ZN", str_pad(row_number(), 3, pad="0")), .before=1)

# Save full mapping info
write.table(remove.germline, "intermediate-files/standardized-znf-sequences-only-somatic-sperm-map.tsv", row.names=F, quote=F, sep="\t")
write.table(znf.sequences, "intermediate-files/standardized-znf-sequences-map.tsv", row.names=F, quote=F, sep="\t")
write.table(bind_rows(znf.sequences, remove.germline), "intermediate-files/standardized-znf-sequences-with-somatic-sperm-map.tsv", row.names=F, quote=F, sep="\t")

# Save just standard name and sequence
write.table(bind_rows(znf.sequences, remove.germline) %>% select(StandardName, Sequence), "intermediate-files/standardized-znf-sequences.tsv", row.names=F, quote=F, sep="\t", col.names=F)
```

---

## Step 3. Compile allele znf content, get list of unique alleles, and give them standardized names
Append author & year to each allele name and save in new file:
```
# exclude files with unnamed znfs
for FILE in intermediate-files/*allele-znf-content.tsv
do
NAME=$(basename ${FILE%-allele*})
awk -v NAME="$NAME" '{print NAME "\t" $1 "\t" $2}' $FILE >> intermediate-files/publication-allele-znf-content.tsv
done
```

Convert publication znf names to standardized names:
```
# done in R

library(tidyr)
library(dplyr)
library(stringr)

pub.znf.map <- read.table("intermediate-files/standardized-znf-sequences-with-somatic-sperm-map.tsv", header=T)

pub.znf.map.long <- pub.znf.map %>%
  select(-Sequence) %>%
  pivot_longer(cols=contains("20"), names_to="Publication", values_to="PubName",
               values_drop_na=T)


pub.allele <- read.table("intermediate-files/publication-allele-znf-content.tsv",
                             header=F, col.names=c("Publication", "AlleleName", "ZnfContent"))

# All unique sequences
allele.znf.sperm <- pub.allele %>%
  mutate(Publication=sub("-", ".", Publication)) %>%
  mutate(ZnfContent=case_when(
    Publication=="alleva.2021" ~ gsub("(.{2})", "\\1_\\2", ZnfContent),
    Publication!="alleva.2021" ~ gsub("(.)", "\\1_\\2", ZnfContent))) %>%
  mutate(ZnfContent=sub("_$", "", ZnfContent)) %>%
  separate_rows(ZnfContent, sep="_") %>%
  rename(PubName=ZnfContent) %>% 
  # ponting.2011 used lowercase version of berg names, so convert join with berg
  mutate(PubName=case_when(
    Publication=="ponting.2011" ~ toupper(PubName),
    TRUE ~ PubName)) %>%
  group_by(Publication, AlleleName) %>%
  # duplicate berg.2010 znfs to make a copy for ponting.2011
  # assume jeffreys.2013 N znf is the same as berg.2010 N znf and add
  left_join(bind_rows(pub.znf.map.long,
                      pub.znf.map.long %>% filter(Publication=="berg.2010") %>%
                        mutate(Publication="ponting.2011"),
                      pub.znf.map.long %>% filter(Publication=="berg.2010", PubName=="N") %>%
                        mutate(Publication="jeffreys.2013"))) %>%
  summarize(StandardZnfContent=str_c(StandardName, collapse="_")) %>%
  # remove ponting.2011 allele L24 since it is incorrect
  filter(!(Publication=="ponting.2011" & AlleleName=="L24")) %>%
  pivot_wider(names_from=Publication, values_from=AlleleName) %>%
  arrange(berg.2010, berg.2011, ponting.2011, borel.2012, jeffreys.2013, hussin.2013, alleva.2021)

# get list of alleles only observed in somatic blood/sperm
remove.germline.allele <- allele.znf.sperm %>%
  filter(across(c(alleva.2021, jeffreys.2013), ~ !is.na(.x))) %>%
  filter(across(c(-alleva.2021, -jeffreys.2013, -StandardZnfContent), is.na)) %>%
  mutate(StandardName=paste0("pr", str_pad(row_number(), 3, pad="0")), .before=1)

# All unique sequences observed in populations
allele.znf <- allele.znf.sperm %>%
  anti_join(remove.germline.allele) %>%
  mutate(StandardName=paste0("PR", str_pad(row_number(), 3, pad="0")), .before=1)

write.table(remove.germline.allele, "intermediate-files/standardized-allele-znf-content-only-somatic-sperm-map.tsv", row.names=F, quote=F, sep="\t")
write.table(allele.znf, "intermediate-files/standardized-allele-znf-content-map.tsv", row.names=F, quote=F, sep="\t")
write.table(bind_rows(allele.znf, remove.germline.allele), "intermediate-files/standardized-allele-znf-content-with-somatic-sperm-map.tsv", row.names=F, quote=F, sep="\t")

# Save just standardized name and znf content
write.table(bind_rows(allele.znf, remove.germline.allele) %>% select(StandardName, StandardZnfContent), "intermediate-files/standardized-allele-znf-content.tsv", row.names=F, quote=F, sep="\t", col.names=F)
```

---

## Step 4. Convert allele DNA sequences to standardized znf names to confirm known alleles and identify ones not in current list
Append author & year to each allele name and save in new file:
```
for FILE in intermediate-files/*allele-sequences.tsv
do
NAME=$(basename ${FILE%-allele*})
awk -v NAME="$NAME" '{print NAME "\t" $1 "\t" $2}' $FILE >> intermediate-files/publication-allele-sequences.tsv
done
```

Replace znf sequences with standardized znf names:
```
cp intermediate-files/publication-allele-sequences.tsv intermediate-files/publication-allele-sequences-standardized.tsv
while read ZNF SEQUENCE
do
sed -i "s/$SEQUENCE/$ZNF\_/g" intermediate-files/publication-allele-sequences-standardized.tsv
done < intermediate-files/standardized-znf-sequences.tsv
sed -i 's/\_$//' intermediate-files/publication-allele-sequences-standardized.tsv
```

See if sequences that still have unnamed nucleotide chunks can be separated into 84bp chunks
```
# done in R

library(tidyr)
library(dplyr)

allele.seqs.znf.converted <- read.table("intermediate-files/publication-allele-sequences-standardized.tsv", 
                                        header=F, col.names=c("Publication", "AlleleName", "ZnfContent"))

unknown <- allele.seqs.znf.converted %>%
  # get only seqs with nucleotides left
  filter(grepl("[ACGT]", ZnfContent)) %>%
  mutate(ZnfContent=gsub("([ACGT])(Z)", "\\1_\\2", ZnfContent, ignore.case=T)) %>%
  separate_rows(ZnfContent, sep="_") %>%
  # remove known znfs
  filter(!grepl("Z", ZnfContent)) %>%
  mutate(ZnfLength=nchar(ZnfContent)) %>%
  # try to split seqs at ending seq for all known znfs
  mutate(ZnfContent=gsub("CAGGGAG", "CAGGGAG_", ZnfContent)) %>%
  separate_rows(ZnfContent, sep="_") %>%
  mutate(ZnfLength=nchar(ZnfContent)) %>%
  filter(ZnfLength!=0)
```
Three sequences from Beyter 2021 have unidentified sequences that don't have the expected lengths of 84bp. 
Ignore alleles `chr5:23527530:FN.0`, `chr5:23527530:FN.1`, and `chr5:23527530:FN.4` for now.

Add any additional allele sequences to main list:
```
# done in R

library(tidyr)
library(dplyr)

allele.seqs.znf.converted <- read.table("intermediate-files/publication-allele-sequences-standardized.tsv", 
                                        header=F, col.names=c("Publication", "AlleleName", "StandardZnfContent"))
allele.seqs.znf.map <- read.table("intermediate-files/standardized-allele-znf-content-with-somatic-sperm-map.tsv", 
                                        header=T)


allele.seqs <- allele.seqs.znf.converted %>%
  # remove seqs with nucleotides left
  filter(!grepl("[ACGT]", StandardZnfContent)) %>%
  mutate(Publication=sub("-", ".", Publication)) %>%
  pivot_wider(names_from=Publication, values_from=AlleleName) %>%
  # join with original list of standardized allele znf content
  full_join(allele.seqs.znf.map) %>%
  relocate(StandardName, .before=StandardZnfContent) %>%
  relocate(beyter.2021, .after=berg.2011) %>%
  filter(!grepl("pr", StandardName)) %>%
  arrange(StandardName) %>%
  # add standardized name to new allele
  mutate(StandardName=ifelse(is.na(StandardName), paste0("PR", str_pad(row_number(), 3, pad="0")), StandardName)) %>%
  bind_rows(allele.seqs.znf.map %>% filter(grepl("pr", StandardName)))

# rewrite files
write.table(allele.seqs %>% filter(!grepl("pr", StandardName)), "intermediate-files/standardized-allele-znf-content-map.tsv", row.names=F, quote=F, sep="\t")
write.table(allele.seqs, "intermediate-files/standardized-allele-znf-content-with-somatic-sperm-map.tsv", row.names=F, quote=F, sep="\t")

write.table(allele.seqs %>% select(StandardName, StandardZnfContent), "intermediate-files/standardized-allele-znf-content.tsv", row.names=F, quote=F, sep="\t", col.names=F)
```
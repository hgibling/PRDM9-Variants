## Step 2. Compile znf sequences, get list of unique znfs, and give them temporary standardized names
Append author & year to each znf name and save in new file:
```
for FILE in intermediate-files/*znf-sequences*
do
NAME=$(basename ${FILE%-znf*})
awk -v NAME="$NAME" '{print NAME "_" $1 "\t" $2}' $FILE >> intermediate-files/publication-znf-sequences.tsv
done
```

Check which sequences are identical/unique and give standardized names:
```
# done in R

library(tidyr)
library(dplyr)
library(stringr)

pub.znf.seqs <- read.table("intermediate-files/publication-znf-sequences.tsv", header=F,
                      col.names=c("PubZnfName", "Sequence"))

# get all unique sequences and give standard names
znf.sequences <- pub.znf.seqs %>%
  # separate author from znf name
  separate(PubZnfName, sep="_", into=c("Publication", "ZnfName"), extra="merge") %>%
  mutate(Publication=sub("-", ".", Publication)) %>%
  mutate(Publication=str_to_title(Publication)) %>%
  # shift Oliver 2009 sequences to match the rest
  mutate(Sequence=ifelse(Publication=="Oliver.2009", sub("TGCAGGGAG", "", Sequence), Sequence)) %>%
  mutate(Sequence=ifelse(Publication=="Oliver.2009", sub("$", "TGCAGGGAG", Sequence), Sequence)) %>%
  # pivot wider to obtain unique znf sequences as rows
  pivot_wider(names_from=Publication, values_from=ZnfName) %>%
  # arrange in publication order (not including Oliver, since Berg is better known)
  arrange(Berg.2010, Berg.2011, Borel.2012, Jeffreys.2013, Hussin.2013, Alleva.2021) %>%
  # give temp standard names
  mutate(StandardName=paste0("z", str_pad(row_number(), 3, pad="0")), .before=1)

# Save full mapping info
write.table(znf.sequences, "intermediate-files/standardized-znf-sequences-map-step2.tsv", row.names=F, quote=F, sep="\t")

# Save just standard name and sequence
write.table(znf.sequences %>% select(StandardName, Sequence), "intermediate-files/standardized-znf-sequences-step2.tsv", row.names=F, quote=F, sep="\t", col.names=F)
```

---

## Step 3. Compile allele znf content, get list of unique alleles, and give them temporary standardized names
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

pub.znf.map <- read.table("intermediate-files/standardized-znf-sequences-map-step2.tsv", header=T)
pub.allele.znfs <- read.table("intermediate-files/publication-allele-znf-content.tsv",
                         header=F, col.names=c("Publication", "AlleleName", "ZnfContent"))

pub.znf.map.long <- pub.znf.map %>%
  select(-Sequence) %>%
  pivot_longer(cols=contains("20"), names_to="Publication", values_to="PubZnfName",
               values_drop_na=T)

# All unique sequences
allele.znfs <- pub.allele.znfs %>%
  mutate(Publication=sub("-", ".", Publication)) %>%
  mutate(Publication=str_to_title(Publication)) %>%
  # separate znfs with spacer
  mutate(ZnfContent=case_when(
    Publication=="Alleva.2021" ~ gsub("(.{2})", "\\1_\\2", ZnfContent),
    Publication!="Alleva.2021" ~ gsub("(.)", "\\1_\\2", ZnfContent))) %>%
  mutate(ZnfContent=sub("_$", "", ZnfContent)) %>%
  separate_rows(ZnfContent, sep="_") %>%
  rename(PubZnfName=ZnfContent) %>% 
  # Ponting.2011 used lowercase version of berg names, so convert to join with berg
  mutate(PubZnfName=case_when(
    Publication=="Ponting.2011" ~ toupper(PubZnfName),
    TRUE ~ PubZnfName)) %>%
  # duplicate Berg.2010 znfs to make a copy for Ponting.2011
  # assume Jeffreys.2013 N znf is the same as Berg.2010 N znf and add
  left_join(bind_rows(pub.znf.map.long,
                      pub.znf.map.long %>% filter(Publication=="Berg.2010") %>%
                        mutate(Publication="Ponting.2011"),
                      pub.znf.map.long %>% filter(Publication=="Berg.2010", PubZnfName=="N") %>%
                        mutate(Publication="Jeffreys.2013"))) %>%
    group_by(Publication, AlleleName) %>%
    summarize(StandardZnfContent=str_c(StandardName, collapse="_")) %>%
  # remove ponting.2011 allele L24 since it is incorrect
  filter(!(Publication=="Ponting.2011" & AlleleName=="L24")) %>%
  pivot_wider(names_from=Publication, values_from=AlleleName) %>%
  # sort in publication order
  arrange(Berg.2010, Berg.2011, Ponting.2011, Borel.2012, Jeffreys.2013, Hussin.2013, Alleva.2021) %>%
  # give temporary standard names
  mutate(StandardName=paste0("p", str_pad(row_number(), 3, pad="0")), .before=1)

# Save full mapping info
write.table(allele.znfs, "intermediate-files/standardized-allele-znf-content-map-step3.tsv", row.names=F, quote=F, sep="\t")

# Save just standardized name and znf content
write.table(allele.znfs %>% select(StandardName, StandardZnfContent), "intermediate-files/standardized-allele-znf-content-step3.tsv", row.names=F, quote=F, sep="\t", col.names=F)
```

---

## Step 4. Convert allele DNA sequences to temporary standardized znf names to confirm known alleles and identify ones not in current list
Append author & year to each allele name and save in new file:
```
for FILE in intermediate-files/*allele-sequences.tsv
do
NAME=$(basename ${FILE%-allele*})
awk -v NAME="$NAME" '{print NAME "\t" $1 "\t" $2}' $FILE >> intermediate-files/publication-allele-sequences.tsv
done
```

Replace znf sequences with temporary standardized znf names:
```
cp intermediate-files/publication-allele-sequences.tsv intermediate-files/publication-allele-sequences-standardized-step4.tsv

while read ZNF SEQUENCE
do
sed -i '' "s/$SEQUENCE/$ZNF\_/g" intermediate-files/publication-allele-sequences-standardized-step4.tsv
done < intermediate-files/standardized-znf-sequences-step2.tsv

# add _ between remaining nucs and named znfs & remove trailing _ after final named znf
sed -i '' -e 's/\([ACGT]\)\([Zz]\)/\1_\2/g' -e 's/\_$//' intermediate-files/publication-allele-sequences-standardized-step4.tsv
```

See if sequences that still have unnamed nucleotide chunks can be separated into 84bp chunks
```
# done in R

library(tidyr)
library(dplyr)

allele.seqs.znf.converted <- read.table("intermediate-files/publication-allele-sequences-standardized-step4.tsv", 
                                        header=F, col.names=c("Publication", "AlleleName", "ZnfContent"))
pub.znf.map <- read.table("intermediate-files/standardized-znf-sequences-map-step2.tsv", header=T)

unknown <- allele.seqs.znf.converted %>%
  # get only seqs with nucleotides left
  filter(grepl("[ACGT]", ZnfContent)) %>%
  separate_rows(ZnfContent, sep="_") %>%
  # remove known znfs
  filter(!grepl("Z", ZnfContent)) %>%
  mutate(ZnfLength=nchar(ZnfContent)) %>%
  # try to split seqs at common ending seq
  mutate(ZnfContent=ifelse(ZnfLength>84, 
                           gsub("CAGGGAG", "CAGGGAG_", ZnfContent), 
                           ZnfContent)) %>%
  separate_rows(ZnfContent, sep="_") %>%
  mutate(ZnfLength=nchar(ZnfContent)) %>%
  filter(ZnfLength!=0)
  # several nuc chunks do not have expected 84bp length

# one allele has 85bp nuc chunk instead of 84, could be 1bp insertion error
check.seq <- unknown %>%
  filter(AlleleName=="chr5:23527530:FN.4") %>%
  select(ZnfContent) %>%
  unname()

# calculate levenshtein distances between 85bp seq and known znf seqs
pub.znf.map %>%
  select(StandardName, Sequence) %>%
  mutate(LevDistanceToUnknown=adist(Sequence, check.seq)) %>%
  filter(LevDistanceToUnknown==min(LevDistanceToUnknown)) %>%
  as_tibble()

# two possible sequences are only off by 1
# StandardName Sequence                                  LevDistanceToUnknow…
# 1 z006         TGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGA…             1
# 2 z048         TGTGGGCGGGGCTTTAGAGATAAGTCACACCTCCTCAGACACCAGA…             1

# too ambiguous to declare, so ignore these three alleles
```
Three sequences from Beyter 2021 have unidentified sequences that don't have the expected lengths of 84bp. 
Ignore alleles `chr5:23527530:FN.0`, `chr5:23527530:FN.1`, and `chr5:23527530:FN.4` for now.

Add any additional allele sequences to main list:
```
# done in R

library(tidyr)
library(dplyr)
library(stringr)

allele.seqs.znf.converted <- read.table("intermediate-files/publication-allele-sequences-standardized-step4.tsv", 
                                        header=F, col.names=c("Publication", "AlleleName", "StandardZnfContent"))
allele.znf.map <- read.table("intermediate-files/standardized-allele-znf-content-map-step3.tsv", 
                                  header=T)


allele.znfs.updated <- allele.seqs.znf.converted %>%
  # remove seqs with nucleotides left
  filter(!grepl("[ACGT]", StandardZnfContent)) %>%
  mutate(Publication=sub("-", ".", Publication)) %>%
  mutate(Publication=str_to_title(Publication)) %>%
  pivot_wider(names_from=Publication, values_from=AlleleName) %>%
  # join with original list of standardized allele znf content
  full_join(allele.znf.map) %>%
  arrange(StandardName) %>%
  # add temp standardized name to new allele
  mutate(StandardName=ifelse(is.na(StandardName), 
                             paste0("P", str_pad(row_number(), 3, pad="0")), 
                             StandardName)) %>%
  relocate(StandardName, .before=StandardZnfContent) %>%
  relocate(Beyter.2021, .after=Berg.2011)

# write updated files
write.table(allele.znfs.updated, "intermediate-files/standardized-allele-znf-content-map-step4.tsv", row.names=F, quote=F, sep="\t")

write.table(allele.znfs.updated %>% select(StandardName, StandardZnfContent), "intermediate-files/standardized-allele-znf-content-step4.tsv", row.names=F, quote=F, sep="\t", col.names=F)
```
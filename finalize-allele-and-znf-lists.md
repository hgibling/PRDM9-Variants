# Finalize _PRDM9_ allele and znf lists
- Add marker if znf/allele has actually been observed in the human population or if it has only been observed in sperm cells or as somatic mutations in blood cells (i.e.; only observed in Jeffreys et al. 2013)
- Arrange alleles in reasonable order:
    - Alphabetical according to common names (Baudat and Berg A-I, L# names)
    - Then in publication order
    - Then in GenBank accession order
- Give permanent standardized names to both znfs and alleles
- Highlight alleles/znfs with inconsistent/repeted names
- Create table of all GenBank accessions, znf content, and standard allele names

---

## Step 10. Check if znfs and alleles have been observed in the human population, arrange and give final standardized names
```
# done in R

library(dplyr)
library(tidyr)
library(stringr)

all.znfs <- read.table("intermediate-files/standardized-znf-sequences-map-step7.tsv",
                       header=T)
all.alleles <- read.table("intermediate-files/standardized-allele-znf-content-map-step9.tsv",
                          header=T)
genbank.accessions <- read.table("genbank-records/standardized-allele-znf-content-step7.tsv",
                                 header=F, col.names=c("GenBankAccession", "StandardZnfContent"))

# get alleles observed in the population
alleles.with.pop.status <- all.alleles %>%
  mutate(InPopulation=case_when(
    !is.na(Baudat.2010) | 
      !is.na(Berg.2010) | 
      !is.na(Berg.2011) | 
      !is.na(Beyter.2021) | 
      !is.na(Borel.2012) | 
      !is.na(Hussin.2013) | 
      !is.na(Oliver.2009) | 
      !is.na(Parvanov.2010) | 
      !is.na(Ponting.2011) | 
      !is.na(GenBank.2021) ~ T,
    grepl("M", Alleva.2021) ~ T,
    TRUE ~ F), .before=2)

# get znfs from the alleles observed in the population
znfs.in.pop <- alleles.with.pop.status %>%
  filter(InPopulation==T) %>%
  select(StandardZnfContent) %>%
  separate_rows(StandardZnfContent, sep="_") %>%
  distinct() %>%
  arrange(StandardZnfContent) %>% 
  pull(StandardZnfContent)

# update znf list
znfs.with.pop.status <- all.znfs %>%
  mutate(InPopulation=case_when(
    StandardZnfName %in% znfs.in.pop ~ T,
    TRUE ~ F), .before=2) %>%
  # already in publication order, since only additions after initial list were GenBank
  # update to final standardized names
  mutate(StandardZnfName=sub("z", "Z", StandardZnfName))

# order for common allele names
baudat.berg.style.names <- c(LETTERS[c(1:6,8,9,11)], paste0("L", 1:nrow(alleles.with.pop.status %>% pivot_longer(cols=contains("20"), names_to="Publication", values_to="PubName", values_drop_na=T) %>% filter(grepl("L[0-9]{1,2}$", PubName)) %>% distinct(PubName))))

# get alleles with common names
alleles.common.names <- alleles.with.pop.status %>%
  pivot_longer(cols=contains("20"), names_to="Publication", values_to="PubName",
               values_drop_na=T) %>%
  filter(PubName %in% baudat.berg.style.names) %>%
  select(-Publication) %>%
  arrange(factor(PubName, levels=baudat.berg.style.names)) %>%
  filter(!duplicated(StandardName)) 

# publication order
pub.order <- c("Oliver.2009", "Parvanov.2010", "Baudat.2010", "Berg.2010", "Ponting.2011", "Berg.2011", "Borel.2012", "Jeffreys.2013", "Hussin.2013", "Beyter.2021", "Alleva.2021","GenBank.2021")

alleva.novel.order <- paste0("M", 1:length(alleles.with.pop.status$Alleva.2021[grepl("M", alleles.with.pop.status$Alleva.2021)]))

# get alleles with non-common names
alleles.not.common.names <- alleles.with.pop.status %>%
  pivot_longer(cols=contains("20"), names_to="Publication", values_to="PubName",
               values_drop_na=T) %>%
  filter(!(StandardName %in% alleles.common$StandardName)) %>%
  left_join(genbank.accessions) %>%
  mutate(PubName=ifelse(Publication=="GenBank.2021", GenBankAccession, PubName)) %>%
  select(-GenBankAccession) %>%
  arrange(factor(Publication, levels=not.common.pub.order), factor(PubName, levels=alleva.novel.order), PubName) %>%
  filter(!duplicated(StandardName))
  
allele.order <- c(alleles.common.names$StandardName, alleles.not.common.names$StandardName)  

# get final list
alleles.with.pop.status.final <- alleles.with.pop.status %>%
  arrange(factor(StandardName, levels=allele.order)) %>%
  # made up names for Jeffreys and Alleva describes them better, so remove
  select(-Jeffreys.2013) %>%
  mutate(StandardName=paste0("P", str_pad(row_number(), 3, pad="0")),
         StandardZnfContent=gsub("z", "Z", StandardZnfContent))

# save to file
write.table(alleles.with.pop.status.final, 
            "standardized-lists/population-and-sperm-somatic-PRDM9-standardized-allele-znf-content-map.tsv",
            row.names=F, quote=F, col.names=T, sep="\t")
write.table(alleles.with.pop.status.final %>% 
              select(StandardName, StandardZnfContent), 
            "standardized-lists/population-and-sperm-somatic-PRDM9-standardized-allele-znf-content.tsv",
            row.names=F, quote=F, col.names=F, sep="\t")
write.table(znfs.with.pop.status, "standardized-lists/population-and-sperm-somatic-PRDM9-standardized-znf-sequences-map.tsv",
            row.names=F, quote=F, col.names=T, sep="\t")
write.table(znfs.with.pop.status %>% 
              select(StandardZnfName, Sequence),
            "standardized-lists/population-and-sperm-somatic-PRDM9-standardized-znf-sequences.tsv",
            row.names=F, quote=F, col.names=F, sep="\t")

# just alleles/znfs found in pop
write.table(alleles.with.pop.status.final %>% 
              filter(InPopulation==T), 
            "standardized-lists/PRDM9-standardized-allele-znf-content-map.tsv",
            row.names=F, quote=F, col.names=T, sep="\t")
write.table(alleles.with.pop.status.final %>% 
              filter(InPopulation==T) %>% 
              select(StandardName, StandardZnfContent), 
            "standardized-lists/PRDM9-standardized-allele-znf-content.tsv",
            row.names=F, quote=F, col.names=F, sep="\t")
write.table(znfs.with.pop.status  %>% 
              filter(InPopulation==T), 
            "standardized-lists/PRDM9-standardized-znf-sequences-map.tsv",
            row.names=F, quote=F, col.names=T, sep="\t")
write.table(znfs.with.pop.status  %>% 
              filter(InPopulation==T) %>% 
              select(StandardZnfName, Sequence),
            "standardized-lists/PRDM9-standardized-znf-sequences.tsv",
            row.names=F, quote=F, col.names=F, sep="\t")
```

## Step 11. Regenerate allele sequences from znf sequences
Replace standard znf names with znf sequences to get full allele sequences for newly organized alleles
```
cp standardized-lists/PRDM9-standardized-allele-znf-content.tsv standardized-lists/PRDM9-standardized-allele-sequences.tsv

while read ZNF SEQUENCE
do
sed -i '' "s/$ZNF/$SEQUENCE/g" standardized-lists/PRDM9-standardized-allele-sequences.tsv
done < standardized-lists/PRDM9-standardized-znf-sequences.tsv

# check if any znfs left
grep "Z" standardized-lists/PRDM9-standardized-allele-sequences.tsv
# no

# remove _ between znfs
sed -i '' -e 's/_//g' standardized-lists/PRDM9-standardized-allele-sequences.tsv
```

Repeat for all alleles (including those only found in sperm/somatic blood)
```
cp standardized-lists/population-and-sperm-somatic-PRDM9-standardized-allele-znf-content.tsv standardized-lists/population-and-sperm-somatic-PRDM9-standardized-allele-sequences.tsv

while read ZNF SEQUENCE
do
sed -i '' "s/$ZNF/$SEQUENCE/g" standardized-lists/population-and-sperm-somatic-PRDM9-standardized-allele-sequences.tsv
done < standardized-lists/population-and-sperm-somatic-PRDM9-standardized-znf-sequences.tsv

# check if any znfs left
grep "Z" standardized-lists/population-and-sperm-somatic-PRDM9-standardized-allele-sequences.tsv
# no

# remove _ between znfs
sed -i '' -e 's/_//g' standardized-lists/population-and-sperm-somatic-PRDM9-standardized-allele-sequences.tsv
```
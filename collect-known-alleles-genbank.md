# Collecting known _PRDM9_ allele DNA sequences from GenBank
Many publications describing _PRDM9_ allele sequences submitted those sequences to GenBank. A broad GenBank searche for _PRDM9_ revealed several additional entries not in publications. Some of these entries were identical to published allleles while others were novel. Some alleles required extra clean-up to clarify the DNA sequences.

Outline of collection process:
- Search GenBank for PRDM9 nucleotide sequences
- Remove flanking sequences and keep only the znf repeat domain
- Replace DNA sequence with znf-temporary name if there is a match
- Examine remaining sequences and determine if they represent novel znfs
- Generate list of all GenBank entries and the allele they represent
- Determine unique list of alleles and give them allele-temporary names

## Step 5. Search GenBank for PRDM9 sequences
Search GenBank for _PRDM9_ nucleotide sequences
- Search in the [NCBI Nucleotide database](https://www.ncbi.nlm.nih.gov/nucleotide/):
    - `(PRDM9) AND "Homo sapiens"[porgn:__txid9606] NOT "patch" NOT "scaffold" NOT "assembly" NOT "chain" NOT "ZCWPW1"`
    - 128 entries (as of 2021-12-19)
- Save results:
    - _Send to: > Complete Record > Choose Destination: File > Format: FASTA > Sort by: Accession > Create File_ and save as `copy-paste-files/genbank-PRDM9-complete-record_2021-12-19.fa`
- Tidy file: `genbank-records/PRDM9-complete-record.fa`
```
# remove blank lines, collapse sequences to single lines
sed '/>/ s/$/NEWLINE/' copy-paste-files/genbank-PRDM9-complete-record_2021-12-19.fa | tr -d '\n' | sed 's/>/\n>/g' | sed 's/NEWLINE/\n/' > genbank-records/PRDM9-complete-record.fa

# turn into tsv
cat genbank-records/PRDM9-complete-record.fa | tr '\n' '\t' | sed 's/\t>/\n/g' | sed 's/ /\t/' > genbank-records/PRDM9-complete-record.tsv

# get list of accession numbers
cut -f1 genbank-records/PRDM9-complete-record.tsv | sort > genbank-records/PRDM9-accessions.txt

# append author & year to each accession and save in new file:
for FILE in genbank-records/*allele-sequence-accessions.txt
do
NAME=$(basename ${FILE%-allele*})
awk -v NAME="$NAME" '{print NAME "\t" $0}' $FILE >> genbank-records/publication-accessions.txt
done

# check that all publication accessions are present
wc -l genbank-records/publication-accessions.txt
grep -f <(cut -f2 genbank-records/publication-accessions.txt) genbank-records/PRDM9-accessions.txt | wc -l
# 95 for each--all present
# 33 additional records from genbank search
```

## Step 6. Replace genbank sequences with standardized znf names
```
cut -f1,3 genbank-records/PRDM9-complete-record.tsv > genbank-records/PRDM9-complete-record-standardized-step6.tsv 

while read ZNF SEQUENCE
do
sed -i '' "s/$SEQUENCE/$ZNF\_/g" genbank-records/PRDM9-complete-record-standardized-step6.tsv
done < intermediate-files/standardized-znf-sequences-step2.tsv

# remove leading and trailing sequences
# all alleles start with z001 (a)
# last part of exon 11 after the last znf is GATGAGTAA; some records only have GATGAG trailing, some have more
sed -i '' -e 's/\t[ACGT]*z001/\tz001/' -e 's/GATGAG\(TAA\)*.*$//' -e 's/\([ACGT]\)\([Zz]\)/\1_\2/g' -e 's/_$//' genbank-records/PRDM9-complete-record-standardized-step6.tsv

egrep "_[ACGT]" genbank-records/PRDM9-complete-record-standardized-step6.tsv | wc -l
# 16 entries have unknown zinc fingers/sequence chunks
```

### Look into sequences without any known znfs
```
# done in R

library(tidyr)
library(dplyr)

genbank.seqs <- read.table("genbank-records/PRDM9-complete-record-standardized-step6.tsv", 
                                        header=F, col.names=c("Accession", "ZnfContent"))

genbank.unknown <- genbank.seqs %>%
  filter(!grepl("z", ZnfContent)) %>%
  # check length, potential number of znfs (84bp)
  mutate(ZnfLength=nchar(ZnfContent)) %>%
  mutate(ZnfNum=ZnfLength/84) %>%
  # check for common znf motifs
  mutate(Start=ifelse(grepl("TGTGG", ZnfContent), T, F),
         End=ifelse(grepl("CAGGGAG", ZnfContent), T, F),
         ShortStart=ifelse(grepl("TGT", ZnfContent), T, F),
         ShortEnd=ifelse(grepl("GGAG", ZnfContent), T, F))
  
as_tibble(genbank.unknown)
#   Accession ZnfContent            ZnfLength ZnfNum Start End   ShortStart ShortEnd
#   <chr>     <chr>                     <int>  <dbl> <lgl> <lgl> <lgl>      <lgl>   
# 1 KT584285… AGGTCTGCAGGGAG               14  0.167 FALSE TRUE  FALSE      TRUE    
# 2 LP837495… ATTGTGAGATGTGTCAGAAC…       275  3.27  TRUE  FALSE TRUE       TRUE    
# 3 MW814869… TGGGAAGCTCCAAATCCTCT…       300  3.57  FALSE FALSE FALSE      FALSE   
# 4 MW814870… TGGGAAGCTCCAAATCCTCT…       300  3.57  FALSE FALSE FALSE      FALSE   
# 5 MW814871… GTAAGTGACACTTTTGGCCA…       208  2.48  FALSE FALSE TRUE       FALSE   
# 6 MW814872… ATTGTGAGATGTGTCAGAAC…       272  3.24  FALSE FALSE TRUE       TRUE    

genbank.unknown %>%
  # try splitting between start and stop motifs
  mutate(ZnfContent=gsub("GAGTGT", "GAG_TGT", ZnfContent)) %>%
  filter(grepl("_", ZnfContent))
  # no sequences
```
These are short fragments not long enough to be full znf regions:
- `KT584285.1` appears to be a short fragment at the end of a znf; GenBank listing says 3' UTR
- `MW814869.1`, `MW814870.1` are intron 1 (as per GenBank; znf region is in intron 11)
- `MW814871.1` is intron 7 (GenBank)
- `LP837495.1`, `MW814872.1` don't seem to be within the znf region at all
Ignore these sequences for further analyses

## Step 7. Look into sequences with some known znfs, some unknown sequences and update list of known znfs
```
# done in R

library(tidyr)
library(dplyr)
library(stringr)

genbank.seqs <- read.table("genbank-records/PRDM9-complete-record-standardized-step6.tsv", 
                           header=F, col.names=c("Accession", "ZnfContent"))
pub.accessions <- read.table("genbank-records/publication-accessions.txt", 
                             header=F, col.names=c("Publication", "Accession"))
known.znfs <- read.table("intermediate-files/standardized-znf-sequences-step2.tsv",
                         header=F, col.names=c("StandardZnfName", "Sequence"))

# look into alleles with some unknown znfs
genbank.partial <- genbank.seqs %>%
  # remove fully unknown sequences
  filter(grepl("z", ZnfContent)) %>%
  # remove fully known sequences
  filter(grepl("[ACGT]", ZnfContent)) %>%
  separate_rows(ZnfContent, sep="_") %>%
  # split znfs longer than 84bp
  mutate(ZnfContent=ifelse(nchar(ZnfContent)>84, 
                           gsub("CAGGGAGTGT", "CAGGGAG_TGT", ZnfContent), 
                           ZnfContent)) %>%
  separate_rows(ZnfContent, sep="_") %>%
  mutate(ZnfLength=ifelse(grepl("z", ZnfContent), 
                          84, nchar(ZnfContent))) %>%
  # give positional ID for each znf
  group_by(Accession) %>%
  mutate(ZnfPosition=row_number())

# look into non-84bp znfs
znfs.non.84bp <- genbank.partial %>%
  filter(ZnfLength!=84) %>%
  left_join(pub.accessions)
# all are from parvanov

# since parvanov had full amino acid seqs of the alleles in their publication, see if sequence can be deduced from those
parvanov.allele.aminos <- read.table("intermediate-files/parvanov-2010-allele-aminos.tsv", 
                                     header=F, col.names=c("Allele", "Amino"))

# compare to parvanov amino acid seqs
parvanov.znf.aminos <- parvanov.allele.aminos %>%
  # split into znf amino sequences
  mutate(Amino=gsub("(.{28})", "\\1_\\2", Amino)) %>%
  mutate(Amino=sub("_$", "", Amino)) %>%
  separate_rows(Amino, sep="_") %>%
  group_by(Allele) %>%
  mutate(ZnfPosition=row_number()) %>%
  mutate(ZnfLength=max(ZnfPosition))

# check if amino sequences match from parvanov publication sequences and genbank sequences
# copied from https://teaching.healthtech.dtu.dk/22110/index.php/Codon_list
aa.codons <- data.frame(
  stringsAsFactors = FALSE,
  Amino.Acid = c("Isoleucine","Leucine",
                 "Valine","Phenylalanine","Methionine","Cysteine",
                 "Alanine","Glycine","Proline","Threonine","Serine",
                 "Tyrosine","Tryptophan","Glutamine","Asparagine",
                 "Histidine","Glutamic acid","Aspartic acid","Lysine",
                 "Arginine","Stop codons"),
  SLC = c("I","L","V","F","M","C","A","G","P","T","S","Y",
          "W","Q","N","H","E","D","K","R","Stop"),
  DNA.codons = c("ATT, ATC, ATA","CTT, CTC, CTA, CTG, TTA, TTG",
                 "GTT, GTC, GTA, GTG","TTT, TTC","ATG","TGT, TGC",
                 "GCT, GCC, GCA, GCG","GGT, GGC, GGA, GGG",
                 "CCT, CCC, CCA, CCG","ACT, ACC, ACA, ACG",
                 "TCT, TCC, TCA, TCG, AGT, AGC","TAT, TAC","TGG",
                 "CAA, CAG","AAT, AAC","CAT, CAC","GAA, GAG",
                 "GAT, GAC","AAA, AAG","CGT, CGC, CGA, CGG, AGA, AGG",
                 "TAA, TAG, TGA")) %>%
  separate_rows(DNA.codons, sep=", ") %>%
  select(DNA.codons, SLC) %>%
  arrange(DNA.codons)

# convert known znf DNA sequences to amino sequences
known.znf.aminos <- known.znfs %>%
  # split into codons
  mutate(Sequence=gsub("(.{3})", "\\1_\\2", Sequence)) %>%
  mutate(Sequence=sub("_$", "", Sequence)) %>%
  separate_rows(Sequence, sep="_") %>%
  left_join(aa.codons, by=c("Sequence"="DNA.codons")) %>%
  select(-Sequence) %>%
  # merge amino sequence for each znf
  group_by(StandardZnfName) %>%
  summarize(Amino=paste0(SLC, collapse=""))

# convert non-84bp znfs to amino acid sequences
znfs.non.84bp.amino <- znfs.non.84bp %>%
  # split into codons
  mutate(ZnfContent=gsub("(.{3})", "\\1_\\2", ZnfContent)) %>%
  mutate(ZnfContent=sub("_$", "", ZnfContent)) %>%
  separate_rows(ZnfContent, sep="_") %>%
  left_join(aa.codons, by=c("ZnfContent"="DNA.codons")) %>%
  # if no AA, put lowercase nucs in place
  mutate(SLC=ifelse(is.na(SLC), tolower(ZnfContent), SLC)) %>%
  # merge amino sequence for each znf
  group_by(Accession, ZnfPosition) %>%
  summarize(SLC=paste0(SLC, collapse="")) %>%
  # move extra nucs to new column
  mutate(SLC=gsub("([A-Z])([acgt])", "\\1_\\2", SLC)) %>%
  separate(SLC, c("SLC","Extra"), sep="_")
# first 5 (56bp) are the same

# check if the 56bp fragment occurs in any parvanov alleles
seq.56 <- znfs.non.84bp.amino %>%
  filter(nchar(SLC)<=56/3) %>%
  ungroup() %>%
  select(SLC) %>%
  distinct() %>%
  as.data.frame() %>%
  unname()

parvanov.znf.aminos %>% filter(grepl(unlist(seq.56), Amino)) %>%
  mutate(Last=(ZnfPosition == ZnfLength))
# occurs in 15/16 alleles, always the last znf
# all the same sequence
# likely the genbank sequence was truncated in znf z010 (j)

# check which znfs the 56bp fragment occurs in 
known.znf.aminos %>% filter(grepl(unlist(seq.56), Amino))
# three possible znfs
# most likely z010  (j) since it usually occurs at the end of alleles and the others have only been observed in sperm/somatic samples

# check if the 85bp fragment occurs in any parvanov alleles
seq.85 <- znfs.non.84bp.amino %>%
  filter(nchar(SLC)>=56/3) %>%
  ungroup() %>%
  select(SLC) %>%
  distinct() %>%
  as.data.frame() %>%
  unname()

parvanov.znf.aminos %>% filter(grepl(seq.85, Amino))
# no
known.znf.aminos %>% filter(grepl(seq.85, Amino))
# no
# could be simple 1bp insertion that is throwing off the reading frame

# calculate levenshtein distances between 85bp seq and known znf seqs
seq.85.nuc <- znfs.non.84bp %>%
  filter(ZnfLength==85) %>%
  ungroup() %>%
  select(ZnfContent) %>%
  as.data.frame() %>%
  unname()

known.znfs %>%
  mutate(LevDistToUnknown=adist(Sequence, seq.85.nuc)) %>%
  filter(LevDistToUnknown==min(LevDistToUnknown))
# unknown is closest to z010 (j)

# summary: 56bp unknown znf most likely truncated z010 & 85bp znf most likely z010 with a 1bp insertion error

# update list of non-84bp znfs with sequence for z010
znfs.non.84bp.fixed <- znfs.non.84bp %>%
  mutate(ZnfContent="z010", ZnfLength=84) %>%
  select(-Publication)

genbank.partial.fixed <- genbank.partial %>%
  # remove non-84bp to replace with correct znf
  filter(ZnfLength==84) %>%
  bind_rows(znfs.non.84bp.fixed) %>%
  arrange(Accession, ZnfPosition)

# create temp standard names for novel znfs
known.znfs.updated <- known.znfs %>%
  bind_rows(genbank.partial.fixed %>%
              filter(grepl("[ACGT]", ZnfContent)) %>%
              ungroup() %>%
              select(ZnfContent) %>%
              rename(Sequence=ZnfContent)) %>%
  distinct() %>%
  mutate(StandardZnfName=ifelse(is.na(StandardZnfName),
                                paste0("z", str_pad(row_number(), 3, pad="0")),
                                StandardZnfName))

# update genbank allele znf content
genbank.seqs.updated <- genbank.partial.fixed %>%
  filter(grepl("[ACGT]", ZnfContent)) %>%
  rename(Sequence=ZnfContent) %>%
  left_join(known.znfs.updated) %>%
  select(-Sequence) %>%
  bind_rows(genbank.partial.fixed %>%
              filter(!grepl("[ACGT]", ZnfContent)) %>%
              rename(StandardZnfName=ZnfContent)) %>%
  arrange(Accession, ZnfPosition) %>%
  select(-ZnfLength) %>%
  # collapse individual znfs into znf content sequence
  summarize(ZnfContent=paste0(StandardZnfName, collapse="_")) %>%
  # add known genbank sequences
  bind_rows(genbank.seqs %>%
              filter(!grepl("[ACGT]", ZnfContent))) %>%
  arrange(Accession)

# save accession znf content to file
write.table(genbank.seqs.updated, 
            "genbank-records/standardized-allele-znf-content-step7.tsv",
            row.names=F, quote=F, sep="\t", col.names=F)

# update list of znfs
write.table(known.znfs.updated, 
            "intermediate-files/standardized-znf-sequences-step7.tsv",
            row.names=F, quote=F, sep="\t", col.names=F)
```

## Step 8. Compare accession sequences with publication sequences and update allele list
```
# done in R

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

genbank.seqs <- read.table("genbank-records/standardized-allele-znf-content-step7.tsv", 
                           header=F, col.names=c("Accession", "StandardZnfContent"))
pub.accessions <- read.table("genbank-records/publication-accessions.txt", 
                             header=F, col.names=c("Publication", "Accession"))
accession.znfs <- read.table("genbank-records/standardized-allele-znf-content-step7.tsv", 
                             header=F, col.names=c("Accession", "StandardZnfContent"))
# for some reason read.table() messes up the last ~23 lines, so use readr package
all.accessions <- read_tsv("genbank-records/PRDM9-complete-record.tsv",
                           col_names=c("Accession", "GenBankName", "Sequence"),
                           col_types = cols())

pub.allele.znfs.map <- read.table("intermediate-files/standardized-allele-znf-content-map-step4.tsv", header=T)
pub.znf.seqs <- read.table("intermediate-files/standardized-znf-sequences-step7.tsv",
                           header=F, col.names=c("StandardZnfName", "Sequence"))

accession.names <- accession.znfs %>%
  left_join(pub.accessions) %>%
  left_join(all.accessions %>%
              select(-Sequence)) %>%
  mutate(Publication=ifelse(!is.na(Publication),
                            str_to_title(Publication),
                            Publication)) %>%
  mutate(Publication=ifelse(!is.na(Publication),
                            sub("-", ".", Publication),
                            Publication))

# get list of unique accession sequences
unique.accession.seqs <- accession.names %>%
  select(-GenBankName) %>%
  arrange(Accession) %>%
  group_by(StandardZnfContent, Publication) %>%
  summarize(Accession=str_c(Accession, collapse="_")) %>%
  # three sequences repeated among genbank entries
  # keep only first accession for simplification
  mutate(Accession=ifelse(grepl("_", Accession), 
                          sub("_.*", "", Accession), Accession)) %>%
  left_join(accession.names) %>%
  # give names based on original fasta header
  mutate(Publication=ifelse(is.na(Publication), "GenBank.2021", Publication),
         GenBankAlleleName=case_when(
           grepl("PRDM9-", GenBankName) ~ sub("PRDM9-", "", 
                                              str_extract(GenBankName, "PRDM9-[A-Z]*[0-9]*")),
           # Berg.2010 misspelled PRDM9 in GenBank entries
           grepl("PRMD9-", GenBankName) ~ sub("PRMD9-", "", 
                                              str_extract(GenBankName, "PRMD9-[A-Z]*[0-9]*")),
           grepl("MW", Accession) ~ sub("_", ".", str_extract(GenBankName, "UP_[0-90]*")),
           grepl("isolate", GenBankName) ~ sub(" ", ".", str_extract(GenBankName, "isolate [0-9]*[a-z]*")),
           # only one seq that doesn't have an identifiable name from GenBank, so use accession
           TRUE ~ Accession
         )) %>%
  arrange(Publication) %>%
  select(-Accession, -GenBankName)

# check that for alleles described in both publications and GenBank, sequences match
genbank.authors <- unique(accession.names$Publication)
pub.authors <- colnames(pub.allele.znfs.map)[grepl("20", colnames(pub.allele.znfs.map))]
authors.pubs.and.genbank <- genbank.authors[genbank.authors %in% pub.authors]

pub.allele.znfs.map %>%
  pivot_longer(cols=contains("20"), names_to="Publication", values_to="PubAlleleName") %>%
  filter(!is.na(PubAlleleName)) %>%
  filter(Publication %in% authors.pubs.and.genbank) %>%
  full_join(unique.accession.seqs %>%
              filter(Publication %in% authors.pubs.and.genbank)) %>%
  filter(!(PubAlleleName==GenBankAlleleName) %in% T)

# StandardName   StandardZnfContent                              Publication PubAlleleName GenBankAlleleNa…
# <chr>          <chr>                                           <chr>       <chr>         <chr>
# 1 p001         z001_z002_z003_z004_z004_z005_z003_z006_z007_z… Berg.2010   A             NA              
# 2 p002         z001_z002_z003_z004_z004_z003_z003_z006_z007_z… Berg.2010   B             NA              
# 3 p003         z001_z002_z003_z004_z004_z003_z003_z006_z011_z… Berg.2010   C             NA              
# 4 p004         z001_z002_z003_z004_z004_z005_z003_z006_z011_z… Berg.2010   D             NA              
# 5 p005         z001_z002_z003_z004_z008_z006_z009_z010         Berg.2010   E             NA              
# 6 p582         z001_z002_z003_z004_z004_z005_z003_z071_z007_z… Hussin.2013 L35           NA              
# 7 p584         z001_z002_z003_z004_z004_z005_z003_z068_z007_z… Hussin.2013 L38           NA              
# 8 NA           z001_z002_z003_z004_z004_z005_z003_z068_z007_z… Hussin.2013 NA            L38             
# 9 NA           z001_z002_z003_z004_z004_z005_z003_z070_z007_z… Hussin.2013 NA            L35       

# Berg.2010 didn't put alleles A-E in GenBank, so those are ok
# Hussin.2013 L35 and L38 do not match publication sequences
# confirmed with Hussin that GenBank sequences are incorrect and publication sequences are correct

# add accession sequences to existing list of alleles
updated.allele.znf.content.map <- pub.allele.znfs.map %>%
  select(-StandardName) %>%
  pivot_longer(cols=contains("20"), names_to="Publication", values_to="PubAlleleName") %>%
  filter(!is.na(PubAlleleName)) %>%
  full_join(unique.accession.seqs %>%
              # remove incorrect Hussin.2013 alleles
              filter(!(Publication=="Hussin.2013" & GenBankAlleleName=="L35") & 
                       !(Publication=="Hussin.2013" & GenBankAlleleName=="L38")), 
            by=c("StandardZnfContent", "Publication", "PubAlleleName"="GenBankAlleleName")) %>%
  pivot_wider(names_from=Publication, values_from=PubAlleleName) %>%
  full_join(pub.allele.znfs.map %>%
              select(StandardName, StandardZnfContent)) %>%
  # rearrange columns
  select(StandardName, StandardZnfContent, 
         Alleva.2021, Baudat.2010, Berg.2010, Berg.2011, 
         Beyter.2021, Borel.2012, Hussin.2013, Jeffreys.2013, 
         Oliver.2009, Parvanov.2010, Ponting.2011, GenBank.2021) %>%
  # arrange in roughly publication order
  arrange(StandardName, Parvanov.2010, Baudat.2010, Berg.2010, 
          Ponting.2011, Berg.2011, Borel.2012, Jeffreys.2013, 
          Hussin.2013, Beyter.2021, Alleva.2021, GenBank.2021) %>%
  # name new alleles and assume all are observed in populations
  mutate(StandardName=ifelse(is.na(StandardName), 
                                   paste0("p", str_pad(row_number(), 3, pad="0")), 
                                   StandardName))

# update list of alleles
write.table(updated.allele.znf.content.map, 
            "intermediate-files/standardized-allele-znf-content-map-step8.tsv",
            row.names=F, quote=F, sep="\t")
```

## Step 9. Determine sequence of Baudat allele K given sequences for other alleles in GenBank
Baudat depicted allele K in their publication figure, but did not submit a sequence for it in GenBank.
```
# done in R
```
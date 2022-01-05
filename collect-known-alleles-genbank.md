# Collecting known _PRDM9_ allele DNA sequences from GenBank
Many publications describing _PRDM9_ allele sequences submitted those sequences to GenBank. A broad GenBank searche for _PRDM9_ revealed several additional entries not in publications. Some of these entries were identical to published allleles while others were novel. Some alleles required extra clean-up to clarify the DNA sequences.

Outline of collection process:
- Search GenBank for PRDM9 nucleotide sequences
- Remove flanking sequences and keep only the znf repeat domain
- Replace DNA sequence with znf-temporary name if there is a match
- Examine remaining sequences and determine if they represent novel znfs
- Generate list of all GenBank entries and the allele they represent
- Determine unique list of alleles and give them allele-temporary names

## Initial GenBank search
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
# 95 for each
# 33 additional records from genbank search
```

## Replace genbank sequences with standardized znf names
```
cut -f1,3 genbank-records/PRDM9-complete-record.tsv > genbank-records/PRDM9-complete-record-temp.tsv 

while read ZNF SEQUENCE
do
sed -i '' "s/$SEQUENCE/$ZNF\_/g" genbank-records/PRDM9-complete-record-temp.tsv
done < intermediate-files/standardized-znf-sequences.tsv

# remove leading and trailing sequences
# all alleles start with Z001
# last part of exon 11 after the last znf is GATGAGTAA; some records only have GATGAG trailing, some have more
sed -i '' -e 's/\t[ACGT]*Z001/\tZ001/' -e 's/GATGAG\(TAA\)*.*$//' -e 's/\([ACGT]\)\([Zz]\)/\1_\2/g' -e 's/_$//' genbank-records/PRDM9-complete-record-temp.tsv
```

## Look into sequences without known znfs
```
# done in R

library(tidyr)
library(dplyr)

genbank.seqs <- read.table("genbank-records/PRDM9-complete-record-temp.tsv", 
                                        header=F, col.names=c("Accession", "ZnfContent"))

genbank.unknown <- genbank.seqs %>%
  filter(!grepl("z", ZnfContent, ignore.case=T)) %>%
  # check length, potential number of znfs (84bp)
  mutate(ZnfLength=nchar(ZnfContent)) %>%
  mutate(ZnfNum=ZnfLength/84) %>%
  # check for common znf motifs
  mutate(Start=ifelse(grepl("TGTGG", ZnfContent), T, F),
         End=ifelse(grepl("CAGGGAG", ZnfContent), T, F),
         ShortStart=ifelse(grepl("TGT", ZnfContent), T, F),
         ShortEnd=ifelse(grepl("GGAG", ZnfContent), T, F))
  
as_tibble(unknown)
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
These are short fragments not long enough to be full znf regions
- `KT584285.1` appears to be a short fragment at the end of a znf; GenBank listing says 3' UTR
- `MW814869.1`, `MW814870.1` are intron 1 (GenBank)
- `MW814871.1` is intron 7
- `LP837495.1`, `MW814872.1` don't seem to be within the znf region at all
Ignore these sequences for further analyses

## Look into sequences with some known znfs, some unknown sequences
```
# done in R

# look alleles with some unknown znfs
genbank.partial <- genbank.seqs %>%
  # remove fully unknown sequences
  filter(grepl("z", ZnfContent, ignore.case=T)) %>%
  # remove fully known sequences
  filter(grepl("[ACGT]", ZnfContent)) %>%
  separate_rows(ZnfContent, sep="_") %>%
  # split znfs longer than 84bp
  mutate(ZnfContent=ifelse(nchar(ZnfContent)>84, 
                           gsub("CAGGGAGTGT", "CAGGGAG_TGT", ZnfContent), 
                           ZnfContent)) %>%
  separate_rows(ZnfContent, sep="_") %>%
  mutate(ZnfLength=ifelse(grepl("z", ZnfContent, ignore.case=T), 
                          84, nchar(ZnfContent))) %>%
  # give positional ID for each znf
  group_by(Accession) %>%
  mutate(ZnfPosition=row_number())

# look into non-84bp znfs
znfs.non.84bp <- genbank.partial %>%
  filter(ZnfLength!=84) %>%
  left_join(pub.accessions)
  # all are from parvanov

# compare to parvanov amino acid seqs
parvanov.allele.aminos <- read.table("intermediate-files/parvanov-2010-allele-aminos.tsv", 
                                  header=F, col.names=c("Allele", "Amino"))

parvanov.znf.aminos <- parvanov.allele.aminos %>%
  mutate(Amino=gsub("(.{28})", "\\1_\\2", Amino)) %>%
  mutate(Amino=sub("_$", "", Amino)) %>%
  separate_rows(Amino, sep="_") %>%
  group_by(Allele) %>%
  mutate(ZnfPosition=row_number()) %>%
  mutate(ZnfLength=max(ZnfPosition))

# check if amino sequences match
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
  mutate(Sequence=gsub("(.{3})", "\\1_\\2", Sequence)) %>%
  mutate(Sequence=sub("_$", "", Sequence)) %>%
  separate_rows(Sequence, sep="_") %>%
  left_join(aa.codons, by=c("Sequence"="DNA.codons")) %>%
  select(-Sequence) %>%
  group_by(StandardName) %>%
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
  mutate(SLC=gsub("([A-Z])([a-z])", "\\1_\\2", SLC)) %>%
  separate(SLC, c("SLC","Extra"), sep="_")
  # first 5 are the same

# 56bp
parvanov.znf.aminos %>% filter(grepl("CGRGFSDRSSLCYHQRTH", Amino)) %>%
  mutate(Last=(ZnfPosition == ZnfLength))
  # occurs in 15/16 alleles, always the last znf
  # all the same sequence
  # likely the genbank sequence was truncated in znf Z010 (j)
known.znf.aminos %>% filter(grepl("CGRGFSDRSSLCYHQRTH", Amino))
# three possible znfs
# most likely Z010  (j) since it occurs at the end of alleles and the others have only been observed in sperm/somatic samples

# 85bp
parvanov.znf.aminos %>% filter(grepl("CGRGFSDRSSLCYHQRTHTRGEALRLQG", Amino))
# no
known.znf.aminos %>% filter(grepl("CGRGFSDRSSLCYHQRTHTRGEALRLQG", Amino))
# no
# could be simple 1bp insertion that is throwing off the reading frame

# calculate levenshtein distances between 85bp seq and known znf seqs
known.znfs %>%
  mutate(LevDistToUnknown=adist(Sequence, "TGTGGGCGGGGCTTTAGCGATAGGTCAAGCCTCTGCTATCACCAGAGGACACACACAAGGGGAGAAGCCCTACGTCTGCAGGGAG")) %>%
  filter(LevDistToUnknown==min(LevDistToUnknown))
# unknown is closest to Z010 (j)

# summary: 56bp unknown znf most likely truncated Z010 & 85bp znf most likely Z010 with a 1bp insertion error


# update list of non-84bp znfs with sequence for Z010
Z010.seq <- known.znfs %>%
  filter(StandardName=="Z010") %>%
  select(Sequence)

znfs.non.84bp.fixed <- znfs.non.84bp %>%
  mutate(ZnfContent="Z010", ZnfLength=84) %>%
  select(-Publication)

genbank.partial.fixed <- genbank.partial %>%
  # remove non-84bp to replace
  filter(ZnfLength==84) %>%
  bind_rows(znfs.non.84bp.fixed) %>%
  arrange(Accession, ZnfPosition)

# create names for novel znfs
known.znfs.updated <- known.znfs %>%
  filter(grepl("Z", StandardName)) %>%
  bind_rows(genbank.partial.fixed %>%
  filter(grepl("[ACGT]", ZnfContent)) %>%
    ungroup() %>%
    select(ZnfContent) %>%
    rename(Sequence=ZnfContent)) %>%
  distinct() %>%
  mutate(StandardName=ifelse(is.na(StandardName),
                             paste0("Z", str_pad(row_number(), 3, pad="0")),
                             StandardName))

# update genbank allele znf content
genbank.seqs.updated <- genbank.partial.fixed %>%
  filter(grepl("[ACGT]", ZnfContent)) %>%
  rename(Sequence=ZnfContent) %>%
  left_join(known.znf.updated) %>%
  select(-Sequence) %>%
  bind_rows(genbank.partial.fixed %>%
              filter(!grepl("[ACGT]", ZnfContent)) %>%
              rename(StandardName=ZnfContent)) %>%
  arrange(Accession, ZnfPosition) %>%
  select(-ZnfLength) %>%
  # collapse individual znfs into znf content sequence
  summarize(ZnfContent=paste0(StandardName, collapse="_")) %>%
  # add known genbank sequences
  bind_rows(genbank.seqs %>%
              filter(!grepl("[ACGT]", ZnfContent))) %>%
  arrange(Accession)

# save to file
write.table(genbank.seqs.updated, 
            "genbank-records/PRDM9-allele-znf-content.tsv",
            row.names=F, quote=F, sep="\t")

write.table(bind_rows(known.znf.updated, 
                      known.znfs %>%
                        filter(grepl("z", StandardName))), 
            "intermediate-files/standardized-znf-sequences.tsv",
            row.names=F, quote=F, sep="\t", col.names=F)
```

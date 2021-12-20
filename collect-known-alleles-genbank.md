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
# all alleles start with ZN001
# last part of exon 11 after the last znf is GATGAGTAA; some records only have GATGAG trailing, some have more
sed -i '' -e 's/\t[ACGT]*ZN001/\tZN001/' -e 's/GATGAG\(TAA\)*.*$//' -e 's/\([ACGT]\)\([Zz]\)/\1_\2/' -e 's/_$//' genbank-records/PRDM9-complete-record-temp.tsv
```

## Look into sequences without known znfs
```
# done in R

library(tidyr)
library(dplyr)

genbank.seqs <- read.table("genbank-records/PRDM9-complete-record-temp.tsv", 
                                        header=F, col.names=c("Accession", "ZnfContent"))

genbank.seqs %>%
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

genbank.seqs %>%
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

## Look into sequences with some known znfs, some unknown sequences
```
# done in R

library(tidyr)
library(dplyr)

genbank.seqs <- read.table("genbank-records/PRDM9-complete-record-temp.tsv", 
                                        header=F, col.names=c("Acession", "ZnfContent"))

```

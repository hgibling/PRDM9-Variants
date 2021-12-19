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

# check that all publication accessions are present
cat genbank-records/*allele-sequence-accessions.txt > genbank-records/publication-accessions.txt
wc -l genbank-records/publication-accessions.txt
# 95
grep -f genbank-records/publication-accessions.txt genbank-records/PRDM9-accessions.txt | wc -l
# 95
# 33 additional records from genbank search
```

## Replace genbank sequences with standardized znf names


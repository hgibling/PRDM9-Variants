# Collecting known _PRDM9_ allele DNA sequences from GenBank
Most publications describing PRDM9 allele sequences had submitted those sequences to GenBank. A broad GenBank searche for _PRDM9_ revealed several additional entries not in publications. Some of these entries were repeats of published allleles while others were novel. Some alleles required extra clean up to clarify the DNA sequences.

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
    - 128 entries (as of 2021-11-24)
- Save results:
    - _Send to: > Complete Record > Choose Destination: File > Format: FASTA > Sort by: Accession > Create File_ and save as `genbank-records/genbank-PRDM9-complete-record_2021-11-24.fa`
- Tidy file
```
# remove blank lines, collapse sequences to single lines
sed '/>/ s/$/NEWLINE/' genbank-records/genbank-PRDM9-complete-record_2021-11-24.fa | tr -d '\n' | sed 's/>/\n>/g' | sed 's/NEWLINE/\n/' > genbank-records/genbank-PRDM9-complete-record-oneline_2021-11-24.fa

# remove non-znf sequences


```

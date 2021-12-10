# Collect known _PRDM9_ alleles and associated zinc fingers (znfs) from publications
Literature searches revealed several publications that describe znf DNA sequences, znf amino acid sequences, allele znf content, allele DNA sequences and/or their accession numbers, allele mutations, and/or allele structural variants (SVs):
- [Oliver 2009](#oliver-et-al-dec-2009)
  - Znf DNA sequences
  - Allele DNA sequence accession numbers
- [Thomas 2009](#thomas-et-al-dec-2009)
  - Znf DNA sequences
- [Parvanov 2010](#parvanov-et-al-feb-2010)
  - Allele DNA sequence accession numbers
- [Baudat 2010](#baudat-et-al-feb-2010)
  - Allele znf content
  - Allele DNA sequence accession numbers
- [Berg 2010](#berg-et-al-oct-2010)
  - Znf DNA sequences
  - Allele znf content
  - Allele DNA sequence accession numbers
- [Kong 2010](#kong-et-al-oct-2010)
  - Allele znf content
- [Ponting 2011](#ponting-may-2011)
  - Allele znf content
- [Berg 2011](#berg-et-al-jul-2011)
  - Znf DNA sequences
  - Allele znf content
- [Borel 2012](#borel-et-al-may-2012)
  - Znf DNA sequences
  - Allele znf content
- [Jeffreys 2013](#jeffreys-et-al-jan-2013)
  - Znf DNA sequences
  - Allele znf content
- [Hussin 2013](#hussin-et-al-mar-2013)
  - Znf DNA sequences
  - Allele znf content
  - Allele DNA sequence accession numbers
- [Beyter 2021](#beyter-et-al-may-2021)
  - Allele structural variants
    - Allele DNA sequences (inferred)
- [Wang 2021](#wang-et-al-jul-2021)
  - Allele mutations
    - Allele DNA sequences (inferred)
- [Alleva 2021](#alleva-et-al-nov-2021)
  - Znf DNA sequences
  - Znf amino acid sequences
  - Allele znf content
  - Allele DNA sequences

When possible, sequences were copy/pasted from publications and saved as `firstauthor-year-type.txt` in the `copy-paste-files` directory. After tidying up, content was saved as `firstauthor-year-type.tsv` in the `intermediate-files` directory. `Type` was one of the following:
- `znf-sequences`
- `znf-aminos`
- `allele-znf-content`
- `allele-sequence-accessions`
- `allele-sequences`
- `allele-mutations`

Genbank accession downloads are described in [additional documentation](/collect-known-alleles-genbank.md).

Analysis steps:
1. [Get allele and znf sequence data from publications](#step-1-get-allele-and-znf-sequence-data-from-publications)
2. [Check if allele content and znf sequences are unique](#step-2-check-if-allele-znf-content-and-znf-sequences-are-unique)
3. [Compile known znf sequences and allele znf content](#step-3-compile-known-znf-sequences-and-allele-znf-content)

---

## Step 1. Get allele and znf sequence data from publications

### Oliver et al. Dec 2009
### Accelerated Evolution of the _Prdm9_ Speciation Gene across Diverse Metazoan Taxa
**PMID: [19997497](https://pubmed.ncbi.nlm.nih.gov/19997497)**\
**GenBank Accession Numbers: [FJ899863.1 - FJ899912.1](https://ncbi.nlm.nih.gov/nuccore/?term=FJ899863.1%3AFJ899912.1%5Baccn%5D)**

Znf DNA sequences:
- Includes znfs `03`-`14`
- Copy/paste from **Supplementary Dataset S1** to: `copy-paste-files/oliver-2009-znf-copy.txt`
  - Sequences are shifted and begin with the last 9 nucleotides of previous zinc finger
- Tidy file: `intermediate-files/oliver-2009-znf-sequences.tsv`
```
# extract human znfs
sed '/>/ s/$/NEWLINE/' copy-paste-files/oliver-2009-znf-copy.txt | tr -d '\n' | sed 's/>/\n/g' | grep "homo_sapien" | sed 's/NEWLINE/\t/' > intermediate-files/oliver-2009-znf-sequences.tsv

# check if all unique
wc -l intermediate-files/oliver-2009-znf-sequences.tsv
# 12
cut -f2 intermediate-files/oliver-2009-znf-sequences.tsv | sort | uniq | wc -l
# 9

# remove duplicate sequences, keeping lowest ID number
sort -u -k2,2 intermediate-files/oliver-2009-znf-sequences.tsv | sort > intermediate-files/TEMP-oliver-2009-znf-sequences.tsv
rm intermediate-files/oliver-2009-znf-sequences.tsv
mv intermediate-files/TEMP-oliver-2009-znf-sequences.tsv intermediate-files/oliver-2009-znf-sequences.tsv
wc -l intermediate-files/oliver-2009-znf-sequences.tsv
# 9
```

Allele DNA sequence accession numbers:
- Save accession numbers to: `genbank-records/oliver-2009-allele-sequence-accessions.txt`
```
# generate sequence of numbers representing genbank accession numbers
for i in $(seq 899863 899912)
do
echo "FJ$i.1" >> genbank-records/oliver-2009-allele-sequence-accessions.txt
done
```

#

### Thomas et al. Dec 2009
### Extraordinary Molecular Evolution in the _PRDM9_ Fertility Gene
**PMID: [20041164](https://pubmed.ncbi.nlm.nih.gov/20041164)**\
**GenBank Accession Numbers: None**

Znf DNA sequences:
- Includes znfs `ZF1`-`ZF12`
- Znf sequences depicted in **Figure 4**
- Not copy/pastable

<img src="https://github.com/hgibling/PRDM9-Variants/blob/main/images/Thomas-2009-Figure4.png?raw=true" width="600">

#

### Parvanov et al. Feb 2010
### _PRDM9_ controls activation of mammalian recombination hotspots
**PMID: [20044538](https://pubmed.ncbi.nlm.nih.gov/20044538)**\
**GenBank Accession Numbers: [GU183914.1 - GU183919.1](https://ncbi.nlm.nih.gov/nuccore/?term=GU216222.1%3AGU216229.1%5Baccn%5D)**

Allele DNA sequence accession numbers:
- Accession numbers listed on **Supplementary Material page 3: Material and Methods, Gene Bank Numbers** include from GU183909-GU183919, but only GU183914-GU183919 are human (the others are mouse)
- Save accession numbers to: `genbank-records/parvanov-2010-allele-sequence-accessions.txt`
```
# generate sequence of numbers representing genbank accession numbers
for i in $(seq 183914 183919)
do
echo "GU$i.1" >> genbank-records/parvanov-2010-allele-sequence-accessions.txt
done
```

#

### Baudat et al. Feb 2010
### _PRDM9_ is a major determinant of meiotic recombination hotspots in humans and mice
**PMID: [20044539](https://pubmed.ncbi.nlm.nih.gov/20044539)**\
**GenBank Accession Numbers: [GU216222.1 - GU216229.1](https://ncbi.nlm.nih.gov/nuccore/?term=GU216222.1%3AGU216229.1%5Baccn%5D)**

Allele znf content:
- Includes alleles `A`-`E`, `F`, `H`-`I`, `K`
- Image in **Figure 2b** depicts znf content as unnamed blocks
  - Appears one block (`-NHR` in the legend) represents more than one znf as it appears in positions 11 & 12 in alleles A and B
  - Additionally, tne block in alleles F and K (`--HR`) and one block in allele H (`-RVS`)that are not also present in A-E
- Type znf content typed out by hand, triple check for accuracy, to: `copy-paste-files/baudat-2010-allele-copy.txt`
  - Name blocks by four-character code in figure legend (orange = `NTOR`; grey = `NTGR`), with `|` to separate blocks
- Tidy file: `intermediate-files/baudat-2010-allele-znf-content-unnamed.tsv`

<img src="https://github.com/hgibling/PRDM9-Variants/blob/main/images/Baudat-2010-Figure2B.png?raw=true" width="600">

```
# tidy file
awk '{print $1 "\t" $3}' copy-paste-files/baudat-2010-allele-copy.txt > intermediate-files/baudat-2010-allele-znf-content-unnamed.tsv
```

Allele DNA sequence accession numbers:
- Includes alleles `A`-`F`, `H`-`I`
- Sequences in **Supplementary Figure S3B** are _screenshots_ instead of copy/pastable text (!)
- Save accession numbers to: `genbank-records/baudat-2010-allele-sequence-accessions.txt`
```
# generate sequence of numbers representing genbank accession numbers
for i in $(seq 216222 216229)
do
echo "GU$i.1" >> genbank-records/baudat-2010-allele-sequence-accessions.txt
done
```

#

### Berg et al. Oct 2010
### _PRDM9_ variation strongly influences recombination hot-spot activity and meiotic instability in humans
**PMID: [20818382](https://pubmed.ncbi.nlm.nih.gov/20818382)**\
**GenBank Accession Numbers: [HM210983.1 – HM211006.1](https://ncbi.nlm.nih.gov/nuccore/?term=HM210983.1%3AHM211006.1%5Baccn%5D)**

Znf DNA sequences:
- Includes znfs `a`-`t`
- Copy/paste from **Supplementary Figure 1a** to: `copy-paste-files/berg-2010-znf-copy.txt`
- Tidy file: `intermediate-files/berg-2010-znf-sequences.tsv`
```
# tidy file
sed -e 's/\s/\t/' -e '$a\' copy-paste-files/berg-2010-znf-copy.txt > intermediate-files/berg-2010-znf-sequences.tsv

# check if all unique
wc -l intermediate-files/berg-2010-znf-sequences.tsv
# 20
cut -f2 intermediate-files/berg-2010-znf-sequences.tsv | sort | uniq | wc -l
# 20
```

Allele znf content:
- Includes alleles `A`-`E`, `L1`-`L24`
- Copy/paste from **Supplementary Figure 1b** to: `copy-paste-files/berg-2010-allele-copy.txt`
- Tidy file: `intermediate-files/berg-2010-allele-znf-content.tsv`

```
# remove extra columns, tidy file, and sort alphabetically
awk '{print $1 "\t" $3}' copy-paste-files/berg-2010-allele-copy.txt | sort -k1,1V > intermediate-files/berg-2010-allele-znf-content.tsv
```

Allele DNA sequence accession numbers:
- Save accession numbers to: `genbank-records/berg-2010-allele-sequence-accessions.txt`
```
# generate sequence of numbers representing genbank accession numbers
for i in $(seq 210983 211006)
do
echo "HM$i.1" >> genbank-records/berg-2010-allele-sequence-accessions.txt
done
```

#

### Kong et al. Oct 2010
### Fine-scale recombination rate differences between sexes, populations and individuals
**PMID: [20981099](https://pubmed.ncbi.nlm.nih.gov/20981099)**\
**GenBank Accession Numbers: None**

Allele znf content:
- Includes alleles `Kong01`-`Kong12` (publication did not provide allele names, so name them here)
- Image in **Supplementary Figure 4** depicts znf content as blocks named with the amino acids at repeat positions -1, 3 and 6 of the alpha helix
- Type znf content typed out by hand, triple check for accuracy, to: `copy-paste-files/kong-2010-allele-copy.txt`
  - Add `-` to separate blocks
- Tidy file: `intermediate-files/kong-2010-allele-znf-content-unnamed.tsv`

<img src="https://github.com/hgibling/PRDM9-Variants/blob/main/images/Kong-2010-SupFigure4.png?raw=true" width="600">

```
# remove gaps, tidy file, give new allele names
awk '{print $2}' copy-paste-files/kong-2010-allele-copy.txt | sed -e 's/-\{2,\}/-/g' | grep "-" | sort | uniq | awk '{printf "Decode%02i\t%s\n", NR, $1}' > intermediate-files/kong-2010-allele-znf-content-unnamed.tsv
```

#

### Ponting May 2011
### What are the genomic drivers of the rapid evolution of _PRDM9_?
**PMID: [21388701](https://pubmed.ncbi.nlm.nih.gov/21388701)**\
**GenBank Accession Numbers: None**

Allele znf content:
- Includes alleles `A`-`L24`
- Review paper
- Image in **Figure 4** depicts znf content as named blocks
  - Cites Berg et al. 2010 for allele znf content
- Type Znf content out by hand, triple check for accuracy, to: `copy-paste-files/ponting-2010-allele-znf-content.txt`
  - Include gaps represented with `_`
- Tidy file: `intermediate-files/ponting-2011-allele-znf-content.tsv`
```
# remove gaps and sort alphabetically
sed -e 's/\s/\t/g' -e 's/_//g' copy-paste-files/ponting-2011-allele-copy.txt | cut -f1,2 | sort -k1,1V > intermediate-files/ponting-2011-allele-znf-content.tsv
```

<img src="https://github.com/hgibling/PRDM9-Variants/blob/main/images/Ponting-2011-Figure3.jpg?raw=true" width="500">

#

### Berg et al. Jul 2011
### Variants of the protein _PRDM9_ differentially regulate a set of human meiotic recombination hotspots highly active in African populations
**PMID: [21750151](https://pubmed.ncbi.nlm.nih.gov/21750151)**\
**GenBank Accession Numbers: None**

Znf DNA sequences
- Includes znfs `a`-`l`, `o`-`v`
- Copy/paste from **Supplementary Figure 1A** to: `copy-paste-files/berg-2011-znf-copy.txt`
- Tidy file: `intermediate-files/berg-2011-znf-sequences.tsv`
```
# tidy file
sed -e 's/\s/\t/'  -e '$a\' copy-paste-files/berg-2011-znf-copy.txt > intermediate-files/berg-2011-znf-sequences.tsv

# check if all unique
wc -l intermediate-files/berg-2011-znf-sequences.tsv
# 20
cut -f2 intermediate-files/berg-2011-znf-sequences.tsv | sort | uniq | wc -l
# 20
```

Allele znf content:
- Includes alleles `A`-`E`, `L1`-`L27`
- Copy/paste from **Supplementary Figure 1B** to: `copy-paste-files/berg-2011-allele-copy.txt`
- Tidy file: `intermediate-files/berg-2010-allele-znf-content.tsv`

```
# remove extra columns, tidy file, and sort alphabetically
awk '{print $1 "\t" $3}' copy-paste-files/berg-2011-allele-copy.txt | sort -k1,1V > intermediate-files/berg-2011-allele-znf-content.tsv
```

#

### Borel et al. May 2012
### Evaluation of _PRDM9_ variation as a risk factor for recurrent genomic disorders and chromosomal non-disjunction
**PMID: [22643917](https://pubmed.ncbi.nlm.nih.gov/22643917)**\
**GenBank Accession Numbers: None**

Znf DNA sequences:
- Includes znfs `a`-`m`, `q`
- Copy/paste from **Supplementary Table S1** to: `copy-paste-files/borel-2012-znf-copy.txt`
- Tidy file: `intermediate-files/borel-2012-znf-sequences.tsv`
```
# tidy file
sed -e 's/\s/\t/' copy-paste-files/borel-2012-znf-copy.txt > intermediate-files/borel-2012-znf-sequences.tsv

# check if all unique
wc -l intermediate-files/borel-2012-znf-sequences.tsv
# 14
cut -f2 intermediate-files/borel-2012-znf-sequences.tsv | sort | uniq | wc -l
# 14
```

Allele znf content:
- Includes alleles `A`-`F`, `I`, `L1`, `L19`, `L28`-`L31`
- Copy/paste from **Supplementary Table S1** to: `copy-paste-files/borel-2012-allele-copy.txt`
- Tidy file: `intermediate-files/borel-2012-allele-znf-content.tsv`
```
# remove extra column
awk '{print $1 "\t" $3}' copy-paste-files/borel-2012-allele-copy.txt > intermediate-files/borel-2012-allele-znf-content.tsv
```

#

### Jeffreys et al. Jan 2013
### Recombination regulator _PRDM9_ influences the instability of its own coding sequence in humans
**PMID: [23267059](https://pubmed.ncbi.nlm.nih.gov/23267059)**\
**GenBank Accession Numbers: None**

Znf DNA sequences:
- Study looks at low-frequency mutations in blood and genotypes of sperm cells; these znfs not necessarily observed as part of a human genotype
- Includes znfs `A`-`L`, `O`-`V`, `a`-`z`, `1-9`, `!`, `@`, `£`, `$`, `%`, `&`, `§`, `*`, `:`, `±`
- Copy/paste from **Supplementary Figure S2** to: `copy-paste-files/jeffreys-2013-znf-copy.txt`
  - Icons `£`, `§`, and `±` may not render properly when copy/pasted (appear as `_`) depending on language settings; these icons are used to name sequences on lines 57, 61, and 64 respectively
- Tidy file: `intermediate-files/jeffreys-2013-znf-sequences.tsv`
```
# remove extra characters and tidy file
awk '{if ($1 ~ /[A-T,V-Z]/) print $1 "\t" $3; else print $1 "\t" $2}' copy-paste-files/jeffreys-2013-znf-copy.txt > intermediate-files/jeffreys-2013-znf-sequences.tsv

# check if all unique
wc -l intermediate-files/jeffreys-2013-znf-sequences.tsv
# 64
cut -f2 intermediate-files/jeffreys-2013-znf-sequences.tsv | sort | uniq | wc -l
# 64
```

Allele znf content:
- Study looks at low-frequency mutations in blood and genotypes of sperm cells; these alleles not necessarily observed as part of a human genotype
- Includes alleles `Jeffreys001`-`Jeffreys559` (publication did not provide allele names, so name them here)
- Copy/paste from **Supplementary Table S1** to: `copy-paste-files/jeffreys-2013-allele-copy.txt`
- Tidy file: `intermediate-files/jeffreys-2013-allele-znf-content.tsv`
```
# convert from 2 'text columns' to one, remove extra columns and duplicated alleles, remove gaps, add temporary allele names Je_###
egrep -vi "fig|type|no|man" copy-paste-files/jeffreys-2013-allele-copy.txt | awk '{if (length($5) >= 4) print $1 "\n" $5; else if (length($6) >= 4) print $1 "\n" $6; else if (NF < 5 && length($3) >= 4) print $1 "\n" $3; else if (NF < 6 && length($3) < 4) print $1}' | sed 's/-//g' | sort | uniq | awk '{printf "Jeffreys%03i\t%s\n", NR, $1}' > intermediate-files/jeffreys-2013-allele-znf-content.tsv
```

#

### Hussin et al. Mar 2013
### Rare allelic forms of _PRDM9_ associated with childhood leukemogenesis
**PMID: [23222848](https://pubmed.ncbi.nlm.nih.gov/23222848)**\
**GenBank Accession Numbers: [JQ044371.1 – JQ044377.1](https://ncbi.nlm.nih.gov/nuccore/?term=JQ044371.1%3AJQ044377.1%5Baccn%5D)**

Znf DNA sequences:
- Includes znfs `a`-`x`
- Copy/paste from **Supplementary Figure S6** to: `copy-paste-files/hussin-2013-znf-copy.txt`
- Tidy file: `intermediate-files/hussin-2013-znf-sequences.tsv`
```
# remove extra lines and tidy file
grep -v "Zinc" copy-paste-files/hussin-2013-znf-copy.txt | sed 's/\s/\t/' > intermediate-files/hussin-2013-znf-sequences.tsv

# check if all unique
wc -l intermediate-files/hussin-2013-znf-sequences.tsv
# 24
cut -f2 intermediate-files/hussin-2013-znf-sequences.tsv | sort | uniq | wc -l
# 24
```

Allele znf content:
- Includes alleles `L32`-`L37`
- Copy/paste from **Supplementary Material Page 7: Supplementary Results, Description of *PRDM9* Alleles and Novel ZnF Types** to `copy-paste-files/hussin-2013-allele-copy.txt`
- Tidy file: `intermediate-files/hussin-2013-allele-sequences.tsv`
```
# tidy file and sort alphabetically
sed -e 's/ is /\t/' -e 's/[,=]/\t/' copy-paste-files/hussin-2013-allele-copy.txt | sort -k1,1V > intermediate-files/hussin-2013-allele-sequences.tsv
```

Allele DNA sequence accession numbers:
- Save accession numbers to: `genbank-records/hussin-2013-allele-sequence-accessions.txt`
```
# generate sequence of numbers representing genbank accession numbers
for i in $(seq 44371 44377)
do
echo "JQ0$i.1" >> genbank-records/hussin-2013-allele-sequence-accessions.txt
done
```

#

### Beyter et al. May 2021
### Long-read sequencing of 3,622 Icelanders provides insight into the role of structural variants in human diseases and other traits
**PMID: [33972781](https://pubmed.ncbi.nlm.nih.gov/33972781)**\
**GenBank Accession Numbers: None**

Allele structural variants:
- Includes SVs `chr5:23526974:DN.1`, `chr5:23526974:DN.2`, `chr5:23526974:DN.4`, `chr5:23526974:XN.5`, `chr5:23527530:FN.0`, `chr5:23527530:FN.1`, `chr5:23527530:FN.4`
  - `chr5:23526974:XN.5` has no alternate allele and is therefore the reference sequence
- Download SV vcf from **[github](https://github.com/DecodeGenetics/LRS_SV_sets/raw/master/ont_sv_high_confidence_SVs.sorted.vcf.gz)** to: `copy-paste-files/beyter-2021-allele-SVs-copy.vcf`
- Tidy file: `intermediate-files/beyter-2021-allele-SVs.vcf`
```
# download vcf
wget https://github.com/DecodeGenetics/LRS_SV_sets/raw/master/ont_sv_high_confidence_SVs.sorted.vcf.gz -O copy-paste-files/beyter-2021-allele-SVs-copy.vcf.gz

# index vcf, subset to PRDM9 znf region
tabix copy-paste-files/beyter-2021-allele-SVs-copy.vcf.gz
tabix -h copy-paste-files/beyter-2021-allele-SVs-copy.vcf.gz chr5:23526673-23527764 > intermediate-files/beyter-2021-allele-SVs.vcf
```

Allele DNA sequences:
- Includes alleles `chr5:23526974:DN.1`, `chr5:23526974:DN.2`, `chr5:23526974:DN.4`, `chr5:23526974:XN.5`, `chr5:23527530:FN.0`, `chr5:23527530:FN.1`, `chr5:23527530:FN.4`
- Modify reference sequence to replace reference sequences with alternate SV sequences from `intermediate-files/beyter-2021-allele-SVs.vcf`
- Tidy file: `intermediate-files/beyter-2021-allele-sequences.tsv`
```
# download GRCh38 chr5 fasta
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr5.fa.gz -O copy-paste-files/GRCh38-chr5.fa.gz
zcat copy-paste-files/GRCh38-chr5.fa | grep -v ">" | tr -d '\n' | sed 's/$/\n/' | gzip > copy-paste-files/GRCh38-chr5.seq.gz

# replace reference allele sequence (znf domain) with alternate SV sequences; but leave sequence as is for variant with no alt sequence (reference sequence)
while read CHR POS ID REF ALT REMAINING
do
if [[ $ALT == "<ALT>" ]]
then awk -v ID=$ID '{print ID"\t"substr($0,23526673,23527764-23526673+1)}' <(zcat copy-paste-files/GRCh38-chr5.seq.gz) >> intermediate-files/beyter-2021-allele-sequences.tsv
else
awk -v ID=$ID -v ALT=$ALT -v POS=$POS -v REF=$REF '{print ID"\t"substr($0,23526673,POS-23526673)""ALT""substr($0,POS+length(REF),23527764-POS-length(REF)+1)}' <(zcat copy-paste-files/GRCh38-chr5.seq.gz) >> intermediate-files/beyter-2021-allele-sequences.tsv
fi
done < <(grep -v "#" intermediate-files/beyter-2021-allele-SVs.vcf)
```

#

### Wang et al. Jul 2021
### Pathogenic variants of meiotic double strand break (DSB) formation genes _PRDM9_ and _ANKRD31_ in premature ovarian insufficiency
**PMID: [34257419](https://pubmed.ncbi.nlm.nih.gov/34257419)**\
**GenBank Accession Numbers: None**

Allele mutations:
- Includes point mutations `c.229C>T:p.Arg77*`, `c.638T>G:p.Ile213Ser`, `c.677A>T:p.Lys226Met` relative to allele B (NM_020227.3)
- Copy/paste from **Table 1** to: `copy-paste-files/wang-2021-allele-mutations-copy.txt`
  - Only copy/pasted from `Patient number` to `E2,pg/mL` for patients 1-4 due to merged cells in publication table
- Tidy file: `intermediate-files/wang-2021-allele-mutations.tsv`
```
# remove extra columns, split position, reference and alternate alleles into separate columns, and sort alphabetically
cut -f3 copy-paste-files/wang-2021-allele-mutations-copy.txt  | sed 's/\s//' | uniq | awk -F: '{print $0 "\t" substr($0, 3,3 ) "\t" substr($0, 6, 1) "\t" substr($0, 8, 1)}' > intermediate-files/wang-2021-allele-mutations.tsv
```

Reference allele:
- Above point mutations relative to allele B ([NM_020227.3](https://ncbi.nlm.nih.gov/nucore/NM_020227.3))
- Copy/paste NCBI fasta record to: `copy-paste-files/NM_020227.3-copy.txt`
- Tidy file: `intermediate-files/NM_020227.3.fa`
```
# collapse reference sequence to single lines
sed '/>/ s/$/NEWLINE/' copy-paste-files/NM_020227.3-copy.txt | tr -d '\n' | sed 's/>/\n>/g' | sed 's/NEWLINE/\n/' > intermediate-files/NM_020227.3.fa
```

Allele sequences:
- Modify reference sequence to incorporate each allele mutation from `intermediate-files/wang-2021-allele-mutations.tsv` and save as allele sequence
- Tidy file: `intermediate-files/wang-2021-allele-sequences.tsv`
```
# generate allele sequence based on point mutation data
while read variant position ref alt
do
grep -v ">" intermediate-files/NM_020227.3.fa |  sed -e "s|\(.\{$position\}\).|\1$alt|" -e "s/^/$variant\t/" >> intermediate-files/wang-2021-allele-sequences.tsv
done < intermediate-files/wang-2021-allele-mutations.tsv
```

#

### Alleva et al. Nov 2021
### Cataloging human _PRDM9_ allelic variation using long-read sequencing reveals _PRDM9_ population specificity and two distinct groupings of related alleles
**PMID: [34805134](https://pubmed.ncbi.nlm.nih.gov/34805134)**\
**GenBank Accession Numbers: None**

Znf DNA sequences:
- Includes `!A`-`!N`, `:A`-`:V`, `|1`-`|9`, `|A`-`|J`, `|a`-`|j`
- Save **Supplementary Data File 2** to: `copy-paste-files/alleva-2021-SD2-znf-copy.tsv`
- Tidy file: `intermediate-files/alleva-2021-znf-sequences.tsv`
```
# remove extra lines and columns
grep "TGT" copy-paste-files/alleva-2021-SD2-znf-copy.tsv | cut -f1,4 > intermediate-files/alleva-2021-znf-sequences.tsv

# check if all unique
wc -l intermediate-files/alleva-2021-znf-sequences.tsv
# 81
cut -f2 intermediate-files/alleva-2021-znf-sequences.tsv | sort | uniq | wc -l
# 81
```

Znf amino acid sequences:
- Includes `!A`-`!N`, `:A`-`:V`, `|1`-`|9`, `|A`-`|J`, `|a`-`|j`
- From the same **Supplementary Data File 2** for znf DNA sequences: `copy-paste-files/alleva-2021-SD2-znf-copy.tsv`
- Tidy file: `intermediate-files/alleva-2021-znf-aminos.tsv`
```
# remove extra lines and columns
grep "TGT" copy-paste-files/alleva-2021-SD2-znf-copy.tsv | cut -f1,5 > intermediate-files/alleva-2021-znf-aminos.tsv
```

Allele znf content:
- Includes alleles `A`-`E`, `F`, `H`-`I`, `L1`-`L27`, `M1`-`M32`
  - Also includes 542 additional alleles from [Jeffreys et al. 2013](#jeffreys-et-al-jan-2013) observed in sperm or blood
- Save **Supplementary Data File 3** to: `copy-paste-files/alleva-2021-SD3-allele-copy.tsv`
- Tidy file: `intermediate-files/alleva-2021-allele-znf-content.tsv`
```
# remove extra lines and columns
grep "TGT" copy-paste-files/alleva-2021-SD3-allele-copy.tsv | cut -f2,3 > intermediate-files/alleva-2021-allele-znf-content.tsv
```

Allele DNA sequences:
- Includes alleles `A`-`E`, `F`, `H`-`I`, `L1`-`L27`, `M1`-`M32`
  - Also includes 542 additional alleles from [Jeffreys et al. 2013](#jeffreys-et-al-jan-2013) observed in sperm or blood
- From the same **Supplementary Data File 3** for allele znf content: `copy-paste-files/alleva-2021-SD3-allele-copy.tsv`
- Tidy file: `intermediate-files/alleva-2021-allele-sequences.tsv`
```
# remove extra lines and columns
grep "TGT" copy-paste-files/alleva-2021-SD3-allele-copy.tsv | cut -f2,6 > intermediate-files/alleva-2021-allele-sequences.tsv
```

---

## Step 2. Compile znf sequences and get list of unique znfs
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

# separate author from znf name
# shift Oliver 2009 sequences to match the rest
# pivot wider to obtain unique znf sequences as rows
# arrange in publication order (not including Oliver)
znf.sequences.sperm <- pub.znf %>%
  separate(PubZnfName, sep="_", into=c("Publication", "ZnfName"), extra="merge") %>%
  mutate(Publication=sub("-", ".", Publication)) %>%
  mutate(Sequence=ifelse(Publication=="oliver.2009", sub("TGCAGGGAG", "", Sequence), Sequence)) %>%
  mutate(Sequence=ifelse(Publication=="oliver.2009", sub("$", "TGCAGGGAG", Sequence), Sequence)) %>%
  pivot_wider(names_from=Publication, values_from=ZnfName) %>%
  arrange(berg.2010, berg.2011, borel.2012, jeffreys.2013, hussin.2013, alleva.2021)

write.table(znf.sequences.sperm, "publication-unique-znf-sequences-with-somatic-and-sperm.tsv", row.names=F, quote=F, sep="\t")

remove.germline <-znf.sequences %>%
  filter(across(c(alleva.2021, jeffreys.2013), ~ !is.na(.x))) %>%
  filter(across(c(-alleva.2021, -jeffreys.2013, -Sequence), is.na))

znf.sequences <- znf.sequences.sperm %>%
  anti_join(remove.germline)

write.table(znf.sequences, "publication-unique-znf-sequences.tsv", row.names=F, quote=F, sep="\t")
```

Give standardized names to unique znf sequences (ZN###) and to somatic/sperm znf sequences (zn$$$):
```
tail -n +2 intermediate-files/publication-unique-znf-sequences.tsv | awk '{printf "ZN%03i\t%s\n", NR, $1}' > intermediate-files/standardized-znf-sequences-list.tsv
cut -f2 intermediate-files/standardized-znf-sequences-list.tsv > intermediate-files/standardized-znf-seqs.txt
grep -vf intermediate-files/standardized-znf-seqs.txt intermediate-files/publication-unique-znf-sequences-with-somatic-and-sperm.tsv | grep -v Sequence | awk '{printf "zn%03i\t%s\n", NR, $1}' > intermediate-files/standardized-znf-sequences-list-somatic-and-sperm.tsv
```

---

## Step 3. Compile known znf sequences and allele znf content

Znf DNA sequences:
- Combine Berg et al. 2010 and Berg et al. 2011, as the first contains znfs m & n and the second contains znfs u & v
- Give Hussin et al. znfs u & v temporary names U & V, since they differ from Berg 2011 znfs u & v
- Add Hussin et al znf sequences to the Berg sequences
- Keep unique sequences
```
# combine lists of znf sequences
cat intermediate-files/berg-2010-znf-sequences.tsv intermediate-files/berg-2011-znf-sequences.tsv <(sed -e 's/u/U/' -e 's/v/V/' intermediate-files/hussin-2013-znf-sequences.tsv) | sort -f | uniq > intermediate-files/publication-znf-sequences.tsv

wc -l intermediate-files/publication-znf-sequences.tsv
# 26 unique znf sequences in total
```

Allele znf content:
- Add Hussin et al. to Berg et al. 2011, as Berg et al. 2011 has all alleles in Berg et al. 2010, plus two additional alleles
- Ignore Ponting because it contains all alleles in Berg et al. 2010, one of which is incorrect
- Ignore Baudat et al. for now because not all znfs can be deduced from the image alone
```
# combine lists of allele znf content
cat intermediate-files/berg-2011-allele-znf-content.tsv intermediate-files/hussin-2013-allele-znf-content.tsv | grep . > intermediate-files/publication-allele-znf-content.tsv

wc -l intermediate-files/publication-allele-znf-content.tsv
# 39 unique alleles in total
```

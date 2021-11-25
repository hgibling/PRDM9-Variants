# Collecting known _PRDM9_ alleles and associated zinc fingers (znfs) from publications
Literature searches revealed several publications that describe znf DNA sequences, allele znf content, allele DNA sequence accession numbers, and/or allele mutations:
- [Oliver 2009](#oliver-et-al-dec-2009)
  - Znf DNA sequences
  - Allele DNA sequence accession numbers
- [Baudat 2010](#baudat-et-al-feb-2010)
  - Allele znf content
  - Allele DNA sequence accession numbers
- [Berg 2010](#berg-et-al-oct-2010)
  - Znf DNA sequences
  - Allele znf content
  - Allele DNA sequence accession numbers
- [Ponting 2011](#ponting-may-2011)
  - Allele znf content
- [Berg 2011](#berg-et-al-jul-2011)
  - Znf DNA sequences
  - Allele znf content
- [Borel 2012](#borel-et-al-may-2012)
  - Znf DNA sequences
  - Allele znf content
- [Hussin 2013](#hussin-et-al-mar-2013)
  - Znf DNA sequences
  - Allele znf content
  - Allele DNA sequence accession numbers
- [Wang 2021](#wang-et-al-jul-2021)
  - Allele mutations

When possible, sequences were copy/pasted from publications and saved as `firstauthor-year-type.txt` in the `copy-paste-files` directory. After tidying up, content was saved as `firstauthor-year-type.tsv` in the `intermediate-files` directory. 

Genbank accession downloads are described in [additional documentation](/collect-known-alleles-genbank.md).

Analysis steps:
1. [Get allele and znf sequence data from publications](#1-get-allele-and-znf-sequence-data-from-publications)
2. [Check if allele content and znf sequences are unique](#2-check-if-allele-content-and-znf-sequences-are-unique)
3. [Compile known znf sequences and allele znf content](#3-compile-known-znf-sequences-and-allele-znf-content)

---

## 1. Get allele and znf sequence data from publications

### Oliver et al. Dec 2009
### Accelerated Evolution of the _Prdm9_ Speciation Gene across Diverse Metazoan Taxa
**PMID: [19997497](https://pubmed.ncbi.nlm.nih.gov/19997497)**\
**GenBank Accession Numbers: [FJ899863.1 - FJ899912.1](https://ncbi.nlm.nih.gov/nuccore/?term=FJ899863.1%3AFJ899912.1%5Baccn%5D)**

Znf DNA sequences:
- Includes znfs 03-14
- Copy/paste from **Supplementary Dataset S1** to: `copy-paste-files/oliver-2009-znf-copy.txt`
- Tidy file: `intermediate-files/oliver-2009-znf-sequences.tsv`
```
# extract human znfs
sed '/>/ s/$/NEWLINE/' copy-paste-files/oliver-2009-znf-copy.txt | tr -d '\n' | sed 's/>/\n/g' | grep "homo_sapien" | sed 's/NEWLINE/\t/' > intermediate-files/oliver-2009-znf-sequences.tsv
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

### Baudat et al. Feb 2010
### _PRDM9_ is a major determinant of meiotic recombination hotspots in humans and mice
**PMID: [20044539](https://pubmed.ncbi.nlm.nih.gov/20044539)**\
**GenBank Accession Numbers: [GU216222.1 - GU216229.1](https://ncbi.nlm.nih.gov/nuccore/?term=GU216222.1%3AGU216229.1%5Baccn%5D)**

Allele znf content:
- Includes alleles A-E, F, H-I, and K
- Image in **Figure 2b** depicts znf content as unnamed blocks
  - Appears one block (`-NHR` in the legend) represents more than one znf as it appears in positions 11 & 12 in alleles A and B
  - Additionally, tne block in alleles F and K (`--HR`) and one block in allele H (`-RVS`)that are not also present in A-E
- Type znf content typed out by hand, triple check for accuracy, to: `copy-paste-files/baudat-2010-allele-copy.txt`
  - Name blocks by four-character code in figure legend (orange = `NTOR`; grey = `NTGR`), with `|` to separate blocks
- Tidy file: `intermediate-files/baudat-2010-allele-copy.txt`

<img src="https://github.com/hgibling/PRDM9-Variants/blob/main/images/Baudat-2010-Figure2B.png?raw=true" width="600">

```
# tidy file
awk '{print $1 "\t" $3}' copy-paste-files/baudat-2010-allele-copy.txt > intermediate-files/baudat-2010-allele-content-unnamed.tsv
```

Allele DNA sequence accession numbers:
- Includes alleles A-F, H-I
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
- Includes znfs a-t
- Copy/paste from **Supplementary Figure 1a** to: `copy-paste-files/berg-2010-znf-copy.txt`
- Tidy file: `intermediate-files/berg-2010-znf-sequences.tsv`
```
# convert znf names to lowercase and tidy file
sed -e 's/^\(.\)/\L\1/' -e 's/\s/\t/'  -e '$a\' copy-paste-files/berg-2010-znf-copy.txt > intermediate-files/berg-2010-znf-sequences.tsv
```

Allele znf content:
- Includes alleles A-E, L1-L24
- Copy/paste from **Supplementary Figure 1b** to: `copy-paste-files/berg-2010-allele-copy.txt`
- Tidy file: `intermediate-files/berg-2010-allele-content.tsv`

```
# remove extra columns, convert znf names to lowercase, tidy file, and sort alphabetically
awk '{print $1 "\t" tolower($3)}' copy-paste-files/berg-2010-allele-copy.txt | sort -k1,1V > intermediate-files/berg-2010-allele-sequences.tsv
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

### Ponting May 2011
### What are the genomic drivers of the rapid evolution of _PRDM9_?
**PMID: [21388701](https://pubmed.ncbi.nlm.nih.gov/21388701)**\
**GenBank Accession Numbers: None**

Allele znf content:
- Includes alleles A-L24
- Review paper
- Image in **Figure 4** depicts znf content as named blocks
  - Cites Berg et al. 2010 for allele znf content
- Type Znf content out by hand, triple check for accuracy, to: `copy-paste-files/ponting-2010-allele-content.txt`
  - Include gaps represented with `_`
- Tidy file: `intermediate-files/ponting-2011-allele-content.tsv`
```
# remove gaps and sort alphabetically
sed -e 's/\s/\t/g' -e 's/_//g' copy-paste-files/ponting-2011-allele-copy.txt | cut -f1,2 | sort -k1,1V > intermediate-files/ponting-2011-allele-content.tsv
```

<img src="https://github.com/hgibling/PRDM9-Variants/blob/main/images/Ponting-2011-Figure3.jpg?raw=true" width="500">

#

### Berg et al. Jul 2011
### Variants of the protein _PRDM9_ differentially regulate a set of human meiotic recombination hotspots highly active in African populations
**PMID: [21750151](https://pubmed.ncbi.nlm.nih.gov/21750151)**\
**GenBank Accession Numbers: None**


Znf DNA sequences
- Includes znfs a-l, o-v
- Copy/paste from **Supplementary Figure 1A** to: `copy-paste-files/berg-2011-znf-copy.txt`
- Tidy file: `intermediate-files/berg-2011-znf-sequences.tsv`
```
# convert znf names to lowercase and tidy file
sed -e 's/^\(.\)/\L\1/' -e 's/\s/\t/'  -e '$a\' copy-paste-files/berg-2011-znf-copy.txt > intermediate-files/berg-2011-znf-sequences.tsv
```

Allele znf content:
- Includes alleles A-E, L1-L27
- Copy/paste from **Supplementary Figure 1B** to: `copy-paste-files/berg-2011-allele-copy.txt`
- Tidy file: `intermediate-files/berg-2010-allele-content.tsv`

```
# remove extra columns, convert znf names to lowercase, tidy file, and sort alphabetically
awk '{print $1 "\t" tolower($3)}' copy-paste-files/berg-2011-allele-copy.txt | sort -k1,1V > intermediate-files/berg-2011-allele-content.tsv
```

#

### Borel et al. May 2012
### Evaluation of _PRDM9_ variation as a risk factor for recurrent genomic disorders and chromosomal non-disjunction
**PMID: [22643917](https://pubmed.ncbi.nlm.nih.gov/22643917)**\
**GenBank Accession Numbers: None**

Znf DNA sequences:
- Includes znfs a-m, q
- Copy/paste from **Supplementary Table S1** to: `copy-paste-files/borel-2012-znf-copy.txt`
- Tidy file: `intermediate-files/borel-2012-znf-sequences.tsv`
```
# convert znf names to lowercase and tidy file
sed -e 's/^\(.\)/\L\1/' -e 's/\s/\t/' copy-paste-files/borel-2012-znf-copy.txt > intermediate-files/borel-2012-znf-sequences.tsv
```

Allele znf content:
- Includes alleles A-F, I, L1, L19, L28-L31
- Copy/paste from **Supplementary Table S1** to: `copy-paste-files/borel-2012-allele-copy.txt`
- Tidy file: `intermediate-files/borel-2012-allele-content.tsv`

```
# remove extra column and convert znf names to lowercase
awk '{print $1 "\t" tolower($3)}' copy-paste-files/borel-2012-allele-copy.txt > intermediate-files/borel-2012-allele-sequences.tsv
```

#

### Hussin et al. Mar 2013
### Rare allelic forms of _PRDM9_ associated with childhood leukemogenesis
**PMID: [23222848](https://pubmed.ncbi.nlm.nih.gov/23222848)**\
**GenBank Accession Numbers: [JQ044371.1 – JQ044377.1](https://ncbi.nlm.nih.gov/nuccore/?term=JQ044371.1%3AJQ044377.1%5Baccn%5D)**

Znf DNA sequences
- Includes znfs a-x
- Copy/paste from **Supplementary Figure S6** to: `copy-paste-files/hussin-2013-znf-copy.txt`
- Tidy file: `intermediate-files/hussin-2013-znf-sequences.tsv`
```
# remove extra lines and tidy file
grep -v "Zinc" copy-paste-files/hussin-2013-znf-copy.txt | sed 's/\s/\t/' > intermediate-files/hussin-2013-znf-sequences.tsv
```

Allele znf content
- Includes alleles L32-L37
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

### Wang et al. Jul 2021
### Pathogenic variants of meiotic double strand break (DSB) formation genes _PRDM9_ and _ANKRD31_ in premature ovarian insufficiency
**PMID: [34257419](https://pubmed.ncbi.nlm.nih.gov/34257419)**\
**GenBank Accession Numbers: None**

Allele mutations:
- Includes point mutations c.229C>T:p.Arg77*, c.638T>G:p.Ile213Ser, c.677A>T:p.Lys226Met relative to allele B (NM_020227.3)
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
- Modify reference sequence to incorporate each allele mutation and save as allele sequence
- Tidy file: `intermediate-files/wang-2021-allele-sequences.tsv`
```
# generate allele sequence based on point mutation data
while read variant position ref alt
do
grep -v ">" intermediate-files/NM_020227.3.fa |  sed -e "s|\(.\{$position\}\).|\1$alt|" -e "s/^/$variant\t/" >> intermediate-files/wang-2021-allele-sequences.tsv
done < intermediate-files/wang-2021-allele-mutations.tsv
```

---

## 2. Check if allele content and znf sequences are unique

---

## 3. Compile known znf sequences and allele znf content

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
cat intermediate-files/berg-2011-allele-content.tsv intermediate-files/hussin-2013-allele-content.tsv | grep . > intermediate-files/publication-allele-content.tsv

wc -l intermediate-files/publication-allele-content.tsv
# 39 unique alleles in total
```

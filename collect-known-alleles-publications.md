# Collecting known _PRDM9_ alleles and associated zinc fingers (znfs) from publications
Literature searches revealed several publications that describe one or more of the following:
- **Znf DNA sequences**
  - Berg 2010
  - Berg 2011
  - Hussin 2013
- **Allele znf content**
  - Berg 2010
  - Baudat 2010
  - Berg 2011
  - Ponting 2011
  - Hussin 2013
- **Allele DNA sequences**
  - Baudat 2010

When possible, sequences were copy/pasted from publications and saved as `firstauthor-year-type.txt` in the `copy-paste-files` directory. After tidying up, content was saved as `firstauthor-year-type.tsv` in the `intermediate-files` directory. 

---

## Publications describing Znf DNA sequences & Allele znf content

### Berg et al. Oct 2010
### _PRDM9_ variation strongly influences recombination hot-spot activity and meiotic instability in humans
**PMID: 20818382**\
**GenBank Accession Numbers: HM210983.1 – HM211006.1**

Znf DNA sequences:
- Includes znfs a-t
- Copy/pasted from **Supplementary Figure 1a** to `copy-paste-files/berg-2010-znf-copy.txt`
- Final file: `intermediate-files/berg-2010-znf-sequences.tsv`
```
# convert znf names to lowercase and tidy file
sed -e 's/^\(.\)/\L\1/' -e 's/\s/\t/'  -e '$a\' copy-paste-files/berg-2010-znf-copy.txt > intermediate-files/berg-2010-znf-sequences.tsv
```

Allele znf content:
- Includes alleles A-E, L1-L24
- Copy/pasted from **Supplementary Figure 1b** to `copy-paste-files/berg-2010-allele-copy.txt`
- Final file: `intermediate-files/berg-2010-allele-content.tsv`

```
# remove extra columns, convert znf names to lowercase, tidy file, and sort alphabetically
awk '{print $1 "\t" tolower($3)}' copy-paste-files/berg-2010-allele-copy.txt | sort -k1,1V > intermediate-files/berg-2010-allele-sequences.tsv
```

#

### Berg et al. Jul 2011
### Variants of the protein _PRDM9_ differentially regulate a set of human meiotic recombination hotspots highly active in African populations
**PMID: 21750151**\
**GenBank Accession Numbers: None**


Znf DNA sequences
- Includes previously described znfs a-l, o-t, plus newly described znfs u-v
- Copy/pasted from **Supplementary Figure 1A** to `copy-paste-files/berg-2011-znf-copy.txt`
- Final file: `intermediate-files/berg-2010-znf-sequences.tsv`
```
# convert znf names to lowercase and tidy file
sed -e 's/^\(.\)/\L\1/' -e 's/\s/\t/'  -e '$a\' copy-paste-files/berg-2011-znf-copy.txt > intermediate-files/berg-2011-znf-sequences.tsv

# check that znfs are identical to those from Berg et al 2010
diff intermediate-files/berg-2010-znf-sequences.tsv intermediate-files/berg-2011-znf-sequences.tsv

# 13,14d12
# < m     TGTGGGCGGGGCTTTAGAGATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAG
# < n     TGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAG
# 20a19,20
# > u     TGTGGGCGGGGCTTTAGCGATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAG
# > v     TGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCAGTCACCAGAGGACACACACAGGGAAGAAGCCCTATGTCTGCAGGGAG

# znfs m-n are missing, znfs u-v new to Berg 2011
```

Allele znf content:
- Includes previously described alleles A-E, L1-L24, plus newly described alleles L25-L27
- Copy/pasted from **Supplementary Figure 1B** to `copy-paste-files/berg-2011-allele-copy.txt`
- Final file: `intermediate-files/berg-2010-allele-content.tsv`

```
# remove extra columns, convert znf names to lowercase, tidy file, and sort alphabetically
awk '{print $1 "\t" tolower($3)}' copy-paste-files/berg-2011-allele-copy.txt | sort -k1,1V > intermediate-files/berg-2011-allele-content.tsv

# check that allele contents are identical to those from Berg 2010
diff intermediate-files/berg-2010-allele-content.tsv intermediate-files/berg-2011-allele-content.tsv

# 29a30,32
# > L25   abcddecfuhfij
# > L26   abcddecfghfqj
# > L27   abcddvcfghfij

# alleles L25-L27 new to Berg 2011
```

#

### Hussin et al. Mar 2013
### Rare allelic forms of _PRDM9_ associated with childhood leukemogenesis
**PMID: 23222848**\
**GenBank Accession Numbers: JQ044371 – JQ044377**

Znf DNA sequences
- Includes previously described znfs a-t, plus newly described znfs u-w
- However, Hussin et al. 2013 znfs u-v **do not** match Berg et al. 2011 znfs u-v
- Copy/pasted from **Supplementary Figure S6** to `copy-paste-files/hussin-2013-znf-copy.txt`
- Final file: `intermediate-files/hussin-2013-znf-sequences.tsv`
```
# remove extra lines and tidy file
grep -v "Zinc" copy-paste-files/hussin-2013-znf-copy.txt | sed 's/\s/\t/' > intermediate-files/hussin-2013-znf-sequences.tsv

# check that znf sequences match those from Berg 2010, 2011
diff intermediate-files/berg-2010-znf-sequences.tsv intermediate-files/hussin-2013-znf-sequences.tsv

# 20a21,24
# > u     TGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGAAGAAGCCCTATGTCTGCAGGGAG
# > v     TGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTTTGTCTGCAGGGAG
# > w     TGTGGGCGGGGCTTTCTCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAG
# > x     TGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAACCCTATGTCTGCAGGGAG

# znfs a-t match those from Berg 2010, and znfs u-x are new to Hussin 2013

diff intermediate-files/berg-2011-znf-sequences.tsv intermediate-files/hussin-2013-znf-sequences.tsv

# 12a13,14
# > m     TGTGGGCGGGGCTTTAGAGATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAG
# > n     TGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAG
# 19,20c21,24
# < u     TGTGGGCGGGGCTTTAGCGATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAG
# < v     TGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCAGTCACCAGAGGACACACACAGGGAAGAAGCCCTATGTCTGCAGGGAG
# ---
# > u     TGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGAAGAAGCCCTATGTCTGCAGGGAG
# > v     TGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTTTGTCTGCAGGGAG
# > w     TGTGGGCGGGGCTTTCTCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAG
# > x     TGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAACCCTATGTCTGCAGGGAG

# znfs m-n not present in Berg 2011; znfs w-x new to Hussin 2013
# IMPORTANT: znfs u-v have different sequences between Berg 2011 and Hussin 2013
```

Allele znf content
- Includes alleles L32-L37
- Copy/pasted from **Supplementary Material Page 7: Supplementary Results, Description of *PRDM9* Alleles and Novel ZnF Types** to `copy-paste-files/hussin-2013-allele-copy.txt`
- Final file: `intermediate-files/hussin-2013-allele-sequences.tsv`
```
# tidy file and sort alphabetically
sed -e 's/ is /\t/' -e 's/[,=]/\t/' copy-paste-files/hussin-2013-allele-copy.txt | sort -k1,1V > intermediate-files/hussin-2013-allele-sequences.tsv
```

---

## Publications describing Allele znf content in alleles and Allele DNA sequences


### Baudat et al. Feb 2010
### _PRDM9_ is a major determinant of meiotic recombination hotspots in humans and mice
**PMID: 20044539**\
**GenBank Accession Numbers: GU216222.1 - GU216229.1**

Allele znf content:
- Includes alleles A-E described in Berg et al. 2010, plus alleles F, H-I, and K
- Image in **Figure 2b** depicts znf content as unnamed blocks
  - Appears one block (`-NHR` in the legend) represents more than one znf as it appears in positions 11 & 12 in alleles A and B, which corresponds to znfs h and i
  - Additionally, tne block in alleles F and K (`--HR`) and one block in allele H (`-RVS`)that are not also present in A-E
    - The znf content for these alleles cannot be deduced from the Berg et al. 2010 znf content
    - Need to obtain DNA sequences for these alleles and determine znf labels from there
  - Nevertheless, Znf content typed out by hand, triple checked for accuracy, to `copy-paste-files/Baudat-2010-allele-copy.txt`
    - Named blocks by four-character code in figure legend (orange = `NTOR`; grey = `NTGR`), with `|` to separate blocks

<img src="https://github.com/hgibling/PRDM9-Variants/blob/main/images/Baudat-2010-Figure2B.png?raw=true" width="600">

```
# tidy file
awk '{print $1 "\t" $3}' copy-paste-files/baudat-2010-allele-copy.txt > intermediate-files/baudat-2010-allele-content-unnamed.tsv

# print znf content for alleles A-E from Berg 2010
egrep "[AB]" intermediate-files/berg-2010-allele-content.tsv 

# A       abcddecfghfij
# B       abcddccfghfij

# alleles A and B have the same final 3 zinc fingers, of which those at positions 11 and 12 are znfs h and i, respectively
```

Allele DNA sequences:
- Includes alleles A-F, H-I
- Sequences are fasta screenshots in **Supplementary Figure S3B** instead of copy/pastable text (!)
- Parse sequences from NCBI accessions instead

---

## Publications describing Allele znf content


### Ponting May 2011
### What are the genomic drivers of the rapid evolution of _PRDM9_?
**PMID: 21388701**\
**GenBank Accession Numbers: None**

Allele znf content:
- Includes alleles A-L24
- Review paper
- Image in **Figure 4** depicts znf content as named blocks
  - Figure cites Berg et al. 2010 for allele znf content
  - However, allele znf content for L24 **does not** match that for Berg et al. 2010
  - Via email communication with Ponting (Aug 2021), confirmed that L24 znf content depicted in Figure 4 is **incorrect**
- Znf content typed out by hand, triple checked for accuracy, to `copy-paste-files/Ponting-2010-allele-content.txt`
  - Included gaps, represented with `_`
- Final file: `intermediate-files/Ponting-2011-allele-content.tsv`
```
# remove gaps and sort alphabetically
sed -e 's/\s/\t/g' -e 's/_//g' copy-paste-files/ponting-2010-allele-copy.txt | cut -f1,2 | sort -k1,1V > intermediate-files/ponting-2010-allele-content.tsv

# check that allele contents are identical to those from Berg 2010, the reference for the image content
diff intermediate-files/berg-2010-allele-content.tsv intermediate-files/ponting-2011-allele-content.tsv

# 29c29
# < L24   abcddecftpfqj
# ---
# > L24   abcddecfthpfqj

# IMPORTANT: content for L24 in Ponting 2011 has an h znf not present in Berg 2010
```

<img src="https://github.com/hgibling/PRDM9-Variants/blob/main/images/Ponting-2011-Figure3.jpg?raw=true" width="500">

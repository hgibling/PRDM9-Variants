# Collecting known _PRDM9_ alleles and associated zinc fingers (znfs)
Literature searches revealed several publications that describe one or more of the following:
- Znf DNA sequences
- Allele znf content
- Allele DNA sequences

---

## Publications describing Znf DNA sequences & Allele znf content

### Berg et al. Oct 2010
### _PRDM9_ variation strongly influences recombination hot-spot activity and meiotic instability in humans
**PMID: 20818382**

Znf sequences:
- Copy/pasted from **Supplementary Figure 1a** to `berg-2010-znf-copy.txt`
- Includes znfs a-t
- Final file: `berg-2010-znf-sequences.tsv`
```
# convert znf names to lowercase and tidy file
sed -e 's/^\(.\)/\L\1/' -e 's/\s/\t/'  -e '$a\' intermediate-files/berg-2010-znf-copy.txt > berg-2010-znf-sequences.tsv
```

Allele znf content:
- Copy/pasted from **Supplementary Figure 1b** to `berg-2010-allele-copy.txt`
- Includes alleles A-E, L1-L24
- Final file: `berg-2010-allele-content.tsv`

```
# remove extra columns, convert znf names to lowercase, tidy file, and sort alphabetically
awk '{print $1 "\t" tolower($3)}' intermediate-files/berg-2010-allele-copy.txt | sort -k1,1V > berg-2010-allele-sequences.tsv
```

#

### Berg et al. Jul 2011
### Variants of the protein _PRDM9_ differentially regulate a set of human meiotic recombination hotspots highly active in African populations
**PMID: 21750151**

Znf sequences
- Copy/pasted from **Supplementary Figure 1A** to `berg-2011-znf-copy.txt`
- Includes previously described znfs a-l, o-t, plus newly described znfs u-v
- Final file: `berg-2010-znf-sequences.tsv`
```
# convert znf names to lowercase and tidy file
sed -e 's/^\(.\)/\L\1/' -e 's/\s/\t/'  -e '$a\' intermediate-files/berg-2011-znf-copy.txt > berg-2011-znf-sequences.tsv

# check that znfs are identical to those from Berg et al 2010
diff berg-2010-znf-sequences.tsv berg-2011-znf-sequences.tsv

# 13,14d12
# < m     TGTGGGCGGGGCTTTAGAGATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAG
# < n     TGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAG
# 20a19,20
# > u     TGTGGGCGGGGCTTTAGCGATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAG
# > v     TGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCAGTCACCAGAGGACACACACAGGGAAGAAGCCCTATGTCTGCAGGGAG

# znfs m-n are missing, znfs u-v added, in berg-2011-znf-sequences.tsv
```

Allele znf content:
- Copy/pasted from **Supplementary Figure 1B** to `berg-2011-allele-copy.txt`
- Includes previously described alleles A-E, L1-L24, plus newly described alleles L25-L27
- Final file: `berg-2010-allele-content.tsv`

```
# remove extra columns, convert znf names to lowercase, tidy file, and sort alphabetically
awk '{print $1 "\t" tolower($3)}' intermediate-files/berg-2011-allele-copy.txt | sort -k1,1V > berg-2011-allele-content.tsv

# check that allelels are identical to those from Berg et al 2010
diff berg-2010-allele-content.tsv berg-2011-allele-content.tsv

# 29a30,32
# > L25   abcddecfuhfij
# > L26   abcddecfghfqj
# > L27   abcddvcfghfij

# alleles L25-L27 added in berg-2011-allele-content.tsv
```

#

### Hussin et al. Mar 2013
### Rare allelic forms of _PRDM9_ associated with childhood leukemogenesis
**PMID: 23222848**

Znf sequences
- Copy/pasted from **Supplementary Figure S6** to `hussin-2013-znf-copy.txt`
- Includes previously described znfs a-t, plus newly described znfs u-w
- However, Hussin et al. 2013 znfs u-v **do not** match Berg et al. 2011 znfs u-v
- Final file: `hussin-2013-znf-sequences.tsv`
```
# remove extra lines and tidy file
grep -v "Zinc" intermediate-files/hussin-2013-znf-copy.txt | sed 's/\s/\t/' > hussin-2013-znf-sequences.tsv

# check that znf sequences match those from Berg et al. 2010, 2011
diff berg-2010-znf-sequences.tsv hussin-2013-znf-sequences.tsv

# 20a21,24
# > u     TGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGAAGAAGCCCTATGTCTGCAGGGAG
# > v     TGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTTTGTCTGCAGGGAG
# > w     TGTGGGCGGGGCTTTCTCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAG
# > x     TGTGGGCGGGGCTTTAGCAATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAACCCTATGTCTGCAGGGAG

# znfs a-t match those from Berg et al 2010, and znfs u-x are new to Hussin 2013

diff berg-2011-znf-sequences.tsv hussin-2013-znf-sequences.tsv

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

# znfs m-n not present in Berg et al 2011; znfs w-x new to Hussin et al 2013
# IMPORTANT: znfs u-v have different sequences between Berg et al 2011 and Hussin et al 2013
```

Allele znf content
- Copy/pasted from **Supplementary Material Page 7: Supplementary Results, Description of *PRDM9* Alleles and Novel ZnF Types** to `hussin-2013-allele-copy.txt`
- Includes alleles L32-L37
- Final file: `hussin-2013-allele-sequences.tsv`
```
# tidy file and sort alphabetically
sed -e 's/ is /\t/' -e 's/[,=]/\t/' intermediate-files/hussin-2013-allele-copy.txt | sort -k1,1V > hussin-2013-allele-sequences.tsv
```

---

## Publications describing Allele znf content in alleles and Allele DNA sequences


### Baudat et al. Feb 2010
### _PRDM9_ is a major determinant of meiotic recombination hotspots in humans and mice
**PMID: 20044539**

Allele znf content:
- Includes alleles A-F, H-I, K
- Image in **Figure 2b** depicts znf content as unnamed blocks
  - Compare colored blocks to znf content in alleles A-E from Berg et al. 2010 to deduce which block color corresponds to which znf
  - Compare znf DNA sequences from inferred colored blocks to allele DNA sequences provided



Allele DNA sequences:
- Includes alleles A-F, H-I
- Sequences are fasta screenshots in **Supplementary Figure S3B** instead of copy/pastable text (!)
- Parse sequences from NCBI accessions instead

#

# Collecting known *PRDM9* alleles and associated zinc fingers (znfs)

## Berg et al. 2010
**PRDM9 variation strongly influences recombination hot-spot activity and meiotic instability in humans**\
**PMID: 20818382**
- Znf sequences: copy/pasted from **Supplementary Figure 1a** to `berg-2010-znf-copy.txt`
  - Includes znfs a-t
  - Final file: `berg-2010-znf-sequences.tsv`
```
# convert znf names to lowercase and tidy file
sed -e 's/^\(.\)/\L\1/' -e 's/\s/\t/'  -e '$a\' intermediate-files/berg-2010-znf-copy.txt > berg-2010-znf-sequences.tsv
```

- Allele sequences: copy/pasted from **Supplementary Figure 1b** to `berg-2010-allele-copy.txt`
  - Includes alleles A-E, L1-L24
  - Final file: `berg-2010-allele-sequences.tsv`

```
# remove extra columns, convert znf names to lowercase, tidy file, and sort alphabetically
awk '{print $1 "\t" tolower($3)}' intermediate-files/berg-2010-allele-copy.txt | sort -k1,1V > berg-2010-allele-sequences.tsv
```

## Berg et al. 2011
**Variants of the protein PRDM9 differentially regulate a set of human meiotic recombination hotspots highly active in African populations**\
**PMID: 21750151**
- Znf sequences: copy/pasted from **Supplementary Figure 1A** to `berg-2011-znf-copy.txt`
  - Includes previously described znfs a-l, o-t, plus znfs u-v
  - Final file: `berg-2010-znf-sequences.tsv`
```
# convert znf names to lowercase and tidy file
sed -e 's/^\(.\)/\L\1/' -e 's/\s/\t/'  -e '$a\' intermediate-files/berg-2011-znf-copy.txt > berg-2011-znf-sequences.tsv

# check that znfs are identical to those in berg-2010-znf-sequences.tsv
diff berg-2010-znf-sequences.tsv berg-2011-znf-sequences.tsv

# 13,14d12
# < m     TGTGGGCGGGGCTTTAGAGATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTACGTCTGCAGGGAG
# < n     TGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCAGTCACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAG
# 20a19,20
# > u     TGTGGGCGGGGCTTTAGCGATAAGTCACACCTCCTCAGACACCAGAGGACACACACAGGGGAGAAGCCCTATGTCTGCAGGGAG
# > v     TGTGGGCGGGGCTTTAGCTGGCAGTCAGTCCTCCTCAGTCACCAGAGGACACACACAGGGAAGAAGCCCTATGTCTGCAGGGAG

# znfs m-n are missing, znfs u-v added, in berg-2011-znf-sequences.tsv
```

- Allele sequences: copy/pasted from **Supplementary Figure 1B** to `berg-2011-allele-copy.txt`
  - Includes alleles A-E, L1-L27
  - Final file: `berg-2010-allele-sequences.tsv`

```
# remove extra columns, convert znf names to lowercase, tidy file, and sort alphabetically
awk '{print $1 "\t" tolower($3)}' intermediate-files/berg-2011-allele-copy.txt | sort -k1,1V > berg-2011-allele-sequences.tsv

# check that allelels are identical to those in berg-2010-allele-sequences.tsv
diff berg-2010-allele-sequences.tsv berg-2011-allele-sequences.tsv

# 29a30,32
# > L25   abcddecfuhfij
# > L26   abcddecfghfqj
# > L27   abcddvcfghfij

# alleles L25-L27 added in berg-2011-allele-sequences.tsv
```

# The *PRDM9* Variant Database
The *PRDM9* fandom is a bit of a mess:
- There is inconsistency in how alleles are named
- Different alleles have the same name
- The same name is used for different alleles
- Publications are claiming to have identified novel alleles that have actually been previously described
- Publications and NCBI Nucleotide entries have typos
- There is no database for all of the alleles uncovered to date

This is a collection of all known/described *PRDM9* allele and zinc finger variants. In addition, they have been renamed in a standardized manner: alleles are P### and zinc finger repeats are Z###. The allele zinc finger content and zinc finger sequence files contain a "map" to connect the standardized names to original names.

## Standardized lists
Variants observed in the population:
- [Allele sequences](standardized-lists/PRDM9-standardized-allele-sequences.tsv)
- [Zinc finger content of alleles](standardized-lists/PRDM9-standardized-allele-znf-content-map.tsv)
- [Zinc finger sequences](standardized-lists/PRDM9-standardized-znf-sequences-map.tsv)

Above lists plus variants observed only in sperm or as somatic mutations in blood (see [Jeffreys et al.](collect-known-alleles-publications.md#jeffreys-et-al-jan-2013))
- [Allele sequences blood/sperm](standardized-lists/population-and-sperm-somatic-PRDM9-standardized-allele-sequences.tsv)
- [Zinc finger content of alleles blood/sperm](standardized-lists/population-and-sperm-somatic-PRDM9-standardized-allele-znf-content-map.tsv)
- [Zinc finger sequences blood/sperm](standardized-lists/population-and-sperm-somatic-PRDM9-standardized-znf-sequences-map.tsv)


## Process
- [Step 1](collect-known-alleles-publications.md): Collect known *PRDM9* alleles and associated zinc fingers from publications
- [Steps 2-4](compile-publication-variants.md): Compile the publication variants
- [Step 5-8](collect-known-alleles-genbank.md): Collect known *PRDM9* allele DNA sequences from NCBI
- [Step 9](infer-ambiguous-alleles-publications.md): Infer allele sequencess mentioned in publications but not explicitly defined
- [Steps 10-11](finalize-allele-and-znf-lists.md): Finalize lists of *PRDM9* alleles and zinc finger repeats

#### TODO
- Add ref https://onlinelibrary.wiley.com/doi/10.1002/mgg3.56 for Genbank alleles L39-L46
- Add ref https://www.nature.com/articles/srep40031 for undefined allele L50
- Add note re: Alleva using allele N for Jeffreys Av:s:0053:M1S:A-A

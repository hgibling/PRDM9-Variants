# Infer _PRDM9_ allele sequencess mentioned in publications not explicitly defined in publications or GenBank
A handful of alleles were only decribed in publications in terms of amino acid sequences for the znfs, or as 3-4 amino acids per znf that bind DNA. The following sequences were not able to be inferred due to amino acid ambiguity:
- Parvanov alleles A3-A10, CH1-CH3
    - Full amino acid sequences for each znf
- Kong alleles Decode01-Decode07, YRI01-YRI09
    - `-1`, `3`, `6` amino acids for znfs `f` and `i` are the same

Allele K from Baudat, however was inferred and confirmed via last author de Massey.

## Step 9. Confirm Baudat figure matches GenBank allele seqeunces and infer znf content of allele K
Baudat depicted alleles A-F, H,I, and K in their publication figure without proper names for znfs. See if color blocks match znfs from GenBank sequences.

In addition, they did not submit a sequence for allele K. Determine the sequence from the other alleles.
```
# done in R

library(dplyr)
library(tidyr)
library(stringr)

baudat.allele.znf.seqs <- read.table("intermediate-files/baudat-2010-allele-znf-content-unnamed.tsv",
                                     header=F, col.names=c("Allele", "ZnfContent"))
pub.genbank.allele.znf.seqs <- read.table("intermediate-files/standardized-allele-znf-content-map-step8.tsv", 
                           header=T)

baudat.pub <- baudat.allele.znf.seqs %>%
  group_by(Allele) %>%
  # clean up blocks
  separate_rows(ZnfContent, sep="\\|") %>%
  mutate(ZnfPosition=row_number())

baudat.genbank <- pub.genbank.allele.znf.seqs %>%
  select(Baudat.2010, StandardZnfContent) %>%
  rename(Allele=Baudat.2010) %>%
  filter(!is.na(Allele)) %>%
  group_by(Allele) %>%
  separate_rows(StandardZnfContent, sep="_") %>%
  mutate(ZnfPosition=row_number())

znf.map <- baudat.genbank %>%
  full_join(baudat.pub %>%
              filter(Allele!="K")) %>%
  arrange(Allele, ZnfPosition) %>%
  ungroup() %>%
  select(-ZnfPosition, -Allele) %>%
  distinct() 

# check each znf has unique mapping to named znf
znf.map %>%
  group_by(ZnfContent) %>%
  summarize(StandardZnfContent=paste0(StandardZnfContent, collapse="_"))

# ZnfContent StandardZnfContent
# <chr>      <chr>             
#  1 --HR       z013              
#  2 --HS       z011              
#  3 --NS       z008              
#  4 -NHR       z006_z009         
#  5 -RVS       z014              
#  6 -RVT       z004              
#  7 -W--       z002              
#  8 -WVR       z012              
#  9 -WVS       z005              
# 10 -WVT       z003              
# 11 NTGR       z010              
# 12 NTOR       z001              
# 13 R-HR       z007

# one znf is ambiguous
# check if ambiguous znf is in allele K

baudat.pub %>%
  filter(Allele=="K", ZnfContent=="-NHR")

# Allele ZnfContent ZnfPosition
# <chr>  <chr>            <int>
# 1 K      -NHR                11
# 2 K      -NHR                14
# 3 K      -NHR                15

# in three places
# confirmed with de Massy (last author) that sequence is most likely z006 (f) at positions 11 and 14 and z009 (i) at position 15

# correct sequence for allele K and add to list of alleles
updated.allele.znf.content <- pub.genbank.allele.znf.seqs %>%
  # do join to check if sequence already in list
  full_join(baudat.pub %>%
              filter(Allele=="K") %>%
              left_join(znf.map) %>%
              filter(!(ZnfPosition==11 & StandardZnfContent=="z009"),
                !(ZnfPosition==14 & StandardZnfContent=="z009"),
                !(ZnfPosition==15 & StandardZnfContent=="z006")) %>%
              arrange(Allele, ZnfPosition) %>%
              group_by(Allele) %>%
              summarize(StandardZnfContent=paste0(StandardZnfContent, collapse="_")) %>%
              rename(Baudat.2010=Allele)) %>%
  mutate(StandardName=ifelse(is.na(StandardName),
                                   paste0("p", str_pad(row_number(), 3, pad="0")), 
                                   StandardName))

# update list of alleles
write.table(updated.allele.znf.content, 
            "intermediate-files/standardized-allele-znf-content-map-step9.tsv",
            row.names=F, quote=F, sep="\t")
```
# in-silico-trypsin
In-silico peptide digestion by trypsin

Peptide digestion can be done in one line of code in R using functions from the stringr package

```
str_subset (as.character(str_split_fixed(peptide,'(K|R)(?!P)',Inf)),'[:alpha:]')
```
Where *peptide* is the input amino acid sequence
The rule here is to cut whenever arginine (R) or lysine (K) are encountered, except when they are followed by proline (P)


## Flanking sequences of a peptide of interest
If we want to do mass spec targetting specific peptides, we will need to include flanking sequences in both directions until the first arginine or lysine. For example, 
- If we are interested in finding this peptide **ANVGAGRHGLYKPE**, 
- and the peptide is part of a bigger protein ALLAMKYTNQ**ANVGAGRHGLYKPE**QLQAIREFN
- Tryspin will digest it into 4 smaller peptides: ALLAM - YTNQ**ANVGAG - HGLYKPE**QLQAI - EFN  
- The peptide of interest is now in the two middle sequences, so we keep these, and delete the other sequences.

I wrote this function for a project involving parallel reaction monitoring (PRM), which is a type of mass spectrometry that is targeted to peptides of specific masses.

### Input dataframe with the following variables: 
- amino acid sequence of peptide of interest (*aa*), 
- longer protein sequence containing the peptide of interest (*longer_seq*)

### Output will be a dataframe with 2 additional variables: 
- *r2r* is aa plus right and left flanking sequences until R, K, or end of *longer_seq*
- *r2r_split* is *r2r* as digested by trypsin


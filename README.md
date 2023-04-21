# in-silico-trypsin
In-silico peptide digestion by trypsin

Peptide digestion is straighforward in R using functions from the stringr package

```
str_subset (as.character(str_split_fixed(peptide,'(K|R)(?!P)',Inf)),'[:alpha:]')
```
Where peptide is the input amino acid sequence
The rule here is to cut whenever arginine (R) or lysine (K) are encountered, except when they are followed by proline(P)


If you want to do mass spec targetting specific peptides, you will want to include flanking sequences in both directions until the first arginine or lysine. I wrote a function for that

For example, if I am interested in finding this peptide ANVGAGRHGLYKPE, and this peptide is part of a bigger protein ALLAMKYTNQANVGAGRHGLYKPEQLQAIREFN
After trypsin digestion, you get this: ALLAM YTNQANVGAG HGLYKPEQLQAI EFN  

The peptide of interest is now in the two middle sequences, so I'd like to keep these, and delete the other sequences. I wrote a function for this task

```
# input dataframe with the following columns: 
# amino acid sequence of peptide of interest (aa), 
# longer protein sequence containing the peptide of interest (longer_seq), and
# position of first amino acid of peptide in the longer protein (aa_position)

trypsin= function () {
  df= data.frame(peptide=character(), epitope=character(), r2r=character(), r2r_split=character())
  t= mutate (t, len=nchar(aa))
  
  for (i in 1:nrow(t)) {
    
    peptide=t$seq[i]
    aa_st=t$Sub.peptide.Position[i]
    aa_nd= aa_st+ t$len[i]-1
    
    cleav=as.data.frame(str_locate_all(peptide,'(K|R)(?!P)'))$start
    st= max(cleav[which(cleav <= aa_st)],0)
    
    nd=min(cleav[which(cleav >aa_nd)], nchar(peptide) ) 
    
    r2r=str_sub(peptide,st, nd)
    #print(r2r)
    
    r2r_split= str_subset (as.character(str_split_fixed(r2r,'(K|R)(?!P)',Inf)),'[:alpha:]')
    r2r_split=str_c (r2r_split, collapse='-')
    #print(r2r_split)
    
    df=add_row(df, peptide=peptide, epitope=t$aa[i], r2r=r2r, r2r_split=r2r_split)
  }
  return (df)
}


```

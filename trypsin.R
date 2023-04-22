####################### 
#Function to pad a smaller peptide of interest by left and right flanking sequences till R, K, or end of longer peptide
# Written by Rayan Najjar
#######################

trypsin_flank= function (aa_df) {
  df= data.frame(longer_seq=character(), aa=character(), r2r=character(), r2r_split=character())
  
  for (i in 1:nrow(aa_df)) {
    peptide=aa_df$longer_seq[i]
    aa= aa_df$aa[i]
    len=nchar(aa)
    aa_st= str_locate(peptide,aa)[,1]
    aa_nd= aa_st+len-1
    
    cleav=as.data.frame(str_locate_all(peptide,'(K|R)(?!P)'))$start
    st= max(cleav[which(cleav <= aa_st)],0)
    nd=min(cleav[which(cleav >aa_nd)], nchar(peptide) ) 
    r2r=str_sub(peptide,st,nd)
    
    r2r_split= str_subset (as.character(str_split_fixed(r2r,'(K|R)(?!P)',Inf)),'[:alpha:]')
    r2r_split=str_c (r2r_split, collapse='-')
    
    df=add_row(df, longer_seq=peptide, aa=aa, r2r=r2r, r2r_split=r2r_split)
  }
  return (df)
}
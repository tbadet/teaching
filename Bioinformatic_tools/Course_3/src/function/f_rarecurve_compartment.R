# function for make rarecurve for each compartment
rarecurve_compartment<-function(compartment_i,type='16S'){
  cat(compartment_i, "\n")
  
  if(type=="16S"){
  compartment_l=df_16S_count %>% 
    rownames_to_column("ASV") %>% 
    pivot_longer(-ASV, names_to = "sample_name", values_to = "count") %>% 
    inner_join(.,df_16S_metadata, by="sample_name") %>% 
    filter(compartment==compartment_i) %>% 
    head(1) %>% 
    pull(compartment_l) %>% 
    as.character()
  
  df_tmp<-df_16S_count %>%t() %>% 
    as.data.frame() %>% 
    rownames_to_column("sample_name") %>% 
    inner_join(df_16S_metadata %>% 
                 dplyr::select(sample_name,compartment),
               ., by="sample_name") %>% 
    filter(compartment==compartment_i) %>% 
    dplyr::select(-compartment) %>% 
    column_to_rownames("sample_name")
  
  rarecurve_data=vegan::rarecurve(df_tmp, step = 100, cex = 0.75, las = 1) 
  
  min_n_seqs=min(
    df_16S_count %>% 
      rownames_to_column("ASV") %>% 
      pivot_longer(-ASV, names_to = "sample_name", values_to = "count") %>% 
      inner_join(.,df_16S_metadata, by="sample_name") %>% 
      filter(compartment==compartment_i) %>% 
      dplyr::group_by(sample_name) %>% 
      dplyr::summarise(total=sum(count)) %>% 
      pull(total)
  )
  
  xintercept_i=16000
  
  }else if(type=="ITS"){
    compartment_l=df_ITS_count %>% 
      rownames_to_column("ASV") %>% 
      pivot_longer(-ASV, names_to = "sample_name", values_to = "count") %>% 
      inner_join(.,df_ITS_metadata, by="sample_name") %>% 
      filter(compartment==compartment_i) %>% 
      head(1) %>% 
      pull(compartment_l) %>% 
      as.character()
    
    df_tmp<-df_ITS_count %>%t() %>% 
      as.data.frame() %>% 
      rownames_to_column("sample_name") %>% 
      inner_join(df_ITS_metadata %>% 
                   dplyr::select(sample_name,compartment),
                 ., by="sample_name") %>% 
      filter(compartment==compartment_i) %>% 
      dplyr::select(-compartment) %>% 
      column_to_rownames("sample_name")
  
    rarecurve_data=vegan::rarecurve(df_tmp, step = 100, cex = 0.75, las = 1) 
    
    min_n_seqs=min(
      df_ITS_count %>% 
        rownames_to_column("ASV") %>% 
        pivot_longer(-ASV, names_to = "sample_name", values_to = "count") %>% 
        inner_join(.,df_ITS_metadata, by="sample_name") %>% 
        filter(compartment==compartment_i) %>% 
        dplyr::group_by(sample_name) %>% 
        dplyr::summarise(total=sum(count)) %>% 
        pull(total)
    )
    
    xintercept_i=83000
    }
  
  alpha_rarefaction<-map_dfr(rarecurve_data,bind_rows) %>% 
    bind_cols(Group=rownames(df_tmp),.) %>% 
    pivot_longer(-Group) %>% 
    drop_na() %>% 
    mutate(n_seqs=as.numeric(str_replace(name,"N",""))) %>% 
    dplyr::select(-name) %>%
    ggplot(aes(x=n_seqs,y=value,group=Group))+
    geom_vline (xintercept=min_n_seqs,color="red")+
    geom_line()+
    theme_classic()+
    geom_vline (xintercept=min_n_seqs,color="red")+
    geom_vline (xintercept=xintercept_i,color="green")+
    #geom_vline (xintercept=1800,color="green")+
    geom_line()+
    theme_classic()+
    labs(title = paste0("Compartment :", compartment_l),
         subtitle = paste0("min sample: ",min_n_seqs)
    )+
    ylab("Taxonomic richness")+ 
    xlab("No. reads")
  
  return(alpha_rarefaction)
}

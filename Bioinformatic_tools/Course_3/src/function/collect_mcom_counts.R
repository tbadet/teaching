# Function to collect bacteria counts and ASV per sample

collect_mcom_counts <- function(storage=data.frame(matrix(data=NA,nrow=0, ncol=0)),df_i, new_column_name,message=F) {
  df_resum=df_i %>% 
    t() %>%
    as.data.frame()%>%
    rownames_to_column("sample_name")%>%
    pivot_longer(-sample_name)%>%
    dplyr::group_by(sample_name)%>%
    dplyr::summarise(Sum=sum(value))%>%
    drop_na() %>% 
    arrange(Sum)
  
  # # Initialize an empty dataframe to store the results if nothing before
   if(length(storage)==0) {
     if(message==T){
       cat(" First input or error \n")
      }else{}
     storage <-df_resum
     colnames(storage)<-c("sample_name",new_column_name)
  }else{
      colnames(df_resum)<-c("sample_name",new_column_name)
      storage<-full_join(storage,df_resum,by="sample_name") 
    }
    
    # Returns of some informations
  if(message==T){
    cat(" Table added \n",
        " Total number of counts : ",sum(df_resum[,2]),
        "\n Number of samples :" ,length(df_i),
        "\n Total number of ASV :" ,nrow(df_i),
        "\n Sample with the least number of counts is : ", df_resum[1,1]%>%pull(), 
        "With : ",df_resum[1,2]%>%pull(), "counts",
        " \n Number of data tables in the backup : ",length(storage)-1,
    "\n"
    )
    print(storage%>% arrange(!!sym(new_column_name)))
    }else{}
    write.csv(file = here::here("data/mcom/output/loos_by_step.csv"),x=storage%>% arrange(!!sym(new_column_name)))  
    return(storage%>% arrange(!!sym(new_column_name)))
   
}

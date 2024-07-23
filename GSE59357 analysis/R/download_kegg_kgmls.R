download_KEGG_pathways_with_searchKey <- function(orgID = "hsa", searchKey = "signaling", 
                                                  outPath = "data/KGMLs/", saveToGMT = T){
  require(KEGGREST)
  require(qdapTools)
  
  # orgID = "hsa"
  # searchKey = "signaling"
  # outPath = "data/KGMLs/"
  dir.create(outPath)
  
  # get list of pathways of an organism
  allPaths = keggList("pathway", orgID)
  allPaths.filt = allPaths %>% 
    .[base::grep(x = ., pattern =  searchKey)] 
  
  gmt_content <- NULL
  for(x in 1: length(allPaths.filt)){
    # get and Save XMLs
    pathID = allPaths.filt[x] %>% names()
    fileContent = keggGet(pathID,"kgml")
    targetFileName = paste0(outPath, pathID, ".xml")
    writeLines(fileContent, targetFileName)
    
    if(saveToGMT){
      # get genes
      path_data = keggGet(pathID)[[1]]
      path_name = path_data$NAME %>% 
        gsub(pattern = " - Homo sapiens (human)", replacement = "", fixed = T) %>%
          gsub(pattern = " ", replacement = "_")
      gene_vector = path_data$GENE
      gene_df <- data.frame(
        Gene.ID = gene_vector[seq(1, length(gene_vector), by = 2)],
        Gene.Symbol = gene_vector[seq(2, length(gene_vector), by = 2)] %>% 
          base::gsub(pattern = ";.*", replacement="")
      ) 
      path_GS = gene_df$Gene.Symbol %>% unique() %>% paste0(collapse = "\t")
      gmt_content <- rbind(gmt_content, c(path_name, pathID, path_GS))
    }
  }
  if(saveToGMT){
    gmt_content %>%
      fwrite(file = paste0(outPath, "KEGG_selected_pathways.csv"), 
             col.names = F, quote = F, sep = "\t", row.names = F)
  }
}
parse_XMLfile<-function(pathway_id, species, database_dir=getwd()) {
  #get pathway gene and their location in the pic
  #except three global maps
  if (pathway_id=="01100" | pathway_id=="01110" | pathway_id=="01120") {
    print(paste("Skip global maps:",pathway_id,sep=""))
    return(NULL)
  }
  require(XML)
  inter1<-xmlTreeParse(paste(database_dir,"/",species,pathway_id,".xml",sep=""),useInternalNodes=TRUE)
  pathName <- getNodeSet(inter1,"//pathway")[[1]] %>% xmlGetAttr(.,  "title")
  inter2<-getNodeSet(inter1,"//entry")
  inter3<-lapply(inter2,  function(xxx) xmlGetAttr(xxx,  "name"))
  inter4<-lapply(inter2,  function(xxx) xmlGetAttr(xxx,  "type"))
  inter5<-sapply(inter2,  function(xxx) getNodeSet(xxx,".//graphics"))
  inter_graphic_type<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "type"))
  inter6<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "x"))
  inter7<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "y"))
  inter8<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "width"))
  inter9<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "height"))
  inter10<-sapply(inter5,  function(xxx) xmlGetAttr(xxx,  "name"))
  result<-NULL
  for (i in 1:length(inter4)) {
    if ((inter4[[i]]=="gene" | inter4[[i]]=="compound") & inter_graphic_type[[i]]!="line") {
      temp<-strsplit(inter3[[i]]," ")[[1]]
      #			temp<-gsub("cpd:","",temp)
      name<-strsplit(inter10[[i]],",")[[1]][1]
      name<-gsub('\\.\\.\\.',"",name)
      for (n in 1:length(temp)) {
        result<-rbind(result,c(pathName, temp[n],inter6[[i]],inter7[[i]],inter8[[i]],inter9[[i]],name))
      }
    }
  }
  # result[,1]<-gsub(paste(species,":",sep=""),'',result[,1])
  return(result)
}

KEGGPathwayXMLs_to_GMT <- function(speciesID = "hsa", xmlPaths = "data/KGMLs/",
                                   outPath = "data/"){
  library(data.table)
  
  pathway_ids = list.files(path = xmlPaths) %>% 
    base::gsub(pattern = ".xml", replacement = "") %>% 
    base::gsub(pattern = speciesID, replacement = "")
  
  result <- NULL
  for(i in 1: length(pathway_ids)){
    pathId = pathway_ids[i]
    path.df = parse_XMLfile(pathway_id = pathId, 
                            species = speciesID, 
                            database_dir = xmlPaths) %>% as.data.frame()
    pathName = path.df[1,1]
    path.gs = path.df %>%  dplyr::select(V7) %>% 
      unlist(use.names = F) %>% unique() %>% paste0(collapse = "\t")
    result <- rbind(result, c(pathName, "dummyText",path.gs))
  }
  result %>%
    fwrite(file = paste0(outPath, "XMLs_to_KEGG_pathways.csv"), col.names = F,
           quote = F,
           sep = "\t", row.names = F)
}

#YanTing hdac 20200309 functions 



readInTidy = function(raw_data_filename,
                      seq = 3,
                      protAcc = 10,
                      des = 13,
                      sp1ctrl = 16,
                      sp2ctrl = 17,
                      sp1low = 18,
                      sp2low = 19,
                      sp1high = 20,
                      sp2high = 21,
                      output_name )
{
  raw = fread(raw_data_filename, stringsAsFactors = F)
  
  sel_col = c(protAcc, seq, des, sp1ctrl, sp2ctrl, sp1low, sp2low, sp1high, sp2high)
  
  sel_data = raw[, ..sel_col]
  
  colnames(sel_data) = c("protein", "seq", "description", "b1_ctrl","b2_ctrl","b1_low","b2_low","b1_high","b2_high")
  
  
  
  
  fil_data = sel_data%>%
    na.omit() 
  
 ### separate proteins 
  
  
  
  
  multi_name = grep(";", fil_data$protein)
  single_name = grep(";", fil_data$protein, invert = T)
  
  
  ### ok here I need to separate the description part also 
  
  
  multi_df = rbindlist(lapply(1:length(multi_name), function(x) {
    
    this_row = multi_name[x]
    
    this_prot = fil_data$protein[this_row]
    this_des = fil_data$description[this_row]
    
    names = unlist(strsplit(this_prot, split = "; "))
    dess = unlist(strsplit(this_des, split = ";"))
    
    
    this_df = fil_data[rep(this_row,length(names)),]
    
    this_df$protein = names  
    this_df$description = dess
    
    return(this_df)  
    
    
  }))
  
  
  combine_df = rbind(multi_df,fil_data[single_name,])%>%
    dplyr::arrange(protein)
  
  
  #### get name fo the gene name 
  
  all_gn = unlist(lapply(1:nrow(combine_df), function(x) {
    
    this_des = combine_df$description[x]
    
    get_des = unlist(strsplit(this_des,split = "GN="))[2]
    
    get_gn = unlist(strsplit(get_des, split = " ", fixed = T))[1]
    
    return(get_gn)
    
  }))
  
  combine_df$description = all_gn
  
  
  pd_combine_df = combine_df%>%
    dplyr::mutate(protein = paste(protein, description, sep = "_"))%>%
    dplyr::select(-description)

  write.table(pd_combine_df, output_name, 
              quote = F, row.names = F, sep = "\t")
  
  
  return(combine_df)
  
  
}


### ok 
### after this I need to make sure the gn col is not lost in the following functions 









##### 


#### I think reliability if need to be evaluated, should be done before agg and before normalization 
npos = function(x)
{
  return(length(which(x>0)))
}


nneg = function(x)
{
  return(length(which(x<=0)))
}


testReliability_afterFC = function(peptide_fc,
                           output_name)
{
#   peptide_fc = fc37
#   output_name = reliability_outputname1
#   
  
  rel_data = peptide_fc%>%
    dplyr::mutate(drl = fc_l_b1* fc_l_b2, drh = fc_h_b1*fc_h_b2)

  
  #### return some lables 
  rel_dr = rel_data%>%
    dplyr::group_by(protein, description)%>%
    dplyr::summarise(npep = n(),
                     same_drl = npos(drl), 
                     diff_drl = nneg(drl),
                     same_drh = npos(drh),
                     diff_drh = nneg(drh))%>%
    dplyr::ungroup()%>%
    dplyr::mutate(rel_low = case_when(diff_drl == 0 ~ T,
                                      diff_drl > 0 ~ F),
                  rel_high = case_when(diff_drh == 0 ~ T,
                                       diff_drh >0 ~ F))
  write.table(rel_dr, output_name, 
              quote = F, row.names = F, sep = "\t")
  
  return(rel_dr)
  
  
}






#### mapDIA process 

#### make it compatible with windows 

mapDIApro = function(working_dir,
                     mapDIA_dir,
                     userSystem,
                     default_input_filename,
                     change_input_filename,
                     peptide_df_filename,
                     analysis_output_name)


{
  
  default_input = readLines(default_input_filename)
  
 
  
  set_SDF = 3
  
  default_input[2] = paste0("FILE= ", peptide_df_filename)
  default_input[17] = paste0("SDF= ", set_SDF)
  
  
  
  
  write.table(default_input, change_input_filename,
              col.names = F, row.names = F, quote = F)
  
  
  
  
  ########## enter the MAPDIA program
  
  
  
 # setwd(workingDir)
  
  linux_command =  paste0( mapDIA_dir, "mapDIA ", change_input_filename)
  #windows_command =  paste0( mapDIA_dir, "mapDIA_win64.exe ", change_input_filename)
  
  
  if(userSystem == "Windows")
  {
    cat("run mapDIA now****************************************","\n")
    
    
    wcmd = paste0(mapDIA_dir, "mapDIA_win64.exe")
    warg = change_input_filename
    
    system2(command = wcmd,
            args = warg)
    

  }else{
    
    cat("run mapDIA now","\n")
    
    system(linux_command)
    
  }

  
  
  this_dir = getwd()
  
  old_output_name = paste0(this_dir, "/analysis_output.txt")
  new_output_name = paste0(working_dir, analysis_output_name)
  
  linux_mv_command = paste0("mv ", old_output_name," ",new_output_name)
  #windows_mv_command = paste0("move ", old_output_name," ",new_output_name)
  
  if(userSystem == "Windows")
  {
    
    cat("$$$$$$$$$$$$$$$$$$$$ ", this_dir,"\n")
    
    w_old_output_name = paste0(this_dir, "/analysis_output.txt")
  #  wo = gsub("/", "\\\\",
   #           w_old_output_name )
    
    w_new_output_name = paste0(working_dir, analysis_output_name)
  #  wn = gsub("/", "\\\\",
   #           w_new_output_name )
    
   windows_mv_command = paste0("move ", w_old_output_name," ",w_new_output_name)
  #  windows_mv_command = paste0("move ", wo," ",wn)
    
    shell(windows_mv_command)
    
    
  }else{
    system(linux_mv_command)
    
  }
}







### another way of fc calculation 
### mean the peptide after normalization rather than aggregate them 

calFoldChange1 = function(peptide_df,
                          output_name)
 {

  
  data_pep = log2(peptide_df[,-c(1,2,3)])
  
  ### normalize 
  

  each_med = apply(data_pep, 2, median)
  
  med_each_med = median(each_med)
  
  norm_pep_log2data = sweep(data_pep, 2, each_med) + med_each_med
  
  ### get fold change
  #### how to keep tracking the gene name of each protein 
  
  norm_log_peptide_df = cbind(peptide_df[,c(1,2,3)], norm_pep_log2data)%>%
    dplyr::group_by(protein, description)%>%
    dplyr::summarise_at(vars(b1_ctrl:b2_high), mean)%>%
    dplyr::ungroup()%>%
    dplyr::mutate(fc_l_b1 = b1_low - b1_ctrl, fc_l_b2 = b2_low -b2_ctrl)%>%
    dplyr::mutate(fc_h_b1 = b1_high - b1_ctrl, fc_h_b2 = b2_high -b2_ctrl)%>%
    dplyr::mutate(fc_l = (fc_l_b1+fc_l_b2)/2,
                  fc_h = (fc_h_b1+fc_h_b2)/2)%>%
    dplyr::mutate(fc_l_diff = fc_l_b2 - fc_l_b1,
                  fc_h_diff = fc_h_b2 - fc_h_b1)
    
  
  write.table(norm_log_peptide_df, 
              output_name,
              quote = F, row.names = F, sep = "\t")
  
    
  return(norm_log_peptide_df)
    
  
  
  
}








joinScatterOutput = function(pepFC1,
                       pepFC2,
                       fc_upper = 1.3,
                       fc_lower = 0.7,
                       mapDIA_output_filename1,
                       mapDIA_output_filename2,
                       mapDIA_sigCutoff,
                       prot_reliable_filename1,
                       prot_reliable_filename2,
                       low_output_name,
                       high_output_name,
                       low_pdf_name,
                       high_pdf_name)

{
  
  join_pep = pepFC1%>%
    dplyr::left_join(pepFC2, by = c("protein","description"))%>%
    na.omit()
  
    # 
  join_pep_low = join_pep%>%
    dplyr::select(protein,description, fc_l.x, fc_l.y, fc_l_diff.x, fc_l_diff.y)#%>%
  #dplyr::mutate(sigFlag = "noneSig", relFlag = "noneRel")

  
  join_pep_high = join_pep%>%
    dplyr::select(protein, description, fc_h.x, fc_h.y, fc_h_diff.x, fc_h_diff.y)#%>%
  # dplyr::mutate(sigFlag = "noneSig", relFlag = "noneRel")
  
  
  
  
  
  ##### get other data parsed 
  if(file.exists(mapDIA_output_filename1) & file.exists(mapDIA_output_filename2))
  {
    mapDIA1 = fread(mapDIA_output_filename1, stringsAsFactors = F)
    mapDIA2 = fread(mapDIA_output_filename2, stringsAsFactors = F)
    
    
    join_fdr = mapDIA1%>%
      dplyr::left_join(mapDIA2, by = c("Protein", "Label2"))%>%
      dplyr::select(Protein, Label2, FDR.x, FDR.y)
    
    prot_gn_df = rbindlist(lapply(1:nrow(join_fdr), function(x)
      {
      this_line  = unlist(strsplit(join_fdr$Protein[x],split = "_"))
      
      this_prot = this_line[1]
      this_des = this_line[2]
      
      this_df = data.frame(protein = this_prot, description = this_des,
                           stringsAsFactors = F)
      
      return(this_df)
      
      
    }))
    
    join_fdr_des = cbind(prot_gn_df, join_fdr[,-1])
    
    
    join_fdr_low = join_fdr_des %>%
      dplyr::filter(Label2 == "low/control")%>%
      na.omit()%>%
      dplyr::mutate(sigFlag = case_when(FDR.x<=mapDIA_sigCutoff & FDR.y <= mapDIA_sigCutoff ~ "bothSig",
                                        FDR.x<=mapDIA_sigCutoff & FDR.y > mapDIA_sigCutoff ~ "firstSig",
                                        FDR.x>mapDIA_sigCutoff & FDR.y <= mapDIA_sigCutoff ~ "secondSig",
                                        FDR.x>mapDIA_sigCutoff & FDR.y > mapDIA_sigCutoff ~ "noneSig"))
    
    join_pep_low = join_pep_low%>%
      dplyr::left_join(join_fdr_low, by = c("protein", "description"))%>%
      dplyr::select(protein, description,fc_l.x, fc_l.y, FDR.x, FDR.y, sigFlag,fc_l_diff.x, fc_l_diff.y )
    
    
    join_fdr_high = join_fdr_des %>%
      dplyr::filter(Label2 == "high/control")%>%
      na.omit()%>%
      dplyr::mutate(sigFlag = case_when(FDR.x<=mapDIA_sigCutoff & FDR.y <= mapDIA_sigCutoff ~ "bothSig",
                                        FDR.x<=mapDIA_sigCutoff & FDR.y > mapDIA_sigCutoff ~ "firstSig",
                                        FDR.x>mapDIA_sigCutoff & FDR.y <= mapDIA_sigCutoff ~ "secondSig",
                                        FDR.x>mapDIA_sigCutoff & FDR.y > mapDIA_sigCutoff ~ "noneSig"))
    
    join_pep_high = join_pep_high%>%
      dplyr::left_join(join_fdr_high, by = c("protein","description"))%>%
      dplyr::select(protein,description, fc_h.x, fc_h.y, FDR.x, FDR.y, sigFlag, fc_h_diff.x, fc_h_diff.y)
    
    
    
  }
  
  ### ok I decide to have use both uniprot accession and gene name as the identifier so that I don't mess up things 
  
  
  
  if(file.exists(prot_reliable_filename1) & file.exists(prot_reliable_filename2))
  {
    
    protR1 = fread(prot_reliable_filename1, stringsAsFactors = F)
    protR2 = fread(prot_reliable_filename2, stringsAsFactors = F)
    
    join_rel = protR1%>%
      dplyr::left_join(protR2, by = c("protein", "description"))%>%
      dplyr::select(protein,description,rel_low.x,rel_low.y, rel_high.x, rel_high.y)
    
    
    join_rel_low = join_rel%>%
      dplyr::select(protein, description, rel_low.x, rel_low.y)%>%
      na.omit()%>%
      dplyr::mutate(relFlag = case_when(rel_low.x  == T & rel_low.y  == T ~ "bothRel",
                                        rel_low.x  == T & rel_low.y  == F ~ "firstRel",
                                        rel_low.x  == F & rel_low.y  == T ~ "secondRel",
                                        rel_low.x  == F & rel_low.y  == F ~ "noneRel"))
    
    join_pep_low = join_pep_low%>%
      dplyr::left_join(join_rel_low, by = c("protein","description"))
    
    
    join_rel_high = join_rel%>%
      dplyr::select(protein,description, rel_high.x, rel_high.y)%>%
      na.omit()%>%
      dplyr::mutate(relFlag = case_when(rel_high.x  == T & rel_high.y  == T ~ "bothRel",
                                        rel_high.x  == T & rel_high.y  == F ~ "firstRel",
                                        rel_high.x  == F & rel_high.y  == T ~ "secondRel",
                                        rel_high.x  == F & rel_high.y  == F ~ "noneRel"))
    
    
    
    join_pep_high = join_pep_high%>%
      dplyr::left_join(join_rel_high, by =c("protein", "description"))
    
  }
    # 
  write.table(join_pep_low, low_output_name, quote = F, row.names = F, sep = "\t")
  write.table(join_pep_high, high_output_name, quote = F, row.names = F, sep = "\t")
  
  #### ok there is some thing with the ploting of names 
 
  
  prot_des_low = join_pep_low%>%
    dplyr::select(protein, description)%>%
    unique()%>%
    dplyr::arrange(description,protein)%>%
    dplyr::mutate(number = 1)%>%
    dplyr::group_by(description)%>%
    dplyr::mutate(ticker = cumsum(number))%>%
    dplyr::mutate(symbol = case_when(ticker>1 ~ paste0(description,"-",ticker),
                                     ticker ==1 ~ description))%>%
    dplyr::select(protein, description, symbol)
  
  
  ### map back 
  symbol_join_pep_low = join_pep_low%>%
    dplyr::left_join(prot_des_low, by = c("protein","description"))
  
  
   
  prot_des_high = join_pep_high%>%
    dplyr::select(protein, description)%>%
    unique()%>%
    dplyr::arrange(description,protein)%>%
    dplyr::mutate(number = 1)%>%
    dplyr::group_by(description)%>%
    dplyr::mutate(ticker = cumsum(number))%>%
    dplyr::mutate(symbol = case_when(ticker>1 ~ paste0(description,"-",ticker),
                                     ticker ==1 ~ description))%>%
    dplyr::select(protein, description, symbol)

  
  ### map back 
  symbol_join_pep_high = join_pep_high%>%
    dplyr::left_join(prot_des_high, by = c("protein","description"))
  
  
  ###  plot scatter 
  ### have the ones marked with gene name 
  
  ### low/control
  pdf(low_pdf_name, useDingbats = F)
  plot(symbol_join_pep_low$fc_l.x, symbol_join_pep_low$fc_l.y,
       type = "p", pch = 16, cex = 0.3,
       xlim = c(-1.5,1.5),
       ylim = c(-1.5,1.5),
       xlab = "log2FC_37C",
       ylab = "log2FC_52C",
       main = "low/ctrl fold change")
  abline(h= log2(fc_upper), lty = 3)
  abline(h= log2(fc_lower), lty = 3)
  abline(v= log2(fc_upper), lty = 3)
  abline(v= log2(fc_lower), lty = 3)
  abline(h = 0, lty = 2)
  abline(v = 0, lty = 2)  
  
  ### add marks to the fc ones 
  
  join_low_hit = symbol_join_pep_low%>%
       dplyr::filter(fc_l.x>log2(fc_upper)| fc_l.x<log2(fc_lower) | fc_l.y>log2(fc_upper)| fc_l.y<log2(fc_lower))

  

  text(jitter(join_low_hit$fc_l.x),jitter(join_low_hit$fc_l.y), 
       labels = join_low_hit$symbol, pos = 1, cex = 0.25)
  
  
  

  if("relFlag"%in%colnames(join_pep_low))
  {

    plot_rel = symbol_join_pep_low%>%
      dplyr::filter(relFlag == "bothRel"| relFlag == "firstRel" | relFlag == "secondRel")
    rel_pch = rep(16, nrow(plot_rel))
    rel_pch[which(plot_rel$relFlag == "bothRel")] = 2
    rel_pch[which(plot_rel$relFlag == "firstRel")] = 0
    rel_pch[which(plot_rel$relFlag == "secondRel")] = 0


    lines(plot_rel$fc_l.x, plot_rel$fc_l.y,
          type = "p", pch = rel_pch, cex = 0.4)



  }


  if("sigFlag"%in%colnames(join_pep_low))
  {
    plot_sig = symbol_join_pep_low%>%
      dplyr::filter(sigFlag == "bothSig" | sigFlag == "firstSig" | sigFlag == "secondSig")

    sig_col = rep("black", nrow(plot_sig))
    sig_col[which(plot_sig$sigFlag == "bothSig")] = "red"
    sig_col[which(plot_sig$sigFlag == "firstSig")] = "green"
    sig_col[which(plot_sig$sigFlag == "secondSig")] = "green"


    lines(plot_sig$fc_l.x, plot_sig$fc_l.y,
          type = "p", pch = 16, cex = 0.5,
          col = sig_col)



  }


  legend(1, 1.5, legend=c("bothSig", "oneSig","bothRel","oneRel"),
         col=c("red", "green","black","black"), pch = c(16,16,2,0), cex=0.8)
  
  
  
  dev.off()
  
  
  #####
  #####
  #####
  
  
  
  
  pdf(high_pdf_name, useDingbats = F)
  plot(symbol_join_pep_high$fc_h.x, symbol_join_pep_high$fc_h.y,
       type = "p", pch = 16, cex = 0.3,
       xlim = c(-1.5,1.5),
       ylim = c(-1.5,1.5),
       xlab = "log2FC_37C",
       ylab = "log2FC_52C",
       main = "high/ctrl fold change")
  abline(h= log2(fc_upper), lty = 3)
  abline(h= log2(fc_lower), lty = 3)
  abline(v= log2(fc_upper), lty = 3)
  abline(v= log2(fc_lower), lty = 3)
  abline(h = 0, lty = 2)
  abline(v = 0, lty = 2)  
  
  
  join_high_hit = symbol_join_pep_high%>%
    dplyr::filter(fc_h.x>log2(fc_upper)| fc_h.x<log2(fc_lower) | fc_h.y>log2(fc_upper)| fc_h.y<log2(fc_lower))
  
  
  text(jitter(join_high_hit$fc_h.x, factor = 5),jitter(join_high_hit$fc_h.y, factor = 5), labels = join_high_hit$symbol, pos = 1, cex = 0.25)
  
  
  

  
  if("relFlag"%in%colnames(join_pep_high))
  {
    
    plot_rel = symbol_join_pep_high%>%
      dplyr::filter(relFlag == "bothRel"| relFlag == "firstRel" | relFlag == "secondRel")
    rel_pch = rep(16, nrow(plot_rel))
    rel_pch[which(plot_rel$relFlag == "bothRel")] = 2
    rel_pch[which(plot_rel$relFlag == "firstRel")] = 0
    rel_pch[which(plot_rel$relFlag == "secondRel")] = 0
    
    
    lines(plot_rel$fc_h.x, plot_rel$fc_h.y,
          type = "p", pch = rel_pch, cex = 0.4)
    
    
    
  }
  
  
  if("sigFlag"%in%colnames(join_pep_high))
  {
    plot_sig = symbol_join_pep_high%>%
      dplyr::filter(sigFlag == "bothSig" | sigFlag == "firstSig" | sigFlag == "secondSig")
    
    sig_col = rep("black", nrow(plot_sig))
    sig_col[which(plot_sig$sigFlag == "bothSig")] = "red"
    sig_col[which(plot_sig$sigFlag == "firstSig")] = "green"
    sig_col[which(plot_sig$sigFlag == "secondSig")] = "green"
    
    
    lines(plot_sig$fc_h.x, plot_sig$fc_h.y,
          type = "p", pch = 16, cex = 0.5,
          col = sig_col)
    
    
    
  }
  
  legend(1, 1.5, legend=c("bothSig", "oneSig","bothRel","oneRel"),
         col=c("red", "green","black","black"), pch = c(16,16,2,0), cex=0.8)
  
  
  
  
  dev.off()
  
  
  
  
  return(join_pep)
  
  
  
}


























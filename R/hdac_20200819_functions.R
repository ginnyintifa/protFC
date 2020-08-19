

### if there are 3 replicates 


### I think I should make it applicable to any number of reps
# 
# library(dplyr)
# library(data.table)
# library(magrittr)
# # 
# rawData_inputfilename1 = "/data/ginny/YanTingProtein/20200813_TF_MEHP37_F_PeptideGroups.txt"
# rawData_inputfilename2 = "/data/ginny/YanTingProtein/20200813_TF_MEHP53_F_PeptideGroups.txt"
# inputSpec_filename =  "/data/ginny/YanTingProtein/20200813_TF_spCol_positions.txt"
# 
# seqCol = 3
# protAccCol = 10
# desCol = 13
# 
# sp1ctrlCol = 16
# sp2ctrlCol = 17
# sp1lowCol = 18
# sp2lowCol = 19
# sp1highCol = 20
# sp2highCol = 21
# workingDir = "/data/ginny/YanTingProtein/foldChange0819/"
# mapDIA_flag = F
# mapDIADir = "/data/ginny/mapDIA-master/"
# userSystem = "Linux"
# defaultInput_parameter = "/data/ginny/YanTingProtein/YTinput"
# cutoff_uFc = 1.3
# cutoff_lFc = 0.7
# cutoff_mapDIA = 0.05
# outputTagS1 = "lowVSctrl"
# outputTagS2 = "highVSctrl"



readInTidy_r = function(raw_data_filename,
                        data_col_filename,
                                        output_name )
{
  
  # raw_data_filename = rawData_inputfilename1
  # data_col_filename = inputSpec_filename
  #  output_name = cleanData_outputname1
  # 
  # 
  raw = fread(raw_data_filename, stringsAsFactors = F, data.table = F)
  spec = fread(data_col_filename, stringsAsFactors = F, data.table = F)
  
 cond_names = grep("cond_", spec$feature, value = T)
  
  #### suppose there are control low high and the number of reps is elastic 
  
  ctrl_cols = spec$position[grep("cond_ctrl", spec$feature)]
  low_cols = spec$position[grep("cond_low", spec$feature)]
  hight_cols = spec$position[grep("cond_high", spec$feature)]
  seq = spec$position[which(spec$feature == "seq")]
  acc = spec$position[which(spec$feature == "acc")]
  des = spec$position[which(spec$feature == "des")]
  
  sel_col = c(acc, seq, des, ctrl_cols, low_cols, hight_cols)
  
  sel_data = raw[,sel_col]
  
  colnames(sel_data) = c("protein", "seq", "description", cond_names )
  
  

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



# 
# npos = function(x)
# {
#   return(length(which(x>0)))
# }
# 
# 
# nneg = function(x)
# {
#   return(length(which(x<0)))
# }

nzero = function(x)
{
  return(length(which(x==0)))
}




test_equal_sign = function(x)
{
  x_sign = sign(x)
  x_sign_equal <- var(x_sign)==0
  return(x_sign_equal)
  
}




calFoldChange1_r = function(peptide_df,
                          output_name)
{
  # 
  # peptide_df = pep37
  # output_name = foldChange_outputname1
  # 
  
  data_pep = log2(peptide_df[,-c(1,2,3)])
  
  ### normalize 
  
  
  each_med = apply(data_pep, 2, median)
  
  med_each_med = median(each_med)
  
  norm_pep_log2data = sweep(data_pep, 2, each_med) + med_each_med
  
  ### get fold change
  #### how to keep tracking the gene name of each protein 
  

  
  
  
  norm_log_peptide_df = cbind(peptide_df[,c(1,2,3)], norm_pep_log2data)%>%
    dplyr::group_by(protein, description)%>%
    dplyr::summarise_if(is.numeric, mean, na.rm = T)%>%
    dplyr::ungroup()
  
  
  norm_log_peptide_df = as.data.frame(norm_log_peptide_df)
  
  cond_names = grep("cond_", colnames(norm_log_peptide_df), value = T)
  num_rep = length(cond_names)/3
  
  ### get correct cols 
  ctrl_cols = rep(0, num_rep)
  low_cols = rep(0, num_rep)
  high_cols = rep(0, num_rep)
  
  
  for(i in 1:num_rep)
  {
    
    ctrl_cols[i] = which(colnames(norm_log_peptide_df) == paste0("cond_ctrl_", i))
    low_cols[i] = which(colnames(norm_log_peptide_df) == paste0("cond_low_", i))
    high_cols[i] = which(colnames(norm_log_peptide_df) == paste0("cond_high_", i))
    
    
  }
  
  
  lc_fc_rep = rep(0, num_rep)
  hc_fc_rep = rep(0, num_rep)
  

  
  data_df = rbindlist(lapply(1:nrow(norm_log_peptide_df), function(x) {
    
    for(i in 1:num_rep)
    {
      lc_fc_rep[i] = norm_log_peptide_df[x, low_cols[i]] - norm_log_peptide_df[x, ctrl_cols[i]]
      hc_fc_rep[i] = norm_log_peptide_df[x, high_cols[i]] - norm_log_peptide_df[x, ctrl_cols[i]]
      
    }
    
    lc_fc_mean = mean(lc_fc_rep, na.rm = T)
    hc_fc_mean = mean(hc_fc_rep, na.rm = T)
    
    lc_fc_sd = sd(lc_fc_rep, na.rm = T)
    hc_fc_sd = sd(hc_fc_rep, na.rm = T)
    
    ### constrauct the order for l and high vs control
    
    lc_fc_order = paste(order(lc_fc_rep), collapse = "_")
    hc_fc_order = paste(order(hc_fc_rep), collapse = "_")
    
    this_df = data.frame(t(lc_fc_rep), t(hc_fc_rep), 
                         lc_fc_mean, hc_fc_mean,
                         lc_fc_sd, hc_fc_sd,
                         lc_fc_order, hc_fc_order,
                         stringsAsFactors = F)
    
    
    colnames(this_df)[1:(2*num_rep)] = c(paste0("lc_fc_", seq(1:num_rep)), paste0("hc_fc_", seq(1:num_rep)))
      
    if(x%%1000 == 0 )
      cat(x,"\n")
    
      return(this_df)
      
    
  }))
  
  
  norm_log_peptide_data_df = cbind(norm_log_peptide_df, data_df)
  
  write.table(norm_log_peptide_data_df, 
              output_name,
              quote = F, row.names = F, sep = "\t")
  
  
  return(norm_log_peptide_data_df)
  
  
  
  
}






testReliability_afterFC_r = function(peptide_fc,
                                   output_name)
{
    # peptide_fc = fc37
    # output_name = reliability_outputname1
    # 
  #   
    # 
    # it = peptide_fc%>%
    #   dplyr::select(protein, description)%>%
    #   unique()
    lc_fc_cols = grep("lc_fc_[0-9]+",colnames(peptide_fc))
    hc_fc_cols = grep("hc_fc_[0-9]+",colnames(peptide_fc))
    
    test_drl_sign = apply(peptide_fc[,lc_fc_cols], 1, test_equal_sign)
    test_drh_sign = apply(peptide_fc[,hc_fc_cols], 1, test_equal_sign)
    
    
    rel_data = peptide_fc%>%
      dplyr::mutate(drl = test_drl_sign, drh = test_drh_sign)
    
  
  #### return some lables 
  rel_dr = rel_data%>%
    dplyr::group_by(protein, description)%>%
    dplyr::summarise(npep = n(),
                     same_drl = npos(drl), 
                     diff_drl = nzero(drl),
                     same_drh = npos(drh),
                     diff_drh = nzero(drh))%>%
    dplyr::ungroup()%>%
    dplyr::mutate(rel_low = case_when(diff_drl == 0 ~ T,
                                      diff_drl > 0 ~ F),
                  rel_high = case_when(diff_drh == 0 ~ T,
                                       diff_drh >0 ~ F))
  
  
  write.table(rel_dr, output_name, 
              quote = F, row.names = F, sep = "\t")
  
  return(rel_dr)
  
  
}






mapDIApro_r = function(num_rep,
                       working_dir,
                     mapDIA_dir,
                     userSystem,
                     
                     default_input_filename,
                     change_input_filename,
                     peptide_df_filename,
                     analysis_output_name)
  
  
{
  # num_rep = rep_number
  # working_dir = workingDir
  # mapDIA_dir = mapDIADir
  # userSystem = userSystem
  # default_input_filename = defaultInput_parameter
  # change_input_filename =   paste0(workingDir, "input1.txt")
  # peptide_df_filename = cleanData_outputname1
  # analysis_output_name = "mapDIA_output1.txt"
  # 
  
  
  default_input = readLines(default_input_filename)
  
  
  
  #set_SDF = 2
  
  default_input[2] = paste0("FILE= ", peptide_df_filename)
  #default_input[17] = paste0("SDF= ", set_SDF)
  default_input[25] = paste0("SIZE= ", num_rep)
  
  
  
  
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







joinScatterOutput_r = function(pepFC1,
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
                             full_output_name,
                             low_pdf_name,
                             high_pdf_name)

{
  # 
  # pepFC1 = fc37
  # pepFC2 = fc52
  # fc_upper = cutoff_uFc
  # fc_lower = cutoff_lFc
  # mapDIA_output_filename1 = mapDIA_out1
  # mapDIA_output_filename2 = mapDIA_out2
  # mapDIA_sigCutoff = cutoff_mapDIA
  # prot_reliable_filename1 = reliability_outputname1
  # prot_reliable_filename2 = reliability_outputname2
  # low_output_name = low_outputname
  # high_output_name = high_outputname
  # low_pdf_name  = low_pdfname
  # high_pdf_name = high_pdfname
  # 
  # 
  
  
  join_pep = pepFC1%>%
    dplyr::left_join(pepFC2, by = c("protein","description"))%>%
    na.omit()
  
  # 
  join_pep_low = join_pep%>%
    dplyr::select(protein,description, lc_fc_mean.x, lc_fc_mean.y, lc_fc_sd.x, lc_fc_sd.y)

  
  join_pep_high = join_pep%>%
    dplyr::select(protein,description, hc_fc_mean.x, hc_fc_mean.y, hc_fc_sd.x, hc_fc_sd.y)
  
  
  
  
  
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
      dplyr::select(protein,lc_fc_mean.x, lc_fc_mean.y, description, FDR.x, FDR.y, sigFlag,lc_fc_sd.x, lc_fc_sd.y)
    
    
    join_fdr_high = join_fdr_des %>%
      dplyr::filter(Label2 == "high/control")%>%
      na.omit()%>%
      dplyr::mutate(sigFlag = case_when(FDR.x<=mapDIA_sigCutoff & FDR.y <= mapDIA_sigCutoff ~ "bothSig",
                                        FDR.x<=mapDIA_sigCutoff & FDR.y > mapDIA_sigCutoff ~ "firstSig",
                                        FDR.x>mapDIA_sigCutoff & FDR.y <= mapDIA_sigCutoff ~ "secondSig",
                                        FDR.x>mapDIA_sigCutoff & FDR.y > mapDIA_sigCutoff ~ "noneSig"))
    
    join_pep_high = join_pep_high%>%
      dplyr::left_join(join_fdr_high, by = c("protein","description"))%>%
      dplyr::select(protein,hc_fc_mean.x, hc_fc_mean.y, description, FDR.x, FDR.y, sigFlag,hc_fc_sd.x, hc_fc_sd.y)
    
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
  plot(symbol_join_pep_low$lc_fc_mean.x, symbol_join_pep_low$lc_fc_mean.y,
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
    dplyr::filter(lc_fc_mean.x>log2(fc_upper)| lc_fc_mean.x<log2(fc_lower) | lc_fc_mean.y>log2(fc_upper)| lc_fc_mean.y<log2(fc_lower))
  
  
  
  text(jitter(join_low_hit$lc_fc_mean.x),jitter(join_low_hit$lc_fc_mean.y), 
       labels = join_low_hit$symbol, pos = 1, cex = 0.25)
  
  
  
  
  if("relFlag"%in%colnames(join_pep_low))
  {
    
    plot_rel = symbol_join_pep_low%>%
      dplyr::filter(relFlag == "bothRel"| relFlag == "firstRel" | relFlag == "secondRel")
    rel_pch = rep(16, nrow(plot_rel))
    rel_pch[which(plot_rel$relFlag == "bothRel")] = 2
    rel_pch[which(plot_rel$relFlag == "firstRel")] = 0
    rel_pch[which(plot_rel$relFlag == "secondRel")] = 0
    
    
    lines(plot_rel$lc_fc_mean.x, plot_rel$lc_fc_mean.y,
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
    
    
    lines(plot_sig$lc_fc_mean.x, plot_sig$lc_fc_mean.y,
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
  plot(symbol_join_pep_high$hc_fc_mean.x, symbol_join_pep_high$hc_fc_mean.y,
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
    dplyr::filter(hc_fc_mean.x>log2(fc_upper)| hc_fc_mean.x<log2(fc_lower) | hc_fc_mean.y>log2(fc_upper)| hc_fc_mean.y<log2(fc_lower))
  
  
  text(jitter(join_high_hit$hc_fc_mean.x, factor = 5),jitter(join_high_hit$hc_fc_mean.y, factor = 5), labels = join_high_hit$symbol, pos = 1, cex = 0.25)
  
  
  
  
  
  if("relFlag"%in%colnames(join_pep_high))
  {
    
    plot_rel = symbol_join_pep_high%>%
      dplyr::filter(relFlag == "bothRel"| relFlag == "firstRel" | relFlag == "secondRel")
    rel_pch = rep(16, nrow(plot_rel))
    rel_pch[which(plot_rel$relFlag == "bothRel")] = 2
    rel_pch[which(plot_rel$relFlag == "firstRel")] = 0
    rel_pch[which(plot_rel$relFlag == "secondRel")] = 0
    
    
    lines(plot_rel$hc_fc_mean.x, plot_rel$hc_fc_mean.y,
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
    
    
    lines(plot_sig$hc_fc_mean.x, plot_sig$hc_fc_mean.y,
          type = "p", pch = 16, cex = 0.5,
          col = sig_col)
    
    
    
  }
  
  
  
  legend(1, 1.5, legend=c("bothSig", "oneSig","bothRel","oneRel"),
         col=c("red", "green","black","black"), pch = c(16,16,2,0), cex=0.8)
  
  
  
  
  dev.off()
  
  
  write.table(join_pep, full_output_name, quote = F, row.names = F, sep = "\t")
  
  return(join_pep)
  
  
  
}
















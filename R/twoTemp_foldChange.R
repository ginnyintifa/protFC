### hdac wrapper function 

#' This function allows you to calculate the fold change of two different temperatures and to visualize.  
#' @param rawDate_inputfilename1 absolute path to your input file of temperature 1 
#' @import dplyr magrittr data.table 
#' @export

twoTemp_foldChange = function(rawData_inputfilename1,
                           rawData_inputfilename2,
                           seqCol = 3,
                           protAccCol = 10,
                           desCol = 13,
                           sp1ctrlCol = 16,
                           sp2ctrlCol = 17,
                           sp1lowCol = 18,
                           sp2lowCol = 19,
                           sp1highCol = 20,
                           sp2highCol = 21,
                           workingDir,
                           mapDIA_flag = F,
                           mapDIADir,
                           userSystem = "Linux",
                           defaultInput_parameter,
                           cutoff_uFc = 1.3,
                           cutoff_lFc = 0.7,
                           cutoff_mapDIA = 0.05,
                           outputTagS1,
                           outputTagS2)

{
   
  
  
  cleanData_outputname1 = paste0(workingDir, "clean1.tsv")
  cleanData_outputname2 = paste0(workingDir, "clean2.tsv")
  
  pep37 = readInTidy(raw_data_filename = rawData_inputfilename1,
                     seq = seqCol,
                     protAcc = protAccCol,
                     des = desCol,
                     sp1ctrl = sp1ctrlCol,
                     sp2ctrl = sp2ctrlCol,
                     sp1low = sp1lowCol,
                     sp2low = sp2lowCol,
                     sp1high = sp1highCol,
                     sp2high = sp2highCol,
                     output_name = cleanData_outputname1)
  
  
  pep52 = readInTidy(raw_data_filename = rawData_inputfilename2,
                     seq = seqCol,
                     protAcc = protAccCol,
                     des = desCol,
                     sp1ctrl = sp1ctrlCol,
                     sp2ctrl = sp2ctrlCol,
                     sp1low = sp1lowCol,
                     sp2low = sp2lowCol,
                     sp1high = sp1highCol,
                     sp2high = sp2highCol,
                     output_name = cleanData_outputname2)
  
  cat ("original files cleaned." , "\n")
  
  
  foldChange_outputname1 = paste0(workingDir, "foldChange1.tsv")
  foldChange_outputname2 = paste0(workingDir, "foldChange2.tsv")
  
  
  fc37 = calFoldChange1(peptide_df = pep37,
                        output_name = foldChange_outputname1)
  
  
  fc52 = calFoldChange1(peptide_df = pep52,
                        output_name = foldChange_outputname2)
  
  
  cat ("fold change calculated." , "\n")
  
  
  
  reliability_outputname1 = paste0(workingDir, "rel1.tsv")
  reliability_outputname2 = paste0(workingDir, "rel2.tsv")
  
  
  rel37 = testReliability_afterFC(peptide_fc = fc37,
                                  output_name = reliability_outputname1)
  
  
  rel52 = testReliability_afterFC(peptide_fc = fc52,
                                  output_name = reliability_outputname2)
  
  
  cat ("reliability assessed." , "\n")
  
  mapDIA_out1 = ""
  mapDIA_out2 = ""
  
  if(mapDIA_flag == T)
  {
    
    cat("start mapDIA","/n")
    
    
    mapDIApro(working_dir = workingDir,
              mapDIA_dir = mapDIADir,
              userSystem = userSystem,
              default_input_filename = defaultInput_parameter,
              change_input_filename =   paste0(workingDir, "input1.txt"),
              peptide_df_filename = cleanData_outputname1,
              analysis_output_name = "mapDIA_output1.txt")
    
    
    
    mapDIApro(working_dir = workingDir,
              mapDIA_dir = mapDIADir,
              userSystem = userSystem,
              default_input_filename =  defaultInput_parameter,
              change_input_filename =   paste0(workingDir, "input2.txt"),
              peptide_df_filename = cleanData_outputname2,
              analysis_output_name =  "mapDIA_output2.txt")
    
    
    mapDIA_out1 = paste0(workingDir, "mapDIA_output1.txt")
    mapDIA_out2 = paste0(workingDir, "mapDIA_output2.txt")
    
    
    
    
    
    cat ("significance calculated by mapDIA." , "\n")
    
    
    
  }
  
  
  
  low_outputname = paste0(workingDir, outputTagS1,".tsv")
  high_outputname = paste0(workingDir, outputTagS2,".tsv")
  
  low_pdfname = paste0(workingDir, outputTagS1,".pdf")
  high_pdfname = paste0(workingDir, outputTagS2,".pdf")
  
  
  join_prot = joinScatterOutput( pepFC1 = fc37,
                                 pepFC2 = fc52,
                                 fc_upper = cutoff_uFc,
                                 fc_lower = cutoff_lFc,
                                 mapDIA_output_filename1 = mapDIA_out1,
                                 mapDIA_output_filename2 = mapDIA_out2,
                                 mapDIA_sigCutoff = cutoff_mapDIA,
                                 prot_reliable_filename1 = reliability_outputname1,
                                 prot_reliable_filename2 = reliability_outputname2,
                                 low_output_name = low_outputname,
                                 high_output_name = high_outputname,
                                 low_pdf_name  = low_pdfname,
                                 high_pdf_name = high_pdfname)
  
  
  cat ("final outputs generated." , "\n")
  cat ("MISSION COMPLETED" , "\n")
  
  
  
}


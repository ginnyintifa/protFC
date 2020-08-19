
# 
# rawData_inputfilename1 = "/data/ginny/YanTingProtein/20200813_TF_MEHP37_F_PeptideGroups.txt"
# rawData_inputfilename2 = "/data/ginny/YanTingProtein/20200813_TF_MEHP53_F_PeptideGroups.txt"
# inputSpec_filename =  "/data/ginny/YanTingProtein/20200813_TF_spCol_positions.txt"
# 
# rep_number = 3
# workingDir = "/data/ginny/YanTingProtein/foldChange0819/"
# mapDIA_flag = T
# mapDIADir = "/data/ginny/mapDIA-master/"
# userSystem = "Linux"
# defaultInput_parameter = "/data/ginny/YanTingProtein/YTinput"
# cutoff_uFc = 1.3
# cutoff_lFc = 0.7
# cutoff_mapDIA = 0.05
# outputTagS1 = "lowVSctrl"
# outputTagS2 = "highVSctrl"
# 


#' @param rawDate_inputfilename1 absolute path to your input file of temperature 1 
#' @import dplyr magrittr data.table 
#' @export

twoTemp_foldChange_r = function(rawData_inputfilename1,
                              rawData_inputfilename2,
                              inputSpec_filename,
                              rep_number,
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



pep37 = readInTidy_r(raw_data_filename = rawData_inputfilename1,
                     data_col_filename = inputSpec_filename,
                      output_name = cleanData_outputname1)
    
    


pep52 = readInTidy_r(raw_data_filename = rawData_inputfilename2,
                     data_col_filename = inputSpec_filename,
                     output_name = cleanData_outputname2)


cat ("original files cleaned." , "\n")

foldChange_outputname1 = paste0(workingDir, "foldChange1.tsv")
foldChange_outputname2 = paste0(workingDir, "foldChange2.tsv")


fc37 = calFoldChange1_r(peptide_df = pep37,
                      output_name = foldChange_outputname1)



fc52 = calFoldChange1_r(peptide_df = pep52,
                        output_name = foldChange_outputname2)


cat ("fold change calculated." , "\n")

reliability_outputname1 = paste0(workingDir, "rel1.tsv")
reliability_outputname2 = paste0(workingDir, "rel2.tsv")



rel37 = testReliability_afterFC_r(peptide_fc = fc37,
                                output_name = reliability_outputname1)




rel52 = testReliability_afterFC_r(peptide_fc = fc52,
                                  output_name = reliability_outputname2)


cat ("reliability assessed." , "\n")

mapDIA_out1 = ""
mapDIA_out2 = ""


if(mapDIA_flag == T)
{
  
  cat("start mapDIA","/n")
  
  
  mapDIApro_r(num_rep = rep_number,
    working_dir = workingDir,
            mapDIA_dir = mapDIADir,
            userSystem = userSystem,
            default_input_filename = defaultInput_parameter,
            change_input_filename =   paste0(workingDir, "input1.txt"),
            peptide_df_filename = cleanData_outputname1,
            analysis_output_name = "mapDIA_output1.txt")
  
  
  
  mapDIApro_r(num_rep = rep_number, 
    working_dir = workingDir,
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
full_outputname = paste0(workingDir, outputTagS1, "_",outputTagS2,".pdf")

low_pdfname = paste0(workingDir, outputTagS1,".pdf")
high_pdfname = paste0(workingDir, outputTagS2,".pdf")




join_prot = joinScatterOutput_r( pepFC1 = fc37,
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
                               full_output_name = full_outputname,
                               low_pdf_name  = low_pdfname,
                               high_pdf_name = high_pdfname)


cat ("final outputs generated." , "\n")
cat ("MISSION COMPLETED" , "\n")




}




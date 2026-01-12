#inputpath=/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2025-10-14_pPb
inputpath=/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2025-10-14_Pbp

root -x -l -b -q plotRawYields2.C'+("", false, '0',"'${inputpath}'")' #, '$ZT_MODE')'


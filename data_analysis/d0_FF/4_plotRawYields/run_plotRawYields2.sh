addpathname=""
# addpathname="_16tagzBins"


inputpath=/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-09-22_pPb_DGaussDefault${addpathname}

root -x -l -b -q plotRawYields2.C'+("", false, '0',"'${inputpath}'")' #, '$ZT_MODE')'

inputpath=/media/niviths/local/analysis_code/data_analysis/d0_FF/2_fitData/D0_FF_DATA_2026-09-22_Pbp_DGaussDefault${addpathname}

root -x -l -b -q plotRawYields2.C'+("", false, '0',"'${inputpath}'")' #, '$ZT_MODE')'


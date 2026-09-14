#infile2=/media/niviths/local/analysis_code/data_analysis/d0_FF/3_makeResponseMatrix/D0_FF_ResMatr_20260210_53_56/D0response_matrices.root
#infile1=/media/niviths/local/analysis_code/data_analysis/d0_FF/3_makeResponseMatrix/D0_FF_ResMatr_20260211_68/D0response_matrices.root
#nameout=mergedVspp15GpthD0
#root -x -l -b -q compareResponseMatrices.'C("'$infile1'","'$infile2'","'$nameout'")'


infile1=/media/niviths/local/analysis_code/data_analysis/d0_FF/3_makeResponseMatrix/D0_FF_ResMatr_20260210_53_56/D0response_matrices.root
infile2=/media/niviths/local/analysis_code/data_analysis/d0_FF/3_makeResponseMatrix/D0_FF_ResMatr_20260223_70_72/D0response_matrices.root
#infile2=/media/niviths/local/analysis_code/data_analysis/d0_FF/3_makeResponseMatrix/D0_FF_ResMatr_20260218_69/D0response_matrices.root
#nameout=69_2018pthVs70plus_2016pth
nameout=53to56_vs_70to72
root -x -l -b -q compareResponseMatrices.'C("'$infile1'","'$infile2'","'$nameout'")'

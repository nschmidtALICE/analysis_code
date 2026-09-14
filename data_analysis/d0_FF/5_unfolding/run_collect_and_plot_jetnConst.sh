#echo filelist for jetnConst unfolding
echo "/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/70_72_filtered.root" > filelist.txt
#append more files if needed
echo "/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/69-72/69_2018sim_D02Kpi_pthgreater15_filtered.root" >> filelist.txt
echo "/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/54_FF_pPb_EPOS_filtered.root" >> filelist.txt

root -l -b -q 'collect_jetnConst.C("filelist.txt","jetnConst_out.root")'


root -l -b -q 'plot_jetnConst.C("jetnConst_out_2026-02-26/jetnConst_out.root","mcProd_compare")'
#if an argument is passed, it is the systematic variation number (0, 1, or 2)
if [ "$#" -eq 1 ]; then
    SYSVAR=$1
else
    SYSVAR=0
fi
isMCswitch="false"

eventfraction=1.0
echo "Running mass fitter with event fraction ${eventfraction}"

#add another variable for massfitter depending on the sysvar input
if [ "$SYSVAR" -eq 0 ]; then
    echo "Running mass fitter for nominal output"
    OUTNAME="DGaussDefault"
elif [ "$SYSVAR" -eq 1 ]; then
    echo "Running mass fitter for systematic variation 1"
    OUTNAME="DGaussTightPID"
elif [ "$SYSVAR" -eq 2 ]; then
    echo "Running mass fitter for systematic variation 2"
    OUTNAME="DGaussLoosePID"
else
    echo "Invalid systematic variation number. Please provide 0, 1, or 2."
    exit 1
fi
# inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/57_FF_pPb_DATA_filtered_sysvar${SYSVAR}.root"
inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/8_Pbp_data_filtered_sysvar${SYSVAR}.root"
root -x -l -b -q MassFitter.C'+("'$inputFileMassFit'",'$isMCswitch',false,true,"DGauss","'${OUTNAME}'", '$eventfraction')'
inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/9_pPb_data_filtered_sysvar${SYSVAR}.root"
# inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/59_Pbp_data_filtered_sysvar${SYSVAR}.root"
root -x -l -b -q MassFitter.C'+("'$inputFileMassFit'",'$isMCswitch',false,true,"DGauss","'${OUTNAME}'", '$eventfraction')'


inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/0_pp_data_full_filtered_sysvar0.root"
# inputFileMassFit="/media/niviths/SSD2/lhcb_analysis_SSD/GANGA/0_pp_Data_filtered_sysvar0.root"
# root -x -l -b -q MassFitter.C'+("'$inputFileMassFit'",'$isMCswitch',false,true,"DGauss","'${OUTNAME}'", '$eventfraction', true)'


pidcalib2.make_eff_hists --sample pATurbo16 --magnet down --particle Pi --pid-cut "MC15TuneV1_ProbNNpi > 0.9 & MC15TuneV1_ProbNNghost < 0.3 & TRCHI2NDOF < 3" --bin-var P --bin-var ETA --binning-file binning.json --output-dir pidcalib_output_pA_09_pEta_pi

pidcalib2.make_eff_hists --sample pATurbo16 --magnet down --particle K --pid-cut "MC15TuneV1_ProbNNk > 0.9 & MC15TuneV1_ProbNNghost < 0.3 & TRCHI2NDOF < 3" --bin-var P --bin-var ETA --binning-file binning.json --output-dir pidcalib_output_pA_09_pEta_k

pidcalib2.make_eff_hists --sample ApTurbo16 --magnet down --particle Pi --pid-cut "MC15TuneV1_ProbNNpi > 0.9 & MC15TuneV1_ProbNNghost < 0.3 & TRCHI2NDOF < 3" --bin-var P --bin-var ETA --binning-file binning.json --output-dir pidcalib_output_Ap_09_pEta_pi

pidcalib2.make_eff_hists --sample ApTurbo16 --magnet down --particle K --pid-cut "MC15TuneV1_ProbNNk > 0.9 & MC15TuneV1_ProbNNghost < 0.3 & TRCHI2NDOF < 3" --bin-var P --bin-var ETA --binning-file binning.json --output-dir pidcalib_output_Ap_09_pEta_k

pidcalib2.plot_calib_distributions --sample pATurbo16 --magnet down --particle Pi --bin-var P --output-dir pidcalib_output --max-files 1 --format pdf --force-uniform --bins 95


pidcalib2.pklhisto2root pidcalib_output_09extfine/effhists-pATurbo16-down-Pi-MC15TuneV1_ProbNNpi\>0.9\&MC15TuneV1_ProbNNghost\<0.3\&TRCHI2NDOF\<3-PT.ETA.pkl 


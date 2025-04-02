#F1  =  ifx -O3  -check all,nouninit -traceback 
#F2  =  ifx -O3 -qopenmp  -check all,nouninit -traceback
F1  =  ifx -O3 
F2  =  ifx -O3 -qopenmp 
LAPACK  = -qmkl -lpthread 
#-- Lapack interface
Lapack    = sub.lapack_inv.f90  $(LAPACK)
#-- Other routines
FILTER    = KIKUCHI/sub.l_filter.f90
TMPSM     = sub.TmpSmWeight.f90 sub.g_Gmatrix_TmpSm.f90
#--------------------------------------------------
all: inv_pre ta ta_d 

inv_pre: ../bin/PDTI_get_knot_info ../bin/PDTI_GreenPointSources ../bin/PDTI_get_ndj_main ../bin/PDTI_pre_inv_new ../bin/PDTI_get_covariance_grn

ta:      ../bin/PDTI_ABIC_TA ../bin/PDTI_get_init_model_para_TmpSm 

ta_d:    ../bin/PDTI_inversionA_TAd 


#--------------------------------------------------
../bin/PDTI_get_knot_info: get_knot_info.f90  sub.c_location_lib.f90 
	$(F1) -o $@ $^

../bin/PDTI_GreenPointSources :  GreenPointSources.f90 KIKUCHI/sub.green_lib1.f90 KIKUCHI/sub.green_lib2.f90 KIKUCHI/sub.l_filter.f90 sub.c_wave_lib.f90 sub.c_cha.f90    sub.taup_lib.f90
	$(F2) -o $@ $^

../bin/PDTI_get_ndj_main : get_ndj_main.f90 sub.c_wave_lib.f90 sub.l_invlibF.f90 sub.r_OBS.f90 sub.g_ndj.f90 sub.g_TR.f90 sub.c_cha.f90 $(FILTER) 
	$(F2) -o $@ $^ $(Lapack)

../bin/PDTI_pre_inv_new : pre_inv_new.f90 sub.c_matrix.f90 sub.c_wave_lib.f90 sub.l_invlibF.f90 sub.r_OBS.f90 sub.r_green_H.f90 sub.g_TR.f90 sub.c_cha.f90 $(FILTER)
	$(F1) -o $@ $^ $(Lapack)

../bin/PDTI_get_covariance_grn :  get_covariance_grn.f90 sub.c_matrix.f90 sub.c_wave_lib.f90 sub.g_TR.f90 sub.get_slip_rate.f90
	$(F2) -o $@ $^ $(Lapack)
#--------------------------------------------------

../bin/PDTI_ABIC_TA : sub.c_matrix.f90 $(TMPSM) sub.l_invlibF.f90 sub.g_TR.f90 sub.l_ABIC_lib.f90 sub.g_covOBS.f90 sub.g_cm.f90 ABIC_TA.f90  sub.r_green_H.f90  sub.g_TR.f90 sub.c_cha.f90  sub.c_wave_lib.f90 sub.l_invlibF.f90 sub.g_dumpS.f90 sub.w_sol.f90  sub.l_disp_lib.f90  KIKUCHI/sub.focalM_lib.f90
	$(F2) -o $@ $^ $(Lapack)

../bin/PDTI_inversionA_TAd : sub.c_matrix.f90 sub.c_wave_lib.f90 $(TMPSM) sub.l_invlibF.f90 KIKUCHI/sub.focalM_lib.f90 sub.g_covOBS.f90 sub.g_cm.f90 sub.r_OBS.f90 $(FILTER) sub.g_TR.f90 sub.c_cha.f90 sub.r_stinfo.f90 sub.g_dumpS.f90 sub.l_disp_lib.f90  sub.r_green_H.f90  sub.w_sol.f90 inversionA_TAd.f90  
	$(F2) -o $@ $^ $(Lapack)

../bin/PDTI_get_init_model_para_TmpSm : get_init_model_para_TmpSm.f90 sub.slip_rate2fort40.f90 sub.l_disp_lib.f90 sub.c_wave_lib.f90 sub.g_TR.f90  KIKUCHI/sub.focalM_lib.f90 
	$(F2) -o $@ $^ $(Lapack)
#--------------------------------------------------

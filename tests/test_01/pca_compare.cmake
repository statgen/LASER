execute_process(COMMAND ${LASER} -g ${GENO} -k 20 -pca 1 -o test_pca RESULT_VARIABLE laser_exit_code)
if(laser_exit_code)
   message(FATAL_ERROR "LASER failed.")
endif()

#execute_process(COMMAND diff -q test_pca.RefPC.coord good/test_pca.RefPC.coord RESULT_VARIABLE diff_exit_code)
#if(diff_exit_code)
#   message(FATAL_ERROR "LASER didn't replicate results.")
#endif()

#execute_process(COMMAND diff -q test_pca.RefPC.var good/test_pca.RefPC.var RESULT_VARIABLE diff_exit_code)
#if(diff_exit_code)
#   message(FATAL_ERROR "LASER didn't replicate results.")
#endif()

#execute_process(COMMAND diff -q test_pca.RefPC.zscore good/test_pca.RefPC.zscore RESULT_VARIABLE diff_exit_code)
#if(diff_exit_code)
#   message(FATAL_ERROR "LASER didn't replicate results.")
#endif()

#execute_process(COMMAND diff -q test_pca.RefPC.grm good/test_pca.RefPC.grm RESULT_VARIABLE diff_exit_code)
#if(diff_exit_code)
#   message(FATAL_ERROR "LASER didn't replicate results.")
#endif()

execute_process(COMMAND ${LASER} -g ${GENO} -k 20 -pca 1 -o test_pca RESULT_VARIABLE laser_exit_code)
if(laser_exit_code)
   message(FATAL_ERROR "LASER failed.")
endif()

execute_process(COMMAND ${TESTLASER} compare_tables good/test_pca.RefPC.coord test_pca.RefPC.coord ssffffffffffffffffffff 0.01 RESULT_VARIABLE test_exit_code)
if(test_exit_code)
   message(FATAL_ERROR "LASER didn't replicate results.")
endif()


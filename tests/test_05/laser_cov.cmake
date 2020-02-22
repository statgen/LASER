execute_process(COMMAND ${LASER} -g ${GENO_REF} -s ${SEQ_STUDY} -cov 2 -o test_laser RESULT_VARIABLE laser_exit_code)
if(laser_exit_code)
   message(FATAL_ERROR "LASER failed.")
endif()

execute_process(COMMAND ${TESTLASER} compare_tables good/test_laser.ind.cov test_laser.ind.cov ssdf 0.001 RESULT_VARIABLE test_exit_code)
if(test_exit_code)
   message(FATAL_ERROR "LASER didn't replicate results.")
endif()

execute_process(COMMAND ${TESTLASER} compare_tables good/test_laser.loc.cov test_laser.loc.cov sdf 0.001 RESULT_VARIABLE test_exit_code)
if(test_exit_code)
   message(FATAL_ERROR "LASER didn't replicate results.")
endif()
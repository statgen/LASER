execute_process(COMMAND ${LASER} -g ${GENO_REF} -c ${COORD_REF} -s ${SEQ_STUDY} -o test_laser RESULT_VARIABLE laser_exit_code)
if(laser_exit_code)
   message(FATAL_ERROR "LASER failed.")
endif()

execute_process(COMMAND ${TESTLASER} compare_tables good/test_laser.SeqPC.coord test_laser.SeqPC.coord ssdfdffff 0.001 RESULT_VARIABLE test_exit_code)
if(test_exit_code)
   message(FATAL_ERROR "LASER didn't replicate results.")
endif()
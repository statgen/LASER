execute_process(COMMAND ${LASER} -g ${GENO_REF} -c ${COORD_REF} -s ${SEQ_STUDY} -o test_laser RESULT_VARIABLE laser_exit_code)
if(laser_exit_code)
   message(FATAL_ERROR "LASER failed.")
endif()



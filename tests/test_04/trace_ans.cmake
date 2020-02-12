execute_process(COMMAND ${TRACE} -g ${GENO_REF} -s ${GENO_STUDY} -k 20 -x 1 -y 700 -o test_trace RESULT_VARIABLE trace_exit_code)
if(trace_exit_code)
   message(FATAL_ERROR "TRACE failed.")
endif()

execute_process(COMMAND ${TESTLASER} compare_tables good/test_trace.ProPC.coord test_trace.ProPC.coord ssddffffffffffffffffffffff 0.1 RESULT_VARIABLE test_exit_code)
if(test_exit_code)
   message(FATAL_ERROR "TRACE didn't replicate results: *.ProPC.coord files differ")
endif()

execute_process(COMMAND ${TESTLASER} compare_tables good/test_trace.RefPC.coord test_trace.RefPC.coord ssffffffffffffffffffff 0.01 RESULT_VARIABLE test_exit_code)
if(test_exit_code)
   message(FATAL_ERROR "TRACE didn't replicate results: *.RefPC.coord files differ")
endif()

execute_process(COMMAND ${TESTLASER} compare_tables good/test_trace.RefPC.var test_trace.RefPC.var df 0.01 RESULT_VARIABLE test_exit_code)
if(test_exit_code)
   message(FATAL_ERROR "TRACE didn't replicate results: *.RefPC.var files differ")
endif()

execute_process(COMMAND ${TRACE} -g ${GENO_REF} -c ${COORD_REF} -s ${GENO_STUDY} -x 1 -y 700 -o test_trace RESULT_VARIABLE trace_exit_code)
if(trace_exit_code)
   message(FATAL_ERROR "TRACE failed.")
endif()

execute_process(COMMAND ${TESTLASER} test_propc_files good/test_trace.ProPC.coord test_trace.ProPC.coord ssddffff RESULT_VARIABLE test_exit_code)
if(test_exit_code)
   message(FATAL_ERROR "TRACE didn't replicate results.")
endif()

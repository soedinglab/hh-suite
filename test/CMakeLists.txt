enable_testing()


add_test(NAME hhblits COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hhblits -h)
add_test(NAME hhsearch COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hhsearch -h)
add_test(NAME hhalign COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hhalign -h)


INSTALL(FILES
        context_data.crf
        context_data.lib
        cs219.lib
        do_not_delete
        do_not_delete.phr
        do_not_delete.pin
        do_not_delete.psq
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
        DESTINATION data
        )

include(ResourceCompiler)

set(COMPILED_RESOURCES
        context_data.crf
        cs219.lib
        )

set(GENERATED_OUTPUT_HEADERS "")
FOREACH(INPUT_FILE ${COMPILED_RESOURCES})
    compile_resource(${INPUT_FILE} OUTPUT_FILE)
    list(APPEND GENERATED_OUTPUT_HEADERS "${OUTPUT_FILE}")
ENDFOREACH()

add_custom_target(generated ALL DEPENDS ${GENERATED_OUTPUT_HEADERS})


# Add the flow library
set(SOURCE_FILES
        adem/adem.h
        adem/biot_savart.h
        adem/signature.h
        data_types.h
        definitions.h
        fit.h
        relations/velocity.h
        relations/stress.h
        relations/spectra.cpp
        relations/spectra.h
        profile.cpp
        profile.h
        io/readers.h
        io/variable_readers.h
        io/variable_writers.h
        utilities/conv.h
        utilities/cumsum.h
        utilities/cumtrapz.h
        utilities/filter.h
        utilities/integration.h
        utilities/interp.h
        utilities/tensors.h
        )
add_library(flow STATIC ${SOURCE_FILES})


# The hello world function which gives an executable target for the library
add_executable(es_flow main.cpp)
target_link_libraries(es_flow ${CERES_LIBRARIES} ${MKL_LIBRARIES} ${TBB_LIBRARIES} ${MATIO_LIBRARIES})


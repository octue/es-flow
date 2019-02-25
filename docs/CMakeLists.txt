
message(INFO "-- Skipping Sphinx documentation build via cmake (using travis instead).")

# Build the documentation, if we can find sphinx
# set(DOCS_INSTALL_DIR "share/doc/es-flow" CACHE STRING "es-flow docs install path relative to CMAKE_INSTALL_PREFIX")
# find_package(Sphinx QUIET)
# if (NOT SPHINX_FOUND)
#     message(WARNING "-- Failed to find Sphinx, disabling build of documentation.")
# else()
#     # Generate the User's Guide (html).
#     # The corresponding target is es_flow_docs, but is included in ALL.
#     message("-- Found Sphinx")
#     add_subdirectory(source)
# endif()

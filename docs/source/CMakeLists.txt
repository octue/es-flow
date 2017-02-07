
find_package(Sphinx REQUIRED)

# HTML output directory
set(SPHINX_HTML_DIR "${CMAKE_BINARY_DIR}/docs/html")

# Install documentation
install(DIRECTORY ${SPHINX_HTML_DIR}
        DESTINATION "${DOCS_INSTALL_DIR}"
        COMPONENT Doc
        PATTERN "${SPHINX_HTML_DIR}/*")

# Building using 'make_docs.py' python script
#add_custom_target(es_flowdoc ALL ${SPHINX_EXECUTABLE} -b html "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/_build" COMMENT Building documentation)

add_custom_target(es_flow_docs ALL
        python
        "${CMAKE_SOURCE_DIR}/scripts/make_docs.py"
        "${CMAKE_SOURCE_DIR}"
        "${CMAKE_BINARY_DIR}/docs"
        "${SPHINX_EXECUTABLE}"
        COMMENT "Building HTML documentation with Sphinx")

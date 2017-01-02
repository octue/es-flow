find_package(Sphinx REQUIRED)

## HTML output directory
#MESSAGE("${CMAKE_BINARY_DIR}")
#set(SPHINX_HTML_DIR "${CMAKE_BINARY_DIR}/docs/html")
#
## Install documentation
#install(DIRECTORY ${SPHINX_HTML_DIR}
#        DESTINATION "${DOCS_INSTALL_DIR}"
#        COMPONENT Doc
#        PATTERN "${SPHINX_HTML_DIR}/*")

## Building using 'make_docs.py' python script
#add_custom_target(es_flow_docs ALL
#                  python
#                  "${CMAKE_SOURCE_DIR}/scripts/make_docs.py"
#                  "${CMAKE_SOURCE_DIR}"
#                  "${CMAKE_BINARY_DIR}/docs"
#                  "${SPHINX_EXECUTABLE}"
#                  COMMENT "Building HTML documentation with Sphinx")

if(NOT DEFINED SPHINX_THEME)
    set(SPHINX_THEME default)
endif()

if(NOT DEFINED SPHINX_THEME_DIR)
    set(SPHINX_THEME_DIR)
endif()

# configured documentation tools and intermediate build results
set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

# HTML output directory
set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

configure_file(
        "../conf.py.in"
        "${BINARY_BUILD_DIR}/conf.py"
        @ONLY)

add_custom_target(es_flow_docs ALL
        ${SPHINX_EXECUTABLE}
        -q -b html
        -c "${BINARY_BUILD_DIR}"
        -d "${SPHINX_CACHE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${SPHINX_HTML_DIR}"
        COMMENT "Building HTML documentation with Sphinx")

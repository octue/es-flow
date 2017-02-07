# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2015 Google Inc. All rights reserved.
# http://ceres-solver.org/
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Google Inc. nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

add_subdirectory(source)
#find_package(Sphinx REQUIRED)
#
### HTML output directory
##MESSAGE("${CMAKE_BINARY_DIR}")
##set(SPHINX_HTML_DIR "${CMAKE_BINARY_DIR}/docs/html")
##
### Install documentation
##install(DIRECTORY ${SPHINX_HTML_DIR}
##        DESTINATION "${DOCS_INSTALL_DIR}"
##        COMPONENT Doc
##        PATTERN "${SPHINX_HTML_DIR}/*")
#
### Building using 'make_docs.py' python script
##add_custom_target(es_flow_docs ALL
##                  python
##                  "${CMAKE_SOURCE_DIR}/scripts/make_docs.py"
##                  "${CMAKE_SOURCE_DIR}"
##                  "${CMAKE_BINARY_DIR}/docs"
##                  "${SPHINX_EXECUTABLE}"
##                  COMMENT "Building HTML documentation with Sphinx")
#
##if(NOT DEFINED SPHINX_THEME)
##    set(SPHINX_THEME default)
##endif()
##
##if(NOT DEFINED SPHINX_THEME_DIR)
##    set(SPHINX_THEME_DIR)
##endif()
#
### configured documentation tools and intermediate build results
##set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")
##
### Sphinx cache with pickled ReST documents
##set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")
#
## HTML output directory
#set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")
##MESSAGE("CMKAKE_CURRENT_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}")
##
##configure_file(
##        "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
##        "${CMAKE_CURRENT_BINARY_DIR}/conf.py"
##        @ONLY)
#
##file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/_static/" DESTINATION "${BINARY_BUILD_DIR}/_static")

#add_custom_target(es_flowdoc ALL ${SPHINX_EXECUTABLE} -b html "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/_build" COMMENT Building documentation)

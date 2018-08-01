project ( LARC )

cmake_minimum_required ( VERSION 2.8 )

include ( CTest )

########################
# Select the LARC Type #
########################

if ( NOT ( (REAL STREQUAL LARC_TYPE) OR (INTEGER STREQUAL LARC_TYPE) OR (COMPLEX STREQUAL LARC_TYPE) ) )

  message ( " " )
  message ( " LARC_TYPE not specified, or not recognized." )
  message ( " To specify a type add one of -DLARC_TYPE=REAL -DLARC_TYPE=COMPLEX -DLARC_TYPE=INTEGER." )
  message ( " " )

  set ( LARC_TYPE REAL )

endif ( NOT ( (REAL STREQUAL LARC_TYPE) OR (INTEGER STREQUAL LARC_TYPE) OR (COMPLEX STREQUAL LARC_TYPE) ) )

message ( " " )
message ( " Using ${LARC_TYPE} for LARC_TYPE" )
message ( " " )

configure_file ( ${PROJECT_SOURCE_DIR}/src/.type${LARC_TYPE}.h ${PROJECT_SOURCE_DIR}/src/type.h )

################################
# Enable testing with Valgrind #
################################

if ( True STREQUAL TEST_WITH_VALGRIND ) 

  message ( " " ) 
  message ( " Running unit tests under valgrind" )
  message ( " " ) 

  set ( VG   "valgrind"           )
  set ( VGO1 "--error-exitcode=1" )
  set ( VGO2 "--leak-check=full"  )

endif ( True STREQUAL TEST_WITH_VALGRIND ) 

##################
# Enable Doxygen #
##################

find_package ( Doxygen )

if ( DOXYGEN_FOUND )

  configure_file ( ${PROJECT_SOURCE_DIR}/src/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/src/Doxyfile @ONLY )

  add_custom_target (
    doc ALL
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/src/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )

endif ( DOXYGEN_FOUND )

install ( DIRECTORY "${CMAKE_BINARY_DIR}/src/html" DESTINATION "${PROJECT_SOURCE_DIR}/install/doc" )
  
##############################
# Compiler specific messages #
##############################

if ( NOT ( GNU STREQUAL ${CMAKE_C_COMPILER_ID} ) )

  message ( " " )
  message ( " The compiler has been set to \"${CMAKE_C_COMPILER_ID}\"." )
  message ( " You may need to set the LD_LIBRARY_PATH environment variable for tests to run." )   
  message ( " " )

endif ( NOT ( GNU STREQUAL ${CMAKE_C_COMPILER_ID} ) )

#####################################
# Add the components of the project #
#####################################

# Build the larc library, and python interface
add_subdirectory ( src ) 

# Add the CPP wrappers.
add_subdirectory ( CPPWrapper )

# Add our main test directory
add_subdirectory ( tests )


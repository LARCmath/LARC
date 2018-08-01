#############################
# Options for building larc #
#############################

# NB The Intel compiler supports the gnu99 extensions.
set ( CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}" )

# Report all warnings.
set ( CMAKE_C_FLAGS "-Wall ${CMAKE_C_FLAGS}" )

#######################
# Building liblarc.so #
#######################

# Build the larc shared library.
add_library ( # Library name
              larc
	      
	      # Library type
	      SHARED 

	      # Source files
	      fft.c
              global.c  
              hash.c  
              info_store.c  
              io.c  
              json.c  
              logging.c  
              matmath.c  
              matrix_store.c  
              op_store.c  
	      organize.c
	      matrix_accessor.c            
	    )
             
# Specify the other libraries on which the larc library depends.  
target_link_libraries ( larc m rt gmp pthread )

########################################
# Building python interfaces with SWIG #
########################################

# NB Searching for the python intepretor package first ensures that we
# will find a compatible version of the python libs.

find_package ( PythonInterp )
find_package ( PythonLibs   )

find_package ( SWIG REQUIRED    )
include      ( ${SWIG_USE_FILE} )

include_directories         ( ${PYTHON_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} ) 
set_source_files_properties ( larc_py.i PROPERTIES C ON                          )
swig_add_module             ( larc_py python larc_py.i                           )
swig_link_libraries         ( larc_py ${PYTHON_LIBRARIES} larc gmp m rt pthread  )

########################
# Build an application #
########################

add_executable ( exampleLARC exampleLARC.c )
target_link_libraries ( exampleLARC larc )

##############################
# Install our build products #
##############################

# Make sure that the larc library gets installed.
install ( TARGETS 

          # The target, in this case liblarc.so, to install	  
          larc 

          DESTINATION ${PROJECT_SOURCE_DIR}/install/lib 
        )

# Also install the header files.
install ( FILES 

          # Header files	  
          global.h 
          hash.h 
          info_store.h 
          io.h 
          json.h 
          matmath.h 
          matrix_store.h 
          organize.h 
          op_store.h 
          version.h 

          DESTINATION ${PROJECT_SOURCE_DIR}/install/include 
        )

# Install the SWIG generated Python interface
install ( FILES 

          ${CMAKE_CURRENT_BINARY_DIR}/_larc_py.so 
          ${CMAKE_CURRENT_BINARY_DIR}/larc_py.py 

          DESTINATION ${PROJECT_SOURCE_DIR}/install/lib 
        )

# Build a sample LARC application
install ( TARGETS

	  exampleLARC

	  DESTINATION ${PROJECT_SOURCE_DIR}/install/bin
        )

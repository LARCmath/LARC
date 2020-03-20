##                   CMakeLists.txt
###################################################################
 #                                                                #
 # Copyright (C) 2014, Institute for Defense Analyses             #
 # 4850 Mark Center Drive, Alexandria, VA; 703-845-2500           #
 # This material may be reproduced by or for the US Government    #
 # pursuant to the copyright license under the clauses at DFARS   #
 # 252.227-7013 and 252.227-7014.                                 #
 #                                                                #
 # LARC : Linear Algebra via Recursive Compression                #
 # Authors:                                                       #
 #   - Steve Cuccaro (IDA-CCS)                                    #
 #   - John Daly (LPS)                                            #
 #   - John Gilbert (UCSB, IDA adjunct)                           #
 #   - Jenny Zito (IDA-CCS)                                       #
 #                                                                #
 # Additional contributors are listed in "LARCcontributors".      #
 #                                                                #
 # Questions: larc@super.org                                      #
 #                                                                #
 # All rights reserved.                                           #
 #                                                                #
 # Redistribution and use in source and binary forms, with or     #
 # without modification, are permitted provided that the          #
 # following conditions are met:                                  #
 #   - Redistribution of source code must retain the above        #
 #     copyright notice, this list of conditions and the          #
 #     following disclaimer.                                      #
 #   - Redistribution in binary form must reproduce the above     #
 #     copyright notice, this list of conditions and the          #
 #     following disclaimer in the documentation and/or other     #
 #     materials provided with the distribution.                  #
 #   - Neither the name of the copyright holder nor the names of  #
 #     its contributors may be used to endorse or promote         #
 #     products derived from this software without specific prior #
 #     written permission.                                        #
 #                                                                #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND         #
 # CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,    #
 # INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF       #
 # MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE       #
 # DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER NOR        #
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   #
 # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT   #
 # NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;   #
 # LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)       #
 # HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN      #
 # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   #
 # OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, #
 # EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             #
 #                                                                #
 ##################################################################


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
set_source_files_properties ( pylarc.i PROPERTIES C ON                          )
swig_add_module             ( pylarc python pylarc.i                           )
swig_link_libraries         ( pylarc ${PYTHON_LIBRARIES} larc gmp m rt pthread  )

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

          ${CMAKE_CURRENT_BINARY_DIR}/_pylarc.so 
          ${CMAKE_CURRENT_BINARY_DIR}/pylarc.py 

          DESTINATION ${PROJECT_SOURCE_DIR}/install/lib 
        )

# Build a sample LARC application
install ( TARGETS

	  exampleLARC

	  DESTINATION ${PROJECT_SOURCE_DIR}/install/bin
        )

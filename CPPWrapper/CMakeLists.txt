# NB as soon as we switch to a compiler supporting C++11, this should
# be changed to -std=c++11
set ( CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}" )
set ( CMAKE_CXX_FLAGS "-Werror ${CMAKE_CXX_FLAGS}" )

if ( GNU STREQUAL ${CMAKE_C_COMPILER_ID} )

   set ( CMAKE_CXX_FLAGS "-DUSE_GNU ${CMAKE_CXX_FLAGS}" )

endif ( GNU STREQUAL ${CMAKE_C_COMPILER_ID} )

include_directories ( ${CMAKE_SOURCE_DIR}/src )

add_library           ( larc++ SHARED larc.cpp indexable.cpp )
target_link_libraries ( larc++ larc )

install ( TARGETS larc++  DESTINATION ${PROJECT_SOURCE_DIR}/install/lib     )
install ( FILES   larc.hh DESTINATION ${PROJECT_SOURCE_DIR}/install/include )

set ( CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}" )
set ( CMAKE_C_FLAGS "-Wall ${CMAKE_C_FLAGS}"      )

# Specify where to find the necessary include directories.
include_directories ( ${PROJECT_SOURCE_DIR}/src )

add_executable        ( C_test1 test1.c         )
target_link_libraries ( C_test1 larc            )
add_test              ( create_a_matrix C_test1 )

add_executable        ( C_test2 test2.c          )
target_link_libraries ( C_test2 larc             )
add_test              ( matrix_accessors C_test2 )

if ( TEST_WITH_VALGRIND ) 

   add_test ( create_a_matrix_valgind  ${VG} ${VG01} ${VG02} ${CMAKE_CURRENT_BINARY_DIR}/C_test1 )
   add_test ( matrix_accessors_valgind ${VG} ${VG01} ${VG02} ${CMAKE_CURRENT_BINARY_DIR}/C_test2 )

endif ( TEST_WITH_VALGRIND ) 

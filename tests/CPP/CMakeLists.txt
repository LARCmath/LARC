# NB as soon as we switch to a compiler supporting C++11, this should
# be changed to -std=c++11
set ( CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}" )
set ( CMAKE_CXX_FLAGS "-Werror ${CMAKE_CXX_FLAGS}" )

include_directories ( ${CMAKE_SOURCE_DIR}/src        )
include_directories ( ${CMAKE_SOURCE_DIR}/CPPWrapper )

add_executable        ( CPP_test1  test1.cpp )
target_link_libraries ( CPP_test1  larc++    )
add_test              ( CPP_Access CPP_test1 )

add_executable        ( CPP_test2 test2.cpp )
target_link_libraries ( CPP_test2 larc++    )
add_test              ( CPP_Add   CPP_test2 )

add_executable        ( CPP_test3 test3.cpp )
target_link_libraries ( CPP_test3 larc++    )
add_test              ( CPP_Mult  CPP_test3 )

add_executable        ( CPP_test4     test4.cpp )
target_link_libraries ( CPP_test4     larc++    )
add_test              ( CPP_Operators CPP_test4 )

add_executable        ( CPP_test5            test5.cpp )
target_link_libraries ( CPP_test5            larc++    )
add_test              ( CPP_KroneckerProduct CPP_test5 )

add_executable        ( CPP_test6   test6.cpp )
target_link_libraries ( CPP_test6   larc++    )
add_test              ( CPP_Adjoint CPP_test6 )

add_executable        ( CPP_test7 test7.cpp )
target_link_libraries ( CPP_test7 larc++    )
add_test              ( CPP_Trace CPP_test7 )

add_executable        ( CPP_test8     test8.cpp )
target_link_libraries ( CPP_test8     larc++    )
add_test              ( CPP_JoinStack CPP_test8 )

if ( TEST_WITH_VALGRIND ) 

   add_test ( CPP_Access_valgrind          ${VG} ${VG01} ${VG02} ${CMAKE_CURRENT_BINARY_DIR}/CPP_test1 )
   add_test ( CPP_Add_valgind              ${VG} ${VG01} ${VG02} ${CMAKE_CURRENT_BINARY_DIR}/CPP_test2 )
   add_test ( CPP_Mult_valgind             ${VG} ${VG01} ${VG02} ${CMAKE_CURRENT_BINARY_DIR}/CPP_test3 )
   add_test ( CPP_Operators_valgind        ${VG} ${VG01} ${VG02} ${CMAKE_CURRENT_BINARY_DIR}/CPP_test4 )
   add_test ( CPP_KroneckerProduct_valgind ${VG} ${VG01} ${VG02} ${CMAKE_CURRENT_BINARY_DIR}/CPP_test5 )
   add_test ( CPP_Adjoint_valgind          ${VG} ${VG01} ${VG02} ${CMAKE_CURRENT_BINARY_DIR}/CPP_test6 )
   add_test ( CPP_Trace_valgind            ${VG} ${VG01} ${VG02} ${CMAKE_CURRENT_BINARY_DIR}/CPP_test7 )
   add_test ( CPP_JoinStack_valgind        ${VG} ${VG01} ${VG02} ${CMAKE_CURRENT_BINARY_DIR}/CPP_test8 )

endif ( TEST_WITH_VALGRIND ) 

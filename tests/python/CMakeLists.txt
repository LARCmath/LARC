add_test ( # Name of test
	   sample_python_unittest 

	   # The command to run
	   python "-B" "${CMAKE_CURRENT_SOURCE_DIR}/test_unittest_example.py"

	   # Specify the path for the modules and libraries
	   "${CMAKE_BINARY_DIR}/src"
         )


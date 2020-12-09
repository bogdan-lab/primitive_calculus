enable_testing()

add_executable(test_newton_1D
		test_newton_1D.c)

target_link_libraries(test_newton_1D PRIVATE newton_method)

add_test(
	NAME TestNewton1D
	COMMAND test_newton_1D)

add_executable(test_other_methods
                test_other_methods.c)

target_link_libraries(test_other_methods PRIVATE newton_method)

add_test(NAME TestOtherMethods
        COMMAND test_other_methods)

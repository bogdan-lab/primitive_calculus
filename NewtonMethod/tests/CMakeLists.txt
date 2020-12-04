enable_testing()

add_executable(test_newton_1D
		test_newton_1D.c)

target_link_libraries(test_newton_1D PRIVATE newton_method)

add_test(
	NAME TestNewton1D
	COMMAND test_newton_1D)

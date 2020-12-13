enable_testing()

add_executable(test_newton_1D
		test_newton_1D.c)

target_link_libraries(test_newton_1D PRIVATE newton_method)

add_test(NAME TestNewton1D
	COMMAND test_newton_1D)

add_executable(bisection_test
                test_bisection.c)

target_link_libraries(bisection_test PRIVATE newton_method)

add_test(NAME BisectionTest
        COMMAND bisection_test)

add_executable(secant_test
                test_secant.c)
target_link_libraries(secant_test PRIVATE newton_method)

add_test(NAME SecantTest
        COMMAND secant_test)

add_executable(third_order_test
        test_third_order.c)
target_link_libraries(third_order_test PRIVATE newton_method)

add_test(NAME ThirdOrderTest
        COMMAND third_order_test)

add_executable(compare_methods
        compare_methods.c)
target_link_libraries(compare_methods PRIVATE newton_method)

add_test(NAME CompareMethods
        COMMAND compare_methods)



add_library(nonlinear_equations_nd
	SHARED
	newton_method.c)


target_include_directories(nonlinear_equations_nd PUBLIC ../include /usr/include/)
target_link_directories(nonlinear_equations_nd PUBLIC /usr/lib/)
target_compile_options(nonlinear_equations_nd PUBLIC -Wall -Wextra -Wshadow -Werror)
set_target_properties(nonlinear_equations_nd PROPERTIES C_STANDARD 99)
target_link_libraries(nonlinear_equations_nd PUBLIC m lapacke gslcblas)

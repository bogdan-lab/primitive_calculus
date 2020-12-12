add_library(newton_method SHARED
                newton_1D.c
                bisection.c
                secant.c
                third_order_newton.c)

target_link_libraries(newton_method PUBLIC m)
target_include_directories(newton_method PUBLIC ../include)
target_compile_options(newton_method PUBLIC -Wall -Wextra -Wshadow -Werror)
set_target_properties(newton_method PROPERTIES C_STANDARD 99)

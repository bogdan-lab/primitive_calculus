add_library(newton_method SHARED
                newton_1D.c
                other_methods.c)

target_link_libraries(newton_method PUBLIC m)
target_include_directories(newton_method PUBLIC ../include)
target_compile_options(newton_method PUBLIC -Wall -Wextra -Wshadow -Werror)

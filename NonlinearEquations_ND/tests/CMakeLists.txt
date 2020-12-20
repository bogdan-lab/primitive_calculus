add_executable(lapack_test
                lapack_tests.c)

target_link_directories(lapack_test PRIVATE /usr/lib/)
target_include_directories(lapack_test PRIVATE /usr/include/)
target_compile_options(lapack_test PRIVATE -Wall -Wextra -Wshadow -Werror)
target_link_libraries(lapack_test PRIVATE lapacke m gslcblas)

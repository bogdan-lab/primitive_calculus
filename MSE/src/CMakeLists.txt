add_executable( testing_mse
         main.c
         mse_utils.c
)

target_compile_options(testing_mse PRIVATE -Wall -Wextra -Wshadow -Werror)
target_link_libraries(testing_mse PUBLIC m lapacke gslcblas)

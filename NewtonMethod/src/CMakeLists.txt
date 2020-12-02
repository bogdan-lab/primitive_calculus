add_library(newton_method SHARED
		newton_1D.c)

target_include_directories(newton_method PUBLIC ../include)

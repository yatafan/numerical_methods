cmake_minimum_required(VERSION 2.8)

project(lab3)

set(INC include/equations.h
		include/vector.h)

set(SRC src/equations.c
		src/vector.c)

add_executable(lab3 main.c ${INC} ${SRC})

target_link_libraries(lab3 m)

cmake_minimum_required(VERSION 2.8)

project(lab2)

set(INC include/matrix.h
		include/vector.h
		include/slau.h)

set(SRC src/matrix.c
		src/vector.c
		src/slau.c)

add_executable(lab2 main.c ${INC} ${SRC})

target_link_libraries(lab2 m)

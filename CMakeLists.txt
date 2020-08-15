cmake_minimum_required(VERSION 3.16.5)

project(linearalgebra C)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror -Wextra -Wstrict-prototypes -pedantic -pedantic-errors")

include_directories("./include")

set(SOURCE_FILES
        ./src/driver.c
        ./src/matrix.c
        ./src/vector.c
        ./src/equations.c
)

add_executable(linearalgebra ${SOURCE_FILES})
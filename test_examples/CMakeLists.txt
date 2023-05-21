cmake_minimum_required(VERSION 3.8)

set(This examples)

set(Source test.cpp)


add_executable(${This} ${Source})

target_link_libraries(${This} PUBLIC gtest_main)

add_test(NAME ${This} COMMAND ${This})

 
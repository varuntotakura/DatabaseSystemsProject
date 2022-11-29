cmake_minimum_required(VERSION 3.0)


set(CMAKE_BUID_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

project(
	DatabaseSystemsProject
	VERSION 1.0.0
	DESCRIPTION "Application"
	LANGUAGES CXX
)

add_executable(
	MainProgram
	main.cpp
)

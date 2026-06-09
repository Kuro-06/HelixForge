cmake_minimum_required(VERSION 3.15)
project(HelixForge)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 1. We added src/AlignmentEngine.cpp to our list of building materials
add_executable(HelixForge src/main.cpp src/AlignmentEngine.cpp src/FastaReader.cpp src/KmerIndex.cpp)

# 2. We tell the contractor: "Whenever our code asks for a blueprint, look in the 'include' folder!"
target_include_directories(HelixForge PRIVATE include)
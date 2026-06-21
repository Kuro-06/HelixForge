cmake_minimum_required(VERSION 3.15)
project(HelixForge)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(HelixForge
  src/main.cpp
  src/AlignmentEngine.cpp 
  
  src/genome.cpp
  

  src/AlignmentWorkflow 
 )

target_include_directories(HelixForge PRIVATE include)
cmake_minimum_required(VERSION 3.22)
project(federated_inference_of_grns)

set(CMAKE_CXX_STANDARD 14)

include_directories(src/cpp_implementation)

add_executable(federated_inference_of_grns
        src/cpp_implementation/global_model.cpp
        src/cpp_implementation/global_model.h
        src/cpp_implementation/train_local.cpp
        src/cpp_implementation/train_local.h
        src/cpp_implementation/utils.cpp
        src/cpp_implementation/utils.h)

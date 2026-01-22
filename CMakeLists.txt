cmake_minimum_required(VERSION 3.16)
project(processRootFiles LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED CONFIG)

add_executable(processRootFiles
    getSensitivityMapsProcesses.cpp
)

target_link_libraries(processRootFiles
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::Hist
    ROOT::RDF
    ROOT::Thread
    ROOT::Gpad
)

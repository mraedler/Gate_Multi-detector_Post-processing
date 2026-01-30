cmake_minimum_required(VERSION 3.16)
project(Gate_Multi-detector_Post-processing LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED CONFIG)

# Source
add_library(core
        source/CASToRTools.cpp
        source/eventSelection.cpp
        source/utils.cpp
)

# Headers
target_include_directories(core
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Add root
target_link_libraries(core
        PUBLIC
        ROOT::Core
        ROOT::RIO
        ROOT::Tree
        ROOT::Hist
        ROOT::Thread
        ROOT::Gpad
)

add_executable(addCastorIDsToGateCoincidencesTree addCastorIDsToGateCoincidencesTree.cpp)
target_link_libraries(addCastorIDsToGateCoincidencesTree PRIVATE core)

add_executable(addTruenessToGateCoincidenceTree addTruenessToGateCoincidenceTree.cpp)
target_link_libraries(addTruenessToGateCoincidenceTree PRIVATE core)

add_executable(runPreselection runPreselection.cpp)
target_link_libraries(runPreselection PRIVATE core)

add_executable(groupCoincidences groupCoincidences.cpp)
target_link_libraries(groupCoincidences PRIVATE core)

add_executable(getSensitivityMap getSensitivityMap.cpp)
target_link_libraries(getSensitivityMap PRIVATE core)

add_executable(processRootFiles getSensitivityMapsProcesses.cpp)
target_link_libraries(processRootFiles PRIVATE core)


cmake_minimum_required(VERSION 3.26)
project(ESTPGeneticAlgorithmWOC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
add_subdirectory(libraries/glfw)


# Source files
set(SOURCES
        src/main.cpp
        src/algorithms/geneticAlgorithm.cpp
        src/algorithms/WOCAlgorithm.cpp
        src/utilities/randomWeightedGraph.cpp
        src/utilities/globalData.cpp
        src/utilities/common.cpp
        src/gui/graphVisualizer.cpp
        # ImGui sources
        libraries/imgui/imgui.cpp
        libraries/imgui/imgui_draw.cpp
        libraries/imgui/imgui_widgets.cpp
        libraries/imgui/imgui_tables.cpp
        libraries/imgui/imgui_demo.cpp
        libraries/imgui/backends/imgui_impl_glfw.cpp
        libraries/imgui/backends/imgui_impl_opengl3.cpp
        # Plot sources
        libraries/implot/implot.cpp
        libraries/implot/implot_items.cpp
        libraries/implot/implot_demo.cpp
)

# Find OpenGL package
find_package(OpenGL REQUIRED)

# Add executable target
add_executable(ESTPGeneticAlgorithmWOC ${SOURCES})

# Include directories for the target
include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}/libraries/imgui
        ${CMAKE_CURRENT_SOURCE_DIR}/libraries/imgui/backends
        ${CMAKE_CURRENT_SOURCE_DIR}/libraries/glfw/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src/utilities
        ${CMAKE_CURRENT_SOURCE_DIR}/libraries/implot
        ${CMAKE_CURRENT_SOURCE_DIR}/src/gui
        ${CMAKE_CURRENT_SOURCE_DIR}/src/algorithms
)

# Link libraries to the target
target_link_libraries(ESTPGeneticAlgorithmWOC
        glfw
        OpenGL::GL
)

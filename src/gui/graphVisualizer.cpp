// graphVisualizer.cpp
// Created by LuisF on 11/14/2024.

#include "graphVisualizer.h"
#include "globalData.h"
#include "geneticAlgorithm.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"
#include <GLFW/glfw3.h>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <numeric>

void graphVisualizer::displayPlot(SteinerTree steinerTree, SteinerTree wocTree) {
    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW.\n";
        return;
    }

    // Configure GLFW
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

    // Create window
    GLFWwindow* window = glfwCreateWindow(1280, 720, "Graph Visualization", NULL, NULL);
    if (!window) {
        std::cerr << "Failed to create GLFW window.\n";
        glfwTerminate();
        return;
    }
    glfwMakeContextCurrent(window);

    // Initialize ImGui and ImPlot
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGuiIO& io = ImGui::GetIO();

    // Setup ImGui style
    ImGui::StyleColorsDark();

    // Initialize ImGui backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    SteinerTree randomPath = steinerTree;

    // Prepare data for plotting vertices
    std::vector<double> vertexX;
    std::vector<double> vertexY;

    for (const auto& vertex : vertices) {
        vertexX.push_back(vertex[1]);
        vertexY.push_back(vertex[2]);
    }

    // Prepare data for plotting edges
    std::vector<std::pair<double, double>> edgePointsA;
    std::vector<std::pair<double, double>> edgePointsB;

    for (const auto& edge : edges) {
        double x1 = vertices[edge.vertexA][1];
        double y1 = vertices[edge.vertexA][2];
        double x2 = vertices[edge.vertexB][1];
        double y2 = vertices[edge.vertexB][2];
        edgePointsA.push_back({ x1, y1 });
        edgePointsB.push_back({ x2, y2 });
    }

    // Prepare data for Steiner Path
    std::vector<std::pair<double, double>> steinerEdgePointsA;
    std::vector<std::pair<double, double>> steinerEdgePointsB;
    for (size_t i = 0; i < randomPath.path.size() - 1; ++i) {
        int from = randomPath.path[i];
        int to = randomPath.path[i + 1];
        double x1 = vertices[from][1];
        double y1 = vertices[from][2];
        double x2 = vertices[to][1];
        double y2 = vertices[to][2];
        steinerEdgePointsA.push_back({ x1, y1 });
        steinerEdgePointsB.push_back({ x2, y2 });
    }

    // Prepare data for WOC Path
    std::vector<std::pair<double, double>> wocEdgePointsA;
    std::vector<std::pair<double, double>> wocEdgePointsB;
    for (size_t i = 0; i < wocTree.path.size() - 1; ++i) {
        int from = wocTree.path[i];
        int to = wocTree.path[i + 1];
        double x1 = vertices[from][1];
        double y1 = vertices[from][2];
        double x2 = vertices[to][1];
        double y2 = vertices[to][2];
        wocEdgePointsA.push_back({ x1, y1 });
        wocEdgePointsB.push_back({ x2, y2 });
    }

    // Prepare data for Average Distances over Generations
    std::vector<double> generations(avgDistanceOverGenerations.size());
    for (size_t i = 0; i < avgDistanceOverGenerations.size(); ++i) {
        generations[i] = static_cast<double>(i + 1); // Start from 1
    }

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        // Poll events
        glfwPollEvents();

        // Start new ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Create ImGui window
        ImGui::Begin("Graph Visualization");

        if (ImGui::BeginTabBar("TabBar")) {
            //-----------------------------------------------------------------------------
            // NEW TAB 1 - Graph
            //-----------------------------------------------------------------------------
            if (ImGui::BeginTabItem("Graph")) {
                // Separate terminal and normal vertices
                std::vector<double> terminalX, terminalY;
                std::vector<double> normalX, normalY;

                for (const auto& vertex : vertices) {
                    if (vertex[3] == 1.0) { // Terminal vertex
                        terminalX.push_back(vertex[1]);
                        terminalY.push_back(vertex[2]);
                    } else { // Normal vertex
                        normalX.push_back(vertex[1]);
                        normalY.push_back(vertex[2]);
                    }
                }

                // Begin plot
                if (ImPlot::BeginPlot("Graph", ImVec2(-1, -1), 0)) {
                    // Compute dynamic axis limits
                    std::vector<double> allX;
                    std::vector<double> allY;
                    for (const auto& point : edgePointsA) {
                        allX.push_back(point.first);
                        allY.push_back(point.second);
                    }
                    for (const auto& point : edgePointsB) {
                        allX.push_back(point.first);
                        allY.push_back(point.second);
                    }

                    double minX = *std::min_element(allX.begin(), allX.end());
                    double maxX = *std::max_element(allX.begin(), allX.end());
                    double minY = *std::min_element(allY.begin(), allY.end());
                    double maxY = *std::max_element(allY.begin(), allY.end());

                    // Add padding
                    double paddingX = (maxX - minX) * 0.05; // 5% padding
                    double paddingY = (maxY - minY) * 0.05; // 5% padding

                    // Set up axes with padding
                    ImPlot::SetupAxes("X Coordinate", "Y Coordinate", 0, 0);
                    ImPlot::SetupAxisLimits(ImAxis_X1, minX - paddingX, maxX + paddingX, ImGuiCond_Always);
                    ImPlot::SetupAxisLimits(ImAxis_Y1, minY - paddingY, maxY + paddingY, ImGuiCond_Always);

                    // Plot edges as lines
                    ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.4f, 0.6f, 0.9f, 1.0f)); // Blue color
                    for (size_t i = 0; i < edgePointsA.size(); ++i) {
                        double xs[2] = { edgePointsA[i].first, edgePointsB[i].first };
                        double ys[2] = { edgePointsA[i].second, edgePointsB[i].second };
                        ImPlot::PlotLine("Edges", xs, ys, 2);
                    }
                    ImPlot::PopStyleColor();

                    // Plot normal vertices as scatter points
                    ImPlot::PushStyleColor(ImPlotCol_MarkerFill, ImVec4(0.4f, 0.6f, 0.9f, 1.0f));  // Blue
                    ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(0.4f, 0.6f, 0.9f, 1.0f));  // Blue
                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 5.0f);
                    ImPlot::PlotScatter("Vertices", normalX.data(), normalY.data(), normalX.size());
                    ImPlot::PopStyleColor(2);

                    // Plot terminal vertices as scatter points
                    ImPlot::PushStyleColor(ImPlotCol_MarkerFill, ImVec4(0.8f, 0.2f, 0.2f, 1.0f));  // Red
                    ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(0.8f, 0.2f, 0.2f, 1.0f));  // Red
                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 5.0f);
                    ImPlot::PlotScatter("Terminal Vertices", terminalX.data(), terminalY.data(), terminalX.size());
                    ImPlot::PopStyleColor(2);

                    // Annotate all vertices with their indices
                    for (size_t i = 0; i < vertices.size(); ++i) {
                        double x = vertices[i][1];
                        double y = vertices[i][2];
                        // Adjust offset for better visibility
                        ImVec2 offset(5, 5); // Adjust as needed
                        // Set text color
                        ImVec4 textColor = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);  // Black
                        // Annotate the vertex with its index
                        ImPlot::Annotation(
                                x, y,                   // Data coordinates
                                textColor,              // Text color
                                offset,                 // Pixel offset from the data point
                                true,                   // Clamp annotation within plot area
                                "%zu", i                // Format string to display the vertex index
                        );
                    }

                    ImPlot::EndPlot();
                }
                ImGui::EndTabItem();
            }

            //-----------------------------------------------------------------------------
            // TAB 2 - Steiner Tree
            //-----------------------------------------------------------------------------
            if (ImGui::BeginTabItem("Steiner Tree")) {
                // Code is similar to previous "Graph" tab, but includes the Steiner Path
                // Separate terminal and normal vertices
                std::vector<double> terminalX, terminalY;
                std::vector<double> normalX, normalY;

                for (const auto& vertex : vertices) {
                    if (vertex[3] == 1.0) { // Terminal vertex
                        terminalX.push_back(vertex[1]);
                        terminalY.push_back(vertex[2]);
                    } else { // Normal vertex
                        normalX.push_back(vertex[1]);
                        normalY.push_back(vertex[2]);
                    }
                }

                // Begin plot
                if (ImPlot::BeginPlot("Steiner Tree", ImVec2(-1, -1), 0)) {
                    // Compute dynamic axis limits
                    std::vector<double> allX;
                    std::vector<double> allY;
                    for (const auto& point : edgePointsA) {
                        allX.push_back(point.first);
                        allY.push_back(point.second);
                    }
                    for (const auto& point : edgePointsB) {
                        allX.push_back(point.first);
                        allY.push_back(point.second);
                    }

                    double minX = *std::min_element(allX.begin(), allX.end());
                    double maxX = *std::max_element(allX.begin(), allX.end());
                    double minY = *std::min_element(allY.begin(), allY.end());
                    double maxY = *std::max_element(allY.begin(), allY.end());

                    // Add padding
                    double paddingX = (maxX - minX) * 0.05; // 5% padding
                    double paddingY = (maxY - minY) * 0.05; // 5% padding

                    // Set up axes with padding
                    ImPlot::SetupAxes("X Coordinate", "Y Coordinate", 0, 0);
                    ImPlot::SetupAxisLimits(ImAxis_X1, minX - paddingX, maxX + paddingX, ImGuiCond_Always);
                    ImPlot::SetupAxisLimits(ImAxis_Y1, minY - paddingY, maxY + paddingY, ImGuiCond_Always);

                    // Plot edges as lines
                    ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.4f, 0.6f, 0.9f, 1.0f));
                    for (size_t i = 0; i < edgePointsA.size(); ++i) {
                        double xs[2] = { edgePointsA[i].first, edgePointsB[i].first };
                        double ys[2] = { edgePointsA[i].second, edgePointsB[i].second };
                        ImPlot::PlotLine("Edges", xs, ys, 2);
                    }
                    ImPlot::PopStyleColor();

                    // Plot normal vertices as scatter points
                    ImPlot::PushStyleColor(ImPlotCol_MarkerFill, ImVec4(0.4f, 0.6f, 0.9f, 1.0f));  // Blue
                    ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(0.4f, 0.6f, 0.9f, 1.0f));  // Blue
                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 5.0f);
                    ImPlot::PlotScatter("Vertices", normalX.data(), normalY.data(), normalX.size());
                    ImPlot::PopStyleColor(2);

                    // Plot terminal vertices as scatter points
                    ImPlot::PushStyleColor(ImPlotCol_MarkerFill, ImVec4(0.8f, 0.2f, 0.2f, 1.0f));  // Red
                    ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(0.8f, 0.2f, 0.2f, 1.0f));  // Red
                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 5.0f);
                    ImPlot::PlotScatter("Terminal Vertices", terminalX.data(), terminalY.data(), terminalX.size());
                    ImPlot::PopStyleColor(2);

                    // Annotate all vertices with their indices
                    for (size_t i = 0; i < vertices.size(); ++i) {
                        double x = vertices[i][1];
                        double y = vertices[i][2];
                        // Adjust offset for better visibility
                        ImVec2 offset(5, 5); // Adjust as needed
                        // Set text color
                        ImVec4 textColor = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);  // Black
                        // Annotate the vertex with its index
                        ImPlot::Annotation(
                                x, y,                   // Data coordinates
                                textColor,              // Text color
                                offset,                 // Pixel offset from the data point
                                true,                   // Clamp annotation within plot area
                                "%zu", i                // Format string to display the vertex index
                        );
                    }

                    // Plot Steiner Path edges as red lines
                    if (!randomPath.path.empty()) {
                        ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.8f, 0.2f, 0.2f, 1.0f)); // Red color
                        for (size_t i = 0; i < steinerEdgePointsA.size(); ++i) {
                            double xs[2] = { steinerEdgePointsA[i].first, steinerEdgePointsB[i].first };
                            double ys[2] = { steinerEdgePointsA[i].second, steinerEdgePointsB[i].second };
                            ImPlot::PlotLine("Steiner Path", xs, ys, 2);
                        }
                        ImPlot::PopStyleColor();
                    }

                    ImPlot::EndPlot();
                }
                ImGui::EndTabItem();
            }
        //-----------------------------------------------------------------------------
        // TAB 3 - WOC Tree
        //-----------------------------------------------------------------------------
            if (ImGui::BeginTabItem("WOC Tree")) {
                // Separate terminal and normal vertices
                std::vector<double> terminalX, terminalY;
                std::vector<double> normalX, normalY;

                for (const auto& vertex : vertices) {
                    if (vertex[3] == 1.0) { // Terminal vertex
                        terminalX.push_back(vertex[1]);
                        terminalY.push_back(vertex[2]);
                    } else { // Normal vertex
                        normalX.push_back(vertex[1]);
                        normalY.push_back(vertex[2]);
                    }
                }

                // Begin plot
                if (ImPlot::BeginPlot("WOC Tree", ImVec2(-1, -1), 0)) {
                    // Compute dynamic axis limits
                    std::vector<double> allX;
                    std::vector<double> allY;
                    for (const auto& point : edgePointsA) {
                        allX.push_back(point.first);
                        allY.push_back(point.second);
                    }
                    for (const auto& point : edgePointsB) {
                        allX.push_back(point.first);
                        allY.push_back(point.second);
                    }

                    double minX = *std::min_element(allX.begin(), allX.end());
                    double maxX = *std::max_element(allX.begin(), allX.end());
                    double minY = *std::min_element(allY.begin(), allY.end());
                    double maxY = *std::max_element(allY.begin(), allY.end());

                    // Add padding
                    double paddingX = (maxX - minX) * 0.05; // 5% padding
                    double paddingY = (maxY - minY) * 0.05; // 5% padding

                    // Set up axes with padding
                    ImPlot::SetupAxes("X Coordinate", "Y Coordinate", 0, 0);
                    ImPlot::SetupAxisLimits(ImAxis_X1, minX - paddingX, maxX + paddingX, ImGuiCond_Always);
                    ImPlot::SetupAxisLimits(ImAxis_Y1, minY - paddingY, maxY + paddingY, ImGuiCond_Always);

                    // Plot edges as lines
                    ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.4f, 0.6f, 0.9f, 1.0f)); // Blue color
                    for (size_t i = 0; i < edgePointsA.size(); ++i) {
                        double xs[2] = { edgePointsA[i].first, edgePointsB[i].first };
                        double ys[2] = { edgePointsA[i].second, edgePointsB[i].second };
                        ImPlot::PlotLine("Edges", xs, ys, 2);
                    }
                    ImPlot::PopStyleColor();

                    // Plot normal vertices as scatter points
                    ImPlot::PushStyleColor(ImPlotCol_MarkerFill, ImVec4(0.4f, 0.6f, 0.9f, 1.0f));  // Blue
                    ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(0.4f, 0.6f, 0.9f, 1.0f));  // Blue
                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 5.0f);
                    ImPlot::PlotScatter("Vertices", normalX.data(), normalY.data(), normalX.size());
                    ImPlot::PopStyleColor(2);

                    // Plot terminal vertices as scatter points
                    ImPlot::PushStyleColor(ImPlotCol_MarkerFill, ImVec4(0.8f, 0.2f, 0.2f, 1.0f));  // Red
                    ImPlot::PushStyleColor(ImPlotCol_MarkerOutline, ImVec4(0.8f, 0.2f, 0.2f, 1.0f));  // Red
                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 5.0f);
                    ImPlot::PlotScatter("Terminal Vertices", terminalX.data(), terminalY.data(), terminalX.size());
                    ImPlot::PopStyleColor(2);

                    // Annotate all vertices with their indices
                    for (size_t i = 0; i < vertices.size(); ++i) {
                        double x = vertices[i][1];
                        double y = vertices[i][2];
                        // Adjust offset for better visibility
                        ImVec2 offset(5, 5); // Adjust as needed
                        // Set text color
                        ImVec4 textColor = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);  // Black
                        // Annotate the vertex with its index
                        ImPlot::Annotation(
                                x, y,                   // Data coordinates
                                textColor,              // Text color
                                offset,                 // Pixel offset from the data point
                                true,                   // Clamp annotation within plot area
                                "%zu", i                // Format string to display the vertex index
                        );
                    }

                    // Plot WOC Path edges as green lines
                    if (!wocTree.path.empty()) {
                        ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(0.2f, 0.8f, 0.2f, 1.0f)); // Green color
                        for (size_t i = 0; i < wocEdgePointsA.size(); ++i) {
                            double xs[2] = { wocEdgePointsA[i].first, wocEdgePointsB[i].first };
                            double ys[2] = { wocEdgePointsA[i].second, wocEdgePointsB[i].second };
                            ImPlot::PlotLine("WOC Path", xs, ys, 2);
                        }
                        ImPlot::PopStyleColor();
                    }

                    ImPlot::EndPlot();
                }
                ImGui::EndTabItem();
            }

            //-----------------------------------------------------------------------------
            // TAB 4 - Average Distances
            //-----------------------------------------------------------------------------
            if (ImGui::BeginTabItem("Average Distances")) {
                // Begin plot
                if (ImPlot::BeginPlot("Average Distance Over Generations", ImVec2(-1, -1), 0)) {
                    // Set up axes
                    ImPlot::SetupAxes("Generation", "Average Total Distance", 0, 0);

                    // Compute axis limits with padding
                    double min_gen = 1.0;
                    double max_gen = static_cast<double>(generations.size());
                    double min_dist = *std::min_element(avgDistanceOverGenerations.begin(), avgDistanceOverGenerations.end());
                    double max_dist = *std::max_element(avgDistanceOverGenerations.begin(), avgDistanceOverGenerations.end());

                    double padding_gen = (max_gen - min_gen) * 0.05; // 5% padding
                    double padding_dist = (max_dist - min_dist) * 0.05; // 5% padding

                    // Set axis limits with padding
                    ImPlot::SetupAxisLimits(ImAxis_X1, min_gen - padding_gen, max_gen + padding_gen, ImGuiCond_Always);
                    ImPlot::SetupAxisLimits(ImAxis_Y1, min_dist - padding_dist, max_dist + padding_dist, ImGuiCond_Always);

                    // Add legend
                    ImPlot::SetupLegend(ImPlotLocation_NorthEast);

                    // Plot the average distance over generations
                    ImPlot::PlotLine("Avg Distance", generations.data(), avgDistanceOverGenerations.data(), avgDistanceOverGenerations.size());

                    ImPlot::EndPlot();
                }
                ImGui::EndTabItem();
            }

            ImGui::EndTabBar();
        }
        ImGui::End();

        // Rendering
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
#ifdef __APPLE__
        glViewport(0, 0, display_w, display_h);
#else
        glViewport(0, 0, display_w, display_h);
#endif
        glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
        glClear(GL_COLOR_BUFFER_BIT);

        // Render ImGui draw data
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        // Swap buffers
        glfwSwapBuffers(window);
    }

    // Cleanup
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    glfwDestroyWindow(window);
    glfwTerminate();
}

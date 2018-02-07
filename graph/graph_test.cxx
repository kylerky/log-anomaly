#include "graph.hxx"
#include "gtest/gtest.h"
#include <iostream>

using namespace GPSOINN;

TEST(basic, init) {
    Digraph<unsigned, int> graph;
    UndirectedGraph<unsigned, int> ugraph;
}

TEST(Digraph, insert_vertex) {
    Digraph<int, int> graph;
    for (int i = 0; i != 4; ++i) {
        graph.insert_vertex(i);
    }
    ASSERT_EQ(4, graph.vertex_count());
    {
        auto iter = graph.cbegin();
        for (unsigned i = 0; i != 4; ++i, ++iter) {
            EXPECT_EQ(i, iter->value());
        }
    }
}

TEST(Digraph, insert_erase_edge) {
    Digraph<int, int> graph;
    for (int i = 0; i != 5; ++i) {
        graph.insert_vertex(i);
    }
    graph.insert_edge(0, 1, 1);
    graph.insert_edge(1, 2, 2);
    graph.insert_edge(1, 3, 4);
    graph.insert_edge(1, 3, 4);
    graph.insert_edge(3, 4, 1);
    graph.insert_edge(4, 2, 1);
    // graph.erase_after_edge(ver, edge);
    ASSERT_EQ(5, graph.vertex_count());
    for (auto ver : graph)
        for (auto edge : ver)
            std::cout << ver.value() << " - " << edge.head << std::endl;
    graph.erase_vertex(1);
    ASSERT_EQ(4, graph.vertex_count());
    std::cout << std::endl;
    for (auto ver : graph)
        for (auto edge : ver)
            std::cout << ver.value() << " - " << edge.head << std::endl;
}

TEST(UndirenctedGraph, insert_erase_vertex) {
    UndirectedGraph<unsigned, int> graph;
    for (int i = 0; i != 5; ++i) {
        graph.insert_vertex(i);
    }
    graph.insert_edge(0, 1, 1);
    graph.insert_edge(1, 2, 2);
    graph.insert_edge(1, 3, 4);
    graph.insert_edge(1, 3, 5);
    graph.insert_edge(3, 4, 1);
    graph.insert_edge(4, 2, 1);
    ASSERT_EQ(5, graph.vertex_count());
    for (auto ver : graph)
        for (auto edge : ver)
            std::cout << ver.value() << " - " << edge.head << " " << edge.weight
                      << std::endl;
    std::cout << std::endl;
    auto v1 = graph.get_vertex_iterator(1);
    for (auto iter = v1->cbegin(), pre = v1->cbefore_begin();
         iter != v1->cend(); ++iter, ++pre) {
        if (iter->head == 3) {
            graph.erase_after_edge(v1, pre);
            break;
        }
    }
    for (auto ver : graph)
        for (auto edge : ver)
            std::cout << ver.value() << " - " << edge.head << " " << edge.weight
                      << std::endl;
}

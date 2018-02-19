cmake_minimum_required (VERSION 2.8.2)
project (infomap)

include_directories(${PROJECT_SOURCE_DIR})

add_library(infomap 
                    infomap/FlowNetwork.cpp
                    infomap/InfomapBase.cpp
                    infomap/InfomapContext.cpp
                    infomap/MemFlowNetwork.cpp
                    infomap/MemNetwork.cpp
                    infomap/MemoryNetworkAdapter.cpp
                    infomap/MultiplexNetwork.cpp
                    infomap/NetworkAdapter.cpp
                    infomap/Network.cpp
                    infomap/Node.cpp
                    infomap/TreeData.cpp

                    io/BipartiteClusterReader.cpp
                    io/ClusterReader.cpp
                    io/HierarchicalNetwork.cpp
                    io/ProgramInterface.cpp
                    io/TreeDataWriter.cpp
                    io/version.cpp

                    utils/FileURI.cpp
                    utils/Logger.cpp

                    Infomap.cpp
)

target_compile_options(infomap PRIVATE "-w")

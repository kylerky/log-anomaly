#ifndef LOG_ANOMALY_GPSOINN_HXX
#define LOG_ANOMALY_GPSOINN_HXX

#include "graph/graph.hxx"

#include "cereal/access.hpp"

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <algorithm>
#include <array>
#include <cmath>
#include <random>

namespace LogAnomaly {

template <unsigned dimension> class GPNet {
  private:
    typedef Eigen::Matrix<double, dimension, 1> NodeVector;
    struct Node {
        NodeVector vector;
        size_t win_count;
    };
    typedef std::array<double, dimension> array;
    typedef UndirectedGraph<Node, unsigned, std::less<NodeVector>,
                            Eigen::aligned_allocator<NodeVector>>
        UGraph;

  public:
    GPNet(unsigned lambda = 20000, unsigned age_max = 50, unsigned k = 1,
          double sigma_2 = 1e-6);
    ~GPNet() {}

    void train(array &data);
    double predict(array &data);

  private:
    UGraph m_graph;
    double m_sigma_2;
    unsigned m_local_opt_coeff = 100;
    unsigned m_age_max;
    unsigned m_lambda;
    unsigned m_k;
    unsigned m_cycles = 0;
    const double m_const_coeff;

    std::pair<double, double> threshold(size_t index,
                                        const NodeVector &x_vector);
}; // class GPNet

} // namespace LogAnomaly

namespace LogAnomaly {

template <unsigned dimension>
GPNet<dimension>::GPNet(unsigned lambda, unsigned age_max, unsigned k,
                        double sigma_2)
    : m_lambda(lambda), m_age_max(age_max), m_k(k), m_sigma_2(sigma_2),
      m_const_coeff(std::sqrt(std::pow(2 * M_PI, dimension))) {
    m_graph.insert_vertex(Node{.vector = NodeVector::Random(), .win_count = 0});
    m_graph.insert_vertex(Node{.vector = NodeVector::Random(), .win_count = 0});
}

template <unsigned dimension> void GPNet<dimension>::train(array &data) {
    using namespace Eigen;
    using NodeM = Matrix<double, dimension, Dynamic>;

    if (m_graph.vertex_count() < 2) {
        m_graph.insert_vertex(
            Node{.vector = NodeVector::Random(), .win_count = 0});
        m_graph.insert_vertex(
            Node{.vector = NodeVector::Random(), .win_count = 0});
    }
    std::random_device rand_dev;
    std::srand(rand_dev());

    typename NodeM::Index min_index, min2_index;
    Map<NodeVector> input(data.data());

    {
        // construct the data matrix
        NodeM nodes(dimension, m_graph.vertex_count());

        size_t col = 0;
        for (auto &node : m_graph) {
            nodes.col(col) = node.value().vector;
            ++col;
        }

        // find the winner and the second winner
        RowVectorXd distances =
            (nodes.colwise() - input).colwise().squaredNorm();
        distances.minCoeff(&min_index);

        if (min_index == 0) {
            distances(0) = distances(1);
            distances.minCoeff(&min2_index);
            if (min2_index == 0)
                min2_index = 1;
        } else {
            distances(min_index) = distances(min_index - 1);
            distances.minCoeff(&min2_index);
            if (min2_index == min_index)
                ++min2_index;
        }
    }

    auto min_iter = m_graph.begin();
    auto min2_iter = m_graph.begin();
    for (size_t i = 0; i != min_index; ++i)
        ++min_iter;
    for (size_t i = 0; i != min2_index; ++i)
        ++min2_iter;

    auto[threshold1, prob1] = threshold(min_index, input);
    auto[threshold2, prob2] = threshold(min2_index, input);

    if ((prob1 > prob2 && prob1 > threshold1) || prob2 > threshold2) {
        m_graph.insert_edge(min_iter, min2_iter, 0);

        auto &win_count = min_iter->value().win_count;
        NodeVector &winner_vec = min_iter->value().vector;
        ++win_count;
        winner_vec.array() += ((input - winner_vec).array() / (win_count + 1));

        for (auto &edge : *min_iter) {
            NodeVector &node_vec = m_graph[edge.head].value().vector;
            node_vec.array() += ((input - node_vec).array() /
                                 m_local_opt_coeff / (win_count + 1));
            ++edge.weight;
        }

        for (auto niter = m_graph.cbegin(); niter != m_graph.cend(); ++niter) {
            for (auto edge = niter->cbegin(), pre = niter->cbefore_begin();
                 edge != niter->cend();) {
                if (edge->weight > m_age_max) {
                    edge = m_graph.erase_after_edge(niter, pre);
                } else {
                    ++edge;
                    ++pre;
                }
            }
        }
    } else {
        m_graph.insert_vertex(Node{.vector = input, .win_count = 0});
    }

    ++m_cycles;
    if (m_cycles == m_lambda) {
        m_cycles = 0;
        for (auto iter = m_graph.cbegin(); iter != m_graph.cend();) {
            if (std::distance(iter->cbegin(), iter->cend()) < m_k)
                iter = m_graph.erase_vertex(iter);
            else
                ++iter;
        }
    }
} // namespace LogAnomaly

template <unsigned dimension>
std::pair<double, double>
GPNet<dimension>::threshold(size_t index, const NodeVector &x_vector) {
    using namespace Eigen;

    Matrix<double, dimension, dimension> local_cov;
    local_cov.setZero();
    auto iter = m_graph.begin();
    for (size_t i = 0; i != index; ++i)
        ++iter;

    NodeVector &winner_vec = iter->value().vector;
    size_t win_sum = 0;
    size_t neighbour_cnt = 0;
    for (auto edge : *iter) {
        NodeVector &node_vec = m_graph[edge.head].value().vector;
        local_cov.array() +=
            ((node_vec - winner_vec) * (node_vec - winner_vec).transpose())
                .array() *
            edge.weight;
        win_sum += edge.weight;
        ++neighbour_cnt;
    }
    if (win_sum)
        local_cov.array() /= win_sum;
    local_cov += MatrixXd::Identity(dimension, dimension) * m_sigma_2;

    // if (neighbour_cnt > dimension) {
    //     Matrix<bool, dimension, 1> comp;
    //     SelfAdjointEigenSolver<Matrix<double, dimension, dimension>>
    //     eigens(local_cov);
    //     comp = ((eigens.eigenvalues()).array() >=
    //     m_sigma_2)
    //                .template cast<bool>();

    //     unsigned count = comp.count();
    //     Matrix<double, dimension, Dynamic> p_components(dimension, count);
    //     VectorXd eigenvals(count);

    //     {
    //         unsigned cnt = 0;
    //         for (unsigned i = 0; i != dimension; ++i) {
    //             if (comp(i)) {
    //                 p_components.col(cnt) =
    //                 eigens.eigenvectors().col(i); eigenvals(cnt) =
    //                 eigens.eigenvalues()(i);
    //                 ++cnt;
    //             }
    //         }
    //     }
    //     local_cov =
    //         p_components * eigenvals.asDiagonal() * p_components.transpose()
    //         + MatrixXd::Identity(dimension, dimension) * m_sigma_2;
    // }

    SelfAdjointEigenSolver<Matrix<double, dimension, dimension>> eigens(
        local_cov);
    double determinant_sqrt = std::sqrt(std::abs(eigens.eigenvalues().prod()));
    double threshold = 1;
    Matrix<double, dimension, dimension> cov_inv =
        eigens.eigenvectors() *
        eigens.eigenvalues().array().inverse().matrix().asDiagonal() *
        eigens.eigenvectors().transpose();
    if (neighbour_cnt) {
        for (auto edge : *iter) {
            NodeVector &node_vec = m_graph[edge.head].value().vector;
            double prob = std::exp(-((x_vector - node_vec).transpose() *
                                     cov_inv * (x_vector - node_vec))(0) /
                                   2) /
                          m_const_coeff / determinant_sqrt;

            if (prob < threshold)
                threshold = prob;
        }
    } else
        threshold = 0.55;

    double xprob = std::exp(-((x_vector - winner_vec).transpose() * cov_inv *
                              (x_vector - winner_vec))(0) /
                            2) /
                   m_const_coeff / determinant_sqrt;
    return {
        threshold,
        xprob,
    };
}
template <unsigned dimension> double GPNet<dimension>::predict(array &data) {
    using namespace Eigen;

    Map<NodeVector> input(data.data());

    double prob = 0;
    size_t wins = 0;
    for (auto &node : m_graph) {
        Matrix<double, dimension, dimension> local_cov;
        local_cov.setZero();
        size_t win_sum = 0;
        for (auto edge : node) {
            NodeVector &edge_vec = m_graph[edge.head].value().vector;
            NodeVector &node_vec = node.value().vector;
            local_cov.array() +=
                ((edge_vec - node_vec) * (edge_vec - node_vec).transpose())
                    .array() *
                edge.weight;
            win_sum += edge.weight;
        }
        if (win_sum)
            local_cov.array() /= win_sum;
        local_cov += MatrixXd::Identity(dimension, dimension) * m_sigma_2;

        SelfAdjointEigenSolver<Matrix<double, dimension, dimension>> eigens(
            local_cov);
        double determinant_sqrt = std::sqrt(eigens.eigenvalues().prod());
        Matrix<double, dimension, dimension> cov_inv =
            eigens.eigenvectors() *
            eigens.eigenvalues().array().inverse().matrix().asDiagonal() *
            eigens.eigenvectors().transpose();

        NodeVector &node_vec = node.value().vector;
        double xprob = std::exp(-((input - node_vec).transpose() * cov_inv *
                                  (input - node_vec))(0) /
                                2) /
                       m_const_coeff / determinant_sqrt;
        prob += xprob * node.value().win_count;
        wins += node.value().win_count;
    }

    if (wins)
        prob /= wins;
    return prob;
}

// specialization
template <> class GPNet<0> {
  public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> NodeVector;
    struct Node {
        NodeVector vector;
        size_t win_count;
    };

  private:
    typedef std::vector<double> vector_d;
    typedef UndirectedGraph<Node, unsigned, std::less<NodeVector>,
                            Eigen::aligned_allocator<NodeVector>>
        UGraph;

    template <typename Archive>
    friend void serialize(Archive &archive, Node &node) {
        archive(node.vector, node.win_count);
    }

    template <typename Archive>
    friend void serialize(Archive &archive, GPNet<0> &net) {
        archive(net.m_graph, net.m_sigma_2, net.m_local_opt_coeff,
                net.m_age_max, net.m_lambda, net.m_k, net.m_cycles,
                net.m_const_coeff, net.m_dimension);
    }

  public:
    GPNet(unsigned dimension = 1, unsigned lambda = 20000,
          unsigned age_max = 50, unsigned k = 1, double sigma_2 = 1e-6);
    ~GPNet() {}

    void train(vector_d &data);
    double predict(vector_d &data);

  private:
    UGraph m_graph;
    double m_sigma_2;
    unsigned m_local_opt_coeff = 100;
    unsigned m_age_max;
    unsigned m_lambda;
    unsigned m_k;
    unsigned m_cycles = 0;
    double m_const_coeff;
    unsigned m_dimension;

    std::pair<double, double> threshold(size_t index,
                                        const NodeVector &x_vector);
}; // class GPNet

GPNet<0>::GPNet(unsigned dimension, unsigned lambda, unsigned age_max,
                unsigned k, double sigma_2)
    : m_lambda(lambda), m_age_max(age_max), m_k(k), m_sigma_2(sigma_2),
      m_const_coeff(std::sqrt(std::pow(2 * M_PI, dimension))),
      m_dimension(dimension) {
    m_graph.insert_vertex(
        Node{.vector = NodeVector::Random(m_dimension, 1), .win_count = 0});
    m_graph.insert_vertex(
        Node{.vector = NodeVector::Random(m_dimension, 1), .win_count = 0});
}

void GPNet<0>::train(vector_d &data) {
    using namespace Eigen;
    using NodeM = Matrix<double, Dynamic, Dynamic>;

    if (m_graph.vertex_count() < 2) {
        m_graph.insert_vertex(
            Node{.vector = NodeVector::Random(m_dimension, 1), .win_count = 0});
        m_graph.insert_vertex(
            Node{.vector = NodeVector::Random(m_dimension, 1), .win_count = 0});
    }
    std::random_device rand_dev;
    std::srand(rand_dev());

    typename NodeM::Index min_index, min2_index;
    Map<NodeVector> input(data.data(), m_dimension);

    {
        // construct the data matrix
        NodeM nodes(m_dimension, m_graph.vertex_count());

        size_t col = 0;
        for (auto &node : m_graph) {
            nodes.col(col) = node.value().vector;
            ++col;
        }

        // find the winner and the second winner
        RowVectorXd distances =
            (nodes.colwise() - input).colwise().squaredNorm();
        distances.minCoeff(&min_index);

        if (min_index == 0) {
            distances(0) = distances(1);
            distances.minCoeff(&min2_index);
            if (min2_index == 0)
                min2_index = 1;
        } else {
            distances(min_index) = distances(min_index - 1);
            distances.minCoeff(&min2_index);
            if (min2_index == min_index)
                ++min2_index;
        }
    }

    auto min_iter = m_graph.begin();
    auto min2_iter = m_graph.begin();
    for (size_t i = 0; i != min_index; ++i)
        ++min_iter;
    for (size_t i = 0; i != min2_index; ++i)
        ++min2_iter;

    auto[threshold1, prob1] = threshold(min_index, input);
    auto[threshold2, prob2] = threshold(min2_index, input);

    if ((prob1 > prob2 && prob1 > threshold1) || prob2 > threshold2) {
        m_graph.insert_edge(min_iter, min2_iter, 0);

        auto &win_count = min_iter->value().win_count;
        NodeVector &winner_vec = min_iter->value().vector;
        ++win_count;
        winner_vec.array() += ((input - winner_vec).array() / (win_count + 1));

        for (auto &edge : *min_iter) {
            NodeVector &node_vec = m_graph[edge.head].value().vector;
            node_vec.array() += ((input - node_vec).array() /
                                 m_local_opt_coeff / (win_count + 1));
            ++edge.weight;
        }

        for (auto niter = m_graph.cbegin(); niter != m_graph.cend(); ++niter) {
            for (auto edge = niter->cbegin(), pre = niter->cbefore_begin();
                 edge != niter->cend();) {
                if (edge->weight > m_age_max) {
                    edge = m_graph.erase_after_edge(niter, pre);
                } else {
                    ++edge;
                    ++pre;
                }
            }
        }
    } else {
        m_graph.insert_vertex(Node{.vector = input, .win_count = 0});
    }

    ++m_cycles;
    if (m_cycles == m_lambda) {
        m_cycles = 0;
        for (auto iter = m_graph.cbegin(); iter != m_graph.cend();) {
            if (std::distance(iter->cbegin(), iter->cend()) < m_k)
                iter = m_graph.erase_vertex(iter);
            else
                ++iter;
        }
    }
} // namespace LogAnomaly

std::pair<double, double> GPNet<0>::threshold(size_t index,
                                              const NodeVector &x_vector) {
    using namespace Eigen;

    Matrix<double, Dynamic, Dynamic> local_cov(m_dimension, m_dimension);
    local_cov.setZero();
    auto iter = m_graph.begin();
    for (size_t i = 0; i != index; ++i)
        ++iter;

    NodeVector &winner_vec = iter->value().vector;
    size_t win_sum = 0;
    size_t neighbour_cnt = 0;
    for (auto edge : *iter) {
        NodeVector &node_vec = m_graph[edge.head].value().vector;
        local_cov.array() +=
            ((node_vec - winner_vec) * (node_vec - winner_vec).transpose())
                .array() *
            edge.weight;
        win_sum += edge.weight;
        ++neighbour_cnt;
    }
    if (win_sum)
        local_cov.array() /= win_sum;
    local_cov += MatrixXd::Identity(m_dimension, m_dimension) * m_sigma_2;

    // if (neighbour_cnt > m_dimension) {
    //     Matrix<bool, Dynamic, 1> comp(m_dimension, 1);
    //     SelfAdjointEigenSolver<Matrix<double, Dynamic, Dynamic>>
    //     eigens(local_cov);
    //     comp = ((eigens.eigenvalues()).array() >=
    //     m_sigma_2)
    //                .template cast<bool>();

    //     unsigned count = comp.count();
    //     Matrix<double, Dynamic, Dynamic> p_components(m_dimension, count);
    //     VectorXd eigenvals(count);

    //     {
    //         unsigned cnt = 0;
    //         for (unsigned i = 0; i != m_dimension; ++i) {
    //             if (comp(i)) {
    //                 p_components.col(cnt) =
    //                 eigens.eigenvectors().col(i); eigenvals(cnt) =
    //                 eigens.eigenvalues()(i);
    //                 ++cnt;
    //             }
    //         }
    //     }
    //     local_cov =
    //         p_components * eigenvals.asDiagonal() * p_components.transpose()
    //         + MatrixXd::Identity(m_dimension, m_dimension) * m_sigma_2;
    // }

    SelfAdjointEigenSolver<Matrix<double, Dynamic, Dynamic>> eigens(local_cov);
    double determinant_sqrt = std::sqrt(eigens.eigenvalues().prod());
    double threshold = 1;
    Matrix<double, Dynamic, Dynamic> cov_inv =
        eigens.eigenvectors() *
        eigens.eigenvalues().array().inverse().matrix().asDiagonal() *
        eigens.eigenvectors().transpose();
    if (neighbour_cnt) {
        for (auto edge : *iter) {
            NodeVector &node_vec = m_graph[edge.head].value().vector;
            double prob = std::exp(-((x_vector - node_vec).transpose() *
                                     cov_inv * (x_vector - node_vec))(0) /
                                   2) /
                          m_const_coeff / determinant_sqrt;

            if (prob < threshold)
                threshold = prob;
        }
    } else
        threshold = 0.55;

    double xprob = std::exp(-((x_vector - winner_vec).transpose() * cov_inv *
                              (x_vector - winner_vec))(0) /
                            2) /
                   m_const_coeff / determinant_sqrt;
    return {
        threshold,
        xprob,
    };
}
double GPNet<0>::predict(vector_d &data) {
    using namespace Eigen;

    Map<NodeVector> input(data.data(), m_dimension);

    double prob = 0;
    size_t wins = 0;
    for (auto &node : m_graph) {
        Matrix<double, Dynamic, Dynamic> local_cov(m_dimension, m_dimension);
        local_cov.setZero();
        size_t win_sum = 0;
        for (auto edge : node) {
            NodeVector &edge_vec = m_graph[edge.head].value().vector;
            NodeVector &node_vec = node.value().vector;
            local_cov.array() +=
                ((edge_vec - node_vec) * (edge_vec - node_vec).transpose())
                    .array() *
                edge.weight;
            win_sum += edge.weight;
        }
        if (win_sum)
            local_cov.array() /= win_sum;
        local_cov += MatrixXd::Identity(m_dimension, m_dimension) * m_sigma_2;

        SelfAdjointEigenSolver<Matrix<double, Dynamic, Dynamic>> eigens(
            local_cov);
        double determinant_sqrt =
            std::sqrt(std::abs(eigens.eigenvalues().prod()));
        Matrix<double, Dynamic, Dynamic> cov_inv =
            eigens.eigenvectors() *
            eigens.eigenvalues().array().inverse().matrix().asDiagonal() *
            eigens.eigenvectors().transpose();

        NodeVector &node_vec = node.value().vector;
        double xprob = std::exp(-((input - node_vec).transpose() * cov_inv *
                                  (input - node_vec))(0) /
                                2) /
                       m_const_coeff / determinant_sqrt;
        prob += xprob * node.value().win_count;
        wins += node.value().win_count;
    }

    if (wins)
        prob /= wins;
    return prob;
}
} // namespace LogAnomaly

namespace cereal {

template <typename Archive>
void save(Archive &archive,
          Eigen::Matrix<double, Eigen::Dynamic, 1> const &vec) {
    archive(vec.size());
    for (unsigned i = 0; i != vec.size(); ++i)
        archive(vec(i));
}
template <typename Archive>
void load(Archive &archive, Eigen::Matrix<double, Eigen::Dynamic, 1> &vec) {
    Eigen::Index sz;
    archive(sz);
    vec.resize(sz);
    for (unsigned i = 0; i != vec.size(); ++i) {
        double elem;
        archive(elem);
        vec(i) = elem;
    }
}

} // namespace cereal

#endif // LOG_ANOMALY_GPSOINN_HXX

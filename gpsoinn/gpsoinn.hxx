#ifndef LOG_ANOMALY_GPSOINN_HXX
#define LOG_ANOMALY_GPSOINN_HXX

#include "graph/graph.hxx"
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <algorithm>
#include <array>
#include <cmath>
#include <random>

namespace LogAnomaly {

template <unsigned dimension> class GPNet {
    static_assert(dimension > 0, "Dimension must be above 0");

  private:
    typedef Eigen::Matrix<float, dimension, 1> NodeVector;
    struct Node {
        NodeVector vector;
        size_t win_count;
    };
    typedef std::array<float, dimension> array;
    typedef UndirectedGraph<Node, unsigned, std::less<NodeVector>,
                            Eigen::aligned_allocator<NodeVector>>
        UGraph;

  public:
    GPNet(unsigned lambda = 20000, unsigned age_max = 50, unsigned k = 1,
          float sigma_2 = 0.001);
    ~GPNet() {}

    void train(array &data);
    float predict(array &data);

  private:
    UGraph m_graph;
    float m_sigma_2;
    unsigned m_local_opt_coeff = 100;
    unsigned m_age_max;
    unsigned m_lambda;
    unsigned m_k;
    unsigned m_cycles = 0;

    std::pair<float, float> threshold(size_t index, const NodeVector &x_vector);
}; // class GPNet

} // namespace LogAnomaly

namespace LogAnomaly {

template <unsigned dimension>
GPNet<dimension>::GPNet(unsigned lambda, unsigned age_max, unsigned k,
                        float sigma_2)
    : m_lambda(lambda), m_age_max(age_max), m_k(k), m_sigma_2(sigma_2) {
    m_graph.insert_vertex(Node{.vector = NodeVector::Random(), .win_count = 0});
    m_graph.insert_vertex(Node{.vector = NodeVector::Random(), .win_count = 0});
}

template <unsigned dimension> void GPNet<dimension>::train(array &data) {
    using namespace Eigen;
    using NodeM = Matrix<float, dimension, Dynamic>;

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
        RowVectorXf distances =
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
            node_vec.array() += ((input - node_vec).array() / m_local_opt_coeff /
                                 (win_count + 1));
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
std::pair<float, float>
GPNet<dimension>::threshold(size_t index, const NodeVector &x_vector) {
    using namespace Eigen;

    Matrix<float, dimension, dimension> local_cov;
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
    local_cov += MatrixXf::Identity(dimension, dimension) * m_sigma_2;

    // if (neighbour_cnt > dimension) {
    //     Matrix<bool, dimension, 1> comp;
    //     SelfAdjointEigenSolver<Matrix<float, dimension, dimension>>
    //     eigens(local_cov);
    //     comp = ((eigens.eigenvalues()).array() >=
    //     sigma_2)
    //                .template cast<bool>();

    //     unsigned count = comp.count();
    //     Matrix<float, dimension, Dynamic> p_components(dimension, count);
    //     VectorXf eigenvals(count);

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
    //         + MatrixXf::Identity(dimension, dimension) * sigma_2;
    // }

    SelfAdjointEigenSolver<Matrix<float, dimension, dimension>> eigens(
        local_cov);
    static const float const_coeff = std::sqrt(std::pow(2 * M_PI, dimension));
    float determinant_sqrt = std::sqrt(eigens.eigenvalues().prod());
    float threshold = 1;
    Matrix<float, dimension, dimension> cov_inv =
        eigens.eigenvectors() *
        eigens.eigenvalues().array().inverse().matrix().asDiagonal() *
        eigens.eigenvectors().transpose();
    if (neighbour_cnt) {
        for (auto edge : *iter) {
            NodeVector &node_vec = m_graph[edge.head].value().vector;
            float prob = std::exp(-((x_vector - node_vec).transpose() *
                                    cov_inv * (x_vector - node_vec))(0) /
                                  2) /
                         const_coeff / determinant_sqrt;

            if (prob < threshold)
                threshold = prob;
        }
    } else
        threshold = 0.55;

    float xprob = std::exp(-((x_vector - winner_vec).transpose() * cov_inv *
                             (x_vector - winner_vec))(0) /
                           2) /
                  const_coeff / determinant_sqrt;
    return {
        threshold,
        xprob,
    };
}
template <unsigned dimension> float GPNet<dimension>::predict(array &data) {
    using namespace Eigen;

    Map<NodeVector> input(data.data());

    float prob = 0;
    size_t wins = 0;
    for (auto &node : m_graph) {
        Matrix<float, dimension, dimension> local_cov;
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
        local_cov += MatrixXf::Identity(dimension, dimension) * m_sigma_2;

        SelfAdjointEigenSolver<Matrix<float, dimension, dimension>> eigens(
            local_cov);
        static const float const_coeff =
            std::sqrt(std::pow(2 * M_PI, dimension));
        float determinant_sqrt = std::sqrt(eigens.eigenvalues().prod());
        Matrix<float, dimension, dimension> cov_inv =
            eigens.eigenvectors() *
            eigens.eigenvalues().array().inverse().matrix().asDiagonal() *
            eigens.eigenvectors().transpose();

        NodeVector &node_vec = node.value().vector;
        float xprob = std::exp(-((input - node_vec).transpose() * cov_inv *
                                 (input - node_vec))(0) /
                               2) /
                      const_coeff / determinant_sqrt;
        prob += xprob * node.value().win_count;
        wins += node.value().win_count;
    }

    if (wins)
        prob /= wins;
    return prob;
}
} // namespace LogAnomaly

#endif // LOG_ANOMALY_GPSOINN_HXX

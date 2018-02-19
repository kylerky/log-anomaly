#ifndef LOG_ANOMALY_GRAPH_HXX
#define LOG_ANOMALY_GRAPH_HXX

#include "multiset.hxx"
#include <algorithm>
#include <forward_list>
#include <functional>
#include <unordered_set>
#include <utility>

namespace LogAnomaly {

template <typename ValueT, typename WeightT,
          typename Compare = std::less<ValueT>,
          typename Allocator = std::allocator<ValueT>>
class Digraph {
  public:
    typedef typename multiset<ValueT, Compare>::size_type index_t;

    friend void swap(Digraph &left, Digraph &right) noexcept {
        using std::swap;
        swap(left.vertices, right.vertices);
    }

    // edge node
    struct EdgeNode {
        WeightT weight;
        index_t head;
    };
    typedef typename std::forward_list<EdgeNode>::iterator edge_iterator;
    typedef typename std::forward_list<EdgeNode>::const_iterator
        const_edge_iterator;

    // class Vertex
    class Vertex {
        friend Digraph;

      public:
        edge_iterator begin() noexcept { return edges.begin(); }

        const_edge_iterator begin() const noexcept { return edges.cbegin(); }

        const_edge_iterator cbegin() const noexcept { return edges.cbegin(); }
        edge_iterator end() noexcept { return edges.end(); }

        const_edge_iterator end() const noexcept { return edges.cend(); }

        const_edge_iterator cend() const noexcept { return edges.cend(); }
        edge_iterator before_begin() noexcept { return edges.before_begin(); }

        const_edge_iterator before_begin() const noexcept {
            return edges.cbefore_begin();
        }

        const_edge_iterator cbefore_begin() const noexcept {
            return edges.cbefore_begin();
        }

        ValueT &value() { return val; }
        const ValueT &value() const { return val; }

      private:
        ValueT val;
        std::forward_list<EdgeNode> edges;
    };

    // loads of typedefs
    typedef multiset<Vertex, Compare, Allocator> set_type;

    typedef typename multiset<ValueT, Compare>::size_type size_type;
    typedef ValueT value_type;
    typedef value_type &reference;
    typedef const value_type &const_reference;
    typedef WeightT weight_type;
    typedef Compare value_compare;
    typedef typename set_type::iterator vertex_iterator;
    typedef typename set_type::const_iterator const_vertex_iterator;

    // constructors
    Digraph() : Digraph(Compare()) {}
    explicit Digraph(const Compare &comp, const Allocator &alloc = Allocator())
        : vertices(comp, alloc) {}
    // copy constructors
    Digraph(const Digraph &other) = default;
    Digraph(Digraph &&other) = default;

    /* operators */
    Digraph &operator=(Digraph other) {
        std::swap(*this, other);
        return *this;
    }

    // destructors
    virtual ~Digraph() = default;

    // vertex related modifiers
    vertex_iterator insert_vertex(const ValueT &value) {
        Vertex vertex;
        vertex.val = value;
        return vertices.insert(vertex);
    }
    vertex_iterator insert_vertex(ValueT &&value) {
        Vertex vertex;
        vertex.val = value;
        return vertices.insert(std::move(vertex));
    }
    vertex_iterator erase_vertex(const_vertex_iterator pos);
    vertex_iterator erase_vertex(const_vertex_iterator beg,
                                 const_vertex_iterator end);
    size_type erase_vertex(const ValueT &key);

    // edge related modifiers
    void insert_edge(const_vertex_iterator tail, const_vertex_iterator head,
                     const WeightT &weight = 1); // tail -> head
    void insert_edge(index_t tail, index_t head, const WeightT &weight = 0) {
        insert_edge(get_vertex_iterator(tail), get_vertex_iterator(head),
                    weight);
    }
    edge_iterator erase_after_edge(const_vertex_iterator tail,
                                   const_edge_iterator edge);
    edge_iterator
    erase_after_edge(const_vertex_iterator tail, const_edge_iterator beg,
                     const_edge_iterator end); // range: (beg, end)

    /* iterators */
    // vertex related
    vertex_iterator begin() noexcept { return vertices.begin(); }
    const_vertex_iterator begin() const noexcept { return vertices.cbegin(); }
    const_vertex_iterator cbegin() const noexcept { return vertices.cbegin(); }
    vertex_iterator end() noexcept { return vertices.end(); }
    const_vertex_iterator end() const noexcept { return vertices.cend(); }
    const_vertex_iterator cend() const noexcept { return vertices.cend(); }

    // index related
    index_t get_index(const_vertex_iterator vertex) const noexcept {
        typename set_type::const_iterator &ver_iter = vertex;
        return vertices.get_index(ver_iter);
    }
    vertex_iterator get_vertex_iterator(index_t index) {
        return vertices.get_iterator(index);
    }
    const_vertex_iterator get_vertex_iterator(index_t index) const {
        return vertices.get_iterator(index);
    }
    Vertex &operator[](index_t index) { return vertices[index]; }
    const Vertex &operator[](index_t index) const { return vertices[index]; }
    bool valid(index_t index) const { return vertices.valid(index); }

    //
    void clear() noexcept { vertices.clear(); }

    /* capacity */
    size_type vertex_count() const noexcept { return vertices.size(); }
    [[nodiscard]] bool empty() const noexcept { return vertex_count() == 0; }

  protected:
    set_type vertices;

  private:
    inline typename set_type::iterator
    iter_remove_c(typename set_type::const_iterator citer) {
        index_t tail_pos = citer.index();
        return vertices.get_iterator(tail_pos);
    }
};
template <typename ValueT, typename WeightT,
          typename Compare = std::less<ValueT>,
          typename Allocator = std::allocator<ValueT>>
class UndirectedGraph : protected Digraph<ValueT, WeightT, Compare, Allocator> {
  private:
    using Digraph = Digraph<ValueT, WeightT, Compare, Allocator>;

  public:
    typedef typename Digraph::EdgeNode EdgeNode;
    typedef typename Digraph::Vertex Vertex;
    typedef typename Digraph::vertex_iterator vertex_iterator;
    typedef typename Digraph::const_vertex_iterator const_vertex_iterator;
    typedef typename Digraph::edge_iterator edge_iterator;
    typedef typename Digraph::const_edge_iterator const_edge_iterator;
    typedef ValueT value_type;
    typedef value_type &reference;
    typedef const value_type &const_reference;
    typedef WeightT weight_type;
    typedef Compare value_compare;
    typedef typename Digraph::index_t index_t;
    typedef typename Digraph::size_type size_type;

    // constructors
    UndirectedGraph() : UndirectedGraph(Compare()) {}
    explicit UndirectedGraph(const Compare &comp) : Digraph(comp) {}
    UndirectedGraph(const UndirectedGraph &other) = default;
    UndirectedGraph(UndirectedGraph &&other) = default;

    // destructors
    ~UndirectedGraph() {}

    // swap
    friend void swap(UndirectedGraph &left, UndirectedGraph &right) noexcept {
        using std::swap;
        Digraph &left_base = left;
        Digraph &right_base = right;
        swap(left_base, right_base);
    }

    // operator=
    UndirectedGraph &operator=(UndirectedGraph other) {
        std::swap(*this, other);
        return *this;
    }

    // vertex realted modifiers
    template <typename Type> vertex_iterator insert_vertex(Type &&value) {
        return Digraph::insert_vertex(std::forward<Type>(value));
    }
    vertex_iterator erase_vertex(const_vertex_iterator pos) {
        return Digraph::erase_vertex(pos);
    }
    vertex_iterator erase_vertex(const_vertex_iterator beg,
                                 const_vertex_iterator end) {
        return Digraph::erase_vertex(beg, end);
    }
    size_type erase_vertex(const ValueT &key) {
        return Digraph::erase_vertex(key);
    }

    // edge related modifiers
    void insert_edge(const_vertex_iterator v1, const_vertex_iterator v2,
                     const WeightT &weight = 0) {
        Digraph::insert_edge(v1, v2, weight);
        Digraph::insert_edge(v2, v1, weight);
    }
    void insert_edge(index_t v1, index_t v2, const WeightT &weight = 1) {
        Digraph::insert_edge(v1, v2, weight);
        Digraph::insert_edge(v2, v1, weight);
    }
    edge_iterator erase_after_edge(const_vertex_iterator v1,
                                   const_edge_iterator edge);

    // iterator related
    vertex_iterator begin() noexcept { return Digraph::begin(); }
    const_vertex_iterator begin() const noexcept { return Digraph::cbegin(); }
    const_vertex_iterator cbegin() const noexcept { return Digraph::cbegin(); }
    vertex_iterator end() noexcept { return Digraph::end(); }
    const_vertex_iterator end() const noexcept { return Digraph::cend(); }
    const_vertex_iterator cend() const noexcept { return Digraph::cend(); }

    // index related
    index_t get_index(const_vertex_iterator vertex) const noexcept {
        return Digraph::get_index(vertex);
    }
    vertex_iterator get_vertex_iterator(index_t index) {
        return Digraph::get_vertex_iterator(index);
    }
    const_vertex_iterator get_vertex_iterator(index_t index) const {
        return Digraph::get_vertex_iterator(index);
    }
    typename Digraph::Vertex &operator[](index_t index) {
        return Digraph::operator[](index);
    }
    typename Digraph::Vertex const &operator[](index_t index) const {
        return Digraph::operator[](index);
    }
    bool valid(index_t index) const { return Digraph::valid(index); }

    //
    void clear() noexcept { Digraph::clear(); }

    /* capacity */
    size_type vertex_count() const noexcept { return Digraph::vertex_count(); }
    [[nodiscard]] bool empty() const noexcept { return Digraph::empty(); }
};

} // namespace LogAnomaly

namespace LogAnomaly {

template <typename ValueT, typename WeightT, typename Compare,
          typename Allocator>
typename Digraph<ValueT, WeightT, Compare, Allocator>::vertex_iterator
Digraph<ValueT, WeightT, Compare, Allocator>::erase_vertex(
    const_vertex_iterator pos) {
    typename set_type::const_iterator &iter = pos;
    auto index = iter.index();
    auto result = vertices.erase(iter);
    for (auto &vertex : vertices) {
        auto pre = vertex.edges.cbefore_begin();
        for (auto edge = vertex.edges.cbegin(); edge != vertex.edges.cend();) {
            if (edge->head == index) {
                edge = vertex.edges.erase_after(pre);
            } else {
                ++edge;
                ++pre;
            }
        }
    }
    return result;
}

template <typename ValueT, typename WeightT, typename Compare,
          typename Allocator>
typename Digraph<ValueT, WeightT, Compare, Allocator>::vertex_iterator
Digraph<ValueT, WeightT, Compare, Allocator>::erase_vertex(
    const_vertex_iterator beg, const_vertex_iterator end) {
    typename set_type::const_iterator &beg_iter = beg;
    typename set_type::const_iterator &end_iter = end;
    std::unordered_set<index_t> indices;
    for (auto iter = beg_iter; iter != end_iter; ++iter) {
        indices.insert(iter.index());
    }
    auto result = vertices.erase(beg_iter, end_iter);
    for (auto &vertex : vertices) {
        auto pre = vertex.edges.cbefore_begin();
        for (auto edge = vertex.edges.cbegin(); edge != vertex.edges.cend();) {
            if (indices.count(edge->head) == 1) {
                edge = vertex.edges.erase_after(pre);
            } else {
                ++edge;
                ++pre;
            }
        }
    }
    return result;
}
template <typename ValueT, typename WeightT, typename Compare,
          typename Allocator>
typename Digraph<ValueT, WeightT, Compare, Allocator>::size_type
Digraph<ValueT, WeightT, Compare, Allocator>::erase_vertex(const ValueT &key) {
    std::unordered_set<index_t> indices;
    for (auto iter = vertices.cbegin(); iter != vertices.cend(); ++iter) {
        if (iter->value() == key)
            indices.insert(iter.index());
    }

    for (auto index : indices)
        vertices.erase(vertices.get_iterator(index));

    for (auto &vertex : vertices) {
        auto pre = vertex.edges.cbefore_begin();
        for (auto edge = vertex.edges.cbegin(); edge != vertex.edges.cend();) {
            if (indices.count(edge->head) == 1) {
                edge = vertex.edges.erase_after(pre);
            } else {
                ++edge;
                ++pre;
            }
        }
    }
    return indices.size();
}
template <typename ValueT, typename WeightT, typename Compare,
          typename Allocator>
void Digraph<ValueT, WeightT, Compare, Allocator>::insert_edge(
    const_vertex_iterator tail, const_vertex_iterator head,
    const WeightT &weight) {
    typename set_type::const_iterator &ctail_iter = tail;
    typename set_type::const_iterator &head_iter = head;
    auto tail_iter = iter_remove_c(ctail_iter);

    tail_iter->edges.push_front({.weight = weight, .head = head_iter.index()});
}

template <typename ValueT, typename WeightT, typename Compare,
          typename Allocator>
typename Digraph<ValueT, WeightT, Compare, Allocator>::edge_iterator
Digraph<ValueT, WeightT, Compare, Allocator>::erase_after_edge(
    const_vertex_iterator tail, const_edge_iterator edge) {
    typename set_type::const_iterator &ctail_iter = tail;
    auto tail_iter = iter_remove_c(ctail_iter);

    auto result = tail_iter->edges.erase_after(edge);
    return result;
}
template <typename ValueT, typename WeightT, typename Compare,
          typename Allocator>
typename Digraph<ValueT, WeightT, Compare, Allocator>::edge_iterator
Digraph<ValueT, WeightT, Compare, Allocator>::erase_after_edge(
    const_vertex_iterator tail, const_edge_iterator beg,
    const_edge_iterator end) {
    typename set_type::const_iterator &ctail_iter = tail;
    auto tail_iter = iter_remove_c(ctail_iter);

    auto before_count = vertices.size();
    auto result = tail_iter->edges.erase_after(beg, end);
    return result;
}

template <typename ValueT, typename WeightT, typename Compare,
          typename Allocator>
typename UndirectedGraph<ValueT, WeightT, Compare, Allocator>::edge_iterator
UndirectedGraph<ValueT, WeightT, Compare, Allocator>::erase_after_edge(
    const_vertex_iterator v1, const_edge_iterator edge) {
    auto target = edge;
    ++target;

    auto result = Digraph::erase_after_edge(v1, edge);
    auto v2 = get_vertex_iterator(target->head);
    if (v1 != v2) {
        index_t index = v1.index();
        for (auto iter = v2->cbegin(), pre = v2->cbefore_begin();
             iter != v2->cend(); ++iter, ++pre) {
            if (iter->head == index && iter->weight == target->weight) {
                Digraph::erase_after_edge(v2, pre);
                break;
            }
        }
    }
    return result;
}

} // namespace LogAnomaly

#endif // LOG_ANOMALY_GRAPH_HXX

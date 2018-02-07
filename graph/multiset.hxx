#ifndef LOG_ANOMALY_MULTISET_HXX
#define LOG_ANOMALY_MULTISET_HXX

#include <cstring>
#include <functional>
#include <iterator>
#include <vector>

namespace LogAnomaly {

template <size_t... args> struct static_max;

template <size_t arg> struct static_max<arg> {
    static const size_t value = arg;
};

template <size_t arg1, size_t arg2, size_t... others>
struct static_max<arg1, arg2, others...> {
    static const size_t value = arg1 >= arg2
                                    ? static_max<arg1, others...>::value
                                    : static_max<arg2, others...>::value;
};

template <typename T, typename T_2> struct multiset_variant {
  private:
    using data_t = typename std::aligned_storage<
        static_max<sizeof(T), sizeof(T_2)>::value,
        static_max<alignof(T), alignof(T_2)>::value>::type;
    bool bit;
    data_t data;

    inline void copy(const data_t *other) {
        if (bit)
            new (&data) T_2(*reinterpret_cast<const T_2 *>(other));
        else
            new (&data) T(*reinterpret_cast<const T *>(other));
    }
    inline void move(const data_t *other) {
        if (bit)
            new (&data) T_2(std::move(*reinterpret_cast<const T_2 *>(other)));
        else
            new (&data) T(std::move(*reinterpret_cast<const T *>(other)));
    }
    inline void destroy() {
        if (bit)
            reinterpret_cast<T_2 *>(&data)->~T_2();
        else
            reinterpret_cast<T *>(&data)->~T();
    }

  public:
    multiset_variant() : bit(false) { new (&data) T(); }
    // copy constructor
    multiset_variant(const multiset_variant &other) : bit(other.bit) {
        copy(&other.data);
    }
    // move constructor
    multiset_variant(multiset_variant &&other) : bit(other.bit) {
        move(&other.data);
    }

    ~multiset_variant() { destroy(); }

    multiset_variant &operator=(multiset_variant other) {
        bit = other.bit;
        copy(&other.data);
        return *this;
    }

    bool second() const { return bit; }

    template <typename Type, typename... Args> void set(Args &&... args) {
        using std::is_same;

        destroy();

        new (&data) Type(std::forward<Args>(args)...);
        bit = is_same<Type, T_2>::value;
    }

    template <typename Type> Type &get() {
        using std::is_same;
        if (bit == is_same<Type, T_2>::value) {
            return *reinterpret_cast<Type *>(&data);
        } else
            throw std::bad_cast();
    }
    template <typename Type> const Type &get() const {
        using std::is_same;
        if (bit == is_same<Type, T_2>::value) {
            return *reinterpret_cast<const Type *>(&data);
        } else
            throw std::bad_cast();
    }
};

template <typename Key, typename Compare = std::less<Key>,
          typename Allocator = std::allocator<multiset_variant<size_t, Key>>>
class multiset
    : protected std::vector<multiset_variant<size_t, Key>, Allocator> {
  public:
    class iterator;
    class const_iterator;
    class reverse_iterator;
    class const_reverse_iterator;

    // member types
    typedef Key key_type;
    typedef Key value_type;
    typedef Compare key_compare;
    typedef Compare value_compare;
    typedef std::vector<multiset_variant<size_t, Key>, Allocator> vector;

    typedef typename vector::allocator_type allocator_type;
    typedef typename vector::size_type size_type;
    typedef std::ptrdiff_t difference_type;
    typedef value_type &reference;
    typedef const value_type &const_reference;

    /* constructors */
    multiset() : multiset(Compare()){};
    explicit multiset(const Compare &comp, const Allocator &alloc = Allocator())
        : less(comp), vector(alloc){};
    // copy constructors
    multiset(const multiset &other) = default;
    multiset(const multiset &other, const Allocator &alloc);
    // move constructors
    multiset(multiset &&other) noexcept = default;
    multiset(multiset &&other, const Allocator &alloc);

    friend void swap(multiset &first, multiset &second) noexcept {
        using std::swap;
        swap(dynamic_cast<vector &>(first), dynamic_cast<vector &>(second));
        swap(first.less, second.less);
        swap(first.count, second.count);
        swap(first.first_free, second.first_free);
    }
    /* destructors */
    virtual ~multiset() = default;

    /* modifiers */
    // insert
    iterator insert(const value_type &value);
    iterator insert(value_type &&value);
    // erase
    iterator erase(const_iterator pos);
    iterator erase(const_iterator beg, const_iterator end);
    size_type erase(const key_type &key);
    // clear
    void clear() noexcept {
        vector::clear();
        first_free = 0;
        count = 0;
    };

    /* capacity */
    [[nodiscard]] bool empty() const noexcept { return count == 0; }
    size_type size() const noexcept { return count; };

    /* operators */
    multiset &operator=(multiset other);
    value_type &operator[](size_type index);
    const value_type &operator[](size_type index) const;

    allocator_type get_allocator() const { return vector::get_allocator(); };

    // iterators
    iterator begin() noexcept;
    const_iterator begin() const noexcept;
    const_iterator cbegin() const noexcept;
    iterator end() noexcept;
    const_iterator end() const noexcept;
    const_iterator cend() const noexcept;

    size_type get_index(const_iterator iter) const noexcept;
    iterator get_iterator(size_type index) noexcept;
    const_iterator get_iterator(size_type index) const noexcept;

  private:
    Compare less;
    size_type count = 0;
    size_type first_free = 0;
};

} // namespace LogAnomaly

namespace LogAnomaly {

template <typename Key, typename Compare, typename Allocator>
class multiset<Key, Compare, Allocator>::iterator {
    friend multiset<Key, Compare, Allocator>;
    friend multiset<Key, Compare, Allocator>::const_iterator;

  public:
    friend bool operator==(const iterator &lhs, const iterator &rhs) {
        return lhs.iter == rhs.iter;
    }
    friend bool operator!=(const iterator &lhs, const iterator &rhs) {
        return !(lhs == rhs);
    }

    typedef std::ptrdiff_t difference_type;
    typedef Key value_type;
    typedef value_type *pointer;
    typedef value_type &reference;
    typedef std::bidirectional_iterator_tag iterator_category;

    iterator() : vector(nullptr) {}
    virtual ~iterator() = default;

    // dereference
    reference operator*() const { return iter->template get<Key>(); }
    pointer operator->() const { return &(operator*()); }
    // increment
    iterator &operator++() {
        ++iter;
        while (iter != vector->end() && !iter->second()) {
            ++iter;
        }
        return *this;
    }
    iterator operator++(int) {
        iterator previous = *this;
        this->operator++();
        return previous;
    }
    // decrement
    iterator &operator--() {
        while (iter != vector->begin() && !iter->second()) {
            --iter;
        }
        return *this;
    }
    iterator operator--(int) {
        iterator previous = *this;
        this->operator--();
        return previous;
    }

  private:
    typename vector::iterator iter;
    const vector *vector;
};

template <typename Key, typename Compare, typename Allocator>
class multiset<Key, Compare, Allocator>::const_iterator {
    friend multiset<Key, Compare, Allocator>;

  public:
    friend bool operator==(const const_iterator &lhs,
                           const const_iterator &rhs) {
        return lhs.iter == rhs.iter;
    }
    friend bool operator!=(const const_iterator &lhs,
                           const const_iterator &rhs) {
        return !(lhs == rhs);
    }
    typedef std::ptrdiff_t difference_type;
    typedef const Key value_type;
    typedef value_type *pointer;
    typedef value_type &reference;
    typedef std::bidirectional_iterator_tag iterator_category;

    const_iterator() : vector(nullptr) {}
    const_iterator(const iterator &other)
        : iter(other.iter), vector(other.vector) {}

    virtual ~const_iterator() = default;

    // dereference
    reference operator*() const { return iter->template get<Key>(); }
    pointer operator->() const { return &(operator*()); }
    // increment
    const_iterator &operator++() {
        ++iter;
        while (iter != vector->cend() && !iter->second()) {
            ++iter;
        }
        return *this;
    }
    const_iterator operator++(int) {
        iterator previous = *this;
        this->operator++();
        return previous;
    }
    // decrement
    const_iterator &operator--() {
        while (iter != vector->cbegin() && !iter->second()) {
            --iter;
        }
        return *this;
    }
    const_iterator operator--(int) {
        iterator previous = *this;
        this->operator--();
        return previous;
    }

  private:
    typename vector::const_iterator iter;
    const vector *vector;
};

template <typename T, typename Compare, typename Allocator>
typename multiset<T, Compare, Allocator>::iterator
multiset<T, Compare, Allocator>::insert(const T &value) {
    if (count == vector::size()) {
        // construct the union to be inserted
        multiset_variant<size_t, T> elem;
        elem.template set<T>(value);
        // call the push back
        vector::push_back(elem);
        // increment count
        ++count;

        iterator result;
        result.iter = vector::begin() + (vector::size() - 1);
        result.vector = this;
        return result;
    }

    // get the elem first
    auto &elem = vector::operator[](first_free);
    // store the position
    size_type pos = first_free;

    // pop front operation on the free list
    first_free = elem.template get<size_t>();
    // assign
    elem.template set<T>(value);

    // increment count
    ++count;

    iterator result;
    result.iter = vector::begin() + pos;
    result.vector = this;
    return result;
}

template <typename T, typename Compare, typename Allocator>
typename multiset<T, Compare, Allocator>::iterator
multiset<T, Compare, Allocator>::insert(T &&value) {
    if (count == vector::size()) {
        // construct the union to be inserted
        multiset_variant<size_t, T> elem;
        elem.template set<T>(value);
        // call the push back
        vector::push_back(std::move(elem));
        // increment count
        ++count;

        iterator result;
        result.iter = vector::begin() + (vector::size() - 1);
        result.vector = this;
        return result;
    }

    // get the elem first
    auto &elem = vector::operator[](first_free);
    // store the position
    size_type pos = first_free;

    // pop front operation on the free list
    first_free = elem.template get<size_t>();
    // assign
    elem.template set<T>(value);

    // increment count
    ++count;

    iterator result;
    result.iter = vector::begin() + pos;
    result.vector = this;
    return result;
}

template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::iterator
multiset<Key, Compare, Allocator>::erase(const_iterator pos) {
    typename vector::difference_type index = pos.iter - vector::cbegin();
    auto iter = vector::begin() + index;

    iterator result;
    result.iter = iter;
    result.vector = pos.vector;
    ++result;

    iter->template set<size_t>(first_free);
    first_free = index;

    // decrement count
    --count;

    return result;
}

template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::iterator
multiset<Key, Compare, Allocator>::erase(const_iterator beg,
                                         const_iterator end) {
    while (beg != end)
        beg = erase(beg);

    iterator result;
    result.iter = vector::begin + (end.iter - vector::cbegin());
    result.vector = this;

    return result;
}

template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::size_type
multiset<Key, Compare, Allocator>::erase(const key_type &key) {
    auto iter = cbegin();
    size_type count = 0;
    while (iter != cend()) {
        if (!less(*iter, key) && !less(key, *iter)) {
            iter = erase(iter);
            ++count;
        } else
            ++iter;
    }
    return count;
}

template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::iterator
multiset<Key, Compare, Allocator>::begin() noexcept {
    iterator iter;
    iter.iter = vector::begin();
    while (iter.iter != vector::end() && !iter.iter->second()) {
        ++iter.iter;
    }
    iter.vector = this;
    return iter;
}

template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::const_iterator
multiset<Key, Compare, Allocator>::begin() const noexcept {
    iterator iter;
    iter.iter = vector::cbegin();
    while (iter.iter != vector::cend() && !iter->second()) {
        ++iter.iter;
    }
    iter.vector = this;
    return iter;
}

template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::const_iterator
multiset<Key, Compare, Allocator>::cbegin() const noexcept {
    const_iterator iter;
    iter.iter = vector::cbegin();
    while (iter.iter != vector::cend() && !iter.iter->second()) {
        ++iter.iter;
    }
    iter.vector = this;
    return iter;
}

template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::iterator
multiset<Key, Compare, Allocator>::end() noexcept {
    iterator iter;
    iter.iter = vector::end();
    iter.vector = this;
    return iter;
}

template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::const_iterator
multiset<Key, Compare, Allocator>::end() const noexcept {
    iterator iter;
    iter.iter = vector::cend();
    iter.vector = this;
    return iter;
}

template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::const_iterator
multiset<Key, Compare, Allocator>::cend() const noexcept {
    const_iterator iter;
    iter.iter = vector::cend();
    iter.vector = this;
    return iter;
}

template <typename Key, typename Compare, typename Allocator>
multiset<Key, Compare, Allocator>::multiset(const multiset &other,
                                            const Allocator &alloc)
    : vector(other, alloc), less(other.less), count(other.count),
      first_free(other.first_free) {}

template <typename Key, typename Compare, typename Allocator>
multiset<Key, Compare, Allocator>::multiset(multiset &&other,
                                            const Allocator &alloc)
    : vector(other, alloc), less(other.less), count(other.count),
      first_free(other.first_free) {}

template <typename Key, typename Compare, typename Allocator>
multiset<Key, Compare, Allocator> &multiset<Key, Compare, Allocator>::
operator=(multiset other) {
    std::swap(*this, other);
    return *this;
}
template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::size_type
multiset<Key, Compare, Allocator>::get_index(const_iterator iter) const
    noexcept {
    return iter.iter - vector::cbegin();
}
template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::iterator
multiset<Key, Compare, Allocator>::get_iterator(size_type index) noexcept {
    iterator iter;
    iter.iter = vector::begin() + index;
    iter.vector = this;
    return iter;
}
template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::const_iterator
multiset<Key, Compare, Allocator>::get_iterator(size_type index) const
    noexcept {
    const_iterator iter;
    iter.iter = vector::cbegin() + index;
    iter.vector = this;
    return iter;
}
template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::value_type &
    multiset<Key, Compare, Allocator>::operator[](size_type index) {
    return vector::operator[](index);
}
template <typename Key, typename Compare, typename Allocator>
typename multiset<Key, Compare, Allocator>::value_type const &
    multiset<Key, Compare, Allocator>::operator[](size_type index) const {
    return vector::operator[](index);
}

} // namespace LogAnomaly

#endif // LOG_ANOMALY_MULTISET_HXX

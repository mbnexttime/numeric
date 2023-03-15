#ifndef LINALG_LINALG_HH
#define LINALG_LINALG_HH

#include <ostream>
#include <vector>

using bool_t = uint16_t;

class Expression
{
};

template <class T, size_t N>
class Array;

template <class E>
Array<typename E::value_type, E::dimensions> evaluate(const E &expr)
{
    using value_type = typename E::value_type;
    const auto n = expr.size();
    const auto dimensions = expr.shape();
    Array<value_type, E::dimensions> result(dimensions);
    for (int i = 0; i < n; ++i)
    {
        result(i) = expr.evaluate(i);
    }
    return result;
}

template <class T, size_t N>
class Array : public Expression
{

public:
    using value_type = T;
    static const size_t dimensions = N;

private:
    std::vector<T> _data;
    std::vector<size_t> _shape;
    std::vector<size_t> _order;

    size_t get_size() const {
        size_t p = 1;
        for (int i = 0; i < dimensions; i++) {
            p *= _shape[i];
        }
        return p;
    }

    size_t get_index_by_coord(std::vector<int>& args) const {
        size_t real_index = 0;
        for (int i = 0; i < dimensions; i++) {
            real_index = real_index * _shape[_order[i]] + args[_order[i]];
        }
        return real_index;
    }

    size_t recalc_index(int i) const {
        std::vector<int> old_index;
        for (int j = dimensions - 1; j >= 0; j--) {
            old_index.push_back(i % _shape[j]);
            i /= _shape[j];
        }
        reverse(old_index.begin(), old_index.end());
        return get_index_by_coord(old_index);
    }

    int levi_chivita_symbol(std::vector<int> & indexes) {
        int ret = 1;
        for (int i = 0; i < indexes.size(); i++) {
            for (int j = i + 1; j < indexes.size(); j++) {
                int a = indexes[j] - indexes[i];
                ret *= ((a > 0) - (a < 0));
            }
        }
        return ret;
    }

    int bpow(int a, int t) {
        if (t == 0) {
            return 1;
        }
        if (t % 2 == 0) {
            return bpow(a, t / 2) * bpow(a, t / 2);
        } else {
            return a * bpow(a, t - 1);
        }
    }

public:
    Array(std::initializer_list<T> rhs) : _data(rhs), _shape({ rhs.size() }) {
        _order = { 0 };
    }
    Array(std::initializer_list<T> rhs, std::initializer_list<size_t> shape)
      : _data(rhs), _shape(shape) {
        for (int i = 0; i < dimensions; i++) {
            _order.push_back(i);
        }
    }
    Array(std::vector<size_t> shape): _shape(std::move(shape)) {
        _data.resize(get_size());
        for (int i = 0; i < dimensions; i++) {
            _order.push_back(i);
        }
    }
    template <class E>
    Array(const E &expr,
          typename std::enable_if<std::is_base_of<Expression, E>::value, E>::type *dummy = nullptr) : Array(::evaluate(expr)) {}

    Array() = default;
    ~Array() = default;
    Array(Array &&rhs) = default;
    Array(const Array &rhs) = default;
    Array &operator=(Array &&rhs) = default;
    Array &operator=(const Array &rhs) = default;

    T &operator()(int i) { return this->_data[recalc_index(i)]; }
    const T &operator()(int i) const { return this->_data[recalc_index(i)]; }
    T evaluate(int i) { return this->_data[recalc_index(i)]; }
    T evaluate(int i) const { return this->_data[recalc_index(i)]; }
    int size() const { return this->_data.size(); }
    std::vector<size_t> shape() const {
        return _shape;
    }

    template <
      typename Arg, typename... Ts,
      typename std::enable_if<std::is_integral<Arg>::value>::type * = nullptr>
    T &operator()(Arg tt, Ts... indexes) {
        std::vector<int> coordinates = {tt, indexes...};
        return this->_data[get_index_by_coord(coordinates)];
    }

    void transpose(std::vector<size_t> new_order) {
        _order = new_order;
    }

    void display(std::ostream &out) const
    {
        out << "Array(";
        const auto n = size();
        for (int i = 0; i < n; i++) {
            out << this->_data[recalc_index(i)] << (i == n - 1 ? "" : ",");
        }
        out << ')';
    }

    int determinant() {
        T ret = 0;
        int rows = _shape[0];
        for (int index = 0; index < bpow(rows, rows); index++) {
            int tempIndex = index;
            std::vector<int> coordinates(rows);
            for (int i = 0; i < rows; i++) {
                coordinates[i] = tempIndex % rows;
                tempIndex /= rows;
            }
            T currentRet = _data[coordinates[0]];
            for (int i = 1; i < rows; i++) {
                currentRet *= _data[i * _shape[1] + coordinates[i]];
            }
            ret += levi_chivita_symbol(coordinates) * currentRet;
        }
        return ret;
    }
};

template <class E>
typename std::enable_if<std::is_base_of<Expression, E>::value, std::ostream &>::type
operator<<(std::ostream &out, const E &expr)
{
    expr.display(out);
    return out;
}

template <class E1, class E2>
class Plus : public Expression
{

public:
    using value_type =
        typename std::common_type<typename E1::value_type, typename E2::value_type>::type;
    static const size_t dimensions = E1::dimensions;

private:
    const E1 &_a;
    const E2 &_b;

public:
    explicit Plus(const E1 &a, const E2 &b) : _a(a), _b(b) {}
    value_type evaluate(int i) { return this->_a.evaluate(i) + this->_b.evaluate(i); }
    value_type evaluate(int i) const { return this->_a.evaluate(i) + this->_b.evaluate(i); }
    int size() const { return this->_a.size(); }
    std::vector<size_t> shape() const {
        return _a.shape();
    }
    void display(std::ostream &out) const
    {
        out << "Plus(" << this->_a << ", " << this->_b << ')';
    }
};

template <class E1, class E2>
class Less : public Expression
{

public:
    using value_type = bool_t;
    static const size_t dimensions = E1::dimensions;

private:
    const E1 &_a;
    const E2 &_b;

public:
    explicit Less(const E1 &a, const E2 &b) : _a(a), _b(b) {}
    value_type evaluate(int i) { return this->_a.evaluate(i) < this->_b.evaluate(i); }
    value_type evaluate(int i) const { return this->_a.evaluate(i) < this->_b.evaluate(i); }
    int size() const { return this->_a.size(); }
    std::vector<size_t> shape() const {
        return _a.shape();
    }
    void display(std::ostream &out) const
    {
        out << "Less(" << this->_a << ", " << this->_b << ')';
    }
};

template <class E1>
class All : public Expression
{

public:
    using value_type = bool_t;
    static const size_t dimensions = 1;

private:
    const E1 &_a;

public:
    explicit All(const E1 &a) : _a(a) {}
    value_type evaluate(int i) { 
        for (int i = 0; i < _a.size(); i++) {
            if (!_a.evaluate(i)) {
                return false;
            }
        }
        return true;
    }
    value_type evaluate(int i) const { 
        for (int i = 0; i < _a.size(); i++) {
            if (!_a.evaluate(i)) {
                return false;
            }
        }
        return true;
    }
    int size() const { return 1; }
    std::vector<size_t> shape() const {
        return { 1 };
    }
    void display(std::ostream &out) const
    {
        out << "All(" << this->_a << ')';
    }
};

template <class E1>
class Any : public Expression
{

public:
    using value_type = bool_t;
    static const size_t dimensions = 1;

private:
    const E1 &_a;

public:
    explicit Any(const E1 &a) : _a(a) {}
    value_type evaluate(int i) { 
        for (int i = 0; i < _a.size(); i++) {
            if (_a.evaluate(i)) {
                return true;
            }
        }
        return false;
    }
    value_type evaluate(int i) const { 
        for (int i = 0; i < _a.size(); i++) {
            if (_a.evaluate(i)) {
                return true;
            }
        }
        return false;
    }
    int size() const { return 1; }
    std::vector<size_t> shape() const {
        return { 1 };
    }
    void display(std::ostream &out) const
    {
        out << "Any(" << this->_a << ')';
    }
};

template <class E1, class E2, class E3>
class Where : public Expression
{

public:
    using value_type = typename std::common_type<typename E2::value_type, typename E3::value_type>::type;;
    static const size_t dimensions = E2::dimensions;

private:
    const E1 &_a;
    const E2 &_b; 
    const E3 &_c;

public:
    explicit Where(const E1 &a, const E2 &b, const E3 &c) : _a(a), _b(b), _c(c) {}
    value_type evaluate(int i) { 
        if (_a.evaluate(i)) {
            return _b.evaluate(i);
        } else {
            return _c.evaluate(i);
        }
    }
    value_type evaluate(int i) const { 
        if (_a.evaluate(i)) {
            return _b.evaluate(i);
        } else {
            return _c.evaluate(i);
        }
    }
    int size() const { return _a.size(); }
    std::vector<size_t> shape() const {
        return _b.shape();
    }
    void display(std::ostream &out) const
    {
        out << "Where(" << this->_a << ',' << this->_b << ',' << this->_c << ')';
    }
};



template <class E1, class E2>
typename std::enable_if<std::is_base_of<Expression, E1>::value &&
                            std::is_base_of<Expression, E2>::value,
                        Plus<E1, E2>>::type
operator+(const E1 &a, const E2 &b)
{
    return Plus<E1, E2>(a, b);
}

template <class E1, class E2>
typename std::enable_if<std::is_base_of<Expression, E1>::value &&
                            std::is_base_of<Expression, E2>::value,
                        Less<E1, E2>>::type
operator<(const E1 &a, const E2 &b)
{
    return Less<E1, E2>(a, b);
}

template <class E1>
static typename std::enable_if<std::is_base_of<Expression, E1>::value, All<E1>>::type all(const E1 &a) {
    return All<E1>(a);
}

template <class E1>
static typename std::enable_if<std::is_base_of<Expression, E1>::value, Any<E1>>::type any(const E1 &a) {
    return Any<E1>(a);
}

template <class E1, class E2, class E3>
static typename std::enable_if<
    std::is_base_of<Expression, E1>::value && 
    std::is_base_of<Expression, E2>::value && 
    std::is_base_of<Expression, E3>::value, Where<E1, E2, E3>>::type where(const E1 &a, const E2 &b, const E3 &c) {
    return Where<E1, E2, E3>(a, b, c);
}

#endif // vim:filetype=cpp

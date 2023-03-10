#include <iostream>

#include <linalg/linalg.hh>

int main(int argc, char* argv[]) {
    using T = float;
    
    std::vector<size_t> shape = {2, 3};
    Array<T, 2> a(shape), b(shape);
    a(0, 2) = 22.0f;
    std::cout << a << std::endl;
    a.transpose({1, 0});
    std::cout << a(2, 0) << std::endl;
    a.transpose({0, 1});
    std::cout << a(0, 2) << std::endl;
    std::cout << a << std::endl;
    std::vector<size_t> sq_shape = {2, 2};
    Array<T, 2> c(sq_shape);
    c(0, 0) = 1;
    c(1, 0) = 2;
    c(0, 1) = 2;
    c(1, 1) = 1;
    std::cout << c.determinant() << std::endl;
    c.transpose({1, 0});
    std::cout << c.determinant() << std::endl;
    std::cout << c << std::endl;
    return 0;
}

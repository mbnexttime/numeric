#include <iostream>

#include <linalg/linalg.hh>

int main(int argc, char* argv[]) {
    using T = float;
    
    std::vector<size_t> shape = {2, 3};
    Array<T, 2> a(shape), b(shape);
    std::cout << evaluate(a < a) << std::endl;
    return 0;
}

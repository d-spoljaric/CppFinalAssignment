#include <cmath>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <map>
#include <stdexcept>
#include <utility>

template <typename T>
class Vector{
    // Your implementation of the Vector class starts here
};

template<typename T, typename U>
typename std::common_type<T,U>::type dot(const Vector<T>& lhs, const Vector<U>& rhs){
    // Your implementation of the dot function starts here
}

template <typename T>
class Matrix{
    // Start your implementation of the matrix class here
};

template<typename T>
int cg(const Matrix<T>& A,
       const Vector<T>& b,
       Vector<T>&       x,
       T                tol     = (T)1e-8,
       int              maxiter = 100){
    // Your implementation of the cg function starts here
}

template <int n, typename T>
class Heat{
    // Your implementation of the heat class starts here
};

int main(int argc, char* argv[]){
    // Your testing of the Heat class starts here
    return 0;
}
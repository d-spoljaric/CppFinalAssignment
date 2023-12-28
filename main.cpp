#include <cmath>
#include <iostream>
#include <initializer_list>
#include <memory>
#include <map>
#include <stdexcept>
#include <utility>

template <typename T>
class Vector {
private:
    int length;
    std::unique_ptr<T[]> data;

public:
    // Constructors
    Vector() : length(0) {}
    Vector(int len) : length(len), data(new T[len]()) {}
    Vector(const Vector& other) : length(other.length), data(new T[other.length]) {
        std::copy(other.data.get(), other.data.get() + length, data.get());
    }
    Vector(Vector&& other) noexcept : length(other.length), data(std::move(other.data)) {
        other.length = 0;
    }
    Vector(std::initializer_list<T> list) : Vector((int)list.size()) {
        std::copy(list.begin(), list.end(), data.get());
    }

    // Assignment Operators
    Vector& operator=(const Vector& other) {
        if (this != &other) {
            length = other.length;
            data.reset(new T[length]); // Adjust memory allocation
            std::copy(other.data.get(), other.data.get() + length, data.get());
        }
        return *this;
    }
    Vector& operator=(Vector&& other) noexcept {
        length = other.length;
        data = std::move(other.data);
        other.length = 0;
        return *this;
    }


    Vector operator-(const Vector& rhs) const {
    Vector result(length);
    for (int i = 0; i < length; ++i) {
        result[i] = data[i] - rhs[i];
    }
    return result;
}

Vector operator*(T scalar) const {
    Vector result(length);
    for (int i = 0; i < length; ++i) {
        result[i] = data[i] * scalar;
    }
    return result;
}
// Inside the Vector class
Vector operator+(const Vector& rhs) const {
    if (length != rhs.length) {
        throw std::invalid_argument("Vectors must be of the same length to add.");
    }
    Vector result(length);
    for (int i = 0; i < length; ++i) {
        result[i] = data[i] + rhs[i];
    }
    return result;
}


   // Access Operators with bounds checking
    T& operator[](int i) {
        if (i < 0 || i >= length) {
            throw std::out_of_range("Vector index out of bounds");
        }
        return data[i];
    }

    const T& operator[](int i) const {
        if (i < 0 || i >= length) {
            throw std::out_of_range("Vector index out of bounds");
        }
        return data[i];
    }

    int len() const { return length; }
};
template<typename T>
Vector<T> operator*(const Vector<T>& vec, const T& scalar) {
    Vector<T> result(vec.len());
    for (int i = 0; i < vec.len(); ++i) {
        result[i] = vec[i] * scalar;
    }
    return result;
}

template<typename T>
Vector<T> operator*(const T& scalar, const Vector<T>& vec) {
    return vec * scalar; // Utilize the above operator
}



// Implement dot function outside of class
template<typename T, typename U>
typename std::common_type<T, U>::type 
dot(const Vector<T>& lhs, const Vector<U>& rhs) {
    if (lhs.len() != rhs.len()) throw std::invalid_argument("Vectors must be of the same length");

    typename std::common_type<T, U>::type result = 0;
    for (int i = 0; i < lhs.len(); ++i) {
        result += lhs[i] * rhs[i];
    }
    return result;
}

template <typename T>
class Matrix {
private:
    std::map<std::pair<int, int>, T> data;
    int rows, cols;

public:
    Matrix(int rows, int cols) : rows(rows), cols(cols) {}

    T& operator[](const std::pair<int, int>& ij) {
        return data[ij];
    }
T operator()(const std::pair<int, int>& ij) const {
    auto it = data.find(ij);
    if (it == data.end()) throw std::runtime_error("Matrix entry not found.");
    return it->second;
}

template<typename U>
Vector<typename std::common_type<T, U>::type> operator*(const Vector<U>& rhs) const {
    if (rhs.len() != cols) {
        throw std::invalid_argument("Vector length must match the number of matrix columns");
    }

    Vector<typename std::common_type<T, U>::type> result(rhs.len());
    for (const auto& kv : data) {
        int i = kv.first.first;
        int j = kv.first.second;
        if (j >= 0 && j < rhs.len()) {
            result[i] += kv.second * rhs[j];
        }
    }
    return result;
}

};

template<typename T>
int cg(const Matrix<T>& A, 
       const Vector<T>& b, 
       Vector<T>& x, 
       T tol = (T)1e-8, 
       int maxiter = 100) {
    Vector<T> r = b - A * x;
    Vector<T> p = r;
    T rsold = dot(r, r);

    for (int i = 0; i < maxiter; ++i) {
        Vector<T> Ap = A * p;
        T alpha = rsold / dot(p, Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        T rsnew = dot(r, r);
        if (sqrt(rsnew) < tol) return i;
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    }
    return -1;
}
template <int n, typename T>
class Heat {
private:
    T alpha;
    int m;
    T dt;
    int total_nodes;
    Matrix<T> M;

public:
    Heat(T alpha, int m, T dt) 
        : alpha(alpha), m(m), dt(dt), total_nodes(std::pow(m, n)), M(total_nodes, total_nodes) {
        T dx = 1.0 / (m + 1);

        // Fill the matrix M
        for (int i = 0; i < total_nodes; ++i) {
            for (int k = 0; k < n; ++k) {
                int stride = std::pow(m, k);
                if (i % stride < m - 1) M[{i, i + stride}] = -alpha * dt / (dx * dx);
                if (i % stride > 0) M[{i, i - stride}] = -alpha * dt / (dx * dx);
                M[{i, i}] += 1;
            }
        }
    }

    Vector<T> exact(T t) const {
        Vector<T> u(total_nodes);
        T factor = std::exp(-n * M_PI * M_PI * alpha * t);

        for (int i = 0; i < total_nodes; ++i) {
            T prod = 1;
            for (int k = 0; k < n; ++k) {
                int index = (i / static_cast<int>(std::pow(m, k))) % m;
                T xk = (index + 1) * (1.0 / (m + 1));
                prod *= std::sin(M_PI * xk);
            }
            u[i] = factor * prod;
        }

        return u;
    }

    Vector<T> solve(T t_final) const {
        Vector<T> u = exact(0); // Initial condition
        Vector<T> u_new(total_nodes);

        int steps = t_final / dt;
        for (int step = 0; step < steps; ++step) {
            cg(M, u, u_new);
            u = u_new;
        }

        return u;
    }
};

int main(int argc, char* argv[]) {
    const double alpha = 0.3125;
    const int m = 3;
    const double dt = 0.1;
    const double t_final = 1.0;

    // Creating a solver for the one-dimensional problem
    Heat<1, double> solver1D(alpha, m, dt);
    auto solution1D = solver1D.solve(t_final);
    auto exactSolution1D = solver1D.exact(t_final);

    // Print results (can be replaced with actual comparison)
    for (int i = 0; i < solution1D.len(); ++i) {
        std::cout << "Numerical: " << solution1D[i] << ", Exact: " << exactSolution1D[i] << std::endl;
    }

    // // Repeat similar steps for the two-dimensional problem
    // Heat<2, double> solver2D(alpha, m, dt);
    // auto solution2D = solver2D.solve(t_final);
    // auto exactSolution2D = solver2D.exact(t_final);

    // // Print results (can be replaced with actual comparison)
    // for (int i = 0; i < solution2D.len(); ++i) {
    //     std::cout << "Numerical: " << solution2D[i] << ", Exact: " << exactSolution2D[i] << std::endl;
    // }   

    return 0;
}

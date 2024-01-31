#include <cmath>
#include <iostream>
#include <initializer_list>
#include <stdexcept>
#include <utility>
#include <map>
#include <chrono>

template <typename T>
class Vector {
public:
    T* elements;
    int length;

    // Constructors

    // Default constructor
    Vector() : elements(nullptr), length(0) {}

    // Size constructor
    explicit Vector(int size) : elements(new T[size]), length(size) {
        std::fill_n(elements, size, T());
    }

    // List constructor
    Vector(std::initializer_list<T> list) : elements(new T[list.size()]), length(static_cast<int>(list.size())) {
        std::copy(list.begin(), list.end(), elements);
    }

    // Copy constructor
    template<typename U>
    Vector(const Vector<U>& other) : elements(new T[other.length]), length(other.length) {
        std::copy(other.elements, other.elements + length, elements);
    }

    Vector(const Vector<T>& other) : elements(new T[other.length]), length(other.length) {
        std::copy(other.elements, other.elements + length, elements);
    }

    // Move constructor
    template<typename U>
    Vector(Vector<U>&& other) noexcept : elements(other.elements), length(other.length) {
        other.elements = nullptr;
        other.length = 0;
    }

    Vector(Vector<T>&& other) noexcept : elements(other.elements), length(other.length) {
        other.elements = nullptr;
        other.length = 0;
    }

    // Destructor
    ~Vector() {
        delete[] elements;
        elements = nullptr;
    }

    // Assignment operators

    // Copy operator
    template<typename U>
    Vector& operator=(const Vector<U>& other) {
        if (this != &other) {
            delete[] elements;
            length = other.length;
            elements = new T[length];
            std::copy(other.elements, other.elements + length, elements);
        }
        return *this;
    }

    Vector& operator=(const Vector<T>& other) {
        if (this != &other) {
            delete[] elements;
            length = other.length;
            elements = new T[length];
            std::copy(other.elements, other.elements + length, elements);
        }
        return *this;
    }

    // Move operator
    template<typename U>
    Vector& operator=(Vector<U>&& other) noexcept {
        if (this != &other) {
            delete[] elements;
            length = other.length;
            elements = other.elements;
            other.elements = nullptr;
            other.length = 0;
        }
        return *this;
    }

    Vector& operator=(Vector<T>&& other) noexcept {
        if (this != &other) {
            delete[] elements;
            length = other.length;
            elements = other.elements;
            other.elements = nullptr;
            other.length = 0;
        }
        return *this;
    }

    // Access operators
    T& operator[](int i) {
        if (i < 0 || i >= length) throw std::out_of_range("Vector index out of bounds");
        return elements[i];
    }
    const T& operator[](int i) const {
        if (i < 0 || i >= length) throw std::out_of_range("Vector index out of bounds");
        return elements[i];
    }

    // Operator overloads
    template<typename U>
    Vector<typename std::common_type<T,U>::type> operator+(const Vector<U>& other) const {
        if (length != other.len()) throw std::invalid_argument("Vectors must be of the same length to add.");
        Vector<typename std::common_type<T,U>::type> result(length);
        for (int i = 0; i < length; i++) {
            result[i] = elements[i] + other[i];
        }
        return result;
    }

    template<typename U>
    Vector<typename std::common_type<T,U>::type> operator-(const Vector<U>& other) const {
        if (length != other.len()) throw std::invalid_argument("Vectors must be of the same length to subtract.");
        Vector<typename std::common_type<T,U>::type> result(length);
        for (int i = 0; i < length; i++) {
            result[i] = elements[i] - other[i];
        }
        return result;
    }

    template<typename U>
    Vector<typename std::common_type<T,U>::type> operator*(const U& scalar) const {
        Vector<typename std::common_type<T,U>::type> result(length);
        for (int i = 0; i < length; i++) {
            result[i] = elements[i] * scalar;
        }
        return result;
    }

    // Utility functions
    int len() const {
        return length;
    }


    void info(const std::string& txt) const {
        std::cout << "Vector " << txt << ": ";
        for (int i = 0; i < length; i++) {
            std::cout << elements[i] << " ";
        }
        std::cout << std::endl;
    }
};

// Scalar multiplication: scalar * vector
template <typename T, typename U>
Vector<typename std::common_type<T, U>::type> operator*(const T& scalar, const Vector<U>& vec) {
    return vec * scalar;
}

// Dot product
template<typename T, typename U>
typename std::common_type<T, U>::type dot(const Vector<T>& lhs, const Vector<U>& rhs) {
    if (lhs.len() != rhs.len()) throw std::invalid_argument("Vectors must be of the same length for dot product.");
    typename std::common_type<T, U>::type result = 0;
    for (int i = 0; i < lhs.len(); i++) {
        result += lhs[i] * rhs[i];
    }
    return result;
}
template <typename T>
class Matrix {
public:
    std::map<std::pair<int, int>, T> data;
    int rows, cols;

    Matrix(){
        rows = 0;
        cols = 0;
        data;
    }
    Matrix(int rows, int cols) : rows(rows), cols(cols) {}
    ~Matrix(){
        data.clear();
        rows = 0;
        cols = 0;
    }

    T& operator[](const std::pair<int, int>& ij) {
        if (ij.first >= rows || ij.second >= cols || ij.first < 0 || ij.second < 0){
            throw std::invalid_argument("Indices out of bound");
        } else {
            return data[ij];
        }
    }
    const T& operator()(const std::pair<int, int>& ij) const {
        auto it = data.find(ij);
        if (it == data.end()) throw std::runtime_error("Matrix entry not found.");
        return it->second;
    }
};
// Matrix, vector multiplication

template<typename T, typename U>
Vector<typename std::common_type<T, U>::type> operator*(const Matrix<T>& lhs, const Vector<U>& rhs){
    if(lhs.cols != rhs.len()){
        throw std::invalid_argument("Vector length must match the number of matrix columns");
    }
    Vector<typename std::common_type<T, U>::type> result(lhs.rows);
    for (auto it = lhs.data.begin(); it != lhs.data.end(); it++){
        int i = it -> first.first;
        int j = it -> first.second;
        T value = it -> second;

        result[i] += value * rhs[j];
    }
    return result;
}




template<typename T>
int cg(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x, T tol = (T)1e-8, int maxiter = 100) {
    Vector<T> r = b - A * x;
    Vector<T> p = r;
    T rsold = dot(r, r);
    Vector<T> Ap(x.len());
    for (int i = 0; i < maxiter; ++i) {
        Ap = A * p;
        T alpha = rsold / dot(p, Ap);
        x = x + p * alpha;
        r = r - Ap * alpha;
        T rsnew = dot(r, r);
        if (rsnew < tol*tol) return i;
        p = r + p * (rsnew / rsold);
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
//    Heat(T alpha, int m, T dt) : alpha(alpha), m(m), dt(dt), total_nodes(std::pow(m, n)), M(total_nodes, total_nodes) {
//        T dx = 1.0 / (m + 1);
//        if (n == 1) {
//            // Assembly for 1D case
//            for (int i = 0; i < total_nodes; ++i) {
//                if (i > 0) M[{i, i - 1}] = -alpha * dt / (dx * dx);
//                M[{i, i}] = 1 + 2 * alpha * dt / (dx * dx);
//                if (i < total_nodes - 1) M[{i, i + 1}] = -alpha * dt / (dx * dx);
//            }
//        } else if (n == 2) {
//            // Assembly for 2D case
//            for (int i = 0; i < total_nodes; ++i) {
//                if (i % m > 0) M[{i, i - 1}] = -alpha * dt / (dx * dx);  // Left
//                if (i % m < m - 1) M[{i, i + 1}] = -alpha * dt / (dx * dx);  // Right
//                if (i >= m) M[{i, i - m}] = -alpha * dt / (dx * dx);  // Top
//                if (i < total_nodes - m) M[{i, i + m}] = -alpha * dt / (dx * dx);  // Bottom
//                M[{i, i}] = 1 + 4 * alpha * dt / (dx * dx);  // Center
//            }
//        }
//    }

    Heat(T alpha, int m, T dt) : alpha(alpha), m(m), dt(dt), total_nodes(std::pow(m, n)), M(total_nodes, total_nodes) {
        T dx = 1.0 / (m + 1);
        T beta = alpha * dt / std::pow(dx, 2);

        for(int i = 0; i < total_nodes; ++i) {
            for(int k = 0; k < n; ++k){
                int j_left = i - std::pow(m, k);
                int j_right = i + std::pow(m, k);

                if(j_left >= 0 ) {
                    M[{i, i - std::pow(m, k)}] = -beta;
                }
                if(j_right < total_nodes){
                    M[{i, i + std::pow(m, k)}] = -beta;
                }
                if (i % int(std::pow(m,k+1)) == 0.0 && i > 0){
                    for (int h = 0; h < std::pow(m,k);h++)
                    {
                        M[{i + h,i - std::pow(m,k) + h}] = 0.0;
                        M[{i + h - std::pow(m,k), i + h}] = 0.0;
                    }
                }

            }
            M[{i, i}] = 1 + beta * 2 * n;
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
    Vector<T> solve(T t_final, T tol = (T)1e-8, int maxiter = 1000) const {
        Vector<T> u = exact(0);
        Vector<T> u_new(total_nodes);
        int steps = static_cast<int>(t_final / dt);
        for (int step = 0; step < steps; ++step) {
            cg(M, u, u_new, tol, maxiter);
            u = u_new;
        }
        return u;
    }
    // Method to get the matrix for verification
    const Matrix<T>& getMatrix() const {
        return M;
    }
};

template<typename T>
void printMatrix(const Matrix<T>& mat, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            try {
                T value = mat({i, j});
                std::cout << value << " ";
            } catch (const std::runtime_error&) {
                std::cout << "0 ";
            }
        }
        std::cout << std::endl;
    }
}

int main(int argc, char* argv[]) {
    bool runSimulation = true;
//    Vector<double> v1d({0.0, 1.0, 2.0});
//    Vector<double> v2d({2.0, 1.0, 0.0});
//    Vector<float> v1f({0.0f, 1.0f, 2.0f});
//    Vector<float> v2f({2.0f, 1.0f, 0.0f});
//    Vector<int> v1i({0, 1, 2});
//    Vector<int> v2i({2, 1, 0});
//
//    Vector<float> v_test = v1d;
//
//    Matrix<double> M(5, 5);
//    printMatrix(M, 4, 4);
//    M[{1,2}] = 10.0;
//    printMatrix(M, 5, 5);

    if(runSimulation) {
        const double alpha = 0.3125;
        const int m = 99;
        const double dt = 0.001;
        const double t_final_1D = 1.0;
        const double t_final_2D = 0.5;

        auto start = std::chrono::high_resolution_clock::now();

        Heat<1, double> solver1D(alpha, m, dt);
        auto solution1D = solver1D.solve(t_final_1D);
        auto exactSolution1D = solver1D.exact(t_final_1D);
//        for (int i = 0; i < solution1D.len(); ++i) {
//            std::cout << "Numerical 1D: " << solution1D[i] << ", Exact 1D: " << exactSolution1D[i] << std::endl;
//        }

//     Print 1D Matrix for verification
//     printMatrix(solver1D.getMatrix(), m, m);

        Heat<2, double> solver2D(alpha, m, dt);
        auto solution2D = solver2D.solve(t_final_2D);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        auto exactSolution2D = solver2D.exact(t_final_2D);
        std::cout << duration.count() << std::endl;
//        std::cout << "Numerical at (0,0): " << numerical[0] << "| Exact at (0.0): " << exact[0] << std::endl;


//        for (int i = 0; i < solution2D.len(); ++i) {
//            std::cout << "Numerical 2D: " << solution2D[i] << ", Exact 2D: " << exactSolution2D[i] << std::endl;
//        }

        // Print 2D Matrix for verification
//     printMatrix(solver2D.getMatrix(), m * m, m * m);
    }
    return 0;
}

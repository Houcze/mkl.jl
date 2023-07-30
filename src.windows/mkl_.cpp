#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <mkl.h>
#include <algorithm>
#include <cmath>

double round_if_close(double num, double epsilon = 1e-5)
{
    double nearest = std::round(num);
    if (std::abs(num - nearest) < epsilon)
    {
        return nearest;
    }
    else
    {
        return num;
    }
}

std::vector<double> simplify(std::vector<double> vec, double precision = 1e-6)
{
    double minElement = *std::min_element(vec.begin(), vec.end(), [precision](double a, double b)
                                          { return std::abs(a) < std::abs(b) && std::abs(a) > precision && a != 0.0; });

    if (minElement == 0.0)
        return vec;

    std::transform(vec.begin(), vec.end(), vec.begin(), [minElement](double a)
                   { return a / minElement; });

    return vec;
}

std::vector<double> genfromtxt(const std::string &filename)
{
    std::ifstream file(filename);
    std::vector<double> matrix;
    double value;
    while (file >> value)
    {
        matrix.push_back(value);
    }
    return matrix;
}

int diagonalize_p(std::vector<double> a_)
{
    int n = sqrt(a_.size());
    std::vector<double> a(n * n);
    std::copy(a_.begin(), a_.end(), a.begin());

    std::vector<double> wr(n), wi(n);
    std::vector<double> vl(n * n), vr(n * n);
    std::vector<int> ipiv(n);

    int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'V', 'V', n, a.data(), n, wr.data(), wi.data(), vl.data(), n, vr.data(), n);

    if (info > 0)
    {
        std::cout << "Failure in solving eigen problem. Info = " << info << std::endl;
        return info;
    }

    std::vector<std::pair<double, std::vector<double>>> eigen;

    for (int i = 0; i < n; i++)
    {
        std::vector<double> vec(n);
        for (int j = 0; j < n; j++)
        {
            vec[j] = vr[j * n + i];
        }
        eigen.push_back({wr[i], vec});
    }

    std::sort(eigen.begin(), eigen.end(),
              [](const auto &a, const auto &b)
              { return a.first < b.first; });

    std::cout << "Sorted eigenvalues and eigenvectors:\n";

    for (const auto &pair : eigen)
    {
        std::cout << "Eigenvalue: " << pair.first << "\n";
        break;
    }

    std::vector<double> P(n * n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            P[i * n + j] = eigen[j].second[i];
        }
    }
    /*
    std::vector<double> P_copy = P;

    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, P.data(), n, ipiv.data());
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, P.data(), n, ipiv.data());
    std::vector<double> a_recon = a_;

    std::vector<double> temp(n * n);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1, P.data(), n, a_recon.data(), n, 0, temp.data(), n);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1, temp.data(), n, P_copy.data(), n, 0, a_recon.data(), n);
    */
    return EXIT_SUCCESS;
}

double *eigvals_p(std::vector<double> a_)
{
    int n = sqrt(a_.size());
    std::vector<double> a(n * n);
    std::copy(a_.begin(), a_.end(), a.begin());

    std::vector<double> wr(n), wi(n);
    std::vector<double> vl(n * n), vr(n * n);
    std::vector<int> ipiv(n);

    int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'V', 'V', n, a.data(), n, wr.data(), wi.data(), vl.data(), n, vr.data(), n);

    if (info > 0)
    {
        std::cout << "Failure in solving eigen problem. Info = " << info << std::endl;
        return nullptr;
    }

    std::vector<std::pair<double, std::vector<double>>> eigen;

    for (int i = 0; i < n; i++)
    {
        std::vector<double> vec(n);
        for (int j = 0; j < n; j++)
        {
            vec[j] = vr[j * n + i];
        }
        eigen.push_back({wr[i], vec});
    }

    std::sort(eigen.begin(), eigen.end(),
              [](const auto &a, const auto &b)
              { return a.first < b.first; });

    std::vector<double> eigvalues_;
    for (const auto &pair : eigen)
    {
        eigvalues_.push_back(pair.first);
    }

    double *eigvalues = new double[eigvalues_.size()];
    std::copy(eigvalues_.begin(), eigvalues_.end(), eigvalues);
    return eigvalues;
}

double *eigvecs_p(std::vector<double> a_)
{
    int n = sqrt(a_.size());
    std::vector<double> a(n * n);
    std::copy(a_.begin(), a_.end(), a.begin());

    std::vector<double> wr(n), wi(n);
    std::vector<double> vl(n * n), vr(n * n);
    std::vector<int> ipiv(n);

    int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'V', 'V', n, a.data(), n, wr.data(), wi.data(), vl.data(), n, vr.data(), n);

    if (info > 0)
    {
        std::cout << "Failure in solving eigen problem. Info = " << info << std::endl;
        return nullptr;
    }

    std::vector<std::pair<double, std::vector<double>>> eigen;

    for (int i = 0; i < n; i++)
    {
        std::vector<double> vec(n);
        for (int j = 0; j < n; j++)
        {
            vec[j] = vr[j * n + i];
        }
        eigen.push_back({wr[i], vec});
    }

    std::sort(eigen.begin(), eigen.end(),
              [](const auto &a, const auto &b)
              { return a.first < b.first; });

    double *eigvecs = new double[n * n];

    for (int i = 0; i < n; i++)
    {
        std::copy(eigen[i].second.begin(), eigen[i].second.end(), &eigvecs[i * n]);
    }

    return eigvecs;
}

extern "C"
{
    __declspec(dllexport) void diagonalize(double *a, int n)
    {
        std::vector<double> a_(a, a + n * n);
        diagonalize_p(a_);
    }

    __declspec(dllexport) double *eigvals(double *a, int n)
    {
        std::vector<double> a_(a, a + n * n);
        double *eigvals = eigvals_p(a_);
        return eigvals;
    }

    __declspec(dllexport) double *eigvecs(double *a, int n)
    {
        std::vector<double> a_(a, a + n * n);
        double *eigvecs = eigvecs_p(a_);
        return eigvecs;
    }

    __declspec(dllexport) double *matmul(double *A, double *B, int m, int k, int n)
    {
        double *C = new double[m * n];
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, A, k, B, n, 0.0, C, n);
        return C;
    }

    __declspec(dllexport) double *cholesky(double *A, int n)
    {
        double *B = new double[n * n];

        std::copy(A, A + n * n, B);

        int info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', n, B, n);
        if (info != 0)
        {
            delete[] B;
            throw std::runtime_error("Cholesky decomposition failed.");
        }
        return B;
    }

    __declspec(dllexport) double *inverse(double *A, int n)
    {
        double *A_ = new double[n * n];
        std::copy(A, A + n * n, A_);

        int *ipiv = new int[n];

        int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A_, n, ipiv);
        if (info == 0)
        {
            info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, A_, n, ipiv);
        }

        delete[] ipiv;

        if (info != 0)
        {
            delete[] A_;
            return nullptr;
        }

        return A_;
    }

    __declspec(dllexport) void kill(double *a)
    {
        delete[] a;
    }
}

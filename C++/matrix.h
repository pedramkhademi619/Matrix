#ifndef MATRIX_H
#define MATRIX_H
#define COLON ,


#include <iostream>
#include "utils.h"
#include <cmath>
#include <vector>


#define ET double

class Matrix{

private:
    int nrows, ncols;

    std::vector<ET> data;

public:

    // default constructor
    inline Matrix() : nrows(0), ncols(0) {}

    // Parameterized constructor
    inline Matrix(int nrows, int ncols){
        minit(nrows, ncols);
    }

    // copy constructor
    inline Matrix(const Matrix & m){
        minit(m.size(1), m.size(2));
        data = m.data;
    }

    // Constructor with a single value
    Matrix(const ET& val){
        minit(1,1);
        csi0(0) = val;
    }

    inline Matrix(const std::vector<ET> &arrdata, int nrows, int ncols){
        minit(nrows, ncols);
        load(arrdata);
    }

    Matrix createSubmatrix(int excludedRow, int excludedCol) const {
        Matrix submatrix(nrows - 1, ncols - 1);
        int subRow = 0, subCol = 0;

        for (int i = 0; i < nrows; i++) {
            if (i != excludedRow) {
                subCol = 0;
                for (int j = 0; j < ncols; j++) {
                    if (j != excludedCol) {
                        submatrix.cij0(subRow, subCol) = cij0(i, j);
                        subCol++;
                    }
                }
                subRow++;
            }
        }

        return submatrix;
    }



    inline void minit(int nr, int nc) {
        nrows = nr;
        ncols = nc;

        // Clear the existing data and resize the vector
        data.clear();
        data.resize(numel(), ET(0));
    }

    inline Matrix& operator=(const Matrix& m) {
        if (this == &m)
            return *this;

        if (numel() == m.numel()) {
            nrows = m.size(1);
            ncols = m.size(2);
            data = m.data;  // انتقال داده با استفاده از std::vector
        } else {
            Matrix tmp(m);
            swap(*this, tmp);
        }

        return *this;
    }


    inline void swap(Matrix& a, Matrix& b) {
        using std::swap;  // Import std::swap into the scope

        swap(a.nrows, b.nrows);
        swap(a.ncols, b.ncols);
        swap(a.data, b.data);
    }




    virtual inline ET& csi0(int i) const {
        return const_cast<ET&>(data[i]);
    }


    virtual inline ET& cij0(int i, int j) const {
        return const_cast<ET&>(data[i + j * nrows]);
    }


    void load(const std::vector<ET> &arrdata) {
        int i, j;
        for (i = 0; i < nrows; i++) {
            for (j = 0; j < ncols; j++) {
                this->cij0(i, j) = arrdata[j * this->nrows + i];
            }
        }
    }

    void print(std::ostream& os) const {
        os << "--------- Matrix("
           << this->nrows << "x"
           << this->ncols <<
           ") --------" << std::endl;
        int i, j;
        for(i = 0; i < nrows; i++){
            for(j = 0; j < ncols; j++){
                os << this->cij0(i,j) << " ";
            }
            os << std::endl;
        }
        os << "------------------------------" << std::endl;
    }

    friend std::ostream & operator<<(std::ostream& os, const Matrix& m){
        m.print(os);
        return os;
    }


    inline int size(int d) const{ return (d==1) ? nrows : ncols; }


    int numel() const{ return nrows*ncols; }

    bool isscalar() const{
        return (numel() == 1) ? 1 : 0;
    }
    bool isempty() const{
        return (numel() == 0) ? 1 : 0;
    }
    inline bool samesize(const Matrix& m) const{
        return nrows == m.size(1) && ncols == m.size(2);
    }

    inline ~Matrix() {


    }

    //operators 😒️
    // unary + operator
    Matrix operator+() const{
        Matrix c(size(1),size(2));
        for (int i = 0 ; i < numel(); i++)
            c.csi0(i) = +csi0(i);
        return c;
    }

    // binary + operator
    friend Matrix operator+(const Matrix & a, const Matrix & b) {
        if(a.isscalar()){
            Matrix c(b.size(1),b.size(2));
            ET v0 = (ET)a.csi0(0);
            for (int i = 0 ; i < b.numel(); i++)
                c.csi0(i) = v0 + b.csi0(i);
            return c;
        }else if(b.isscalar()){
            Matrix c(a.size(1),a.size(2));
            ET v0 = (ET)b.csi0(0);
            for (int i = 0 ; i < a.numel(); i++)
                c.csi0(i) = a.csi0(i) + v0;
            return c;
        }else if(a.samesize(b)){
            Matrix c(a.size(1),a.size(2));
            for (int i = 0 ; i < a.numel(); i++)
                c.csi0(i) = a.csi0(i) + b.csi0(i);
            return c;
        }else
            ERROR("Matrix dimension mismatch in operation '+'.");
        return Matrix::empty(); // never is reached
    }

    friend inline Matrix operator+(const Matrix & a, const double& b)
    { return a + Matrix(b); }
    friend inline Matrix operator+(const double& a, const Matrix & b)
    { return Matrix(a) + b; }


    // unary - operator
    Matrix operator-() const{
        Matrix c(size(1),size(2));
        for (int i = 0 ; i < numel(); i++)
            c.csi0(i) = -csi0(i);
        return c;
    }


    // binary - operator
    friend Matrix operator-(const Matrix & a, const Matrix & b) {
        if(a.isscalar()){
            Matrix c(b.size(1),b.size(2));
            ET v0 = (ET)a.csi0(0);
            for (int i = 0 ; i < b.numel(); i++)
                c.csi0(i) = v0 - b.csi0(i);
            return c;
        }else if(b.isscalar()){
            Matrix c(a.size(1),a.size(2));
            ET v0 = (ET)b.csi0(0);
            for (int i = 0 ; i < a.numel(); i++)
                c.csi0(i) = a.csi0(i) - v0;
            return c;
        }else if(a.samesize(b)){
            Matrix c(a.size(1),a.size(2));
            for (int i = 0 ; i < a.numel(); i++)
                c.csi0(i) = a.csi0(i) - b.csi0(i);
            return c;
        }else
        ERROR("Matrix dimension mismatch in operation '-'.");
        return Matrix::empty(); // never is reached
    }

    friend inline Matrix operator-(const Matrix & a, const double& b){

        return (a - Matrix(b));
    }
    friend inline Matrix operator-(const double& a, const Matrix & b){

        return (Matrix(a) - b);
    }



    // binary * operator
    friend Matrix operator*(const Matrix& a, const Matrix& b) {
        Matrix result;

        if (a.isscalar()) {
            result = Matrix(b.size(1), b.size(2));
            ET v1 = static_cast<ET>(a.csi0(0));
            for (int i = 0; i < b.numel(); i++)
                result.csi0(i) = v1 * b.csi0(i);
        } else if (b.isscalar()) {
            result = Matrix(a.size(1), a.size(2));
            ET v2 = static_cast<ET>(b.csi0(0));
            for (int i = 0; i < a.numel(); i++)
                result.csi0(i) = a.csi0(i) * v2;
        } else if (a.samesize(b)) {
            result = Matrix(a.size(1), a.size(2));
            for (int i = 0; i < a.numel(); i++)
                result.csi0(i) = a.csi0(i) * b.csi0(i);
        } else {
            ERROR("Matrix dimension mismatch in operation '*'.");
            return Matrix::empty();
        }

        return result;
    }

    friend inline Matrix operator*(const Matrix & a, const double& b){
    Matrix c = Matrix(a.size(1), a.size(2));
    for (int i = 0; i  < a.numel(); i++){
        c.csi0(i)  = a.csi0(i) * b;
    }
    return c;
    }

    friend inline Matrix operator*(const double& a, const Matrix & b){
        Matrix c = Matrix(b.size(1), b.size(2));
        for (int i = 0; i  < b.numel(); i++){
            c.csi0(i)  = b.csi0(i) * a;
        }
        return c;
    }

    // Binary / operator for matrix division
    friend Matrix operator/(const Matrix& a, const Matrix& b) {
        ASSERT(b.nrows == b.ncols, "Matrix B must be square for division.");
        Matrix result;

        if (b.numel() == 1 && b.csi0(0) == 0) {
            ERROR("Matrix division by zero.");
            return Matrix::empty();
        }

        result = a * b.inverse();
        return result;
    }
// Binary / operator for scalar division
    friend Matrix operator/(const Matrix& a, const ET& scalar) {
        if (scalar == 0) {
            ERROR("Scalar division by zero.");
            return Matrix::empty();
        }

        Matrix result(a.nrows, a.ncols);
        for (int i = 0; i < a.numel(); i++) {
            result.csi0(i) = a.csi0(i) / scalar;
        }
        return result;
    }

// Binary / operator for scalar division with scalar on the right
    friend Matrix operator/(const ET& scalar, const Matrix& b) {
        if (b.numel() == 1 && b.csi0(0) == 0) {
            ERROR("Matrix division by zero.");
            return Matrix::empty();
        }

        Matrix result(b.nrows, b.ncols);
        for (int i = 0; i < b.numel(); i++) {
            result.csi0(i) = scalar / b.csi0(i);
        }
        return result;
    }


    // Function to return equal larger numbers
    Matrix operator>=(const Matrix &other) const{
        Matrix result(nrows, ncols);
        for (int i =0; i < numel(); i++){
            result.csi0(i) = (csi0(i) >= other.csi0(i)) ? 1 : 0;
        }
        return result;
    }

    // Function to return equal smaller numbers
    Matrix operator<=(const Matrix &other) const{
        Matrix result(nrows, ncols);
        for (int i =0; i < numel(); i++){
            result.csi0(i) = (csi0(i) <= other.csi0(i)) ? 1 : 0;
        }
        return result;
    }

    // Equality return function
    Matrix operator==(const Matrix &other) const {
        Matrix result(nrows, ncols);
        for (int i = 0; i < numel(); i++) {
            result.csi0(i) = (csi0(i) == other.csi0(i)) ? 1 : 0;
        }
        return result;
    }

    Matrix operator!=(const Matrix &other) const {
        Matrix result(nrows, ncols);
        for (int i = 0; i < numel(); i++) {
            result.csi0(i) = (csi0(i) != other.csi0(i)) ? 1 : 0;
        }
        return result;
    }

    // function to return indices that meet the combination conditions
    Matrix operator&(const Matrix &other) const {
        ASSERT(samesize(other), "Matrix dimensions must be the same for elementwise AND operation.");
        Matrix result(nrows, ncols);
        for (int i = 0; i < numel(); i++) {
            result.csi0(i) = (csi0(i) && other.csi0(i)) ? 1 : 0;
        }
        return result;
    }
    // function to return indices that meet at least one of the conditions
    Matrix operator|(const Matrix &other) const {
        ASSERT(samesize(other), "Matrix dimensions must be the same for elementwise OR operation.");
        Matrix result(nrows, ncols);
        for (int i = 0; i < numel(); i++) {
            result.csi0(i) = (csi0(i) || other.csi0(i)) ? 1 : 0;
        }
        return result;
    }







    // axillury functions
    void range(ET start, ET end, ET step = 1){
        int ub = (end - start) / step;
        for(int i = 0; i < numel() && i < ub; i++){
            this->csi0(i) = start;
            start+=step;
        }
    }
    Matrix power(int exponent) const{
        if (exponent == 0){
            return eye(nrows);
        }
        else if (exponent > 0){
            Matrix result = *this;
            for (int i = 0; i < exponent; i++){
                result = result * (*this);
            }

        }
        else{
            inverse().power(-exponent);
        }
    }
    friend Matrix elementwisePower(const Matrix& m, const Matrix& exponent) {
        ASSERT(m.samesize(exponent), "Matrix dimension mismatch in elementwise power.");

        Matrix result(m.nrows, m.ncols);
        for (int i = 0; i < m.numel(); i++) {
            result.csi0(i) = std::pow(m.csi0(i), exponent.csi0(i));
        }

        return result;
    }

    Matrix cos() const {
        Matrix result(nrows, ncols);
        for (int i = 0; i < numel(); i++) {
            result.csi0(i) = std::cos(csi0(i));
        }
        return result;
    }
    Matrix tan() const {
        Matrix result(nrows, ncols);
        for (int i = 0; i < numel(); i++) {
            result.csi0(i) = std::tan(csi0(i));
        }
        return result;
    }

    Matrix exp() const {
        Matrix result(nrows, ncols);
        for (int i = 0; i < numel(); i++) {
            result.csi0(i) = std::exp(csi0(i));
        }
        return result;
    }

    //Statistical functions
    ET min() const {
        ET minValue = csi0(0);
        for (int i = 1; i < numel(); i++) {
            minValue = std::min(minValue, csi0(i));
        }
        return minValue;
    }

    ET max() const {
        ET maxValue = csi0(0);
        for (int i = 1; i < numel(); i++) {
            maxValue = std::max(maxValue, csi0(i));
        }
        return maxValue;
    }

    ET mean() const {
        ET sum = 0;
        for (int i = 0; i < numel(); i++) {
            sum += csi0(i);
        }
        return sum / numel();
    }

    Matrix sumRows() const {
        Matrix result(1, ncols);
        for (int j = 0; j < ncols; j++) {
            ET sum = 0;
            for (int i = 0; i < nrows; i++) {
                sum += cij0(i, j);
            }
            result.cij0(0, j) = sum;
        }
        return result;
    }

    // تابع برگشت مجموع اعداد ستونی
    Matrix sumCols() const {
        Matrix result(nrows, 1);
        for (int i = 0; i < nrows; i++) {
            ET sum = 0;
            for (int j = 0; j < ncols; j++) {
                sum += cij0(i, j);
            }
            result.cij0(i, 0) = sum;
        }
        return result;
    }

    // تابع برگشت ضرب اعداد سطری
    Matrix productRows() const {
        Matrix result(1, ncols);
        for (int j = 0; j < ncols; j++) {
            ET product = 1;
            for (int i = 0; i < nrows; i++) {
                product *= cij0(i, j);
            }
            result.cij0(0, j) = product;
        }
        return result;
    }

    // تابع برگشت ضرب اعداد ستونی
    Matrix productCols() const {
        Matrix result(nrows, 1);
        for (int i = 0; i < nrows; i++) {
            ET product = 1;
            for (int j = 0; j < ncols; j++) {
                product *= cij0(i, j);
            }
            result.cij0(i, 0) = product;
        }
        return result;
    }

    ET variance() const {
        ET meanValue = mean();
        ET sumSquaredDeviations = 0;
        for (int i = 0; i < numel(); i++) {
            ET deviation = csi0(i) - meanValue;
            sumSquaredDeviations += deviation * deviation;
        }
        return sumSquaredDeviations / numel();
    }

    ET std() const {
        return std::sqrt(variance());
    }
    Matrix columnVar() const {
        Matrix result(1, ncols);
        for (int j = 0; j < ncols; j++) {
            ET meanValue = sumCols().csi0(j) / nrows;
            ET sumSquaredDeviations = 0;
            for (int i = 0; i < nrows; i++) {
                ET deviation = cij0(i, j) - meanValue;
                sumSquaredDeviations += deviation * deviation;
            }
            result.cij0(0, j) = sumSquaredDeviations / nrows;
        }
        return result;
    }

    Matrix columnStdDev() const {
        Matrix result(1, ncols);
        for (int j = 0; j < ncols; j++) {
            ET meanValue = sumCols().csi0(j) / nrows;
            ET sumSquaredDeviations = 0;
            for (int i = 0; i < nrows; i++) {
                ET deviation = cij0(i, j) - meanValue;
                sumSquaredDeviations += deviation * deviation;
            }
            result.cij0(0, j) = std::sqrt(sumSquaredDeviations / nrows);
        }
        return result;
    }

    Matrix rowVar() const {
        Matrix result(nrows, 1);
        for (int i = 0; i < nrows; i++) {
            ET meanValue = sumRows().csi0(i) / ncols;
            ET sumSquaredDeviations = 0;
            for (int j = 0; j < ncols; j++) {
                ET deviation = cij0(i, j) - meanValue;
                sumSquaredDeviations += deviation * deviation;
            }
            result.cij0(i, 0) = sumSquaredDeviations / ncols;
        }
        return result;
    }

    // تابع برگشت انحراف معیار هر سطر
    Matrix rowStdDev() const {
        Matrix result(nrows, 1);
        for (int i = 0; i < nrows; i++) {
            ET meanValue = sumRows().csi0(i) / ncols;
            ET sumSquaredDeviations = 0;
            for (int j = 0; j < ncols; j++) {
                ET deviation = cij0(i, j) - meanValue;
                sumSquaredDeviations += deviation * deviation;
            }
            result.cij0(i, 0) = std::sqrt(sumSquaredDeviations / ncols);
        }
        return result;
    }



    // unimplemented Auxiliary functions
//    int sum(){}


    /**********************/
    /*** static methods ***/
    /**********************/
    static inline const Matrix empty() {
        return Matrix(0,0);
    }

    static Matrix eye(const int& e){
        Matrix retmat(e,e);
        int i, j;
        for(i = 0; i < e; i++){
            for(j = 0; j < e; j++){
                retmat.cij0(i,j) = 0;
                if (i == j)
                    retmat.cij0(i,j) = 1;
            }
        }
        return retmat;
    }

    static Matrix rand(int nrows, int ncols){
        Matrix retmat(nrows, ncols);
        int i;
        for(i = 0; i < retmat.numel(); i++){
            retmat.csi0(i) = std::rand() / ET(RAND_MAX);
        }
        return retmat;
    }

    static Matrix abs(const Matrix& m){
        Matrix retmat(m.size(1), m.size(2));
        int i;
        for(i = 0; i < retmat.numel(); i++){
            retmat.csi0(i) = std::abs(m.csi0(i));
        }
        return retmat;
    }
    static Matrix randn(int nrows, int ncols) {
        Matrix retmat(nrows, ncols);
        int i;
        for (i = 0; i < retmat.numel(); i++) {
            double u1 = (std::rand() + 1.0) / (RAND_MAX + 1.0);  // To avoid log(0)
            double u2 = (std::rand() + 1.0) / (RAND_MAX + 1.0);
            double rand_std_normal = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
            retmat.csi0(i) = rand_std_normal;
        }
        return retmat;
    }
    ET determinant() const {
        ASSERT(nrows == ncols, "Matrix must be square for determinant calculation.");

        if (nrows == 1) {
            return csi0(0);
        } else if (nrows == 2) {
            return cij0(0, 0) * cij0(1, 1) - cij0(0, 1) * cij0(1, 0);
        } else {
            ET det = 0;
            for (int i = 0; i < ncols; i++) {
                Matrix submatrix = createSubmatrix(0, i);
                det += (i % 2 == 0 ? 1 : -1) * cij0(0, i) * submatrix.determinant();
            }
            return det;
        }
    }

    // Function to find the inverse of the matrix
    Matrix inverse()  const {
        ASSERT(nrows == ncols, "Matrix must be square for inverse calculation.");

        ET det = determinant();
        ASSERT(det != 0, "Matrix is singular, and its inverse does not exist.");

        Matrix invMatrix(nrows, ncols);

        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                Matrix submatrix = createSubmatrix(i, j);
                invMatrix.cij0(i, j) = ((i + j) % 2 == 0 ? 1 : -1) * submatrix.determinant() / det;
            }
        }

        return invMatrix.transpose();
    }

    Matrix Pinv()  {
        Matrix transposed = transpose();
        Matrix result = transposed * (*this);
        return result.inverse() * transposed;
    }

    Matrix transpose()  {
        Matrix transposed(ncols, nrows);
        for (int i = 0; i < nrows; i++) {
            for (int j = 0; j < ncols; j++) {
                transposed.cij0(j, i) = cij0(i, j);
            }
        }
        return transposed;
    }

    //indices

    //num indexing
    inline ET& operator()(int i) const {
        ASSERT(i < numel(), "Index out of bounds.");
        return csi0(i);
    }

    //vector indexing
    Matrix operator()(const std::vector<int>& indices) const {
        Matrix result(indices.size(), 1);
        for (size_t i = 0; i < indices.size(); ++i) {
            int index = indices[i];
            ASSERT(index < numel(), "Index out of bounds.");
            result.csi0(i) = csi0(index);
        }
        return result;
    }

//    //range indexing
//    Matrix operator()(const std::string& indices) const {
//        size_t colonPos = indices.find(':');
//        if (colonPos == std::string::npos) {
//            // No colon found, treat it as a single index
//            int index = std::stoi(indices);
//            ASSERT(index < numel(), "Index out of bounds.");
//            return Matrix(csi0(index));
//        } else {
//            // Colon found, treat it as a range
//            int start = std::stoi(indices.substr(0, colonPos));
//            int end = std::stoi(indices.substr(colonPos + 1));
//            ASSERT(start <= end, "Invalid range: start should be less than or equal to end.");
//            ASSERT(end < numel(), "Index out of bounds.");
//
//            int rangeSize = end - start + 1;
//            Matrix result(rangeSize, 1);
//
//            for (int i = 0; i < rangeSize; ++i) {
//                result.csi0(i) = csi0(start + i);
//            }
//
//            return result;
//        }
//    }

    Matrix operator()(int start, int end) const {
        ASSERT(start >= 1 && end <= numel() && start <= end, "Invalid range indices.");

        int size = end - start + 1;
        Matrix result(size, 1);

        for (int i = 0; i < size; ++i) {
            result.csi0(i) = csi0(start + i - 1);
        }

        return result;
    }






};







#endif // MATRIX_H

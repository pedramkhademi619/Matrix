#ifndef MATRIX_H
#define MATRIX_H

#include <cstdlib>
#include <iostream>
#include "utils.h"

using namespace std;

#define ET float

class Matrix{

private:
    int nrows, ncols;
    unsigned char attribs;
    ET *data;

public:
    // default constructor
    inline Matrix(){
        this->nrows = 0;
        this->ncols = 0;
        this->data = NULL;
    }

    inline Matrix(int nrows, int ncols){
        minit(nrows, ncols);
    }

    // copy constructor
    inline Matrix(const Matrix & m){
        minit(m.size(1), m.size(2));
        for(int i = 0; i < numel(); i++){
            csi0(i) = m.csi0(i);
        }
    }

    Matrix(const ET& val){
        minit(1,1);
        csi0(0) = val;
    }

    inline void minit(int nr, int nc){
        nrows = nr;
        ncols = nc;
        if(!isempty()){
            if (numel()){
                data = (ET*) malloc(sizeof(ET)*numel());
                ASSERT(data != NULL, "minit Cannot allocate memory! ");
            } else {
                data = (ET*) NULL;
            }
        }else{
            data = (ET*)(void *)0;
        }
    }


    inline Matrix& operator=(const Matrix& m){
        if(this == &m)
            return *this;
        if(numel() == m.numel()){
            nrows = m.size(1);
            ncols = m.size(2);
            attribs = m.attribs;
            for(int i = 0; i < m.numel(); i++){
                csi0(i) = m.csi0(i);
            }
        }else{
            Matrix tmp(m);
            swap(this,&tmp);
        }
        return *this;
    }

    inline void swap(Matrix* a, Matrix* b){
        int tmpi;
        tmpi = a->nrows; a->nrows = b->nrows; b->nrows = tmpi;
        tmpi = a->ncols; a->ncols = b->ncols; b->ncols = tmpi;

        ET* tmpr;
        tmpr = a->data; a->data = b->data; b->data = tmpr;

        unsigned char tmpa;
        tmpa = a->attribs; a->attribs = b->attribs; b->attribs = tmpa;
    }



    virtual inline ET& csi0(int i) const{
        return data[i];
    }

    virtual inline ET& cij0(int i , int j) const{
        return data[i+j*nrows];
    }

    void load(ET *arrdata){
        int i, j;
        for(i = 0; i < nrows; i++){
            for(j = 0; j < ncols; j++){
                this->cij0(i,j) = arrdata[j*this->nrows+i];
            }
        }
    }

    void print(ostream& os) const {
        os << "--------- Matrix("
           << this->nrows << "x"
           << this->ncols <<
           ") --------" << std::endl;
        int i, j;
        for(i = 0; i < nrows; i++){
            for(j = 0; j < ncols; j++){
                os << this->cij0(i,j) << " ";
            }
            os << endl;
        }
        os << "------------------------------" << std::endl;
    }

    friend ostream & operator<<(ostream& os, const Matrix& m){
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

    inline ~Matrix(){
        free(this->data);
    }

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
    friend Matrix operator*(const Matrix & a, const Matrix & b) {
        if (a.isscalar()){
            Matrix c = Matrix(b.size(1), b.size(2));
            ET v1 = (ET)a.csi0(0);
            for (int i = 0; i < b.numel(); i++){
                b.csi0(i) *= v1;
            }
        }
        else if (b.isscalar()) {
            Matrix c = Matrix(a.size(1), a.size(2));
            ET v2 = (ET) b.csi0(0);
            for (int i = 0; i < a.numel(); i++) {
                a.csi0(i) *= v2;
            }
        }
        else {
            ASSERT(a.size(2) == b.size(1), "Matrix dimension mismatch in operation '*'.");
            Matrix c = Matrix(a.size(1), b.size(2));
            for (int i = 0; i < a.size(1); i++) {
                for (int j = 0; j < b.size(2); j++) {
                    c.cij0(i, j) = 0;
                    for (int k = 0; k < b.size(1); k++) {
                        c.cij0(i, j) +=  a.cij0(i, k) * b.cij0(k, j);


                    }
                }
            }
            return c;
        }
        return Matrix::empty();
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


    // axillury functions
    void range(ET start, ET end, ET step = 1){
        int ub = (end - start) / step;
        for(int i = 0; i < numel() && i < ub; i++){
            this->csi0(i) = start;
            start+=step;
        }
    }

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
};


#endif // MATRIX_H


#include "matrix.h"
#include "utils.h"

int test_matrix(){
    std::cout << __func__ << " Started!" << std::endl << std::endl;

    // test default matrix constructor
    Matrix mat0;
    std::cout << mat0 << std::endl;

    // test matrix construction
    Matrix mat(3, 4);
    mat.range(0, 12);
    std::cout << mat << std::endl;

    printvar(mat.cij0(1,1))

    // unary + operator
    std::cout << +mat << std::endl;

#define printvarnl(x) std::cout << #x << " \n" \
    << (x) << std::endl<< std::endl;

    // add value to matrix
    Matrix meye = Matrix::eye(3);
    printvarnl(meye)
    printvarnl(meye + 1)
    printvarnl(2 + meye)

    // add matrix to matrix
    Matrix mat3 = mat + mat;
    std::cout << mat3 << std::endl;

    // unary - operator
    std::cout << -mat3 << std::endl;

    // binary - operator
    std::cout << mat3-mat3 << std::endl;
    std::cout << 3-mat3 << std::endl;
    std::cout << mat3-3 << std::endl;

    // value to matrix multiplication
    std::cout << 3*mat << std::endl;
    std::cout << mat*3 << std::endl;

    // matrix to matrix multiplication
    Matrix matt(4,3);
    matt.range(0,12);
    std::cout << mat*matt << std::endl;
    std::cout << __func__ << " Ended!" << std::endl << std::endl;
    return 0;
}


int main(){
     test_matrix();
     return 0;
}
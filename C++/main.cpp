
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
    std::cout<<mat.transpose()<<std::endl;

    // matrix to matrix multiplication
    Matrix matt(4,3);
    matt.range(0,12);
    std::cout << mat*matt << std::endl;

//    std::vector<ET> arr = {1.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0, 1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0};
//    std::vector<ET> arr10 = {1.0, .45, .67, .67, 1., .33, .33, .55, 1.};
    std::vector<ET> arr10 = {1.,0.,0.,1.};
// Create a Matrix using the new constructor
    Matrix mat45(arr10, 2,2);
    std::cout<<mat45 * mat45;

    // Print the matrix
//    cout << mat45.transpose() << "\n\n\n\n\n\n\n"<< endl;
//    cout << mat45.inverse() << "\n\n\n\n\n\n\n"<< endl;
//    cout<<mat45(8)<<endl;

    Matrix myMatrix(4, 4);
    myMatrix.range(1, 16); // مقداردهی اولیه با اعداد از 1 تا 16

    std::vector<int> indices = {0, 4, 5, 8,14}; // اندیس‌های مورد نظر
    Matrix result = myMatrix({0,1,2,3,8,5,4});

    std::cout << result << std::endl; // چاپ ماتریس حاصل


    std::cout<<myMatrix(3 COLON 9);
    std::cout << __func__ << " Ended!" << std::endl << std::endl;




    return 0;
}


int main(){
     test_matrix();
     return 0;
}
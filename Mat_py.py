import random
import math


class Matrix:
    def __init__(self, data=list()):
        self.data = data
        self.m = len(data)
        try:
            self.n = len(data[0])
        except IndexError:
            self.n = 0

    def ones(self, m: int, n: int, dtype="float"):
        return Matrix([[1 for i in range(n)] for j in range(m)])

    def zeros(self, m: int, n: int, dtype="float"):
        return Matrix([[0 for i in range(n)] for j in range(m)])

    def get_minor(self, row, col):
        # Get the minor of the matrix by removing the specified row and column
        minor_data = [
            row[:col] + row[col + 1 :]
            for row in (self.data[:row] + self.data[row + 1 :])
        ]
        return Matrix(minor_data)

    @property
    def determinant(self):
        # Calculate the determinant of a 2x2 matrix
        if self.m == 2 and self.n == 2:
            return self.data[0][0] * self.data[1][1] - self.data[0][1] * self.data[1][0]
        # Calculate the determinant of a larger matrix using cofactor expansion
        elif self.m == self.n:
            det = 0
            for i in range(self.n):
                minor = self.get_minor(0, i)
                det += self.data[0][i] * ((-1) ** i) * minor.determinant()
            return det
        else:
            raise ValueError("Matrix should be square for determinant calculation.")

    @property
    def transpose(self):
        # Transpose the matrix (swap rows and columns)
        transposed_data = [
            [self.data[j][i] for j in range(self.m)] for i in range(self.n)
        ]
        return Matrix(transposed_data)

    def eye(self, n):
        identity = [[0] * n for _ in range(n)]
        for i in range(n):
            identity[i][i] = 1
        return Matrix(identity)

    @property
    def pinv(self):
        pseudo_inverse = self.zeros(self.m, self.n)
        transpose_matrix = [[row[i] for row in self.data] for i in range(len(self.n))]
        for i in range(len(transpose_matrix)):
            for j in range(self.n):
                for k in range(self.m):
                    pseudo_inverse[i][j] += transpose_matrix[i][k] * self.data[k][j]
        return Matrix(pseudo_inverse)

    @property
    def invert_matrix(self):
        n = self.m
        identity = self.eye(n).data

        for col in range(n):
            for row in range(col, n):
                if self.data[row][col] != 0:
                    self.data[col], self.data[row] = self.data[row], self.data[col]
                    identity[col], identity[row] = identity[row], identity[col]
                    break

            pivot = self.data[col][col]
            for j in range(col, n):
                self.data[col][j] /= pivot
            for j in range(n):
                identity[col][j] /= pivot
            for i in range(n):
                if i == col:
                    continue
                factor = self.data[i][col]
                for j in range(col, n):
                    self.data[i][j] -= factor * self.data[col][j]
                for j in range(n):
                    identity[i][j] -= factor * identity[col][j]
        return Matrix(identity)

    @property
    def show(self):
        for i in range(len(self.data)):
            print(self.data[i], end="\n")

    def generate_random_matrix(self, m, n, start, end):
        if start > end:
            raise ValueError("شروع بازه باید کمتر یا مساوی پایان بازه باشد")

        # تولید ماتریس خالی با ابعاد m × n
        random_matrix = []
        for i in range(m):
            row = []
            for j in range(n):
                # تولید عدد تصادفی در بازه [start, end]
                random_number = random.uniform(start, end)
                row.append(random_number)
            random_matrix.append(row)

        return Matrix(random_matrix)

    def __str__(self):
        for row in self.data:
            for col in row:
                print(col, end = '\t')

    def __add__(self, B):
        C = Matrix.zeros(self, self.m, self.n)
        for i in range(self.m):
            for j in range(self.n):
                C.data[i][j] = self.data[i][j] + B.data[i][j]
        return C

    def __mul__(self, B):
        rows1, cols1, cols2 = self.m, self.n, B.n
        C = [[0 for i in range(cols2)] for j in range(rows1)]
        for i in range(rows1):
            for j in range(cols2):
                for k in range(cols1):
                    C[i][j] += self.data[i][k] * B.data[k][j]

        return Matrix(C)

    def __truediv__(self, B):
        return self * B.invert_matrix

    def __floordiv__(self, B):
        return self.invert_matrix * B

    def __eq__(self, B) -> bool:
        if self.m != B.m:
            return False
        elif self.n != B.n:
            return False
        else:
            for i in range(self.m):
                for j in range(self.n):
                    if self.data[i][j] != B.data[i][j]:
                        return False
        return True

    def __ne__(self, B) -> bool:
        if self.m != B.m:
            return True
        elif self.n != B.n:
            return True
        else:
            for i in range(self.m):
                for j in range(self.n):
                    if self.data[i][j] != B.data[i][j]:
                        return True
        return False

    def __gt__(self, B: "Matrix"):
        if not isinstance(B, Matrix):
            raise ValueError("B is not a Matrix.")
        if (self - B).Matrix:
            return True
        return False

    def sin(self):
        B = Matrix(self.data[:])
        for i in range(self.m):
            for j in range(self.n):
                B.data[i][j] = math.sin(self.data[i][j])
        return B

    def cos(self):
        B = Matrix(self.data[:])
        for i in range(self.m):
            for j in range(self.n):
                B.data[i][j] = math.cos(self.data[i][j])
        return B

    def exp(self):
        B = Matrix(self.data[:])
        for i in range(self.m):
            for j in range(self.n):
                B.data[i][j] = math.exp(self.data[i][j])
        return B

    def tan(self):
        B = Matrix(self.data[:])
        for i in range(self.m):
            for j in range(self.n):
                B.data[i][j] = math.tan(self.data[i][j])
        return B


A = Matrix()
pi = 3.14
A = Matrix([[pi / 2, pi], [pi / 6, pi / 4]])
B = A.sin()
print(B > A)


# A = Matrix([[1, 2, 3], [3, 4, 3], [5, 9, 2]])
# B = A.invert_matrix
# # B.show
# C = B * A
# # C.show
# (A // B).show
# print(A != A)
# A.show

# B = Matrix()
# B = B.eye(5)
# B.show
# B = Matrix([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]])
# (A * B).show

import numpy as np


class Matrix:
    def __init__(self, data=np.ndarray([]), m: int = 0, n: int = 0):
        if m == 0:
            self.m = len(data) if type(data) == list else data.size
        else:
            self.m = m

        if n == 0:
            try:
                self.n = len(data[0])
            except:
                # If data[0] is not np.ndarray or list , then it is a one-dimensional array of numbers.
                self.n = 1
        else:
            self.n = n

        self.data = data if n == 0 else np.reshape(data, (self.m, self.n))

    def __str__(self):
        s = ""
        for row in self.data:
            for col in row:
                s += f"{col}\t"
            s += "\n"
        return s

    @property
    def shape(self):
        return self.data.shape

    @property
    def dtype(self):
        return self.data.dtype

    @staticmethod
    def ones(m, n):
        return Matrix(np.ones((m, n)))

    @staticmethod
    def zeros(m, n):
        return Matrix(np.zeros((m, n)))

    @staticmethod
    def eye(n):
        return Matrix(np.eye(n))

    @staticmethod
    def rand(m, n):
        return Matrix(np.random.rand(m, n))

    @staticmethod
    def randn(m, n):
        return Matrix(np.random.randn(m, n))

    @property
    def transpose(self):
        return Matrix(self.data.T)

    @property
    def inv(self):
        try:
            inv_matrix = np.linalg.inv(self.data)
            return Matrix(inv_matrix)
        except np.linalg.LinAlgError:
            return "ماتریس وارون‌پذیر نیست."

    @property
    def pinv(self):
        pinv_matrix = np.linalg.pinv(self.data)
        return Matrix(pinv_matrix)

    def power(self, exponent):
        return Matrix(np.power(self.data, exponent))

    def elementwise_power(self, exponent):
        return Matrix(self.data**exponent)

    def dot(self, other):
        if self.shape[1] == other.shape[0]:
            dot_product = np.dot(self.data, other.data)
            return Matrix(dot_product)
        else:
            return (
                "تعداد ستون‌های ماتریس اول باید با تعداد سطر‌های ماتریس دوم برابر باشد."
            )

    def __add__(self, other):
        if self.shape == other.shape:
            return Matrix(self.data + other.data)
        else:
            return "سایز ماتریس‌ها یکسان نیست."

    def __sub__(self, other):
        if self.shape == other.shape:
            return Matrix(self.data - other.data)
        else:
            return "سایز ماتریس‌ها یکسان نیست."

    def __mul__(self, other):
        return self.data @ other.data

    def __truediv__(self, other):
        try:
            inv_other = other.inv()
            return self.dot(inv_other)
        except AttributeError:
            return "ماتریس وارون‌پذیر نیست."

    def __floordiv__(self, other):
        try:
            inv_self = self.inv()
            return inv_self.dot(other)
        except AttributeError:
            return "ماتریس وارون‌پذیر نیست."

    def __pow__(self, other):
        return self.power(other)

    def __eq__(self, other):
        return np.array_equal(self.data, other.data)

    def __ne__(self, other):
        return not np.array_equal(self.data, other.data)

    def __lt__(self, other):
        return self.data < other.data

    def __le__(self, other):
        return self.data <= other.data

    def __gt__(self, other):
        return self.data > other.data

    def __ge__(self, other):
        return self.data >= other.data

    def __and__(self, other):
        return np.logical_and(self.data, other.data)

    def __or__(self, other):
        return np.logical_or(self.data, other.data)

    @property
    def __invert__(self):
        return np.logical_not(self.data)

    def __getitem__(self, key, skey=None):
        if isinstance(key, tuple) and len(key) == 2:
            row, col = key
            return self.data[row, col]
        elif isinstance(key, slice):
            return Matrix(self.data[key])
        elif isinstance(key, int):
            return self.data[key]
        else:
            raise IndexError(
                "Invalid: Indexing should be in the form (row, column) or a slice."
            )

    def __setitem__(self, key, value):
        if isinstance(key, tuple) and len(key) == 2:
            row, col = key
            self.data[row, col] = value
        else:
            raise IndexError("Invalid: Indexing should be in the form (row, column).")

    def __call__(self, *args):
        return self.data[args]

    @property
    def min(self):
        return np.min(self.data)

    @property
    def max(self):
        return np.max(self.data)

    @property
    def sin(self):
        return np.sin(self.data)

    @property
    def cos(self):
        return np.cos(self.data)

    @property
    def tan(self):
        return np.tan(self.data)

    @property
    def exp(self):
        return np.exp(self.data)

    def mean(self, axis=None):
        return np.mean(self.data, axis=axis)

    def var(self, axis=None):
        return np.var(self.data, axis=axis)

    def std(self, axis=None):
        return np.std(self.data, axis=axis)

    def sum(self, axis=None):
        return np.sum(self.data, axis=axis)

    def product(self, axis=None):
        return np.prod(self.data, axis=axis)


A = Matrix.randn(3, 3)
print(A, end="\n\n\n")
print(A(slice(1, 3)))

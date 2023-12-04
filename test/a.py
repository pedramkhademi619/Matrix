a = ones(1000)
b = rand(1000, 1000)
c = a + b
d = a * b
a.dot() * b #elementwise
f =  a.dot() / b


g = a ** 2
h = a.dot() ** b

a - a = zeroes(size(a))
a + b - a == b

all(all(b >= 0)) == True->bool
all(any(eye(7, 7))) == True
sum(sum(ones(10, 10)) - 10)




sum(sum(abs(b * inv(b)  - eye(size(b))) )< eps)

and(any(b(b > 1))) == False
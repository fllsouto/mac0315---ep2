#! /bin/octave -qf
m = scanf("%d", 1);
n = scanf("%d", 1);
A = rand(m, n);
b = xb = zero(m, 1);
c = zero(1, n);
global base = zero(m, 1);
global is_base = zero(n, 1);
for j = 1:n
    c(j) = scanf
endfor
for i = 1:m
    for j = 1:n
        A(i, j) = scanf("%f", 1);
    endfor
endfor
for i = 1:m
endfor
for i = 1:m
    b(i) = scanf("%f", 1);
endfor
for i = 1:m
    base(i) = scanf("%d", 1);
    xb(i) = scanf("%f", 1);
    is_base(base(i)) = 1;
endfor
global invB = inv(A(:, base)); 

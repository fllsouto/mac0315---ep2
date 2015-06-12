#! /bin/octave -qf
m = scanf("%d", 1);
n = scanf("%d", 1);
A = zeros(m, n);
b = xb = zeros(m, 1);
c = zeros(1, n);
global base = zeros(m, 1);
global is_base = zeros(n, 1);
for j = 1:n
    c(j) = scanf("%f", 1);
endfor
for i = 1:m
    for j = 1:n
        A(i, j) = scanf("%f", 1);
    endfor
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
disp(base);
disp(invB);

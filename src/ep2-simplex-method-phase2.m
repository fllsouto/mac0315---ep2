#! /bin/octave -qf

function c_bar = calculate_c_bar (c, j, p_line, A)
  c_bar = c(j) - p_line*A(:, j);
endfunction

function new_B_inv = calculate_new_B_inv(imin, B_inv, u, m)
  for i=1:m
    if(i != imin)
      B_inv(i, :) = B_inv(i , :) - (u(i)/u(imin))*B_inv(imin , :);
    endif
  endfor
  B_inv(imin, :) = B_inv(imin, :)/u(imin);
  new_B_inv = B_inv;
endfunction

function [theta, imin] = calculate_theta(u, m, x, base)
  theta = Inf;
  imin = 0;
  for i=1:m
    if(u(i) > -10^-10)
      if(x(i)/u(i) < theta)
        theta = x(i)/u(i);
        imin = i;
      endif
    endif
  endfor
endfunction

function [ind v] = simplex(A, b, c, m, n, x)
  
  base = zeros(1, m);
  l = 1;
  for k = 1:n
    if(x(k) > 0)
      base(l) = k;
      l = l + 1;
    endif
  endfor

  xb = x(base);
  B_inv = inv(A(:, base));
  count = 0;
  while (true)
    printf("\n-------------------------------------------------------------------------\n");
    printf("Iterando ; %d\n", count);
    for k=1:m
      printf("%d %f\n",base(k),xb(k));
    endfor
    printf("\nValor da funcao objetivo : %f\n\n", c(base)*xb);

    p_line = c(base)*B_inv;
    is_base = zeros(1:n);
    is_base(base) = 1;

    indx_j = 0;
    achou = 0;
    c_bar = zeros(n, 1);
    printf("\nCustos reduzidos\n")
    for j = 1:n
      if(is_base(j) == 0)
        c_bar(j) = calculate_c_bar(c , j, p_line, A);
        printf("%d : %f\n",j, c_bar(j));
        if (c_bar(j) < -10^-10 && !achou)
            indx_j = j;
            achou = 1;
        endif
      endif
    endfor
    
    if (indx_j == 0)
      ind = 0;
      v = zeros(n, 1);
      v(base) = xb;
      break;
    endif

    printf("\nEntra na base : %d\n", indx_j);
    u = B_inv*A(:, indx_j);
    
    printf("\nDirecao\n");
    for k=1:m
      printf("%d %f\n",base(k), u(k));
    endfor

    [theta, imin] = calculate_theta(u, m, xb, base);
    
    if(theta == Inf)
      ind = -1;
      v = zeros(n , 1);
      v(base) = u;
      break;
    else
      old_base_imin = base(imin);
      base(imin) = indx_j;
      xb = xb - theta*u;
      xb(imin) = theta;
      B_inv = calculate_new_B_inv(imin, B_inv, u, m);
    endif

    printf("\nTheta*\n%f\n",theta);
    printf("\nSai da base: %d\n",old_base_imin);
    count = count +1;
  endwhile

endfunction

#m = 2; 
#n = 4;
#A = [2,1,1,0;1,2,0,1];
#b = [4;3];
#base = [1,4];
#c = [-1, -1, 0, 0];
#m = 3;
#n = 6;
#A = [ 1, 2, 2, 1, 0, 0; 
#        2, 1, 2, 0, 1, 0;
#        2, 2, 1, 0, 0, 1;];
#b = [20; 20; 20];
#base = [4, 5, 6];
#c = [-10, -12, -12, 0, 0, 0];
A = [4.73, -4.73,  10.3,  -10.3, -1,  0, 0, 0, 0;
    -0.83,  0.83,  1.05,  -1.05,  0, -1, 0, 0, 0;
    -1.56,  1.56, -0.02,   0.02,  0,  0, 1, 0, 0;
    -2.11,  2.11, 10.22, -10.22,  0,  0, 0, 1, 0;
     8.32, -8.32, 12.38, -12.38,  0,  0, 0, 0, 1];
b = [4.88; -3.61; 1.78; 77.52; 120.26];
c = [1.5; -1.5; -0.25; 0.25; 0; 0; 0; 0; 0]';
n = 9;
m = 5;
base = [1, 4, 7, 8, 9];
x = zeros(n , 1);
x(base) = A(:,base) \ b;
[ind v] = simplex(A, b, c, m, n, x);

if(ind == 0)
  best_c = c*v;
  printf("\nSolucao otima encontrada com custo : %f\n",best_c);
  for k=1:n
    printf("%d %f\n", k, v(k));
  endfor
else
  printf("\nCusto otimo menos infinito\n");
  printf("\nDirecao viavel : \n");
  for i=1:n
    printf("%d %f\n",i,v(i));
  endfor
endif


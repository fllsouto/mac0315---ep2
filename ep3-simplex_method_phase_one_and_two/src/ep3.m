#! /bin/octave -qf
#Universidade de Sao Paulo - USP
#Instituto de Matematica e Estatistica
#
#MAC0315 - Programacao linear
#Exercicio Programa III - Implementacao do Metodo Simplex
#Prof : Ernesto G. Birgin
#
#Fellipe Souto Sampaio - 7990422
#Lucas Romão Silva - 8536214  

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
    if(u(i) > 10^-10)
      if(x(i)/u(i) < theta)
        theta = x(i)/u(i);
        imin = i;
      endif
    endif
  endfor
endfunction

function [ind, v, last_base, last_B_inverse] = simplex_phase_two(A, b, c, m, n, x, B_inv, base)
  
  xb = x(base);
  count = 0;
  while (true)
    printf("\n---------------------------------------------------\n");
    printf("Iterando : %d\n", count);
    for k=1:m
      printf("%d %f\n",base(k),xb(k));
    endfor
    printf("\nValor da funcao objetivo : %f\n\n", c(base)*xb);

    p_line = c(base)*B_inv;
    is_base = zeros(1, n);
    is_base(base) = 1;
    indx_j = 0;
    achou = 0;
    c_bar = zeros(n, 1);
    printf("\nCustos reduzidos\n");
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
      last_base = base;
      last_B_inverse = B_inv;
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
      last_base = base;
      last_B_inverse = B_inv;
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


function [base, new_b, new_B_inv, new_A, new_m] = prepare_to_phase_two(A, b, B_inv, m, n, last_base)
  new_B_inv = B_inv;
  new_b = b;
  new_m = m;
  result = 0;
  new_A = A;
  base = last_base;
  for i=1: m
    if(last_base(i) > n)
      for j=1:n
        result = new_B_inv(i,:)*new_A(:,j);
        if(result != 0)
          base(i) = j;
          new_B_inv = calculate_new_B_inv(i, new_B_inv, new_B_inv*new_A(:,j), m);
          break;
        endif
      endfor

      if(j == n && result == 0)
        new_m = new_m - 1;
        new_A(i,:) = [];
        base(i) = [];
        new_B_inv(i,:) = [];
        new_B_inv(:,i) = [];
        new_b(i) = [];
      endif
    endif
  endfor
endfunction

function [ind, x, d] = simplex(A, b, c, m, n)
  base = (zeros(1,m));
  for j = 1:m+n
    if(j <= m && b(j) < 0)
      b(j)    = b(j)*(-1);
      A(j,:)  = A(j,:)*(-1);
    endif

    if(j > n)
      base(j - n) = j;
    endif
  endfor

  A_line = [A, eye(m)];
  x_line = [zeros(n, 1); b];
  c_line = [zeros(1, n), ones(1, m)];
  printf("\n---------------------------------------------------\n");
  printf("Fase 1\n");
  printf("\n---------------------------------------------------\n");
  [ind, v, last_base, B_inv] = simplex_phase_two(A_line, b, c_line, m, m + n, x_line, eye(m), base);
  if(c_line*v != 0 )
    ind = 1;
    x = v;
    d = 0;
  else
    [base, new_b, new_B_inv, new_A, new_m] = prepare_to_phase_two(A, b, B_inv, m, n, last_base);
    printf("\n---------------------------------------------------\n");
    printf("Fase 2\n");
    printf("\n---------------------------------------------------\n");
    [ind, v, last_base, B_inv] = simplex_phase_two(new_A, new_b, c, new_m, n, v, new_B_inv, base);
    if(ind == 0)
      x = v;
      d = 0;
    else
      x = 0;
      d = v;
    endif
  endif
  

endfunction



##########################################################
#Defina aqui as variaveis (A, b, c, m, n)#################
##########################################################

A = [4.73, -4.73,  10.3,  -10.3, -1,  0, 0, 0, 0;
    -0.83,  0.83,  1.05,  -1.05,  0, -1, 0, 0, 0;
    -1.56,  1.56, -0.02,   0.02,  0,  0, 1, 0, 0;
    -2.11,  2.11, 10.22, -10.22,  0,  0, 0, 1, 0;
     8.32, -8.32, 12.38, -12.38,  0,  0, 0, 0, 1];
b = [4.88; -3.61; 1.78; 77.52; 120.26];
c = [1.5; -1.5; -0.25; 0.25; 0; 0; 0; 0; 0]';
n = 9;
m = 5;
[ind, x, d] = simplex(A, b, c, m, n);

switch(ind)
  case -1
    printf("\nCusto otimo menos infinito\n");
    printf("\nDirecao viavel : \n");
    for i=1:n
      printf("%d %f\n",i,d(i));
    endfor
  case 0
    best_c = c*x;
    printf("\nSolucao otima encontrada com custo : %f\n",best_c);
    for k=1:n
      printf("%d %f\n", k, x(k));
    endfor
  case 1
    printf("\nProblema inviável\n");
endswitch

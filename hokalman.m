function [A, B, C, D] = hokalman(Gbl)
m = size(Gbl,2);

D = squeeze(Gbl(1,:, :));
%constroi a matriz de blocos de Hankel
H = hankel(Gbl);

%performa o SVD na matriz de Hankel
[U, S, V] = svd(H);
%reduz as matrizes U, S, V para a ordem escolhida do sistema
[U1, S1, V1]= modelReduction(U, S, V);

%dividir blocos
Obs = U1*S1;
Ctrl = S1*V1';

k3 = size(H,1);

A = Obs(1 : k3 - m, :) \ Obs(m + 1 : k3, :);

B = Ctrl(:,1:m);

C = Obs(1 : m, :);

end

function [U1, S1, V1] = modelReduction(U, S, V)
  d = diag(S);
  %0.0008 e' o menor valor que nao e' zero
  nx = find(d>0.0008); 
  nx = nx(end);
 
  % Reduz para nx
  U1 = U(:, 1:nx);
  S1 = sqrtm(S(1 : nx, 1 : nx));
  V1 = V(:, 1:nx);
end


function H = hankel(Gbl)
  m = size(Gbl,2);

  nIR = 448;
  kl = 1;
  kc = 1;
  k3 = (nIR - 1)/2;
  for i = 1:k3
    for j = 1:k3
        H(kl: kl + m - 1, kc : kc + m - 1) = squeeze(Gbl(i+j,:,:));
        kc = kc + m;
    end
    kl = kl + m;
    kc = 1;
  end
end
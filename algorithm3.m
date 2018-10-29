function [Gbl, G] = algorithm3(U, Y)
%m e' o numero de entradas e saidas
 m = size(U,1);
[Ubl, Ybl] = mountMatrix(U, Y);

%encontra a pseudoinversa
U_dagger = pinv(transpose(Ubl));

%encontra as respostas ao impulso
G = transpose(U_dagger*transpose(Ybl));

%apos corte dos blocos - o numero nIR foi encontrado experimentalmente
nIR = 448;
G = G(:, 1:nIR*m);
k = 1;

for i = 1:nIR
    Gbl(i, 1:m, 1:m) = G(1:m, k:k + m - 1); 
    k = k + m;
end

end

function [Ubl, Ybl] = mountMatrix(U, Y)
    N = size(U,2); %numero experimentos
    m = size(U,1); %numero entradas
    
    k = floor(N/2);
    kk = k;
    kf = k;
    %loop que vai construir a matriz de blocos Ubl
    %anda de m em m blocos
    for i = 1:m:m*k
        for j = 1:(N-kf+1)
            Ubl(i:i+m-1,j) = U(1:m, kk);
            kk = kk + 1;
        end
        k = k - 1;
        kk = k;
    end
    
    %loop que constroi a matriz Ybl
    for i = kf:N
        for j = 1:(N-kf+1)
            Ybl(1:m, j) = Y(:,i);
        end
    end
end
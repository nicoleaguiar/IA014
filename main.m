%%
clear all;
load('MS_Harm_h3_N1568_RMS70_P2P350.mat')

U = u_m';
Y = y_m';
%%
%encontra as respostas ao impulso do sistema
[Gbl, G] = algorithm3(U, Y);
%%
%encontra A, B, C, D
[A, B, C, D] = hokalman(Gbl);

%%

h = ss(A, B, C, D, .1);
impulse(h)
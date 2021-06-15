E = 210e9; I = 0.0833; q = 0; P = 20000; M = 0; density = 7850;

le = 10;

%Element stiffness matrix and load vector
k1 = Beam_stiff_matrix(E,I,[0,10]);
f = Load_vec_matrix(q,[0,10]);

f(3) = f(3) + P;
f(4) = f(4) + M;

%Solving equations for deflection ( Essential Boundary conditions : w1 = 0
%and (dw/dx)1 = 0)
W = zeros(4,1);

Kreduced = k1([3,4],[3,4]);
Freduced = f([3,4]);
wreduced = inv(Kreduced)*Freduced;

W(3:4) = wreduced(1:2);
Reaction = k1*W;

%Natural Frequency
Mass_mat = Beam_Mass_matrix(1,density,10);


%Plotting the solution
w_final = zeros(length([-1,0.1,1]),1);
w_actual = zeros(length([-1,0.1,1]),1);
i = 1;
for e = -1 : 0.1 : 1
    N1 = 0.25*(1-e)*(1-e)*(2+e);
    N2 = 0.25*(1-e)*(1-e)*(1+e);
    N3 = 0.25*(1+e)*(1+e)*(2-e);
    N4 = 0.25*(1+e)*(1+e)*(e-1);
    w_final(i) = N1*W(1) + le*0.5*N2*W(2) + N3*W(3) + le*0.5*N4*W(4);
    X_fin = (1+e)*0.5*10;
    w_actual(i) = (P*X_fin*X_fin/(6*E*I))*(3*le - X_fin);
    i = i+1;
end

x = linspace(0,10,21);
%plot(x,-1*w_final,'blue');
plot(x,-1*w_actual,'red');
title('Analytical Solution'); xlabel('x(m)');ylabel('Deflection(m)');




% Solución por el método de los momentos de la distribución de carga 
% sobre una placa cuadrada de lado L, conectada a un potencial V
% Ver libro de Gibson Ej:3.2.
% 
% Solution using MoM of the charge distribution on a square metal plate of
% side length L, connected to a potential V, see Gibson's book.
%
% Germán Augusto Ramírez Arroyave
% CMUN - 2017

%% Physical parameters of the plate and discretization:
L = 1; % side length in meters
N = 18; % side partition (there will be N^2 cells)
DeltaL = L/N; 
DeltaS = DeltaL^2; 
a = DeltaL/2; % With the proposed discretization DeltaL = 2a, according to the notation used in the text

%% Define the boundary condition
eps = 8.85e-12;
B = ones(N*N,1); % Potential on the conductor's surface is set to 1V.

%% Discretization 
[xx, yy] = ndgrid(linspace(DeltaL,L,N),linspace(DeltaL,L,N));
xobs = xx - a;
yobs = yy - a;

%% Compute Z matrix
% Self terms (diagonal terms) are given by:
Z = 2*a*log(1+sqrt(2))*eye(N*N)/(pi*eps);
% The remaining terms are calculated by the integral:
for m = 1:N*N
    for n = 1:N*N
        if m~=n
           Z(m,n) = 1/sqrt( (xobs(m)-xobs(n))^2 + (yobs(m)-yobs(n))^2 );
        end
    end
end
Z = Z*DeltaS;

%% Invert the system to find the required charge distribution:
% cond(Z) % check matrix condition number (this problem is well posed)
rho = Z\B;
plot(rho); % This vector must be converted to a 2D array to be meaningful

%% Array the charge distribution in a matrix form in order to plot on the surface
rhom = zeros(N,N);
for k = 1:N
    rhom(k,:) = rho((k-1)*N+1:k*N);
end
figure;
surf(xx-a,yy-a,rhom);
title('Distribución de carga para un plano conductor');
xlabel('Longitud (m)');
ylabel('Longitud (m)');
zlabel('Densidad de carga C/m');
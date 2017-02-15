% Programa para calcular la distribución de carga en un hilo de alambre,
% empleando el método de los momentos, según la descripción del libro de
% Gibson.
%
% Simple program to calculate charge distribution on a wire using a MoM 
% implementation according to Gibson's book.
% The equations system: Za = b, is to be solved, being 'bi = 4pieps', Z the
% coefficients matrix from the integral of the Green's función and 'a' the
% vector with the charge distribution to be found
% 
% Germán Augusto Ramírez Arroyave
% CMUN - 2017

%% Physical parameters involved
a = 0.001; % Radius of the wire in meters
long = 1; % Length of the wire in meters 
%% Parameters for the N points longitudinal discretization of the wire:
N = 100; % Number of points for the discretization 
deltaX = long/N; % Discretiation size
Xm = deltaX/2;
X = linspace(deltaX,long,N); 
Xobs = X - deltaX/2; % Observation points are taken in the segments' centroids:

%% Define boundary condition
eps = 8.85e-12;
B = 4*pi*eps*ones(N,1); % potential on the conductor's surface is set to 1V.

%% Calculation of Z matrix
tic
Z = zeros(N,N); % Allocate memory for the coefficients matrix
for m = 1:N
    for n = 1:N
        Xb = n*deltaX; % ending point of the segment 
        Xa = (n-1)*deltaX; % initial point of the segment 
        Z(m,n) = log( ((Xb-Xm) + sqrt((Xb-Xm)^2 - a^2)) / ((Xa-Xm) + sqrt((Xa-Xm)^2 - a^2)) );
    end
    Xm = Xm + deltaX;
end
toc

%% Improved version
% The first approach is really slow. It is convenient to take advantage of 
% Matlab's vectorization to fit each row of the matrix using a single cycle
% Matrix symmetry could be used to optimize computation time even further
tic
M = zeros(N,N);
for cont=1:N
	M(cont,:) = log( ((X(cont)-Xobs) + sqrt((X(cont)-Xobs).^2 - a^2)) ./ ((X(cont)-deltaX-Xobs) + sqrt((X(cont)-deltaX-Xobs).^2 - a^2)) );   
end
toc

%% Invert the system to find the required charge distribution:
% cond(Z) % check matrix condition number (this problem is well posed)
rho1 = Z\B; % X = A\B solves the system A*X = B 
rho2 = M\B; 

%% Plot charge distribution:
plot(Xobs,abs(rho1),'b', Xobs,abs(rho2),'r'); legend('cycles','vectorized');
title('Distribución de carga para un alambre lineal');
xlabel('Longitud (m)');
ylabel('Densidad de carga C/m');
% Programa que emplea el M�todo de los momentos para calcular la 
% capacitancia de un arreglo arbitrario de objetos conductores a partir 
% de un archivo .dat que describe la geometr�a.
% Se emplean las expresiones enunciadas en el art�culo de Wilton Rao y
% Glisson para calcular las integralas de las funciones de potencial.

% =======================================================================
clc; clear all; close all; % ADVERTENCIA... Guardar antes de ejecutar esto...
% =======================================================================
tic


% Se carga la geometr�a a analizar:

% ------------------ Placas paralelas (lxwxh) = 20x20x1 cm --------------
 geom = readDatMesh('Pl_par_Finest_20x20x1.dat');   
% geom = readDatMesh('Pl_par_Fine_20x20x1.dat'); 
% geom = readDatMesh('Pl_par_Good_20x20x1.dat'); 
% geom = readDatMesh('Pl_par_Moderate_20x20x1.dat'); 
% geom = readDatMesh('Pl_par_Coarsest_20x20x1.dat'); 
% ------------------ Placas paralelas diferentes ------------------------
% geom = readDatMesh('DosPlacas_Coarse.dat');
% geom = readDatMesh('Placas_par_20x20x5_Coarse.dat');
% geom = readDatMesh('Placas_par_20x20x5_Fine.dat');
% ------------------ Esferas diferents ------------------
% geom = readDatMesh('Esferas_15x15x5_Coarse.dat');
% geom = readDatMesh('Esferas_10_10_1_Moderate.dat'); 
% ------------------ Cilindros paralelos (lxrxd) 50x5x2 cm -------------
% geom = readDatMesh('Cilindros_par_5_1_Fine.dat'); 
% geom = readDatMesh('Cilindros_par_5_1_Coarse.dat'); 

% Se buscan los centroides y los vectores normales a cada uno de los
% tri�ngulos, y de paso se muestra la geometr�a
[centros normales areas] = Geomet(geom);

% N�mero total de celdas de la geometr�a
Ncells = size(geom.cells);
% Se define que lo que est� sobre el plano xy es un objeto y lo que 
% est� por debajo es otro objeto
Nobj1 = sum( centros(:,3) > 0 );
Nobj2 = sum( centros(:,3) < 0 );

% Se reserva espacio en memoria para la matriz de impedancia
Z = zeros(Ncells(1),Ncells(1));

% Se forma la matriz llamando la rutina correspondiente
for i = 1:Ncells(1)
    % Arreglo con las coordenadas de los v�rtices del tri�ngulo
    v(1,:) = geom.verts(geom.cells(i,3),:);
    v(2,:) = geom.verts(geom.cells(i,2),:);
    v(3,:) = geom.verts(geom.cells(i,4),:);

    Z(:,i) = MatrizZ(centros(:,:), normales(i,:), v);
end

% Se definen las condiciones de frontera
eps0 = 8.85e-12;
V1 = 1;
V2 = 1;

% B1 = 4*pi*eps0*[V1*ones(1,Nobj1) zeros(1,Nobj2)]';
% B2 = 4*pi*eps0*[zeros(1,Nobj1) V2*ones(1,Nobj2)]';

%exciting top plate
B1 = zeros(Nobj1+Nobj2,1);
B1(centros(:,3)>0)= 4*pi*eps0*V1;

%exciting bottom plate
B2 = zeros(Nobj1+Nobj2,1);
B2(centros(:,3)<0)= 4*pi*eps0*V2;


% Se invierte para obtener la distribuci�n superficial de carga en cada
% tri�ngulo
rho1 = Z\B1;
rho2 = Z\B2;

% Se multiplica por el �rea para hallar la carga total de cada tri�ngulo
qi1 = rho1.*areas;
qi2 = rho2.*areas;

figure;
plot(rho1); hold on;
plot(rho2, 'r');
title('Distribuci�n de carga en las celdas');
xlabel('N�mero de celda ');
ylabel('Sigma (C/m^2)');
legend('Sup1 a 1V, Sup2 a 0V', 'Sup1 a 0V, Sup2 a 1V');

figure;
trisurf(geom.cells(:,2:4),geom.verts(:,1), geom.verts(:,2), geom.verts(:,3), rho1');
axis equal;
title('Distribuci�n de carga en las celdas (Sup1 a 1V, Sup2 a 0V)');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
colorbar;
figure;
trisurf(geom.cells(:,2:4),geom.verts(:,1), geom.verts(:,2), geom.verts(:,3), rho2');
axis equal;
title('Distribuci�n de carga en las celdas (Sup1 a 0V, Sup2 a 1V)');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
colorbar;

% Carga en cada superficie para la condici�n Sup1 a 1V, Sup2 a 0V
q1 = sum(qi1(1:Nobj1));
q2 = sum(qi1(Nobj1+1:size(qi1,1)));
% Carga en cada superficie para la condici�n Sup1 a 0V, Sup2 a 1V
q3 = sum(qi2(1:Nobj1));
q4 = sum(qi2(Nobj1+1:size(qi2,1)));

V = [1 0; 0 1];
Q = [q1 q2; q3 q4];

C = Q/V;
cap = (C(1,1)*C(2,2)-C(1,2)*C(2,1)) / (C(1,1)+C(2,2)+C(1,2)+C(2,1));
toc
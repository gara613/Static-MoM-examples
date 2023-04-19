% Función encargada de cargar y mostrar la geometría del problema, además
% de calcular parámetros relevantes relacionados con la misma.
% Salidas:
%       centro [Mx3]: coordenadas de los centroídes de cada triángulo
%       nuni [Mx3]: vector unitario normal a cada triángulo
%       Area [Mx1]: área de cada triángulo
% Entradas:
%       geom: estructura de vértices y coordenadas tal como es devuelta por
%       la función readDatMesh()

function [centro nuni Area] = Geomet(geom)

Nverts = size(geom.verts);
Ncells = size(geom.cells);

% Reservar espacio para los vectores que se entregan con las coordenadas
% de los centros y el vector normal a cada triángulo
centro = zeros(Ncells(1),3);
nuni = zeros(Ncells(1),3);
Area = zeros(Ncells(1),1);

for i = 1:Ncells(1)
    % Xcomp es un arreglo con las componentes X de los vértices del
    % triángulo
    
    xcomp = geom.verts(geom.cells(i,2:4),1); % Componentes en x
    ycomp = geom.verts(geom.cells(i,2:4),2); % Componentes en y
    zcomp = geom.verts(geom.cells(i,2:4),3); % Componentes en z
    
    % Determina el centroíde de cada triángulo
    centro(i,1) = sum(xcomp(1:3)/3); % En x
    centro(i,2) = sum(ycomp(1:3)/3); % En y
    centro(i,3) = sum(zcomp(1:3)/3); % En z
    
    % Pi es la coordenada (x,y,z) del vértice i del triángulo
    p1 = geom.verts(geom.cells(i,2),:);
    p2 = geom.verts(geom.cells(i,3),:);
    p3 = geom.verts(geom.cells(i,4),:);
       
    % Se determina y normaliza el vector normal al triángulo dado
    n = cross(p2-p1, p2-p3);
    Area(i) = 0.5*norm(cross(p1-p2, p1-p3));
    nuni(i,:) = n/norm(n);
    
    % Se refiere el vector normal al centro de cada triángulo (Por si se quiere visualizar)
    Nuni = [centro(i,1), centro(i,2), centro(i,3); centro(i,1)+10*Area(i)*nuni(i,1),centro(i,2)+10*Area(i)*nuni(i,2),centro(i,3)+10*Area(i)*nuni(i,3)];
    plot3(Nuni(:,1),Nuni(:,2),Nuni(:,3),'r');
    hold on;
end  

% Agrega el punto del centroíde al gráfico
hold on;
plot3(centro(:,1),centro(:,2),centro(:,3),'rx');

trisurf(geom.cells(:,2:4),geom.verts(:,1), geom.verts(:,2), geom.verts(:,3));
hold on;

axis equal;
title('Geometría del problema');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
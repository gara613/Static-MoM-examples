% Cálculo de las integrales de potencial para distribuciones superficiales 
% de carga en polígonos.
% zmn = MatrizZ(r, nuni, vert)
% Salida:
%       zmn [Mx1]= Elemento [mn] de la matriz de "impedancia" calculado por
%       medio de las integrales reportadas en el artículo de RWG.
% Entradas: 
%       r [Mx3]: coordenadas del punto de observación.
%       nuni [Mx1]: vector unitario normal al plano del polígono.
%       vert [3x3]: vértices del polígono.

function z = MatrizZ(r, nuni, vert)

z = zeros(size(r,1),1); % Inicializar el valor de retorno
vert(4,:)= vert(1,:); % Vértice adicional para cerrar el triángulo

% Para cada uno de los segmentos del polígono
for i = 1:3 
    
    rplus = vert(i+1,:); % Vector que describe el extremo superior del segmento
    rminus = vert(i,:); % Vector que describe el extremo inferior del segmento
    % Otra opción puede ser: rplus = repmat(vert(i+1,:), size(r,1),1);
    
    % Altura del punto de observación sobre el plano de la fuente
    % d = dot(nuni, r-rplus);
     r1 = bsxfun(@minus, r, rplus);
     d = r1*nuni';
     
    % Proyección del vector del punto de prueba sobre el plano de la fuente
    % rho = r - nuni*dot(nuni, r);
     aux = (r*nuni')*nuni;
     rho = r - aux;
                
    % Proyecciones del vector de fuente sobre el plano de la fuente
    rhoplus = rplus - nuni * dot(nuni, rplus); 
    rhominus = rminus - nuni * dot(nuni, rminus);

    % Vectores y magnitudes definidas en el artículo de RWG
    % luni = (rhoplus - rhominus) / norm(rhoplus - rhominus); %
    % Inicialmente, luego lo cambian por:
    luni = (rplus - rminus) / norm(rplus - rminus); 
    Uuni = cross(luni, nuni);
    
    % lplus = dot(rhoplus - rho, luni); 
    % lminus = dot(rhominus - rho, luni);
    aux = bsxfun(@minus,rhominus,rho);
    lminus = aux*luni'; 
    aux = bsxfun(@minus,rhoplus,rho);
    lplus = aux*luni';

    % Po = abs( dot (rhoplus - rho, Uuni));
    Po = abs(aux*Uuni');
    % Pouni = ((rhoplus-rho) - lplus*luni) / Po;
    Pouni = bsxfun(@rdivide, (aux - lplus*luni), Po); 
     
    Pplus = sqrt(Po.^2 + lplus.^2);
    Pminus = sqrt(Po.^2 + lminus.^2);
    
    Ro = sqrt(Po.^2 + d.^2);

    Rplus = sqrt(Pplus.^2 + d.^2);
    Rminus = sqrt(Pminus.^2 + d.^2);

    % z = z + dot(Pouni,Uuni) * (Po * log((Rplus + lplus)/(Rminus + lminus)) - abs(d) * (atan (Po*lplus/(Ro^2 + abs(d)*Rplus)) - atan(Po*lminus/(Ro^2 + abs(d)*Rminus)) ) );
    z = z + (Pouni*Uuni') .* (Po .* log((Rplus + lplus)./(Rminus + lminus))...
        - abs(d) .* (atan (Po.*lplus./(Ro.^2 + abs(d).*Rplus)) - atan(Po.*lminus./(Ro.^2 + abs(d).*Rminus)) ) );
end
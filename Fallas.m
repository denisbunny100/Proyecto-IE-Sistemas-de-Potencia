function [Vfalla, Ifalla, Scorto] = Fallas(y_barra0, y_barra1, conexiones, datos_prefalla, Recsub0, Recsub1)
%Recsub0 -> barra rec0
%conexiones -> barra1 barra2 conP conS isTri Trirepdel
y_sec0 = y_barra0;
y_sec1 = y_barra1;
if ~isempty(conexiones)
    for i = 1:length(conexiones(:,1))
        Barra_delta = conexiones(i,conexiones(i,3:4) == 2);
        Barra_y_solid = conexiones(i,conexiones(i,3:4) == 1);
        if ~isempty(Barra_delta)
            if isempty(Barra_y_solid) 
                yx = y_sec0(Barra_delta(0),Barra_delta(1));
                y_sec0(Barra_delta(0),Barra_delta(1)) = y_sec0(Barra_delta(0),Barra_delta(1)) - yx;
                y_sec0(Barra_delta(1),Barra_delta(0)) = y_sec0(Barra_delta(1),Barra_delta(0)) - yx;
                if conexiones(i,6) ~= 0
                    y_sec0(Barra_delta(0),Barra_delta(0)) = y_sec0(Barra_delta(0),Barra_delta(0)) + yx;
                    y_sec0(Barra_delta(1),Barra_delta(1)) = y_sec0(Barra_delta(1),Barra_delta(1)) + yx;
                end
            else
                yx = y_sec0(Barra_y_solid,Barra_delta);
                if conexiones(i,5) == 1
                    y_sec0(Barra_delta,Barra_y_solid) = yx;
                end
                y_sec0(Barra_y_solid,Barra_delta) = y_sec0(Barra_y_solid,Barra_delta) - yx;
                y_sec0(Barra_delta,Barra_y_solid) = y_sec0(Barra_delta,Barra_y_solid) - yx;
                if conexiones(i,6) ~= 0
                    y_sec0(Barra_delta,Barra_delta) = y_sec0(Barra_delta,Barra_delta) + yx;
                end
            end
        end
    end
end
y_sec0(Recsub0(:,1),Recsub0(:,1)) = diag(diag(y_sec0(Recsub0(:,1),Recsub0(:,1))) + 1./Recsub0(:,2));
y_sec1(Recsub1(:,1),Recsub1(:,1)) = diag(diag(y_sec1(Recsub1(:,1),Recsub1(:,1))) + 1./Recsub1(:,2));
zsec0 = inv(y_sec0);
zsec1 = inv(y_sec1);
%Corrientes de falla
%datos prefalla -> Vf nfalla tipofalla Zf
%Tipos falla
%ASIMÉTRICAS
%1: Linea Tierra
%2: Linea Linea
%3: Doble linea a Tierra
% 4: Trifásicas
a = 1*exp(1j*deg2rad(120));
A = [1 1 1; 1 a^2 a; 1 a a^2];
Vf = datos_prefalla(1);
Zf = datos_prefalla(4);
znn1 = zsec1(datos_prefalla(2),datos_prefalla(2));
znn0 = zsec0(datos_prefalla(2),datos_prefalla(2));
switch datos_prefalla(3)
    case 1
        I = Vf/(1j*(znn0 + znn1 + 3*Zf));
        I = [I; I; I];
    case 2
        I = [0;0;0];
        I(1:2:end) = Vf/(1j*(znn1+Zf));
    case 3
        I1 = Vf/(1j*znn1); I2 = - I1; I0 = 0;
        I = [I1; I2;I0];
    case 4
        I=Vf/(inv(y_barra1(datos_prefalla(2),datos_prefalla(2))));
        I = [0; I; 0];
end
Ifalla = A*I;
zsecs = zeros(3); zsecs(1,1) = znn0; zsecs(2,2) = znn1;
Vaux = [0; Vf; 0];
Vfalla = Vaux - 1j*zsecs*I;
Scorto = conj(Vfalla)'*conj(Ifalla);
%Corrientes del generador 
Ig = [];
for i = 1:length(Recsub0)
    zs = zeros(3); zs(1,1) = zsec0(datos_prefalla(2),Recsub0(i,1)); zs(2,2) = zsec1(datos_prefalla(2),Recsub1(i,1));
    Vg = Vaux - 1j*zs*I;
    Ig(i) = (1 - Vg)/Recsub1(i,2);
end

end


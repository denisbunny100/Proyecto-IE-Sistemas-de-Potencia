function [PS_out, flujos, perdidas] = Flujos(V, delta, Y_barra, datos_potencia, datos_linea, datos_trafo, r2p, r3p, Sb)
    v = V.*exp(1j*delta);
    Sbarra = conj(v).*sum(conj(v').*Y_barra,2);
    Sg = conj(Sbarra) + (datos_potencia(:,5)+1j*datos_potencia(:,6))/Sb;
    I = zeros(length(datos_linea(:,1)),2);
    S = zeros(length(datos_linea(:,1)),2);
    Sp = zeros(length(datos_linea(:,1)),1);
    for i = 1:length(datos_linea(:,1))
        I(i,1) = (v(datos_linea(i,1)) - v(datos_linea(i,2)))./sum(datos_linea(i,3:4),2) + v(datos_linea(i,1)).*datos_linea(i,5); 
        I(i,2) = (v(datos_linea(i,2)) - v(datos_linea(i,1)))./sum(datos_linea(i,3:4),2) + v(datos_linea(i,2)).*datos_linea(i,5);
        S(i,1) = v(datos_linea(i,1)).*conj(I(i,1));
        S(i,2) = v(datos_linea(i,2)).*conj(I(i,2));
        Sp(i) = S(i,2)-S(i,1);
    end
    count = length(datos_linea(:,1))+1;
    c2 = 1; c3 = 1;
    S3 = [];
    if ~isempty(datos_trafo)
        for i = 1:length(datos_trafo(:,1))
            if datos_trafo(i,3)~=0   
                %I1 I2 I3
                I3p = r3p(:,:,c3)*[v(datos_trafo(i,1));v(datos_trafo(i,2));v(datos_trafo(i,3))];
                I3p(2:end) = -I3p(2:end);
                S3(:,c3) = v(datos_trafo(i,1)).*conj(I3p);
                c3 = c3 + 1;
            else
                I2p = r2p(:,:,c2)*[v(datos_trafo(i,1));v(datos_trafo(i,2))];
                S(count,1) = v(datos_trafo(i,1))*conj(I2p(1));
                S(count,2) = v(datos_trafo(i,2))*conj(-I2p(2));
                Sp(count) = S(count,2)-S(count,1);
                count = count + 1;
                c2 = c2 + 1;
            end
        end
    end
    LS = [datos_linea(:,1:2); datos_trafo(datos_trafo(:,3)==0,1:2)];
    S = reshape(S, 2*length(LS(:,1)), []);
    Barras = datos_potencia(:,1);
    del_grados = rad2deg(delta);
    Pg = real(Sg); Qg = imag(Sg); Pd = datos_potencia(:,5)/Sb; Qd = datos_potencia(:,6)/Sb; Pp = real(Sp); Qp = imag(Sp);  
    PS_out = table(Barras, V, del_grados, Pg, Qg, Pd, Qd);
    L1 = LS(:,1);L2 = LS(:,2);
    perdidas = table(L1, L2, Pp, Qp);
    f1 = [LS(:,1);LS(:,2)];f2 = [LS(:,2);LS(:,1)]; 
    P = real(S); Q = imag(S); Sabs = abs(S); Sang = rad2deg(angle(S));
    flujos = table(f1, f2, P, Q, Sabs, Sang);
end


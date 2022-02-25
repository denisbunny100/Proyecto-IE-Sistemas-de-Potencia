function [S] = Flujos(V, delta, datos_linea, datos_trafo, datos_potencia, Sb)
    v = V.*delta;
    I = [];
    Sx = [];
    for i = 1:length(datos_linea(:,1))
        I(i) = (v(datos_linea(i,1))-v(datos_linea(i,2)))/sum(datos_linea(i,3:4)) + V(datos_linea(i,1))*datos_linea(i,5);
    end
    for i = 1:length(datos_linea(:,1))
        if 
    end
    Sd = (datos_potencia(:,5)+1j*datos_potencia(:,6))/Sb;
    
end


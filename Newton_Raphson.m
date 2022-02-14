function [Jacob] = Newton_Raphson( Y_barra, datos_potencia, init, Sb)
    datos_potencia(datos_potencia(:,9) == false,9) = init;
    V = datos_potencia(:,9);
    delta = datos_potencia(:,10);
    Pcalc = zeros(length(Y_barra(2:end,1)),1); Qcalc = zeros(length(Y_barra(2:end,1)),1);    
    for i = 1:length(Y_barra(2:end,1))
        for j = 1:length(Y_barra(1,:))
            Pcalc(i) = Pcalc(i) + abs(Y_barra(i+1,j)*V(i+1)*V(j))*cos(angle(Y_barra(i+1,j))+delta(j)-delta(i+1));
            Qcalc(i) = Qcalc(i) + abs(Y_barra(i+1,j)*V(i+1)*V(j))*sin(angle(Y_barra(i+1,j))+delta(j)-delta(i+1));
        end
    end
    Qcalc = -Qcalc;
    H = []; N = []; L = []; M = [];
    for i = 2:length(Y_barra(:,1))
        for j = 2:length(Y_barra(1,:))
            H(i-1,j-1) = -abs(Y_barra(i,j)*V(i)*V(j))*sin(angle(Y_barra(i,j))+delta(j)-delta(i));
            L(i-1,j-1) = H(i-1,j-1);
            N(i-1,j-1) = abs(Y_barra(i,j)*V(i)*V(j))*cos(angle(Y_barra(i,j))+delta(j)-delta(i));
            M(i-1,j-1) = -N(i-1,j-1);
        end
    end
    datos_potencia(:,3:8) = datos_potencia(:,3:8)./Sb; 
    Pprog = datos_potencia(2:end,3)-datos_potencia(2:end,5);
    Qprog = datos_potencia(2:end,4)-datos_potencia(2:end,6);
    deltaP = Pprog - Pcalc;
    deltaQ = Qprog - Qcalc;
    deltaQ1 = deltaQ(~(datos_potencia(2:end,2) == 1));
    delta_del = deltaQ(datos_potencia(2:end,2) == 1); 
    deltas = [deltaP;deltaQ];
    pos = find((deltas == delta_del)==1);
    Jacob = [H N;M L];
    for i = 1:length(pos)
        try
            Jacob(pos(i),:) = 0;
            Jacob(:,pos(i)) = 0;
        catch ME
            Jacob(pos,:) = 0;
            Jacob(:,pos) = 0;
        end
    end
    deltas_pot = [deltaP;deltaQ1];
    Jacob = reshape(Jacob(Jacob~=0), [], length(deltas_pot));
    count = 1;
    while count <=1
        deltas_v = inv(Jacob)*deltas_pot;
        if count == 1
            deltav_v = deltas_v(length(deltaP)+1:end);
            del_v = deltav_v.*datos_potencia(datos_potencia(2:end,2)~=1,9);
            vs = datos_potencia(datos_potencia(2:end,2)~=1,9)+del_v;
        end
        V(datos_potencia(2:end,2)~=1) = vs;
        delta(2:end) = deltas_v(2:length(deltaP));
        for i = 1:length(Y_barra(2:end,1))
            for j = 1:length(Y_barra(1,:))
                Pcalc(i) = Pcalc(i) + abs(Y_barra(i+1,j)*V(i+1)*V(j))*cos(angle(Y_barra(i+1,j))+delta(j)-delta(i+1));
                Qcalc(i) = Qcalc(i) + abs(Y_barra(i+1,j)*V(i+1)*V(j))*sin(angle(Y_barra(i+1,j))+delta(j)-delta(i+1));
            end
        end
        deltaP = Pprog - Pcalc;
        deltaQ = Qprog - Qcalc;
        count = count+1;
    end
end


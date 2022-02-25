function [V, delta, count] = Newton_Raphson( Y_barra, datos_potencia, init, Sb, metodo)
    datos_potencia(datos_potencia(:,9) == false,9) = init;
    V = datos_potencia(:,9);
    delta = datos_potencia(:,10);
    datos_potencia(:,3:8) = datos_potencia(:,3:8)./Sb;
    Pprog = datos_potencia(2:end,3)-datos_potencia(2:end,5);
    Qprog = datos_potencia(2:end,4)-datos_potencia(2:end,6);
    error = 1e-3;
    count = 1;
    while true
        Pcalc = zeros(length(Y_barra(2:end,1)),1); Qcalc = zeros(length(Y_barra(2:end,1)),1);
        for i = 1:length(Y_barra(2:end,1))
            for j = 1:length(Y_barra(1,:))
                Pcalc(i) = Pcalc(i) + abs(Y_barra(i+1,j)*V(i+1)*V(j))*cos(angle(Y_barra(i+1,j))+delta(j)-delta(i+1));
                Qcalc(i) = Qcalc(i) + abs(Y_barra(i+1,j)*V(i+1)*V(j))*sin(angle(Y_barra(i+1,j))+delta(j)-delta(i+1));
            end
        end
        Qcalc = -Qcalc;
        deltaP = Pprog - Pcalc;
        deltaQ = Qprog - Qcalc;
        %delta_del = deltaQ(datos_potencia(2:end,2) == 1); 
        Qcalc_PV = Qcalc(datos_potencia(2:end,2) == 1);
        limits = datos_potencia(datos_potencia(:,2)==1,7:8);
        posiciones_limits = find(datos_potencia(:,2) == 1);
        Qd_PV = datos_potencia(datos_potencia(:,2)==1,4);
        if count > 1
            for i = 1:length(limits(:,1))
                try
                    if Qcalc_PV(i) + Qd_PV(i) < limits(i,1)
                        datos_potencia(posiciones_limits(i), 2) = 2;
                        pos = find((Qcalc == Qcalc_PV(i))==1);
                        datos_potencia(pos,4) = limits(i,1);
                        Qprog(pos) = limits(i,1);
                        deltaQ(pos) = Qprog(pos) - Qcalc(pos);
                    elseif Qcalc_PV(i) + Qd_PV(i) > limits(i,2)
                        datos_potencia(posiciones_limits(i), 2) = 2;
                        pos = find((Qcalc == Qcalc_PV(i))==1);
                        Qprog(pos) = limits(i,2);
                        deltaQ(pos) = Qprog(pos) - Qcalc(pos);
                    end
                catch ME
                    if Qcalc_PV + Qd_PV < limits(i,1)
                        datos_potencia(posiciones_limits, 2) = 2;
                        pos = find((Qcalc == Qcalc_PV)==1);
                        datos_potencia(pos,4) = limits(i,1);
                        Qprog(pos) = limits(i,1);
                        deltaQ(pos) = Qprog(pos) - Qcalc(pos);
                    elseif Qcalc_PV(i) + Qd_PV > limits(i,2)
                        datos_potencia(posiciones_limits(i), 2) = 2;
                        pos = find((deltaQ == Qcalc_PV)==1);
                        datos_potencia(pos,4) = limits(i,2);
                        Qprog(pos) = limits(i,2);
                        deltaQ(pos) = Qprog(pos) - Qcalc(pos);
                    end
                end
            end
        end
        delta_del = deltaQ(datos_potencia(2:end,2) == 1);
        deltaQ1 = deltaQ(~(datos_potencia(2:end,2) == 1));
        H = []; N = []; L = []; M = [];
        for i = 2:length(Y_barra(:,1))
            for j = 2:length(Y_barra(1,:))
                if i ~= j
                    H(i-1,j-1) = -abs(Y_barra(i,j)*V(i)*V(j))*sin(angle(Y_barra(i,j))+delta(j)-delta(i));
                    L(i-1,j-1) = H(i-1,j-1);
                    N(i-1,j-1) = abs(Y_barra(i,j)*V(i)*V(j))*cos(angle(Y_barra(i,j))+delta(j)-delta(i));
                    M(i-1,j-1) = -N(i-1,j-1);
                else
                    H(i-1,j-1) = -Qcalc(i-1)-imag(Y_barra(i,j))*V(i)^2;
                    L(i-1,j-1) = Qcalc(i-1)-imag(Y_barra(i,j))*V(i)^2;
                    N(i-1,j-1) = Pcalc(i-1)+real(Y_barra(i,j))*V(i)^2;
                    M(i-1,j-1) = Pcalc(i-1)-real(Y_barra(i,j))*V(i)^2;
                end
            end
        end
        deltas_pot = [deltaP;deltaQ1];
        switch metodo
            case "normal"
                Jacob = [H N;M L];
                if ~isempty(delta_del)
                    iaux2 = logical([zeros(size(deltaP));datos_potencia(2:end,2)==1]);
                    Jacob = Jacob(~iaux2,~iaux2); 
                end
                deltas_v = inv(Jacob)*deltas_pot;
                delta_deltas = deltas_v(1:length(deltaP));
                deltav_v = deltas_v(length(deltaP)+1:end);
            case "NRD"
                delta_deltas = inv(H)*deltaP;
                deltav_v = inv(L(~(datos_potencia(2:end,2) == 1), ~(datos_potencia(2:end,2) == 1)))*deltaQ1;
            case "NRDR"
                deltas_pot = deltas_pot./V;
                delta_deltas = inv(-imag(deltasY_barra(2:end,2:end)))*deltas_pot(1:length(deltaP));
                deltav_v = (-imag(x(~(datos_potencia(2:end,2)==1),~(datos_potencia(2:end,2)==1))))*deltas_pot(length(deltaP)+1:end);
        end
        del_v = deltav_v.*datos_potencia(datos_potencia(2:end,2)~=1,9);
        V(datos_potencia(1:end,2)== 2) = V(datos_potencia(1:end,2)== 2)+del_v;
        delta(2:end) = delta_deltas+delta(2:end);
        if abs(deltas_pot) < error
            break
        elseif count > 20
            break
        else
            deltas_pot
        end
        count = count + 1;
    end
    count = count-1;
end



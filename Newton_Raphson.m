function [V, delta] = Newton_Raphson( Y_barra, datos_potencia, init, Sb)
    datos_potencia(datos_potencia(:,9) == false,9) = init;
    V = datos_potencia(:,9);
    delta = datos_potencia(:,10);
    Pcalc = zeros(length(Y_barra(2:end,1)),1); Qcalc = zeros(length(Y_barra(2:end,1)),1);
    datos_potencia(:,3:8) = datos_potencia(:,3:8)./Sb;
    Pprog = datos_potencia(2:end,3)-datos_potencia(2:end,5);
    Qprog = datos_potencia(2:end,4)-datos_potencia(2:end,6);
    error = 1e-3;
    count = 1;
    while true
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
        deltas = [deltaP;deltaQ];
        if count > 1
            for i = 1:length(limits(:,1))
                try
                    if Qcalc_PV(i) < limits(i,1)
                        datos_potencia(posiciones_limits(i), 2) = 2;
                        pos = find((Qcalc == Qcalc_PV(i))==1);
                        datos_potencia(pos,4) = limits(i,1);
                        Qprog(pos) = limits(i,1);
                        deltaQ(pos) = Qprog(pos) - Qcalc(pos);
                    elseif Qcalc_PV(i) > limits(i,2)
                        datos_potencia(posiciones_limits(i), 2) = 2;
                        pos = find((Qcalc == Qcalc_PV(i))==1);
                        Qprog(pos) = limits(i,2);
                        deltaQ(pos) = Qprog(pos) - Qcalc(pos);
                    end
                catch ME
                    if Qcalc_PV < limits(i,1)
                        datos_potencia(posiciones_limits, 2) = 2;
                        pos = find((Qcalc == Qcalc_PV)==1);
                        datos_potencia(pos,4) = limits(i,1);
                        Qprog(pos) = limits(i,1);
                        deltaQ(pos) = Qprog(pos) - Qcalc(pos);
                    elseif Qcalc_PV(i) > limits(i,2)
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
        if ~isempty(delta_del)
            pos_del = find((deltas == delta_del)==1);
        else
            pos_del = [];
        end
        H = []; N = []; L = []; M = [];
        for i = 2:length(Y_barra(:,1))
            for j = 2:length(Y_barra(1,:))
                H(i-1,j-1) = -abs(Y_barra(i,j)*V(i)*V(j))*sin(angle(Y_barra(i,j))+delta(j)-delta(i));
                L(i-1,j-1) = H(i-1,j-1);
                N(i-1,j-1) = abs(Y_barra(i,j)*V(i)*V(j))*cos(angle(Y_barra(i,j))+delta(j)-delta(i));
                M(i-1,j-1) = -N(i-1,j-1);
            end
        end
        Jacob = [H N;M L];
        if ~isempty(pos_del)
            for i = 1:length(pos_del)
                try
                    Jacob(pos_del(i),:) = 0;
                    Jacob(:,pos_del(i)) = 0;
                catch ME
                    Jacob(pos_del,:) = 0;
                    Jacob(:,pos_del) = 0;
                end
            end
        end
        deltas_pot = [deltaP;deltaQ1];
        if deltas_pot < error
            break
        end
        if ~isempty(pos_del)
            Jacob = reshape(Jacob(Jacob~=0), [], length(deltas_pot));
        end
        disp(Jacob)
        if count > 5
            break
        end
        deltas_v = inv(Jacob)*deltas_pot;
        deltav_v = deltas_v(length(deltaP)+1:end);
        del_v = deltav_v.*datos_potencia(datos_potencia(2:end,2)~=1,9);
        vs = datos_potencia(datos_potencia(2:end,2)~=1,9)+del_v;
        V(datos_potencia(1:end,2)== 2) = vs;
        delta(2:end) = deltas_v(1:length(deltaP));  
        count = count + 1;
    end 
end



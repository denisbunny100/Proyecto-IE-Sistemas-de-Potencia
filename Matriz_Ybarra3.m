function [y_barra, r2p, r3p] = Matriz_Ybarra3(datos_linea,n,datos_trafo,pasos)
    %numero de barras
    %datos_linea de la matriz - tiene dos columnas y n numero de filas, puede
    %tener una columna de más cuando existen compensadores 
    y_barra = zeros(n);
    Y = 1./sum(datos_linea(:,3:4),2);
    datos_linea(:,3) = real(Y);
    datos_linea(:,4) = 1j*imag(Y);
    r2p = [];r3p = [];
    for i = 1:n 
        for j = 1:n
            if i == j
                yes_it_is = datos_linea(:,1) == i | datos_linea(:,2) == i;
                there_is = datos_linea(yes_it_is, :);
                y_barra(i,j) = sum(sum(there_is(:,3:end)));
            else
                yes_it_is = datos_linea(:,1) == i & datos_linea(:,2) == j;
                there_is = datos_linea(yes_it_is,:);
                if ~isempty(there_is)
                    y_barra(i,j) = -sum(there_is(3:4));
                    y_barra(j,i) = y_barra(i,j);
                end
            end    
        end
    end
    count = 1;
    if ~isempty(datos_trafo)
        if sum(datos_trafo(:,3)~=0) ~= 0
            dat_trafos_tri = datos_trafo(datos_trafo(:,3)~=0,:);
            NDDH = [];%Matriz tridevanado
            for i = 1:sum(datos_trafo(:,3)~=0)
                if dat_trafos_tri(i,7) ~= 0
                    c1 = (1+pasos(count)*dat_trafos_tri(i,8))*exp(1j*deg2rad(dat_trafos_tri(i,9)));
                    count = count+1;
                else
                    c1 = exp(1j*deg2rad(dat_trafos_tri(i,9)));
                end
                c2 = exp(1j*deg2rad(dat_trafos_tri(i,10)));
                yp = 1/dat_trafos_tri(i,4);
                ys = 1/dat_trafos_tri(i,5);
                yt = 1/dat_trafos_tri(i,6);
                
                NDDH(1,1) = (yp*(yt+ys))/(yp+ys+yt);
                NDDH(1,2) = -(yp*ys)/(c1*(yp+ys+yt));
                NDDH(1,3) = -(yp*yt)/(c2*(yp+ys+yt));
                NDDH(2,1) = -(c1*yp*ys)/(abs(c1)^2*(yp+ys+yt));
                NDDH(2,2) = (ys*(yp+yt))/(abs(c1)^2*(yp+ys+yt));
                NDDH(2,3) = -(yt*ys*c1)/(c2*abs(c1)^2*(yp+ys+yt));
                NDDH(3,1) = -(yp*yt*c2)/(abs(c2)^2*(yp+ys+yt));
                NDDH(3,2) = -(ys*yt*c2)/(c1*abs(c2)^2*(yp+ys+yt));
                NDDH(3,3) = (yt*(yp+ys))/(abs(c2)^2*(yp+ys+yt));
                
                r3p(:,:,i) = [(yp*(yt+ys))/(yp+ys+yt) -(yp*ys)/(c1*(yp+ys+yt)) -(yp*yt)/(c2*(yp+ys+yt)); 
                             -(c1*yp*ys)/(abs(c1)^2*(yp+ys+yt)) (ys*(yp+yt))/(abs(c1)^2*(yp+ys+yt)) -(yt*ys*c1)/(c2*abs(c1)^2*(yp+ys+yt));
                             -(yp*yt*c2)/(abs(c2)^2*(yp+ys+yt)) -(ys*yt*c2)/(c1*abs(c2)^2*(yp+ys+yt)) (yt*(yp+ys))/(abs(c2)^2*(yp+ys+yt))];
                for z = 1:3
                    for x = 1:3
                        y_barra(dat_trafos_tri(i,z),dat_trafos_tri(i,x)) = y_barra(dat_trafos_tri(i,z),dat_trafos_tri(i,x))+NDDH(z,x); 
                    end
                end 
            end
        end
        x = sum(datos_trafo(:,3)~=0)+1;
        c2p = 1;
        for i = x:length(datos_trafo(:,1))
            if datos_trafo(i,7) ~= 0
                c = (1+pasos(count)*datos_trafo(i,8))*exp(1j*deg2rad(datos_trafo(i,9)));
                count = count+1;
            else
                c = exp(1j*deg2rad(datos_trafo(i,9)));
            end
            y = 1/datos_trafo(i,4); 
            y_barra(datos_trafo(i,1),datos_trafo(i,1)) = y_barra(datos_trafo(i,1),datos_trafo(i,1))+y;
            y_barra(datos_trafo(i,1),datos_trafo(i,2)) = y_barra(datos_trafo(i,1),datos_trafo(i,2))-c*y;
            y_barra(datos_trafo(i,2),datos_trafo(i,1)) = y_barra(datos_trafo(i,2),datos_trafo(i,1))-y*conj(c);
            y_barra(datos_trafo(i,2),datos_trafo(i,2)) = y_barra(datos_trafo(i,2),datos_trafo(i,2))+y*abs(c)^2;
            r2p(:,:,c2p) = [y -c*y; -y*conj(c) y*abs(c)^2];
            c2p = c2p + 1;
        end
    end
end


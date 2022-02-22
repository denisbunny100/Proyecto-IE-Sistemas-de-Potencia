function y_barra = Matriz_Ybarra(datos_linea,n,contin,datos_trafo,pasos)
    %numero de barras
    %datos_linea de la matriz - tiene dos columnas y n numero de filas, puede
    %tener una columna de más cuando existen compensadores 
    y_barra = zeros(n);
    Y = 1./sum(datos_linea(:,3:4),2);
    datos_linea(:,3) = real(Y);
    datos_linea(:,4) = 1j*imag(Y);
    if contin
        datos_linea(contin,3:5) = 0;
    end
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
        for i = 1:length(datos_trafo(:,1))
            if datos_trafo(i,7) ~= 0
                if datos_trafo(i,3) ~= 0
                    c = (1+pasos(count)*datos_trafo(i,8))*exp(1j*deg2rad(datos_trafo(i,9:11)));
                    y = datos_trafo(i,4:6);
                    B = -c.*y;
                    C = -y.*conj(c);
                    D = y.*abs(c).^2; 
                    for a = 1:length(c)-1
                        for u = a+1:length(c)
                            A1 = y(u) + y(a); 
                            A2 = B(u) + B(a);
                            A3 = C(u) + C(a); 
                            A4 = D(u) + D(a);
                            y_barra(datos_trafo(i,1),datos_trafo(i,1)) = y_barra(datos_trafo(i,1),datos_trafo(i,1))+A1;
                            y_barra(datos_trafo(i,1),datos_trafo(i,2)) = y_barra(datos_trafo(i,1),datos_trafo(i,2))+A2;
                            y_barra(datos_trafo(i,2),datos_trafo(i,1)) = y_barra(datos_trafo(i,2),datos_trafo(i,1))+A3;
                            y_barra(datos_trafo(i,2),datos_trafo(i,2)) = y_barra(datos_trafo(i,2),datos_trafo(i,2))+A4;
                        end
                    end
                else
                    y = 1/(1j*datos_trafo(i,4));
                    c = (1+pasos(count)*datos_trafo(i,8))*exp(1j*deg2rad(datos_trafo(i,9)));
                    y_barra(datos_trafo(i,1),datos_trafo(i,1)) = y_barra(datos_trafo(i,1),datos_trafo(i,1))+y;
                    y_barra(datos_trafo(i,1),datos_trafo(i,2)) = y_barra(datos_trafo(i,1),datos_trafo(i,2))-c*y;
                    y_barra(datos_trafo(i,2),datos_trafo(i,1)) = y_barra(datos_trafo(i,2),datos_trafo(i,1))-y*conj(c);
                    y_barra(datos_trafo(i,2),datos_trafo(i,2)) = y_barra(datos_trafo(i,2),datos_trafo(i,2))+y*abs(c)^2;
                end
                count = count+1;
            end
        end
    end
end
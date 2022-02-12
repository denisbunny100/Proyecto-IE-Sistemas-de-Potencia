function y_barra = Matriz_Ybarra(datos,n,contin)
    %numero de barras
    %datos de la matriz - tiene dos columnas y n numero de filas, puede
    %tener una columna de más cuando existen compensadores 
    y_barra = zeros(n);
    Y = 1./sum(datos(:,3:4),2);
    datos(:,3) = real(Y);
    datos(:,4) = 1j*imag(Y);
    if contin
        datos(contin,3:5) = 0;
    end
    for i = 1:n 
        for j = 1:n
            if i == j
                yes_it_is = datos(:,1) == i | datos(:,2) == i;
                there_is = datos(yes_it_is, :);
                y_barra(i,j) = sum(sum(there_is(:,3:end)));
            else
                yes_it_is = datos(:,1) == i & datos(:,2) == j;
                there_is = datos(yes_it_is,:);
                if ~isempty(there_is)
                    y_barra(i,j) = -sum(there_is(3:4));
                    y_barra(j,i) = y_barra(i,j);
                end
            end    
        end
    end
end
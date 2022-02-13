function Y = Impedancias_tridevanado(X_tridevanado, Sb)
%CAMBIO_B Summary of this function goes here
%   Detailed explanation goes here
    X_new = (X_tridevanado(1,:)./X_tridevanado(2,:))*Sb;
    Xp = (1/2)*(sum(X_new)-2*X_new(2));
    Xs = (1/2)*(sum(X_new)-2*X_new(3));
    Xt = (1/2)*(sum(X_new)-2*X_new(1));
    X = [Xp Xs Xt];
    Y = 1./(1j*X);
end


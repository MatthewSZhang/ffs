% Fisher's Iris
clear
clc

load fisheriris

X = meas([1, 2, 51, 52, 101, 102, 103], :);
Y = [1, 0;
     1, 0;
     0, 1;
     0, 1;
     0, 0;
     0, 0;
     0, 0];
 
n = size(X, 2);
 
%% Step 1
Xc = X - mean(X);
Yc = Y - mean(Y); 

%% Step 2
[U, S, Rv] = svd([Xc, Yc], 'econ');
F = S*Rv';
Fx = F(:, 1:n);
Fy = F(:, n+1:end);


V1 = Fy(:, 1);
V2 = Fy(:, 2) - Fy(:, 2)'*V1./(V1'*V1)*V1;


%% Step 3
Xs = [];
Fxr = X;

Fxs = [];
Fxr = Fx;

Wr = Fxr;

%% Step 4
t11 = V1'*Wr(:, 1)*Wr(:, 1)'*V1./(Wr(:, 1)'*Wr(:, 1)*(V1'*V1));
t12 = V2'*Wr(:, 1)*Wr(:, 1)'*V2./(Wr(:, 1)'*Wr(:, 1)*(V2'*V2));
t21 = V1'*Wr(:, 2)*Wr(:, 2)'*V1./(Wr(:, 2)'*Wr(:, 2)*(V1'*V1));
t22 = V2'*Wr(:, 2)*Wr(:, 2)'*V2./(Wr(:, 2)'*Wr(:, 2)*(V2'*V2));
t31 = V1'*Wr(:, 3)*Wr(:, 3)'*V1./(Wr(:, 3)'*Wr(:, 3)*(V1'*V1));
t32 = V2'*Wr(:, 3)*Wr(:, 3)'*V2./(Wr(:, 3)'*Wr(:, 3)*(V2'*V2));
t41 = V1'*Wr(:, 4)*Wr(:, 4)'*V1./(Wr(:, 4)'*Wr(:, 4)*(V1'*V1));
t42 = V2'*Wr(:, 4)*Wr(:, 4)'*V2./(Wr(:, 4)'*Wr(:, 4)*(V2'*V2));

R2xr1Y = t11 + t12;
R2xr2Y = t21 + t22;
R2xr3Y = t31 + t32;
R2xr4V = t41 + t42;

%% Step 3
Xs = X(:, 3);
Fxr = X(:, [1, 2, 4]);

Fxs = Fx(:, 3);
Fxr = Fx(:, [1, 2, 4]);

Ws = Fxs;
Wr = [];
Wr(:, 1) = Fxr(:, 1) - Fxr(:, 1)'*Ws./(Ws'*Ws)*Ws;
Wr(:, 2) = Fxr(:, 2) - Fxr(:, 2)'*Ws./(Ws'*Ws)*Ws;
Wr(:, 3) = Fxr(:, 3) - Fxr(:, 3)'*Ws./(Ws'*Ws)*Ws;

%% Step 4
t11 = V1'*Wr(:, 1)*Wr(:, 1)'*V1./(Wr(:, 1)'*Wr(:, 1)*(V1'*V1));
t12 = V2'*Wr(:, 1)*Wr(:, 1)'*V2./(Wr(:, 1)'*Wr(:, 1)*(V2'*V2));
t21 = V1'*Wr(:, 2)*Wr(:, 2)'*V1./(Wr(:, 2)'*Wr(:, 2)*(V1'*V1));
t22 = V2'*Wr(:, 2)*Wr(:, 2)'*V2./(Wr(:, 2)'*Wr(:, 2)*(V2'*V2));
t31 = V1'*Wr(:, 3)*Wr(:, 3)'*V1./(Wr(:, 3)'*Wr(:, 3)*(V1'*V1));
t32 = V2'*Wr(:, 3)*Wr(:, 3)'*V2./(Wr(:, 3)'*Wr(:, 3)*(V2'*V2));



R2xr1Y = t11 + t12;
R2xr2Y = t21 + t22;
R2xr3Y = t31 + t32;

%% Step 3
Xs = X(:, [3, 4]);
Fxr = X(:, [1, 2]);

Fxs = Fx(:, [3, 4]);
Fxr = Fx(:, [1, 2]);

Ws = [];
Ws(:, 1) = Fxs(:, 1);
Ws(:, 2) = Fxs(:, 2) - Fxs(:, 2)'*Ws(:, 1)./(Ws(:, 1)'*Ws(:, 1))*Ws(:, 1);

Wr = [];
Wr(:, 1) = Fxr(:, 1) - Fxr(:, 1)'*Ws(:, 1)./(Ws(:, 1)'*Ws(:, 1))*Ws(:, 1) - Fxr(:, 1)'*Ws(:, 2)./(Ws(:, 2)'*Ws(:, 2))*Ws(:, 2);
Wr(:, 2) = Fxr(:, 2) - Fxr(:, 2)'*Ws(:, 1)./(Ws(:, 1)'*Ws(:, 1))*Ws(:, 1) - Fxr(:, 2)'*Ws(:, 2)./(Ws(:, 2)'*Ws(:, 2))*Ws(:, 2);

%% Step 4
t11 = V1'*Wr(:, 1)*Wr(:, 1)'*V1./(Wr(:, 1)'*Wr(:, 1)*(V1'*V1));
t12 = V2'*Wr(:, 1)*Wr(:, 1)'*V2./(Wr(:, 1)'*Wr(:, 1)*(V2'*V2));
t21 = V1'*Wr(:, 2)*Wr(:, 2)'*V1./(Wr(:, 2)'*Wr(:, 2)*(V1'*V1));
t22 = V2'*Wr(:, 2)*Wr(:, 2)'*V2./(Wr(:, 2)'*Wr(:, 2)*(V2'*V2));

R2xr1Y = t11 + t12;
R2xr2Y = t21 + t22;

%% CCA
[~, ~, Rx342Y] = canoncorr(X(:, [3, 4, 2]),Y);
R2x342Y = Rx342Y.^2;


[IND_ols, criteria] = ffs(X, Y, 3);

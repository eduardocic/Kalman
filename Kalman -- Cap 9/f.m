function xd = f(x, t)

% Gravidade.
g = 32.2;

% Problema no âmbito contínuo.
A = [0 1 0 0;
     0 0 0 0;
     0 0 0 1;
     0 0 0 0];
     
B = [ 0;
      0;
      0;
     -g];

xd = A*x + B;

end
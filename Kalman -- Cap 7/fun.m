function x_d = fun(beta, x, t)
x_d = zeros(2,1);

g = 32.2;
    
x_d(1) = x(2);
x_d(2) = 0.0034*g*exp(-x(1)/22000)*(x(2)^2)/(2*beta) - g; 
    
end
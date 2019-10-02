function x_out = RK2(f, x_in, t, h)

x_out = x_in + (0.5)*h*(f(x_in, t) + f(x_in, t + h) );
   
end
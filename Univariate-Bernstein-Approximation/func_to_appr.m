function y = func_to_appr(x)
%y = x.^2.2;
y = x.^.45;
%y = (sin(4*x)+2)/4;
%y = (sin(8*x)+2)/4;
%y = ((x+0.055)/1.055).^2.4;
%y = x.^(1/2.4) * 1.055 - 0.055;
%y = cos(x) .* cos(sin(x));
%y = (1 + 0.3*cos(5*x)).*sin(x);
%y = (1 + 0.3*cos(5*x)).*cos(x)/2;
%y = exp(-3*x);
%y = log(x+1.5);
%y = tanh(2*x);
% end of function func_to_appr

function dxdt = plant(t,x,K,M,L,F,g)

dxdt = zeros(4,1); u = K*x;

dxdt(1) = x(2);
dxdt(2) = 1/M*u-F/M*x(2);
dxdt(3) = x(4);
dxdt(4) = g/L*sin(x(3))-1/(L*M)*u*cos(x(3))+F/(L*M)*x(2)*cos(x(3));

end
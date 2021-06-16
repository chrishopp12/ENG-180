%Chris Hopp
%915866326
%ENG-180 Project 10: Eigenvalues
%12/17/2020


clc
clf
clear all


Problem1a()
Problem1b()
Problem1c()


%% Problem 1a
function Problem1a()
E = 1;
I = 1;
element = [10,20,40];

for i=1:3
[lambda,u] = InversePower(element(i));
y = linspace(0,1,element(i)+1);
figure(1)
subplot(1,3,i)
plot(u,y)
xlim([0,.5])
xlabel('u')
ylabel('y')
title(['Buckling: ', num2str(element(i)), ' Elements'])


eigenvector(i) = lambda;
P(i) = lambda*E*I;
end

display(eigenvector)
display(P)

end

%% Problem 1b
function Problem1b()

[lambdaMax, uMax] = Power(40);
[lambdaMin, uMin] = InversePower(40);

x = linspace(0,1,41);
figure(2)
subplot(2,2,1)
plot(x,uMax)
xlabel('x')
ylabel('T_0')
title(['Heat Equation: \lambda = ', num2str(lambdaMax)])

subplot(2,2,2)
plot(x,uMin)
xlabel('x')
ylabel('T_0')
title(['Heat Equation: \lambda = ', num2str(lambdaMin)])


k = .1;
alphaMin = k*lambdaMin;
alphaMax = k*lambdaMax;

t = x;
t2 = linspace(0,.01,41);            % Shorter t vector for exponential decay


for i= 1:41                             % Builds temperature plots
TMax(i,:) = exp(-alphaMax*t2(i)).*uMax;
TMin(i,:) = exp(-alphaMin*t(i)).*uMin;
end


figure(2)
subplot(2,2,3)
pcolor(x,t2,TMax)
shading interp
xlabel('x')
ylabel('t')
cb = colorbar;
cb.Label.String = "T";
title(['Heat Equation: \lambda = ', num2str(lambdaMax)])


subplot(2,2,4)
pcolor(x,t,TMin)
xlabel('x')
ylabel('t')
cb = colorbar;
cb.Label.String = "T";
title(['Heat Equation: \lambda = ', num2str(lambdaMin)])
shading interp


end


%% Problem 1c

function Problem1c()
[lambda, u] = InversePower(50);
x = linspace(0,1,51);
t = x;
c = 10;
omega = c*sqrt(lambda);
frequency = omega/(2*pi);
display(frequency)

for i= 1:51                           
U(i,:) = sin(omega*t(i)).*sin(sqrt(lambda)*x);
end

figure(3)
subplot(1,2,1)
plot(x,u)
xlabel('x')
ylabel('u_0')
title('String vibration: Mode 1')

subplot(1,2,2)
pcolor(x,t,U)
xlabel('x')
ylabel('t')
title(['String Vibration'])
cb = colorbar;
cb.Label.String = "Amplitude";
shading interp

end
%% Inverse power method
function [lambda, u] = InversePower(elements)
cor = 1;
iteration = 0;
tol = 1e-10;    
dy = 1/elements;
n = elements + 1;
uNew = ones(n,1);                   % Initial guess
uNew(1) = 0;
uNew(end) = 0;
uNew = uNew/(sqrt(uNew'*uNew));     % Normalizes guess


while cor > tol 
    iteration = iteration + 1;
    uOld = uNew;
    
    a = -ones(n-1,1)*1/dy^2;      % Forms tridiagonal system from central difference
    b = ones(n,1)*(2/dy^2);
    c = a;
    
    a(end) = 0;                   % Boundary conditions
    b(1) = 1;
    b(end) = 1;
    c(1) = 0;
    
    d = uOld;
    d(1) = 0;
    d(end) = 0;
    
    w = Thomas3(a,b,c,d);       
    w = w';
    A = diag(b) + diag(a,-1) + diag(c,1);
    uNew = w/(sqrt(w'*w));          % Normalizes solution
    cor = max(abs(uNew-uOld));
 
end

lambda = uNew'*A*uNew;
u = uNew;

end

%% Power method
function [lambda, u] = Power(elements)
cor = 1;
iteration = 0;
tol = 1e-10;    
dy = 1/elements;
n = elements + 1;
uNew = ones(n,1);                   % Initial guess
uNew(1) = 0;
uNew(end) = 0;
uNew = uNew/(sqrt(uNew'*uNew));     % Normalizes guess


while cor > tol 
    iteration = iteration + 1;
    uOld = uNew;
    
    a = -ones(n-1,1)*1/dy^2;      % Forms tridiagonal system from central difference
    b = ones(n,1)*(2/dy^2);
    c = a;
    
    a(end) = 0;                   % Boundary conditions
    b(1) = 1;
    b(end) = 1;
    c(1) = 0;
    
    A = diag(b) + diag(a,-1) + diag(c,1);  % Forms A matrix
    
    w = A*uOld;       

    uNew = w/(sqrt(w'*w));          % Normalizes solution
    cor = max(abs(uNew-uOld));
 
end

lambda = uNew'*A*uNew;
u = uNew;

end




%% Tridiagonal Thomas Algorithm Function
function xBar = Thomas3(a,b,c,d)
    n = length(b);             % Length of total diagonal
    aBar = [0;a];       % Forms vectors of n length for manipulated values
    bBar = b;
    cBar = [c;0];
    dBar = d;

 % "Zip-down" eliminates subdiagonal by subtracting previous row scaled to 'a' term
    for i=2:n     
        bBar(i) = b(i) - aBar(i)*cBar(i-1)/bBar(i-1);
        dBar(i) = d(i) - aBar(i)*dBar(i-1)/bBar(i-1);
    end
    
    % "Zip-up" solves for x from last row single variable equation to first row
    xBar(n) = dBar(n)/ bBar(n);         
    for i = n-1:-1:1
        xBar(i) = ((dBar(i)-(cBar(i)*xBar(i+1)))/bBar(i));
    end
end

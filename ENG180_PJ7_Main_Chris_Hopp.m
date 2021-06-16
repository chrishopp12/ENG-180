%Chris Hopp
%915866326
%ENG-180 Project 7: Ordinary Differential Equations: IVP
%11/18/2020


clc
clf
clear all


Problem1a()
Problem1b()
Problem1ci()
Problem1cii()
Problem2()
Problem3()
Problem4()
%% Problem 1a
function Problem1a()
dt = 0.01;
t = [0:dt:10];
n = length(t);
v = zeros(n,1);
x = v;
v(1) = 1;
x(1) = 0;
m = 1;
g = 2;
c1 = m*g;
c2 = m;

% Analytical
xAnalytical = (m/c2)*((c1/c2)-1)*exp(-c2*t/m)-((m/c2)*((c1/c2)-1))+(c1/c2)*t;
vAnalytical = -((c1/c2)-1)*exp(-c2*t/m) + c1/c2;

% Forward Euler
vForward = v;
xForward = x;
for i = 1:n-1
    vForward(i+1) = ((dt/m)*(c1-c2*vForward(i)))+vForward(i);
    xForward(i+1) = dt*vForward(i)+xForward(i);
end

% Backward Euler
vBackward = v;
xBackward = x;
for i = 1:n-1
    vBackward(i+1) = (m/(m+(dt*c2)))*(((dt*c1)/m)+vBackward(i));
    xBackward(i+1) = dt*vBackward(i+1)+xBackward(i);
end

% Leap Frog
vLeap = v;
xLeap = x;
vLeap(2) = vLeap(1)+(dt/m)*(c1-c2*vLeap(1));
xLeap(2) = vLeap(1)*dt + xLeap(1);
for i = 2:n-1
    vLeap(i+1) = (((2*dt)/m)*(c1-c2*vLeap(i)))+vLeap(i-1);
    xLeap(i+1) = 2*dt*vLeap(i)+xLeap(i-1);
end

fV = @(V) (1/m)*(c1-c2*V);   % I should have done this sooner
fX = @(V) V;

% Crank Nicolson
vCrank = v;
xCrank = x;

for i = 1:n-1
    vCrank(i+1) = (1/(1+(dt*c2/(2*m)))) * ( (dt*c1/(2*m))+ (dt/2)*fV(vCrank(i))+vCrank(i));
    xCrank(i+1) = ((dt/2) * (fX(vCrank(i+1)) + fX(vCrank(i)))) + xCrank(i);
end


fV = @(V,t) (1/m)*(c1-c2*V) + 0*t;   % Had to add the term to make more general RK functions
fX = @(V,t) V + 0*t;

% Runge-Kutta 2
vRK2 = RK2(fV,1,0,dt,n);
xRK2 = RK2(fX,0,0,dt,n,vRK2);

% Runge-Kutta 3
vRK3 = RK3(fV,1,0,dt,n);
xRK3 = RK3(fX,0,0,dt,n,vRK3);

% Runge-Kutta 4
vRK4 = RK4(fV,1,0,dt,n);
xRK4 = RK4(fX,0,0,dt,n,vRK4);

% Quadratic Explicit
fV = @(V) (1/m)*(c1-c2*V);   
fX = @(V) V;
vQuad = v;
xQuad = x;

for i = 1:n-1
    vStar = vQuad(i) + (dt*fV(vQuad(i)));
    b = fV(vQuad(i));
    c = (1/(2*dt))*(fV(vStar)-fV(vQuad(i)));
    vQuad(i+1) = vQuad(i) + b*dt + c*dt^2;
    
    bx = fX(vQuad(i));                          % x is function of v, so uses v*
    cx = (1/(2*dt))*(fX(vStar)-fX(vQuad(i)));
    xQuad(i+1) = xQuad(i) + bx*dt + cx*dt^2;
end

% Quadratic Explicit
fV = @(V) (1/m)*(c1-c2*V);   
fX = @(V) V;
vQuad = v;
xQuad = x;

for i = 1:n-1
    vStar = vQuad(i) + (dt*fV(vQuad(i)));
    b = fV(vQuad(i));
    c = (1/(2*dt))*(fV(vStar)-fV(vQuad(i)));
    vQuad(i+1) = vQuad(i) + b*dt + c*dt^2;
    
    bx = fX(vQuad(i));                          % x is function of v, so uses v*
    cx = (1/(2*dt))*(fX(vStar)-fX(vQuad(i)));
    xQuad(i+1) = xQuad(i) + bx*dt + cx*dt^2;
end


% Cubic Explicit
fV = @(V) (1/m)*(c1-c2*V);   
fX = @(V) V;
vCubic = v;
xCubic = x;

for i = 1:n-1
    vStar = vCubic(i) + (dt*fV(vCubic(i)));
    vStar2 = vCubic(i) + (1/2)*dt*(fV(vStar)+fV(vCubic(i)));
    b = fV(vCubic(i));
    c = (1/(2*dt))*(fV(vStar)-fV(vCubic(i)));
    vStar3 = vCubic(i) + (b*dt/2) + (c*dt^2)/4;

    vCubic(i+1) = vCubic(i) + (1/6)*dt*(fV(vCubic(i)) + 4*fV(vStar3) + fV(vStar2));
  
    xStar = vStar;          %x is a function of v, so uses v*, v**, v***
    xStar2 = vStar2;
    xStar3 = vStar3;
    bx = fX(vCubic(i));
    cx = (1/(2*dt))*(fX(xStar)-fX(vCubic(i)));

    
    xCubic(i+1) = xCubic(i) + (1/6)*dt*(fX(vCubic(i)) + 4*fX(xStar3) + fX(xStar2));
    
end




figure(1)
subplot(1,2,1)
plot(t,xAnalytical,'DisplayName', 'Analytical')
hold on
plot(t,xForward,'DisplayName', 'Forward Euler')
plot(t,xBackward,'DisplayName', 'Backward Euler')
plot(t,xLeap,'DisplayName', 'Leap Frog')
plot(t,xCrank,'DisplayName', 'Crank Nicolson')
plot(t,xRK2,'DisplayName', 'Runge-Kutta 2')
plot(t,xRK3,'DisplayName', 'Runge-Kutta 3')
plot(t,xRK4,'DisplayName', 'Runge-Kutta 4')
plot(t,xQuad,'DisplayName', 'Quadratic Collocation')
plot(t,xCubic,'DisplayName', 'Cubic Collocation')
title('Position')
xlabel('Time')
ylabel('x')
legend('Location','best')

subplot(1,2,2)
plot(t,vAnalytical,'DisplayName', 'Analytical')
hold on
plot(t,vForward,'DisplayName', 'Forward Euler')
plot(t,vBackward,'DisplayName', 'Backward Euler')
plot(t,vLeap,'DisplayName', 'Leap Frog')
plot(t,vCrank,'DisplayName', 'Crank Nicolson')
plot(t,vRK2,'DisplayName', 'Runge-Kutta 2')
plot(t,vRK3,'DisplayName', 'Runge-Kutta 3')
plot(t,vRK4,'DisplayName', 'Runge-Kutta 4')
plot(t,vQuad,'DisplayName', 'Quadratic Collocation')
plot(t,vCubic,'DisplayName', 'Cubic Collocation')
title('Velocity')
xlabel('Time')
ylabel('v')
legend('Location','best')
hold off
end

%% Problem 1b
function Problem1b()
dt = 0.01;
t = [0:dt:50];
n = length(t);
x = zeros(n,1);
vx = x;
y=x;
vy = x;
x(1) = 1;           % Initial conditions
vx(1) = 0;
y(1) = 0;
vy(1) = 1;
mu1 = [.4,.5,.6,1];

for j = 1:4   % Iterates through mu values
   mu = mu1(j); 

    yVec = [vx(1);vy(1);x(1);y(1)];
    
    w =  mu/((x(1)^2+y(1)^2)^(3/2));     % w = w^2, initital

    for i = 1:n-1                       % Solves via RK4
        f = @(yVec)[yVec(3)*(-w);yVec(4)*(-w);yVec(1);yVec(2)];
        yVec = RK4Step(dt,f,yVec);      
        vx(i+1) = yVec(1);
        vy(i+1) = yVec(2);
        x(i+1) = yVec(3);
        y(i+1) = yVec(4);
        w = mu/((x(i+1)^2+y(i+1)^2)^(3/2));         % Updates with steps in x,y
    
    end
    
    H = ((1/2)*((vx.^2).*(vy.^2)))-(x.^2 + y.^2).^(-1/2);
    M = x.*vy - y.*vx;
    
    figure(2)
    subplot(2,2,j)
    plot(x,y)
    xlabel('x')
    ylabel('y')
    title({'Orbital Mechanics',['\mu = ', num2str(mu)]})
    
    figure(3)
    subplot(2,2,j)
    plot(t,H,'DisplayName', 'H')
    hold on
    plot(t,M,'DisplayName', 'M')
    xlabel('t')
    ylabel('H, M')
    legend('location','best')
    ylim([-1,1.5])
    xlim([0,50])
    title({'Orbital Mechanics',['\mu = ', num2str(mu)]})
    
end
    
end

%% Problem 1c i
function Problem1ci()
% Duffings
dt = .01;
time = [0:dt:50];
n = length(time);
x = zeros(n,1);
v = x;
x(1) = 0;           % Initial conditions
v(1) = 1;

delta = .22;
w = 1;
beta = 1;
a = 1;
fVar = .3;

yVec = [v(1);x(1)];


    for i = 1:n-1                       % Solves via RK4
        f = @(yVec,t) [fVar*cos(w*t) - delta*yVec(1) + beta*(yVec(2)) - a*yVec(2)^3; yVec(1)];
        yVec = RK4Step(dt,f,yVec,time(i));     
        v(i+1) = yVec(1);
        x(i+1) = yVec(2);
    end
    figure(4)
        subplot(1,2,1)
    plot(time,x,'DisplayName','x')
    hold on
    plot(time,v,'DisplayName','v')
    xlabel('Time')
    ylabel('x,v')
    title("Duffing's Equation")
    legend('location','best')
    xlim([0,50])
    hold off
    subplot(1,2,2)
    plot(x,v)
    xlabel('x')
    ylabel('v')
    title(["Duffing's Equation","v vs x"])

end
%% Problem 1c ii
function Problem1cii()
% Van der Pol
dt = .01;
time = [0:dt:50];
n = length(time);
x = zeros(n,1);
v = x;
x(1) = 0;           % Initial conditions
v(1) = 1;

epsilon = 10;

yVec = [v(1);x(1)];


    for i = 1:n-1                       % Solves via RK4
        f = @(yVec) [-epsilon*((yVec(2)^2 -1)*yVec(1))-yVec(2); yVec(1)];
        yVec = RK4Step(dt,f,yVec);     
        v(i+1) = yVec(1);
        x(i+1) = yVec(2);   
    end
    figure(5)
        subplot(1,2,1)
    plot(time,x,'DisplayName','x')
    hold on
    plot(time,v,'DisplayName','v')
    xlabel('Time')
    ylabel('x,v')
    title("Van del Pol Equation")
    legend('location','best')
    xlim([0,50])
    hold off
    subplot(1,2,2)
    plot(x,v)
    xlabel('x')
    ylabel('v')
    title(["Van der Pol Equation","v vs x"])

end


%% Problem 2 
function Problem2()
% Part c) Heat Radiation
dt = .01;
time = [0:dt:15];
n = length(time);
T = zeros(n,1);

T(1) = 1;           % Initial conditions
Tinf = .1;
beta = 1;



    for i = 1:n-1                       % Solves via RK4
        f = @(T) -beta*(T^4-Tinf^4);
        T(i+1) = RK4Step(dt,f,T(i));       
    end
    figure(6)
    plot(time,T)
    xlabel('Time')
    ylabel('T')
    title("Heat Radiation")


end

%% Problem 3 
function Problem3()
% Part a) Landau Equation
dt = .01;
time = [0:dt:10];
n = length(time);
T = zeros(n,1);

u(1) = 1;           % Initial conditions
a=1;
b=.04;



    for i = 1:n-1                       % Solves via RK4
        f = @(u) a*u - b*u^3;
        u(i+1) = RK4Step(dt,f,u(i));       
    end
    figure(7)
    plot(time,u)
    xlabel('Time')
    ylabel('u')
    title("Landau Equation")


end

%% Problem 4 
function Problem4()
% Part a) Jerk Equation
dt = .01;
time = [0:dt:30];
n = length(time);
x = zeros(n,1);
v = x;
a = x;

x(1) = .02;           % Initial conditions
v(1) = 0;
a(1) = 0;

A = 2.017;


    for i = 1:n-1                       % Solves via RK4
        y = [a(i);v(i);x(i)];
        f = @(y)[-A*y(1)+y(2)^2-y(3);y(1);y(2)];
        y = RK4Step(dt,f,y);
        a(i+1)=y(1);
        v(i+1)=y(2);
        x(i+1)=y(3);
    end
    figure(8)
    plot(x,v)
    xlabel('x')
    ylabel('v')
    xlim([-.25,.4])
    title("Jerk Equation")


end


%% Runge Kutta 2
function y = RK2(f,f0,t0,dt,n,x)
switch nargin
    case 5

    k1 = dt*f(f0,t0);
    k2 = dt*f(f0+k1/2,t0+.5*dt);
    y = zeros(n,1);
    y(1) = f0;
    for i = 1:n-1                               %RK2 algorithm that updates y
        y(i+1) = y(i) + (1/2)*(k1+k2);
        k1 = dt*f(y(i+1),(i*1)*dt);
        k2 = dt*f(y(i+1)+k1/2,(i+1.5)*dt);
    end

    case 6            %Used when y is function of another variable like x_dot = f(v)
    k1 = dt*f(x(1),t0);
    k2 = dt*f(x(1)+k1/2,t0+.5*dt);
    y = zeros(n,1);
    y(1) = f0;
    for i = 1:n-1
        y(i+1) = y(i) + (1/2)*(k1+k2);
        k1 = dt*f(x(i+1),(i*1)*dt);
        k2 = dt*f(x(i+1)+k1/2,(i+1.5)*dt);
    end 
end
end

%% Runge Kutta 3
function y = RK3(f,f0,t0,dt,n,x)
switch nargin
    case 5
                                %RK3 algorithm that updates y
    k1 = dt*f(f0,t0);
    k2 = dt*f(f0+k1/2,t0+.5*dt);
    k3 = dt*f(f0-k1+2*k2,(t0 + 1)*dt);
    y = zeros(n,1);
    y(1) = f0;
    for i = 1:n-1
        y(i+1) = y(i) + (1/6)*(k1+4*k2+k3);
        k1 = dt*f(y(i+1),(i*1)*dt);
        k2 = dt*f(y(i+1)+k1/2,(i+1.5)*dt);
        k3 = dt*f(y(i+1)-k1+2*k2,(i + 2)*dt);
    end

    case 6            %Used when y is function of another variable like x_dot = f(v)
    k1 = dt*f(x(1),t0);
    k2 = dt*f(x(1)+k1/2,t0+.5*dt);
    k3 = dt*f(x(1)-k1+2*k2,(t0 + 1)*dt);
    y = zeros(n,1);
    y(1) = f0;
    for i = 1:n-1
        y(i+1) = y(i) + (1/6)*(k1+4*k2+k3);
        k1 = dt*f(x(i+1),(i*1)*dt);
        k2 = dt*f(x(i+1)+k1/2,(i+1.5)*dt);
        k3 = dt*f(x(i+1)-k1+2*k2,(i + 2)*dt);
    end 
end
end

%% Runge Kutta 4
function y = RK4(f,f0,t0,dt,n,x)
switch nargin
    case 5
                                %RK4 algorithm that updates y
    k1 = dt*f(f0,t0);
    k2 = dt*f(f0+k1/2,t0+.5*dt);
    k3 = dt*f(f0+k2/2,t0+.5*dt);
    k4 = dt*f(f0+k3,(t0 + 1)*dt);
    y = zeros(n,1);
    y(1) = f0;
    for i = 1:n-1
        y(i+1) = y(i) + (1/6)*(k1+2*k2+2*k3+k4);
        k1 = dt*f(y(i+1),(i*1)*dt);
        k2 = dt*f(y(i+1)+k1/2,(i+1.5)*dt);
        k3 = dt*f(y(i+1)+k2/2,(i+1.5)*dt);
        k4 = dt*f(y(i+1)+k3,(i + 2)*dt);
    end

    case 6            %Used when y is function of another variable like x_dot = f(v)
    k1 = dt*f(x(1),t0);
    k2 = dt*f(x(1)+k1/2,t0+.5*dt);
    k3 = dt*f(x(1)+k2/2,t0+.5*dt);
    k4 = dt*f(x(1)+k3,(t0 + 1)*dt);
    y = zeros(n,1);
    y(1) = f0;
    for i = 1:n-1
        y(i+1) = y(i) + (1/6)*(k1+2*k2+2*k3+k4);
        k1 = dt*f(x(i+1),(i*1)*dt);
        k2 = dt*f(x(i+1)+k1/2,(i+1.5)*dt);
        k3 = dt*f(x(i+1)+k2/2,(i+1.5)*dt);
        k4 = dt*f(x(i+1)+k3,(i + 2)*dt);
    end 
end
end


%% Runge Kutta 4 - Stepwise
function y = RK4Step(dt,f,yn,tn)
switch nargin
    case 3
        k1 = dt*f(yn);
        k2 = dt*f(yn+k1/2);
        k3 = dt*f(yn+k2/2);
        k4 = dt*f(yn+k3);
        y = yn + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    case 4
        k1 = dt*f(yn,tn);
        k2 = dt*f(yn+k1/2,tn+.5*dt);
        k3 = dt*f(yn+k2/2,tn+.5*dt);
        k4 = dt*f(yn+k3, tn + dt);
        y = yn + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end
        
end
%Chris Hopp
%915866326
%ENG-180 Project 9: Partial Differential Equations
%12/6/2020


clc
clf
clear all




Problem1a()
Problem1b()
Problem1c()
Problem2a()
Problem2b()

%% Problem 1a
function Problem1a()
nx = 101;
nt = 1001;
x = linspace(-1,1,nx);
t = linspace(0,10,nt);

dx = x(2) - x(1);
dt = t(2) - t(1);

u = zeros(nt,nx);   % u(t,x), seems backwards but this puts x in columns (horizontal)
 
u(1,1:end) = sin(pi.*x);

% c = -1
u1 = u;
for j = 1:nt-1                  % Solves with stencil at j,i+1 level
    u1(j+1,nx) =  u1(j,nx) + ((dt/dx) * ( u1(j,2) - u1(j,nx)));  % Periodic B.C.
    for i=nx-1:-1:1
        u1(j+1,i) = u1(j,i) + ((dt/dx) * ( u1(j,i+1) - u1(j,i)));
    end
end
ux1 = u1(nt,:);

% c = 1
u2 = u;

for j = 1:nt-1                  % Solves with stencil at j, i-1 level
    u2(j+1,1) = u2(j,1) - ((dt/dx) * ( u2(j,1) - u2(j,nx-1)));
    for i=2:nx
        u2(j+1,i) = u2(j,i) - ((dt/dx) * ( u2(j,i) - u2(j,i-1)));
    end
end
ux2 = u2(nt,:);

% c = u 
u3 = u;

for j = 1:nt-1
    for i=1:nx
        if u3(j,i) >= 0                 % Checks if c pos or neg
            if i == 1
                u3(j+1,i) = u3(j,i+1);  % Enforces boundary at endpoints
            elseif i == nx
                u3(j+1,i) = u3(j,i-1);
            else
                u3(j+1,i) = u3(j,i) + u3(j,i)*((dt/dx) * (u3(j,i-1) - u3(j,i)));
            end
        else
            if i == nx                  % Enforces boundary at endpoints
                u3(j+1,i) = u3(j,i-1);
            elseif i == 1
                u3(j+1,i) = u3(j,i+1);
            else
                u3(j+1,i) = u3(j,i) + u3(j,i)*((dt/dx) * ( u3(j,i) - u3(j,i+1))); 
            end
        end
    end
end
ux3 = u3(nt,:);


figure(1)
subplot(3,2,1)
imagesc(x,t,u1)
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('t')
title({'Transport Equation', 'u_0 = sin(\pi x), c = -1'})

subplot(3,2,2)
plot(x,ux1)
xlabel('x')
ylabel('u')
title({'Transport equation t = 10', 'u_0 = sin(\pi x), c = -1'})

subplot(3,2,3)
imagesc(x,t,u2)
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('t')
title({'Transport Equation', 'u_0 = sin(\pi x), c = 1'})

subplot(3,2,4)
plot(x,ux2)
xlabel('x')
ylabel('u')
title({'Transport equation t = 10', 'u_0 = sin(\pi x), c = 1'})


subplot(3,2,5)
imagesc(x,t,u3)
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('t')
title({'Transport Equation', 'u_0 = sin(\pi x), c = u'})

subplot(3,2,6)
plot(x,ux3)
xlabel('x')
ylabel('u')
title({'Transport equation t = 10', 'u_0 = sin(\pi x), c = u'})



end

%% Problem 1b
function Problem1b()
dx = .01;
dt = dx;
x = [-1:dx:1];
t = [0:dt:10];

nx = length(x);
nt = length(t);

u = zeros(nt,nx);   % u(t,x), seems backwards but this puts x in columns (horizontal)
 
u(1,1:end) = -sin(pi.*x);

% c = -1
u1 = u;
for j = 1:nt-1                  % Solves with stencil at j,i+1 level
    u1(j+1,nx) =  u1(j,nx) + ((dt/dx) * ( u1(j,2) - u1(j,nx)));  % Periodic B.C.
    for i=nx-1:-1:1
        u1(j+1,i) = u1(j,i) + ((dt/dx) * ( u1(j,i+1) - u1(j,i)));
    end
end
ux1 = u1(nt,:);

% c = 1
u2 = u;

for j = 1:nt-1                  % Solves with stencil at j, i-1 level
    u2(j+1,1) = u2(j,1) - ((dt/dx) * ( u2(j,1) - u2(j,nx-1)));
    for i=2:nx
        u2(j+1,i) = u2(j,i) - ((dt/dx) * ( u2(j,i) - u2(j,i-1)));
    end
end
ux2 = u2(nt,:);


% c = u
u3 = u;
epsilon = .5*dx;
tol = 1e-8;

for j = 1:nt-1 
    cor = 1;
    uNew = u3(j+1,:);
    uOld = zeros(size(uNew));
    uOld = uOld';
    iteration = 0;
    
    while cor>tol                       % Iterates, lagging u_Old term
        iteration = iteration+1;
        uOld = uNew(:);
                                        % B.C. set with ficticous point and
                                        % enforced at u_1, u_nx
        a = ((-dt/(2*dx))*uOld(2:nx)) - epsilon*dt/dx^2;
        a(end) = (-epsilon*2*dt)/dx^2;
        b = ones(nx,1)*(1 + ((epsilon*2*dt)/dx^2));
        c = ((dt/(2*dx))*uOld(1:nx-1)) - epsilon*dt/dx^2;
        c(1) = (-epsilon*2*dt)/dx^2;
        d = u3(j,:);

        uNew = Thomas3(a,b,c,d);
        uNew = uNew';
        cor = max(abs(uNew-uOld));
    end
    
    u3(j+1,:) = uNew(:);                 % Entire row solved at same time
end

ux3 = u3(nt,:);


figure(2)
subplot(3,2,1)
imagesc(x,t,u1)
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('t')
title({'Transport Equation', 'u_0 = -sin(\pi x), c = -1'})

subplot(3,2,2)
plot(x,ux1)
xlabel('x')
ylabel('u')
title({'Transport equation t = 10', 'u_0 = -sin(\pi x), c = -1'})

subplot(3,2,3)
imagesc(x,t,u2)
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('t')
title({'Transport Equation', 'u_0 = -sin(\pi x), c = 1'})

subplot(3,2,4)
plot(x,ux2)
xlabel('x')
ylabel('u')
title({'Transport equation t = 10', 'u_0 = -sin(\pi x), c = 1'})


subplot(3,2,5)
imagesc(x,t,u3)
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('t')
title({'Burgers Equation', 'u_0 = -sin(\pi x), c = u'})

subplot(3,2,6)
plot(x,ux3)
xlabel('x')
ylabel('u')
title({'Burgers equation t = 10', 'u_0 = -sin(\pi x), c = u'})



end

%% Problem 1c
function Problem1c
dx = .01;
dy = .01;

x = [-5:dx:5];
y = [0:dy:10];

nx = length(x);
ny = length(y);

u = zeros(ny,nx);
v = zeros(ny,nx);

u(1,1:nx) = sin(pi.*x);

% Leap frog method

u1 = u;
v1 = v;
Q = u + v;                  % Uses method of characteristics to get 2nd row
R = u - v;

R(2,1) = u1(1,nx-1) - v1(1,nx-1);
R(2,nx) = R(2,1);
Q(2,1) = u1(1,2) + v1(1,2);
Q(2,nx) = Q(2,1);
    for i = 2:nx-1
        R(2,i) = u1(1,i-1) - v1(1,i-1);
        Q(2,i) = u1(1,i+1) + v1(1,i+1);
    end
    for i = 1:nx
        u1(2,i) = (Q(2,i)+R(2,i))/2;
        v1(2,i) = (Q(2,i)-R(2,i))/2;
    end
        
for j=2:ny-1                    % Periodic BC at x=0
        u1(j+1,1) = ((dy/dx)*(v1(j,2)-v1(j,nx-1)))+u1(j-1,1);
        v1(j+1,1) = ((dy/dx)*(u1(j,2)-u1(j,nx-1)))+v1(j-1,1);
    for i = 2:nx-1
        u1(j+1,i) = ((dy/dx)*(v1(j,i+1)-v1(j,i-1)))+u1(j-1,i);
        v1(j+1,i) = ((dy/dx)*(u1(j,i+1)-u1(j,i-1)))+v1(j-1,i);
    end                         % Periodic BC at x = nx
        u1(j+1,nx) = ((dy/dx)*(v1(j,2)-v1(j,nx-1)))+u1(j-1,nx);
        v1(j+1,nx) = ((dy/dx)*(u1(j,2)-u1(j,nx-1)))+v1(j-1,nx);
end

    
% Method of characteristics
u2 = u; 
v2 = v;
Q = u + v;            
R = u - v;  
    
for j = 1:ny-1                              % Gets endpoints using BC
    R(j+1,1) = u2(j,nx-1) - v2(j,nx-1);
    R(j+1,nx) = R(j+1,1);
    Q(j+1,1) = u2(j,2) + v2(j,2);
    Q(j+1,nx) = Q(j+1,1);
    
    for i = 2:nx-1                          % Solves R,Q
        R(j+1,i) = u2(j,i-1) - v2(j,i-1);
        Q(j+1,i) = u2(j,i+1) + v2(j,i+1);
    end
    
    for i = 1:nx                            % Recovers u,v
        u2(j+1,i) = (Q(j+1,i)+R(j+1,i))/2;
        v2(j+1,i) = (Q(j+1,i)-R(j+1,i))/2;
    end
end


% Augmented system

u3 = u;
v3 = v;
epsilon = dy/2;

Q = u + v;                  % Uses method of characteristics to get 2nd row
R = u - v;

R(2,1) = u3(1,nx-1) - v3(1,nx-1);
R(2,nx) = R(2,1);
Q(2,1) = u3(1,2) + v3(1,2);
Q(2,nx) = Q(2,1);
    for i = 2:nx-1
        R(2,i) = u3(1,i-1) - v3(1,i-1);
        Q(2,i) = u3(1,i+1) + v3(1,i+1);
    end
    for i = 1:nx
        u3(2,i) = (Q(2,i)+R(2,i))/2;
        v3(2,i) = (Q(2,i)-R(2,i))/2;
    end

for j = 2:ny-1          % Periodic BC at x = 0
    u3(j+1,1) = (1/((epsilon/dy^2) + 1/(2*dy)))* (   ((v3(j,2)-v3(j,nx-1))/(2*dx)) + ((epsilon/dx^2)*(u3(j,2)-2*u3(j,1)+u3(j,nx-1))) + ((epsilon/dy^2)*((2*u3(j,1))-u3(j-1,1))) + u3(j-1,1)/(2*dy));
    v3(j+1,1) = (1/((epsilon/dy^2) + 1/(2*dy)))* (   ((u3(j,2)-u3(j,nx-1))/(2*dx)) + ((epsilon/dx^2)*(v3(j,2)-2*v3(j,1)+v3(j,nx-1))) + ((epsilon/dy^2)*((2*v3(j,1))-v3(j-1,1))) + v3(j-1,1)/(2*dy));
    for i = 2:nx-1
    u3(j+1,i) = (1/((epsilon/dy^2) + 1/(2*dy)))* (   ((v3(j,i+1)-v3(j,i-1))/(2*dx)) + ((epsilon/dx^2)*(u3(j,i+1)-2*u3(j,i)+u3(j,i-1))) + ((epsilon/dy^2)*((2*u3(j,i))-u3(j-1,i))) + u3(j-1,i)/(2*dy));
    v3(j+1,i) = (1/((epsilon/dy^2) + 1/(2*dy)))* (   ((u3(j,i+1)-u3(j,i-1))/(2*dx)) + ((epsilon/dx^2)*(v3(j,i+1)-2*v3(j,i)+v3(j,i-1))) + ((epsilon/dy^2)*((2*v3(j,i))-v3(j-1,i))) + v3(j-1,i)/(2*dy));
    end                 % Periodic BC at x = nx
    u3(j+1,nx) = (1/((epsilon/dy^2) + 1/(2*dy)))* (   ((v3(j,2)-v3(j,nx-1))/(2*dx)) + ((epsilon/dx^2)*(u3(j,2)-2*u3(j,nx)+u3(j,nx-1))) + ((epsilon/dy^2)*((2*u3(j,nx))-u3(j-1,nx))) + u3(j-1,nx)/(2*dy));
    v3(j+1,nx) = (1/((epsilon/dy^2) + 1/(2*dy)))* (   ((u3(j,2)-u3(j,nx-1))/(2*dx)) + ((epsilon/dx^2)*(v3(j,2)-2*v3(j,nx)+v3(j,nx-1))) + ((epsilon/dy^2)*((2*v3(j,nx))-v3(j-1,nx))) + v3(j-1,nx)/(2*dy));
end






figure(3)

[X,Y] = meshgrid(x,y);

subplot(1,3,1)
contourf(X,Y,u1)
axis square
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('y')
title('Leap Frog')

subplot(1,3,2)
contourf(X,Y,u2)
axis square
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('y')
title('Method of Characteristics')

subplot(1,3,3)
contourf(X,Y,u3)
axis square
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('y')
title('Augmented System')
    
end

%% Problem 2a
function Problem2a
dx = .1;
dy = .1;

x = [0:dx:pi];
y = [0:dy:5*pi];

nx = length(x);
ny = length(y);

u = zeros(ny,nx);
u(2:end-1,2:end-1) = 1;
u(1,:) = sin(x);
cor = 1;
tol = 1e-8;
w=1.4;

% Vertical line SOR

uNew = u;
while cor > tol
    for i = 2:nx-1                        % forms tridiagonal system at each column
        uOld = uNew;
        a = ones(ny-1,1)/(dy^2);
        a(end) = 0;
        b = ones(ny,1)*((-2/dy^2)-(2/dx^2));
        b(1) = 1;
        b(end) = 1;
        c = a;
        c(1) = 0;
        d = zeros(ny,1);                % enforces BC
        d(1) = u(1,i);
        d(end) = u(ny,i);
        for j = 2:ny-1
            d(j)= (-uOld(j,i+1)-uNew(j,i-1))/dx^2;
        end
        new = Thomas3(a,b,c,d);
        new = uOld(:,i) + w*(new' - uOld(:,i));             % SOR
        uNew(:,i)=new;
        correction(i) = max(abs(uNew(:,i) - uOld(:,i)));    % Cor at each column
    end
        cor = max(correction);
end
u1 = uNew;


% Horizontal Line SOR
dx = .1;
dy = .1;

x = [0:dx:pi];
y = [0:dy:5*pi];

nx = length(x);
ny = length(y);

u = zeros(ny,nx);
u(2:end-1,2:end-1) = 1;
u(1,:) = sin(x);
cor = 1;
tol = 1e-8;
w=1.2;


uNew = u;

correction = zeros(nx,1);
while cor > tol
    for j = 2:ny-1                        % forms tridiagonal system at each row
        uOld = uNew;
        a = ones(nx-1,1)/(dx^2);
        a(end) = 0;
        b = ones(nx,1)*((-2/dy^2)-(2/dx^2));
        b(1) = 1;
        b(end) = 1;
        c = a;
        c(1) = 0;
        d = zeros(nx,1);                % enforces BC
        d(1) = u(j,1);
        d(end) = u(j,nx);
        for i = 2:nx-1
            d(i)= (-uOld(j+1,i)-uNew(j-1,i))/dy^2;
        end
        new = Thomas3(a,b,c,d);
        new = uOld(j,:) + w*(new - uOld(j,:));             % SOR
        uNew(j,:)=new;
        correction(j) = max(abs(uNew(j,:) - uOld(j,:)));    % Cor at each row
    end
        cor = max(correction);
end
u2 = uNew;




figure(4)

subplot(1,2,1)

pcolor(x,y,u1)             
shading interp
set(gca,'YDir','normal')
colorbar
xlabel('x')
ylabel('y')
title('Vertical SOR')

subplot(1,2,2)
pcolor(x,y,u2)               
set(gca,'YDir','normal')
colorbar
shading interp
xlabel('x')
ylabel('y')
title('Horizontal SOR')

end

%% Problem 2b
function Problem2b()
dx = .025;
dy = .025;

x = [-1:dx:2];
y = [0:dy:2];

[u1,res1] = Subsonic(0);
[u2,res2] = Subsonic(0.4);

cp1 = 2*(1-u1(1,:));
cp2 = 2*(1-u2(1,:));

figure(5)

subplot(2,3,1)

contourf(x,y,u1)             
shading interp
set(gca,'YDir','normal')
colorbar
xlabel('x/c')
ylabel('y')
title('Subsonic Flow: M = 0')

subplot(2,3,2)
loglog(res1)
xlabel('Iterations')
ylabel('Correction')
title('Convergence: M = 0')

subplot(2,3,3)
plot(x,cp1)
xlabel('x/c')
ylabel('C_p')
title('Coefficient of Pressure: M = 0')

subplot(2,3,4)
contourf(x,y,u2)             
shading interp
set(gca,'YDir','normal')
colorbar
xlabel('x/c')
ylabel('y')
title('Subsonic Flow: M = 0.4')

subplot(2,3,5)
loglog(res2)
xlabel('Iterations')
ylabel('Correction')
title('Convergence: M = 0.4')

subplot(2,3,6)
plot(x,cp2)
xlabel('x/c')
ylabel('C_p')
title('Coefficient of Pressure: M = 0.4')


end

function [u1,residual] = Subsonic(M)
dx = .025;
dy = .025;

x = [-1:dx:2];
y = [0:dy:2];

nx = length(x);
ny = length(y);

u = ones(ny,nx);
cor = 1;
tol = 1e-5;
tau = .1;
w = 1.3;

uNew = u;
vVec = zeros(1,nx);                     % Forms v vector along x axis
yB = zeros(1,nx);
foilX = [0:dx:1];
foil = tau * sin(pi*foilX);             % airfoil geometry
yB(41:81) = foil(:);
iteration = 0;
for k = 2:nx-1                          % dyB/dx = v
    vVec(k) = (yB(k+1) - yB(k-1))/(2*dx);
end

while cor > tol
    iteration = iteration + 1;
    
    for i = 2:nx-1
        uOld = uNew;
        a = ones(ny-1,1)/(dy^2);
        b = ones(ny,1)*((-2/dy^2)-((2*(1-M^2))/dx^2));
        b(end) = 1;
        c = a;
        a(end) = 0;
        c(1) = 2/dy^2;
        d = zeros(ny,1);                % enforces BC
        d(1) = ((1-M^2)*(-uOld(1,i+1)-uNew(1,i-1))/dx^2)+(vVec(i+1)-vVec(i-1))/(dx*dy);
        d(end) = u(ny,i);
        for j = 2:ny-1
            d(j)= (1-M^2)*(-uOld(j,i+1)-uNew(j,i-1))/dx^2;
        end
        new = Thomas3(a,b,c,d);
        new = uOld(:,i) + w *(new' - uOld(:,i));             % SOR
        uNew(:,i)=new;
        correction(i) = max(abs(uNew(:,i) - uOld(:,i)));  % Cor at each column
    end
    
        cor = max(correction);
        residual(iteration) = cor;
end
u1 = uNew;

end

%% Tridiagonal Thomas Algorithm Function
function xBar = Thomas3(a,b,c,d)
    n = length(b);                      % Length of total diagonal
    aBar = [0;a];                       % Forms vectors of n length for manipulated values, padded as necessary
    bBar = b;
    cBar = [c;0];
    dBar = d;

    for i=2:n                           % "Zip-down" eliminates subdiagonal by subtracting the previous row scaled to the 'a' term
        bBar(i) = b(i) - aBar(i)*cBar(i-1)/bBar(i-1);
        dBar(i) = d(i) - aBar(i)*dBar(i-1)/bBar(i-1);
    end
    
    xBar(n) = dBar(n)/ bBar(n);         % "Zip-up" solves for x from last row single variable equation to first row
    for i = n-1:-1:1
        xBar(i) = ((dBar(i)-(cBar(i)*xBar(i+1)))/bBar(i));
    end
end

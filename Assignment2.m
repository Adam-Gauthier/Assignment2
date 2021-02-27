close all; clear all; clc;
%% Assignment 2 ELEC4700
% Adam Gauthier 100947233
% 
% This report summarizes and discusses the numerical results obtained from
% numerically computing the electrostatic potential of an area of specified
% conductivity.
% The technique for computing the results is finite difference.

%% Part 1 Solving Laplace's Equation
% In part one of this assignment, Laplace's equation was solved using
% finite difference. A rectangular area was setup with uniform(constant)
% conductivity an solved using the matrix equation G*F=V.
% 
% $$ Laplace's Eq: \nabla^2 V = 0 $$
% 
% Laplace equation can be turned into a difference equation of the following
% form:
% 
% $$\frac{V_{x-1,y}-2V_{x,y}+V_{x+1,y}}{(\Delta x)^2} +
% \frac{V_{x,y-1}-2V_{x,y}+V_{x,y+1}}{(\Delta y)^2}=0$$
% 
% For part one, the boundary conditions (BC) were setup as specified in the
% instructions, with one side having a constant voltage and the other 0.

% Initialization Part 1
W = 2; % width
L = 3; % length 
dx=0.1; % step x
dy=0.1; % step y 
nx=L/dx; % points in x
ny=W/dy; % points in y
V0=1; % BC voltage

%Terms to populate G from difference eq
Vxy = -2*(1/dx^2 +1/dy^2); % middle term
Vy = 1/(dy^2); % side term
Vx = 1/(dx^2); % side term
%%
% Question 1 a) solve Laplace's equation using GV=F

% Build G matrix
G_a = sparse(nx*ny,nx*ny);
for i = 2:(nx-1) 
    for j=2:(ny-1) 
         p=mapping_equation(i,j,nx);
         G_a(p,p)=Vxy;
         G_a(p,mapping_equation(i,j-1,nx))=Vy;
         G_a(p,mapping_equation(i,j+1,nx))=Vy; 
         G_a(p,mapping_equation(i-1,j,nx))=Vx;
         G_a(p,mapping_equation(i+1,j,nx))=Vx; 
    end
end

% Build F vector + BCs
F_a = sparse(nx*ny,1);
for j=1:ny
    p=mapping_equation(1,j,nx);
    F_a(p)=V0;
    G_a(p,p)=1; % BC
    p=mapping_equation(nx,j,nx);
    G_a(p,p)=1; % BC
end
for i=2:(nx-1)
    p=mapping_equation(i,1,nx);
    G_a(p,p)=1;
    G_a(p,mapping_equation(i,2,nx))=-1; % BC
    p=mapping_equation(i,ny,nx);
    G_a(p,p)=1;
    G_a(p,mapping_equation(i,ny-1,nx))=-1; % BC
end

result_a = G_a\F_a;
result_a = reshape(result_a,[],ny)';

figure(1);
[X,Y]=meshgrid(linspace(0,L,nx),linspace(0,L,ny));
surf(X,Y,result_a)
xlabel('Length');
ylabel('Width');
title('Laplace Equation Finite Difference');
view(45,45);
colorbar
%%
% The voltage for the rectangular geometry was solved using Laplace's
% equation and finite difference. The above figure shows the result in 3D.
% As we can see, the voltage is linear from one side to the other as
% expected since one side has potential and the other is ground and the
% conductivity is constant throughout.
%
% Question 1 b) Numerical vs Analytic Solution
% 
% In part b) the boundary conditions were changed so that curvature was
% introduced. The voltage was solved the same as before. The numerical
% result was compared to an analytic solution which is in the form of a
% series. 

% build new G and F matrix/vector with new BCs
G_b = sparse(nx*ny,nx*ny);
for i = 2:(nx-1) 
    for j=2:(ny-1) 
         p=mapping_equation(i,j,nx);
         G_b(p,p)=Vxy;
         G_b(p,mapping_equation(i-1,j,nx))=Vx;
         G_b(p,mapping_equation(i+1,j,nx))=Vx;
         G_b(p,mapping_equation(i,j-1,nx))=Vy;
         G_b(p,mapping_equation(i,j+1,nx))=Vy;
    end
end 

F_b=sparse(nx*ny,1);
for i=1:nx
    p = mapping_equation(i,1,nx);
    G_b(p,p) = 1;
    p = mapping_equation(i,ny,nx);
    G_b(p,p) = 1;
end
for j=1:ny
    p = mapping_equation(1,j,nx);
    F_b(p) = V0;
    G_b(p,p) = 1;
    p = mapping_equation(nx,j,nx);
    G_b(p,p) = 1;
    F_b(p) = V0;
end
F_b(mapping_equation(1,1,nx)) = 0;
F_b(mapping_equation(1,ny,nx)) = 0;
F_b(mapping_equation(nx,1,nx)) = 0;
F_b(mapping_equation(nx,ny,nx)) = 0;

result_b = G_b\F_b;
result_b = reshape(result_b,[],ny)'; 

figure(2);
[X,Y]=meshgrid(linspace(0,L,nx),linspace(0,L,ny));
surf(X,Y,result_b)
xlabel('Length');
ylabel('Width');
title('Numerical Solution');
colorbar

% analytic solution using series
analytic=zeros(ny,nx);
x=repmat(linspace(-L/2,L/2,nx),ny,1);
y=repmat(linspace(0,W,ny),nx,1)';
for term = 1:2:250
    analytic=analytic + (1./term).*cosh(term.*pi.*x./W).*sin(term.*pi.*y./W)./cosh(term.*pi.*(L./2)./W);
end
analytic=analytic.*4.*V0./pi;

figure(3);
surf(linspace(0,L,nx),linspace(0,W,ny),analytic);
xlabel('Length');
ylabel('Width');
title('Analytic Solution')
colorbar
%%
% The two figures above show the numerically computed solution using
% finite difference and the analytic solution computed using a series
% equation. As we can see, the solutions are the same or very close. The
% curvature is present in both figures from the boundary conditions. From
% subtracting one solution from the other, we found that the diffence was
% very small. As more terms are computed in the sum for the analytic
% solution, it converges to the true value and shows less error with
% respect to the numerical solution. The same would be expected by using
% finer meshing in the numerical solution.
%
% Both approaches work well and provide valid results with some error. The
% analytical approach may be preferred as it requires less computational
% power and gives insight into the result from it's equation. The analytical
% solution series can be truncated once the terms you are adding are below
% some threshold (small enough). For the numerical solution, it is harder to
% gain precision because adding mesh density greatly increases computation
% time. One benefit of the numerical solution is that it is flexible and
% does not require you to understand or solve the problem
% theoretically/analytically. 

%% Part 2 Added Conductivities and Current
% In part 2, a varying conductivity was introduced to change Laplace's
% equation. Two boxes which 'bottle neck' the flow from one side to the
% other were added into the rectangular space. The two boxes were given a
% low conductivity as specified in the manual. The remaining region had a
% conductivity of 1. Parameters like the box heights, the mesh density and
% the conductivity of the boxes were changed while recording the current.
% The plots of these parameters versus the current is shown below.
% 
% New equation: $$\nabla(\sigma_{xy} \nabla V) = 0 $$
%
% In part 2, a general function "compute_Q2.m" was used to compute the
% voltage result from the finite differences equations as shown in part 1.


clear all;clc;
nx = 200; %x points
ny = 100; %y points
cond = 1e-2; % low condunctivity in boxes

%compute results from general Q2 problem function
[cond_map,voltage,Ex,Ey,Jx,Jy,current] = compute_Q2(nx,ny,cond,1);

%%
% Part 2 a) plots of conductivity, voltage, field and current.

figure(4) % conductivity plot
surf(cond_map,'edgecolor','none')
colorbar
xlabel('x')
ylabel('y')
title('Conductivity Plot')
view(0,90)

%%
% The plot above shows the conductivity after adding our two low
% conductivity boxes. As we can see the conductivity in these box regions
% is much less (1/100) than the surrounding region. We expect these regions
% to have very little current flow as the lower resistance path is to go
% through the higher conductivity region.

figure(5)
contourf(voltage,50)
xlabel('x')
ylabel('y')
title('Voltage Plot')
colorbar
%%
% The plot above shows the result for the voltage compute using finite
% difference from the modified Laplace equation which accounts for the
% conductivities of the added boxes. The plot is a contour plot which shows
% the lines of equal potential. As we can see the potential appears
% somewhat linear on the edges and deforms as it approaches the middle
% where the conductivity changes. 

figure(6)
quiver(Ex,Ey)
ylim([0,100])
title('Electric Field Plot')
xlabel('x')
ylabel('y')
%%
% The gradient function was used and negated to get the electric field from
% the voltage. The plot aboce shows the electric field vectors in both
% dimensions. As we can see, the field is fairly straight except in the
% middle where it is distorted by the conductivity changes. The field is
% strongest inside the low conductivity zone due to the relationship
% between field, voltage and conductivity. The field between the boxes is
% slightly stronger than the surrounding area due to the proximity to the
% boxes.

figure(7)
quiver(Jx,Jy)
ylim([0,ny])
xlim([0,nx])
title('Current Plot')
xlabel('x')
ylabel('y')

%%
% The figure above shows the current map of the rectangular area. It is
% complementary to the previous electric field plot as the current is
% calculated from the electric field and conductivity. As we can see, there
% is very little current going through the low conductivity regions as
% expected. Most of the current is going through the low resistance center
% path between the boxes. Therefore, the boxes act as a bottle neck to
% funnel the current through an area, thereby increasing the resistance.

fprintf('Current is calculated as %d A\n',current)

%%
% From the analysis above, the parameters were changed in order to see the
% effect on the current. The mesh density, bottle neck size, and
% conductivity of the boxes was change and the current calculated in order
% to see the effect of each parameter on the current. 
%
% Question 2 b) current vs mesh size.
index=1;
multipliers =[0.25,0.5,0.75,1];
for m = multipliers
    [cond_map,voltage,Ex,Ey,Jx,Jy,current] = compute_Q2(nx*m,ny*m,cond,1);
    mesh_array(index) = current;
    index=index+1;
end
figure(8)
hold on;
plot(multipliers,mesh_array)
plot(multipliers,mesh_array,'o')
title('Current VS Mesh Density')
xlabel('Mesh Density Multiplier')
ylabel('Current')
hold off;
%%
% The plot above shows the result of changing the mesh size by some
% multiplier. The multiplier was used to scale the mesh size below the
% current value and see the effecton the current. As we can see in the
% plot, as the mesh size increases the value increases asymptotically and
% approaches the true value. At low mesh densities the mesh is too coarse
% to get good precision but increasing the mesh density increases the
% precision and better represents the true value.
% 
%  Question 2 c) current vs bottle neck
narrow = [5:5:40];
index=1;
for k = narrow
    [cond_map,voltage,Ex,Ey,Jx,Jy,current] = compute_Q2(nx,ny,cond,k);
    narrow_array(index) = current;
    index=index+1;
end
figure(9)
hold on;
plot(narrow,narrow_array)
plot(narrow,narrow_array,'o')
title('Current VS Bottle Neck')
xlabel('Box Height (Narrowing Bottle Neck)')
ylabel('Current')
hold off;
%%
% The plot above shows the result of current as you increase the heigh of
% the boxes, thereby reduce the bottle neck gap and allowing less current
% through the high conductivity region. As expected, the current drops off
% as the bottle neck increases due to the low conductivity regions being
% larger and not allowing current to pass. As the bottle neck begins to
% close we expect a large drop in current since there will be no high
% conductivity path for the current to take through the area. 
%
% Question 2 c) current vs conductivity
cond_mult = [0.5*10:0.25*10:2.5*10];
index=1;
for q = cond_mult
    [cond_map,voltage,Ex,Ey,Jx,Jy,current] = compute_Q2(nx,ny,cond*q,1);
    cond_array(index) = current;
    index=index+1;
end
figure(10)
hold on;
plot(cond_mult,cond_array)
plot(cond_mult,cond_array,'o')
title('Current VS Conductivity')
xlabel('Conductivity Multiplier')
ylabel('Current')
hold off;
%%
% Lastly, the conductivity in the boxes was changed by multiplying the
% current conductivity by a factor that scaled the value. As expected, as
% the conductivity went up in the boxes the overall current increases
% somewhat linearly. This is expected as there is a linear relationship
% between current and conductivity. 
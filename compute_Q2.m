% Function to compute results for Q2

function [cond_map,voltage, Ex,Ey,Jx,Jy,current] = compute_Q2(nx,ny,cond,narrowing)
    cond_map = ones(ny,nx); % all initialize to 1 at start
    if narrowing~=1
        box(1,1)=nx/2;
        box(1,2)=round(narrowing/2);
        box(1,3)=narrowing+1;
        box(1,4)=0.2*nx;
        box(2,1)=nx/2;
        box(2,2)=100-round(narrowing/2);
        box(2,3)=narrowing+1;
        box(2,4)=0.2*nx;
    else
        box(1,1)=nx/2;
        box(1,2)=0.1*nx;
        box(1,3)=0.2*nx+1;
        box(1,4)=0.2*nx;
        box(2,1)=nx/2;
        box(2,2)=0.4*nx;
        box(2,3)=0.2*nx+1;
        box(2,4)=0.2*nx;
    end
    
    % add conductivity from boxes
    cond_map=add_boxes(cond_map,box(1,:),cond,nx,ny);
    cond_map=add_boxes(cond_map,box(2,:),cond,nx,ny);
    
    %build G matrix for 
    G_2 = sparse(nx*ny);
    F_2 = sparse(1,nx*ny);
    % Populate G according to example given in class and BCs
    for i = 1:ny 
        for j = 1:nx 
            n = i + (j - 1) * ny; % mapping eq
            if j==1 
                G_2(n,:) = 0;
                G_2(n,n)=1;
                F_2(n)=1;
            elseif j==nx
                G_2(n,:)=0;
                G_2(n,n)=1;
            elseif i==1 
                xm = i + (j - 2) * ny;
                xp = i + (j) * ny;
                yp = i + 1 + (j - 1) * ny;
                x_minus = (cond_map(i,j)+cond_map(i,j-1))/2;
                x_plus = (cond_map(i,j)+cond_map(i,j+1))/2;
                y_plus = (cond_map(i,j)+cond_map(i+1,j))/2;
                G_2(n, n) = -(x_minus+x_plus+y_plus);
                G_2(n, xm) = x_minus;
                G_2(n, xp) = x_plus;
                G_2(n, yp) = y_plus;
            elseif i==ny
                xm = i + (j - 2) * ny;
                xp = i + (j) * ny;
                ym = i - 1 + (j - 1) * ny;
                x_minus = (cond_map(i, j) + cond_map(i,j-1))/2;
                x_plus = (cond_map(i, j) + cond_map(i,j+1))/2;
                y_minus = (cond_map(i, j) + cond_map(i-1,j))/2;
                G_2(n,n) = -(x_minus+x_plus+y_minus);
                G_2(n,xm) = x_minus;
                G_2(n,xp) = x_plus;
                G_2(n,ym) = y_minus;
            else
                xm = i + (j-2)*ny;
                xp = i + (j)*ny;
                ym = i-1 + (j-1)*ny;
                yp = i+1 + (j-1)*ny;

                x_minus = (cond_map(i, j) + cond_map(i , j-1))/2;
                x_plus = (cond_map(i, j) + cond_map(i , j+1))/2;
                y_minus = (cond_map(i, j) + cond_map(i-1, j ))/2;
                y_plus = (cond_map(i, j) + cond_map(i+1, j))/2;

                G_2(n,n) = -(x_minus+x_plus+y_minus+y_plus);
                G_2(n,xm) = x_minus;
                G_2(n,xp) = x_plus;
                G_2(n,ym) = y_minus;
                G_2(n,yp) = y_plus;
            end
        end
    end

    result_2 = G_2\F_2';

    voltage = zeros(ny,nx);
    for i = 1:ny
        for j = 1:nx
            n = i + (j - 1) * ny;
            voltage(i, j) = result_2(n);
        end
    end
    [Ex,Ey]=gradient(voltage); % Calculates the gradient in one line
    Ex=-Ex; 
    Ey=-Ey;
    Jx=Ex.*cond_map;
    Jy=Ey.*cond_map;
    current = sum(Jx(1:ny,1));


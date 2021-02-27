function [matrix] = add_boxes(matrix,boxes,cond,nx,ny)
    for j=1:nx
        for i=1:ny     
            if(abs(i-boxes(1,2)) < boxes(1,3)/2 && abs(j-boxes(1,1)) < boxes(1,4)/2)  
                matrix(i,j)=cond;   % change cond if insie box
            end
        end
    end
end
function [xcoords,ycoords]=nodesmax(im,threshold,s)

%     [ycoords,xcoords]=f_dynadetection(handles.var.img,threshold);
   
    padDim = max([3*round(threshold),s+2]);
    imgPadded = padarray(im,[padDim,padDim],'replicate');
    if threshold==0
        pstruct.x=length(imgPadded(:,1))/2;
        pstruct.y=length(imgPadded(:,1))/2;
    else 
        [pstruct, ~, ~, ~] = pointSourceDetection(imgPadded, threshold, 'mode','xy');
    end
    xcoords=[];
    ycoords=[];
    if ~isempty(pstruct)
%         if (length(pstruct.x)>20)
            id=1;
            for i=1:length(pstruct.x)
               %if (round(pstruct.y(i)-s)>0)&&(round(pstruct.y(i)+s)<length(im(:,1,:)))&&(round(pstruct.x(i)-s)>0)&&(round(pstruct.x(i)+s)<length(im(1,:,:)))
               
               A=imgPadded(round(pstruct.y(i))-s:round(pstruct.y(i))+s,round(pstruct.x(i))-s:round(pstruct.x(i))+s);
               [~,I]=max(A(:));
               [I_row, I_col] = ind2sub(size(A),I);
               xcurrent=I_col+pstruct.x(i)-s-1-padDim;
               ycurrent=I_row+pstruct.y(i)-s-1-padDim;
               if xcurrent>0&&ycurrent>0
                    xcoords(id)=xcurrent;
                    ycoords(id)=ycurrent;
                    id=id+1;
               end
            end
%         else
%             xcoords=[];
%             ycoords=[];
%         end
    end
    if isempty(xcoords)
               pstruct.x=length(imgPadded(:,1))/2;
               pstruct.y=length(imgPadded(:,1))/2;
               xcoords = pstruct.x' - padDim;
               ycoords = pstruct.y' - padDim;
               A=imgPadded(round(pstruct.y)-s:round(pstruct.y)+s,round(pstruct.x)-s:round(pstruct.x)+s);
               [~,I]=max(A(:));
               [I_row, I_col] = ind2sub(size(A),I);
               xcoords=I_col+xcoords-s-1;
               ycoords=I_row+ycoords-s-1;
    end
end
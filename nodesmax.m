function [xcoords,ycoords]=nodesmax(im,threshold,s)

%     [ycoords,xcoords]=f_dynadetection(handles.var.img,threshold);
    padDim = 3*round(threshold);
    imgPadded = padarray(im,[padDim,padDim],'replicate');
    [pstruct, ~, ~, ~] = pointSourceDetection(imgPadded, threshold, 'mode','xy');
    if ~isempty(pstruct)
%         if (length(pstruct.x)>20)
            xcoords = pstruct.x' - padDim;
            ycoords = pstruct.y' - padDim;
            for i=1:length(xcoords)
               A=imgPadded(round(pstruct.y(i))-s:round(pstruct.y(i))+s,round(pstruct.x(i))-s:round(pstruct.x(i))+s);
               [~,I]=max(A(:));
               [I_row, I_col] = ind2sub(size(A),I);
               xcoords(i)=I_col+xcoords(i)-s-1;
               ycoords(i)=I_row+ycoords(i)-s-1;
            end
%         else
%             xcoords=[];
%             ycoords=[];
%         end
    else
        xcoords=[];
        ycoords=[];
    end
end
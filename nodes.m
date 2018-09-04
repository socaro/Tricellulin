function [xcoords,ycoords]=nodes(im,threshold)

%     [ycoords,xcoords]=f_dynadetection(handles.var.img,threshold);
    padDim = 3*round(threshold);
    imgPadded = padarray(im,[padDim,padDim],'replicate');
    [pstruct, ~, ~, ~] = pointSourceDetection(imgPadded, threshold, 'mode','xy');
    if ~isempty(pstruct)
%         if (length(pstruct.x)>20)
            xcoords = pstruct.x' - padDim;
            ycoords = pstruct.y' - padDim;
%         else
%             xcoords=[];
%             ycoords=[];
%         end
    else
        xcoords=[];
        ycoords=[];
    end
end
function rarray = f_drawrectangles(data, maxval)
% data already ordered 4x4 mapping to electrode array
    pos = [70, 510, 960, 1370];
    rarray = gobjects(4,4);
    for i = 1:4
        for j = 1:4
            red = max(1-data(i,j)/maxval, 0);
            green = max(1-data(i,j)/maxval,0);
            rarray(i,j) = rectangle('Position',[pos(j),pos(i),90,90], 'FaceColor', [red, green, 1], 'EdgeColor', 'none');
        end
    end
    
end
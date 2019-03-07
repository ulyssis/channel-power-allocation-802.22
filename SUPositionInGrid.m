% Creat a grid, and randomly locate the SUs on the grid points. The grid is maximal 10 x 10 area
% without location schufful, the suqence of user 1 to user n is up-right.
function [pos, posET] = SUPositionInGrid(n, nET, lengthSide, coverage)

    for sth = 1: 10
       if(n >= sth^2 && n <= (sth+1)^2)
           scale1 = sth+1;
           break;
       end
    end
    scale = scale1^2;           % scale is minimum x^2 which is bigger than n

    pos = zeros(2, scale);      % (x, y) coordinates
    for i = 1:scale1
       for j = 1:scale1
           pos(1, (i-1)*scale1 + j) = 0.5*lengthSide/scale1 + (i-1)*lengthSide/scale1;
           pos(2, (i-1)*scale1 + j) = 0.5*lengthSide/scale1 + (j-1)*lengthSide/scale1;

%              pos(1, (i-1)*scale1 + j) = lengthSide/(scale1+1) + (i-1)*lengthSide/(scale1+1);        % concentrate a little bit
%              pos(2, (i-1)*scale1 + j) = lengthSide/(scale1+1) + (j-1)*lengthSide/(scale1+1);
       end
    end
    
% when n ~= scale1^2, the following section is executed
% choose scale1^2 - n columns to delete.
    for i=1:scale - n
        column2delete = floor(rand(1)*size(pos,2))+1;
        pos(:, column2delete) =[];
    end


% random distribution of ETs within each cell.
posET = [];
for i=1:n
    x1 = pos(1,i);
    y1 = pos(2,i);
    posETs = [];
    for t = 1: nET % loop until doing nET points inside the circle
        [x, y] = cirrdnPJ(x1,y1, coverage);
    %     plot(x,y,'x');
        posETs = [posETs, [x;y]];
    end % generate ETs for one WBS cell
    
    posET = [posET, posETs];
    
end

% % random distribution of ETs in the whole simulation area
% posET=zeros(2, n*nET);
% for i= 1: n*nET
%     posET(1, i) = rand*lengthSide;
%     posET(2, i) = rand*lengthSide;
% end
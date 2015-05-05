function [posET] = ETlocation(n, nET, lengthSide, coverage)

    for sth = 1: 10
       if(n > sth^2 && n <= (sth+1)^2)
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
%     shuffledpPos = pos(:,randperm(size(pos, 2)));
%     pos = shuffledpPos(:, 1:n);

posET = [];
for i=1:n
    x1 = pos(1,i);
    y1 = pos(2,i);
    posETs = [];
    for t = 1: nET % loop until doing nET points inside the circle
        [x y] = cirrdnPJ(x1,y1, coverage);
    %     plot(x,y,'x');
        posETs = [posETs, [x;y]];
    end % generate ETs for one WBS cell
    
    posET = [posET, posETs];
    
end
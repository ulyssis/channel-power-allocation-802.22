function [pos] = SUPositionRandom(n, lengthSide)


    pos = zeros(2, n);
    for i = 1:n
        pos(1, i) = rand * lengthSide;
        pos(2, i) = 2 + rand * lengthSide * 0.8;
    end
    pos = pos(:,randperm(size(pos, 2)));

function [pos] = PUPositionRandom(m, lengthSide)

    pos = zeros(2, m);

    fourCornors = [0, lengthSide, 0, lengthSide; 
                   0, 0, lengthSide, lengthSide]; 
    fourCornors = fourCornors(:,randperm(4));
    
    for i = 1: m 
        pos(:, i) = fourCornors (:, i);
    end
    
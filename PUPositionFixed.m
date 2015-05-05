function [pos] = PUPositionFixed(m, lengthSide)

    pos = zeros(2, m);

    fourCornors = [lengthSide/2, lengthSide+100, lengthSide/2, -100, -100, -150, -110; 
                   -100, lengthSide/2,           lengthSide+100, lengthSide/2, lengthSide+100, lengthSide+150, -50]; 
%    fourCornors = fourCornors(:,randperm(4));
    
    for i = 1: m 
        pos(:, i) = fourCornors (:, i);
    end
    
function [pos] = PUContorPosition(m, lengthSide)

    pos = zeros(2, m);
    
%     fourCornors = [lengthSide/2,    lengthSide*1.13,    lengthSide/2,      -lengthSide*0.3; 
%                    0,               lengthSide/2,   lengthSide*1.2,   lengthSide/2];

% the 4 contours randomly locate within one ring around the SUs
%     for i = 1: m 
%         if (rand > 0.5) % upper/down side
%             pos(1, i) = rand*lengthSide;
%             if(rand > 0.5)
%             pos(2, i) = lengthSide+rand*lengthSide/2;
%             else
%             pos(2, i) = -rand*lengthSide/2;
%             end
%         else
%             pos(2, i) = rand*lengthSide;
%             if (rand > 0.5)
%                 pos(1, i) = lengthSide+rand*lengthSide/2;
%             else
%                 pos(1, i) = -rand*lengthSide/2;
%             end
%         end
%         
%     end
    

% m contours randomly locate within one ring around the SUs
    seqq=randperm(m);
    for i=1:m
        
        if (seqq(i)==1)
            pos(1, i) = rand*lengthSide;
            pos(2, i) = lengthSide + rand*lengthSide/2;
        end
        
        if (seqq(i)==2)            
            pos(1, i) = lengthSide + rand*lengthSide/2;
            pos(2, i) = rand*lengthSide;  
        end
        
        if (seqq(i)==3)
            pos(1, i) = rand*lengthSide;
            pos(2, i) = -rand*lengthSide/2;            
        end
        
        if (seqq(i)==4)
            pos(1, i) = -rand*lengthSide/2;
            pos(2, i) = rand*lengthSide;               
        end
        
        if (seqq(i)==5)
            if (rand > 0.5) % upper/down side
                pos(1, i) = rand*lengthSide;
                if(rand > 0.5)
                pos(2, i) = lengthSide+rand*lengthSide/2;
                else
                pos(2, i) = -rand*lengthSide/2;
                end
            else
                pos(2, i) = rand*lengthSide;
                if (rand > 0.5)
                    pos(1, i) = lengthSide+rand*lengthSide/2;
                else
                    pos(1, i) = -rand*lengthSide/2;
                end
            end             
        end

        if (seqq(i)==6)
            if (rand > 0.5) % upper/down side
                pos(1, i) = rand*lengthSide;
                if(rand > 0.5)
                pos(2, i) = lengthSide+rand*lengthSide/2;
                else
                pos(2, i) = -rand*lengthSide/2;
                end
            else
                pos(2, i) = rand*lengthSide;
                if (rand > 0.5)
                    pos(1, i) = lengthSide+rand*lengthSide/2;
                else
                    pos(1, i) = -rand*lengthSide/2;
                end
            end               
        end    
        
        if (seqq(i)==7)
            if (rand > 0.5) % upper/down side
                pos(1, i) = rand*lengthSide;
                if(rand > 0.5)
                pos(2, i) = lengthSide+rand*lengthSide/2;
                else
                pos(2, i) = -rand*lengthSide/2;
                end
            else
                pos(2, i) = rand*lengthSide;
                if (rand > 0.5)
                    pos(1, i) = lengthSide+rand*lengthSide/2;
                else
                    pos(1, i) = -rand*lengthSide/2;
                end
            end               
        end          
        
    end
               
    
%     fourCornors = [lengthSide*1.1,        lengthSide*1.12,     lengthSide*1.16,       lengthSide*1.2; 
%                    lengthSide/2,        lengthSide/2,       lengthSide/2,       lengthSide/2];
    
% % performance is not so good:
%     fourCornors = [lengthSide/2,    lengthSide*1.1,    lengthSide/2,      -lengthSide*0.3; 
%                    -lengthSide*0.06,               lengthSide/2,   lengthSide*1.2,   lengthSide/2];

%     fourCornors = [lengthSide/2,    lengthSide*1.08,    lengthSide/2,      -lengthSide*0.09; 
%                    -lengthSide*0.06,               lengthSide/2,   lengthSide*1.07,   lengthSide/2];
               
%    fourCornors = fourCornors(:,randperm(4));
    
%     % same place:
%     fourCornors = [lengthSide/2,    lengthSide/2,     lengthSide/2,       lengthSide/2; 
%                    -0.1*lengthSide, -0.1*lengthSide,       -0.1*lengthSide,     -0.1*lengthSide];
    
%     for i = 1: m 
%         pos(:, i) = fourCornors (:, i);
%     end
    
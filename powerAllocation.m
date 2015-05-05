% adjust power after choosing channel!
% newUtility is proposed to adjust transmission power.
% max: E*R./p .* (1-exp(-0.5.*W*h*p/(R*f))).^L, from '98 paper

function  [B_powerAllocation] = powerAllocation(B, n, m, Gtilde, GtildeAll, TVpower, delta, SUcellRadius, pathlossfactor)

    F = (B* B' ~= 0);   % F illustrates the interferce relations
    F = F - eye(n);

    interfTV = sum(GtildeAll(n+m+1:n+m+m, 1:n)'.* (B~=0) *TVpower, 2);  % interference on every users

    InterferenceonAll = Gtilde.* F * sum(B, 2) + interfTV + delta;
 
    
    listQuasiSINR = InterferenceonAll./ (sum(B, 2)*SUcellRadius^(-pathlossfactor));    % f/p
    
%     [maxQuasiSINR, nodeIndexToStartFirst] =  max (listQuasiSINR);
    
    % start reduce power
    
    % maximize new utility: u = E*R./p .* (1-exp(-0.5.*gamma)).^L
    % when transmission power is not smaller than that has been decided,
    % and the quasiSINR is not smaller than 80% of the previous value.
    
    % discount of quasiSINR
% %     % discount = 0.8;
    

    % W: spectrum bandwidth
    % R: rate, bits/sec
    % f: interferece 
    % h: path gain
    % E: engerny
    % L: length of packet
    %

    W = 6e+6;
    R = 1e+6;
    f = 1e-6; %% renew
    h = 1e-6; %% renew
    E=1;
    L=80;
    %t = 1;

    p_max_list = sum(B, 2);
    convergence = 0;
    preB = B;
    
    while(convergence ==0)
    
        % update from 1 to N
        for i = 1:n

        %     % maximal power
        %     p_max = sum(B, 2)(nodeIndexToStartFirst);

            f = InterferenceonAll(n);
            h = SUcellRadius^(-pathlossfactor);

            % utility function
            func =@(p)(-E*R./p .* (1-exp(-0.5.*W*h*p/(R*f))).^L);
            %func =@(p)(-h*p/(f));
            
            p = fminbnd(func, 0, p_max_list(i));  

            % update matrix B
            index_power = find(B(i, :));
            B(i, index_power) = p;
            
            if(i ==n)
                if(preB == B)
                    convergence =1;
                else
                    preB = B;
                end
            end
        end

    end

B_powerAllocation = B;


    

        
        



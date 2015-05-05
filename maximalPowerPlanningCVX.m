    function [P, max_cvx_statusMsg] = maximalPowerPlanningCVX(n, m, infBound, GtildeAll, miniP, maxP)
        c = m;
        P = zeros(n*c, c);
        cvx_statusMsg = zeros(1, m);
        
        for channel = 1: m
            %------------ cvx start -----------
            cvx_precision default 
            cvx_begin gp %quiet % no output on the screa
            cvx_problem

            variables pMax(n)

            % objective function is convec
            minimize( -sum(log(pMax)) )  
            subject to

            miniP*ones(n,1)<=pMax<=maxP*ones(n,1); % 2n constraints 
            GtildeAll(n+channel, 1:n)*pMax<=infBound;          % 1 constraint

            cvx_end
            %------------ cvx end -----------

            cvx_statusMsg(channel) = cvx_optval;
            if (cvx_optval < 0)
                for i = 1: n
                P((i-1)*c + channel, channel) = pMax(i);
                end
            end
        end
        max_cvx_statusMsg = max(cvx_statusMsg);
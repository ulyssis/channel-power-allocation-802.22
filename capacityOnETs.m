function [averageShannonCPerCell] = capacityOnETs(B, n, w, GtildeETsSUs, nET)
    
    % n*nET x n*w, the power received from every SU
    rcvdPower_ET_SU = repmat(GtildeETsSUs(1:n*nET, n*nET+1:n*nET+n), 1, w).*repmat((ones(n*nET,1)*sum(B, 2)'), 1,w); 
    
    % for each cell, easy to decide on the transmitter, not easy to find
    % interfering SUs
    
    T = zeros(n*nET, w);
    
    B= repmat(B, w, 1);
    F = (B* B' ~= 0); % n*w x n*w, interfering relationship between SUs
    F = F - eye(size(B, 1));
    
    F_IntFe = zeros(n*nET, n*w);
    
    % power of good signals received
    p_signal = zeros(n*nET, 1);
    for j = 1: w
        for i = 1: n*nET
            F_IntFe(i,:) = F(((j-1)*n + ceil(i/nET)), :);
            p_signal(i,:) = rcvdPower_ET_SU(i, (j-1)*n + ceil(i/nET));
        end
        
        % bad interference caused by one set of co-address WBSs
        %n*nET x n*w
        bad_inf_ET_SU = rcvdPower_ET_SU.* F_IntFe; 
        Interference = sum(bad_inf_ET_SU, 2);

        % Shannon capacity
        T(:, j) = 6000*10*log10(Interference./p_signal);
    end
    
    % The aggregated Shannon capacity on the ETs.
    SC = sum(T, 2);
    
    averageShannonCPerCell = zeros(1, n);
    for i = 1:n
        averageShannonCPerCell(i) = sum(SC((i-1)*n+1:(i-1)*n+n));
    end
    
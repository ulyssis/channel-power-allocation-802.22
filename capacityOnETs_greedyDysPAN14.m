function [averageShannonCPerCell] = capacityOnETs_greedyDysPAN14(B, n, w, GtildeETsSUs, nET, delta)
    
    c = size(B, 2);
    w = c;
    % n*nET x n*w, the power received from every SU
    rcvdPower_ET_SU = repmat(GtildeETsSUs(1:n*nET, n*nET+1:n*nET+n), 1, w).*repmat((ones(n*nET,1)*sum(B, 2)'), 1,w); 
    
    % for each cell, easy to decide on the transmitter, not easy to find
    % interfering SUs
    
    T = zeros(n*nET, w);
    if(nnz(B) < size(B, 1)*w)
       originalBAbnormal =1; 
    end
    
    disp(B);
    
    %B = expandB(B, w);
    expandedB = expandB_greedyDysPAN14(B);
    
    disp(expandedB);
        
    %B= repmat(B, w, 1);
    F = (expandedB* expandedB' ~= 0); % n*w x n*w, interfering relationship between SUs
    F = F - eye(size(expandedB, 1));
    

    
    % the interfering WBSs to each ET.
    F_IntFe = zeros(n*nET, n*w);
    
    p_signal = zeros(n*nET, n*w);
    
    for i = 1: n*nET
        
        % power level of interferences
        for j = 1:w
            %--debug
            disp(size(F((j-1)*n + ceil(i/nET), :), 2));
            disp(size(rcvdPower_ET_SU(i, :), 2));
            if(size(F((j-1)*n + ceil(i/nET), :), 2) ~= size(rcvdPower_ET_SU(i, :), 2))
                stop =1;
                disp(w);
            end
            % --debug
            temp = F((j-1)*n + ceil(i/nET), :) .* rcvdPower_ET_SU(i, :);
            F_IntFe(i,:) = F_IntFe(i,:) + temp;
        end

        % find out the serving WBSs:
        indexWBS = ceil(i/nET);
        
        % power level of the received signals
        indexServingWBSsBase = zeros(1, n*w);
        for j = 1:w
            indexServingWBSsBase((j-1)*n+indexWBS) =1;
        end
        p_signal(i,:) = rcvdPower_ET_SU(i, :).*indexServingWBSsBase;

    end
        
    % interference on each channel on an ET
    B_workingChannelBinary = (expandedB~=0);
    interferenceOnETs = F_IntFe * B_workingChannelBinary; % n*nET X c
    % signal level on each channel on an ET
    signalOnETs = p_signal * B_workingChannelBinary;      % n*nET X c
    

    T = 6000*10*log2(1+signalOnETs./(interferenceOnETs + delta));

    % ---debug, existense of Inf
    T2 = isinf(T);
    if(nnz(T2) > 0)
        stop =1;
    end
    % ---debug, existense of Inf
    
%     % -------
%     for j = 1: w
%         for i = 1: n*nET
%             F_IntFe(i,:) = F_IntFe(i,:)+ F(((j-1)*n + ceil(i/nET)), :);
%             
%             % power of good signals received
%             p_signal(i,:) = rcvdPower_ET_SU(i, (j-1)*n + ceil(i/nET));
%         end
%         
%         % bad interference caused by one set of co-address WBSs
%         %n*nET x n*w
%         bad_inf_ET_SU = rcvdPower_ET_SU.* F_IntFe; 
%         Interference = sum(bad_inf_ET_SU, 2);
% 
%         % Shannon capacity
%         T(:, j) = 6000*10*log2(1+Interference./p_signal);
%     end
%     %-------
    % The aggregated Shannon capacity on the ETs.
    T(isnan(T))=0;
    SC = sum(T, 2);
    
    % average Shannon capacity per cell:
    averageShannonCPerCell = zeros(1, n);
    for i = 1:n
        averageShannonCPerCell(i) = mean(SC((i-1)*nET+1 : (i-1)*nET+nET));
    end
    
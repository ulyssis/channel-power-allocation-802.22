function  [SINR_ETs] = SINR_ETs(posSU, posET, B, n, m, nET, TVpower, delta, SUcellRadius, coverage, pathlossfactor, s)

posETs = [];    % the location of endusers of one WBS
posETsSUs = [];
SINR_ETs_oneWBS = [];
SINR_ETs = [];

posETcell = mat2cell(posET, [2], nET*ones(1, n));

for i = 1: n
    
    % generate pathlos matrix 
    posETsSUs = [posETcell{i}, posSU];
    d = dist(posETsSUs);    % get distance matrix for SUs
    G = ones(n+nET, n+nET)./ (d  + diag(ones(n+nET, 1))).^pathlossfactor;
    shadow = normrnd (0, s, n+nET, n+nET);
    G = G.*10.^(shadow/10);
    Gtilde = G - diag(diag(G));     % path loss between any ET of one WBS and any other WBS
    
    g_signal = Gtilde(1: nET, nET+i); % pathloss between all ETs and mother WBS, nET x 1
    
    currentChannel = find(B(i, :)); % which channel is being used:
%     listCochannelWBS = find(:, currentChannel); % list of co-channel WBSs
%     g_interference = Gtilde(1: nET, listCochannelWBS); % pathloss between all ETs and interfering WBS, nTE x m
    
    g_interference = Gtilde(1: nET, nET+1: nET+n); % nET  x n
    
    
    WBScochannel = B(:, currentChannel); % the power of WBS which interfers ETs, o means there is no interference caused.
    MumWBPower = WBScochannel(i); % power of mother WBS
    WBScochannel(i) = 0; % set the power of ith WB to 0
    
    SINR_ETs_oneWBS = 10*log10((g_signal*MumWBPower)./(g_interference* WBScochannel + delta)); % calculate SINR on all of them
    SINR_ETs = [SINR_ETs, SINR_ETs_oneWBS']; % store
   
end % every WBS
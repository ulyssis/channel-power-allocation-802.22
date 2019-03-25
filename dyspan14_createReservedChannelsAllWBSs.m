function [reservedChannelsAllWBSs] = dyspan14_createReservedChannelsAllWBSs(n, c)

reservedChannelsAllWBSs = cell(1, n);

for i = 1:n
    reservedChannelsAllWBSs(i) = {ones(1, c)};
    %contains the channel indice
    %reservedChannelsAllWBSs(i) = {1: c};
end
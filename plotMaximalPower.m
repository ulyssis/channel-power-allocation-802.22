function [] = plotMaximalPower(P, n, c)

%mp=cell(1,c);
mp = zeros(n, c);


for i = 1: c
	tem = P(:,i);
    tem(P(:, i)==0) = [];
	mp(:, i) = tem;
end

figure(8)



bar(mp);

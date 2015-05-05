% sum up the utilites over all nodes
% Obtain a n x 1 utility function 

    B=[0.000000 	1.000000 	0.000000; 
0.000000 	1.000000 	0.000000; 
0.000000 	0.000000 	1.000000; 
-0.000000 	1.000000 	-0.000000; 
0.000000 	0.000000 	1.000000; 
-0.000000 	-0.000000 	1.000000; 
1.000000 	0.000000 	0.000000; 
0.000000 	1.000000 	0.000000; 
1.000000 	-0.000000 	-0.000000; 
1.000000 	0.000000 	0.000000]; 

P = PHistory(:, :, 1);
Gtilde = GtildeHistory(:, :, 1);

sortedList = 1: c: n*c;
for i=1:n
	B0(i,:) = sum(P(sortedList(i): sortedList(i)+c-1, :), 1);
end
B = B0.* B;



[sumUtility_lindo, averageI_lindo, averageP_lindo, averageSINR, averageSINR2] = obtainPerformance(B, n, Gtilde, delta);
% sprintf('Lindo: all utility, averaged Iinterference, averaged P, averaged PsudoSINR, %f, %f, %f, %f', averageUtility_lindo*n, averageI_lindo,  averageP_lindo, averageP_lindo/(averageI_lindo+delta))
recordruns3 = [sumUtility_lindo, averageI_lindo, averageP_lindo, averageSINR, averageSINR2];
disp([recordruns0; recordruns1; recordruns2; recordruns3]);
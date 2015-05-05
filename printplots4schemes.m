% function printplots(powerHistory, averageSinrHistory, averageStdHistory, snrRatio_random, snrRatio_lindo, snrRatio_whitecat, snrRatio_whitecase, snrRatio_noregret)          
function printplots4schemes(powerHistory, averageSinrHistory, averageStdHistory, snrRatio_random, snrRatio_whitecat, snrRatio_whitecase, snrRatio_noregret)
% the latter snrRatio_random, snrRatio_whitecat, snrRatio_whitecase,
% snrRatio_noregret are sinr on each WBS in the last run, which is a
% intersection for sinr in one run


figure(2);
%            linetype =['b*-'; 'cx-'; 'g+-'; 'r<-'; 'k.-'];
           linetype =['k*-'; 'gx-'; 'r+-'; 'b<-'];

        
            subplot(2,2,1) % 
            for i=1:4
               hold on;
               semilogy (powerHistory(i, :), linetype(i,:));
            end
%             legend('random', 'lindo', 'whitecat','whitecaseish','noregret')
            legend('random', 'whitecat','whitecaseish','noregret')
%             set(gca,'YLim',[90 110],'Layer','top') 
            grid on
            title('average Power')
            xlabel('Topology index')
            ylabel([' Average power over all Secondary',sprintf('\n'),'Base Stations (w)']);


            subplot(2,2,2)
            for i=1:4
               hold on;
               semilogy (averageSinrHistory(i, :), linetype(i,:));
            end
%             set(gca,'YLim',[40 60],'Layer','top')  
            grid on
            title('average QuasiSINR')
            xlabel('Topology index')
            ylabel([' Average SINR over all Secondary',sprintf('\n'),'Base Stations']);


            subplot(2,2,3)
            for i=1:4
               hold on;
               semilogy (averageStdHistory(i, :), linetype(i,:));
            end
            grid on
            title('average Std of QuasiSINR')
            xlabel('Topology index')
            ylabel('Standard deviation')


            subplot(2,2,4)
%               sinr = bar ([snrRatio_random'; snrRatio_lindo'; snrRatio_whitecat'; snrRatio_whitecase'; snrRatio_noregret']');
            sinr = bar ([snrRatio_random'; snrRatio_whitecat'; snrRatio_whitecase'; snrRatio_noregret']');
            set(sinr(1),'Displayname','random');
            set(sinr(2),'Displayname','Lindo');
            set(sinr(2),'Displayname','whitecat');
            set(sinr(3),'Displayname','whitecaseish');
            set(sinr(4),'Displayname','noregret');
            legend('Location','northeast');
            grid on   
            title('SINR of SUs in one sinario')
            xlabel('Secondary base station index')
            ylabel('SINR of each SBS')

            
            
        
figure(3);     
%   average Power
%   average QuasiSINR

%   average qusai SINR on reference position
            subplot(1,2,1)
%             labels = {'rand'  'lindo' 'whitecat' 'whitecase'  'noreg'};
            labels = {'rand' 'whitecat' 'whitecase'  'noreg'};
            bar((1:4), mean(powerHistory,2)');
            set(gca, 'XTick', 1:4, 'XTickLabel', labels);
%             set(gca,'YLim',[80 110],'Layer','top') 
            h = gca;
            H = rotateticklabel(h, 90); 
            grid on
            title(['average power consumption of ',sprintf('\n'),'all WBs over different simulations'])
            ylabel(['Average value of the averaged power of all SBs',sprintf('\n'),'over all scenarios (w)']);
            
            
            
            subplot(1,2,2)
%             labels = {'rand'  'lindo' 'whitecat' 'whitecase'  'noreg'};
            labels = {'rand' 'whitecat' 'whitecase'  'noreg'};
            bar((1:4), mean(averageSinrHistory,2)'); 
            set(gca, 'XTick', 1:4, 'XTickLabel', labels);            
%             set(gca,'YLim',[40 60],'Layer','top') 
            h = gca;
            H = rotateticklabel(h, 90); 
            grid on
            title(['average qusai SINR on reference position ',sprintf('\n'),' of all WBs over different simulations'])
            ylabel(['Average value of the averaged qusai SINR (db) of all',sprintf('\n'),' SBs over all scenarios']);
   




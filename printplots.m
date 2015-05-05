% function printplots(powerHistory, averageSinrHistory, averageStdHistory, snrRatio_random, snrRatio_lindo, snrRatio_WhiteCat, snrRatio_self, snrRatio_Noregret)          
function printplots(powerHistory, averageSinrHistory, averageStdHistory, snrRatio_random, snrRatio_WhiteCat, snrRatio_self, snrRatio_lindo, snrRatio_Noregret)          


figure(2);
            linetype =['b*-'; 'cx-'; 'g+-'; 'r<-'; 'k.-'];
%            linetype =['k*-'; 'gx-'; 'r+-'; 'b<-'];
        
            subplot(2,2,1)
            for i=1:5
               hold on;
               semilogy (powerHistory(i, :), linetype(i,:));
            end
%             set(gca,'YLim',[90 110],'Layer','top') 
            grid on
            title('average Power')
            xlabel('Topology index')
            ylabel([' Average power over all Secondary',sprintf('\n'),'Base Stations (w)']);


            subplot(2,2,2)
            for i=1:5
               hold on;
               semilogy (averageSinrHistory(i, :), linetype(i,:));
            end
%             set(gca,'YLim',[40 60],'Layer','top')  
            grid on
            title('average QuasiSINR')
            xlabel('Topology index')
%             ylabel([' Average SINR over all Secondary',sprintf('\n'),'Base Stations (dB)']);
            ylabel([' Average SINR',sprintf('\n'),'on Base Stations']);

            
            %subplot(2,2,3);
            legend('random', 'lindo', 'WhiteCat','WhiteCase','Noregret');
            axes('Position',[0.13 0.11 0.17 0.341162790697674]);
            for i=1:5
               hold on;
               semilogy (averageStdHistory(i, :), linetype(i,:));
            end
            axis tight;
            title('average Std of QuasiSINR')
            xlabel('Topology index')
            ylabel('Standard deviation')


            
            
            subplot(2,2,4)
%              sinr = bar ([log(snrRatio_random)'; log(snrRatio_WhiteCat)'; log(snrRatio_self)'; log(snrRatio_Noregret)']');
%            sinr = bar ([snrRatio_random'; snrRatio_lindo'; snrRatio_WhiteCat'; snrRatio_self'; snrRatio_Noregret']');
            sinr = bar ([snrRatio_random'; snrRatio_WhiteCat'; snrRatio_self'; snrRatio_lindo'; snrRatio_Noregret']');
%             set(sinr(1),'Displayname','random');
%             set(sinr(2),'Displayname','Lindo');
%             set(sinr(3),'Displayname','WhiteCat');
%             set(sinr(4),'Displayname','WhiteCase');
%             set(sinr(5),'Displayname','Noregret');
%             legend('Location','northeast');
            axis tight;
            title('SINR of SUs in one sinario')
            xlabel('Secondary base station index')
            ylabel('SINR of each SBS(dB)')
            
            % set the legend where I really want it.            
            set(legend,'Position',[0.34 0.11 0.147741935483871 0.341162790697674]);
        
            
figure(3);        
            subplot(2,2,1)
            labels = {'rand'  'lindo' 'WhiteCat' 'WhiteCase'  'noreg'};
%             labels = {'rand' 'WhiteCat' 'self'  'noreg'};

            bar((1:5), mean(powerHistory,2)');
            set(gca, 'XTick', 1:5, 'XTickLabel', labels);
%             set(gca,'YLim',[80 110],'Layer','top') 
            h = gca;
            H = rotateticklabel(h, 90); 
            grid on
            title('average power consumption')
            ylabel([' Averaged averaged power of all SBs',sprintf('\n'),'over all scenarios (w)']);
            
            
            
            subplot(2,2,2)
            labels = {'rand'  'lindo' 'WhiteCat' 'WhiteCase'  'noreg'};
%             labels = {'rand' 'WhiteCat' 'self'  'noreg'};
            %H=get(gcf,'CurrentAxes');
            bar((1:5), mean(averageSinrHistory,2)'); 
            set(gca, 'XTick', 1:5, 'XTickLabel', labels);            
%             set(gca,'YLim',[40 60],'Layer','top') 
            h = gca;
            H = rotateticklabel(h, 90); 
            grid on
            title('average SINR')
            ylabel([' Averaged averaged SINR of all',sprintf('\n'),' SBs over all scenarios(dB)']);
            
                        
            subplot(2,1,2)
            labels = {'random' 'lindo' 'WhiteCat' 'WhiteCase' 'Noregret'};
%             labels = {'rand' 'WhiteCat' 'self'  'noreg'};
            bar((1:5), mean(averageStdHistory,2)');
            set(gca, 'XTick', 1:5, 'XTickLabel', labels);
            grid on
            title('average standard divation')
            ylabel(['Average standard deviation of SINR',sprintf('\n'),'over all scenarios'])    

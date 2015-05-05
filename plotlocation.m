function []= plotlocation(n, m, lengthSide, posSU, posET, posPUContor)
        % plot the topology of uesrs 
        hold off;
        figure(1);
        plot(posSU(1,:), posSU(2,:), 'r+');
        hold on;
        plot(posET(1,:), posET(2,:), 'bo');
        indexSU = (1:n)';
        hold on;
        plot(posPUContor(1,:), posPUContor(2,:), 'b<');
        indexPU = (1:m)';
        text(posSU(1,:), posSU(2,:), num2str(indexSU));
        text(posPUContor(1,:), posPUContor(2,:), num2str(indexPU));
        leg3=legend('Secondary users');
        axis([-lengthSide lengthSide*2 -lengthSide lengthSide*2]);
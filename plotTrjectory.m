function [  ] = plotTrjectory( bd,fd,rawData ,note)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i = 2:size(rawData,1)
    tnt=find(fd{i-1,1}(end-5:end,11),1,'first');
    h=figure();
    subplot(2,1,1);
    hist=rawData{i,2};
    hist=hist(hist<=fd{i-1,1}(end-6,1));
    plot(hist/12+min(rawData{i,4})/365,rawData{i,3}(1,1:length(hist)),'*b')
    hold on
   % plot(fd{i-1,1}(:,1)./12+min(rawData{i,4})/365,fd{i-1,1}(:,2),'or')
    plot(fd{i-1,1}(1:end-6,1)./12+min(rawData{i,4})/365,fd{i-1,1}(1:end-6,5),'s-k')
    
    errorbar(fd{i-1,1}(end-5:end,1)./12+min(rawData{i,4})/365,fd{i-1,1}(end-5:end,5),fd{i-1,1}(end-5:end,8).*1.96,'^-k')
    plot(fd{i-1,1}(end-6:end-5,1)./12+min(rawData{i,4})/365,fd{i-1,1}(end-6:end-5,5),'-k')
    legend('Observed','Filtered','Predicted w/ 95% Conf.','Location','best')
    xlabel('Age (year)')
    ylabel('MD')
   % xlim([min(rawData{i,4})/365,ceil(max(fd{i-1,1}(:,1))./12+min(rawData{i,4})/365)+1])
   % ylim([-20 -5])
    grid on
    grid minor
   % text(fd{i-1,1}(end-5:end,1)./12+min(rawData{i,4})/365-0.25,fd{i-1,1}(end-5:end,5)+4,(num2cell(round(fd{i-1,1}(end-5:end,12).*100,0)./100)))
    title(['Time to Next Test = ' num2str(tnt*6) ' mo'] )
    set(gca,'FontSize',14)
    
    subplot(2,1,2);
    hist=rawData{i,2};
    hist=hist(hist<=fd{i-1,1}(end-6,1));
    plot(hist/12+min(rawData{i,4})/365,rawData{i,3}(2,1:length(hist)),'*b')
    hold on
  %  plot(fd{i-1,1}(:,1)./12+min(rawData{i,4})/365,fd{i-1,1}(:,3),'or')
    plot(fd{i-1,1}(1:end-6,1)./12+min(rawData{i,4})/365,fd{i-1,1}(1:end-6,6),'s-k')
    
    errorbar(fd{i-1,1}(end-5:end,1)./12+min(rawData{i,4})/365,fd{i-1,1}(end-5:end,6),fd{i-1,1}(end-5:end,9).*1.96,'^-k')
    plot(fd{i-1,1}(end-6:end-5,1)./12+min(rawData{i,4})/365,fd{i-1,1}(end-6:end-5,6),'-k')
    legend('Observed','Filtered','Predicted w/ 95% Conf.','Location','best')
    xlabel('Age (year)')
    ylabel('IOP')
  %  xlim([min(rawData{i,4})/365,ceil(max(fd{i-1,1}(:,1))./12+min(rawData{i,4})/365)+1])
   % ylim([0 25])
    grid on
    grid minor
    %text(fd{i-1,1}(end-5:end,1)./12+min(rawData{i,4})/365-0.25,fd{i-1,1}(end-5:end,5)+8,(num2cell(round(fd{i-1,1}(end-5:end,12).*100,0)./100)))
   % title(['IOP Plot for Patient ',rawData{i,1}, note] )
   set(gca,'FontSize',14)
    set(h,'PaperUnits','Inches');
    set(h, 'PaperSize', [8.5 11*0.6]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition', [0 0 8.5 11*0.6]);
    print(h,[rawData{i,1} note '.pdf'],'-dpdf','-r0');
    close all;
end
end


%% Program to produce figures from particles

x=PARTICLE(:,1:20);
y=LL;
z=PARTICLE(:,21:end);

%plot maximum likelihood value
figure
[maximum,maxlike]=max(y);
disp(maximum)
% hold on
x0=PlotPARTICLES( maxlike, x, z,1);
drawnow

%%  Plot matrix of correlations and posteriors
EstIndex=[2:5,9:17,19:20];
figure
[S,AX,BigAx,H,HAx] =plotmatrix(PARTICLE(:,EstIndex),'k.');  h=get(gcf,'children'); set(h,'FontSize',10)
Correlations=corrcoef(PARTICLE(:,EstIndex));
colors=hot(11);
Correlations=abs(round(Correlations,1));
colormapvec=linspace(0,1,11);

for i =1:length(EstIndex)
    for j =1:length(EstIndex)
    colorind=find((colormapvec<Correlations(i,j)+0.05 & colormapvec>Correlations(i,j)-0.05));
    S(i,j).Color=colors(12-colorind,:);
    end
end   

%% Plot sample particles
k=20;

figure

sampleIndices=randsample(length(y),k);    %collect a vector of indices for 100 samples of chain

for i=1:k
    hold on
    modelOutput(:,:,i)=PlotPARTICLES( sampleIndices(i), x, z,0);
    drawnow
    h=get(gcf,'children');
    
    for j=1:7
        AgeGroupCounts(i,:,j)=sum(modelOutput(33+(j-1)*52+1:33+j*52,:,i));
    end
end
AgeGroupCounts=reshape(AgeGroupCounts,[],7);

load('GermanCaseNotification.mat')
StratifiedCases= AgeStratify( GermanCaseNotification, ageGroupBreaks );
medmod=median(modelOutput,3);
for index=1:7
subplot(4,2,index) ; hold on
plot(medmod(:,index).'/PopulationSize(index),'k.','LineWidth',1);
end

set(h,'xlim',[0,418],'xgrid','on','FontSize',16, ...
    'xticklabel',{'2008','2009','2010','2011','2012','2013','2014','2015','2016'},'XTick',[33,85,137,189,241,293,345,397])

figure
boxplot((AgeGroupCounts),'Boxstyle','filled','colorgroup',[1:7],'colors',winter(7),'labels',{'0-8','8-18','18-26','26-37','37-50','50-70','70+'}); 
hold on;
for i=1:7; plot(sum(StratifiedCases(33+(i-1)*52+1:33+i*52,:)),'ko','markersize',10,'LineWidth',3);end
ylabel('Number of cases per complete season')
ylim([0,inf])
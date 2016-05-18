lambda=[0.5 1 2 3 4 5 6 7];
t=0:0.01:5;
stable=figure;
hold on;
for i=1:length(lambda)
  plot(t,e.^(-lambda(i)*t)+0.2/lambda(i),'-r','LineWidth',2);
end
xlabel("t");
ylabel("y");
grid on;
set(gca,'FontSize',30,'fontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold');
saveas(stable,['stable','.eps'],'eps') 

unstable=figure;
hold on;
for i=1:length(lambda)
  plot(t,e.^(lambda(i)/max(lambda)*t)+0.2/lambda(i),'-r','LineWidth',2);
end
xlabel("t");
ylabel("y");
grid on;
set(gca,'FontSize',30,'fontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold');
saveas(unstable,['unstable','.eps'],'eps')

asyn_stable=figure;
hold on;
for i=1:length(lambda)
  plot(t,t+lambda(i),'-r','LineWidth',2);
end
xlabel("t");
ylabel("y");
grid on;
set(gca,'FontSize',30,'fontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold');
saveas(asyn_stable,['asyn_stable','.eps'],'eps')
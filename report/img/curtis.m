h_im=0.25;
h_ei=0.25;
h_fin=1/1000;
lambda=-50;
t_im=0:h_im:2;
t_ei=0:h_ei:2;
t_fin=0:h_fin:2;
y_ei=zeros(1,length(t_ei));
y_im=zeros(1,length(t_im));
y_ei(1)=0.15;
y_im(1)=0.15;
for i=2:length(t_im)
  y_im(i)=1/(1-h_im*0.5*lambda)*((1+0.5*h_im*lambda)*y_im(i-1)-h_im*lambda*cos(t_im(i-1)+h_im/2));
end
for i=2:length(t_ei)
  y_ei(i)=1/(1-h_ei*lambda)*(y_ei(i-1)-h_ei*lambda*cos(t_ei(i)));
end  
imp_mid=figure
plot(t_im,y_im,'-k','LineWidth',2)
hold on
plot(t_fin,lambda*lambda/(lambda*lambda+1)*cos(t_fin)-lambda/(lambda*lambda+1)*sin(t_fin)-lambda*lambda/(lambda*lambda+1)*e.^(lambda*t_fin),'--k','LineWidth',2)
xlabel("t");
ylabel("y");
axis([0 2 -0.5 2]);
grid on;
set(gca,'FontSize',30,'fontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold');
saveas(imp_mid,['imp_mid','.eps'],'eps') 
eul_imp=figure
plot(t_ei,y_ei,'-k','LineWidth',2);
hold on
plot(t_fin,lambda*lambda/(lambda*lambda+1)*cos(t_fin)-lambda/(lambda*lambda+1)*sin(t_fin)-lambda*lambda/(lambda*lambda+1)*e.^(lambda*t_fin),'--k','LineWidth',2)
xlabel("t");
ylabel("y");
axis([0 2 -0.5 2]);
grid on;
set(gca,'FontSize',30,'fontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold');
saveas(eul_imp,['eul_imp','.eps'],'eps') 
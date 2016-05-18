%% Lotka-Volterra model
h=0.1;
t=0:h:6.4;
u_split_1=zeros(1,length(t));
v_split_1=zeros(1,length(t));
u_split_2=zeros(1,length(t));
v_split_2=zeros(1,length(t));
u_split_2_supp=zeros(1,length(t));
u_split_1(1)=6;
v_split_1(1)=2;
u_split_2(1)=6;
v_split_2(1)=2;
for i=2:length(t)
    u_split_1(i)=exp(h*(v_split_1(i-1)-2))*u_split_1(i-1);
    v_split_1(i)=exp(h*(1-u_split_1(i)))*v_split_1(i-1);
    u_split_2_supp(i)=exp(h/2*(v_split_2(i-1)-2))*u_split_2(i-1);
    v_split_2(i)=exp(h*(1-u_split_2_supp(i)))*v_split_2(i-1);
    u_split_2(i)=exp(h/2*(v_split_2(i)-2))*u_split_2_supp(i);    
end
plot(u_split_1,v_split_1,'ob');
hold on;
plot(u_split_2,v_split_2,'xr');

options = odeset('RelTol',1e-7,'AbsTol',[1e-7 1e-7]);
[T,Y] = ode45(@rigid,[0 6.4],[6 2],options);

plot(Y(:,1),Y(:,2),'-k');
legend('Lie-Trotter splitting','Strang splitting','Exact solution')
xlabel('Predators')
ylabel('Preys')
grid minor
set(gca,'FontSize',30,'fontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold');

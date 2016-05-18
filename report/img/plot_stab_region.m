[x,y] = meshgrid(linspace(-4,0),linspace(-2,2));
hlambda = x + i*y;
gamma=1/3;
StabilityFunction = (1+hlambda*(1-2*gamma)+hlambda*hlambda*(0.5-2*gamma+gamma*gamma))./(1-gamma*hlambda).^2;
stab_domain=contourf(x,y,-abs(StabilityFunction),-[1 1])
axis equal, axis([-4 2 -2 2]), grid on
set(gca,'FontSize',30,'fontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold');

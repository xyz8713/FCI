function comph_c(coefs,shift,r,c)
x = linspace(c+r-1,c+r+7,200);
y = fh_c(coefs, shift, x);
plot(x,y,'b-','linewidth',2);
hold on
x1 = linspace(r+c,r+c+7, 160);
plot(x1, 1 ./ x1,'r-.','linewidth',2)
%ylim([-0.2,1/(r+c)+1.2])
ylim([0,1])
hold on;
grid on;
plot([c,c],[0,0],'b*')
plot([c+r,c+r],[0,1/(c+r)+1],'g--')
lngd = legend('8 poles', 'Exact 1/x');
set(lngd, 'Location', 'northeast');
outerposition = get(lngd, 'OuterPosition');
set(lngd, 'interpreter','latex', 'fontsize', 20);
tightInset = get(gca, 'TightInset');
position(1) = tightInset(1);
position(2) = tightInset(2);
position(3) = 1 - tightInset(1) - tightInset(3);
position(4) = 1 - tightInset(2) - tightInset(4);
set(gca, 'Position', position);
pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)
end





function plotnodes_c(shift,r,c)
x1 = linspace(-pi,pi,100);
y1 = r*exp(1i*x1)+c;
plot(real(y1),imag(y1))
hold on;
axis equal
plot(real(shift),-imag(shift),'ro','MarkerSize', 12);
hold on;
plot([c,c],[0,0],'b*')
plot([0,0],[-r,r],'b--')
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
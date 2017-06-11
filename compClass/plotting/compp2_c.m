function compp2_c(coefs,shift,r,c)
x = linspace(c-r-1,c+r+1,200);
y = fp2_c(coefs, shift, x);
plot(x,y,'b-','linewidth',2);
hold on;
grid on;
plot([c,c],[0,0],'b*')
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

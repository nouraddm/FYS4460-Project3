function UIcall(slidehandle)
% a=get(slidehandle,'Value')
% x=linspace(0,4*pi,100);
% y=sin(a*x);
% plot(x,y,'-r') 
p=get(slidehandle,'Value')
L=100;
r = rand(L,L);
z = r < p;
[lw,num] = bwlabel(z,4);

img = label2rgb(lw,'jet','k','shuffle');
image(img);


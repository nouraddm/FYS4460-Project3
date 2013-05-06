%-------------------------------------------------------------------------
%                          
%                           Project 3
%
%                          Exercise (a)
%
%-------------------------------------------------------------------------

close all;

L = 100;
p = 0.1;

fighandle=figure(1)
slidehandle = uicontrol(fighandle,'Style','Slider')
set (slidehandle,'Value',1,'Min',0,'Max',1,'Callback','UIcall(slidehandle)')
set (slidehandle,'SliderStep',[0.1 1])

r = rand(L,L);
z = r < p;
[lw,num] = bwlabel(z,4);

img = label2rgb(lw,'jet','k','shuffle');
image(img);


% 
% s1 = regionprops(lw,'Area');
% s2 = regionprops(lw,'BoundingBox');
% 
% area = cat(1,s1.Area)
% bbox = cat(1,s2.BoundingBox);
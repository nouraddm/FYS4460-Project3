%-------------------------------------------------------------------------
%                          
%                           Project 3
%                           
%                         Exercise (e)
%
%                    Cluster number density 
%-------------------------------------------------------------------------
close all;
abovePc = 0;

pc = 0.59275;
numberOfSamples = 1;

if(abovePc == 0)
    p = [0.5 0.54 0.57 pc];
else
    p = [pc 0.6 0.64 0.68];
end

sizeOfP = size(p,2);
tau = 187/91;

Lx = 256;
Ly = Lx;
A = Lx*Ly;

logbinsize = 2; 
logbinmax = A;
 
for i = 1 : sizeOfP
    for j = 1 : numberOfSamples
        z = rand(Lx,Ly);
        zz = z < p(i);
        [lw,num] = bwlabel(zz,4);
        perc_y = intersect(lw(:,1),lw(:,Ly)); 
        perc_x = intersect(lw(1,:),lw(Lx,:)); 
        perc_xy = union(perc_x,perc_y);
        perc = find(perc_xy >0);
        s = regionprops(lw,'Area'); 
        clusterareas = cat(1,s.Area);

        % Find the cluster number density, get rid of percolating clusters
        ind = (1: num );
        indnoP = setxor(ind,perc_xy(perc));
        % Do statistics on area(indnoP)
        clusta = clusterareas(indnoP);
        [x,dx,n] = logbin(clusta,logbinsize,logbinmax); 
        if (j == 1)
            nnsp = n/A;
            nnsp = nnsp'./dx;
            nsp = nnsp;
        else
            nnsp = n/A;
            nnsp = nnsp'./dx;
            nsp = nsp + nnsp; 
        end
    end
    nsp = nsp/numberOfSamples;
    ind2 = find(nsp > 0);
 
    % n(s,p) plot
    
    linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',...
    ':','-',':','-.','--','-',':','-.','--','-',':','-.'));
    MarkerEdgeColors=jet(21);  
    Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.',...
        '+','*','o','x','^','<','h','.','>','p','s','d','v',...
        'o','x','+','*','s','d','v','^','<','>','p','h','.'];

    figure(1);
    plot(log10(x(ind2)),log10(nsp(ind2)),[linestyles{i} Markers(i)],'Color',MarkerEdgeColors(i,:));
    xlabel('log10(s)');
    ylabel('log10(n(s,p_c))');

    drawnow;
    hold on;
end
hold off;
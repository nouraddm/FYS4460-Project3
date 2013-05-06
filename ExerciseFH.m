%-------------------------------------------------------------------------
%                          
%                           Project 3
%                           
%                       Exercise (f) & (h)
%
%                     Cluster number density 
%                               &
%                 Mass scaling of percolating cluster
%-------------------------------------------------------------------------
close all;

pc = 0.59275;
numberOfSamples = 10;
tau = 187/91;

lengthL = 6;

P = zeros(lengthL,1);
M = zeros(lengthL,1);

for i = 1:lengthL
    Lx = 2^(i+3);
    Ly = Lx;
    A = Lx*Ly;
    L(i) = Lx;

    logbinsize = 2; 
    logbinmax = A;
 
    for j = 1 : numberOfSamples
        z = rand(Lx,Ly);
        zz = z < pc;
        [lw,num] = bwlabel(zz,4);
        perc_y = intersect(lw(:,1),lw(:,Ly)); 
        perc_x = intersect(lw(1,:),lw(Lx,:)); 
        perc_xy = union(perc_x,perc_y);
        perc = find(perc_xy >0);
        s = regionprops(lw,'Area'); 
        clusterareas = cat(1,s.Area);
        if (length(perc)>0)
            ar = sum(clusterareas(perc_xy(perc)));
            P(i) = P(i) + ar/A; 
        end

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
    P(i) = P(i)/numberOfSamples;
    M(i) = Lx^(2)*P(i);
    nsp = nsp/numberOfSamples;
    ind2 = find(nsp > 0);
 
    % n(s, pc; L) plot
    
    linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',...
    ':','-',':','-.','--','-',':','-.','--','-',':','-.'));
    MarkerEdgeColors=jet(21);  
    Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.',...
        '+','*','o','x','^','<','h','.','>','p','s','d','v',...
        'o','x','+','*','s','d','v','^','<','>','p','h','.'];

    figure(1);
    plot(log10(x(ind2)),log10(nsp(ind2)),[linestyles{i} Markers(i)],'Color',MarkerEdgeColors(i,:));
    xlabel('log10(s)');
    ylabel('log10(n(s,p_c;L))');

    drawnow;
    hold on;
end

figure(2);
plot(log10(L),log10(M))
xlabel('log10(L)');
ylabel('log10(M)');
hold off;
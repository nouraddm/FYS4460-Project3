%-------------------------------------------------------------------------
%                          
%                           Project 3
%
%                          Exercise (c)
%
%-------------------------------------------------------------------------


Pcumsum = zeros(1e5,1);
Pcum = zeros(1e5,1);
xsum = zeros(1e5,1);

for k = 1:10
    zc        = rand(1e5,1).^(-2);
    unqZC     = unique(zc);
    xsum = xsum + unqZC; 
    countZC = histc(zc,unqZC);
    relFreq   = countZC/numel(zc);
    
    for i = 1:length(unqZC)
        Pcum(i) = 1.0 - sum(relFreq(1:i));
        Pcumsum(i) = Pcumsum(i) + Pcum(i);
    end
end

Pcumsum = Pcumsum/10;
xsum = xsum/10;

%i = 1;
%while(Pcumsum(i) > 0)
 %   y(i) = Pcumsum(i);
%    i = i+1;
%end

figure(1);
plot(log10(xsum),log10(Pcumsum))


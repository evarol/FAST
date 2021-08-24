function h=radial_histogram(im,location,maxradius)
%2D image
prevblock=zeros(size(im));
for i=1:maxradius
    currentblock=zeros(size(im));
    currentblock(location(1)-(i-1):location(1)+(i-1),location(2)-(i-1):location(2)+(i-1))=1;
    currentmask=currentblock-prevblock;
    prevblock=currentblock;
    h(i)=mean(im(currentmask==1));
end
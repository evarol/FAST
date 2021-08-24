function newpatchradius=patch_bound(location,patchradius,imagesize)
%given a prescribed patch size, truncate the patch (keeping it centered) if
%it goes over image bounds
newpatchradius=zeros(size(patchradius));
for i=1:length(location)
newpatchradius(i)=min(min(location(i)-1,imagesize(i)-location(i)),patchradius(i));
end

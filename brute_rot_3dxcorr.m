function [bestxcorr,Rbest]=brute_rot_3dxcorr(X,Y,sigma,anglelimit,gridsize)

angle_1=anglelimit(1,1):gridsize:anglelimit(1,2);
angle_2=anglelimit(2,1):gridsize:anglelimit(2,2);
angle_3=anglelimit(3,1):gridsize:anglelimit(3,2);
bestxcorr=0;
for i=1:length(angle_1)
    for j=1:length(angle_2)
        for k=1:length(angle_3)
            [i j k]
            R=eul2rotm(deg2rad([angle_1(i),angle_2(j),angle_3(k)]));
            Yrot=Y*R;
            xc=sum(sum(exp(-pdist2(X,Yrot).^2/2*sigma^2)));
            if bestxcorr<xc
                bestxcorr=xc;
                Rbest=R;
            end
        end
    end
    
end
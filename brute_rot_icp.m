function [bestmatch,Rbest]=brute_rot_icp(X,Y,indelcost,anglelimit,gridsize)

angle_1=anglelimit(1,1):gridsize:anglelimit(1,2);
angle_2=anglelimit(2,1):gridsize:anglelimit(2,2);
angle_3=anglelimit(3,1):gridsize:anglelimit(3,2);
bestmatch=0;
for i=1:length(angle_1)
    for j=1:length(angle_2)
        for k=1:length(angle_3)
            R=eul2rotm([angle_1(i),angle_2(j),angle_3(k)]);
            Yrot=Y*R;
            [~,D] = knnsearch(Yrot,X,'K',1);
            if bestmatch<sum(D<indelcost)
                bestmatch=sum(D<indelcost);
                Rbest=R;
            end
        end
    end
    
end
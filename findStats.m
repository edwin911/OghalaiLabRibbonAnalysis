function [Radii,Centers,DCircles,DRadii] = findStats(app,type,slice)
    if type==1
        Centers=[];
        Radii=[];
        counter=1;
        for i = 1:size(app.NCenters,1)
            if app.NRadii(i)^2>(app.NCenters(i,3)-slice*app.voxel(3))^2

                Radii(counter)=((((app.NRadii(i))^2)-((app.NCenters(i,3)-slice*app.voxel(3))^2))^(.5))/app.voxel(1);
                Centers(counter,:)=[app.NCenters(i,1)/app.voxel(1),app.NCenters(i,2)/app.voxel(2),app.NCenters(i,3)/app.voxel(1)];
                counter=counter+1;
            end
        end
        A=app.discardedN(:,3)==slice;
        DCircles=app.discardedN(A,1:2);
        DRadii=app.discardedN(A,4);

    elseif type==2
        A=app.discardedRPre(:,3)==slice;
        Centers=app.presynapticPositions(slice).points;
        Radii=str2double(app.RibbonRadiusEditField.Value);
        DCircles=app.discardedRPre(A,:);
        DRadii=app.RibbonRadiusEditField;
    elseif type==3
        A=app.discardedRPost(:,3)==slice;
        Centers=app.postsynapticPositions(slice).points;
        Radii=str2double(app.RibbonRadiusEditField.Value);
        DCircles=app.discardedRPost(A,:);
        DRadii=app.RibbonRadiusEditField;
    end
end
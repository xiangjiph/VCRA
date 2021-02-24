function myplot3(pts,args)
if iscell(args)
    plot3(pts(:,1),pts(:,2),pts(:,3),args{:})
else
    plot3(pts(:,1),pts(:,2),pts(:,3),args)
end
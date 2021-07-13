% colors
mycol=[1 .4 0; 1 0 1; 0 0.7 0.7];

% loop over shapes
for k=2:3
    
    % select shape
    if k==1
        mytulip
    elseif k==2
        myellipse3
    elseif k==3
        myellipse4
    end
    
    % clear folder
	system('rm TEMP/*');

    % set up features and speedmatrix
    [X,Y]=meshgrid(-1.5:0.01:3,-2.5:0.01:2.5);
    myinpol=inpolygon(X(:),Y(:),yvals,xvals);
    speedmatrix=1.0-reshape(myinpol,size(X,1),size(X,2));

    h=patch(yvals,xvals,mycol(k,:)); hold on;
    h.FaceAlpha=0.1;
    h.EdgeColor='none';
    axis equal;
    plot(yvals,xvals,'-k');
    
    % FMM
    dlmwrite('TEMP/myX.txt',X,' ');
    dlmwrite('TEMP/myY.txt',Y,' ');
    dlmwrite('TEMP/myspeedmatrix.txt',speedmatrix,' ');
    system('source ~/.bash_profile; source activate skfmm; python singlefeaturefmm.py');
    
    % identify individual fronts
    tempfronts=dlmread('TEMP/fronts.dat');
    fronts=tempfronts;
    counter=1; % counter for fronts
    fronts(1,1)=counter;
    for i=2:size(tempfronts,1)
        % identify whether this is a new front now
        % break up front if long stretch in direction of front propagation
        % (such as through obstacle)
        if ~isequal(fronts(i,2),fronts(i-1,2)) || abs(fronts(i,4)-fronts(i-1,4))>0.1
            counter=counter+1;
        end
        % write the front id into fronts
        fronts(i,1)=counter;
    end
    
    % loop over front ids
    for i=(unique(fronts(:,1)))'
    %for i=7:7
        mylocfront=fronts(fronts(:,1)==i,3:4);
        % plot front as line
        pl=plot(mylocfront(:,1),mylocfront(:,2),'.','Color',mycol(k,:),'LineWidth',2);
        pl.Color(4)=0.4;
    end

end

% plot analytical solution for long distances
halfcirc=[];
for phi=-pi/2:0.01:pi/2
    halfcirc=cat(1,halfcirc,[1.75*cos(phi),1.75*sin(phi)]);
end
plot(halfcirc(:,1),halfcirc(:,2)+1,'--','Color','black','LineWidth',1.5);
plot(halfcirc(:,1),halfcirc(:,2)-1,'--','Color','black','LineWidth',1.5);
plot(0,1,'o','Color','black','MarkerSize',7,'LineWidth',2);
plot(0,-1,'o','Color','black','MarkerSize',7,'LineWidth',2);

axis([-1.5 3 -2.5 2.5])
axis off
hold off
saveas(gcf,'OUTPUT/obstacles_ellipses.png')

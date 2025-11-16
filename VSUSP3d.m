% Script to find suspension parameters of a 3D quarter car model through
% bump and steer

% Input values are carParams structure and maxDroop, maxJounce, maxLeft,
% and maxRight values. carParams structure must include members inboardF,
% outboardF, outerTire, tireOutboard, tireInboard, and tireContactPoint.
% These members are arrays of x, y, and z coordinates. 

% Coordinates are measured in mm. Positive x is forward, positive y is
% left, and positive z is up. outboardF and inboardF members must be a
% matrix in the format:
% [upper fore x, upper fore y, upper fore z;
%  upper aft x,  upper aft y,  upper aft z;
%  lower fore x, lower fore y, lower fore z;
%  lower aft x,  lower aft y,  lower aft z;
%  pushrod x,    pushrod y,    pushrod z;
%  tie rod x,    tie rod y,    tie rod z]

function VSUSP3d(carParams, maxDroop, maxJounce, maxLeft, maxRight)
    increment = 0.1;                                                      % Increment of travel from 0 to maximum values
     
    if maxDroop == 0 && maxJounce == 0
        midBump = 1;                                                       % Middle index of bump
    else
        % midBump = round(((maxJounce-maxDroop)/increment)/2+1)
        midBump = round(abs(maxDroop)/increment)+1;
    end

    if maxLeft == 0 && maxRight == 0
        midSteer = 1;
    else
        midSteer = round(abs(maxLeft)/increment+1);                            % Middle index of steer
    end

    inboard = carParams.inboardF([1:4,6],:);                               % Defines beginning inboard coordinates
    outboard = carParams.outboardF([1:4,6],:);                             % Defines beginning outboard coordinates
    tireOutboard = carParams.tireOutboardF;                                % Defines beginning outer point of tire normal vector
    tireInboard = carParams.tireInboardF;                                  % Defines beginning outer point of tire normal vector
    tireContactPt = carParams.tireContactPtF(1,:);                         % Defines beginning tire contact point
    lengths = zeros(1,13);
    lengths(1) = norm(outboard(3,[2,3])-inboard(3,[2,3]));                 % Defines known lengths between different points
    lengths(2) = norm(outboard(1,:)-inboard(1,:));
    lengths(3) = norm(outboard(1,:)-inboard(2,:));
    lengths(4) = norm(outboard(5,:)-inboard(5,:));
    lengths(5) = norm(outboard(1,:)-outboard(3,:));
    lengths(6) = norm(outboard(1,:)-outboard(5,:));
    lengths(7) = norm(outboard(3,:)-outboard(5,:));
    lengths(8) = norm(tireInboard-outboard(1,:));
    lengths(9) = norm(tireInboard-outboard(3,:));
    lengths(10) = norm(tireInboard-outboard(5,:));
    lengths(11) = norm(tireOutboard-outboard(1,:));
    lengths(12) = norm(tireOutboard-outboard(3,:));
    lengths(13) = norm(tireOutboard-outboard(5,:));
    lengths(14) = norm(tireContactPt-tireInboard);
    lengths(15) = norm(tireContactPt-tireOutboard);
    lengths(16) = norm(tireContactPt-outboard(5,:));

    params = zeros(6,round(1+(maxJounce-maxDroop)/increment),round(1+(maxRight-maxLeft)/increment));    % Creates array to store parameters (toe, camber, caster etc.)
    coords = zeros(7,3,round(1+(maxJounce-maxDroop)/increment),round(1+(maxRight-maxLeft)/increment));  % Creates array to store suspension coordinates
    tire = zeros(101,12,round(1+(maxJounce-maxDroop))/increment,round(1+(maxRight-maxLeft)/increment)); % Creates array to store tire coordinates
    
    coordsNew = [outboard(1,:); outboard(3,:); outboard(5,:); tireInboard; tireOutboard; tireContactPt; inboard(5,:)];  % Places starting coordinates in easier to call variable

    [coords(:,:,midBump,midSteer), params(2:6,midBump,midSteer), tire(:,:,midBump,midSteer)] = solveBump(coordsNew, inboard, lengths, 0);   % Finds parameters at static

    for i = midSteer-1:-1:1                                                                                     % Works from zero to max right steer
        [coordsNew, params(2:6,midBump,i), tire(:,:,midBump,i)] = solveSteer(coordsNew, lengths, -increment);   % Finds parameters at steer value
        coords(:,:,midBump,i) = coordsNew;
    end

    coordsNew = [outboard(1,:); outboard(3,:); outboard(5,:); tireInboard; tireOutboard; tireContactPt; inboard(5,:)];  % Resets starting point

    for i = midSteer+1:size(coords,4)                                                                           % Works from zero to max left steer
        [coordsNew, params(2:6,midBump,i), tire(:,:,midBump,i)] = solveSteer(coordsNew, lengths, increment);    % Finds parameters at steer value
        coords(:,:,midBump,i) = coordsNew;
    end
    
    coordsNew = [outboard(1,:); outboard(3,:); outboard(5,:); tireInboard; tireOutboard; tireContactPt; inboard(5,:)];  % Resets starting point

    for i = midBump-1:-1:1
        [coordsNew, params(2:6,i,midSteer), tire(:,:,i,midSteer)] = solveBump(coordsNew, inboard, lengths, -increment);
        coords(:,:,i,midSteer) = coordsNew;
        params(1,i,:) = (midBump-i)*-increment;

        for j = midSteer-1:-1:1
            [coordsNew, params(2:6,i,j), tire(:,:,i,j)] = solveSteer(coordsNew, lengths, -increment);
            coords(:,:,i,j) = coordsNew;
        end

        coordsNew(:,:) = coords(:,:,i,midSteer);

        for j = midSteer+1:size(coords,4)
            [coordsNew, params(2:6,i,j), tire(:,:,i,j)] = solveSteer(coordsNew, lengths, increment);
            coords(:,:,i,j) = coordsNew;
        end
        
        coordsNew(:,:) = coords(:,:,i,midSteer);
    end
    
    coordsNew = [outboard(1,:); outboard(3,:); outboard(5,:); tireInboard; tireOutboard; tireContactPt; inboard(5,:)];
    
    for i = midBump+1:size(coords,3)
        [coordsNew, params(2:6,i,midSteer), tire(:,:,i,midSteer)] = solveBump(coordsNew, inboard, lengths, increment);
        coords(:,:,i,midSteer) = coordsNew;
        params(1,i,:) = (i-midBump)*increment;
        
        for j = midSteer-1:-1:1
            [coordsNew, params(2:6,i,j), tire(:,:,i,j)] = solveSteer(coordsNew, lengths, -increment);
            coords(:,:,i,j) = coordsNew;
        end

        coordsNew(:,:) = coords(:,:,i,midSteer);

        for j = midSteer+1:size(coords,4)
            [coordsNew, params(2:6,i,j), tire(:,:,i,j)] = solveSteer(coordsNew, lengths, increment);
            coords(:,:,i,j) = coordsNew;
        end
        
        coordsNew(:,:) = coords(:,:,i,midSteer);
    end

    for i = 1:size(coords, 3)
        for j = 1:size(coords,4)
            % j = 1+size(coords, 3)-i;
            % j = midSteer;
            coordsNew(:,:) = coords(:,:,i,j);
            clf
            hold on
            % constantplane("z", tireContactPt(1,3))
            set(gca, 'YDir','reverse')
            xlim([floor(inboard(2,1)/50)*50-200,ceil(inboard(1,1)/50)*50+200])
            ylim([0,floor(outboard(3,2)/50)*50+300])
            zlim([floor(tireContactPt(3)/50)*50-200,ceil(outboard(1,3)/50)*50+200])
            view(60, 15)
            plot3([coordsNew(1,1), inboard(1,1)], [coordsNew(1,2), inboard(1,2)], [coordsNew(1,3), inboard(1,3)], '-o', 'Color', 'b', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 4)
            plot3([coordsNew(1,1), inboard(2,1)], [coordsNew(1,2), inboard(2,2)], [coordsNew(1,3), inboard(2,3)], '-o', 'Color', 'b', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 4)
            plot3([coordsNew(2,1), inboard(3,1)], [coordsNew(2,2), inboard(3,2)], [coordsNew(2,3), inboard(3,3)], '-o', 'Color', 'b', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 4)
            plot3([coordsNew(2,1), inboard(4,1)], [coordsNew(2,2), inboard(4,2)], [coordsNew(2,3), inboard(4,3)], '-o', 'Color', 'b', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 4)
            plot3([coordsNew(3,1), inboard(5,1)], [coordsNew(3,2), coordsNew(7,2)], [coordsNew(3,3), inboard(5,3)], '-o', 'Color', 'r', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 4)
            plot3([coordsNew(1,1), coordsNew(2,1)], [coordsNew(1,2), coordsNew(2,2)], [coordsNew(1,3), coordsNew(2,3)], 'Color', 'k', 'LineWidth', 4)
            plot3([coordsNew(3,1), coordsNew(2,1)], [coordsNew(3,2), coordsNew(2,2)], [coordsNew(3,3), coordsNew(2,3)], 'Color', 'k', 'LineWidth', 4)
            plot3([coordsNew(1,1), coordsNew(3,1)], [coordsNew(1,2), coordsNew(3,2)], [coordsNew(1,3), coordsNew(3,3)], 'Color', 'k', 'LineWidth', 4)
            plot3([coordsNew(4,1), coordsNew(5,1)], [coordsNew(4,2), coordsNew(5,2)], [coordsNew(4,3), coordsNew(5,3)], 'Color', 'k', 'LineWidth', 4)
            surf([tire(:,1,i,j), tire(:,7,i,j)], [tire(:,2,i,j), tire(:,8,i,j)], [tire(:,3,i,j), tire(:,9,i,j)], 'FaceColor', 'k', 'FaceAlpha', 0.5)
            surf([tire(:,1,i,j), tire(:,4,i,j)], [tire(:,2,i,j), tire(:,5,i,j)], [tire(:,3,i,j), tire(:,6,i,j)], 'FaceColor', 'k', 'FaceAlpha', 0.5)
            surf([tire(:,7,i,j), tire(:,10,i,j)], [tire(:,8,i,j), tire(:,11,i,j)], [tire(:,9,i,j), tire(:,12,i,j)], 'FaceColor', 'k', 'FaceAlpha', 0.5)
            surf([tire(:,4,i,j), tire(:,10,i,j)], [tire(:,5,i,j), tire(:,11,i,j)], [tire(:,6,i,j), tire(:,12,i,j)], 'FaceColor', 'k', 'FaceAlpha', 0.5)
            constantplane("z", coordsNew(6,3))
            % midTire = (coordsNew(4,:)+coordsNew(5,:))/2;
            % v = coordsNew(2,:)-coordsNew(1,:);
            % t = (coordsNew(6,3)-coordsNew(2,3))/v(3);
            % x = coordsNew(2,1)+t*v(1);
            % y = coordsNew(2,2)+t*v(2);
            % plot3(midTire(1), midTire(2), coordsNew(6,3), '-o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
            % plot3(x, y, coordsNew(6,3), '-o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
            % plot3([x, coordsNew(1,1)],[y, coordsNew(1,2)],[coordsNew(6,3), coordsNew(1,3)])
            % v = coordsNew(5,:)-coordsNew(4,:);
            % d1 = v(1)*coordsNew(4,1)+v(2)*coordsNew(4,2)+v(3)*coordsNew(4,3);
            % d2 = v(1)*coordsNew(5,1)+v(2)*coordsNew(5,2)+v(3)*coordsNew(5,3);
            % constantplane(v, d1)
            % constantplane(v, d2)
            drawnow
            hold off
        end
    end

    subplot(2,3,1);
    hold on
    plot(params(1,:,midSteer), params(2,:,midSteer), 'LineWidth', 2,'Color',[0, 0.4470, 0.7410],'displayName','Toe vs. Bump')
    grid on
    xlabel('Front Left Bump (in)', 'FontSize', 12)
    ylabel('Front Left Toe (degrees)', 'FontSize', 12)
    title('Toe vs. Bump', 'FontSize', 20)
    xlim([params(1,1,midSteer),params(1,size(params,2),midSteer)])
    % ylim([floor(params(2,1,midSteer)*10)/10,ceil(params(2,size(params,2),midSteer)*10)/10])
    plot(params(1,midBump,midSteer),params(2,midBump,midSteer),'-o','MarkerFaceColor','k','MarkerEdgeColor','k','Color',[0, 0.4470, 0.7410],'LineWidth',2)
    legend('','Static Position','Location','northwest')
    hold off
    subplot(2,3,2);
    hold on
    plot(params(1,:,midSteer), params(3,:,midSteer), 'LineWidth', 2,'Color',[0, 0.4470, 0.7410],'displayName','Camber vs. Bump')
    grid on
    xlabel('Front Left Bump (in)', 'FontSize', 12)
    ylabel('Front Left Camber (degrees)', 'FontSize', 12)
    title('Camber vs. Bump', 'FontSize', 20)
    xlim([params(1,1,midSteer),params(1,size(params,2),midSteer)])
    ylim([floor(params(3,size(params,2),midSteer)*10)/10,ceil(params(3,1,midSteer)*10)/10])
    plot(params(1,midBump,midSteer),params(3,midBump,midSteer),'-o','MarkerFaceColor','k','MarkerEdgeColor','k','Color',[0, 0.4470, 0.7410],'LineWidth',2)
    legend('','Static Position')
    hold off
    subplot(2,3,4);
    hold on
    plot(params(1,:,midSteer), params(4,:,midSteer), 'LineWidth', 2,'Color',[0, 0.4470, 0.7410],'displayName','Caster vs. Bump')
    grid on
    xlabel('Front Left Bump (in)', 'FontSize', 12)
    ylabel('Front Left Caster (degrees)', 'FontSize', 12)
    title('Caster vs. Bump', 'FontSize', 20)
    xlim([params(1,1,midSteer),params(1,size(params,2),midSteer)])
    ylim([floor(params(4,1,midSteer)*10)/10,ceil(params(4,size(params,2),midSteer)*10)/10])
    plot(params(1,midBump,midSteer),params(4,midBump,midSteer),'-o','MarkerFaceColor','k','MarkerEdgeColor','k','Color',[0, 0.4470, 0.7410],'LineWidth',2)
    legend('','Static Position','Location','northwest')
    hold off
    subplot(2,3,5);
    hold on
    plot(params(1,:,midSteer), params(5,:,midSteer), 'LineWidth', 2,'Color',[0, 0.4470, 0.7410],'displayName','Mechanical Trail vs. Bump')
    grid on
    xlabel('Front Left Bump (in)', 'FontSize', 12)
    ylabel('Front Left Mechanical Trail (in)', 'FontSize', 12)
    title('Mechanical Trail vs. Bump', 'FontSize', 20)
    xlim([params(1,1,midSteer),params(1,size(params,2),midSteer)])
    ylim([floor(params(5,1,midSteer)*20)/20,ceil(params(5,size(params,2),midSteer)*20)/20])
    plot(params(1,midBump,midSteer),params(5,midBump,midSteer),'-o','MarkerFaceColor','k','MarkerEdgeColor','k','Color',[0, 0.4470, 0.7410],'LineWidth',2)
    legend('','Static Position','Location','northwest')
    hold off
    subplot(2,3,3)
    hold on
    plot(squeeze(params(2,midBump,:)-params(2,midBump,midSteer)),squeeze(params(3,midBump,:)), 'LineWidth', 2,'Color',[0, 0.4470, 0.7410],'displayName','Camber vs. Steer')
    grid on
    xlabel('Front Left Steer (degrees)', 'FontSize', 12)
    ylabel('Front Left Camber (degrees)', 'FontSize', 12)
    title('Camber vs. Steer', 'FontSize', 20)
    xlim([params(2,midBump,1)-params(2,midBump,midSteer),params(2,midBump,size(params,3))-params(2,midBump,midSteer)])
    plot(0,params(3,midBump,midSteer),'-o','MarkerFaceColor','k','MarkerEdgeColor','k','Color',[0, 0.4470, 0.7410],'LineWidth',2)
    legend('','Static Position')
    hold off
    subplot(2,3,6)
    hold on
    surf(squeeze(params(1,:,:)), squeeze(params(2,:,:)-params(2,midBump,midSteer)), squeeze(params(3,:,:)))
    xlabel('Front Left Bump (in)', 'FontSize', 12)
    ylabel('Front Left Steer (degrees)', 'FontSize', 12)
    zlabel('Front Left Camber (degrees)', 'FontSize', 12)
    title('Camber vs. Steer vs. Bump', 'FontSize', 20)
    plot3(params(1,midBump,midSteer),0,params(3,midBump,midSteer),'-o','MarkerFaceColor','k','MarkerEdgeColor','k','Color',[0, 0.4470, 0.7410],'LineWidth',2)
    legend('','Static Position')
    hold off
    % ylim([floor(params(3,midBump,size(params,3))*10)/10,ceil(params(3,midBump,1)*10)/10])
end

function coords = solveSpheres(p1, p2, p3, v1, v2, v3, p0)
    coords = zeros(1,3);
    syms x y z
    a = (2*p2(1)-2*p1(1));
    b = (2*p2(2)-2*p1(2));
    c = (2*p2(3)-2*p1(3));
    d = v1^2-v2^2-p1(1)^2+p2(1)^2-p1(2)^2+p2(2)^2-p1(3)^2+p2(3)^2;
    A = (2*p3(1)-2*p1(1));
    B = (2*p3(2)-2*p1(2));
    C = (2*p3(3)-2*p1(3));
    D = v1^2-v3^2-p1(1)^2+p3(1)^2-p1(2)^2+p3(2)^2-p1(3)^2+p3(3)^2;
    A = [a b c d; A B C D];
    A = rref(A);
    eq1 = x*A(1,1)+y*A(1,2)+z*A(1,3) == A(1,4);
    eq2 = x*A(2,1)+y*A(2,2)+z*A(2,3) == A(2,4);
    x = solve(eq1, x);
    y = solve(eq2, y);
    eq3 = (x-p1(1))^2+(y-p1(2))^2+(z-p1(3))^2 == v1^2;
    zcoord = solve(eq3, z);
    z1 = double(zcoord(1));
    z2 = double(zcoord(2));
    x1 = double((A(1,4)-y*A(1,2)-z1*A(1,3))/A(1,1));
    x2 = double((A(1,4)-y*A(1,2)-z2*A(1,3))/A(1,1));
    y1 = double((A(2,4)-x*A(2,1)-z1*A(2,3))/A(2,2));
    y2 = double((A(2,4)-x*A(2,1)-z2*A(2,3))/A(2,2));
    if abs(norm([x1, y1, z1]-p0))<abs(norm([x2, y2, z2]-p0))
        coords(1) = x1;
        coords(2) = y1;
        coords(3) = z1;
    else
        coords(1) = x2;
        coords(2) = y2;
        coords(3) = z2;
    end
end

function [tireInner, tireOuter] = findTire(p1, p2, toeAngle, camberAngle, r)
    t = (0:(2*pi/100):2*pi)';
    quat = quaternion([toeAngle,camberAngle,0],"eulerd","ZXY","point");

    x1 = r*25.4*cos(t);
    y1 = 0*t;
    z1 = r*25.4*sin(t);

    x2 = r*25.4*cos(t);
    y2 = 0*t;
    z2 = r*25.4*sin(t);
    tireInner = [x1 y1 z1];
    tireOuter = [x2 y2 z2];
    
    for i = 1:length(t)
        tireInner(i,:) = rotatepoint(quat, tireInner(i,:));
        tireOuter(i,:) = rotatepoint(quat, tireOuter(i,:));
    end
    tireInner = tireInner+p1;
    tireOuter = tireOuter+p2;
end

function [coordsNew, params, tire] = solveBump(coordsOld, inboard, lengths, bump)
    coordsNew = zeros(7,3);
    params = zeros(1,5);
    tire = zeros(101,12);

    coordsNew(7,:) = coordsOld(7,:);
    coordsNew(2,1) = coordsOld(2,1);
    coordsNew(2,3) = coordsOld(2,3)+bump*25.4;                % Adds bump to lower control arm z coordinate
    coordsNew(2,2) = inboard(3,2)+sqrt(lengths(1)^2-(coordsNew(2,3)-inboard(3,3))^2);   % Finds new lower control arm y coordinate
    coordsNew(1,:) = solveSpheres(inboard(1,:), inboard(2,:), coordsNew(2,:), lengths(2), lengths(3), lengths(5), coordsOld(1,:));  % Finds new upper control arm coordinates
    coordsNew(3,:) = solveSpheres(coordsNew(1,:), coordsNew(2,:), coordsNew(7,:), lengths(6), lengths(7), lengths(4), coordsOld(3,:));
    coordsNew(4,:) = solveSpheres(coordsNew(1,:), coordsNew(2,:), coordsNew(3,:), lengths(8), lengths(9), lengths(10), coordsOld(4,:));
    coordsNew(5,:) = solveSpheres(coordsNew(1,:), coordsNew(2,:), coordsNew(3,:), lengths(11), lengths(12), lengths(13), coordsOld(5,:));
    coordsNew(6,:) = solveSpheres(coordsNew(4,:), coordsNew(5,:), coordsNew(3,:), lengths(14), lengths(15), lengths(16), coordsOld(6,:));

    params(1) = asind((coordsNew(5,1)-coordsNew(4,1))/norm((coordsNew(5,1:2)-coordsNew(4,1:2)))); % Finds toe value
    params(2) = -asind((coordsNew(5,3)-coordsNew(4,3))/norm((coordsNew(5,:)-coordsNew(4,:)))); % Finds camber value
    params(3) = asind((coordsNew(2,1)-coordsNew(1,1))/norm(coordsNew(1,[1,3])-coordsNew(2,[1,3])));      % Finds caster value

    midTire = (coordsNew(4,:)+coordsNew(5,:))/2;
    v = coordsNew(2,:)-coordsNew(1,:);
    t = (coordsNew(6,3)-coordsNew(2,3))/v(3);
    x = coordsNew(2,1)+t*v(1);
    y = coordsNew(2,2)+t*v(2);
    params(4) = (x-midTire(1))/25.4;
    params(5) = (y-midTire(2))/25.4;
    [tire(:,1:3), tire(:,7:9)] = findTire(coordsNew(4,:), coordsNew(5,:), -params(1), -params(2),8);
    [tire(:,4:6), tire(:,10:12)] = findTire(coordsNew(4,:), coordsNew(5,:), -params(1), -params(2),5);
end

function [coordsNew, params, tire] = solveSteer(coordsOld, lengths, steer)
    coordsNew = zeros(7,3);
    params = zeros(1,5);
    tire = zeros(101,12);

    coordsNew(1:2,:) = coordsOld(1:2,:);
    coordsNew(7,[1,3]) = coordsOld(7,[1,3]);
    coordsNew(7,2) = coordsOld(7,2)+steer*25.4;
    coordsNew(3,:) = solveSpheres(coordsNew(1,:), coordsNew(2,:), coordsNew(7,:), lengths(6), lengths(7), lengths(4), coordsOld(3,:));
    coordsNew(4,:) = solveSpheres(coordsNew(1,:), coordsNew(2,:), coordsNew(3,:), lengths(8), lengths(9), lengths(10), coordsOld(4,:));
    coordsNew(5,:) = solveSpheres(coordsNew(1,:), coordsNew(2,:), coordsNew(3,:), lengths(11), lengths(12), lengths(13), coordsOld(5,:));  
    coordsNew(6,:) = solveSpheres(coordsNew(4,:), coordsNew(5,:), coordsNew(3,:), lengths(14), lengths(15), lengths(16), coordsOld(6,:));

    params(1) = asind((coordsNew(5,1)-coordsNew(4,1))/norm((coordsNew(5,1:2)-coordsNew(4,1:2)))); % Finds steer angle
    params(2) = -asind((coordsNew(5,3)-coordsNew(4,3))/norm((coordsNew(5,:)-coordsNew(4,:)))); % Finds camber value
    params(3) = asind((coordsNew(2,1)-coordsNew(1,1))/norm(coordsNew(1,[1,3])-coordsNew(2,[1,3])));      % Finds caster value

    midTire = (coordsNew(4,:)+coordsNew(5,:))/2;
    v = coordsNew(2,:)-coordsNew(1,:);
    t = (coordsNew(6,3)-coordsNew(2,3))/v(3);
    x = coordsNew(2,1)+t*v(1);
    y = coordsNew(2,2)+t*v(2);
    params(4) = (x-midTire(1))/25.4;
    params(5) = (y-midTire(2))/25.4;

    [tire(:,1:3), tire(:,7:9)] = findTire(coordsNew(4,:), coordsNew(5,:), -params(1), -params(2),8);
    [tire(:,4:6), tire(:,10:12)] = findTire(coordsNew(4,:), coordsNew(5,:), -params(1), -params(2),5);
end
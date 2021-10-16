lengths = readmatrix('lengths.txt');
path = readmatrix('Path.txt');
a = [0,lengths(2),lengths(3)];
alpha = [pi/2,0,0];
d = [lengths(1),0,0];

theta = readmatrix('angles.txt');
size_theta = size(theta);
pos = zeros(3,4);
total_pos = zeros(3,4);
t_pause = 0.01;

l_axis = sum(a) + sum(d);

figure(1)
Plot = plot3(path(:,1),path(:,2),path(:,3),'r*','LineWidth',3);
axis([-2 4 -4 1 0 6])
xlabel('X axis')
ylabel('Y axis')
zlabel('Z axis')
grid on;
hold on;


for i = 1:size_theta(1)
    A1 = compute_dh_matrix(a(1),alpha(1),d(1),theta(i,1));
    A2 = compute_dh_matrix(a(2),alpha(2),d(2),theta(i,2));
    A3 = compute_dh_matrix(a(3),alpha(3),d(3),theta(i,3));

    % Transformation from link0 to link1
    T0_1 = A1;
    pos(1,2) = T0_1(1,4);
    pos(2,2) = T0_1(2,4);
    pos(3,2) = T0_1(3,4);
    
    % Transformation from link0 to link2
    T0_2 = A1*A2;
    pos(1,3) = T0_2(1,4);
    pos(2,3) = T0_2(2,4);
    pos(3,3) = T0_2(3,4);
    
    % Transformation from link0 to link3
    T0_3 = A1*A2*A3;
    pos(1,4) = T0_3(1,4);
    pos(2,4) = T0_3(2,4);
    pos(3,4) = T0_3(3,4);
    
    total_pos(:,i) = pos(:,4); 
    
    figure(1)
    Plot = plot3(pos(1,:),pos(2,:),pos(3,:),'b-*','LineWidth',2);
    %axis([-l_axis l_axis -l_axis l_axis 0 l_axis])
    %axis([-3 4 -4 1 0 5])
    xlabel('X axis')
    ylabel('Y axis')
    zlabel('Z axis')
    grid on;
    %hold on;
    pause(t_pause)
    if i ~= size_theta(1)
        delete(Plot);
    end
end

%figure(1)
%hold on
%plot3(total_pos(1,:),total_pos(2,:),total_pos(3,:),'-*','LineWidth',1);
%axis([-l_axis l_axis -l_axis l_axis 0 l_axis])
%grid on;


%plot(pos(1,:),pos(2,:))
for i = 1:size_theta(1)
    A1 = compute_dh_matrix(a(1),alpha(1),d(1),theta(i,1));
    A2 = compute_dh_matrix(a(2),alpha(2),d(2),theta(i,2));
    A3 = compute_dh_matrix(a(3),alpha(3),d(3),theta(i,3));

    % Transformation from link0 to link1
    T0_1 = A1;
    pos(1,2) = T0_1(1,4);
    pos(2,2) = T0_1(2,4);
    pos(3,2) = T0_1(3,4);
    
    % Transformation from link0 to link2
    T0_2 = A1*A2;
    pos(1,3) = T0_2(1,4);
    pos(2,3) = T0_2(2,4);
    pos(3,3) = T0_2(3,4);
    
    % Transformation from link0 to link3
    T0_3 = A1*A2*A3;
    pos(1,4) = T0_3(1,4);
    pos(2,4) = T0_3(2,4);
    pos(3,4) = T0_3(3,4);
    
    figure(2)
    xlabel('X axis')
    ylabel('Y axis')
    zlabel('Z axis')
    Plot = plot3(pos(1,:),pos(2,:),pos(3,:),'-*','LineWidth',2);
    %Plot = plot3(pos(1,4),pos(2,4),pos(3,4),'','LineWidth',2);
    axis([-2 4 -4 1 0 6])
    grid on;
    hold on;
    pause(t_pause)
end
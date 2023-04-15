

nodes = [
    -1,0;
    -1,1;
    1,1;
    1,0;
];

traj = CardinalSpline(nodes, 0.5, 1e-3);
traj.get_distance_sum();

hold on
point_num = 400;
l = linspace(0, size(nodes,1)-1, point_num);


p = traj.get_positions(l);
disp(traj.get_distance_sum())

axis equal
plot(p(:,1), p(:,2));
%plot(p(:,1), p(:,2), ".");
plot(nodes(:,1), nodes(:,2), "o");


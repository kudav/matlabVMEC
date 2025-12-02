function intersections = calculateIntersections(lens, axis, phi)
% numLines = size(lens, 2);
% intersections = zeros(3,numLines); % Prepare for (x, y, z) points
% % X value for all intersections based on fixed phi
% x_fixed = radius' .* cos(phi);



t = (tan(phi).*lens(1,:)-lens(2,:))./(axis(2,:)-axis(1,:).*tan(phi));
intersections = lens + t .* axis;
% for i = 1:numLines
%     t = (x_fixed - lens(1, i)) / axis(1, i);
%     intersections(:, i) = lens(:, i) + t(i) .* axis(:, i);
% end
intersections(1,:)=sqrt(intersections(1,:).^2+intersections(2,:).^2);
intersections=intersections([1,3],:); %Only R,Z coordinates
% Filter out any intersections that do not make physical sense, e.g., if t < 0
intersections = intersections(:,t >= 0);
end
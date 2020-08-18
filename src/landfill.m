function [ landfilled_hist ] = landfill( joint_dist)

landfilled_hist = joint_dist;
x = size(landfilled_hist, 2); 
y = size(landfilled_hist, 1);
for i = 1:y 
  for j = 1:x 
    landfilled_hist(i,j) = max(max(joint_dist(i:y, j:x)));
  end 
end


% Function to view slices of three dimensional data sets
%
% Usage: view(dat,S)
%
% dat: 3d data set (any shape) of total size prod(S)=S(1)*S(2)*S(3)
% S: dimensions of dat in a 3-vector

function view(dat,S)

fprintf('\nRemember to hit <enter> or <spacebar> after each plot!\n\n');

for k=1:3
  if k==1
    fprintf('m1=0 slice (m3, m2 along left-, right- axes):\n');
  elseif k==2
    fprintf('m2=0 slice (m3, m1 along left-, right- axes):\n');
  elseif k==3
    fprintf('m3=0 slice (m2, m1 along left-, right- axes):\n');
  end

  mesh(slice(dat,S,1,k)); pause;
end

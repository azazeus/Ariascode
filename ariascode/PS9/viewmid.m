% Function to view slices through center point of three dimensional data sets
% 
% Usage: viewmid(dat,S)
%
% dat: 3d data set (any shape) of total size prod(S)=S(1)*S(2)*S(3)
% S: dimensions of dat in a 3-vector

function viewmid(dat,S)

  fprintf('\nRemember to hit <enter> or <spacebar> after each \
	  plot!\n\n');

  for k=1:3
    if k==1
      fprintf('m1=%d slice (m2, m3 along left, right \
			    axes):\n',S(1)/2-1);
    elseif k==2
      fprintf('m2=%d slice (m1, m3 along left, right \
			    axes):\n',S(2)/2-1);
    elseif k==3
      fprintf('m3=%d slice (m1, m2 along left, right \
			    axes):\n',S(3)/2-1);
    end
    
    mesh( real( slice(dat,S,S(k)/2,k) )' ); pause;
  end
endfunction


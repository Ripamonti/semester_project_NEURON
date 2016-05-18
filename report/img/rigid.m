function dy = rigid(t,y)
dy = zeros(2,1);    % a column vector
dy(1) = y(1) * ( y(2)-2 );
dy(2) = y(2) * ( 1-y(1) );
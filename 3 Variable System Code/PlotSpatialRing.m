%function to generate visual plot for spatial ring case
function PlotSpatialRing(n,t,a_mat)
num_seg_pts = 100; 
theta = linspace(0,2*pi, num_seg_pts*n); 

%using sine waves: 
x = cos(theta);
y = sin(theta);
z = zeros(1,length(theta)); 
t_idx_inc = 80; %used for plotting. increase t_idx_inc to plot at fewer time steps and decrease function run time.


    for i = 1:t_idx_inc:round(length(a_mat),-2)
        A = a_mat(i,:); 

        for j = 0:(n-1)
            idx1 = j*num_seg_pts+1;
            idx2 = num_seg_pts/2 + idx1 - 1;
            %z is a matrix of superimposed sine waves corresponding to the artist densities in each region at each time step:
            z(idx1:idx2) = A(j+1)*sin(n*theta(idx1:idx2));    
        end
        figure(2)
        view(3)
        plot3(x,y,z)
        zlim([0 1])
    end
end
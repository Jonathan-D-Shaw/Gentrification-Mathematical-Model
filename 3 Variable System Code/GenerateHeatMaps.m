%% Heat Map Generation
%housekeeping: 
close all; clear; clc; 
load("heat_map_data.mat")

n_vec = (4:1:10)'; 
z = 0.01; 
epsilon = 0.1; 
tau_n = 0.1:0.001:20;

%generating heat maps: 
for it = 1:length(n_vec)
    n = n_vec(it); 
    sig_p = sigma_prime(1/n,z,epsilon);
    tau_n_trunc = tau_n(tau_n>(1/sig_p)); 
    tau_c = 1./(sig_p - 1./tau_n_trunc); 

    %plotting heat map:
    figure(2*it-1)
    hold on      
    tau_n_map = heat_map_data{it}(:,4); 
    tau_c_map = heat_map_data{it}(:,5); 
    var = heat_map_data{it}(:,6); var(var<0) = NaN;
    dt = delaunayTriangulation(tau_n_map,tau_c_map) ;
    tri = dt.ConnectivityList ;
    xi = dt.Points(:,1) ; 
    yi = dt.Points(:,2) ; 
    F = scatteredInterpolant(tau_n_map,tau_c_map,var);
    zi = F(xi,yi) ;
    trisurf(tri,xi,yi,zi) 
    view(2)
    shading flat
    c = colorbar; 
    clim([1e-10 max(var)])
    c.Label.String = '\lambda_R';
    c.Label.FontSize = 12;
    zvec = ones(length(tau_n_trunc),1);
    plot3(tau_n_trunc,tau_c,zvec, 'r','linewidth',2) %plotting stability condition over surfac
    ylim([0 20])
    xlabel('\tau_N', 'fontsize', 12)
    ylabel('\tau_C', 'fontsize', 12)
    l = strcat('Heat Map, n = ', num2str(n), ', z = ', num2str(z), ', \epsilon = ', num2str(epsilon));
    title(l, 'fontweight', 'bold', 'fontsize', 12)
    set(gca,'Color','k')
    hold off  
    F = getframe(gcf);
    im = frame2im(F);
    imwrite(im,sprintf('neq%d_heat_map.png',n));
   

    %plotting binary map: 
    figure(2*it)
    hold on
    increase = heat_map_data{it}(:,7);         
    F = scatteredInterpolant(tau_n_map,tau_c_map,increase);
    zi = F(xi,yi) ; 
    trisurf(tri,xi,yi,zi) 
    view(2)
    shading interp  
    plot3(tau_n_trunc,tau_c,zvec, 'r','linewidth',2) %plotting stability condition over surface
    ylim([0 20])
    %xline(1/sig_p,'r','linewidth',2)
    xlabel('\tau_N', 'fontsize', 12)
    ylabel('\tau_C', 'fontsize', 12)
    l = strcat('Binary Heat Map, n = ', num2str(n), ', z = ', num2str(z), ', \epsilon = ', num2str(epsilon));
    title(l, 'fontweight', 'bold', 'fontsize', 12)
    hold off
    saveas(gcf,sprintf('neq%d_binary.png',n))
    
end
   
%function to determine whether interesting behavior is present:
function [peak_data,ind_indicators,grp_indicators, fixed_pt, base_osc, grp_osc, sep_grp_osc, chaos] = InterestingBehavior(n,t,a_mat,params)

fixed_pt = 0; 
base_osc = 0; 
grp_osc = 0;
sep_grp_osc = 0; 
chaos = 0; 
num_peaks_in_reg = zeros(n,1); 
peak_data = [];
ind_indicators = [];


%finding peaks and performing individual region analysis:
    for i = 1:n
        [peak_vals{i},peak_idxs{i}] = findpeaks(a_mat(:,i),'MinPeakHeight',1/n); %pulling peak values and peak indices for each region
        
        num_peaks_in_reg(i) = length(peak_vals{i});                         %number of peaks in region i 
        
        t_vec = t(peak_idxs{i}); 
        A_vec = peak_vals{i};
        
        %discarding fixed-point case: 
        if isempty(peak_vals{i}) || mean(peak_vals{i})<(params.A_abs_tol + 1/n)
            ind_indicators = 0; 
            fixed_pt = 1; 
            break
        end
        
        peak_data = [peak_data; t_vec A_vec];                               %concatenating peak data across all regions
                
        %individual region period analysis: 
        T_theory = (t_vec(end)-t_vec(1))/(num_peaks_in_reg(i) -1);          %theoretical period of oscillations if base case is occuring
        t_tol = 0.1*T_theory;                                              %period proximity tolerance proportional to theoretical period
        T_actual = diff(t_vec);                                             %actual difference between time locations of peaks          
        T_dif = abs(T_actual - T_theory);                                   %difference between actual and theoretical periods
        T_ind_vec = T_dif(T_dif>t_tol);                                     %indicator vector. 1 indicates the absolute value of the difference between two actual points and the theoretical period has been exceed, per t_tol.
                                                 %test vector becomes logical array with 1 entries if difference has exceeded tolerance

        %individual region amplitude analysis: 
        A_ref = A_vec(1);                                            %reference peak value for each region
        A_dif = abs(A_vec - A_ref);                                  %difference between all peak values and reference value
        A_ind_vec = A_dif(A_dif>params.A_rel_tol);                              %same test vector idea as before      


        ind_indicators = [ind_indicators; sum(T_ind_vec)~=0 sum(A_ind_vec)~=0]; %matrix of individual region indicators in the form (period violation, peak violation). 1 will be entered if a tolerance has been exceeded.
    end

%group analysis: 
    if fixed_pt == 1                                                        %discarding fixed-point case:
        grp_indicators = zeros(1,2); 
    else
        sorted_peak_data = sortrows(peak_data);                             %putting peaks in chronological order
       
        t_vec  = sorted_peak_data(:,1);                                     %time vector for all peaks in chronological order
        A_vec = sorted_peak_data(:,2);                                      %amplitude vector for all peaks in chronological order
        tot_num_peaks = length(sorted_peak_data);                               %total number of peaks
        
        %group period analysis: 
        T_theory = (t_vec(end) - t_vec(1))/(tot_num_peaks - 1); %theoretical time spacing between peaks should be linear
        T_actual = diff(t_vec);                                             %actual difference between time locations of peaks          
        T_dif = abs(T_actual - T_theory);                                   %difference between actual and theoretical periods
        T_ind_vec = T_dif(T_dif>t_tol);                                     %indicator vector. 1 indicates the absolute value of the difference between two actual points and the theoretical period has been exceed, per t_tol.
                                                 

        %individual region amplitude analysis: 
        A_ref = A_vec(1);                                                   %reference peak value for each region
        A_dif = abs(A_vec - A_ref);                                         %difference between all peak values and reference value
        A_ind_vec = A_dif(A_dif>params.A_rel_tol);                           

        grp_indicators = [sum(T_ind_vec)~=0, sum(A_ind_vec)~=0];

   
    
        %categorizing interesting oscillatory behaviors: 

        groups = dbscan(sorted_peak_data, t_tol,2);                                 %looking for groups of size 2 inside a radius of t_tol (theory explained in the notes document)
        unique_groups = unique(groups);                                             %vector in the form 1:u where u is the number of unique groups
        ind_sep_grps = zeros(length(unique_groups),1);                              %indicators for whether amplitudes in a group are constant
        num_peaks_not_in_grp = tot_num_peaks; 

        %determining if the amplitude of oscillations in each group is constant: 
            for i = 1:length(unique_groups)

                if unique_groups(i) ~= -1                                           %discarding degenerative case, means a peaks doesn't belong to a group
                group_peaks = nonzeros(sorted_peak_data(:,2).*(groups == unique_groups(i))); %pulls peak amplitudes in group i
                num_peaks_not_in_grp = num_peaks_not_in_grp - length(group_peaks); 

                if abs(max(group_peaks)-min(group_peaks))>params.A_rel_tol
                    ind_sep_grps(i) = 1;                                            %indicator set to 1 if the difference 
                end
                end
            end
    end
  
%categorizing behavior: 
    if fixed_pt==1
        
    elseif (sum(ind_indicators,'all') + sum(grp_indicators))==0
        base_osc = 1; 
        
    elseif (num_peaks_not_in_grp<floor(0.1*tot_num_peaks)) && (sum(grp_indicators)>0) && (sum(ind_sep_grps(i))==0) %arguement is that if the number of peaks not in a group is less than 10% of the total number of peaks, there is grouping
        grp_osc = 1; 
        
    elseif (num_peaks_not_in_grp<floor(0.1*tot_num_peaks)) && (sum(grp_indicators)>0) && (sum(ind_indicators,'all')>0)  
        sep_grp_osc = 1; 
        
    else
        chaos = 1;    
    end
end
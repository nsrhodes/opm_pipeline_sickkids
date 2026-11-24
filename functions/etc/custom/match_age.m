function idx_closest = match_age(age, matching_group, n)

age_bin_min = floor(age);
age_bin_max = ceil(age);
 
% Sort subjects based on age
age = age*ones(size(matching_group,1), 1);
age_diff = abs(age - matching_group.age);
[~, idx_sorted] = sort(age_diff);
matching_group_sorted = matching_group(idx_sorted, :);

% Ensure they are in the same age bin
idx_sorted = idx_sorted(age_bin_min <= matching_group_sorted.age & matching_group_sorted.age < age_bin_max);

if n <= length(idx_sorted)
    idx_closest = idx_sorted(1:n);
else
    warning('N = %i but only %i potential matches found',n,length(idx_sorted));
    idx_closest = idx_sorted;
end
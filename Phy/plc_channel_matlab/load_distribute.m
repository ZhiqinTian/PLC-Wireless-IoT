function level_load_type = load_distribute(load_type)
% This function generate the type of load for each cables

%%%%%%%%% the following is considered as a uniformly distributed variable
%%%%%%%%% which means that the load type is between 1 to load_type
level_load_type = randperm(load_type);
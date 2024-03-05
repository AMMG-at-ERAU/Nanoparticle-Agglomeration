clc
clear all

A = [0 2:2:20]
load('CNT_Hole_Cond_2D_Vf_8p_agg_0p_Hole_rad_Lover20_number_(2_2_20)_5000.mat')
result = (result1-result1(:,1))./result1(:,1);
Result = mean(result);
plot(A,Result,'*-g')
hold on
load('CNT_Hole_Cond_2D_Vf_10p_agg_0p_Hole_rad_Lover20_number_(2_2_20)_5000.mat')
result = (result1-result1(:,1))./result1(:,1);
Result = mean(result);
plot(A,Result,'*-k')
hold on
load('CNT_Hole_Cond_2D_Vf_15p_agg_0p_Hole_rad_Lover20_number_(2_2_20)_5000.mat')
result = (result1-result1(:,1))./result1(:,1);
Result = mean(result);
plot(A,Result,'*-b')
hold on
load('CNT_Hole_Cond_2D_Vf_20p_agg_0p_Hole_rad_Lover20_number_(2_2_20)_5000.mat')
result = (result1-result1(:,1))./result1(:,1);
Result = mean(result);
plot(A,Result,'*-r')
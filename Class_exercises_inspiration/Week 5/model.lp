maximize
obj: 20 p_d + 6.3 p_up_1_ - 12.15 p_dw_1_ + 7.700000000000001 p_up_2_ - 14.850000000000001 p_dw_2_
subject to
max_production: 1 p_d <= 100
deviation_1_: -1 p_d - 1 p_up_1_ + 1 p_dw_1_ = -82.5
deviation_2_: -1 p_d - 1 p_up_2_ + 1 p_dw_2_ = -24
Bounds
p_d >= 0
p_up_1_ >= 0
p_up_2_ >= 0
p_dw_1_ >= 0
p_dw_2_ >= 0
End

V30 :0x14 grid_setup_mod
18 grid_setup_mod.f90 S624 0
07/30/2017  20:53:12
use pgi_acc_common private
use iso_c_binding private
enduse
D 56 24 646 8 645 7
D 65 24 649 8 648 7
D 74 24 646 8 645 7
D 95 24 722 8 721 7
D 4359 21 9 2 5718 5905 0 0 0 0 0
 0 13 3 3 13 13
 0 5799 13 3 5799 5799
D 4362 21 9 2 5718 5905 0 0 0 0 0
 0 13 3 3 13 13
 0 5799 13 3 5799 5799
D 4365 21 9 2 5718 5799 0 0 0 0 0
 0 13 3 3 13 13
 0 5710 13 3 5710 5710
D 4368 21 9 1 3 6599 0 0 0 0 0
 0 6599 3 3 6599 6599
D 4371 21 9 1 3 6599 0 0 0 0 0
 0 6599 3 3 6599 6599
D 4374 21 9 1 3 5773 0 0 0 0 0
 0 5773 3 3 5773 5773
D 4377 21 9 1 3 6599 0 0 0 0 0
 0 6599 3 3 6599 6599
D 4380 21 9 1 3 6599 0 0 0 0 0
 0 6599 3 3 6599 6599
D 4383 21 9 1 3 5773 0 0 0 0 0
 0 5773 3 3 5773 5773
D 4386 21 9 1 3 13 0 0 0 0 0
 0 13 3 3 13 13
D 4389 21 9 2 5718 7012 0 0 1 0 0
 0 13 3 3 13 13
 0 7010 13 3 7011 7011
D 4392 21 9 1 3 7014 0 0 1 0 0
 0 7013 3 3 7014 7014
S 624 24 0 0 0 9 1 0 5011 10015 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 23 0 0 0 0 0 0 grid_setup_mod
S 629 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 645 25 6 iso_c_binding c_ptr
R 646 5 7 iso_c_binding val c_ptr
R 648 25 9 iso_c_binding c_funptr
R 649 5 10 iso_c_binding val c_funptr
R 683 6 44 iso_c_binding c_null_ptr$ac
R 685 6 46 iso_c_binding c_null_funptr$ac
R 686 26 47 iso_c_binding ==
R 688 26 49 iso_c_binding !=
R 721 25 6 pgi_acc_common c_devptr
R 722 5 7 pgi_acc_common cptr c_devptr
R 725 6 10 pgi_acc_common c_null_devptr$ac
R 729 26 14 pgi_acc_common =
S 5578 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 5585 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 5594 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 49 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 5607 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 64 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 5624 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 256 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 10784 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 193 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 11219 27 0 0 0 6 11230 624 73575 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 make_grid
S 11220 27 0 0 0 9 11232 624 73585 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 get_visualization_grid
S 11221 7 4 0 4 4359 11222 624 73608 800014 100 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 11228 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 global_x1d
S 11222 7 4 0 4 4362 11223 624 73619 800014 100 A 0 0 0 0 B 0 0 0 0 0 2048 0 0 0 0 0 0 11228 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 global_y1d
S 11223 7 4 0 4 4365 11224 624 73630 800014 100 A 0 0 0 0 B 0 0 0 0 0 4096 0 0 0 0 0 0 11228 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 global_z1d
S 11224 7 4 0 4 4368 11225 624 73641 800014 100 A 0 0 0 0 B 0 0 0 0 0 4608 0 0 0 0 0 0 11228 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 viz_x1d
S 11225 7 4 0 4 4371 11226 624 73649 800014 100 A 0 0 0 0 B 0 0 0 0 0 6160 0 0 0 0 0 0 11228 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 viz_y1d
S 11226 7 4 0 4 4374 1 624 73657 800014 100 A 0 0 0 0 B 0 0 0 0 0 7712 0 0 0 0 0 0 11228 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 viz_z1d
S 11227 6 4 0 0 16 1 624 73665 80001c 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 11229 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 initialize_1d_grids_called
S 11228 11 0 0 4 9 11206 624 73692 40800010 805000 A 0 0 0 0 B 0 0 0 0 0 8104 0 0 11221 11226 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _grid_setup_mod$6
S 11229 11 0 0 0 9 11228 624 73710 40800010 805000 A 0 0 0 0 B 0 0 0 0 0 4 0 0 11227 11227 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _grid_setup_mod$12
S 11230 23 5 0 0 0 11231 624 73575 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 make_grid
S 11231 14 5 0 0 0 1 11230 73575 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4837 0 0 0 0 0 0 0 0 0 0 0 0 0 48 0 624 0 0 0 0 make_grid
F 11231 0
S 11232 23 5 0 0 0 11236 624 73585 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_visualization_grid
S 11233 7 3 2 0 4377 1 11232 73729 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ax
S 11234 7 3 2 0 4380 1 11232 73732 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ay
S 11235 7 3 2 0 4383 1 11232 73735 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 az
S 11236 14 5 0 0 0 1 11232 73585 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4838 3 0 0 0 0 0 0 0 0 0 0 0 0 86 0 624 0 0 0 0 get_visualization_grid
F 11236 3 11233 11234 11235
S 11237 23 5 0 0 0 11238 624 73738 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 initialize_1d_grids
S 11238 14 5 0 0 0 1 11237 73738 10 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4842 0 0 0 0 0 0 0 0 0 0 0 0 0 101 0 624 0 0 0 0 initialize_1d_grids
F 11238 0
S 11239 23 5 0 0 0 11248 624 73758 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 grid_1dwork
S 11240 1 3 1 0 9 1 11239 73770 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 smin
S 11241 1 3 1 0 9 1 11239 73775 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 smax
S 11242 6 3 1 0 6 1 11239 73780 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nes
S 11243 6 3 1 0 6 1 11239 73784 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nus
S 11244 1 3 2 0 9 1 11239 73788 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dels
S 11245 7 3 2 0 4389 1 11239 73793 800214 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 global_s1d
S 11246 7 3 2 0 4392 1 11239 73804 800214 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 viz_s1d
S 11247 7 3 1 0 4386 1 11239 70166 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gllp
S 11248 14 5 0 0 0 1 11239 73758 210 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4843 8 0 0 0 0 0 0 0 0 0 0 0 0 137 0 624 0 0 0 0 grid_1dwork
F 11248 8 11240 11241 11242 11243 11247 11244 11245 11246
S 11249 6 1 0 0 6 1 11239 73812 40800016 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_7010
S 11250 6 1 0 0 6 1 11239 73821 40800016 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_7012
S 11251 6 1 0 0 6 1 11239 73830 40800016 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_7014
A 13 2 0 0 0 6 629 0 0 0 13 0 0 0 0 0 0 0 0 0 0
A 67 1 0 0 0 56 683 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 70 1 0 0 0 65 685 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 86 1 0 0 0 95 725 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 5710 2 0 0 4579 6 5585 0 0 0 5710 0 0 0 0 0 0 0 0 0 0
A 5718 2 0 0 4572 6 5578 0 0 0 5718 0 0 0 0 0 0 0 0 0 0
A 5773 2 0 0 5244 6 5594 0 0 0 5773 0 0 0 0 0 0 0 0 0 0
A 5799 2 0 0 5174 6 5607 0 0 0 5799 0 0 0 0 0 0 0 0 0 0
A 5905 2 0 0 5851 6 5624 0 0 0 5905 0 0 0 0 0 0 0 0 0 0
A 6599 2 0 0 4809 6 10784 0 0 0 6599 0 0 0 0 0 0 0 0 0 0
A 7010 1 0 0 5818 6 11242 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 7011 1 0 0 5831 6 11249 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 7012 1 0 0 6988 6 11250 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 7013 1 0 0 6147 6 11243 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 7014 1 0 0 6149 6 11251 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
J 149 1 1
V 67 56 7 0
S 0 56 0 0 0
A 0 6 0 0 1 2 0
J 150 1 1
V 70 65 7 0
S 0 65 0 0 0
A 0 6 0 0 1 2 0
J 31 1 1
V 86 95 7 0
S 0 95 0 0 0
A 0 74 0 0 1 67 0
Z

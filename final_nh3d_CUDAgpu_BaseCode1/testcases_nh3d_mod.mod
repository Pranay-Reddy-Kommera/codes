V30 :0x14 testcases_nh3d_mod
22 testcases_nh3d_mod.f90 S624 0
07/30/2017  20:53:12
use pgi_acc_common private
use iso_c_binding private
enduse
D 56 24 647 8 646 7
D 65 24 650 8 649 7
D 74 24 647 8 646 7
D 95 24 723 8 722 7
D 4335 21 9 3 88 5799 0 0 0 0 0
 0 13 3 3 13 13
 0 13 13 3 13 13
 0 13 5710 3 13 13
D 4338 21 9 3 88 5799 0 0 0 0 0
 0 13 3 3 13 13
 0 13 13 3 13 13
 0 13 5710 3 13 13
D 4341 21 9 3 88 5799 0 0 0 0 0
 0 13 3 3 13 13
 0 13 13 3 13 13
 0 13 5710 3 13 13
D 4344 21 9 3 88 5799 0 0 0 0 0
 0 13 3 3 13 13
 0 13 13 3 13 13
 0 13 5710 3 13 13
D 4347 21 9 3 88 5799 0 0 0 0 0
 0 13 3 3 13 13
 0 13 13 3 13 13
 0 13 5710 3 13 13
D 4350 21 9 3 88 5799 0 0 0 0 0
 0 13 3 3 13 13
 0 13 13 3 13 13
 0 13 5710 3 13 13
D 4353 21 9 3 88 5799 0 0 0 0 0
 0 13 3 3 13 13
 0 13 13 3 13 13
 0 13 5710 3 13 13
S 624 24 0 0 0 9 1 0 5011 10015 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 12 0 0 0 0 0 0 testcases_nh3d_mod
S 630 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 646 25 6 iso_c_binding c_ptr
R 647 5 7 iso_c_binding val c_ptr
R 649 25 9 iso_c_binding c_funptr
R 650 5 10 iso_c_binding val c_funptr
R 684 6 44 iso_c_binding c_null_ptr$ac
R 686 6 46 iso_c_binding c_null_funptr$ac
R 687 26 47 iso_c_binding ==
R 689 26 49 iso_c_binding !=
R 722 25 6 pgi_acc_common c_devptr
R 723 5 7 pgi_acc_common cptr c_devptr
R 726 6 10 pgi_acc_common c_null_devptr$ac
R 730 26 14 pgi_acc_common =
S 781 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 21 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 5585 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 5607 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 64 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 11223 27 0 0 0 6 11237 624 73559 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 load_params
S 11224 27 0 0 0 6 11239 624 73571 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 initialize_testcase
S 11225 27 0 0 0 9 11235 624 73591 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 print_vital_stat
S 11226 27 0 0 0 9 11287 624 73608 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 compute_pressure
S 11227 27 0 0 0 9 11282 624 73625 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 compute_soundspeed
S 11228 27 0 0 0 9 11297 624 73644 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 pressure
S 11229 27 0 0 0 9 11292 624 73653 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 soundspeed
S 11230 27 0 0 0 6 11257 624 73664 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 nh3d_warm_bubble
S 11231 27 0 0 0 6 11241 624 73681 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 nh3d_igw
S 11232 27 0 0 0 6 11259 624 73690 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 load_nhvars
S 11233 27 0 0 0 6 11245 624 73702 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 iterative_z2eta
S 11234 27 0 0 0 6 11243 624 73718 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 nh3d_steady_state
S 11235 23 5 0 0 0 11236 624 73591 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 print_vital_stat
S 11236 14 5 0 0 0 1 11235 73591 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4821 0 0 0 0 0 0 0 0 0 0 0 0 0 29 0 624 0 0 0 0 print_vital_stat
F 11236 0
S 11237 23 5 0 0 0 11238 624 73559 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 load_params
S 11238 14 5 0 0 0 1 11237 73559 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4822 0 0 0 0 0 0 0 0 0 0 0 0 0 58 0 624 0 0 0 0 load_params
F 11238 0
S 11239 23 5 0 0 0 11240 624 73571 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 initialize_testcase
S 11240 14 5 0 0 0 1 11239 73571 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4823 0 0 0 0 0 0 0 0 0 0 0 0 0 94 0 624 0 0 0 0 initialize_testcase
F 11240 0
S 11241 23 5 0 0 0 11242 624 73681 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nh3d_igw
S 11242 14 5 0 0 0 1 11241 73681 10 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4824 0 0 0 0 0 0 0 0 0 0 0 0 0 116 0 624 0 0 0 0 nh3d_igw
F 11242 0
S 11243 23 5 0 0 0 11244 624 73718 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nh3d_steady_state
S 11244 14 5 0 0 0 1 11243 73718 10 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4825 0 0 0 0 0 0 0 0 0 0 0 0 0 188 0 624 0 0 0 0 nh3d_steady_state
F 11244 0
S 11245 23 5 0 0 0 11256 624 73702 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 iterative_z2eta
S 11246 1 3 1 0 6 1 11245 73736 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ie
S 11247 1 3 1 0 6 1 11245 73739 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 je
S 11248 1 3 1 0 6 1 11245 73742 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ke
S 11249 7 3 2 0 4335 1 11245 71373 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 u
S 11250 7 3 2 0 4338 1 11245 71353 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rho
S 11251 7 3 2 0 4341 1 11245 71357 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 the
S 11252 7 3 2 0 4350 1 11245 71413 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rhb
S 11253 7 3 2 0 4353 1 11245 71417 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 thb
S 11254 7 3 2 0 4347 1 11245 71425 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 prb
S 11255 7 3 2 0 4344 1 11245 71429 800014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 fcori
S 11256 14 5 0 0 0 1 11245 73702 10 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4826 10 0 0 0 0 0 0 0 0 0 0 0 0 249 0 624 0 0 0 0 iterative_z2eta
F 11256 10 11246 11247 11248 11249 11250 11251 11252 11253 11254 11255
S 11257 23 5 0 0 0 11258 624 73664 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nh3d_warm_bubble
S 11258 14 5 0 0 0 1 11257 73664 10 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4837 0 0 0 0 0 0 0 0 0 0 0 0 0 340 0 624 0 0 0 0 nh3d_warm_bubble
F 11258 0
S 11259 23 5 0 0 0 11281 624 73690 10 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 load_nhvars
S 11260 1 3 1 0 6 1 11259 73736 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ie
S 11261 1 3 1 0 6 1 11259 73739 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 je
S 11262 1 3 1 0 6 1 11259 73742 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ke
S 11263 1 3 1 0 6 1 11259 6070 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 i
S 11264 1 3 1 0 6 1 11259 73745 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 j
S 11265 1 3 1 0 6 1 11259 73747 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 k
S 11266 1 3 1 0 9 1 11259 71373 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 u
S 11267 1 3 1 0 9 1 11259 71375 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 v
S 11268 1 3 1 0 9 1 11259 51750 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 w
S 11269 1 3 1 0 9 1 11259 71357 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 the
S 11270 1 3 1 0 9 1 11259 71417 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 thb
S 11271 1 3 1 0 9 1 11259 71388 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 thp
S 11272 1 3 1 0 9 1 11259 71361 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 prs
S 11273 1 3 1 0 9 1 11259 71425 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 prb
S 11274 1 3 1 0 9 1 11259 71365 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 spd
S 11275 1 3 1 0 9 1 11259 71353 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rho
S 11276 1 3 1 0 9 1 11259 71413 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rhb
S 11277 1 3 1 0 9 1 11259 71384 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rhp
S 11278 1 3 1 0 9 1 11259 71421 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rtb
S 11279 1 3 1 0 9 1 11259 73749 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rtp
S 11280 1 3 1 0 9 1 11259 73753 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cori
S 11281 14 5 0 0 0 1 11259 73690 10 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4838 21 0 0 0 0 0 0 0 0 0 0 0 0 412 0 624 0 0 0 0 load_nhvars
F 11281 21 11260 11261 11262 11263 11264 11265 11266 11267 11268 11269 11270 11271 11272 11273 11274 11275 11276 11277 11278 11279 11280
S 11282 23 5 0 0 9 11286 624 73625 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 compute_soundspeed
S 11283 1 3 0 0 9 1 11282 71361 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 prs
S 11284 1 3 0 0 9 1 11282 71353 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rho
S 11285 1 3 0 0 9 1 11282 73758 14 1003000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 speed
S 11286 14 5 0 0 9 1 11282 73625 4 1400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4860 2 0 0 11285 0 0 0 0 0 0 0 0 0 449 0 624 0 0 0 0 compute_soundspeed
F 11286 2 11283 11284
S 11287 23 5 0 0 9 11291 624 73608 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 compute_pressure
S 11288 1 3 0 0 9 1 11287 71353 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rho
S 11289 1 3 0 0 9 1 11287 71357 14 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 the
S 11290 1 3 0 0 9 1 11287 71361 14 1003000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 prs
S 11291 14 5 0 0 9 1 11287 73608 4 1400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4863 2 0 0 11290 0 0 0 0 0 0 0 0 0 458 0 624 0 0 0 0 compute_pressure
F 11291 2 11288 11289
S 11292 23 5 0 2 0 11296 624 73653 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 soundspeed
S 11293 1 3 1 0 9 1 11292 71361 8014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 prs
S 11294 1 3 1 0 9 1 11292 71353 8014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rho
S 11295 1 3 2 0 9 1 11292 73764 8014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 speed
S 11296 14 5 0 2 0 1 11292 73653 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4866 3 0 0 0 0 0 0 0 0 0 0 0 0 475 0 624 0 0 0 0 soundspeed
F 11296 3 11293 11294 11295
S 11297 23 5 0 2 0 11301 624 73644 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 pressure
S 11298 1 3 1 0 9 1 11297 71353 8014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rho
S 11299 1 3 1 0 9 1 11297 71357 8014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 the
S 11300 1 3 2 0 9 1 11297 71361 8014 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 prs
S 11301 14 5 0 2 0 1 11297 73644 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4870 3 0 0 0 0 0 0 0 0 0 0 0 0 487 0 624 0 0 0 0 pressure
F 11301 3 11298 11299 11300
A 13 2 0 0 0 6 630 0 0 0 13 0 0 0 0 0 0 0 0 0 0
A 67 1 0 0 0 56 684 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 70 1 0 0 0 65 686 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 86 1 0 0 0 95 726 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 88 2 0 0 0 6 781 0 0 0 88 0 0 0 0 0 0 0 0 0 0
A 5710 2 0 0 4578 6 5585 0 0 0 5710 0 0 0 0 0 0 0 0 0 0
A 5799 2 0 0 5174 6 5607 0 0 0 5799 0 0 0 0 0 0 0 0 0 0
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
V30 :0x14 physical_const_mod
22 physical_const_mod.f90 S624 0
07/30/2017  20:53:09
use iso_c_binding public 0 indirect
use pgi_acc_common public 0 indirect
use cudafor_lib public 0 indirect
use cudafor public 0 indirect
use basic_mod public 0 direct
enduse
D 56 24 644 8 643 7
D 65 24 647 8 646 7
D 74 24 644 8 643 7
D 95 24 720 8 719 7
S 624 24 0 0 0 9 1 0 5011 10005 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 29 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 physical_const_mod
R 643 25 6 iso_c_binding c_ptr
R 644 5 7 iso_c_binding val c_ptr
R 646 25 9 iso_c_binding c_funptr
R 647 5 10 iso_c_binding val c_funptr
R 681 6 44 iso_c_binding c_null_ptr$ac
R 683 6 46 iso_c_binding c_null_funptr$ac
R 684 26 47 iso_c_binding ==
R 686 26 49 iso_c_binding !=
R 719 25 6 pgi_acc_common c_devptr
R 720 5 7 pgi_acc_common cptr c_devptr
R 723 6 10 pgi_acc_common c_null_devptr$ac
R 727 26 14 pgi_acc_common =
S 10931 16 0 0 0 9 10933 624 71312 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10932 6731 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 grv
S 10932 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1076076216 1374389535 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 10933 16 1 0 0 9 10935 624 71316 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10934 6733 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 c_p
S 10934 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1083138048 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 10935 16 1 0 0 9 10937 624 71320 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10936 6735 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 c_v
S 10936 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1082550272 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 10937 16 1 0 0 9 10939 624 71324 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10938 6737 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 r_d
S 10938 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1081208832 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 10939 16 0 0 0 9 10941 624 71328 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10940 6739 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 cp_cv
S 10940 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1073112970 -473225127 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 10941 16 0 0 0 9 10943 624 71334 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10942 6741 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rd_cp
S 10942 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1070746489 1591362385 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 10943 16 0 0 0 9 10945 624 71340 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10944 6743 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 cp_rd
S 10944 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1074527342 -1691049841 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 10945 16 0 0 0 9 10947 624 71346 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10946 6745 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 c_0
S 10946 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1077645340 -397417591 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 10947 16 0 0 0 9 10949 624 71350 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10948 6747 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 t_0
S 10948 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1081262080 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 10949 16 0 0 0 9 10951 624 71354 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10950 6749 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 n_0
S 10950 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1065646817 1202590843 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 10951 16 0 0 0 9 10953 624 71358 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10952 6751 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 p_0
S 10952 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1090021888 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 10953 16 0 0 0 9 10955 624 71362 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10954 6753 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 k_nu
S 10954 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1079164928 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 10955 16 0 0 0 9 10957 624 71367 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10956 6755 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 omg_earth
S 10956 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1058217364 261551826 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
S 10957 16 0 0 0 9 0 624 71377 800004 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10958 6757 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 r_earth
S 10958 3 0 0 0 9 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 1096306151 1073741824 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 9
A 67 1 0 0 0 56 681 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 70 1 0 0 0 65 683 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 86 1 0 0 0 95 723 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6731 2 0 0 6667 9 10932 0 0 0 6731 0 0 0 0 0 0 0 0 0 0
A 6733 2 0 0 6443 9 10934 0 0 0 6733 0 0 0 0 0 0 0 0 0 0
A 6735 2 0 0 5604 9 10936 0 0 0 6735 0 0 0 0 0 0 0 0 0 0
A 6737 2 0 0 5606 9 10938 0 0 0 6737 0 0 0 0 0 0 0 0 0 0
A 6739 2 0 0 6540 9 10940 0 0 0 6739 0 0 0 0 0 0 0 0 0 0
A 6741 2 0 0 5930 9 10942 0 0 0 6741 0 0 0 0 0 0 0 0 0 0
A 6743 2 0 0 4986 9 10944 0 0 0 6743 0 0 0 0 0 0 0 0 0 0
A 6745 2 0 0 6664 9 10946 0 0 0 6745 0 0 0 0 0 0 0 0 0 0
A 6747 2 0 0 6121 9 10948 0 0 0 6747 0 0 0 0 0 0 0 0 0 0
A 6749 2 0 0 6447 9 10950 0 0 0 6749 0 0 0 0 0 0 0 0 0 0
A 6751 2 0 0 5682 9 10952 0 0 0 6751 0 0 0 0 0 0 0 0 0 0
A 6753 2 0 0 6218 9 10954 0 0 0 6753 0 0 0 0 0 0 0 0 0 0
A 6755 2 0 0 5619 9 10956 0 0 0 6755 0 0 0 0 0 0 0 0 0 0
A 6757 2 0 0 6558 9 10958 0 0 0 6757 0 0 0 0 0 0 0 0 0 0
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
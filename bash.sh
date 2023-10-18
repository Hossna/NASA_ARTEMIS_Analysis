#!/bin/bash
# this bash code gives the original state from ARTEMIS (download directly from cdaweb ) and by deleting nonuseful info gives us a clean file "THB_L1_STATE_27227_finalEdited.txt1" which first coloumn is UT (real universal time and then x,y,z of ARTEMIS).
sed '/Year/d'  THB_L1_STATE_3519.txt1 > test1_temp.txt      # from "THB_L1_STATE_3519.txt1" file, delete line starting(containing # ) and save the result in test1_tmp.txt
sed '/UT/d' test1_temp.txt > test2_temp.txt              # from "test1_tmp.txt" file, delete line starting(containing UT ) and save the result in test2_tmp.txt
sed '/#/d' test2_temp.txt > test3_temp.txt              # from "test2_tmp.txt" file, delete line starting(containing dd ) and save the result in test3_tmp.txt
sed '/ (@/d' test3_temp.txt > test4_temp.txt             # from "test3_tmp.txt" file, delete line starting(containing @ ) and save the result in test3_tmp.txt
tr -s ' ' <test4_temp.txt | tr ' ' ',' > test5_temp.txt  #find all spaces from "test4_temp.txt" file and replace by comma(,) and save the result in test5_tmp.txt
awk -F "\"*,\"*" '{print $2,$4,$5,$6}' test5_temp.txt >> THB_L1_STATE_27227_finalEdited.txt1 # in file "test5_tmp.txt" find comma between colomn and copy coloumn $2,4,5 and 6 and put in a new file: THB_L1_STATE_27227_finalEdited.txt1. (But you substract  "0.001168" from coloumn $2 and then devide the results by 1000000 and then put in the file. )
rm *_temp.txt                                           # remove all files contain (_temp) in their name

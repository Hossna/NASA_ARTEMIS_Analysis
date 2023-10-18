set term postscript eps enhanced  "Helvetica" 11
set palette rgbformulae 22,13,10
set size 0.5,0.5

#set fit quiet

#********************************************************************************************************************************
#******************    cap  **************************************************************************************************
#********************************************************************************************************************************
set key at 20,3200

set xlabel "hhmm 2013 March 13" font "Times,12"
set yl "Km" font "Times,12"

set yr [-1000.0:3800.0]
set xr [5.0:50.0]

set xtics 5.0,10,55.0
set xtics ("8:40" 5.0, "8:50" 15.0, "9:00" 25.0, "9:12" 37.0, "9:25" 50.0)


set out 'cap_30.eps'
plot "ARTEMIS_normalized_rotated_positions.txt" using 1:5 title "R_c" with points pt 6 ps 0.5,"ARTEMIS_normalized_rotated_positions.txt" using 1:($6*-1.0) title "X_{sc}" with points pt 2 ps 0.5,\
"ARTEMIS_normalized_rotated_positions.txt" using 1:7 title "Y_{sc}" with points pt 8 ps 0.5
#********************************************************************************************************************************
#******************    cap **************************************************************************************************
#********************************************************************************************************************************



set yl "Density (cm^-3)" font "Times-Bold,12"
set xlabel "Universal Time"  font "Times-Bold,12"



##### density_measurment.out is the file that contains UT and density of Artemis I creat by copy and paste data from website of artemis
##### char_2.out is the density of Hutchinson model I wrote a code (calculating CapRadius) to get it



#set out "comparing_densityWihArtemis_full.eps"
#plot "Results_in_ART_Pos1.out" using 1:($11*3.04) title "Characteristic"  with points pt 7 ps 0.4 ,\
#"Results_in_ART_Pos0.out" using 1:($11*3.04) title "Analytic"  with points pt 6 ps 0.9, \
#"THB_Density_full_16851.out" using 1:2 title "Artemis"  with points pt 5 ps 0.5


#********************************************************************************************************************************
#******************    Curve  **************************************************************************************************
#********************************************************************************************************************************
unset key
set size ratio 0.5
set yr [-2.2:2.2]
set xr [-0.8:8.0]
set xlabel "X(R_{c})" font "Times,12"
set yl "Y(R_{c})" font "Times,12"

     
set out 'curve_art_30.eps'


set label "(dy/dx)^{+}" at 6.8,1.8 font "Times,10"
set label "(dy/dx)^{-}" at 6.8,-1.8 font "Times,10"
set label "8:42" at 6.0,0.0 font "Times,11"
#set label "9:00" at 6.1,0.0 font "Times,11"
set label "9:15" at 1.8,0.9 font "Times,11"
set label "9:00" at 0.9,0.45 font "Times,11"
set label "9:20" at 2.9,1.4 font "Times,11"
set label "9:23" at 6.0,1.6 font "Times,11"
set object  circle at 0.0,0.0 size first 1.0 fc rgb "white" front 

plot "Results_in_ART_Pos.out" using 3:4 with points pt 6 ps 0.5,'char_solarWind.out' u 1:2  w l lt 1 lw 2 ,\
'char_solarWind.out' u 1:3  w l lt 1 lw 2 ,\
#********************************************************************************************************************************
#********************************************************************************************************************************
#********************************************************************************************************************************









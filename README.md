# NASA_ARTEMIS_Analysis
1-I have downloaded "THB_L2_ESA_5259_allQuantities_Full.txt1" and THB_L2_ESA_5672_allQuantities_Reduced.txt1 from "cdfweb".
Then with bash.sh I created : "THB_Density_full.out" or "THB_Density_reduced.out" which has two coloumns (UT,Density) and will be used to 
compare with my Characteristic results (in plot_measurement.out). (my Bash.sh will make downloaded files from "cdfweb" of ARTEMIS, clean.)

2-I have downloaded "THB_L1_STATE_3519.txt1" from "cdfweb" website.It has position of ARTEMIS (x,y,z) in terms of date and UT time.
"THB_L1_STATE_27227_finalEdited.txt1" this file is the clean file which is: UT,x,y,z that can be used in my code for calculations.
3- I read UT,x,y,z from the file "THB_L1_STATE_27227_finalEdited.txt1" in my code and do calculations(cap radious stuff and rotations and then interpolations)
4- I calculate and save the density and velocities of my characteristic code in the positions of ARTEMIS in the file "Results_in_ART_Pos.out".
5- at the end I will compare "THB_Density_full.out" using 1:2" and ""Results_in_ART_Pos.out" using 1:($11*7.0)". coloumn $1 in both is UT (universal time) and the other coloumn is the total density.


charPDE_new_Analytic.f90 :is the file which made comparison between Analytic and ARTEMIS.
charPDE_new.f90          :is the file which made comparison between Characteristic and ARTEMIS.


In file "version1" for 1.1<x<1.9, I used Analytic model because my char model does not work in these regions 

#pos_sse.dat            : this file contains data of state of the ARTEMIS from the website (without any
                          changes) 
#cap_radius.f90         : this code read positions of ARTEMIS from pos_sse.dat, then calculate the cap
                          radius and normalized distance. put them in a file "normalized.out"
                          
                          
#ART_Hutch_90           : read positions from "pos_sse.dat" calculate the density from analytic
                          model(my modified Hutchinson model) on these positions from pos_sse.dat and write them in the file:char_2.out . 
                          
                          


#density_measurment.out : is the file that contains UT and density of Artemis I creat by
                          copy and paste data from website of artemis 
                          

#plot_measurment.gnu    : is a script to plot density of ARTEMIS (density_measurment.out)vs time in a same frame of density of 
                          Analytic model (char_2.out) vs time to compare. creat eps files

#plot_measurment_old.gnu : the same as "plot_measurment.gnu" but creat jpg files.



bash.sh                 : is the bash (command line language)file which is get the positions from ARTEMIS and delet everything except data.
                          
make the bash.sh file executable:         chmod +x bash.sh
remove the bash.sh from executable mode:         chmod +x bash.sh
ll                :                  shows that if it is executable if it has "x" means executable



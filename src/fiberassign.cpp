#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <exception>
#include <sys/time.h>
#include <map>
#include <stdexcept>
#include "modules/htmTree.h"
#include "modules/kdTree.h"
#include "misc.h"
#include "feat.h"
#include "structs.h"
#include "collision.h"
#include "global.h"


int main (int argc, char ** argv) {
    // argv[1] is the features file
    //// Initializations ---------------------------------------------
    //
    //  M is collection of targets, sky fibers, and standard stars
    //  each target has a priority provided by the mtl file
    //  all targets with the same priority are collected into a class
    // check_args(argc);
    // t for global, time for local
    Time t, time;
    init_time(t);
    Feat F;
    MTL M;
    Plates P;
    //for tabulation of results
    int total_used_SS = 0;
    int total_used_SF = 0;   
    //accumulate tiles actually used
    std::vector <int> tiles_really_used;

    //for tex 
    FILE * ftex;
    ftex=fopen("fortex.tex","w");
    
    F.parseCommandLine(argc, argv);

    // Read input files for standards, skys and targets.
    // Try to read SS and SF before targets to avoid wasting time if these
    // smaller files can't be read.
    init_time_at(time, "# Read target, SS, SF files", t);
    MTL SStars = read_MTLfile(F.SStarsfile, F, F.StarMask, 0);
    MTL SkyF   = read_MTLfile(F.SkyFfile,   F, 0, 1);
    MTL Targ   = read_MTLfile(F.Targfile,   F, 0, 0);
    print_time(time, "# ...read targets  took :");

    // combine the three input files
    M = Targ;
    printf(" Target size %lu \n", M.size() );
    std::cout.flush();

    // need to be able to match immutable target id to position in list
    // check for duplicates on mtl only to allow duplication with SS
    std::map <long long, int> invert_target;
    std::map <long long, int>::iterator targetid_to_idx;
    std::pair <std::map <long long, int>::iterator, bool> ret;
    for (unsigned i = 0; i < M.size(); ++i) {
        ret = invert_target.insert(std::make_pair(M[i].id, i) );
        // check for duplicates (std::map.insert only created keys, fails on
        // duplicate keys)
        if (ret.second == false ) {
            std::ostringstream o;
            o << "Duplicate targetid " << M[i].id << " in MTL";
            throw std::logic_error(o.str().c_str() );
        }
    }
    M.insert(M.end(), SStars.begin(), SStars.end() );
    printf(" Standard Star size %lu \n", M.size() );
    M.insert(M.end(), SkyF.begin(), SkyF.end() );
    printf(" Sky Fiber size %lu \n", M.size() );
    init_time_at(time, "# map position in target list to immutable targetid",
                 t);
    init_time_at(time, "# assign priority classes", t);
    F.Ngal = M.size();
    assign_priority_class(M);
    std::vector <int> total_used_by_class(M.priority_list.size(), 0);
    std::vector <int> count_class(M.priority_list.size(), 0);
    int Lya_count=0;
    int LRG2_count=0;
    for (size_t i = 0; i < M.size(); ++i) {
        if (!M[i].SS && !M[i].SF) {
            count_class[M[i].priority_class] += 1;
	    //keep track of Lya QSOs using their NUMOBS_NEEDED>3
	    if(M[i].nobs_remain>3){Lya_count +=1;}
	    //keep track of LRG2 using their NUMOBS_NEEDED==2
	    if(M[i].nobs_remain==2){LRG2_count +=1;}
        }
    }
    //numbers needed for results in percentages
    int ELG_start=count_class[0];
    int LRG1_start=count_class[1]-LRG2_count;
    int QSOt_start=count_class[2]-Lya_count;
    int Lya_start=Lya_count;
    int LRG2_start=LRG2_count;


    for (size_t i = 0; i < M.priority_list.size(); ++i) {
        printf("  class  %lu  number  %d\n", i, count_class[i]);
    }
    printf("   Lya QSO number %d\n",Lya_count);
    printf("   LRG2  number %d \n",LRG2_count);
    print_time(time, "# ...priority list took :");
    init_time_at(time, "# Start positioners", t);
    std::cout<<std::flush;    
    // fiber positioners
    F.Npetal = 10;  // spectrometers run 0 to 9 unless pacman
    FP pp = read_fiber_positions(F);
    read_fiber_status(pp, F);

    // order the fibers by their fiber number (fib_num) not simply order in
    // list
    // need to fix spectrom (List) and fp
    // each fiber has two co-ordinates so divide by two
    F.Nfiber = pp.size();

    // fibers per petal = 500
    F.Nfbp = F.Nfiber / F.Npetal;
    print_time(time, "# ..posiioners  took :");
    std::cout<<std::flush;   
    init_time_at(time, "# Start plates", t);

    // P is original list of plates
    P = read_plate_centers(F);
    F.Nplate = P.size();
    printf("# Read %s plate centers from %s and %d fibers from %s\n",
           f(F.Nplate).c_str(), F.tileFile.c_str(), F.Nfiber,
           F.fibFile.c_str() );
    print_time(time, "# ..plates   took :");

    // Computes geometries of cb and fh: pieces of positioner - used to
    // determine possible collisions
    F.cb = create_cb();  // cb=central body
    F.fh = create_fh();  // fh=fiber holder

    //// Collect available galaxies <-> tilefibers --------------------
    // HTM Tree of galaxies
    const double MinTreeSize = 0.01;
    init_time_at(time, "# Start building HTM tree", t);
    htmTree <struct target> T(M, MinTreeSize);
    print_time(time, "# Doing kd-tree... took :");  // T.stats();
    init_time_at(time, "# collect galaxies at ", t);

    // For plates/fibers, collect available galaxies; done in parallel  P[plate
    // j].av_gal[k]=[g1,g2,..]
    collect_galaxies_for_all(M, T, P, pp, F);
    print_time(time, "# ... took :");  // T.stats();
    init_time_at(time, "# collect available tile-fibers at", t);

    // For each galaxy, computes available tilefibers  G[i].av_tfs =
    // [(j1,k1),(j2,k2),..]
    collect_available_tilefibers(M, P, F);
    //set up results 
    printf(" epochs %d \n",F.num_epoch);
    std::cout<<std::flush;
    double QSOt_list[F.num_epoch];
    double LRG1_list[F.num_epoch];
    double ELG_list[F.num_epoch];
    double LRG2_list[2][F.num_epoch];
    double Lya_list[5][F.num_epoch];
 
    //// Assignment |||||||||||||||||||||||||||||||||||||||||||||||||||
    printf(" Nplate %d  Ngal %d   Nfiber %d \n", F.Nplate, F.Ngal, F.Nfiber);
    Assignment A(M, F);
    //loop over epochs
    for(int epoch=0;epoch<F.num_epoch;epoch++){
    // Make a plan ----------------------------------------------------
    print_time(t, "# Start assignment at : ");
    simple_assign(M, P, pp, F, A, epoch);

    // check to see if there are tiles with no galaxies
    // need to keep mapping of old tile list to new tile list and inverse map
    A.inv_order = initList(F.Nplate, -1);
    //need to clear out A.suborder
    A.suborder.clear();
    int inv_count = 0;
    for (int j = F.epoch_list[epoch]; j < F.Nplate; ++j) {
        bool not_done = true;
	int last_in_list=F.Nplate;
	if(epoch<F.num_epoch-1){last_in_list=F.epoch_list[epoch+1]-1;}
        for (int k = 0; k < F.Nfiber && not_done; ++k) {

            if (A.TF[j][k] != -1) {
                // suborder[jused] is jused-th used plate
                A.suborder.push_back(j);
                not_done = false;
                // inv_order[j] is -1 unless used
                A.inv_order[j] = inv_count;
                inv_count++;
		if(j<last_in_list){
		  tiles_really_used.push_back(j);
		}

            }
        }
    }
    F.NUsedplate = A.suborder.size();
    printf(" Plates actually used %d \n", F.NUsedplate);

    // Smooth out distribution of free fibers, and increase the number of
    // assignments
    // probably should not hard wire the limits i<1, i<3 in redistribute and
    // improve
    // more iterations will improve performance slightly
    for (int i = 0; i < 1; i++) {
        redistribute_tf(M, P, pp, F, A, 0);
    }
    for (int i = 0; i < 1; i++) {
        improve(M, P, pp, F, A, 0);
        redistribute_tf(M, P, pp, F, A, 0);
    }
    init_time_at(time, "# assign SS and SF ", t);
   
        // try assigning SF and SS before real time assignment
    for (int jused =0 ; jused < F.NUsedplate; ++jused) {
      
        int j = A.suborder[jused];

        // Assign SS and SF for each tile
        assign_sf_ss(j, M, P, pp, F, A);
        assign_unused(j, M, P, pp, F, A);
    }
    
       // fill all unused fibers with sky fibers
    for (int jused = 0; jused < F.NUsedplate; ++jused) {
        int j = A.suborder[jused];
        fill_unused_with_sf(j, M, P, pp, F, A);
    }
     

    // Results -------------------------------------------------------*/
    clean_up(M, P, pp, F,A,epoch);

    
    for (int jused = 0; jused < F.NUsedplate; ++jused) {
    std::vector <int> used_by_class(M.priority_list.size(), 0);
    int used_SS = 0;
    int used_SF = 0;           
    int j = A.suborder[jused];
        for (int k = 0; k < F.Nfiber; ++k) {
            int g = A.TF[j][k];
            if (g != -1) {
                if (M[g].SS) {
                    total_used_SS++;
                    used_SS++;
                } else if (M[g].SF) {
                    used_SF++;
                    total_used_SF++;
                } 
		  else {
                    used_by_class[M[g].priority_class]++;
                    total_used_by_class[M[g].priority_class]++;
		  }
	    }
	}
    }
  
    init_time_at(time, "# count SS and SF ", t);
    printf(" Totals SS   %4d    SF   %4d", total_used_SS, total_used_SF);
    for (size_t pr = 0; pr < M.priority_list.size(); ++pr) {
        printf(" class %2lu   %5d \n", pr, total_used_by_class[pr]);
    }
    std::cout<<std::flush;
    printf("\n done with epoch %d\n",epoch);
    int ELG_obs=0;
    int LRG1_obs=0;
    std::vector <int> LRG2_obs(2,0);
    int QSOt_obs=0;
    std::vector <int> QSOLya_obs(5,0);    
    int lastj=F.Nplate;
    if(epoch!=F.num_epoch-1){
      lastj=F.epoch_list[epoch+1]-1;
    }
    printf("** tiles really used %d",tiles_really_used.size());
    for (int jp = 0; jp < tiles_really_used.size(); ++jp) {
        int j=tiles_really_used[jp];  

        for (int k = 0; k < F.Nfiber; ++k) {
            int g = A.TF[j][k];
            if (g != -1) {
		    if(M[g].t_priority==3000){
		      ELG_obs++;
		    }
		    else if(M[g].t_priority==3200){
			LRG1_obs++;
		    }
		    else if(M[g].t_priority==3300){
		      LRG2_obs[M[g].nobs_done-1]++;
		    }   
		    else if(M[g].t_priority==3400){
			QSOt_obs++;
		    }

		    else if(M[g].t_priority==3500){
		      QSOLya_obs[M[g].nobs_done-1]++;
		    }
	    }
	}
    }

    double ELGpct=float(ELG_obs)/ELG_start;
    double LRG1pct=float(LRG1_obs)/LRG1_start;
    double QSOtpct=float(QSOt_obs)/QSOt_start;
 
    std::vector <double> Lya(5,0);
    std::vector <double> LRG2(2,0);
    for(int i=0;i<5;i++){
      Lya[i]=float(QSOLya_obs[i])/(i+1)/Lya_start;
    }
    for(int i=0;i<2;i++){
      LRG2[i]=float(LRG2_obs[i])/(i+1)/LRG2_start;
    }
	  
    printf("XX ELGs  %6.3f   \n",ELGpct);
    printf("XX LRG1s %6.3f   \n",LRG1pct);
    printf("XX LRG2s  (1)  %6.3f    (2)   %6.3f \n",LRG2[0],LRG2[1]);
    printf("XX QSO-t %6.3f   \n",QSOtpct);
    printf("XX QSO-Lya  (1)  %6.3f  (2)  %6.3f  (3) %6.3f  (4)  %6.3f  (5)  %6.3f\n",Lya[0],Lya[1],Lya[2],Lya[3],Lya[4]);
    printf("QSOt  %d   \n",QSOt_obs);
    
    ELG_list[epoch]=ELGpct;
    LRG1_list[epoch]=LRG1pct;
    QSOt_list[epoch]=QSOtpct;
    for (int i=0;i<2;i++){LRG2_list[i][epoch]=LRG2[i];}
    for (int i=0;i<5;i++){Lya_list[i][epoch]=Lya[i];}   
   
    }//end of epoch loop

    //print summary table
    printf(" type    epochs\n");

    
    fprintf(ftex,"\\begin{tabular}{l rrrrr}\\\\ \n \\hline\n");
    fprintf(ftex,"Type  ");
    for (int ep=0;ep<F.num_epoch;ep++){
      fprintf(ftex,"& epoch  %d ",ep);
    }
    
    fprintf(ftex,"\\\\ \n"); 
    fprintf(ftex,"  \\hline \n");
    fprintf(ftex," ELG  ");
    for (int ep=0;ep<F.num_epoch;ep++){
         fprintf(ftex,"& %6.3f  ",ELG_list[ep]);
    }  
    fprintf(ftex,"\\\\ \n"); 
    fprintf(ftex," LRG1  ");

    for (int ep=0;ep<F.num_epoch;ep++){
         fprintf(ftex,"& %6.3f  ",LRG1_list[ep]);
    }  
    fprintf(ftex,"\\\\ \n");
    fprintf(ftex,"LRG2 1  ");
    for (int ep=0;ep<F.num_epoch;ep++){
      fprintf(ftex," & %6.3f  ",LRG2_list[0][ep]);
    }  
    fprintf(ftex,"\\\\ \n");
  
    fprintf(ftex,"LRG2 2 ");
    for (int ep=0;ep<F.num_epoch;ep++){
      fprintf(ftex,"&  %6.3f ",LRG2_list[1][ep]);
    }  
    fprintf(ftex,"\\\\ \n");

    fprintf(ftex," QSO-t  ");
    for (int ep=0;ep<F.num_epoch;ep++){
         fprintf(ftex,"& %6.3f  ",QSOt_list[ep]);
    }  
    fprintf(ftex,"\\\\ \n"); 


    fprintf(ftex,"Lya 1 ");
    for (int ep=0;ep<F.num_epoch;ep++){
      fprintf(ftex,"& %6.3f   ",Lya_list[0][ep]);
    }  
    fprintf(ftex,"\\\\ \n");

    fprintf(ftex,"Lya 2 ");
    for (int ep=0;ep<F.num_epoch;ep++){
      fprintf(ftex,"& %6.3f  ",Lya_list[1][ep]);
    }  
    fprintf(ftex,"\\\\ \n");

    fprintf(ftex,"Lya  3");
    for (int ep=0;ep<F.num_epoch;ep++){
      fprintf(ftex,"& %6.3f  ",Lya_list[2][ep]);
    }  
    fprintf(ftex,"\\\\ \n");

    fprintf(ftex,"Lya 4  ");
    for (int ep=0;ep<F.num_epoch;ep++){
      fprintf(ftex," & %6.3f  ",Lya_list[3][ep]);
    }  
    fprintf(ftex,"\\\\ \n");

    fprintf(ftex,"Lya 5 ");
    for (int ep=0;ep<F.num_epoch;ep++){
      fprintf(ftex," & %6.3f  ",Lya_list[4][ep]);
    }  
    fprintf(ftex,"\\\\ \n");
    fprintf(ftex,"\\hline\n \\end{tabular}\n");

    init_time_at(time, "# print fits files ", t);
    /*
    for (int jused = 0; jused < F.NUsedplate; jused++) {
        int j = A.suborder[jused];
        fa_write(j, F.outDir, M, P, pp, F, A);  // Write output
    }
    */
    print_time(t, "# Finished !... in");
    return (0);
    
}

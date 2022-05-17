

#ifndef DEBYEWALLER_DIST
#define DEBYEWALLER_DIST

#include "analysis.h"
#include <string>

namespace std {

class DebyeWaller_Dist: public Analysis
{
    float * distribution;	//stores distribution in bins
    int n_bins;			//number of distribution bins
    float binsize;		//size of distribution bins
    int time_index;		//index of time-separation for which to calculate msd
    float time;			//actual time-separation corresponding to msd
    int atomcount;		//number of atoms used to produce this distribution (used for normalization)
    int timeweighting;		//normalization factor reflecting the number of timegaps used to produce distribution
    
    float mean;
    float variance;
    float square_term;
    
    int currenttime, nexttime;
   
   
  public:
    DebyeWaller_Dist(System* sys, int bins, float maxvalue, float t);
    
    Analysis_Type what_are_you(){Analysis_Type type = debyewaller_dist; return type;};		//virtual method to report the type of analysis
    
    void write(string)const;			//write distribution to file
    void write(ofstream&)const;
    
    void analyze(Trajectory_List* t_list);
    void list_displacementkernel(int timegapii, int thisii, int nextii);
    void listkernel(Trajectory * current_trajectory);
    void postprocess_list();
    
};

}
#endif
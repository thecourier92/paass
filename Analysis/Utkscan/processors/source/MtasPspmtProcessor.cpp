///@file MtasPspmtProcessor.cpp
///@brief Processes information from a series of SiPMs to determine position. Based on PspmtProcessor.cpp
///@author D.McKinnon, A. Keeler, S. Go, S. V. Paulauskas
///@date February 21, 2019



#include "DammPlotIds.hpp"
#include "DetectorDriver.hpp"
#include "Globals.hpp"
#include "Messenger.hpp"
#include "MtasPspmtProcessor.hpp"
#include "GeProcessor.hpp"
#include "HelperFunctions.hpp"

using namespace std;
using namespace dammIds::mtaspspmt;
vector<ChanEvent *> geEvents;

namespace dammIds {
namespace mtaspspmt {
const int DD_POS_DIAG = 1;
const int DD_POS_IMPL = 2;
// added by Cooper
const int D_RATE_IMPL = 3;
const int D_MULTIPL = 4;
const int DD_ENERGY = 5;
const int D_ENERGY = 6;
const int D_RAW_ENERGY = 7;
const int D_Beta_Gate_HPGe = 8;
const int D_Beta_Gamma_Time = 9;
const int DD_Beta_Gamma_Time2D = 10;
const int DD_Beta_Tdiff2D = 11;
const int DD_gamma_cycle = 12;
const int DD_Gated_Gamma_2D = 13;
const int DD_Left_Gated_2D = 14;
const int DD_Right_Gated_2D = 15;
const int DD_Both_Gated_2D = 16;
const int DD_Energy_A_Fire_A_Spec = 17;
const int DD_Energy_A_Fire_B_Spec = 18;
const int DD_Energy_B_Fire_B_Spec = 19;
const int DD_Energy_B_Fire_A_Spec = 20;
const int DD_Energy_AB_Fire_AB_Spec = 21;
const int D_Left_Gated_1D = 22;
const int D_Right_Gated_1D = 23;
const int D_Both_Gated_1D = 24;
const int D_DIAG_HITS = 25;
const int DD_DIAG_MULTIPLICITY_PER_CHANNEL = 26;
const int D_PSPMT_CHAN_1 = 27;
const int D_PSPMT_CHAN_2 = 28;
const int D_PSPMT_CHAN_3 = 29;
const int D_PSPMT_CHAN_4 = 30;
const int D_SUM_PSPMT = 31;

}  // namespace mtaspspmt
}  // namespace dammIds

void MtasPspmtProcessor::DeclarePlots() {
    DeclareHistogram2D(DD_POS_DIAG, SE, SE, "Diagnostic Detector Positions");
    DeclareHistogram2D(DD_POS_IMPL, SB, SB, "Implant Detector Positions");
    // Added by Cooper
    DeclareHistogram1D(D_RATE_IMPL, SE, "Implant Rate vs Current Time check");
    DeclareHistogram1D(D_MULTIPL, S4, "Left, Right, Total Hits");
    DeclareHistogram2D(DD_ENERGY, SD, S6, "Calibrated Ge Singles");
    DeclareHistogram1D(D_ENERGY,SD,"1D gamma spectrum");
    DeclareHistogram1D(D_RAW_ENERGY,SE,"1D Uncalibrated Spectrum");
    DeclareHistogram1D(D_Beta_Gate_HPGe,SD, "Beta gated gamma spectrum");
    DeclareHistogram1D(D_Beta_Gamma_Time,SC,"Time difference beta gamma");
    DeclareHistogram2D(DD_Beta_Gamma_Time2D,SC,SC,"Gamma energy versus time difference");
    DeclareHistogram2D(DD_Beta_Tdiff2D,SC,SC,"Beta energy versus time difference");
    DeclareHistogram2D(DD_gamma_cycle,SC,SC,"Ge energy time in cycle");
    DeclareHistogram2D(DD_Gated_Gamma_2D,SC,SD,"Gated Ge energy vs time 2D");
    DeclareHistogram1D(D_DIAG_HITS,S4,"2x2 Detector Hits Per Channel");
    DeclareHistogram2D(DD_DIAG_MULTIPLICITY_PER_CHANNEL,S4,S6,"2x2 Detector Multiplicity Per Channel");
    DeclareHistogram1D(D_PSPMT_CHAN_1,SE,"PSPMT 1D Energy Spectra, Channel 1");
    DeclareHistogram1D(D_PSPMT_CHAN_2,SE,"PSPMT 1D Energy Spectra, Channel 2");
    DeclareHistogram1D(D_PSPMT_CHAN_3,SE,"PSPMT 1D Energy Spectra, Channel 3");
    DeclareHistogram1D(D_PSPMT_CHAN_4,SE,"PSPMT 1D Energy Spectra, Channel 4");
    DeclareHistogram1D(D_SUM_PSPMT,SE,"All 4 Channels of PSPMT Summed");
    /*
    DeclareHistogram2D(DD_Left_Gated_2D,SC,SE,"Left SiPM Gated Ge energy vs time 2D");
    DeclareHistogram2D(DD_Right_Gated_2D,SC,SE,"Right SiPM Gated Ge energy vs time 2D");
    DeclareHistogram2D(DD_Both_Gated_2D,SC,SE,"Both SiPM Gated Ge energy vs time 2D");
    DeclareHistogram2D(DD_Energy_A_Fire_A_Spec,SC,SE,"Left SiPM fire, Left SiPM Spectrum 2D");
    DeclareHistogram2D(DD_Energy_A_Fire_B_Spec,SC,SE,"Left SiPM fire, Right SiPM Spectrum 2D");
    DeclareHistogram2D(DD_Energy_B_Fire_B_Spec,SC,SE,"Right SiPM fire, Right SiPM Spectrum 2D");
    DeclareHistogram2D(DD_Energy_B_Fire_A_Spec,SC,SE,"Right SiPM fire, Left SiPM Spectrum 2D");
    DeclareHistogram2D(DD_Energy_AB_Fire_AB_Spec,SC,SE,"Both SiPM fire, Both SiPM Spectrum 2D");
    */
    // DeclareHistogram1D(StartHist,SE,"Start times of beam")
}


MtasPspmtProcessor::MtasPspmtProcessor(const double &impl_scale, const unsigned int &impl_offset,
                                       const double &impl_threshold, const double &diag_scale,
                                       const unsigned int &diag_offset, const double &diag_threshold)
    : EventProcessor(OFFSET, RANGE, "MtasPspmtProcessor") {
    implPosScale_ = impl_scale;
    implPosOffset_ = impl_offset;
    implThreshold_ = 5;
    diagPosScale_ = diag_scale;
    diagPosOffset_ = diag_offset;
    diagThreshold_ = 30;

    associatedTypes.insert("mtaspspmt");
    // cout << "1";
    associatedTypes.insert("ge");
}

bool MtasPspmtProcessor::PreProcess(RawEvent &event) {
    if (!EventProcessor::PreProcess(event)) {
        return false;
    }
    // std::ofstream beta_t;
    // beta_t.open ("beta_times.csv",std::ios_base::app);
    // std::ofstream beta_e;
    // beta_e.open ("beta_energies.csv",std::ios_base::app);
    // std::ofstream gamma_t;
    // gamma_t.open ("gamma_times.csv",std::ios_base::app);
    // std::ofstream gamma_e;
    // gamma_e.open ("gamma_energies.csv",std::ios_base::app);

    if (evtNum_ == 0) {
        first_time = DetectorDriver::get()->GetFirstEventTimeinNs();
        first_time_beta = first_time;
    }
    map<string,int> group_map = 
    { 
        {"xa",0},
        {"xb",1},
        {"ya",2},
        {"yb",3},
    };
    double time1 = 0;
    double bunch_time = 0.0;
    double counter_true = 0.0;
    beta_ev = false;
    bool left_beta_ev = false;
    bool right_beta_ev = false;
    bool both_beta_ev = false; 
    bool gamma_ev = 0;
    bool dual_ev = 0;
    bool flag_on = 0;
    double sumImpl = 0;
    double beta_time,gamma_time, time_dif;
    evtNum_ = DetectorDriver::get()->GetEventNumber();
    // bool beta_ev, gamma_ev,has_beta;
    vector<ChanEvent *> impl = event.GetSummary("mtaspspmt:implant")->GetList();
    vector<ChanEvent *> diag = event.GetSummary("mtaspspmt:diag")->GetList();
    vector<ChanEvent *> geEvents = event.GetSummary("ge")->GetList();
    // cout << "check" << 2 << endl;
    // cout << "GE event size is " << (event.GetSummary("ge")->GetList()).size() << endl;

    // cout << 3 << endl;
    // static const vector<ChanEvent *> diag = event.GetSummary("mtaspspmt:diagnostic")->GetList();
    // if (!impl.empty()){beta_ev = true;}
    if (event.GetSummary("ge")->GetList().size() != 0){gamma_ev = true;}
    if (beta_ev && gamma_ev){dual_ev = true;}
    if (!impl.empty()){
    // cout << impl.size() << "5check" << endl;
    currentTime = impl.front()->GetTimeSansCfd() * 10; // removed 1.e9
    beta_time = currentTime;
    // cout << "check4" << endl;
    double delta_T = currentTime - first_time;
    double counter_bunch = 0.0;
    bunch_time = bunch_time + delta_T;
    if (delta_T > pow(10, 9)) {
        bin_num += 1;
        first_time = currentTime;
    }


    double energy = 0;
    double xa = 0, xb = 0, ya = 0, yb = 0;
    bool counter_xa = false, counter_xb = false;
    // cout << "checktime" << endl;

    // IMPLANT DETECTOR Calculations

    int flag1 = 0;
    for (auto it = impl.begin(); it != impl.end(); it++) {
        energy = (*it)->GetCalibratedEnergy();
        if (energy < implThreshold_) {
            continue;
        }
        beta_ev = true;
        double currentTime2 = impl.front()->GetTimeSansCfd() * 10;
        double lastTape = TreeCorrelator::get()->place("Beam_On")->last().time * Globals::get()->GetClockInSeconds();
//      cout <<"lastTape= "<< lastTape<<endl;
        double gTime = (*it)->GetTimeSansCfd() * Globals::get()->GetClockInSeconds();
        double  CycleTdiff =  gTime - lastTape ;


        std::string group = (*it)->GetChanID().GetGroup();
        if (group == "xa" && xa == 0) {
            left_beta_ev = true;
            xa = energy;
            sumImpl += energy;
            counter_xa = true;
            plot(D_MULTIPL,1);
            /* This is commented out solely to free up histogram space (scanor) for MTAS processor
            plot(DD_Energy_A_Fire_A_Spec,xa,CycleTdiff);
            */
        }
        //cout << "checkagain" << endl;
        if (group == "xb" && xb == 0) {
            right_beta_ev = true;
            xb = energy;
            sumImpl += energy;
            counter_xb = true;
            plot(D_MULTIPL,2);
            /* This is commented out solely to free up histogram space (scanor) for MTAS processor
            plot(DD_Energy_B_Fire_B_Spec,xb,CycleTdiff); */
        }

        if (xa > 0 && xb > 0) {
            double x = 0., y = 0.;
            both_beta_ev = true;
            x = (xa - xb) / (xa + xb);
            y = 0.5;  // shift it up off the axis
            plot(DD_POS_IMPL, x * implPosScale_ + implPosOffset_, y * implPosScale_ + implPosOffset_);
            plot(D_MULTIPL,3);
            plot(D_RATE_IMPL, bin_num);
            /* This is commented out solely to free up histogram space (scanor) for MTAS processor
            plot(DD_Energy_AB_Fire_AB_Spec,sumImpl,CycleTdiff);
            plot(DD_Energy_B_Fire_A_Spec,xa,CycleTdiff);
            plot(DD_Energy_A_Fire_B_Spec,xb,CycleTdiff);
            */

            // beta_t << std::to_string(currentTime) << "\n";
            // beta_e << std::to_string(xa+xb) << "\n";
        }

        if (DetectorDriver::get()->GetSysRootOutput()) {
            MTASPSstruct = processor_struct::MTASPSPMT_DEFAULT_STRUCT;
            MTASPSstruct.energy = energy;
            MTASPSstruct.time = currentTime2;
            // MTASPSstruct.subtype = "mtaspspmt:implant";
            // MTASPSstruct.group = group;
            MTASPSstruct.modNum = (*it)->GetModuleNumber();
            MTASPSstruct.chanNum = (*it)->GetChannelNumber();

        }
    }

    }

    // DIAGNOSTIC DETECTOR Calculations
    if (event.GetSummary("ge")->GetList().size() != 0){
    vector<ChanEvent *> geEvents = event.GetSummary("ge")->GetList();
    gamma_time = geEvents.front()->GetTimeSansCfd() * 10;
    double sumDiag = 0;
    // cout << geEvents.size();
    for (vector<ChanEvent *>::const_iterator ge = geEvents.begin(); ge != geEvents.end(); ge++) {
        // cout << "check";
        double energy = (*ge)->GetCalibratedEnergy();
        if (energy < geThreshold) {
            continue;
        }
        // cout << "yes";

        //if ((*ge)->IsSaturated() || (*ge)->IsPileup())
         //   continue;

        plot(DD_ENERGY, (*ge)->GetCalibratedEnergy(), (*ge)->GetChanID().GetLocation());
        if (TreeCorrelator::get()->place("Beam_On")->status()){
            double lastTape = TreeCorrelator::get()->place("Beam_On")->last().time * Globals::get()->GetClockInSeconds();
//            cout <<"lastTape= "<< lastTape<<endl;
            double gTime = (*ge)->GetTimeSansCfd() * Globals::get()->GetClockInSeconds();
            double  CycleTdiff =  gTime - lastTape ;
            plot(DD_gamma_cycle,CycleTdiff , energy); 
            plot(D_ENERGY,(*ge)->GetCalibratedEnergy());
            plot(DD_Gated_Gamma_2D,energy,CycleTdiff*1000);
            // This is commented out solely to free up histogram space (scanor) for MTAS processor
            if (beta_ev){
                plot(D_Beta_Gate_HPGe,(*ge)->GetCalibratedEnergy());
                /*time_dif = (gamma_time - beta_time) + 250; 
                plot(D_Beta_Gamma_Time,time_dif);
                plot(DD_Beta_Gamma_Time2D,time_dif,energy);
                plot(DD_Beta_Tdiff2D,time_dif,sumImpl); 

                if (left_beta_ev){
                plot(DD_Left_Gated_2D,energy,CycleTdiff);
                plot(D_Left_Gated_1D,energy);
                }
                if (right_beta_ev){
                plot(DD_Right_Gated_2D,energy,CycleTdiff);
                plot(D_Right_Gated_1D,energy);
                }
                if (both_beta_ev){
                plot(DD_Both_Gated_2D,energy,CycleTdiff);
                plot(D_Both_Gated_1D,energy);
                }
                */
            }
            
            plot(D_RAW_ENERGY,(*ge)->GetEnergy());
        }
        

        // gamma_t << std::to_string(currentTime) << "\n";
        // gamma_e << std::to_string(energy) << "\n";

    }
    }
    /*
    int flag_11 = 0;
    static const vector<ChanEvent *> &events = sumMap["logic"]->GetList();
    for (vector<ChanEvent *>::const_iterator it = events.begin(); it != events.end(); it++) {
        ChanEvent *chan = *it;
        string place = (*it)->GetChanID().GetPlaceName();
        string subtype = chan->GetChanID().GetSubtype();
        unsigned int loc = chan->GetChanID().GetLocation();

        if (flag_11 < 1){
        double t_0 = chan->GetTimeSansCfd()*10;
        flag_11 = 10;
        }
        double t_cur = chan->GetTimeSansCfd()*10 - t_0;

        if (subtype == "start") {
            if (!std::isnan(lastStartTime.at(loc))) {
                plot(StartHist,t_0);
            }

    }
    */



    // if(dual_ev){
    //     // cout << "has a dual event" << endl;
    //         vector<ChanEvent *> geEvents = event.GetSummary("ge")->GetList();
    // // cout << geEvents.size();
    //     for (vector<ChanEvent *>::const_iterator ge = geEvents.begin(); ge != geEvents.end(); ge++) {
    //         // cout << "check";
    //         double energy = (*ge)->GetCalibratedEnergy();
    //         if (energy < geThreshold) {
    //             continue;
    //         }
    //     // cout << "yes";

    //     //if ((*ge)->IsSaturated() || (*ge)->IsPileup())
    //      //   continue;

    //         plot(D_Beta_Gate_HPGe,(*ge)->GetCalibratedEnergy());
    // }
        
    // }
    
    //     if(beta_ev){
    //     // cout << "has a dual event" << endl;
    //     vector<vector<ChanEvent *>>::iterator row;
    //     vector<ChanEvent *>::iterator col;
    //     for (row = beta_gammas.begin(); row != beta_gammas.end(); row++) {
    //         for (col = row->begin(); col !=row.end();col++){
    //         // cout << "check";
    //             double energy = (*col)->GetCalibratedEnergy();
    //             if (energy < geThreshold) {
    //                 continue;
    //             }
   //             plot(D_Beta_Gate_HPGe,energy);
    //         }

    //     // cout << "yes";

    //     //if ((*ge)->IsSaturated() || (*ge)->IsPileup())
    //      //   continue;

            
    //     }      
    //     vector<vector<ChanEvent *>> beta_gammas;
    //     }
    // }

    
    // cout << "10";
    double evtnum0 = DetectorDriver::get()->GetEventNumber();
    
    // DIAGNOSTIC DETECTOR Calculations
    double energy = 0,xa = 0,xb = 0,ya = 0,yb = 0;
    double sumDiag = 0;
    double xa_hit = 0, xb_hit = 0, ya_hit = 0, yb_hit = 0;
    if (DetectorDriver::get()->GetSysRootOutput()){
        MTASPSstruct = processor_struct::MTASPSPMT_DEFAULT_STRUCT;
        for(auto it = diag.begin(); it != diag.end(); it++){
            energy = (*it)->GetCalibratedEnergy();
            MTASPSstruct.energy = energy;
            MTASPSstruct.time = diag.front()->GetTimeSansCfd() * 10;
            // MTASPSstruct.subtype = "mtaspspmt:diag";
            MTASPSstruct.group = (group_map.find((*it)->GetChanID().GetGroup())->second);
            MTASPSstruct.modNum = (*it)->GetModuleNumber();
            MTASPSstruct.chanNum = (*it)->GetChannelNumber();
            MTASPSstruct.evtnum1 = evtnum0;
            // cout << "the type of the group is " << typeid(group_map.find(MTASPSstruct.group)->second).name() << endl;
            pixie_tree_event_->mtaspspmt_vec_.emplace_back(MTASPSstruct);
            MTASPSstruct = processor_struct::MTASPSPMT_DEFAULT_STRUCT;
        }
    } 
    for (auto it = diag.begin(); it != diag.end(); it++) {


        // cout << "had a diag event" << endl;
        energy = (*it)->GetCalibratedEnergy();
        if ((energy < diagThreshold_) || (energy > 60000)){
            // cout << "below threshold" << endl;
            continue;
        }
        std::string group = (*it)->GetChanID().GetGroup();

        if (group == "xa"){
            plot(D_DIAG_HITS,1);
            xa_hit = xa_hit + 1;
        }
        
        if (group == "xb"){
            plot(D_DIAG_HITS,2);
            xb_hit = xb_hit + 1;
        }
        
        if (group == "ya"){
            ya_hit = ya_hit + 1;
            plot(D_DIAG_HITS,3);
        }
        
        if (group == "yb"){
            yb_hit = yb_hit + 1;
            plot(D_DIAG_HITS,4);
        }
        // cout << group << endl;
        if (group == "xa" && xa == 0) {
            xa = energy;
            sumDiag += energy;
            // cout << "inside xa" << endl;
        }
        if (group == "xb" && xb == 0) {
            xb = energy;
            sumDiag += energy;
        }
        if (group == "ya" && ya == 0) {
            ya = energy;
            sumDiag += energy;
        }
        if (group == "yb" && yb == 0) {
            yb = energy;
            sumDiag += energy;
        }
        // check how many times each group is being triggered per event





        if (xa > 0.0 && xb > 0.0 && ya > 0.0 && yb > 0.0) {
            double x = 0., y = 0.;
            double x1, x2,x3,y1,y2,y3, x4, y4, x_val,y_val;
            x = ((xa + ya) - (xb + yb)) / (xa + xb + ya + yb);
            y = ((xa + xb) - (ya + yb)) / (xa + xb + ya + yb);
            // xa = xa*10;
            // xb = xb*10;
            // ya = ya*10;
            // yb = yb*10;
            double tot1 = (xa + xb + ya + yb);
            // x = 1*((xa/(tot1))*1.0 + (xb/tot1)*1000.0 + (ya/tot1)*1.0 + (yb/tot1)*1000.0);
            // y = 1*((xa/(tot1))*1000.0 + (xb/tot1)*1000.0 + (ya/tot1)*1.0 + (yb/tot1)*1.0);
            // x1 = (xa/(xa+xb))*xa + (xb/(xa+xb))*xb;
            // y1 = (xa/(xa+ya))*xa + (ya/(xa+ya))*ya;

            // x2 = (xa/(xa+xb))*xa + (xb/(xa+xb))*xb;
            // y2 = (xb/(xb+yb))*xb + (yb/(xb+yb))*yb;

            // x3 = (ya/(ya+yb))*ya + (yb/(ya+yb))*yb;
            // y3 = (xa/(xa+ya))*xa + (ya/(xa+ya))*ya;

            // x4 = (ya/(ya+yb))*ya + (yb/(ya+yb))*yb;
            // y4 = (xb/(xb+ya))*xb + (ya/(xb+ya))*ya;            

            x1 = (xa/(xa+xb))*0.5;
            y1 = (xa/(xa+ya))*0.5;

            x2 = (xb/(xa+xb))*0.5;
            y2 = (xb/(xb+yb))*0.5;

            x3 = (ya/(ya+yb))*0.5;
            y3 = (ya/(xa+ya))*0.5;

            x4 = (yb/(ya+yb))*0.5;
            y4 = (yb/(xb+yb))*0.5;  

            double val1 = (x1+y1);
            double val2 = (x2+y2);  
            double val3 = (x3+y3);
            double val4 = (x4+y4);

            // x = (val1*1.0 + val2*2 + val3*1.0 + val4*2)*150;
            // y = (val1*2 + val2*2 + val3*1.0 + val4*1.0)*150;

            // x1 = (xa+ya)/2;
            // x2 = (xb+yb)/2;
            // x3 = (x1+x2)/2;
            // y1 = (xa+xb)/2;
            // y2 = (ya+yb)/2;
            // y3 = (y1 +y2)/2;

            // x1 = (xa + ya)/2;
            // x2 = (xa + xb)/2;
            // x3 = (xb + yb)/2;
            // x4 = (ya + yb)/2;
            
            // x_val = x3 - x1;
            // y_val = x2 - x4;
            // plot(DD_POS_DIAG, x3 * diagPosScale_ + diagPosOffset_, y3 * diagPosScale_ + diagPosOffset_);
            plot(DD_POS_DIAG, x*100 + 1000, y*100 + 1000);
	    //            if (group == "xa"){
                plot(D_PSPMT_CHAN_1,xa);
		//}
		//if (group == "xb"){
                plot(D_PSPMT_CHAN_2,xb);
		//}
		//if (group == "ya"){
                plot(D_PSPMT_CHAN_3,ya);
		//}
		//            if (group == "yb"){
                plot(D_PSPMT_CHAN_4,yb);
		//            }

            double temp_var1 = xa+xb+ya+yb;
            plot(D_SUM_PSPMT,temp_var1);
        }

    }

    plot(DD_DIAG_MULTIPLICITY_PER_CHANNEL,0,xa_hit);
    plot(DD_DIAG_MULTIPLICITY_PER_CHANNEL,1,xb_hit);
    plot(DD_DIAG_MULTIPLICITY_PER_CHANNEL,2,ya_hit);
    plot(DD_DIAG_MULTIPLICITY_PER_CHANNEL,3,yb_hit);
    
    // cout << "10";
    EndProcess();

    return (true);
}

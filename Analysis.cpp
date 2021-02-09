//g++ -Wall `root-config --glibs --cflags` Analysis.cpp -o Analysis && ./Analysis ../wavedump-3.10.0/wave0.dat

#include <fstream>
#include <iostream>
#include <vector>
#include <bitset>
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TPaveLabel.h"
#include "TROOT.h"
#include "TLine.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

using namespace std;

int main(int argc, char **argv){
  if (argc==1){
    cout << "give a filename"  << endl;
    return 0;
  }
  string filename = argv[1];
  cout << "Read file : " << filename << endl;
  TApplication app("ROOT Application", 0, 0);  

  //read file
  ifstream file(filename, ios::binary | ios::in);
  char buffer[4];
  int event_size, boardID, channel, counter, wf_length, wf_bin;
  unsigned long timer;
  double baseline, baselinerms;
  double tzero, tend, slope, slope_a, slope_b, pz;
  double energy, energy2;
  vector<double> x , y, y2, y3, y4;
  string outputname;
  long file_pos, file_size;
  
  TH1D* h_slope = new TH1D("h_slope","h_slope",1000,0,0.1);
  TH1D* h_energy = new TH1D("h_energy","h_energy",3000,0,3000);
  TH1D* h_energy2 = new TH1D("h_energy2","h_energy2",3000,0,3000);
  TH2D* h_ET = new TH2D("h_ET","h_ET",3000,0,3000,100,0,10);
  file.seekg (0, ios::end);
  file_size = file.tellg();
  cout << "file size: " << file_size << endl;
  file_pos = 0;
  file.seekg(0);
  
  
  if(file.is_open()){
    TCanvas*  mycanvas = new TCanvas("mycanvas","mycanvas",800,600);
    TCanvas*  c1 = new TCanvas("c1","c1",800,600);
    TCanvas*  c2 = new TCanvas("c2","c2",800,600);
    TCanvas*  c3 = new TCanvas("c3","c3",800,600);

    
    while (!file.eof()) {    
      file.read(buffer, 4);  //header0 4x8bit == 32bit
      event_size = 0;
      for (int i = 0;i<4;i++){
        event_size += ((buffer[i] & 0xFF)<< 8*i);
      }
      file_pos += event_size;
      char buffer_event[event_size-4];

      file.read(buffer_event,event_size-4);
      boardID = 0;
      channel = 0;
      counter = 0;
      timer = 0;
      
      for (int i = 0;i<4;i++){
        boardID += ((buffer_event[i] & 0xFF) << 8*i);
        channel += ((buffer_event[i+8] & 0xFF) << 8*i);
        counter += ((buffer_event[i+12] & 0xFF) << 8*i);        
        timer   += ((buffer_event[i+16] & 0xFF) << 8*i);
//        std::bitset<8> x(buffer_event[i+16]);
//        cout << i << " " << (buffer_event[i+16] & 0xFF) << " " << x << endl;
      }
     
      wf_length = (event_size - 6*4)/2;
      x.clear();
      y.clear();
      y2.clear();
      y3.clear();
      y4.clear();
      x.reserve(wf_length);
      y.reserve(wf_length);   //wf
      y2.reserve(wf_length);  //wf -baseline
      y3.reserve(wf_length);  //pz
      y4.reserve(wf_length);  //pz int
      
      baseline = 0;
      baselinerms = 0;
      tzero = 0;
      tend = 0;
      slope_a = 0;
      slope_b = 0;
      energy =0;
      energy2 =0;
      
      for (int i = 0;i<wf_length;i++){
        wf_bin = (buffer_event[2*i+20] & 0xFF);
        wf_bin += (buffer_event[2*i+21] & 0xFF) << 8;
        x[i] = i*4./1000;
        y[i] = wf_bin;
        if(i<200) baseline+= wf_bin/200.;
      }
      
      
      for (int i = 0;i<wf_length;i++){
        y2[i] = -1*(y[i] - baseline);
        y3[i] = y2[i]; //for pz;
        if(i<200){
          baselinerms += y2[i]*y2[i];
        }
        if(i==200) baselinerms = sqrt(baselinerms/200);
        if ( (i>200) && (tzero==0)) {  ///assume tzero > 200 otherwise extra loop
          if ( (y2[i] > 5*baselinerms) && (y2[i-1] > 5*baselinerms) ){
            tzero = i*4./1000; 
          }
        }
        if ( (i>wf_length-400) && (i<wf_length-200)) slope_a += y2[i]/200.;
        if ( (i>wf_length-200) && (i<wf_length)) slope_b += y2[i]/200.;
      }
     
      slope = log(slope_a/slope_b) / (200*4/1000.);
      h_slope->Fill(slope);
      if (slope<1E-3) continue;
      
      for (int i = 0;i<wf_length;i++){
        for(int j  = i+1;j<wf_length;j++){
          y3[j] -= y3[i] * exp(-0.0214171*(j-i)*4/1000.); 
        }
        if(i==0) y4[i]=y3[i];
        else y4[i]= y4[i-1] + y3[i];
        if( i > wf_length-200) energy += y4[i]/200;
        if ((tend==0) && (i> tzero*1000/4.)){
          slope_a = 0;
          slope_b = 0;
          for(int j = 0;j<20;j++){
            slope_a += y4[i-j-20];
            slope_b += y4[i-j];
          }
          //cout << i*4/1000. << " " << slope_a - slope_b << " " << baselinerms << endl;
          if (slope_b - slope_a < 3*baselinerms) tend = i*4/1000.;     
        }
        y3[i] = y3[i]*10;        
      }
      // David 1us ~ 0.11 keV
      //1460: 0.1s == 10keV : 1s == 78.77channel
      //2614: 0.2s == 20keV : 2s == 157.54channel
      //-2 ~ offset
      energy2 = energy + (tend-tzero-2)*78.77;       

      energy = energy*0.7877 -2.7836;
      energy2 = energy2*0.7877 -2.7836;
      
      h_energy->Fill(energy);
      h_energy2->Fill(energy2);
      h_ET->Fill(energy,tend-tzero);
      
      if(counter%50 == 0){
        cout << "-------" << endl;
        cout << "Size  : " << event_size << endl;
        cout << "where : " << 1.0*file_pos/file_size << endl;
        cout << "ID    : " << boardID << endl;
        cout << "Ch    : " << channel << endl;
        cout << "#     : " << counter << endl;
        cout << "time  : " << timer << endl;
        cout << "length: " << wf_length << endl;
        cout << "baseline : " << baseline << "+-" << baselinerms <<  endl;
        cout << "T0       : " << tzero << endl;
        cout << "TEND     : " << tend << endl;
        cout << "slope    : " << slope << endl;
        cout << "energy   : " << energy << endl;
        cout << "energy2   : " << energy2 << endl;
        
        c1->cd();
        h_slope->Draw();
        c1->Update();
        c2->cd();
        h_energy->SetLineColor(kBlack);
        h_energy->Draw();
        h_energy2->SetLineColor(kRed);
        h_energy2->Draw("SAME");
        c2->Update();
        c3->cd();    
        h_ET->Draw("COLZ");
        c3->Update();
        gSystem->ProcessEvents();  
      }
     
      if((tend==0)){
        TLine *line1 = new TLine(tzero,0,tzero,energy);
        line1->SetLineColor(kBlue);
        TLine *line2 = new TLine(tend,0,tend,energy);
        line2->SetLineColor(kBlue);     
      
        TGraph* gr = new TGraph(wf_length,&x[0],&y4[0]);
        gr->SetLineColor(kBlue);
        TGraph* gr2 = new TGraph(wf_length,&x[0],&y2[0]);
        gr->SetLineColor(kBlack);
        TGraph* gr3 = new TGraph(wf_length,&x[0],&y3[0]);
        gr3->SetLineColor(kRed);
      
        TMultiGraph *mg = new TMultiGraph();
        mg->Add(gr);
        mg->Add(gr2);
        mg->Add(gr3);      
      
        mycanvas->cd();
        mg->Draw("AL");
        line1->Draw("SAME");
        line2->Draw("SAME");
        mycanvas->Modified();
        mycanvas->Update();
        gSystem->ProcessEvents();  
        /*
        cin.get();
        mycanvas->Modified();
        mycanvas->Update();
        gSystem->ProcessEvents();  
        cin.get();
        */
        outputname = "wf_"+ to_string(counter) +".pdf";
        mycanvas->SaveAs(outputname.c_str());
        delete mg, gr, gr2, gr3 ,line1, line2;
      }
      
      //if(counter>1000) break;

    }
    file.close();
    
    c1->cd();
    h_slope->Draw();
    c1->Update();
    c2->cd();
    h_energy->SetLineColor(kBlack);
    h_energy->Draw();
    h_energy2->SetLineColor(kRed);
    h_energy2->Draw("SAME");
    c2->Update();
    c3->cd();    
    h_ET->Draw("COLZ");
    c3->Update();

    TFile* out = new TFile("output.root","Recreate");
    h_slope->Write();
    h_energy->Write();
    h_energy2->Write();
    h_ET->Write();
    out->Close();
  }

  app.Run();
  return 0;

}





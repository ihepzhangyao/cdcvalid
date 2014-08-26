int getMatch(TFile* f, char* cut){
 TTree* t = f->Get("evt");
 return t->Draw("nMatchMdc",cut,"goff");
}
drawEffDiff(){
  const int flag_eff_found =0;//0, found nTk>=1; 1,found nTk==1 as found
  const int N = 6;
  string pt_x[N] ={"50", "100","150","200","250","300"};
  char name1[N][200];
  for(int ii=0;ii<6;ii++){
    //sprintf(name1[ii],"/ihepbatch/besd13/zhangyao/fixpt/official/xpt_muon_%sMeV.root",pt_x[ii].c_str());//combine pat+tsf 663p01 official
    //sprintf(name1[ii],"/ihepbatch/besd13/zhangyao/fixpt/663p01/pat/fixpt_muon_%sMeV.root",pt_x[ii].c_str());//663p01 pat only
    sprintf(name1[ii],"/ihepbatch/besd13/zhangyao/fixpt/660/pat/fixpt_muonplus_%sMeV.root",pt_x[ii].c_str());//663p01 pat only
    std::cout<<name1[ii]<<"   "<<std::endl;
  }
  TFile* f1[N];
  f1[0] = new TFile(name1[0]);
  f1[1] = new TFile(name1[1]);
  f1[2] = new TFile(name1[2]);
  f1[3] = new TFile(name1[3]);
  f1[4] = new TFile(name1[4]);
  f1[5] = new TFile(name1[5]);

  double eff663p01[N];
  double efferr663p01[N];
  double multi663p01[N];
  double multierr663p01[N];
  //std::cout<<"663p01    "<<std::endl;
  for( int i=0; i<N; i++){
    float tot = getMatch(f1[i],""); 
    float found = getMatch(f1[i],"nMatchMdc>=1");
    float found0 = getMatch(f1[i],"nMatchMdc==0");
    float found1 = getMatch(f1[i],"nMatchMdc==1");
    float foundMulti = getMatch(f1[i],"nMatchMdc>1");
    float foundMulti1 = getMatch(f1[i],"nTrk>1");
    std::cout<<pt_x[i]<<" found "<<found<<" found0="<<found0<<" found1="<<found1
      <<" foundMulti="<<foundMulti<<" nTot "<<tot<<" missing="<<found0/tot*100.<<" eff="<<found1/tot*100.<<" multi="<<foundMulti/tot*100.<<" nTrk>1 "<<foundMulti1<<std::endl;
    if( flag_eff_found == 0 ) {
      eff663p01[i]=found/tot;
    }else if(flag_eff_found==1) {
      eff663p01[i]=found1/tot;
    }
    efferr663p01[i]=sqrt(eff663p01[i]*(1-eff663p01)/tot);
    multi663p01[i]=foundMulti/tot;
    multierr663p01[i]=sqrt(multi663p01[i]*(1-multi663p01)/tot);
  }

  TFile* f2[N];
  char name2[N][200];
  for(int ii=0;ii<6;ii++){
    //sprintf(name2[ii],"/ihepbatch/besd13/zhangyao/fixpt/official_noPtCut/xpt_muon_%sMeV.root",pt_x[ii].c_str());// 663p01 no pt cut 
    //sprintf(name2[ii],"/ihepbatch/besd13/zhangyao/fixpt/official_noPtCut_gapCut/xpt_muon_%sMeV.root",pt_x[ii].c_str());// 663p01 no pt cut with gap cut v1
    //sprintf(name2[ii],"/ihepbatch/besd13/zhangyao/fixpt/official_noPtCut_gapCutV2/xpt_muon_%sMeV.root",pt_x[ii].c_str());// 663p01 no pt cut with gap cut v2
    //sprintf(name2[ii],"/ihepbatch/besd13/zhangyao/fixpt/663p01/tsf/fixpt_muon_%sMeV.root",pt_x[ii].c_str());//663p01 tsf only
    sprintf(name2[ii],"/ihepbatch/besd13/zhangyao/fixpt/660/tsf/fixpt_muonplus_%sMeV.root",pt_x[ii].c_str());//663p01 pat only
    ///ihepbatch/besd13/zhangyao/fixpt/660/pat/fixpt_muonplus_50MeV.root
    std::cout<<name2[ii]<<"   "<<std::endl;
  }
  f2[0] = new TFile(name2[0]);
  f2[1] = new TFile(name2[1]);
  f2[2] = new TFile(name2[2]);
  f2[3] = new TFile(name2[3]);
  f2[4] = new TFile(name2[4]);
  f2[5] = new TFile(name2[5]);

  //std::cout<<"663p01 noPtCut   "<<std::endl;
  double effNoPtCut[N];
  double efferrNoPtCut[N];
  double multiNoPtCut[N];
  double multierrNoPtCut[N];
  for( int i=0; i<N; i++){
    float tot = getMatch(f2[i],""); 
    float found = getMatch(f2[i],"nMatchMdc>=1");
    float found0 = getMatch(f2[i],"nMatchMdc==0");
    float found1 = getMatch(f2[i],"nMatchMdc==1");
    float foundMulti = getMatch(f2[i],"nMatchMdc>1");
    float foundMulti1 = getMatch(f2[i],"nTrk>1");
    std::cout<<pt_x[i]<<" found "<<found<<" found0="<<found0<<" found1="<<found1
      <<" foundMulti="<<foundMulti<<" nTot "<<tot<<" missing="<<found0/tot*100.<<" eff="<<found1/tot*100.<<" multi="<<foundMulti/tot*100.<<" nTrk>1 "<<foundMulti1<<std::endl;
    if( flag_eff_found == 0 ) {
      effNoPtCut[i]=found/tot;
    }else if(flag_eff_found==1) {
      effNoPtCut[i]=found1/tot;
    }
    efferrNoPtCut[i]=sqrt(effNoPtCut[i]*(1-effNoPtCut)/tot);
    multiNoPtCut[i]=foundMulti/tot;
    multierrNoPtCut[i]=sqrt(multiNoPtCut[i]*(1-multiNoPtCut)/tot);
  }

  double x[N] ={50,100,150,200,250,300};
  double xerr[N] ={0,0,0,0,0,0};
  TCanvas *c1 = new TCanvas("c1","single track efficiency",600,900);  //ww,wh
  c1->Divide(1,2);
  c1->cd(1);

  //===================================================================================
  TGraphErrors* gr663p01= new TGraphErrors(N, x, eff663p01, xerr, efferr663p01);
  TGraphErrors* grNoPtCut= new TGraphErrors(N, x, effNoPtCut, xerr, efferrNoPtCut);

  gr663p01->SetMarkerStyle(kFullCircle);
  gr663p01->SetMarkerColor(kBlue);
  gr663p01->SetMarkerSize(1.5);
  gr663p01->SetLineColor(kBlue);

  grNoPtCut->SetMarkerStyle(kFullCircle);
  grNoPtCut->SetMarkerColor(kRed);
  grNoPtCut->SetMarkerSize(1.5);
  grNoPtCut->SetLineColor(kRed);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr663p01);
  mg->Add(grNoPtCut);
  mg->Draw("APL");
  mg->SetTitle("sigle track efficiency 660;p_{T}(MeV/c);Tracking Efficiency");

  TLegend *leg = new TLegend(0.66,0.63,0.9,0.82);  // x1,y1,x2,y2
  //leg->SetTextSize(0.04);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  //leg->AddEntry(gr663p01,"663p01","p");
  //leg->AddEntry(grNoPtCut,"noPtCut","p");
  leg->AddEntry(gr663p01,"pat","p");
  leg->AddEntry(grNoPtCut,"tsf","p");
  //leg->SetTextFont(32);  // font= times bold italic
  leg->Draw();

  //===================================================================================
  c1->cd(2);
  TGraphErrors* grMulti663p01= new TGraphErrors(N, x, multi663p01, xerr, multierr663p01);
  TGraphErrors* grMultiNoPtCut= new TGraphErrors(N, x, multiNoPtCut, xerr, multierrNoPtCut);

  grMulti663p01->SetMarkerStyle(kFullCircle);
  grMulti663p01->SetMarkerColor(kBlue);
  grMulti663p01->SetMarkerSize(1.5);
  grMulti663p01->SetLineColor(kBlue);

  grMultiNoPtCut->SetMarkerStyle(kFullCircle);
  grMultiNoPtCut->SetMarkerColor(kRed);
  grMultiNoPtCut->SetMarkerSize(1.5);
  grMultiNoPtCut->SetLineColor(kRed);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(grMulti663p01);
  mg->Add(grMultiNoPtCut);
  mg->Draw("APL");
  mg->SetTitle("Multi-track rate 660;p_{T}(MeV/c);Multi-track rate");

  TLegend *leg = new TLegend(0.66,0.63,0.9,0.82);  // x1,y1,x2,y2
  //leg->SetTextSize(0.04);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  //leg->AddEntry(grMulti663p01,"663p01","p");
  //leg->AddEntry(grMultiNoPtCut,"noPtCut","p");
  leg->AddEntry(grMulti663p01,"pat","p");
  leg->AddEntry(grMultiNoPtCut,"tsf","p");
  //leg->SetTextFont(32);  // font= times bold italic
  leg->Draw();

  /*
     c1->cd(2);
     double diffEff[N];
     double diffEffErr[N];
     for(int i=0;i<N; i++){
     diffEff[i]=patEff[i]-tsfEff[i];
     diffEffErr[i]=sqrt(patEffErr[i]*patEffErr[i]+tsfEffErr[i]*tsfEffErr[i]);
     }
     TGraphErrors* grdiff= new TGraphErrors (N, x, diffEff, xerr, diffEffErr);
     grdiff->Draw("AP");
   */
}

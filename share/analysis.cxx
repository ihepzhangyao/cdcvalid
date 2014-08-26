#define ALGNAME "all"
#define MAG "1.0"
#define MAGINDEX 2
//		     0.1   0.5   0.9   1.0  m1.0 m1.0n   2.0
float PRESCUT[7] = {0.25, 0.10, 0.05, 0.05, 0.05, 0.05, 0.05};
#define VRCUT 0.1
#define VZCUT 5

void analysis()
{
	char treename[100], histname[100];
	sprintf(treename, "/panfs/panfs.ihep.ac.cn/home/data/huangl/fixpt-%s/%s_tree.root", MAG, ALGNAME);
	sprintf(histname, "/panfs/panfs.ihep.ac.cn/home/data/huangl/fixpt-%s/%s_hist.root", MAG, ALGNAME);

	TFile *treefile = new TFile(treename);
	TTree *rec = treefile->Get("rec");
	TTree *truth = treefile->Get("truth");

	TFile *histfile = new TFile(histname);
	TFolder *folderTrk = (TFolder *)histfile->FindKey("track")->ReadObj();

	//Results
	TH1F *h_resi = (TH1F *)folderTrk->FindObject("resi_all");
	int nRecTrk = 0, nTruthTrk;
	TH1F *h_nHit = new TH1F("h_nHit", "nHit", 110, 0, 110);
	TH1F *h_pres = new TH1F("h_pres", "P RES", 200, -(PRESCUT[MAGINDEX] * 1.25), PRESCUT[MAGINDEX] * 1.25);
	TH1F *h_vr = new TH1F("h_vr", "Vr", 100, -(VRCUT * 1.25), VRCUT * 1.25);
	TH1F *h_vz = new TH1F("h_vz", "Vz", 100, -(VZCUT * 1.25), VZCUT * 1.25);
	TH1F *h_cn = new TH1F("h_cn", "chi2ndof", 100, 0, 3);

	int iTrk, nTrk, nFailed = 0;
	int trkId;
	double rateTruth, bestrate = -1;
	int nHit, bestnhit;
	double p, pTruth, bestpres;
	double vr, bestvr;
	double vz, bestvz;
	double chi2, bestcn;
	int ndof;
	int filter, bestfilter;

	rec->SetBranchAddress("rateTruth", &rateTruth);
	rec->SetBranchAddress("trkId", &trkId);
	rec->SetBranchAddress("nHit", &nHit);
	rec->SetBranchAddress("p", &p);
	rec->SetBranchAddress("vr", &vr);
	rec->SetBranchAddress("vz", &vz);
	rec->SetBranchAddress("chi2", &chi2);
	rec->SetBranchAddress("ndof", &ndof);
	rec->SetBranchAddress("filterStat_pi", &filter);
	rec->SetBranchAddress("pTruth", &pTruth);

	nTrk = rec->GetEntries();
	nTruthTrk = truth->GetEntries();
	for (iTrk = 0; iTrk < nTrk; iTrk++) {
		rec->GetEntry(iTrk);
		if (rateTruth < -4) continue;
		if (trkId == 0) {
			if (bestrate > 0) {
				nRecTrk++;
				h_nHit->Fill(bestnhit);
				if (fabs(bestpres) < PRESCUT && bestfilter == 0) h_pres->Fill(bestpres);
				if (fabs(bestvr) < VRCUT) h_vr->Fill(bestvr);
				if (fabs(bestvz) < VZCUT) h_vz->Fill(bestvz);
				if (bestfilter == 1) nFailed++;
				h_cn->Fill(bestcn);
			}	//Do something
			bestrate = rateTruth;
			bestnhit = nHit;
			bestpres = p - pTruth;
			bestvr = vr;
			bestvz = vz;
			bestfilter = filter;
			bestcn = chi2 / ndof;
		}
		else if (rateTruth > bestrate) {
			bestrate = rateTruth;
			bestnhit = nHit;
			bestpres = p - pTruth;
			bestvr = vr;
			bestvz = vz;
			bestfilter = filter;
			bestcn = chi2 / ndof;
		}
	}
	if (bestrate > 0) {
		nRecTrk++;
		h_nHit->Fill(bestnhit);
		if (fabs(bestpres) < PRESCUT && bestfilter == 0) h_pres->Fill(bestpres);
		if (fabs(bestvr) < VRCUT) h_vr->Fill(bestvr);
		if (fabs(bestvz) < VZCUT) h_vz->Fill(bestvz);
		if (bestfilter == 1) nFailed++;
		h_cn->Fill(bestcn);
	}	//Do something

	char epsname[100];
	TCanvas *canvas = new TCanvas("canvas", "");

	sprintf(epsname, "%s_%s_pres.eps", MAG, ALGNAME);
	h_pres->Draw();
	TF1 *ff = new TF1("ff", "gaus", -0.02, 0.02);
	h_pres->Fit("ff", "R");
	canvas->Update();
	canvas->Print(epsname);

	sprintf(epsname, "%s_%s_nhit.eps", MAG, ALGNAME);
	h_nHit->Draw();
	canvas->Update();
	canvas->Print(epsname);

	sprintf(epsname, "%s_%s_cn.eps", MAG, ALGNAME);
	h_cn->Draw();
	canvas->Update();
	canvas->Print(epsname);

	sprintf(epsname, "%s_%s_resi.eps", MAG, ALGNAME);
	h_resi->Draw();
	canvas->Update();
	canvas->Print(epsname);

	sprintf(epsname, "%s_%s_vr.eps", MAG, ALGNAME);
	h_vr->Draw();
	canvas->Update();
	canvas->Print(epsname);

	sprintf(epsname, "%s_%s_vz.eps", MAG, ALGNAME);
	h_vz->Draw();
	canvas->Update();
	canvas->Print(epsname);

	cout << "P Res " << h_pres->GetRMS() << endl;
	cout << "Resi " << h_resi->GetRMS() << endl;
	cout << "Vr " << h_vr->GetRMS() << endl;
	cout << "Vz " << h_vz->GetRMS() << endl;
	cout << "hHit Mean " << h_nHit->GetMean() << endl;
	cout << "Failed number " << nFailed << endl;
	cout << "All NTrack " << nTruthTrk << "\tRec NTrack " << nRecTrk << endl;
}

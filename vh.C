#include "SpacePoint.h"

void Helper1DHistDraw(TCanvas* c, TH1D* h, std::string pdf){
    c->Clear();
    h->Draw("E1p");
    c->Print(pdf.c_str());
}
void Helper2DHistDraw(TCanvas* c, TH2D* h, std::string pdf){
    c->Clear();
    h->Draw("colz");
    c->Print(pdf.c_str());
}

void heavyneutrino(const char* pathname) {
    //pathname is the directory where standard root intput is saved, the output HEPevt file will be save under the same directory	
      
	std::ofstream fout;
	int file_size = 50000; //each file contain 50k events info
 	std::string path(pathname);


	//set up detector geometry
	SpacePoint sp; 
	sp.SetX(-1.55, 254.8);  //x, y, z boundary of the active TPC, in cm
	sp.SetY(-115.53, 117.47);
	sp.SetZ(0.1, 1036.8); 
	sp.SetTime(3125, 4725); //time of interaction, in ns


        //open the standard root input
	TFile* tfin = new TFile((path+"/vh.root").c_str(), "read");
	TTree* tin = (TTree*)tfin->Get("t");

	std::string coherent_photon_cut = "nfem==1 && pdgf[0]==22";

	int Nparticle = 4;	
	double Ef[Nparticle];
	double pxf[Nparticle];
	double pyf[Nparticle];
	double pzf[Nparticle];
	int pdgf[Nparticle];
	int nfem;
  	int nf;
	int neu;
	double Q2;
	double Ev;
	
	tin->SetBranchAddress("neu", &neu);
	tin->SetBranchAddress("Q2", &Q2);
	tin->SetBranchAddress("Ev", &Ev);
	tin->SetBranchAddress("nfem", &nfem);
	tin->SetBranchAddress("nf", &nf);
	tin->SetBranchAddress("pdgf", &pdgf[0]);
	tin->SetBranchAddress("Ef", &Ef[0]);
	tin->SetBranchAddress("pxf", &pxf[0]);
	tin->SetBranchAddress("pyf", &pyf[0]);
	tin->SetBranchAddress("pzf", &pzf[0]);

	//helper histograms
	TH2D* h_xy_dist = new TH2D("xydistribution", "X vs Y coordinates; x coordinate (cm); y coordinate (cm)", 50,-1.55, 254.8, 50, -115.53, 117.47);
	TH2D* h_xz_dist = new TH2D("xzdistribution", "X vs Z coordinates; x coordinate (cm); z coordinate (cm)", 50, -1.55, 254.8, 50, 0.1, 1036.8);
	TH2D* h_yz_dist = new TH2D("yzdistribution", "Y vs Z coordinates; y coordinate (cm); z coordinate (cm)", 50, -115.3, 117.47, 50, 0.1, 1036.8); 
	TH1D* h_x_coords = new TH1D("xdistribution", "Distribution of x Coordinates; x coordinate (cm); Entries [area normalized]", 50, -1.55, 254.8);
	TH1D* h_y_coords = new TH1D("ydistribution", "Distribution of y Coordinates; y coordinate (cm); Entries [area normalized]", 50, -115.3, 117.47);
	TH1D* h_z_coords = new TH1D("zdistribution", "Distribution of z Coordinates; z coordinate (cm); Entries [area normalized]", 50, 0.1, 1036.8);	

	for (int i = 0; i < tin->GetEntries(); ++i) {
	//for (int i = 0; i < 5; ++i) {

		// handle txt file writeout
		if(i% file_size == 0){
		    std::cout << "On Entry: " << i<<"/" << tin->GetEntries() << std::endl;

		    if(fout.is_open()) fout.close(); //close previous file
		    fout.open((path+"/gntp_HEPevt_" + std::to_string(int(i/file_size))+".txt").c_str(), std::ofstream::out|std::ofstream::trunc);
		}


		tin->GetEntry(i);
		
		//locate index of outgoing photon
		int photon_index = -1;
		if(pdgf[0] == 22) photon_index = 0;
		else if(pdgf[1] == 22) photon_index = 1;

		int num_particles = 2;
		int status_code = 1;
		int first_mother = 0;
		int second_mother = 0;
		int first_daughter = 0;
		int second_daughter = 0;
		int pdg_code = pdgf[photon_index];    
		double momentum_x = pxf[photon_index];
		double momentum_y = pyf[photon_index];
		double momentum_z = pzf[photon_index];
		double energy = Ef[photon_index];
		double mass = 0.0;
		std::vector<double> random_position = sp.Sample3D(); //in cm
		double time = sp.SampleTime(); //in ns

		//std::cout << "nf: " << nf << std::endl;
		//for(int j =0; j<nf; ++j){
		//	std::cout << j << " : " <<pdgf[j] << ", " << pxf[j] << ", " << pyf[j] << ", " <<pzf[j] << ", " << Ef[j] << std::endl;
		//}

		fout << i << " " << num_particles << std::endl;
		fout << status_code << " " << pdg_code << " " << first_mother << " " << second_mother << " " << first_daughter << " "<< second_daughter << " " << momentum_x << " "<< momentum_y << " "<< momentum_z << " "<< energy << " "<< mass << " "<< random_position[0] << " "<< random_position[1] << " "<< random_position[2] << " " << time << std::endl;
		fout << 0 << " " << neu << " 0 0 0 0 0 0 0 " << Ev << " " <<Q2 << " 0 0 0 0" << std::endl;

		//fill in helper histograms
		h_x_coords->Fill(random_position[0]);
		h_y_coords->Fill(random_position[1]);
		h_z_coords->Fill(random_position[2]);
		h_xy_dist->Fill(random_position[0], random_position[1]);
		h_yz_dist->Fill(random_position[1], random_position[2]);
		h_xz_dist->Fill(random_position[0], random_position[2]);
	}


    	fout.close();


	//draw and save helper histograms
	std::string pdfout = path+"/HelperHistograms.pdf";
	TCanvas* c = new TCanvas();
  	gStyle->SetOptStat(0);
	c->Print((pdfout+"[").c_str());
	h_x_coords->Scale(1/h_x_coords->Integral());
	h_y_coords->Scale(1/h_y_coords->Integral());
	h_z_coords->Scale(1/h_z_coords->Integral());
	h_xy_dist->Scale(1/h_xy_dist->Integral());
	h_yz_dist->Scale(1/h_yz_dist->Integral());
	h_xz_dist->Scale(1/h_xz_dist->Integral());
	Helper1DHistDraw(c, h_x_coords, pdfout);
	Helper1DHistDraw(c, h_y_coords, pdfout);
	Helper1DHistDraw(c, h_z_coords, pdfout);
	Helper2DHistDraw(c, h_xy_dist, pdfout);
	Helper2DHistDraw(c, h_yz_dist, pdfout);
	Helper2DHistDraw(c, h_xz_dist, pdfout);
	c->Print((pdfout+"]").c_str());
}

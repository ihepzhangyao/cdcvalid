#include "MdcValidSingle/MdcValidSingle.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "EventModel/EventHeader.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "HepPDT/ParticleDataTable.hh"
#include "HepPDT/ParticleData.hh"
#include "GaudiKernel/IPartPropSvc.h"

#include "EventModel/Event.h"
#include "EvTimeEvent/RecEsTime.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
#include "MdcRecEvent/RecMdcHit.h"
#include "MdcRecEvent/RecMdcDedx.h"  
#include "MdcRecEvent/RecMdcDedxHit.h"
#include "MdcRawEvent/MdcDigi.h"
#include "Identifier/Identifier.h"
#include "Identifier/MdcID.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "McTruth/McParticle.h"
#include "McTruth/MdcMcHit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
//#include "VertexFit/Helix.h"
#include "TrackUtil/Helix.h"
#include "DstEvent/DstMdcKalTrack.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>

using namespace std;
using namespace Event;


MdcValidSingle::MdcValidSingle(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator)
{
  declareProperty("checkRecTrack",	m_checkRecTrack = true);
  declareProperty("checkKalTrack",	m_checkKalTrack = false);
  declareProperty("debug",		m_debug = 0);
  declareProperty("pid",		m_pid = 0);
  declareProperty("cutCharge",		m_cutCharge = true);
  declareProperty("cutPhi0",		m_cutPhi0 = -999.);
  declareProperty("cutCosTa",		m_cutCosTa= -999.);
  declareProperty("cutPt",		m_cutPt = -999.);
  declareProperty("cutDr",		m_cutDr = -999.);
  declareProperty("cutDz",		m_cutDz = -999.);
  declareProperty("cutHitRate",		m_cutHitRate = 0.8);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode MdcValidSingle::initialize(){

  MsgStream log(msgSvc(), name());
  StatusCode sc;
  log << MSG::INFO << "in initialize()" << endreq;
  m_totevt = 0;
  m_nNotMatch = 0;
  m_nMatchSingle = 0;
  m_nMatchMulti = 0;

  //sc = service("MdcGeomSvc", m_mdcGeomSvc);
  //if(!sc.isSuccess()){
  //  log << MSG::ERROR << " Could not initialize MdcGeomSvc" << endreq;
  //  return StatusCode::FAILURE;
  //}

  //Initailize MdcPrintSvc
  IMdcPrintSvc* imdcPrintSvc;
  sc = service ("MdcPrintSvc", imdcPrintSvc);
  m_mdcPrintSvc = dynamic_cast<MdcPrintSvc*> (imdcPrintSvc);
  if ( sc.isFailure() ){
    log << MSG::FATAL << "Could not load MdcPrintSvc!" << endreq;
    return StatusCode::FAILURE;
  }

  //Initailize magnetic filed
  sc = service ("MagneticFieldSvc",m_pIMF);
  if(sc!=StatusCode::SUCCESS) {
    log << MSG::ERROR << name() << " Unable to open magnetic field service"<<endreq;
  }
  m_bfield = new BField(m_pIMF);
  log << MSG::INFO <<name() << " field z = "<<m_bfield->bFieldNominal()<< endreq;

  Bz = m_bfield->bFieldNominal();

  sc = service("VertexDbSvc", m_vtxsvc);
  if(!sc.isSuccess()){
    log << MSG::ERROR << " Could not initialize VertexDbSvc" << endreq;
    return StatusCode::FAILURE;
  }

  NTuplePtr nt1(ntupleSvc(), "MdcValidSingle/diffMc");
  if ( nt1 ){
    m_ntupleDiffMc = nt1;
  } else {
    m_ntupleDiffMc = ntupleSvc()->book("MdcValidSingle/diffMc", CLID_ColumnWiseTuple, "diffMc");
    if(m_ntupleDiffMc){
      m_ntupleDiffMc->addItem("run",         m_runDiff);
      m_ntupleDiffMc->addItem("evt",         m_evtDiff);
      m_ntupleDiffMc->addItem("tkId",        m_tkId);
      m_ntupleDiffMc->addItem("nTk",         m_nTk);
      m_ntupleDiffMc->addItem("dPhi0",       m_dPhi0);
      m_ntupleDiffMc->addItem("dCosTa",      m_dCosTa);
      m_ntupleDiffMc->addItem("dPt",         m_dPt);
      m_ntupleDiffMc->addItem("dP",          m_dP);
      m_ntupleDiffMc->addItem("dPz",         m_dPz);
      m_ntupleDiffMc->addItem("dDr",         m_dDr);
      m_ntupleDiffMc->addItem("dDz",         m_dDz);
      m_ntupleDiffMc->addItem("dPhi0Kal",    m_dPhi0Kal);
      m_ntupleDiffMc->addItem("dCosTaKal",   m_dCosTaKal);
      m_ntupleDiffMc->addItem("dPtKal",      m_dPtKal);
      m_ntupleDiffMc->addItem("dDrKal",      m_dDrKal);
      m_ntupleDiffMc->addItem("dDzKal",      m_dDzKal);
    } else {   
      log << MSG::ERROR << "Cannot book tuple MdcValidSingle/diffMc" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  NTuplePtr nt(ntupleSvc(), "MdcValidSingle/evt");
  if ( nt ){
    m_ntupleEvt = nt;
  } else {
    m_ntupleEvt = ntupleSvc()->book("MdcValidSingle/evt", CLID_ColumnWiseTuple, "evt");
    if(m_ntupleEvt){
      m_ntupleEvt->addItem("run",        m_run);
      m_ntupleEvt->addItem("evt",        m_evt);
      m_ntupleEvt->addItem("t0",         m_t0);
      m_ntupleEvt->addItem("t0Stat",     m_t0Stat);
      m_ntupleEvt->addItem("t0Truth",    m_t0Truth);
      m_ntupleEvt->addItem("nTk",        m_nTrk);
      m_ntupleEvt->addItem("nKalTk",     m_nKalTrk);
      m_ntupleEvt->addItem("nLayerMc",   m_nLayerMc);
      m_ntupleEvt->addItem("nDigi",      m_nDigi);
      m_ntupleEvt->addItem("chi2",       m_chi2);
      m_ntupleEvt->addItem("ndof",       m_ndof);
      m_ntupleEvt->addItem("nMatchMdc",  m_nMatchMdc);
      m_ntupleEvt->addItem("nMatchKal",  m_nMatchKal);

      m_ntupleEvt->addItem("firstHitZ",  m_firstHitZ);
      m_ntupleEvt->addItem("lastHitZ",   m_lastHitZ);
      m_ntupleEvt->addItem("pidTruth",   m_pidTruth);
      m_ntupleEvt->addItem("costaTruth", m_costaTruth);
      m_ntupleEvt->addItem("phiTruth",   m_phi0Truth);
      m_ntupleEvt->addItem("drTruth",    m_drTruth);
      m_ntupleEvt->addItem("dzTruth",    m_dzTruth);
      m_ntupleEvt->addItem("ptTruth",    m_ptTruth);
      m_ntupleEvt->addItem("pTruth",     m_pTruth);
      m_ntupleEvt->addItem("fpTruth",    m_fpTruth);
      m_ntupleEvt->addItem("maxDp",      m_maxDp);
      m_ntupleEvt->addItem("qTruth",     m_qTruth);
      m_ntupleEvt->addItem("vertexTruth", 3, m_vertexTruth);

      m_ntupleEvt->addItem("charge",     m_charge);
      m_ntupleEvt->addItem("trkId",      m_trkId);
      m_ntupleEvt->addItem("dr",         m_dr);
      m_ntupleEvt->addItem("phi0",       m_phi0);
      m_ntupleEvt->addItem("kap",        m_kap);
      m_ntupleEvt->addItem("dz",         m_dz);
      m_ntupleEvt->addItem("tanl",       m_tanl);
      m_ntupleEvt->addItem("p",          m_p);
      m_ntupleEvt->addItem("pt",         m_pt);
      m_ntupleEvt->addItem("px",         m_px);
      m_ntupleEvt->addItem("py",         m_py);
      m_ntupleEvt->addItem("pz",         m_pz);
      m_ntupleEvt->addItem("fp",         m_fp);
      m_ntupleEvt->addItem("fpt",        m_fpt);
      m_ntupleEvt->addItem("fpx",        m_fpx);
      m_ntupleEvt->addItem("fpy",        m_fpy);
      m_ntupleEvt->addItem("fpz",        m_fpz);
      m_ntupleEvt->addItem("costa",      m_cosTheta);
      m_ntupleEvt->addItem("phi",        m_phi);
      m_ntupleEvt->addItem("r",          m_r);
      m_ntupleEvt->addItem("vr",         m_vr);
      m_ntupleEvt->addItem("vz",         m_vz);
      m_ntupleEvt->addItem("match",      m_match);
      m_ntupleEvt->addItem("stat",       m_tkStat);
      m_ntupleEvt->addItem("fiTerm",     m_fiTerm);
      m_ntupleEvt->addItem("multiRate",  m_multiRate);
      m_ntupleEvt->addItem("lastlay",    m_lastLayer);
      m_ntupleEvt->addItem("nhit",       m_nhit,0,7000);

      m_ntupleEvt->addIndexedItem("hitid",   m_nhit,     m_hitId);
      m_ntupleEvt->addIndexedItem("layer",   m_nhit,     m_hitLayer);
      m_ntupleEvt->addIndexedItem("wire",    m_nhit,     m_hitWire);
      m_ntupleEvt->addIndexedItem("hitstat", m_nhit,     m_hitStat);
      m_ntupleEvt->addIndexedItem("truthx",  m_nhit,     m_hitTruthX);
      m_ntupleEvt->addIndexedItem("truthy",  m_nhit,     m_hitTruthY);
      m_ntupleEvt->addIndexedItem("truthz",  m_nhit,     m_hitTruthZ);
      m_ntupleEvt->addIndexedItem("truthp",  m_nhit,     m_hitTruthP);
      m_ntupleEvt->addIndexedItem("driftDMc",m_nhit,     m_driftDMc);
      //m_ntupleEvt->addIndexedItem("measured",m_nHit,     m_hitMeasureD);
      //m_ntupleEvt->addIndexedItem("truthd",  m_nHit,     m_hitTruthD);
      
      //m_ntupleEvt->addItem("particleType", m_particleType, 0, 5);
      m_ntupleEvt->addItem("id_kal",     5, m_id_kal);
      m_ntupleEvt->addItem("stat_kal",   5,2, m_stat_kal);
      m_ntupleEvt->addItem("ndof_kal",   5,2, m_ndof_kal);
      m_ntupleEvt->addItem("chi2_kal",   5,2, m_chi2_kal);
      m_ntupleEvt->addItem("length_kal", 5, m_length_kal);
      m_ntupleEvt->addItem("tof_kal",    5, m_tof_kal);
      m_ntupleEvt->addItem("charge_kal", 5, m_charge_kal);
      m_ntupleEvt->addItem("dr_kal",     5, m_dr_kal);
      m_ntupleEvt->addItem("phi0_kal",   5, m_phi0_kal);
      m_ntupleEvt->addItem("kap_kal",    5, m_kap_kal);
      m_ntupleEvt->addItem("dz_kal",     5, m_dz_kal);
      m_ntupleEvt->addItem("tanl_kal",   5, m_tanl_kal);
      m_ntupleEvt->addItem("p_kal",      5, m_p_kal);
      m_ntupleEvt->addItem("pt_kal",     5, m_pt_kal);
      m_ntupleEvt->addItem("px_kal",     5, m_px_kal);
      m_ntupleEvt->addItem("py_kal",     5, m_py_kal);
      m_ntupleEvt->addItem("pz_kal",     5, m_pz_kal);
      m_ntupleEvt->addItem("theta_kal",  5, m_theta_kal);
      m_ntupleEvt->addItem("phi_kal",    5, m_phi_kal);
      m_ntupleEvt->addItem("vr_kal",     5, m_vr_kal);
      m_ntupleEvt->addItem("vz_kal",     5, m_vz_kal);
      m_ntupleEvt->addItem("nhit_kal",   5, m_nhit_kal);
      m_ntupleEvt->addItem("zp_kal",     5, m_zp_kal);
      m_ntupleEvt->addItem("zpt_kal",    5, m_zpt_kal);
      m_ntupleEvt->addItem("fp_kal",     5, m_fp_kal);
      m_ntupleEvt->addItem("fpt_kal",    5, m_fpt_kal);
      m_ntupleEvt->addItem("lp_kal",     5, m_lp_kal);
      m_ntupleEvt->addItem("lpt_kal",    5, m_lpt_kal);
      m_ntupleEvt->addItem("match_kal",  5, m_match_kal);
      m_ntupleEvt->addItem("lastlay_kal",5, m_lastlay_kal);//not implemented
    } else {   
      log << MSG::ERROR << "Cannot book tuple MdcValidSingle/evt" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode MdcValidSingle::execute() {

  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

  t_trigger=0;

  SmartDataPtr<MdcDigiCol> mdcDigiCol(eventSvc(),"/Event/Digi/MdcDigiCol");
  if (!mdcDigiCol) {
    log << MSG::FATAL << "Could not find event" << endreq;
    return -999;
  }
  if(mdcDigiCol->size()<5) return StatusCode::SUCCESS;
  m_nDigi = mdcDigiCol->size();

  if(0 == (m_totevt % 10000)) cout << name()<< " execute event " << m_totevt << endl;
  m_totevt++;

  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
  if (!eventHeader) {
    log << MSG::FATAL << "Could not find EventHeader" << endreq;
    return StatusCode::FAILURE;
  }
  t_run = eventHeader->runNumber();
  t_evt = eventHeader->eventNumber();
  m_run = t_run;
  m_evt = t_evt;

  std::cout<<name()<<" run "<<t_run<<" evt "<<t_evt<<"   "<<std::endl;
  int nMatchMdc = -999;
  int nMatchKal = -999;

  if(m_debug>2) m_mdcPrintSvc->printDigi(); 

  //-- mc particle and hit
  SmartDataPtr<McParticleCol> mcParticleCol(eventSvc(),"/Event/MC/McParticleCol");
  if (!mcParticleCol) {
    log << MSG::ERROR << "Could not find McParticle" << endreq;
    return StatusCode::FAILURE; 
  }
  fillMcTruth();
 

  //-- rec track
  if(m_checkRecTrack){
    SmartDataPtr<RecMdcTrackCol> recMdcTrackCol(eventSvc(), "/Event/Recon/RecMdcTrackCol");
    if(!recMdcTrackCol){
      log << MSG::FATAL<< "Could not get RecMdcTrackCol" << endreq;
      return StatusCode::FAILURE;
    }
    m_nTrk = recMdcTrackCol->size();
    RecMdcTrackCol::iterator it = recMdcTrackCol->begin();
    for (;it!= recMdcTrackCol->end();it++){
      fillMdcTrack(*it); 
    }
  }//end of fill Mdc track

  if(m_checkKalTrack){
    //-- kalman track
    SmartDataPtr<RecMdcKalTrackCol> recMdcKalTrackCol(eventSvc(),"/Event/Recon/RecMdcKalTrackCol");
    if(!recMdcKalTrackCol){
      log << MSG::FATAL<< "Could not retrieve RecMdcKalTrackCol" << endreq;
      return StatusCode::FAILURE;
    }
    m_nKalTrk = recMdcKalTrackCol->size();
    std::cout<<name()<<" evt  "<< t_evt <<" nKal=  "<<recMdcKalTrackCol->size()<<std::endl;
    RecMdcKalTrackCol::iterator itKal = recMdcKalTrackCol->begin();
    for (;itKal!= recMdcKalTrackCol->end();itKal++){
      fillKalTrack(*itKal); 
    }
    if(m_debug>2) m_mdcPrintSvc->printRecMdcTrackCol(); 
  }//end of fill Kalman track

  if(!(m_checkRecTrack && m_checkKalTrack)){
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
    if (!evtRecEvent) {
      log << MSG::FATAL << "Could not find EvtRecEvent" << endreq;
      return StatusCode::FAILURE;
    }
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
    if (!evtRecTrkCol) {
      log << MSG::FATAL << "Could not find EvtRecTrackCol" << endreq;
      return StatusCode::FAILURE;
    }

    m_nTrk = evtRecEvent->totalCharged();
    //int tkIdMatched = matchRecTrack(nMatchMdc,nMatchKal);
    //if( tkIdMatched!=-999 ){ 
    //  for(int i = 0; i < evtRecEvent->totalCharged(); i++) {
    //    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    //    RecMdcTrack *mdctrk = (*itTrk)->mdcTrack();
    //    if(mdctrk->trackId() != tkIdMatched) continue;
    //    RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
    //    fillRecMdcTrack(mdctrk, mdcKalTrk); 
    //  }
    //}
    //RecMdcTrack* mdcTrk = NULL;
    //RecMdcKalTrack* mdcKalTrk = NULL;

    //int tkIdMatched = matchRecTrack(nMatchMdc,nMatchKal);
    //if( tkIdMatched!=-999 ){
    //  //-- mdc track
    //  RecMdcTrackCol::iterator it = recMdcTrackCol->begin();
    //  for (;it!= recMdcTrackCol->end();it++){
    //    if(tkIdMatched==(*it)->trackId()) mdcTrk = (*it);
    //  }
    //}
  }
  //int tkIdMatched = matchRecTrack(nMatchMdc,nMatchKal);
  //m_nMatchMdc = nMatchMdc;
  //m_nMatchKal = nMatchKal;



  m_ntupleEvt->write();
  return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode MdcValidSingle::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in finalize()" << endreq;
  cout <<name() << " total event: " << m_totevt <<" matchSingle "<<m_nMatchSingle<<" notMatch "<<m_nNotMatch <<" matchMulti "<<m_nMatchMulti<< endl;

  return StatusCode::SUCCESS;
}


void MdcValidSingle::getVertex(HepVector a, HepSymMatrix Ea, double& vr, double& vz){
  Hep3Vector xorigin(0,0,0);
  double* dbv = m_vtxsvc->PrimaryVertex();
  double* vv = m_vtxsvc->SigmaPrimaryVertex();
  xorigin.setX(dbv[0]);
  xorigin.setY(dbv[1]);
  xorigin.setZ(dbv[2]);

  /*
     HepPoint3D point0(0.,0.,0.);
     VFHelix helixip(point0,a,Ea);
     HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
     helixip.pivot(IP);
     HepVector vecipa = helixip.a();
     vr = vecipa[0];
     vz = vecipa[3];
   */
}

//------------------------------------------------------------------
//given pid of single primary particle
//fill truth info, match track by param cut
//count matched mdc and kalman track number
// if MC track not valid return -999, else return nMatchMdc
//------------------------------------------------------------------
int MdcValidSingle::matchRecTrack(int& nMatchMdc, int& nMatchKal){
  nMatchMdc = -999;
  nMatchKal = -999;
  MsgStream log(msgSvc(), name());
  StatusCode sc;
  std::vector<RecMdcTrack*> mdcTrackVec;
  std::vector<RecMdcKalTrack*> mdcKalTrackVec;

  if(!(m_checkRecTrack||m_checkKalTrack)){
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
    if (!evtRecEvent) {
      log << MSG::FATAL << "Could not find EvtRecEvent" << endreq;
      return -999; 
    }
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
    if (!evtRecTrkCol) {
      log << MSG::FATAL << "Could not find EvtRecTrackCol" << endreq;
      return -999; 
    }


    t_nTk = evtRecEvent->totalCharged();
    for(int i = 0; i < evtRecEvent->totalCharged(); i++){
      EvtRecTrackIterator itTrk =evtRecTrkCol->begin() + i;
      if(!(*itTrk)->isMdcTrackValid()) continue;
      RecMdcTrack* tk = (*itTrk)->mdcTrack();
      mdcTrackVec.push_back(tk);
      RecMdcKalTrack *kalTk = (*itTrk)->mdcKalTrack();
      if(kalTk!=NULL) mdcKalTrackVec.push_back(kalTk);
    }
  }else{
    //-- mdc track
    SmartDataPtr<RecMdcTrackCol> recMdcTrackCol(eventSvc(), "/Event/Recon/RecMdcTrackCol");
    if(!recMdcTrackCol && m_checkRecTrack){
      log << MSG::FATAL<< "Could not get RecMdcTrackCol" << endreq;
      return -999;
    }
    RecMdcTrackCol::iterator it = recMdcTrackCol->begin();
    for (;it!= recMdcTrackCol->end();it++){
      mdcTrackVec.push_back(*it);
    }
    //-- kalman track
    SmartDataPtr<RecMdcKalTrackCol> recMdcKalTrackCol(eventSvc(),"/Event/Recon/RecMdcKalTrackCol");
    if(!recMdcKalTrackCol && m_checkKalTrack){
      log << MSG::FATAL<< "Could not retrieve RecMdcKalTrackCol" << endreq;
      return -999;
    }
    RecMdcKalTrackCol::iterator itKal = recMdcKalTrackCol->begin();
    for (;itKal!= recMdcKalTrackCol->end();itKal++){
      mdcKalTrackVec.push_back(*itKal);
    }
  }

  //-----------------------------------
  //loop over mdc and kalman tracks, match by track param
  //count matched track number
  //-----------------------------------
  nMatchMdc=0;
  nMatchKal=0;
  //m_nTrk = mdcTrackVec.size();
  for (int iTk = 0; iTk< mdcTrackVec.size(); iTk++){
    bool mdcMatched = true;
    bool kalMatched = true;
    //-- MDC track
    //RecMdcTrack* tk = (*itTrk)->mdcTrack();
    RecMdcTrack* tk = mdcTrackVec[iTk];
    int tkId = tk->trackId();
    double dr    = tk->helix(0);
    double phi0  = tk->helix(1);
    double cpa   = tk->helix(2);
    double dz    = tk->helix(3);
    double tanl  = tk->helix(4);
    int    q     = tk->charge();

    if(m_debug>2) std::cout<<__FILE__<<" "<<__LINE__<<" MdcTrack helix "<<tk->helix()<<std::endl;
    double dPhi0,dCosTa,dPt,dP,dPz,dDr,dDz= -999.;
    dPhi0 = phi0 - t_mcHelix[1];//mchelix->phi0();
    dCosTa= cos(tk->theta()) - t_costaTruth;
    dP  = tk->p() - t_pTruth;
    dPt = tk->pxy() - t_ptTruth;
    dPz = tk->pz() - t_pzTruth;
    dDr = dr - t_mcHelix[0];//mchelix->dr();
    dDz = dz - t_mcHelix[3];//mchelix->dz();

    //std::cout<<__FILE__<<" "<<__LINE__<<" ptMdc  "<<tk->pxy()<<" mc pt "<<mchelix->pt()<<std::endl;
    m_runDiff = t_run;
    m_evtDiff = t_evt;
    m_tkId = tkId;
    m_nTk  = t_nTk;

    m_dPhi0 = dPhi0;
    m_dCosTa= dCosTa;
    m_dPt = dPt;
    m_dPz = dPz;
    m_dP = dP;
    m_dDr = dDr;
    m_dDz = dDz;

    //std::cout<<__FILE__<<" "<<__LINE__<<" recHits  "<<std::endl;
    int nMatchHit = 0;
    HitRefVec hitList = tk->getVecHits();
    HitRefVec::iterator hitIt = hitList.begin();
    for(;hitIt!=hitList.end();hitIt++){
      RecMdcHit *h = *hitIt;
      if(!h) continue;
      Identifier id = h->getMdcId();
      int layer = MdcID::layer(id);
      int wire = MdcID::wire(id);
      //std::cout<<"("<<layer<<","<<wire<<") t_McTkId  "<<t_McTkId<<" truthmap "<<truthmap[layer][wire]<<std::endl;
      if(t_truthmap[layer][wire]==t_McTkId) nMatchHit++;
    }

    if(nMatchHit< m_cutHitRate * hitList.size()) {
      if(m_debug>1) std::cout<<t_run<<" "<<t_evt<<" cut by nMatchHit! nMatchHit="<<nMatchHit<<" nHit "<<hitList.size()<<" rate= "<<(float) nMatchHit/((float) hitList.size())<<" cut "<<m_cutHitRate<<" nMcHit "<<t_nMcHit<<std::endl;
      mdcMatched = false;
    }

    RecMdcKalTrack* kalTk = NULL;//TEMP yzhang
    //-- kalman track
    double dPhi0Kal,dCosTaKal,dPtKal,dDrKal,dDzKal = -999;
    int qKal=0;
    //RecMdcKalTrack *kalTk = (*itTrk)->mdcKalTrack();
    if(mdcKalTrackVec.size()>0){
      RecMdcKalTrack* kalTk = mdcKalTrackVec[iTk];
      if(NULL != kalTk){
	//enum PidType { null = -1, electron = 0, muon = 1, pion = 2, kaon = 3, proton = 4 };
	int ipid = 0;
	if(fabs(m_pid) == 13) ipid = 1; 
	DstMdcKalTrack::PidType pidType = DstMdcKalTrack::PidType(ipid);
	kalTk->setPidType(pidType);
	dPhi0Kal = fabs(kalTk->fi0() - t_mcHelix[1]/*mchelix->phi0()*/ + CLHEP::pi*2.);
	while(dPhi0Kal>CLHEP::pi) dPhi0Kal -= CLHEP::twopi;
	dCosTaKal= cos(kalTk->theta()) - t_mcHelix[4];//mchelix->cosTheta(); 
	dPtKal = kalTk->pxy() - t_ptTruth;
	//std::cout<<__FILE__<<" "<<__LINE__<<" ptKal  "<<kalTk->pxy()<<" mc pt "<<mchelix->pt()<<std::endl;
	dDrKal = kalTk->dr() - t_mcHelix[0];//mchelix->dr();
	dDzKal = kalTk->dz() - t_mcHelix[3];//mchelix->dz();

	qKal = kalTk->charge();
	m_dPhi0Kal = dPhi0Kal;
	m_dCosTaKal= dCosTaKal;
	m_dPtKal = dPtKal;
	m_dDrKal = dDrKal;
	m_dDzKal = dDzKal;
      }
    }

    m_ntupleDiffMc->write();

    //--cut mdc
    if(m_cutCharge && t_qTruth != q) {
      if(m_debug>1) std::cout<<t_run<<" "<<t_evt<<" cut by charge! q_rec "<<q<<" q_mc "<<t_qTruth <<std::endl;
      mdcMatched = false;
    }
    if(m_cutPhi0>-999. && fabs(dPhi0)>m_cutPhi0){
      if(m_debug>1) std::cout<<t_run<<" "<<t_evt<<" cut by phi0! dPhi0 "<<dPhi0<<" cutPhi0  "<<m_cutPhi0<<std::endl;
      mdcMatched = false;
    }
    if(m_cutCosTa> -999. && fabs(dCosTa)>m_cutCosTa){
      if(m_debug>1) std::cout<<t_run<<" "<<t_evt<<" cut by cos theta! dCosTa "<<dCosTa<<" cutCosTa "<<m_cutCosTa<<std::endl;
      mdcMatched = false;
    }
    if(m_cutPt > -999. && fabs(dPt)>m_cutPt){
      if(m_debug>1) std::cout<<t_run<<" "<<t_evt<<" cut by pt! dPt "<<dPt<<" cutPt "<<m_cutPt<<std::endl;
      mdcMatched = false;
    }
    if(m_cutDr > -999. && fabs(dDr)>m_cutDr){
      if(m_debug>1) std::cout<<t_run<<" "<<t_evt<<" cut by dr! dDr "<<dDr<<" cutDr "<<m_cutDr<<std::endl;
      mdcMatched = false;
    }
    if(m_cutDz > -999. && fabs(dDz)>m_cutDz){
      if(m_debug>1) std::cout<<t_run<<" "<<t_evt<<" cut by dz! dDz "<<dDz<<" cutDz "<<m_cutDz<<std::endl;
      mdcMatched = false;
    }

    if(mdcKalTrackVec.size()<=0) {kalMatched=false;}
    //if(mdcKalTrackVec.size()>0){
    //  //-- cut kalman
    //  if(t_qTruth != qKal) {
    //    if(m_debug>1) std::cout<<t_run<<" "<<t_evt<<" cut by charge kalman! q_rec "<<qKal<<" q_mc"<<t_qTruth<<std::endl;
    //    kalMatched = false;
    //  }
    //  if(m_cutPhi0>-999. && fabs(dPhi0Kal)>m_cutPhi0){
    //    if(m_debug>1) std::cout<<t_run<<" "<<t_evt<<" cut by phi0 kalman! dPhi "<<dPhi0Kal<<" cutPhi0  "<<m_cutPhi0<<std::endl;
    //    kalMatched = false;
    //  }
    //  if(m_cutCosTa> -999. && fabs(dCosTaKal)>m_cutCosTa){
    //    if(m_debug>1) std::cout<<t_run<<" "<<t_evt<<" cut by cos theta kalman! dCosTa "<<dCosTaKal<<" cutCosTa "<<m_cutCosTa<<std::endl;
    //    kalMatched = false;
    //  }
    //  if(m_cutPt > -999. && fabs(dPtKal)>m_cutPt){
    //    if(m_debug>1) std::cout<<t_run<<" "<<t_evt<<" cut by pt kalman! dPt "<<dPtKal<<" cutPt "<<m_cutPt<<std::endl;
    //    kalMatched = false;
    //  }
    //  if(m_cutDr > -999. && fabs(dDrKal)>m_cutDr){
    //    if(m_debug>1) std::cout<<t_run<<" "<<t_evt<<" cut by dr kalman! dDr "<<dDrKal<<" cutDr "<<m_cutDr<<std::endl;
    //    kalMatched = false;
    //  }
    //  if(m_cutDz > -999. && fabs(dDzKal)>m_cutDz){
    //    if(m_debug>1) std::cout<<t_run<<" "<<t_evt<<" cut by dz kalman! dDz "<<dDzKal<<" cutDz "<<m_cutDz<<std::endl;
    //    kalMatched = false;
    //  }
    //}
    if(mdcMatched) {
      nMatchMdc++;
      //matchedTkId = tkId;
    }
    if(kalMatched) {
      nMatchKal++;
    }

    fillRecMdcTrack(tk, mdcMatched, kalTk, kalMatched); 


  }//end of loop mdcTrackVec
  if(m_debug>2) std::cout<<"nMatchMdc="<<nMatchMdc<<" nMatchKal="<<nMatchKal<<"   "<<std::endl;

  //if(NULL != mchelix) delete mchelix;

  if(nMatchMdc==1) m_nMatchSingle++;
  if(nMatchMdc==0) m_nNotMatch++;
  if(nMatchMdc>1) m_nMatchMulti++;
  if(nMatchMdc!=1 && m_debug>0) {
    std::cout<<" run "<<t_run<<" evt "<<t_evt<<" nMatchMdc="<<nMatchMdc<<std::endl;
    return -999;
  }
  if(m_checkKalTrack && nMatchKal!=1 && m_debug>0) {
    std::cout<<" run "<<t_run<<" evt "<<t_evt<<" nMatchKal="<<nMatchKal<<std::endl;
  }

  return nMatchMdc;
}//end of matchRecTrack

void MdcValidSingle::fillRecMdcTrack(RecMdcTrack* tk, bool mdcMatched, RecMdcKalTrack* kalTk, bool kalMatched){
  m_charge = tk->charge();
  m_trkId  = tk->trackId();
  m_dr     = tk->helix(0);
  m_phi0   = tk->helix(1);
  m_kap    = tk->helix(2);
  m_dz     = tk->helix(3);
  m_tanl   = tk->helix(4);
  m_chi2   = tk->chi2();
  m_ndof   = tk->ndof();
  m_p      = tk->p();
  m_pt     = tk->pxy();
  m_px     = tk->px();
  m_py     = tk->py();
  m_pz     = tk->pz();
  m_cosTheta  = cos(tk->theta());
  m_phi    = tk->phi();
  //double Bz = m_bfield->bFieldNominal();
  m_r      = -333.567/(Bz*(tk->helix(2)));
  double vr,vz;
  //getVertex(tk->helix(),tk->err(),vr,vz);
  m_vr     = vr;
  m_vz     = vz;
  m_tkStat = tk->stat();
  m_fiTerm = tk->getFiTerm();
  m_match  = mdcMatched;
  m_lastLayer = tk->lastLayer();
  m_nhit = tk->getVecHits().size();


  if(NULL != kalTk){
    for (int i=0; i<5; i++) {
      //enum PidType { null = -1, electron = 0, muon = 1, pion = 2, kaon = 3, proton = 4 };
      DstMdcKalTrack::PidType pidType = DstMdcKalTrack::PidType(i);
      kalTk->setPidType(pidType);
      for(int jj=0;jj<2;jj++){
	for(int kk=0;kk<5;kk++){
	  m_stat_kal[jj][kk] = kalTk->getStat(jj,kk);
	  m_ndof_kal[jj][kk] = kalTk->getNdf(jj,kk);
	  m_chi2_kal[jj][kk] = kalTk->getChisq(jj,kk);
	}
      } 
      m_charge_kal[i]     = kalTk->charge();
      m_dr_kal[i]         = kalTk->dr();
      m_phi0_kal[i]       = kalTk->fi0();
      m_kap_kal[i]        = kalTk->kappa();
      m_dz_kal[i]         = kalTk->dz();
      m_tanl_kal[i]       = kalTk->tanl();
      m_zpt_kal[i]          = kalTk->p();
      m_pt_kal[i]         = kalTk->pxy();
      m_px_kal[i]         = kalTk->px();
      m_py_kal[i]         = kalTk->py();
      m_pz_kal[i]         = kalTk->pz();
      m_theta_kal[i]      = kalTk->theta();
      m_phi_kal[i]        = kalTk->phi();
      double vrKal,vzKal;
      //getVertex(kalTk->helix(),kalTk->err(),vrKal,vzKal);
      m_vr_kal[i]         = vrKal;
      m_vz_kal[i]         = vzKal;
      m_match_kal[i]      = kalMatched;
      m_zp_kal[i]         = sqrt(1+kalTk->getZHelixE()[i]*kalTk->getZHelixE()[i])/kalTk->getZHelixE()[i];
      m_zpt_kal[i]        = 1/kalTk->getZHelixE()[i];
      HepVector fHelix = kalTk->getFHelixE();
      double kappa = fHelix[2];
      double tanl = fHelix[4];
      double pxy = 0;
      if(kappa!=0) pxy = 1.0/fabs(kappa);
      double pz = pxy * tanl;
      double ptot = sqrt(pxy*pxy+pz*pz);
      m_fp_kal[i]         = ptot;
      //std::cout<<__FILE__<<"   "<<__LINE__<<" fp  "<<ptot<<" fp_kal "<<m_fp_kal[i]<<std::endl;
      m_fpt_kal[i]        = pxy;
      m_lastlay_kal[i]    = kalTk->getLastLayer(0,i);//not implemented
    }
  }

}//end of fillRecMdcTrack

void MdcValidSingle::fillMdcTrack(RecMdcTrack* tk){
  m_charge = tk->charge();
  m_trkId  = tk->trackId();
  m_dr     = tk->helix(0);
  m_phi0   = tk->helix(1);
  m_kap    = tk->helix(2);
  m_dz     = tk->helix(3);
  m_tanl   = tk->helix(4);
  m_chi2   = tk->chi2();
  m_ndof   = tk->ndof();
  m_p      = tk->p();
  m_pt     = tk->pxy();
  m_px     = tk->px();
  m_py     = tk->py();
  m_pz     = tk->pz();
  m_cosTheta  = cos(tk->theta());
  m_phi    = tk->phi();
  //double Bz = m_bfield->bFieldNominal();
  m_r      = -333.567/(Bz*(tk->helix(2)));
  double vr,vz;
  //getVertex(tk->helix(),tk->err(),vr,vz);
  m_vr     = vr;
  m_vz     = vz;
  m_tkStat = tk->stat();
  m_fiTerm = tk->getFiTerm();
  m_lastLayer = tk->lastLayer();
  m_nhit = tk->getVecHits().size();
}

void MdcValidSingle::fillKalTrack(RecMdcKalTrack* kalTk){
  if(NULL == kalTk) return;
  for (int i=0; i<5; i++) {
    //enum PidType { null = -1, electron = 0, muon = 1, pion = 2, kaon = 3, proton = 4 };
    DstMdcKalTrack::PidType pidType = DstMdcKalTrack::PidType(i);
    kalTk->setPidType(pidType);
    for(int jj=0;jj<2;jj++){
      for(int kk=0;kk<5;kk++){
	m_stat_kal[jj][kk] = kalTk->getStat(jj,kk);
	m_ndof_kal[jj][kk] = kalTk->getNdf(jj,kk);
	m_chi2_kal[jj][kk] = kalTk->getChisq(jj,kk);
      }
    } 
    m_charge_kal[i]     = kalTk->charge();
    m_dr_kal[i]         = kalTk->dr();
    m_phi0_kal[i]       = kalTk->fi0();
    m_kap_kal[i]        = kalTk->kappa();
    m_dz_kal[i]         = kalTk->dz();
    m_tanl_kal[i]       = kalTk->tanl();
    m_p_kal[i]          = kalTk->p();
    m_pt_kal[i]         = kalTk->pxy();
    m_px_kal[i]         = kalTk->px();
    m_py_kal[i]         = kalTk->py();
    m_pz_kal[i]         = kalTk->pz();
    m_theta_kal[i]      = kalTk->theta();
    m_phi_kal[i]        = kalTk->phi();
    m_nhit_kal[i]       = kalTk->ndof()+5;
    double vrKal,vzKal;
    //getVertex(kalTk->helix(),kalTk->err(),vrKal,vzKal);
    m_vr_kal[i]         = vrKal;
    m_vz_kal[i]         = vzKal;
    HepVector zHelix = kalTk->getZHelixE();
    HepVector fHelix = kalTk->getFHelixE();
    HepVector lHelix = kalTk->getLHelixE();

    double kappa = fHelix[2];
    double tanl = fHelix[4];
    double pxy = 0;
    if(kappa!=0) pxy = 1.0/fabs(kappa);
    double pz = pxy * tanl;
    double ptot = sqrt(pxy*pxy+pz*pz);

    m_fp_kal[i]         = ptot;
      //std::cout<<__FILE__<<"   "<<__LINE__<<" fp  "<<ptot<<" fp_kal "<<m_fp_kal[i]<<std::endl;
    m_fpt_kal[i]        = pxy;
  }
}//end of fillKalTrack



void MdcValidSingle::fillMcTruth(){
  SmartDataPtr<McParticleCol> mcParticleCol(eventSvc(),"/Event/MC/McParticleCol");
  if (!mcParticleCol) {
    cout <<  "Could not find McParticle" << endl;
    return; 
  }
  t_t0Truth=-1.;
  t_qTruth = 0;
  int itk = 0;
  int t_pidTruth = -999;
  t_McTkId = -999;
  t_nTrkMC=0;
  //Helix* mchelix;
  //Hep3Vector true_mom;
  t_mcHelix[0]=t_mcHelix[1]=t_mcHelix[2]=t_mcHelix[3]=t_mcHelix[4]=-999;
  HepPoint3D vertexTruth;

  McParticleCol::iterator iter_mc = mcParticleCol->begin();
  for (;iter_mc != mcParticleCol->end(); iter_mc++ ) {
    if(!(*iter_mc)->primaryParticle()) continue;
    t_pidTruth = (*iter_mc)->particleProperty();

    if(m_debug>2) cout<< "tk "<<itk<< " pid="<< t_pidTruth<< endl;
    if((m_pid!=0) && (t_pidTruth != m_pid)){ continue; }

    t_t0Truth=(*iter_mc)->initialPosition().t();

    if(itk>=10) break;

    int pid = t_pidTruth; 
    if( pid == 0 ) { 
      //cout << "Wrong particle id" <<endl;
      //continue;
      t_qTruth = 1;//for changedgeantino
      if(m_debug>2) std::cout << " particle changedgeantino charge= "<<t_qTruth<<std::endl;
    }else{
      IPartPropSvc* p_PartPropSvc;
      static const bool CREATEIFNOTTHERE(true);
      StatusCode PartPropStatus = service("PartPropSvc", p_PartPropSvc, CREATEIFNOTTHERE);
      if (!PartPropStatus.isSuccess() || 0 == p_PartPropSvc) {
	std::cout<< " Could not initialize Particle Properties Service" << std::endl;
      }
      HepPDT::ParticleDataTable* p_particleTable = p_PartPropSvc->PDT(); 
      std::string name;
      if( p_particleTable->particle(pid) ){
	name = p_particleTable->particle(pid)->name();
	t_qTruth = p_particleTable->particle(pid)->charge();
      }else if( p_particleTable->particle(-pid) ){
	name = "anti " + p_particleTable->particle(-pid)->name();
	t_qTruth = (-1)*p_particleTable->particle(-pid)->charge();
      }
      if(m_debug>2) std::cout << " particle "<<name<<" charge= "<<t_qTruth<<std::endl;
    }

    t_McTkId = (*iter_mc)->trackIndex();

    t_true_mom = (*iter_mc)->initialFourMomentum().vect();
    const CLHEP::HepLorentzVector& posIni = (*iter_mc)->initialPosition()/10.;
    vertexTruth.setX(posIni.x());
    vertexTruth.setY(posIni.y());
    vertexTruth.setZ(posIni.z());
    Helix mchelix(vertexTruth,t_true_mom,t_qTruth);
    if(m_debug>1)std::cout<<"\n --Helix before: "<<mchelix.dr()<<","<<mchelix.phi0()<<","<<mchelix.kappa()<<","<<mchelix.dz()<<","<<mchelix.tanl()<<std::endl;
    mchelix.pivot(HepPoint3D(0.,0.,0.));
    if(m_debug>1)std::cout<<"\n --Helix after:"<<mchelix.dr()<<","<<mchelix.phi0()<<","<<mchelix.kappa()<<","<<mchelix.dz()<<","<<mchelix.tanl()<<std::endl;
    t_mcHelix[0]=mchelix.dr();
    t_mcHelix[1]=mchelix.phi0();
    t_mcHelix[2]=mchelix.kappa();
    t_mcHelix[3]=mchelix.dz();
    t_mcHelix[4]=mchelix.tanl();

    t_ptTruth    = fabs(1./t_mcHelix[2]);
    t_pTruth     = sqrt(1+t_mcHelix[4]*t_mcHelix[4])/t_mcHelix[2];
    t_pzTruth    = t_pTruth * t_mcHelix[4];

    if(m_debug>2){
      std::cout<<"Truth tk "<<itk<<" pid "<<pid<<" charge "<<t_qTruth
	<< " helix "<< mchelix.a()<<" t_pTruth "<<t_pTruth<<" mom_p "<<t_true_mom.mag()
	<<" t_ptTruth "<<t_ptTruth<<" mom_pt "<<t_true_mom.perp()<<std::endl;
    }

    t_nTrkMC++;
    itk++;
  }//end of loop over mcParticleCol

  if(t_nTrkMC!=1) {
    std::cout<<"WARNING run "<<t_run<<" evt "<<t_evt<<" not single event. nTrkMc="<<t_nTrkMC<<std::endl;
    return;
  }else{
    if(m_debug>2) std::cout<<"nTrkMc=1"<<std::endl;
  }

  for(int i=0;i<NLAYER;i++)for(int j=0;j<MAXCELL;j++) {
    t_mcDrift[i][j]=-999;
    t_mcX[i][j] =t_mcY[i][j]= t_mcZ[i][j]=-999;
    t_hitTruthP[i][j]=-999;
  }
  double firstHitZ = 9999.;
  double lastHitZ = 9999.;
  double maxDp = 0.;
  //------------------Retrieve MC truth MdcMcHit------------
  SmartDataPtr<Event::MdcMcHitCol> mcMdcMcHitCol(eventSvc(),"/Event/MC/MdcMcHitCol"); 
  if (!mcMdcMcHitCol) {
    cout <<  "Could not find MdcMcHit" << endl; 
  }else{
    Event::MdcMcHitCol::iterator iter_mchit = mcMdcMcHitCol->begin();
    int ihit = 0;
    for (;iter_mchit != mcMdcMcHitCol->end(); iter_mchit++,ihit++ ) {
      if((*iter_mchit)->getPositionFlag()==-999) continue;
      const Identifier id= (*iter_mchit)->identify();
      int layer = MdcID::layer(id);
      int wire = MdcID::wire(id);
      t_mcDrift[layer][wire] = (*iter_mchit)->getDriftDistance()/10.;  //drift in MC.
      t_mcX[layer][wire] = (*iter_mchit)->getPositionX()/10.;
      t_mcY[layer][wire] = (*iter_mchit)->getPositionY()/10.;
      t_mcZ[layer][wire] = (*iter_mchit)->getPositionZ()/10.;
      if(layer==999 && (wire=998||wire==999)) {
	t_trigger++;
	//if(m_debug>0)
	std::cout<<"trigger   "<<t_trigger<<" pos  ("<<t_mcX<<","<<t_mcY<<","<<t_mcZ<<")"<<std::endl;
      }
      //double mcX = (*iter_mchit)->getPositionX()/10.;
      //double mcY = (*iter_mchit)->getPositionY()/10.;
      double mcZ = (*iter_mchit)->getPositionZ()/10.;
      if(ihit == 0 && layer ==0){
	firstHitZ = mcZ;
      }

      if(ihit+1 == mcMdcMcHitCol->size() && layer ==0){
	lastHitZ = mcZ;
      }
    }

    ihit=0;
    //get track mom at this hit position
    iter_mchit = mcMdcMcHitCol->begin();
    for (;iter_mchit != mcMdcMcHitCol->end(); iter_mchit++ ) {
      //std::cout<<__FILE__<<"   "<<__LINE__<<"   "<<(*iter_mchit)->getPositionFlag()<<std::endl;
      if((*iter_mchit)->getPositionFlag()!=-999) continue;
      Identifier mdcid = (*iter_mchit)->identify();
      int layerId = MdcID::layer(mdcid);
      int wireId = MdcID::wire(mdcid);
      double fpx = (*iter_mchit)->getPositionX()/1000.;//MeV 2 GeV
      double fpy = (*iter_mchit)->getPositionY()/1000.;//MeV 2 GeV
      double fpz = (*iter_mchit)->getPositionZ()/1000.;//MeV 2 GeV
      double p = sqrt(fpx*fpx + fpy*fpy + fpz*fpz);
      t_hitTruthP[layerId][wireId] = p;
      if(ihit==0&&layerId==0) t_fpTruth = p;
      if(fabs(-1.*p - t_pTruth)> fabs(maxDp)) maxDp = -1.*p - t_pTruth;
      if(m_debug>2)std::cout<<"("<<layerId<<","<<wireId<<") p  "<<p<<" fpTruth "<<t_fpTruth<<" maxDp "<<maxDp<<std::endl;
      ihit++;
    }
  }//end of loop over MdcMcHitCol

  int nHits=0; int nSt=0;
  double maxTrackLength = -999;
  Hep3Vector trueMomFinal;
  HepPoint3D positionFinal;
  HitRefVec  hitRefVec;
  m_firstHitZ = firstHitZ;
  m_lastHitZ = lastHitZ;

  m_pidTruth   = t_pidTruth;
  //m_drTruth    = mchelix[0];
  //m_phi0Truth  = mchelix[1];
  //m_dzTruth    = mchelix[3];
  //double tanlTruth = mchelix[4];
  m_drTruth    = t_mcHelix[0];//mchelix->dr();
  m_phi0Truth  = t_mcHelix[1];//mchelix->phi0();;
  m_dzTruth    = t_mcHelix[3];//mchelix->dz();
  double tanlTruth = t_mcHelix[4];//mchelix->tanl();
  //m_tanlTruth  = tanlTruth;
  t_costaTruth = tanlTruth / sqrt(1+tanlTruth*tanlTruth);
  m_costaTruth = t_costaTruth;

  m_ptTruth    = t_true_mom.perp();//t_ptTruth
  m_pTruth     = t_true_mom.mag();//t_pTruth;
  m_fpTruth    = t_fpTruth;
  m_maxDp      = maxDp;
  m_vertexTruth[0] = vertexTruth.x();
  m_vertexTruth[1] = vertexTruth.y();
  m_vertexTruth[2] = vertexTruth.z();
  m_qTruth     = t_qTruth;


  //-- event start time
  SmartDataPtr<RecEsTimeCol> aevtimeCol(eventSvc(),"/Event/Recon/RecEsTimeCol");
  if (aevtimeCol) {
    RecEsTimeCol::iterator iter_evt = aevtimeCol->begin();
    for (;iter_evt!=aevtimeCol->end();iter_evt++){
      m_t0 = (*iter_evt)->getTest();
      m_t0Stat = (*iter_evt)->getStat();
      m_t0Truth = t_t0Truth;
    }
  }else{
    cout <<  "Could not find RecEsTimeCol" << endl;
    m_t0 = 0;
  } 

  t_nMcHit=0;

  //mdc digi 
  for(int i=0;i<NLAYER;i++) for(int j=0;j<MAXCELL;j++) t_truthmap[i][j]=-1;
  SmartDataPtr<MdcDigiCol> mdcDigiCol(eventSvc(),"/Event/Digi/MdcDigiCol");
  if (!mdcDigiCol) {
    cout <<  "Could not find event" << endl;
    return;
  }else{
    MdcDigiCol::iterator iter_digi = mdcDigiCol->begin();
    for(; iter_digi != mdcDigiCol->end(); iter_digi++) {
      int layer = MdcID::layer((*iter_digi)->identify());
      int wire = MdcID::wire((*iter_digi)->identify());
      int tkId = (*iter_digi)->getTrackIndex();
      if(tkId>=1000) tkId-=1000;
      t_truthmap[layer][wire]= tkId;
    }
  }
  bool nLayerMcPassed[NLAYER];
  int nLayerMc = 0;
  for (int i=0;i<NLAYER;i++){ nLayerMcPassed[i]=false;}
  for(int i=0;i<NLAYER;i++) for(int j=0;j<MAXCELL;j++) {
    if(t_truthmap[i][j]==t_McTkId) {
      t_nMcHit++;
      nLayerMcPassed[i]=true;
      //std::cout<<"("<<i<<","<<j<<") "<<t_McTkId<<std::endl;
    }
  }
  for (int i=0;i<NLAYER;i++){ 
    if(nLayerMcPassed[i]) nLayerMc++;
  }
  m_nLayerMc = nLayerMc;
}


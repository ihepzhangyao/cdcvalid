#ifndef MDCVALIDSINGLE_H
#define MDCVALIDSINGLE_H

#include "GaudiKernel/Algorithm.h"
//#include "MdcGeomSvc/IMdcGeomSvc.h"
//#include "MdcGeomSvc/MdcGeomSvc.h"
//#include "MdcGeomSvc/MdcGeoWire.h"
//#include "MdcGeomSvc/MdcGeoLayer.h"
#include "McTruth/McParticle.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/INTupleSvc.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "MagneticField/IMagneticFieldSvc.h"
#include "BField/BField.h"
#include <string>
#include "MdcPrintSvc/MdcPrintSvc.h"

class MdcValidSingle:public Algorithm {
public:
     MdcValidSingle(const std::string& name, ISvcLocator* pSvcLocator);
     StatusCode initialize();
     StatusCode execute();
     StatusCode finalize();

     void getVertex(HepVector a, HepSymMatrix Ea, double& vr, double& vz);
     int matchRecTrack(int& nMatchMdc, int& nMatchKal);
     void fillRecMdcTrack(RecMdcTrack* tk, bool mdcMatched, RecMdcKalTrack* kalTk, bool kalMatched);
     void fillMdcTrack(RecMdcTrack* tk);
     void fillKalTrack(RecMdcKalTrack* kalTk);
     HepVector calcHelix(HepPoint3D position, Hep3Vector momentum, int q , double Bz);
     void fillMcTruth();
  
private:
     int    t_trigger;
     double t_mcHelix[5];
     static const int NLAYER=18;//43 besIII
     static const int MAXCELL=300;//288 besIII
     static const int MAXWIERNUMBER=8000;//6796 besIII
     int m_totevt;
     int t_nTotTk;
     int m_nNotMatch;
     int m_nMatchSingle;
     int m_nMatchMulti;

     bool   m_cutCharge;
     double m_cutPhi0;
     double m_cutCosTa;
     double m_cutPt;
     double m_cutDr;
     double m_cutDz;
     double m_cutHitRate;

     int    t_nTrkMC,t_McTkId,t_nTk;
     double t_t0Truth;
     double m_checkRecTrack;
     double m_checkKalTrack;
     int    m_pid;
     int    m_debug;
     double t_fpTruth;
     long   t_run,t_evt;
     int    t_nMcHit,t_qTruth;
     double t_costaTruth;
     double t_mcDrift[NLAYER][MAXCELL];
     int    t_truthmap[NLAYER][MAXCELL];
     double t_mcX[NLAYER][MAXCELL];
     double t_mcY[NLAYER][MAXCELL];
     double t_mcZ[NLAYER][MAXCELL];
     double t_hitTruthP[NLAYER][MAXCELL];
     Hep3Vector t_true_mom;
     double t_pTruth,t_ptTruth,t_pzTruth;

     //IMdcGeomSvc*         m_mdcGeomSvc;
     MdcPrintSvc*          m_mdcPrintSvc;
     IVertexDbSvc*        m_vtxsvc;
     BField*	    	  m_bfield;
     IMagneticFieldSvc*   m_pIMF;
     double		  Bz;

     NTuple::Tuple*       m_ntupleDiffMc;
     NTuple::Item<int>   m_runDiff;
     NTuple::Item<int>   m_evtDiff;
     NTuple::Item<int>   m_tkId;
     NTuple::Item<int>   m_nTk;
     NTuple::Item<double> m_dPhi0;
     NTuple::Item<double> m_dCosTa;
     NTuple::Item<double> m_dPt;
     NTuple::Item<double> m_dPz;
     NTuple::Item<double> m_dP;
     NTuple::Item<double> m_dDr;
     NTuple::Item<double> m_dDz;
     NTuple::Item<double> m_dPhi0Kal;
     NTuple::Item<double> m_dCosTaKal;
     NTuple::Item<double> m_dPtKal;
     NTuple::Item<double> m_dPzKal;
     NTuple::Item<double> m_dPKal;
     NTuple::Item<double> m_dDrKal;
     NTuple::Item<double> m_dDzKal;

     NTuple::Tuple*       m_ntupleEvt;
     NTuple::Item<int>   m_run;
     NTuple::Item<int>   m_evt;
     NTuple::Item<int>   m_nTrk;
     NTuple::Item<int>   m_nKalTrk;
     NTuple::Item<int>   m_nLayerMc;
     NTuple::Item<int>   m_nDigi;
     NTuple::Item<double> m_t0;
     NTuple::Item<int>   m_t0Stat;
     NTuple::Item<double> m_t0Truth;

     NTuple::Item<double> m_chi2;
     NTuple::Item<int>   m_ndof;
     NTuple::Item<int>   m_nMatchMdc;
     NTuple::Item<int>   m_nMatchKal;

     NTuple::Item<double> m_firstHitZ;
     NTuple::Item<double> m_lastHitZ;
     NTuple::Item<double> m_pidTruth;
     NTuple::Item<double> m_costaTruth;
     NTuple::Item<double> m_tanlTruth;
     NTuple::Item<double> m_phi0Truth;
     NTuple::Item<double> m_drTruth;
     NTuple::Item<double> m_dzTruth;
     NTuple::Item<double> m_ptTruth;
     NTuple::Item<double> m_pTruth,m_fpTruth,m_pzTruth,m_qTruth;
     NTuple::Item<double> m_maxDp;
     NTuple::Array<double> m_vertexTruth;

     NTuple::Item<int>     m_trkId,m_charge,m_tkStat;
     NTuple::Item<double>  m_dr,m_phi0,m_kap,m_dz,m_tanl;
     NTuple::Item<double>  m_p,m_pt,m_px,m_py,m_pz;
     NTuple::Item<double>  m_fp,m_fpt,m_fpx,m_fpy,m_fpz;
     NTuple::Item<double>  m_cosTheta,m_phi,m_r,m_vr,m_vz,m_fiTerm;
     NTuple::Item<int>     m_nhit,m_lastLayer,m_multiRate,m_match;
     NTuple::Array<int>    m_hitId,m_hitLayer,m_hitWire,m_hitStat;
     NTuple::Array<double> m_hitTruthX,m_hitTruthY,m_hitTruthZ,m_hitTruthP;

     //mdc hits
     NTuple::Array<int>     m_layer,m_wire,m_ambig,m_flaglr,m_Tdc,m_Adc;
     NTuple::Array<double>  m_resid,m_driftD,m_driftDMc,m_driftT,m_sigma,m_doca;
     NTuple::Array<double>  m_fltLen,m_entra,m_zhit,m_chi2add;

     //kalman track
     NTuple::Matrix<int>    m_stat_kal,m_ndof_kal;
     NTuple::Matrix<double> m_chi2_kal;
     NTuple::Array<int>     m_id_kal,m_length_kal,m_charge_kal,m_nhit_kal,m_match_kal,m_lastlay_kal;
     NTuple::Array<double>  m_tof_kal, m_dr_kal, m_phi0_kal, m_kap_kal, m_dz_kal, m_tanl_kal;
     NTuple::Array<double>  m_p_kal, m_pt_kal, m_px_kal, m_py_kal, m_pz_kal;
     NTuple::Array<double>  m_theta_kal, m_phi_kal, m_vr_kal, m_vz_kal;
     NTuple::Array<double>  m_zp_kal,m_zpt_kal,m_fp_kal, m_fpt_kal, m_lp_kal, m_lpt_kal;

};
#endif /* MDCVALIDSINGLE_H */

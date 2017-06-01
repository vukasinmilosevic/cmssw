#include "RecoTauTag/HLTProducers/interface/L1JetsTMatching.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/JetReco/interface/PFJet.h"

//
// class decleration
//
using namespace reco   ;
using namespace std    ;
using namespace edm    ;
using namespace trigger;


template< typename T>
std::pair<std::vector<T>,std::vector<T>> categorise(const std::vector<T>& pfMatchedJets, double pt1, double pt2, double Mjj)
{
    std::pair<std::vector<T>,std::vector<T>> output;
    unsigned int i1 = 0;
    unsigned int i2 = 0;
    double mjj = 0;
    if (pfMatchedJets.size()>1){
        for (unsigned int i = 0; i < pfMatchedJets.size()-1; i++)
            for (unsigned int j = i+1; j < pfMatchedJets.size(); j++)
            {
                const T &  myJet1 = (pfMatchedJets)[i];
                const T &  myJet2 = (pfMatchedJets)[j];
                
                const double mjj_test = (myJet1.p4()+myJet2.p4()).M();
                
                if (mjj_test > mjj){
                    
                    mjj =mjj_test;
                    i1 = i;
                    i2 = j;
                }
            }
        
        const T &  myJet1 = (pfMatchedJets)[i1];
        const T &  myJet2 = (pfMatchedJets)[i2];
        
        if ((mjj > Mjj) && (myJet1.pt() >= pt1) && (myJet2.pt() > pt2) )
        {
            
            output.first.push_back(myJet1);
            output.first.push_back(myJet2);
            
        }
        
        if ((mjj > Mjj) && (myJet1.pt() < pt1) && (myJet1.pt() > pt2) && (myJet2.pt() > pt2))
        {
            
            const T &  myJetTest = (pfMatchedJets)[0];
            if (myJetTest.pt()>pt1){
                output.second.push_back(myJet1);
                output.second.push_back(myJet2);
                output.second.push_back(myJetTest);
                
            }
        }
        
    }
    
    return output;
    
}
template< typename T>
L1JetsTMatching<T>::L1JetsTMatching(const edm::ParameterSet& iConfig):
  jetSrc_    ( consumes<std::vector<T>>                     (iConfig.getParameter<InputTag>("JetSrc"      ) ) ),
  jetTrigger_( consumes<trigger::TriggerFilterObjectWithRefs>(iConfig.getParameter<InputTag>("L1JetTrigger") ) ),
  pt1_Min   ( iConfig.getParameter<double>("pt1_Min")),
  pt2_Min   ( iConfig.getParameter<double>("pt2_Min")),
  mjj_Min   ( iConfig.getParameter<double>("mjj_Min"))
{  
  produces<std::vector<T>>("TwoJets");
  produces<std::vector<T>>("ThreeJets");
}
template< typename T>
L1JetsTMatching<T>::~L1JetsTMatching(){ }

template< typename T>
void L1JetsTMatching<T>::produce(edm::StreamID iSId, edm::Event& iEvent, const edm::EventSetup& iES) const
{
    
  unique_ptr<std::vector<T>> pfMatchedJets(new std::vector<T>);
    std::pair<std::vector<T>,std::vector<T>> output;
    
  double deltaR    = 1.0;
  double matchingR = 0.5;
  
  // Getting HLT jets to be matched
  edm::Handle<std::vector<T> > pfJets;
  iEvent.getByToken( jetSrc_, pfJets );
        
  edm::Handle<trigger::TriggerFilterObjectWithRefs> l1TriggeredJets;
  iEvent.getByToken(jetTrigger_,l1TriggeredJets);
                
  //l1t::TauVectorRef jetCandRefVec;
  l1t::JetVectorRef jetCandRefVec;
  l1TriggeredJets->getObjects( trigger::TriggerL1Jet,jetCandRefVec);

  math::XYZPoint a(0.,0.,0.);
        
 //std::cout<<"PFsize= "<<pfJets->size()<<endl<<" L1size= "<<jetCandRefVec.size()<<std::endl;
 for(unsigned int iJet = 0; iJet < pfJets->size(); iJet++){
    for(unsigned int iL1Jet = 0; iL1Jet < jetCandRefVec.size(); iL1Jet++){
      // Find the relative L2pfJets, to see if it has been reconstructed
      const T &  myJet = (*pfJets)[iJet];
    //  if ((iJet<3) && (iL1Jet==0))  std::cout<<myJet.p4().Pt()<<" ";
      deltaR = ROOT::Math::VectorUtil::DeltaR2(myJet.p4().Vect(), (jetCandRefVec[iL1Jet]->p4()).Vect());
      if(deltaR < matchingR ) {
        pfMatchedJets->push_back(myJet);
        break;
      }
    }
  }  
   
    output= categorise(*pfMatchedJets,pt1_Min,pt2_Min, mjj_Min);
    unique_ptr<std::vector<T>> output1(new std::vector<T>(output.first));
    unique_ptr<std::vector<T>> output2(new std::vector<T>(output.second));
    
    iEvent.put(std::move(output1),"TwoJets");
    iEvent.put(std::move(output2),"ThreeJets");

}
/*template< typename T>
void L1JetsTMatching<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("L1JetTrigger", edm::InputTag("hltL1DiJetVBF"))->setComment("Name of trigger filter"    );
  desc.add<edm::InputTag>("JetSrc"      , edm::InputTag("hltAK4PFJetsTightIDCorrected"))->setComment("Input collection of PFJets");
  desc.add<double>       ("pt1_Min",95.0)->setComment("Minimal pT1 of PFJets to match");
  desc.add<double>       ("pt2_Min",35.0)->setComment("Minimal pT2 of PFJets to match");
  desc.add<double>       ("mjj_Min",650.0)->setComment("Minimal mjj of matched PFjets");
  descriptions.setComment("This module produces collection of PFJetss matched to L1 Taus / Jets passing a HLT filter (Only p4 and vertex of returned PFJetss are set).");
  descriptions.add       ("L1JetsTMatching",desc);
}*/

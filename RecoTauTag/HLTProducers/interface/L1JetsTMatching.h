#ifndef RecoTauTag_HLTProducers_L1JetsTMatching_h
#define RecoTauTag_HLTProducers_L1JetsTMatching_h

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"


#include <map>
#include <vector>

template< typename T>
class L1JetsTMatching: public edm::global::EDProducer<> {
 public:
  explicit L1JetsTMatching(const edm::ParameterSet&);
  ~L1JetsTMatching();
  virtual void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  //static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
    
  const edm::EDGetTokenT<std::vector<T>> jetSrc_;
  const edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> jetTrigger_;
  const double pt1_Min;
  const double pt2_Min;
  const double mjj_Min;
    

};
#endif

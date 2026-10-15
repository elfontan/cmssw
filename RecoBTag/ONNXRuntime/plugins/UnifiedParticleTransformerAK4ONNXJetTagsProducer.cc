#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/BTauReco/interface/UnifiedParticleTransformerAK4TagInfo.h"
#include "DataFormats/BTauReco/interface/UnifiedParticleTransformerAK4Features.h"

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"
#include <algorithm>
#include <iostream>
#include <cassert>

#include "RecoBTag/ONNXRuntime/interface/tensor_fillers.h"
#include "RecoBTag/ONNXRuntime/interface/tensor_configs.h"
#include "RecoBTag/ONNXRuntime/interface/UParTNoPID.h"

using namespace cms::Ort;

class UnifiedParticleTransformerAK4ONNXJetTagsProducer : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
public:
  explicit UnifiedParticleTransformerAK4ONNXJetTagsProducer(const edm::ParameterSet&, const ONNXRuntime*);
  ~UnifiedParticleTransformerAK4ONNXJetTagsProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

  static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const ONNXRuntime*);

private:
  typedef std::vector<reco::UnifiedParticleTransformerAK4TagInfo> TagInfoCollection;
  typedef reco::JetTagCollection JetTagCollection;

  void produce(edm::Event&, const edm::EventSetup&) override;

  void make_inputs(btagbtvdeep::UnifiedParticleTransformerAK4Features features);
  void make_test_inputs();
  void print_inputs();
  void get_input_sizes(const reco::FeaturesTagInfo<btagbtvdeep::UnifiedParticleTransformerAK4Features> taginfo);

  const edm::EDGetTokenT<TagInfoCollection> src_;
  std::vector<std::string> flav_names_;
  std::vector<std::string> input_names_;
  bool use_dynamic_axes_ = false;
  bool usePID_ = true; // default for UParT v00 and v01 models
  bool verbose_ = false;
  std::vector<std::string> output_names_;

  unsigned n_cpf_;
  unsigned n_features_cpf_ = UparT::N_InputFeatures.at(UparT::kChargedCandidates);
  unsigned n_lt_;
  unsigned n_npf_;
  unsigned n_sv_;
  std::vector<unsigned> input_sizes_;
  std::vector<std::vector<int64_t>> input_shapes_;  // shapes of each input group (-1 for dynamic axis)

  // hold the input data
  FloatArrays data_;

  // test mode config: "off" | "N" | "N+M"
  std::string testMode_ = "off"; // For ONNX & dyn. shapes integration tests
  unsigned testN_ = 1;
  unsigned testM_ = 0;
};

UnifiedParticleTransformerAK4ONNXJetTagsProducer::UnifiedParticleTransformerAK4ONNXJetTagsProducer(
    const edm::ParameterSet& iConfig, const ONNXRuntime* cache)
    : src_(consumes<TagInfoCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      flav_names_(iConfig.getParameter<std::vector<std::string>>("flav_names")),
      input_names_(iConfig.getParameter<std::vector<std::string>>("input_names")),
        use_dynamic_axes_(
          (iConfig.getParameter<edm::FileInPath>("model_path").fullPath().find("V01") != std::string::npos)
                ||
          (iConfig.getParameter<edm::FileInPath>("model_path").fullPath().find("V02") != std::string::npos)
              ),
        usePID_(iConfig.getParameter<bool>("usePID")),
        verbose_(iConfig.getParameter<bool>("verbose")),
      output_names_(iConfig.getParameter<std::vector<std::string>>("output_names")),
      testMode_(iConfig.getParameter<std::string>("testMode")),
      testN_(iConfig.getParameter<unsigned>("testN")),
      testM_(iConfig.getParameter<unsigned>("testM")) {
  // get output names from flav_names
  n_features_cpf_ = usePID_ ? UparT::N_InputFeatures.at(UparT::kChargedCandidates) : UparTNoPID::n_charged_features;
  for (const auto& flav_name : flav_names_) {
    produces<JetTagCollection>(flav_name);
  }

  if (testMode_ != "off" && testMode_ != "N" && testMode_ != "N+M" && testMode_ != "print") {
    throw cms::Exception("Configuration") << "Invalid value for testMode: '" << testMode_
                                          << "'. Allowed values are: off, N, N+M, print";
  }
}

void UnifiedParticleTransformerAK4ONNXJetTagsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // pfUnifiedParticleTransformerAK4JetTags
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("pfUnifiedParticleTransformerAK4TagInfos"));
  desc.add<std::vector<std::string>>(
      "input_names", {"input_1", "input_2", "input_3", "input_4", "input_5", "input_6", "input_7", "input_8"});
  desc.add<edm::FileInPath>("model_path",
                            edm::FileInPath("RecoBTag/Combined/data/UParTAK4/PUPPI/V01/modelfile/model.onnx"));
  desc.add<std::vector<std::string>>("output_names", {"softmax"});
  desc.add<bool>("usePID", true);
  desc.add<bool>("verbose", false);
  desc.add<std::string>("testMode", "off");  // "off" | "N" | "N+M"
  desc.add<unsigned>("testN", 1);
  desc.add<unsigned>("testM", 0);
  desc.add<std::vector<std::string>>(
      "flav_names",
      std::vector<std::string>{"probb",        "probbb",       "problepb",     "probc",         "probs",
                               "probu",        "probd",        "probg",        "probele",       "probmu",
                               "probtaup1h0p", "probtaup1h1p", "probtaup1h2p", "probtaup3h0p",  "probtaup3h1p",
                               "probtaum1h0p", "probtaum1h1p", "probtaum1h2p", "probtaum3h0p",  "probtaum3h1p",
                               "ptcorr",       "ptreshigh",    "ptreslow",     "ptnu",          "probemudata",
                               "probemumc",    "probdimudata", "probdimumc",   "probmutaudata", "probmutaumc"});

  descriptions.add("pfUnifiedParticleTransformerAK4JetTags", desc);
}

std::unique_ptr<ONNXRuntime> UnifiedParticleTransformerAK4ONNXJetTagsProducer::initializeGlobalCache(
    const edm::ParameterSet& iConfig) {
  return std::make_unique<ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("model_path").fullPath());
}

void UnifiedParticleTransformerAK4ONNXJetTagsProducer::globalEndJob(const ONNXRuntime* cache) {}

void UnifiedParticleTransformerAK4ONNXJetTagsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<TagInfoCollection> tag_infos;
  iEvent.getByToken(src_, tag_infos);

  // initialize output collection
  std::vector<std::unique_ptr<JetTagCollection>> output_tags;
  if (!tag_infos->empty()) {
    auto jet_ref = tag_infos->begin()->jet();
    auto ref2prod = edm::makeRefToBaseProdFrom(jet_ref, iEvent);
    for (std::size_t i = 0; i < flav_names_.size(); i++) {
      output_tags.emplace_back(std::make_unique<JetTagCollection>(ref2prod));
    }
  } else {
    for (std::size_t i = 0; i < flav_names_.size(); i++) {
      output_tags.emplace_back(std::make_unique<JetTagCollection>());
    }
  }

  for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
    const auto& taginfo = (*tag_infos)[jet_n];
    std::vector<float> outputs(flav_names_.size(), -1.0);
    if (taginfo.features().is_filled) {
      get_input_sizes(taginfo);

      // run prediction with dynamic batch size per event
      input_shapes_ = {{(int64_t)1, (int64_t)n_cpf_, (int64_t)n_features_cpf_},
                       {(int64_t)1, (int64_t)n_lt_, (int64_t)UparT::N_InputFeatures.at(UparT::kLostTracks)},
                       {(int64_t)1, (int64_t)n_npf_, (int64_t)UparT::N_InputFeatures.at(UparT::kNeutralCandidates)},
                       {(int64_t)1, (int64_t)n_sv_, (int64_t)UparT::N_InputFeatures.at(UparT::kVertices)},
                       {(int64_t)1, (int64_t)n_cpf_, (int64_t)UparT::N_InputFeatures.at(UparT::kChargedCandidates4Vec)},
                       {(int64_t)1, (int64_t)n_lt_, (int64_t)UparT::N_InputFeatures.at(UparT::kLostTracks4Vec)},
                       {(int64_t)1, (int64_t)n_npf_, (int64_t)UparT::N_InputFeatures.at(UparT::kNeutralCandidates4Vec)},
                       {(int64_t)1, (int64_t)n_sv_, (int64_t)UparT::N_InputFeatures.at(UparT::kVertices4Vec)}};

      if (verbose_){
        // print input names
        std::cout << "Input names: ";
        for (const auto& name : input_names_) {
          std::cout << name << " ";
        }
        std::cout << std::endl;
        // print output names
        std::cout << "Output names: ";
        for (const auto& name : output_names_) {
          std::cout << name << " ";
        }      std::cout << std::endl;

      }

      if (testMode_ == "print") {
        print_inputs();
      }
      // Always run the model to produce real scores (also for 'print' mode)
      outputs = globalCache()->run(input_names_, data_, input_shapes_, output_names_, 1)[0];
      assert(outputs.size() == flav_names_.size());
    }

    const auto& jet_ref = tag_infos->at(jet_n).jet();
    for (std::size_t flav_n = 0; flav_n < flav_names_.size(); flav_n++) {
      // print flav name and output value
      if (verbose_) {
        std::cout << "Jet " << jet_ref.key() << ": " << flav_names_[flav_n] << " = " << outputs[flav_n] << std::endl;
      }
      (*(output_tags[flav_n]))[jet_ref] = outputs[flav_n];
    }
  }

  // put into the event
  for (std::size_t flav_n = 0; flav_n < flav_names_.size(); ++flav_n) {
    iEvent.put(std::move(output_tags[flav_n]), flav_names_[flav_n]);
  }
}

void UnifiedParticleTransformerAK4ONNXJetTagsProducer::get_input_sizes(
    const reco::FeaturesTagInfo<btagbtvdeep::UnifiedParticleTransformerAK4Features> taginfo) {
  const auto& features = taginfo.features();

  const unsigned max_cpf = use_dynamic_axes_ && !usePID_ ? 25u : UparT::n_cpf_accept;
  const unsigned max_lt = UparT::n_lt_accept;
  const unsigned max_npf = use_dynamic_axes_ && !usePID_ ? 20u : UparT::n_npf_accept;
  const unsigned max_sv = use_dynamic_axes_ && !usePID_ ? 6u : UparT::n_sv_accept;

  if (testMode_ == "N" || testMode_ == "N+M") {
    const unsigned fakeLen = (testMode_ == "N") ? testN_ : (testN_ + testM_);
    n_cpf_ = std::min(std::max<unsigned>(1u, fakeLen), max_cpf);
    n_lt_ = std::min(std::max<unsigned>(1u, fakeLen), max_lt);
    n_npf_ = std::min(std::max<unsigned>(1u, fakeLen), max_npf);
    n_sv_ = std::min(std::max<unsigned>(1u, fakeLen), max_sv);

    if (verbose_) {
      std::cout << "[TEST MODE] mode=" << testMode_ << " N=" << testN_ << " M=" << testM_
                << " -> n_cpf=" << n_cpf_ << ", n_lt=" << n_lt_ << ", n_npf=" << n_npf_ << ", n_sv=" << n_sv_
                << std::endl;
    }
  } else if (use_dynamic_axes_) {
    // Use actual sizes for dynamic axes version
    n_cpf_ = std::clamp<unsigned>(features.c_pf_features.size(), 1u, max_cpf);
    n_lt_ = std::clamp<unsigned>(features.lt_features.size(), 1u, max_lt);
    n_npf_ = std::clamp<unsigned>(features.n_pf_features.size(), 1u, max_npf);
    n_sv_ = std::clamp<unsigned>(features.sv_features.size(), 1u, max_sv);

  } else {
    // Use fixed sizes for original version
    n_cpf_ = (unsigned int)UparT::n_cpf_accept;
    n_lt_ = (unsigned int)UparT::n_lt_accept;
    n_npf_ = (unsigned int)UparT::n_npf_accept;
    n_sv_ = (unsigned int)UparT::n_sv_accept;
  }

    if (verbose_) {
                std::cout << "Jet " << taginfo.jet().key() << ": n_cpf = " << n_cpf_ << ", n_lt = " << n_lt_
                          << ", n_npf = " << n_npf_ << ", n_sv = " << n_sv_ << std::endl;
              }

  input_sizes_ = {
      n_cpf_ * UparT::N_InputFeatures.at(UparT::kChargedCandidates),
      n_lt_ * UparT::N_InputFeatures.at(UparT::kLostTracks),
      n_npf_ * UparT::N_InputFeatures.at(UparT::kNeutralCandidates),
      n_sv_ * UparT::N_InputFeatures.at(UparT::kVertices),
      n_cpf_ * UparT::N_InputFeatures.at(UparT::kChargedCandidates4Vec),
      n_lt_ * UparT::N_InputFeatures.at(UparT::kLostTracks4Vec),
      n_npf_ * UparT::N_InputFeatures.at(UparT::kNeutralCandidates4Vec),
      n_sv_ * UparT::N_InputFeatures.at(UparT::kVertices4Vec),
  };

  if (verbose_) {
    std::cout << "Input sizes: ";
    for (const auto& size : input_sizes_) {
      std::cout << size << " ";
    }
    std::cout << std::endl;
  }

  // init data storage
  data_.clear();
  for (const auto& len : input_sizes_) {
    data_.emplace_back(1 * len, 0);
  }
  if (testMode_ == "off" || testMode_ == "print") {
    make_inputs(features);
  } else {
    make_test_inputs();
  }
  // Retain the release tensor fillers; adapt only the custom NoPID charged input.
  if (!usePID_) {
    UparTNoPID::removeChargedPIDFeatures(data_[UparT::kChargedCandidates]);
    input_sizes_[UparT::kChargedCandidates] = data_[UparT::kChargedCandidates].size();
  }
}

void UnifiedParticleTransformerAK4ONNXJetTagsProducer::make_inputs(
    btagbtvdeep::UnifiedParticleTransformerAK4Features features) {
  const float* start = nullptr;
  unsigned offset = 0;

  auto max_c_pf_n = std::min(features.c_pf_features.size(), (std::size_t)n_cpf_);
  auto max_lt_n = std::min(features.lt_features.size(), (std::size_t)n_lt_);
  auto max_n_pf_n = std::min(features.n_pf_features.size(), (std::size_t)n_npf_);
  auto max_sv_n = std::min(features.sv_features.size(), (std::size_t)n_sv_);
  // c_pf candidates
  UparT_tensor_filler(data_, UparT::kChargedCandidates, features.c_pf_features, max_c_pf_n, start, offset);
  // lt candidates
  UparT_tensor_filler(data_, UparT::kLostTracks, features.lt_features, max_lt_n, start, offset);
  // n_pf candidates
  UparT_tensor_filler(data_, UparT::kNeutralCandidates, features.n_pf_features, max_n_pf_n, start, offset);
  // sv candidates
  UparT_tensor_filler(data_, UparT::kVertices, features.sv_features, max_sv_n, start, offset);
  // cpf pairwise features (4-vectors)
  UparT_tensor_filler(data_, UparT::kChargedCandidates4Vec, features.c_pf_features, max_c_pf_n, start, offset);
  // lt pairwise features (4-vectors) specific case requiring (pt,eta,phi,e)
  UparT_tensor_filler(data_, UparT::kLostTracks4Vec, features.lt_features, max_lt_n, start, offset);
  // npf pairwise features (4-vectors)
  UparT_tensor_filler(data_, UparT::kNeutralCandidates4Vec, features.n_pf_features, max_n_pf_n, start, offset);
  // sv pairwise features (4-vectors)
  UparT_tensor_filler(data_, UparT::kVertices4Vec, features.sv_features, max_sv_n, start, offset);
}

void UnifiedParticleTransformerAK4ONNXJetTagsProducer::make_test_inputs() {
  // data_ already zero-initialized in get_input_sizes()

  const unsigned ones_cpf = (testMode_ == "N") ? n_cpf_ : std::min(testN_, n_cpf_);
  const unsigned ones_lt = (testMode_ == "N") ? n_lt_ : std::min(testN_, n_lt_);
  const unsigned ones_npf = (testMode_ == "N") ? n_npf_ : std::min(testN_, n_npf_);
  const unsigned ones_sv = (testMode_ == "N") ? n_sv_ : std::min(testN_, n_sv_);

  auto fillLeadingOnes = [&](unsigned idx, unsigned nObj, unsigned nFeat, unsigned nOnesObj) {
    for (unsigned i = 0; i < std::min(nObj, nOnesObj); ++i) {
      const unsigned base = i * nFeat;
      for (unsigned f = 0; f < nFeat; ++f) {
        data_[idx][base + f] = 1.f;
      }
    }
  };

  fillLeadingOnes(UparT::kChargedCandidates, n_cpf_, UparT::N_InputFeatures.at(UparT::kChargedCandidates), ones_cpf);
  fillLeadingOnes(UparT::kLostTracks, n_lt_, UparT::N_InputFeatures.at(UparT::kLostTracks), ones_lt);
  fillLeadingOnes(UparT::kNeutralCandidates, n_npf_, UparT::N_InputFeatures.at(UparT::kNeutralCandidates), ones_npf);
  fillLeadingOnes(UparT::kVertices, n_sv_, UparT::N_InputFeatures.at(UparT::kVertices), ones_sv);

  fillLeadingOnes(UparT::kChargedCandidates4Vec, n_cpf_, UparT::N_InputFeatures.at(UparT::kChargedCandidates4Vec), ones_cpf);
  fillLeadingOnes(UparT::kLostTracks4Vec, n_lt_, UparT::N_InputFeatures.at(UparT::kLostTracks4Vec), ones_lt);
  fillLeadingOnes(UparT::kNeutralCandidates4Vec, n_npf_, UparT::N_InputFeatures.at(UparT::kNeutralCandidates4Vec), ones_npf);
  fillLeadingOnes(UparT::kVertices4Vec, n_sv_, UparT::N_InputFeatures.at(UparT::kVertices4Vec), ones_sv);

  if (verbose_) {
    std::cout << "[TEST MODE] Filled synthetic inputs: mode=" << testMode_ << " N=" << testN_ << " M=" << testM_
              << std::endl;
  }
}

void UnifiedParticleTransformerAK4ONNXJetTagsProducer::print_inputs() {
  // Print each input group as a Python nested list with shape info.
  for (std::size_t idx = 0; idx < input_names_.size(); ++idx) {
    const auto& name = input_names_[idx];
    const auto& shape = input_shapes_[idx];
    const int64_t batch = shape[0];
    const int64_t nObj = shape[1];
    const int64_t nFeat = shape[2];

    std::cout << "# Input '" << name << "' shape=(" << batch << "," << nObj << "," << nFeat << ")\n";
    std::cout << name << " = [" << std::endl;
    // batch dim (expected 1)
    for (int64_t b = 0; b < batch; ++b) {
      std::cout << "  [" << std::endl;
      // each object
      const auto& flat = data_[idx];
      for (int64_t i = 0; i < nObj; ++i) {
        std::cout << "    [";
        for (int64_t f = 0; f < nFeat; ++f) {
          const std::size_t pos = static_cast<std::size_t>(i * nFeat + f);
          // ensure within vector bounds
          float v = 0.f;
          if (pos < flat.size()) v = flat[pos];
          std::cout << v;
          if (f + 1 < nFeat) std::cout << ", ";
        }
        std::cout << "]";
        if (i + 1 < nObj) std::cout << ",";
        std::cout << std::endl;
      }
      std::cout << "  ]" << std::endl;
    }
    std::cout << "]" << std::endl << std::endl;
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(UnifiedParticleTransformerAK4ONNXJetTagsProducer);

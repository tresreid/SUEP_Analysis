#define isrTagger_cxx

std::tuple<Jet, Jet, Jet, Jet> getLeadingJets(std::vector<Jet> jets) {
  Jet lead_pt_jet;
  Jet second_pt_jet;

  Jet lead_y_jet;
  Jet second_y_jet;
  Jet lead_y_p_jet;
  Jet second_y_p_jet;
  Jet lead_y_m_jet;
  Jet second_y_m_jet;

  float lead_pt = 0;
  float lead_y_p = 0;
  float lead_y_m = 0;
  for (auto jet : jets) {
    if (jet.p4.Pt() > lead_pt) {
      lead_pt = jet.p4.Pt();
      second_pt_jet = lead_pt_jet;
      lead_pt_jet = jet;
    }
    if (jet.p4.Rapidity() > 0){
      if (jet.p4.Rapidity() > lead_y_p) {
        lead_y_p = jet.p4.Rapidity();
        second_y_p_jet = lead_y_p_jet;
        lead_y_p_jet = jet;
      }
    }else {
      if (jet.p4.Rapidity() < lead_y_m) {
        lead_y_m = jet.p4.Rapidity();
        second_y_m_jet = lead_y_m_jet;
        lead_y_m_jet = jet;
      }
    }
  }

  if (abs(lead_y_p_jet.p4.Rapidity()-second_y_p_jet.p4.Rapidity()) > abs(lead_y_m_jet.p4.Rapidity()-second_y_m_jet.p4.Rapidity())) {
    lead_y_jet = lead_y_p_jet;
    second_y_jet = second_y_p_jet;
  }else {
    lead_y_jet = lead_y_m_jet;
    second_y_jet = second_y_m_jet;
  }
  return std::make_tuple(lead_pt_jet,second_pt_jet,lead_y_jet,second_y_jet);
}

Jet isrTagger(std::vector<Jet> jets) {
  std::tuple<Jet, Jet, Jet, Jet> leadingJets = getLeadingJets(std::vector<Jet> jets);

  float r_pt = std::get<0>(leadingJets).p4.Pt()/std::get<1>(leadingJets).p4.Pt();
  float dy = abs(std::get<2>(leadingJets).p4.Rapidity()-std::get<3>(leadingJets).p4.Rapidity());
  plotter.Plot1D(Form("%s_offline_r_pt",s_sample.c_str()),";r_pt", r_pt, 50,0,3 );
  plotter.Plot1D(Form("%s_offline_dy",s_sample.c_str()),";dy", dy, 50,0,4 );

  bool isPtLeading = false;
  bool isYLeading = false;

  if (r_pt > 2) {
    if (abs(std::get<0>(leadingJets).p4.Rapidity()) > 1 &&
        abs(std::get<0>(leadingJets).p4.Rapidity()-std::get<1>(leadingJets).p4.Rapidity()) > 0.5) {
      isPtLeading = true;
    }
  }
  if (dy > 1.5){
    if (abs(std::get<2>(leadingJets).p4.Rapidity()) > 1 &&
        abs(std::get<2>(leadingJets).p4.Rapidity()-std::get<3>(leadingJets).p4.Rapidity()) > 0.5) {
      isYLeading = true;
    }
  }

  if (isPtLeading && isYLeading){
    std::cout << "Both jets look promising!" << "\n";
    return 
  }else if (isPtLeading && !isYLeading){
    return std::get<0>(leadingJets);
  }else if (!isPtLeading && isYLeading) {
    return std::get<2>(leadingJets);
  }else {
    std::cout << "No jet looks promissing!" << "\n";
  }
  Jet nullJet;
  return nullJet;
}

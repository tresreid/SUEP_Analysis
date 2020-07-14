#define isrTagger_cxx

std::pair<Jet, Jet> isrTagger(std::vector<Jet> jets) {
  Jet isr_jet;
  Jet lead_pt_jet;
  Jet second_pt_jet;

  Jet lead_y_jet;
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
    if (jet.p4.Rapidity()>0){
      if (jet.p4.Rapidity()>lead_y_p) {
        lead_y_p = jet.p4.Rapidity();
        second_y_p_jet = lead_y_p_jet;
        lead_y_p_jet = jet;
      }
    }else {
      if (jet.p4.Rapidity()<lead_y_m) {
        lead_y_m = jet.p4.Rapidity();
        second_y_m_jet = lead_y_m_jet;
        lead_y_m_jet = jet;
      }
    }
  }

  if (abs(lead_y_p_jet.p4.Rapidity()-second_y_p_jet).p4.Rapidity()>abs(lead_y_m_jet.p4.Rapidity()-second_y_m_jet.p4.Rapidity())) {
    lead_y_jet = lead_y_p_jet;
  } else {
    lead_y_jet = lead_y_m_jet;
  }
  return std::make_pair(lead_pt_jet,lead_y_jet);
}

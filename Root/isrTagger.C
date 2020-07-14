#define isrTagger_cxx

Jet isrTagger(std::vector<Jet> jets) {
  Jet isr_jet;
  Jet lead_pt_jet;
  Jet second_pt_jet;

  float lead_pt = 0;
  float lead_y = 0;

  for (auto jet : jets) {
    if (jet.p4.Pt() > lead_pt) {
      lead_pt = jet.p4.Pt();
      second_pt_jet = lead_pt_jet;
      lead_pt_jet = jet;
      //std::cout << lead_pt << "\n";
    }
  }
  return lead_pt_jet;
  //return isr_jet;
}

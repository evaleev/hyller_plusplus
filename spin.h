
#ifndef _hyller_spin_h_
#define _hyller_spin_h_

namespace hyller {

  typedef enum {SpinAlpha = 1, SpinBeta = -1} SpinCase1;
  typedef enum {SpinAlphaAlpha = 1, SpinAlphaBeta = 0} SpinCase2;
  inline SpinCase1 spin1(SpinCase2 sc2) { return SpinAlpha; }
  inline SpinCase1 spin2(SpinCase2 sc2) { return sc2==SpinAlphaBeta ? SpinBeta : SpinAlpha; }
  inline bool samespin(SpinCase2 sc2) { return sc2 != SpinAlphaBeta; }
  typedef enum {SpinSinglet = 1, SpinTriplet = 3} Spin2;

};

#endif

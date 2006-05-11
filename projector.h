
#ifndef _hyller_projector_h_
#define _hyller_projector_h_

#include "determinant.h"
#include "csf.h"
#include "hylleraas.h"

namespace hyller {

  /// Projector from Slater Determinant to Hylleraas basis
  class SD_2_Hylleraas {
  public:
    SD_2_Hylleraas(const SDBasisSet& sbs, const HylleraasBasisSet& hbs);
    ~SD_2_Hylleraas();

    /// Returns Slater Determinants (in rows) expanded in terms of Hylleraas functions
    double** X() const { return X_; }
  private:
    const SDBasisSet& sbs_;
    const HylleraasBasisSet& hbs_;
    double** X_;
  };

  /// Projector from CSF to Hylleraas basis
  class CSF_2_Hylleraas {
  public:
    CSF_2_Hylleraas(const CSFBasisSet& cbs, const HylleraasBasisSet& hbs);
    ~CSF_2_Hylleraas();

    /// Returns Hylleraas functions (with m=0; in rows) expanded in CSFs
    double** X() const { return X_; }
    /// Returns CSFs expaned in Hylleraas basis
    double** Xinv() const { return Xinv_; }

  private:
    const CSFBasisSet& cbs_;
    const HylleraasBasisSet& hbs_;
    double** X_;
    double** Xinv_;
  };

  /// Projector from SD to CSF
  class SD_2_CSF {
  public:
    SD_2_CSF(const SDBasisSet& sbs, const CSFBasisSet& cbs);
    ~SD_2_CSF();

    /// Returns CSFs expanded in Slater determinants
    double** X() const { return X_; }

  private:
    const SDBasisSet& sbs_;
    const CSFBasisSet& cbs_;
    double** X_;
  };

};

#endif

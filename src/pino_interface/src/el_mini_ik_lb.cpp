#include "pino_interface/ik_fast_comm.h"

class IKSolver_lb {
public:
  IkReal j0, cj0, sj0, htj0, j0mul, j1, cj1, sj1, htj1, j1mul, j2, cj2, sj2,
      htj2, j2mul, new_px, px, npx, new_py, py, npy, new_pz, pz, npz, pp;
  unsigned char _ij0[2], _nj0, _ij1[2], _nj1, _ij2[2], _nj2;

  IkReal j100, cj100, sj100;
  unsigned char _ij100[2], _nj100;
  bool ComputeIk(const IkReal *eetrans, const IkReal *eerot,
                 const IkReal *pfree, IkSolutionListBase<IkReal> &solutions) {
    // printf("IKSolver_lb");

    j0 = numeric_limits<IkReal>::quiet_NaN();
    _ij0[0] = -1;
    _ij0[1] = -1;
    _nj0 = -1;
    j1 = numeric_limits<IkReal>::quiet_NaN();
    _ij1[0] = -1;
    _ij1[1] = -1;
    _nj1 = -1;
    j2 = numeric_limits<IkReal>::quiet_NaN();
    _ij2[0] = -1;
    _ij2[1] = -1;
    _nj2 = -1;
    for (int dummyiter = 0; dummyiter < 1; ++dummyiter) {
      solutions.Clear();
      px = eetrans[0];
      py = eetrans[1];
      pz = eetrans[2];

      new_px = px;
      new_py = ((-1.0) * py);
      new_pz = ((-0.045056) + (((-1.0) * pz)));
      px = new_px;
      py = new_py;
      pz = new_pz;
      pp = ((px * px) + (py * py) + (pz * pz));
      {
        IkReal j0eval[2];
        j0eval[0] = ((px * px) + (py * py));
        j0eval[1] = ((IKabs(px)) + (IKabs(py)));
        if (IKabs(j0eval[0]) < 0.0000010000000000 ||
            IKabs(j0eval[1]) < 0.0000010000000000) {
          {
            IkReal evalcond[1];
            bool bgotonextstatement = true;
            do {
              evalcond[0] = ((px * px) + (py * py));
              if (IKabs(evalcond[0]) < 0.0000050000000000) {
                bgotonextstatement = false;
                {
                  IkReal j2array[2], cj2array[2], sj2array[2];
                  bool j2valid[2] = {false};
                  _nj2 = 2;
                  if ((((0.84697032327434) +
                        (((-16.5782493363067) * (pz * pz))))) <
                          -1 - IKFAST_SINCOS_THRESH ||
                      (((0.84697032327434) +
                        (((-16.5782493363067) * (pz * pz))))) >
                          1 + IKFAST_SINCOS_THRESH)
                    continue;
                  IkReal x15 = IKasin(((0.84697032327434) +
                                       (((-16.5782493363067) * (pz * pz)))));
                  j2array[0] = ((1.57080004761718) + (((-1.0) * x15)));
                  sj2array[0] = IKsin(j2array[0]);
                  cj2array[0] = IKcos(j2array[0]);
                  j2array[1] = ((4.71239270120697) + x15);
                  sj2array[1] = IKsin(j2array[1]);
                  cj2array[1] = IKcos(j2array[1]);
                  if (j2array[0] > IKPI) {
                    j2array[0] -= IK2PI;
                  } else if (j2array[0] < -IKPI) {
                    j2array[0] += IK2PI;
                  }
                  j2valid[0] = true;
                  if (j2array[1] > IKPI) {
                    j2array[1] -= IK2PI;
                  } else if (j2array[1] < -IKPI) {
                    j2array[1] += IK2PI;
                  }
                  j2valid[1] = true;
                  for (int ij2 = 0; ij2 < 2; ++ij2) {
                    if (!j2valid[ij2]) {
                      continue;
                    }
                    _ij2[0] = ij2;
                    _ij2[1] = -1;
                    for (int iij2 = ij2 + 1; iij2 < 2; ++iij2) {
                      if (j2valid[iij2] &&
                          IKabs(cj2array[ij2] - cj2array[iij2]) <
                              IKFAST_SOLUTION_THRESH &&
                          IKabs(sj2array[ij2] - sj2array[iij2]) <
                              IKFAST_SOLUTION_THRESH) {
                        j2valid[iij2] = false;
                        _ij2[1] = iij2;
                        break;
                      }
                    }
                    j2 = j2array[ij2];
                    cj2 = cj2array[ij2];
                    sj2 = sj2array[ij2];

                    {
                      IkReal j1eval[3];
                      px = 0;
                      py = 0;
                      pp = pz * pz;
                      IkReal x16 = pz * pz;
                      IkReal x17 = (pz * sj2);
                      IkReal x18 = (cj2 * pz);
                      j1eval[0] = ((1.0) + (((44.4503709630156) * x16)));
                      j1eval[1] =
                          ((IKabs(((1439.904) + (((-1274.915) * cj2)) +
                                   (((-8500.0) * x17)) +
                                   (((2320000000.0) * x18)) +
                                   (((-347976800.0) * sj2)) +
                                   (((-1300000000.0) * pz))))) +
                           (IKabs(((-194987000.0) + (((-1274.915) * sj2)) +
                                   (((2320000000.0) * x17)) +
                                   (((-9600.0) * pz)) + (((8500.0) * x18)) +
                                   (((347976800.0) * cj2))))));
                      j1eval[2] =
                          IKsign(((224970001.0) + (((10000000000.0) * x16))));
                      if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                          IKabs(j1eval[1]) < 0.0000010000000000 ||
                          IKabs(j1eval[2]) < 0.0000010000000000) {
                        {
                          IkReal j1eval[3];
                          px = 0;
                          py = 0;
                          pp = pz * pz;
                          IkReal x19 = pz * pz;
                          IkReal x20 = (cj2 * pz);
                          IkReal x21 = (pz * sj2);
                          j1eval[0] = ((1.0) + (((44.4503709630156) * x19)));
                          j1eval[1] =
                              IKsign(((0.0043194240192) + (((0.192) * x19))));
                          j1eval[2] =
                              ((IKabs((
                                   (-766.287659247114) + (((-0.02496) * pz)) +
                                   (((-0.003314779) * sj2)) +
                                   (((904.73968) * cj2)) + (((14999.0) * x19)) +
                                   (((-1.632e-7) * x21)) +
                                   (((0.044544) * x20))))) +
                               (IKabs(((-0.0037437504) +
                                       (((5108.92498998009) * pz)) +
                                       (((0.0221) * x21)) +
                                       (((-100000.0) * (pz * pz * pz))) +
                                       (((0.00668115456) * cj2)) +
                                       (((-2.4478368e-8) * sj2)) +
                                       (((-6032.0) * x20))))));
                          if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                              IKabs(j1eval[1]) < 0.0000010000000000 ||
                              IKabs(j1eval[2]) < 0.0000010000000000) {
                            {
                              IkReal j1eval[1];
                              px = 0;
                              py = 0;
                              pp = pz * pz;
                              j1eval[0] = ((1.12941176470588) +
                                           (((-272941.176470588) * sj2)) +
                                           (((6.66711114074272) * pz * sj2)) +
                                           (((1019675.82152536) * pz)) +
                                           (((-1.0) * cj2)) +
                                           (((-1819729.15841448) * cj2 * pz)));
                              if (IKabs(j1eval[0]) < 0.0000010000000000) {
                                continue; // no branches [j0, j1]

                              } else {
                                {
                                  IkReal j1array[1], cj1array[1], sj1array[1];
                                  bool j1valid[1] = {false};
                                  _nj1 = 1;
                                  IkReal x22 = cj2 * cj2;
                                  IkReal x23 = (cj2 * pz);
                                  IkReal x24 = (pz * sj2);
                                  IkReal x25 = (cj2 * sj2);
                                  CheckValue<IkReal> x26 =
                                      IKPowWithIntegerCheck(
                                          ((0.01439904) +
                                           (((-3479.768) * sj2)) +
                                           (((13000.0) * pz)) +
                                           (((-23200.0) * x23)) +
                                           (((-0.01274915) * cj2)) +
                                           (((0.085) * x24))),
                                          -1);
                                  if (!x26.valid) {
                                    continue;
                                  }
                                  CheckValue<IkReal> x27 =
                                      IKPowWithIntegerCheck(
                                          ((1439.904) + (((-1274.915) * cj2)) +
                                           (((-347976800.0) * sj2)) +
                                           (((8500.0) * x24)) +
                                           (((1300000000.0) * pz)) +
                                           (((-2320000000.0) * x23))),
                                          -1);
                                  if (!x27.valid) {
                                    continue;
                                  }
                                  if (IKabs(
                                          ((x26.value) *
                                           (((0.00724) + (((-0.03944) * x22)) +
                                             (((-5382.39999992775) * x25)) +
                                             (((3015.9999999184) * sj2)) +
                                             (((-14999.0) * pz)) +
                                             (((0.033322) * cj2)))))) <
                                          IKFAST_ATAN2_MAGTHRESH &&
                                      IKabs(((x27.value) *
                                             (((55970000.992775) +
                                               (((-2210.0) * sj2)) +
                                               (((603200000.0) * cj2)) +
                                               (((-538239999.992775) * x22)) +
                                               (((3944.0) * x25)))))) <
                                          IKFAST_ATAN2_MAGTHRESH &&
                                      IKabs(
                                          IKsqr(
                                              ((x26.value) *
                                               (((0.00724) +
                                                 (((-0.03944) * x22)) +
                                                 (((-5382.39999992775) * x25)) +
                                                 (((3015.9999999184) * sj2)) +
                                                 (((-14999.0) * pz)) +
                                                 (((0.033322) * cj2)))))) +
                                          IKsqr(
                                              ((x27.value) *
                                               (((55970000.992775) +
                                                 (((-2210.0) * sj2)) +
                                                 (((603200000.0) * cj2)) +
                                                 (((-538239999.992775) * x22)) +
                                                 (((3944.0) * x25)))))) -
                                          1) <= IKFAST_SINCOS_THRESH)
                                    continue;
                                  j1array[0] = IKatan2(
                                      ((x26.value) *
                                       (((0.00724) + (((-0.03944) * x22)) +
                                         (((-5382.39999992775) * x25)) +
                                         (((3015.9999999184) * sj2)) +
                                         (((-14999.0) * pz)) +
                                         (((0.033322) * cj2))))),
                                      ((x27.value) *
                                       (((55970000.992775) +
                                         (((-2210.0) * sj2)) +
                                         (((603200000.0) * cj2)) +
                                         (((-538239999.992775) * x22)) +
                                         (((3944.0) * x25))))));
                                  sj1array[0] = IKsin(j1array[0]);
                                  cj1array[0] = IKcos(j1array[0]);
                                  if (j1array[0] > IKPI) {
                                    j1array[0] -= IK2PI;
                                  } else if (j1array[0] < -IKPI) {
                                    j1array[0] += IK2PI;
                                  }
                                  j1valid[0] = true;
                                  for (int ij1 = 0; ij1 < 1; ++ij1) {
                                    if (!j1valid[ij1]) {
                                      continue;
                                    }
                                    _ij1[0] = ij1;
                                    _ij1[1] = -1;
                                    for (int iij1 = ij1 + 1; iij1 < 1; ++iij1) {
                                      if (j1valid[iij1] &&
                                          IKabs(cj1array[ij1] -
                                                cj1array[iij1]) <
                                              IKFAST_SOLUTION_THRESH &&
                                          IKabs(sj1array[ij1] -
                                                sj1array[iij1]) <
                                              IKFAST_SOLUTION_THRESH) {
                                        j1valid[iij1] = false;
                                        _ij1[1] = iij1;
                                        break;
                                      }
                                    }
                                    j1 = j1array[ij1];
                                    cj1 = cj1array[ij1];
                                    sj1 = sj1array[ij1];
                                    {
                                      IkReal evalcond[5];
                                      IkReal x28 = IKcos(j1);
                                      IkReal x29 = IKsin(j1);
                                      IkReal x30 = ((0.232) * cj2);
                                      IkReal x31 = ((1.0) * pz);
                                      IkReal x32 = ((0.232) * sj2);
                                      IkReal x33 = ((8.5e-7) * cj2);
                                      IkReal x34 = ((8.5e-7) * sj2);
                                      IkReal x35 = ((8.5e-7) * x28);
                                      evalcond[0] =
                                          ((-0.13) + (((-1.0) * x28 * x31)) +
                                           (((-1.0) * x34)) +
                                           (((-0.14999) * x29)) + x30);
                                      evalcond[1] =
                                          ((-9.6e-7) + (((-1.0) * x29 * x31)) +
                                           x33 + x32 + (((0.14999) * x28)));
                                      evalcond[2] = ((0.0172892498998009) +
                                                     (((-1.0) * pz * x31)) +
                                                     (((-0.26) * pz * x28)) +
                                                     (((2.879808e-7) * x28)) +
                                                     (((-1.92e-6) * pz * x29)) +
                                                     (((-0.0389974) * x29)));
                                      evalcond[3] =
                                          ((-0.14999) + (((-1.0) * x29 * x34)) +
                                           ((x29 * x30)) +
                                           (((-1.0) * x28 * x32)) +
                                           (((-1.0) * x28 * x33)) +
                                           (((-0.13) * x29)) +
                                           (((9.6e-7) * x28)));
                                      evalcond[4] =
                                          ((((-9.6e-7) * x29)) + ((x29 * x33)) +
                                           ((x29 * x32)) +
                                           (((-1.0) * x28 * x34)) +
                                           ((x28 * x30)) + (((-1.0) * x31)) +
                                           (((-0.13) * x28)));
                                      if (IKabs(evalcond[0]) >
                                              IKFAST_EVALCOND_THRESH ||
                                          IKabs(evalcond[1]) >
                                              IKFAST_EVALCOND_THRESH ||
                                          IKabs(evalcond[2]) >
                                              IKFAST_EVALCOND_THRESH ||
                                          IKabs(evalcond[3]) >
                                              IKFAST_EVALCOND_THRESH ||
                                          IKabs(evalcond[4]) >
                                              IKFAST_EVALCOND_THRESH) {
                                        continue;
                                      }
                                    }

                                    {
                                      IkReal j0array[1], cj0array[1],
                                          sj0array[1];
                                      bool j0valid[1] = {false};
                                      _nj0 = 1;
                                      j0array[0] = 0;
                                      sj0array[0] = IKsin(j0array[0]);
                                      cj0array[0] = IKcos(j0array[0]);
                                      if (j0array[0] > IKPI) {
                                        j0array[0] -= IK2PI;
                                      } else if (j0array[0] < -IKPI) {
                                        j0array[0] += IK2PI;
                                      }
                                      j0valid[0] = true;
                                      for (int ij0 = 0; ij0 < 1; ++ij0) {
                                        if (!j0valid[ij0]) {
                                          continue;
                                        }
                                        _ij0[0] = ij0;
                                        _ij0[1] = -1;
                                        for (int iij0 = ij0 + 1; iij0 < 1;
                                             ++iij0) {
                                          if (j0valid[iij0] &&
                                              IKabs(cj0array[ij0] -
                                                    cj0array[iij0]) <
                                                  IKFAST_SOLUTION_THRESH &&
                                              IKabs(sj0array[ij0] -
                                                    sj0array[iij0]) <
                                                  IKFAST_SOLUTION_THRESH) {
                                            j0valid[iij0] = false;
                                            _ij0[1] = iij0;
                                            break;
                                          }
                                        }
                                        j0 = j0array[ij0];
                                        cj0 = cj0array[ij0];
                                        sj0 = sj0array[ij0];

                                        {
                                          std::vector<
                                              IkSingleDOFSolutionBase<IkReal>>
                                              vinfos(3);
                                          vinfos[0].jointtype = 1;
                                          vinfos[0].foffset = j0;
                                          vinfos[0].indices[0] = _ij0[0];
                                          vinfos[0].indices[1] = _ij0[1];
                                          vinfos[0].maxsolutions = _nj0;
                                          vinfos[1].jointtype = 1;
                                          vinfos[1].foffset = j1;
                                          vinfos[1].indices[0] = _ij1[0];
                                          vinfos[1].indices[1] = _ij1[1];
                                          vinfos[1].maxsolutions = _nj1;
                                          vinfos[2].jointtype = 1;
                                          vinfos[2].foffset = j2;
                                          vinfos[2].indices[0] = _ij2[0];
                                          vinfos[2].indices[1] = _ij2[1];
                                          vinfos[2].maxsolutions = _nj2;
                                          std::vector<int> vfree(0);
                                          solutions.AddSolution(vinfos, vfree);
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }

                          } else {
                            {
                              IkReal j1array[1], cj1array[1], sj1array[1];
                              bool j1valid[1] = {false};
                              _nj1 = 1;
                              IkReal x36 = pz * pz;
                              IkReal x37 = (cj2 * pz);
                              IkReal x38 = (pz * sj2);
                              CheckValue<IkReal> x39 = IKPowWithIntegerCheck(
                                  IKsign(
                                      ((0.0043194240192) + (((0.192) * x36)))),
                                  -1);
                              if (!x39.valid) {
                                continue;
                              }
                              CheckValue<IkReal> x40 = IKatan2WithCheck(
                                  IkReal(((-0.0037437504) +
                                          (((5108.92498998009) * pz)) +
                                          (((-6032.0) * x37)) +
                                          (((-100000.0) * (pz * pz * pz))) +
                                          (((0.0221) * x38)) +
                                          (((0.00668115456) * cj2)) +
                                          (((-2.4478368e-8) * sj2)))),
                                  IkReal(((-766.287659247114) +
                                          (((14999.0) * x36)) +
                                          (((-0.02496) * pz)) +
                                          (((-0.003314779) * sj2)) +
                                          (((904.73968) * cj2)) +
                                          (((0.044544) * x37)) +
                                          (((-1.632e-7) * x38)))),
                                  IKFAST_ATAN2_MAGTHRESH);
                              if (!x40.valid) {
                                continue;
                              }
                              j1array[0] =
                                  ((-1.5707963267949) +
                                   (((1.5707963267949) * (x39.value))) +
                                   (x40.value));
                              sj1array[0] = IKsin(j1array[0]);
                              cj1array[0] = IKcos(j1array[0]);
                              if (j1array[0] > IKPI) {
                                j1array[0] -= IK2PI;
                              } else if (j1array[0] < -IKPI) {
                                j1array[0] += IK2PI;
                              }
                              j1valid[0] = true;
                              for (int ij1 = 0; ij1 < 1; ++ij1) {
                                if (!j1valid[ij1]) {
                                  continue;
                                }
                                _ij1[0] = ij1;
                                _ij1[1] = -1;
                                for (int iij1 = ij1 + 1; iij1 < 1; ++iij1) {
                                  if (j1valid[iij1] &&
                                      IKabs(cj1array[ij1] - cj1array[iij1]) <
                                          IKFAST_SOLUTION_THRESH &&
                                      IKabs(sj1array[ij1] - sj1array[iij1]) <
                                          IKFAST_SOLUTION_THRESH) {
                                    j1valid[iij1] = false;
                                    _ij1[1] = iij1;
                                    break;
                                  }
                                }
                                j1 = j1array[ij1];
                                cj1 = cj1array[ij1];
                                sj1 = sj1array[ij1];
                                {
                                  IkReal evalcond[5];
                                  IkReal x41 = IKcos(j1);
                                  IkReal x42 = IKsin(j1);
                                  IkReal x43 = ((0.232) * cj2);
                                  IkReal x44 = ((1.0) * pz);
                                  IkReal x45 = ((0.232) * sj2);
                                  IkReal x46 = ((8.5e-7) * cj2);
                                  IkReal x47 = ((8.5e-7) * sj2);
                                  IkReal x48 = ((8.5e-7) * x41);
                                  evalcond[0] = ((-0.13) + (((-1.0) * x47)) +
                                                 (((-1.0) * x41 * x44)) +
                                                 (((-0.14999) * x42)) + x43);
                                  evalcond[1] =
                                      ((-9.6e-7) + (((-1.0) * x42 * x44)) +
                                       x46 + x45 + (((0.14999) * x41)));
                                  evalcond[2] = ((0.0172892498998009) +
                                                 (((-1.0) * pz * x44)) +
                                                 (((2.879808e-7) * x41)) +
                                                 (((-0.26) * pz * x41)) +
                                                 (((-0.0389974) * x42)) +
                                                 (((-1.92e-6) * pz * x42)));
                                  evalcond[3] =
                                      ((-0.14999) + (((-1.0) * x42 * x47)) +
                                       (((-1.0) * x41 * x46)) +
                                       (((-1.0) * x41 * x45)) +
                                       (((9.6e-7) * x41)) + ((x42 * x43)) +
                                       (((-0.13) * x42)));
                                  evalcond[4] =
                                      (((x41 * x43)) + (((-1.0) * x44)) +
                                       (((-1.0) * x41 * x47)) +
                                       (((-9.6e-7) * x42)) + ((x42 * x45)) +
                                       ((x42 * x46)) + (((-0.13) * x41)));
                                  if (IKabs(evalcond[0]) >
                                          IKFAST_EVALCOND_THRESH ||
                                      IKabs(evalcond[1]) >
                                          IKFAST_EVALCOND_THRESH ||
                                      IKabs(evalcond[2]) >
                                          IKFAST_EVALCOND_THRESH ||
                                      IKabs(evalcond[3]) >
                                          IKFAST_EVALCOND_THRESH ||
                                      IKabs(evalcond[4]) >
                                          IKFAST_EVALCOND_THRESH) {
                                    continue;
                                  }
                                }

                                {
                                  IkReal j0array[1], cj0array[1], sj0array[1];
                                  bool j0valid[1] = {false};
                                  _nj0 = 1;
                                  j0array[0] = 0;
                                  sj0array[0] = IKsin(j0array[0]);
                                  cj0array[0] = IKcos(j0array[0]);
                                  if (j0array[0] > IKPI) {
                                    j0array[0] -= IK2PI;
                                  } else if (j0array[0] < -IKPI) {
                                    j0array[0] += IK2PI;
                                  }
                                  j0valid[0] = true;
                                  for (int ij0 = 0; ij0 < 1; ++ij0) {
                                    if (!j0valid[ij0]) {
                                      continue;
                                    }
                                    _ij0[0] = ij0;
                                    _ij0[1] = -1;
                                    for (int iij0 = ij0 + 1; iij0 < 1; ++iij0) {
                                      if (j0valid[iij0] &&
                                          IKabs(cj0array[ij0] -
                                                cj0array[iij0]) <
                                              IKFAST_SOLUTION_THRESH &&
                                          IKabs(sj0array[ij0] -
                                                sj0array[iij0]) <
                                              IKFAST_SOLUTION_THRESH) {
                                        j0valid[iij0] = false;
                                        _ij0[1] = iij0;
                                        break;
                                      }
                                    }
                                    j0 = j0array[ij0];
                                    cj0 = cj0array[ij0];
                                    sj0 = sj0array[ij0];

                                    {
                                      std::vector<
                                          IkSingleDOFSolutionBase<IkReal>>
                                          vinfos(3);
                                      vinfos[0].jointtype = 1;
                                      vinfos[0].foffset = j0;
                                      vinfos[0].indices[0] = _ij0[0];
                                      vinfos[0].indices[1] = _ij0[1];
                                      vinfos[0].maxsolutions = _nj0;
                                      vinfos[1].jointtype = 1;
                                      vinfos[1].foffset = j1;
                                      vinfos[1].indices[0] = _ij1[0];
                                      vinfos[1].indices[1] = _ij1[1];
                                      vinfos[1].maxsolutions = _nj1;
                                      vinfos[2].jointtype = 1;
                                      vinfos[2].foffset = j2;
                                      vinfos[2].indices[0] = _ij2[0];
                                      vinfos[2].indices[1] = _ij2[1];
                                      vinfos[2].maxsolutions = _nj2;
                                      std::vector<int> vfree(0);
                                      solutions.AddSolution(vinfos, vfree);
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }

                      } else {
                        {
                          IkReal j1array[1], cj1array[1], sj1array[1];
                          bool j1valid[1] = {false};
                          _nj1 = 1;
                          IkReal x49 = (pz * sj2);
                          IkReal x50 = (cj2 * pz);
                          CheckValue<IkReal> x51 = IKPowWithIntegerCheck(
                              IKsign(((224970001.0) +
                                      (((10000000000.0) * (pz * pz))))),
                              -1);
                          if (!x51.valid) {
                            continue;
                          }
                          CheckValue<IkReal> x52 = IKatan2WithCheck(
                              IkReal(((-194987000.0) + (((-1274.915) * sj2)) +
                                      (((2320000000.0) * x49)) +
                                      (((-9600.0) * pz)) + (((8500.0) * x50)) +
                                      (((347976800.0) * cj2)))),
                              IkReal(((1439.904) + (((-8500.0) * x49)) +
                                      (((-1274.915) * cj2)) +
                                      (((-347976800.0) * sj2)) +
                                      (((-1300000000.0) * pz)) +
                                      (((2320000000.0) * x50)))),
                              IKFAST_ATAN2_MAGTHRESH);
                          if (!x52.valid) {
                            continue;
                          }
                          j1array[0] = ((-1.5707963267949) +
                                        (((1.5707963267949) * (x51.value))) +
                                        (x52.value));
                          sj1array[0] = IKsin(j1array[0]);
                          cj1array[0] = IKcos(j1array[0]);
                          if (j1array[0] > IKPI) {
                            j1array[0] -= IK2PI;
                          } else if (j1array[0] < -IKPI) {
                            j1array[0] += IK2PI;
                          }
                          j1valid[0] = true;
                          for (int ij1 = 0; ij1 < 1; ++ij1) {
                            if (!j1valid[ij1]) {
                              continue;
                            }
                            _ij1[0] = ij1;
                            _ij1[1] = -1;
                            for (int iij1 = ij1 + 1; iij1 < 1; ++iij1) {
                              if (j1valid[iij1] &&
                                  IKabs(cj1array[ij1] - cj1array[iij1]) <
                                      IKFAST_SOLUTION_THRESH &&
                                  IKabs(sj1array[ij1] - sj1array[iij1]) <
                                      IKFAST_SOLUTION_THRESH) {
                                j1valid[iij1] = false;
                                _ij1[1] = iij1;
                                break;
                              }
                            }
                            j1 = j1array[ij1];
                            cj1 = cj1array[ij1];
                            sj1 = sj1array[ij1];
                            {
                              IkReal evalcond[5];
                              IkReal x53 = IKcos(j1);
                              IkReal x54 = IKsin(j1);
                              IkReal x55 = ((0.232) * cj2);
                              IkReal x56 = ((1.0) * pz);
                              IkReal x57 = ((0.232) * sj2);
                              IkReal x58 = ((8.5e-7) * cj2);
                              IkReal x59 = ((8.5e-7) * sj2);
                              IkReal x60 = ((8.5e-7) * x53);
                              evalcond[0] =
                                  ((-0.13) + (((-1.0) * x53 * x56)) + x55 +
                                   (((-0.14999) * x54)) + (((-1.0) * x59)));
                              evalcond[1] =
                                  ((-9.6e-7) + (((-1.0) * x54 * x56)) +
                                   (((0.14999) * x53)) + x58 + x57);
                              evalcond[2] = ((0.0172892498998009) +
                                             (((-1.92e-6) * pz * x54)) +
                                             (((-0.0389974) * x54)) +
                                             (((-0.26) * pz * x53)) +
                                             (((-1.0) * pz * x56)) +
                                             (((2.879808e-7) * x53)));
                              evalcond[3] =
                                  ((-0.14999) + (((-1.0) * x53 * x58)) +
                                   (((-1.0) * x53 * x57)) +
                                   (((-1.0) * x54 * x59)) + (((9.6e-7) * x53)) +
                                   (((-0.13) * x54)) + ((x54 * x55)));
                              evalcond[4] =
                                  ((((-1.0) * x53 * x59)) + ((x53 * x55)) +
                                   (((-9.6e-7) * x54)) + (((-1.0) * x56)) +
                                   (((-0.13) * x53)) + ((x54 * x57)) +
                                   ((x54 * x58)));
                              if (IKabs(evalcond[0]) > IKFAST_EVALCOND_THRESH ||
                                  IKabs(evalcond[1]) > IKFAST_EVALCOND_THRESH ||
                                  IKabs(evalcond[2]) > IKFAST_EVALCOND_THRESH ||
                                  IKabs(evalcond[3]) > IKFAST_EVALCOND_THRESH ||
                                  IKabs(evalcond[4]) > IKFAST_EVALCOND_THRESH) {
                                continue;
                              }
                            }

                            {
                              IkReal j0array[1], cj0array[1], sj0array[1];
                              bool j0valid[1] = {false};
                              _nj0 = 1;
                              j0array[0] = 0;
                              sj0array[0] = IKsin(j0array[0]);
                              cj0array[0] = IKcos(j0array[0]);
                              if (j0array[0] > IKPI) {
                                j0array[0] -= IK2PI;
                              } else if (j0array[0] < -IKPI) {
                                j0array[0] += IK2PI;
                              }
                              j0valid[0] = true;
                              for (int ij0 = 0; ij0 < 1; ++ij0) {
                                if (!j0valid[ij0]) {
                                  continue;
                                }
                                _ij0[0] = ij0;
                                _ij0[1] = -1;
                                for (int iij0 = ij0 + 1; iij0 < 1; ++iij0) {
                                  if (j0valid[iij0] &&
                                      IKabs(cj0array[ij0] - cj0array[iij0]) <
                                          IKFAST_SOLUTION_THRESH &&
                                      IKabs(sj0array[ij0] - sj0array[iij0]) <
                                          IKFAST_SOLUTION_THRESH) {
                                    j0valid[iij0] = false;
                                    _ij0[1] = iij0;
                                    break;
                                  }
                                }
                                j0 = j0array[ij0];
                                cj0 = cj0array[ij0];
                                sj0 = sj0array[ij0];

                                {
                                  std::vector<IkSingleDOFSolutionBase<IkReal>>
                                      vinfos(3);
                                  vinfos[0].jointtype = 1;
                                  vinfos[0].foffset = j0;
                                  vinfos[0].indices[0] = _ij0[0];
                                  vinfos[0].indices[1] = _ij0[1];
                                  vinfos[0].maxsolutions = _nj0;
                                  vinfos[1].jointtype = 1;
                                  vinfos[1].foffset = j1;
                                  vinfos[1].indices[0] = _ij1[0];
                                  vinfos[1].indices[1] = _ij1[1];
                                  vinfos[1].maxsolutions = _nj1;
                                  vinfos[2].jointtype = 1;
                                  vinfos[2].foffset = j2;
                                  vinfos[2].indices[0] = _ij2[0];
                                  vinfos[2].indices[1] = _ij2[1];
                                  vinfos[2].maxsolutions = _nj2;
                                  std::vector<int> vfree(0);
                                  solutions.AddSolution(vinfos, vfree);
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            } while (0);
            if (bgotonextstatement) {
              bool bgotonextstatement = true;
              do {
                evalcond[0] = ((IKabs(px)) + (IKabs(py)));
                if (IKabs(evalcond[0]) < 0.0000050000000000) {
                  bgotonextstatement = false;
                  {
                    IkReal j2array[2], cj2array[2], sj2array[2];
                    bool j2valid[2] = {false};
                    _nj2 = 2;
                    if ((((0.84697032327434) +
                          (((-16.5782493363067) * (pz * pz))))) <
                            -1 - IKFAST_SINCOS_THRESH ||
                        (((0.84697032327434) +
                          (((-16.5782493363067) * (pz * pz))))) >
                            1 + IKFAST_SINCOS_THRESH)
                      continue;
                    IkReal x61 = IKasin(((0.84697032327434) +
                                         (((-16.5782493363067) * (pz * pz)))));
                    j2array[0] = ((1.57080004761718) + (((-1.0) * x61)));
                    sj2array[0] = IKsin(j2array[0]);
                    cj2array[0] = IKcos(j2array[0]);
                    j2array[1] = ((4.71239270120697) + x61);
                    sj2array[1] = IKsin(j2array[1]);
                    cj2array[1] = IKcos(j2array[1]);
                    if (j2array[0] > IKPI) {
                      j2array[0] -= IK2PI;
                    } else if (j2array[0] < -IKPI) {
                      j2array[0] += IK2PI;
                    }
                    j2valid[0] = true;
                    if (j2array[1] > IKPI) {
                      j2array[1] -= IK2PI;
                    } else if (j2array[1] < -IKPI) {
                      j2array[1] += IK2PI;
                    }
                    j2valid[1] = true;
                    for (int ij2 = 0; ij2 < 2; ++ij2) {
                      if (!j2valid[ij2]) {
                        continue;
                      }
                      _ij2[0] = ij2;
                      _ij2[1] = -1;
                      for (int iij2 = ij2 + 1; iij2 < 2; ++iij2) {
                        if (j2valid[iij2] &&
                            IKabs(cj2array[ij2] - cj2array[iij2]) <
                                IKFAST_SOLUTION_THRESH &&
                            IKabs(sj2array[ij2] - sj2array[iij2]) <
                                IKFAST_SOLUTION_THRESH) {
                          j2valid[iij2] = false;
                          _ij2[1] = iij2;
                          break;
                        }
                      }
                      j2 = j2array[ij2];
                      cj2 = cj2array[ij2];
                      sj2 = sj2array[ij2];

                      {
                        IkReal j1eval[3];
                        px = 0;
                        py = 0;
                        pp = pz * pz;
                        IkReal x62 = pz * pz;
                        IkReal x63 = (pz * sj2);
                        IkReal x64 = (cj2 * pz);
                        j1eval[0] = ((1.0) + (((44.4503709630156) * x62)));
                        j1eval[1] =
                            ((IKabs(((1439.904) + (((-1274.915) * cj2)) +
                                     (((-347976800.0) * sj2)) +
                                     (((-1300000000.0) * pz)) +
                                     (((-8500.0) * x63)) +
                                     (((2320000000.0) * x64))))) +
                             (IKabs(((-194987000.0) + (((-1274.915) * sj2)) +
                                     (((8500.0) * x64)) + (((-9600.0) * pz)) +
                                     (((2320000000.0) * x63)) +
                                     (((347976800.0) * cj2))))));
                        j1eval[2] =
                            IKsign(((224970001.0) + (((10000000000.0) * x62))));
                        if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                            IKabs(j1eval[1]) < 0.0000010000000000 ||
                            IKabs(j1eval[2]) < 0.0000010000000000) {
                          {
                            IkReal j1eval[3];
                            px = 0;
                            py = 0;
                            pp = pz * pz;
                            IkReal x65 = pz * pz;
                            IkReal x66 = (cj2 * pz);
                            IkReal x67 = (pz * sj2);
                            j1eval[0] = ((1.0) + (((44.4503709630156) * x65)));
                            j1eval[1] =
                                IKsign(((0.0043194240192) + (((0.192) * x65))));
                            j1eval[2] =
                                ((IKabs(((-766.287659247114) +
                                         (((0.044544) * x66)) +
                                         (((-0.02496) * pz)) +
                                         (((-0.003314779) * sj2)) +
                                         (((14999.0) * x65)) +
                                         (((-1.632e-7) * x67)) +
                                         (((904.73968) * cj2))))) +
                                 (IKabs(((-0.0037437504) +
                                         (((5108.92498998009) * pz)) +
                                         (((-100000.0) * (pz * pz * pz))) +
                                         (((0.00668115456) * cj2)) +
                                         (((-2.4478368e-8) * sj2)) +
                                         (((-6032.0) * x66)) +
                                         (((0.0221) * x67))))));
                            if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                                IKabs(j1eval[1]) < 0.0000010000000000 ||
                                IKabs(j1eval[2]) < 0.0000010000000000) {
                              {
                                IkReal j1eval[1];
                                px = 0;
                                py = 0;
                                pp = pz * pz;
                                j1eval[0] =
                                    ((1.12941176470588) +
                                     (((-272941.176470588) * sj2)) +
                                     (((6.66711114074272) * pz * sj2)) +
                                     (((1019675.82152536) * pz)) +
                                     (((-1.0) * cj2)) +
                                     (((-1819729.15841448) * cj2 * pz)));
                                if (IKabs(j1eval[0]) < 0.0000010000000000) {
                                  {
                                    IkReal evalcond[1];
                                    bool bgotonextstatement = true;
                                    do {
                                      IkReal x68 =
                                          ((-1.0) +
                                           (((-1819729.15841448) * pz)));
                                      IkReal x69 =
                                          ((-272941.176470588) +
                                           (((6.66711114074272) * pz)));
                                      IkReal x70 =
                                          ((1.12941176470588) +
                                           (((1019675.82152536) * pz)));
                                      IkReal x71 = ((x68 * x68) + (x69 * x69));
                                      if ((((74496885814.1488) +
                                            (((3311414210028.33) *
                                              (pz * pz))))) < -0.00001)
                                        continue;
                                      IkReal x72 = IKabs(IKsqrt((
                                          (74496885814.1488) +
                                          (((3311414210028.33) * (pz * pz))))));
                                      CheckValue<IkReal> x78 =
                                          IKPowWithIntegerCheck(x72, -1);
                                      if (!x78.valid) {
                                        continue;
                                      }
                                      IkReal x73 = x78.value;
                                      IkReal x79 = x71;
                                      if (IKabs(x79) == 0) {
                                        continue;
                                      }
                                      IkReal x74 = pow(x79, -0.5);
                                      if ((((1.0) + (((-1.0) * (x70 * x70) *
                                                      (x73 * x73))))) <
                                          -0.00001)
                                        continue;
                                      IkReal x75 = IKsqrt(
                                          ((1.0) + (((-1.0) * (x70 * x70) *
                                                     (x73 * x73)))));
                                      IkReal x76 = ((1.0) * x70 * x73 * x74);
                                      IkReal x77 = (x74 * x75);
                                      CheckValue<IkReal> x80 = IKatan2WithCheck(
                                          IkReal(x68), IkReal(x69),
                                          IKFAST_ATAN2_MAGTHRESH);
                                      if (!x80.valid) {
                                        continue;
                                      }
                                      if ((x71) < -0.00001)
                                        continue;
                                      CheckValue<IkReal> x81 =
                                          IKPowWithIntegerCheck(
                                              IKabs(IKsqrt(x71)), -1);
                                      if (!x81.valid) {
                                        continue;
                                      }
                                      if (((x70 * (x81.value))) <
                                              -1 - IKFAST_SINCOS_THRESH ||
                                          ((x70 * (x81.value))) >
                                              1 + IKFAST_SINCOS_THRESH)
                                        continue;
                                      IkReal gconst0 =
                                          ((((-1.0) * (x80.value))) +
                                           (((-1.0) *
                                             (IKasin((x70 * (x81.value)))))));
                                      IkReal gconst1 = ((((-1.0) * x68 * x77)) +
                                                        (((-1.0) * x69 * x76)));
                                      IkReal gconst2 = (((x69 * x77)) +
                                                        (((-1.0) * x68 * x76)));
                                      IkReal x82 = j2;
                                      CheckValue<IkReal> x88 = IKatan2WithCheck(
                                          IkReal(
                                              ((-1.0) +
                                               (((-1819729.15841448) * pz)))),
                                          IkReal(((-272941.176470588) +
                                                  (((6.66711114074272) * pz)))),
                                          IKFAST_ATAN2_MAGTHRESH);
                                      if (!x88.valid) {
                                        continue;
                                      }
                                      IkReal x83 = x88.value;
                                      IkReal x84 = x83;
                                      if ((((74496885814.1488) +
                                            (((3311414210028.33) *
                                              (pz * pz))))) < -0.00001)
                                        continue;
                                      CheckValue<IkReal> x89 =
                                          IKPowWithIntegerCheck(
                                              IKabs(
                                                  IKsqrt(((74496885814.1488) +
                                                          (((3311414210028.33) *
                                                            (pz * pz)))))),
                                              -1);
                                      if (!x89.valid) {
                                        continue;
                                      }
                                      IkReal x85 = x89.value;
                                      if ((((((1.12941176470588) * x85)) +
                                            (((1019675.82152536) * pz *
                                              x85)))) <
                                              -1 - IKFAST_SINCOS_THRESH ||
                                          (((((1.12941176470588) * x85)) +
                                            (((1019675.82152536) * pz *
                                              x85)))) >
                                              1 + IKFAST_SINCOS_THRESH)
                                        continue;
                                      IkReal x86 = IKasin(
                                          ((((1.12941176470588) * x85)) +
                                           (((1019675.82152536) * pz * x85))));
                                      IkReal x87 = x86;
                                      if (((((x84 * x86)) + ((x82 * x86)) +
                                            ((x82 * x83)) + ((x86 * x87)) +
                                            ((j2 * x82)) + ((j2 * x84)) +
                                            ((j2 * x87)) + ((x83 * x87)) +
                                            ((x83 * x84)))) < -0.00001)
                                        continue;
                                      evalcond[0] =
                                          ((-3.14159265358979) +
                                           (IKfmod(((3.14159265358979) +
                                                    (IKsqrt((((x84 * x86)) +
                                                             ((x82 * x86)) +
                                                             ((x82 * x83)) +
                                                             ((x86 * x87)) +
                                                             ((j2 * x82)) +
                                                             ((j2 * x84)) +
                                                             ((j2 * x87)) +
                                                             ((x83 * x87)) +
                                                             ((x83 * x84)))))),
                                                   6.28318530717959)));
                                      if (IKabs(evalcond[0]) <
                                          0.0000050000000000) {
                                        bgotonextstatement = false;
                                        {
                                          IkReal j1eval[2];
                                          IkReal x90 =
                                              ((6.66711114074272) * pz);
                                          IkReal x91 = x68;
                                          IkReal x92 = pz * pz;
                                          IkReal x93 =
                                              ((-272941.176470588) + x90);
                                          IkReal x94 =
                                              ((1.12941176470588) +
                                               (((1019675.82152536) * pz)));
                                          IkReal x95 =
                                              ((x93 * x93) + (x91 * x91));
                                          IkReal x96 = x72;
                                          CheckValue<IkReal> x104 =
                                              IKPowWithIntegerCheck(x96, -1);
                                          if (!x104.valid) {
                                            continue;
                                          }
                                          IkReal x97 = x104.value;
                                          IkReal x105 = x95;
                                          if (IKabs(x105) == 0) {
                                            continue;
                                          }
                                          IkReal x98 = pow(x105, -0.5);
                                          if ((x95) < -0.00001)
                                            continue;
                                          CheckValue<IkReal> x106 =
                                              IKPowWithIntegerCheck(
                                                  IKabs(IKsqrt(x95)), -1);
                                          if (!x106.valid) {
                                            continue;
                                          }
                                          IkReal x99 = x106.value;
                                          IkReal x100 = ((1.0) * x91 * x98);
                                          if ((((1.0) + (((-1.0) * (x94 * x94) *
                                                          (x97 * x97))))) <
                                              -0.00001)
                                            continue;
                                          IkReal x101 = IKsqrt(
                                              ((1.0) + (((-1.0) * (x94 * x94) *
                                                         (x97 * x97)))));
                                          IkReal x102 =
                                              ((1.0) * x94 * x97 * x98);
                                          IkReal x103 = (x101 * x98);
                                          px = 0;
                                          py = 0;
                                          pp = x92;
                                          sj2 = gconst1;
                                          cj2 = gconst2;
                                          CheckValue<IkReal> x107 =
                                              IKatan2WithCheck(
                                                  IkReal(
                                                      ((-1.0) +
                                                       (((-1818181.81818182) *
                                                         pz)))),
                                                  IkReal(((-273224.043715847) +
                                                          x90)),
                                                  IKFAST_ATAN2_MAGTHRESH);
                                          if (!x107.valid) {
                                            continue;
                                          }
                                          if (((x99 * (((1.12941176) +
                                                        (((1020408.16326531) *
                                                          pz)))))) <
                                                  -1 - IKFAST_SINCOS_THRESH ||
                                              ((x99 * (((1.12941176) +
                                                        (((1020408.16326531) *
                                                          pz)))))) >
                                                  1 + IKFAST_SINCOS_THRESH)
                                            continue;
                                          j2 = ((((-1.0) * (x107.value))) +
                                                (((-1.0) *
                                                  (IKasin(
                                                      (x99 *
                                                       (((1.12941176) +
                                                         (((1020408.16326531) *
                                                           pz))))))))));
                                          CheckValue<IkReal> x108 =
                                              IKatan2WithCheck(
                                                  IkReal(x91), IkReal(x93),
                                                  IKFAST_ATAN2_MAGTHRESH);
                                          if (!x108.valid) {
                                            continue;
                                          }
                                          if (((x94 * x99)) <
                                                  -1 - IKFAST_SINCOS_THRESH ||
                                              ((x94 * x99)) >
                                                  1 + IKFAST_SINCOS_THRESH)
                                            continue;
                                          IkReal gconst0 =
                                              ((((-1.0) * (x108.value))) +
                                               (((-1.0) *
                                                 (IKasin((x94 * x99))))));
                                          IkReal gconst1 =
                                              ((((-1.0) * x102 * x93)) +
                                               (((-1.0) * x100 * x101)));
                                          IkReal gconst2 =
                                              ((((-1.0) * x100 * x94 * x97)) +
                                               ((x103 * x93)));
                                          IkReal x109 = pz * pz;
                                          j1eval[0] =
                                              ((1.0) +
                                               (((44.4503709630156) * x109)));
                                          j1eval[1] = IKsign(
                                              ((224970001.0) +
                                               (((10000000000.0) * x109))));
                                          if (IKabs(j1eval[0]) <
                                                  0.0000010000000000 ||
                                              IKabs(j1eval[1]) <
                                                  0.0000010000000000) {
                                            {
                                              IkReal j1array[1], cj1array[1],
                                                  sj1array[1];
                                              bool j1valid[1] = {false};
                                              _nj1 = 1;
                                              IkReal x110 = gconst1 * gconst1;
                                              IkReal x111 = gconst2 * gconst2;
                                              IkReal x112 = (gconst2 * pz);
                                              IkReal x113 = (gconst1 * pz);
                                              IkReal x114 = (gconst1 * gconst2);
                                              CheckValue<IkReal> x115 =
                                                  IKPowWithIntegerCheck(
                                                      ((-194987000.0) +
                                                       (((9600.0) * pz)) +
                                                       (((-8500.0) * x112)) +
                                                       (((347976800.0) *
                                                         gconst2)) +
                                                       (((-2320000000.0) *
                                                         x113)) +
                                                       (((-1274.915) *
                                                         gconst1))),
                                                      -1);
                                              if (!x115.valid) {
                                                continue;
                                              }
                                              CheckValue<IkReal> x116 =
                                                  IKPowWithIntegerCheck(
                                                      ((-1949.87) +
                                                       (((-23200.0) * x113)) +
                                                       (((3479.768) *
                                                         gconst2)) +
                                                       (((-0.01274915) *
                                                         gconst1)) +
                                                       (((0.096) * pz)) +
                                                       (((-0.085) * x112))),
                                                      -1);
                                              if (!x116.valid) {
                                                continue;
                                              }
                                              if (IKabs((
                                                      (x115.value) *
                                                      (((224970000.990784) +
                                                        (((-3944.0) * x114)) +
                                                        (((4454.4) * gconst1)) +
                                                        (((0.01632) *
                                                          gconst2)) +
                                                        (((-0.007225) * x111)) +
                                                        (((-538240000.0) *
                                                          x110)))))) <
                                                      IKFAST_ATAN2_MAGTHRESH &&
                                                  IKabs((
                                                      (x116.value) *
                                                      (((-0.01248) +
                                                        (((-0.01972) * x111)) +
                                                        (((-5382.39999992775) *
                                                          x114)) +
                                                        (((0.01972) * x110)) +
                                                        (((3015.9999999184) *
                                                          gconst1)) +
                                                        (((0.033322) *
                                                          gconst2)) +
                                                        (((14999.0) * pz)))))) <
                                                      IKFAST_ATAN2_MAGTHRESH &&
                                                  IKabs(
                                                      IKsqr((
                                                          (x115.value) *
                                                          (((224970000.990784) +
                                                            (((-3944.0) *
                                                              x114)) +
                                                            (((4454.4) *
                                                              gconst1)) +
                                                            (((0.01632) *
                                                              gconst2)) +
                                                            (((-0.007225) *
                                                              x111)) +
                                                            (((-538240000.0) *
                                                              x110)))))) +
                                                      IKsqr((
                                                          (x116.value) *
                                                          (((-0.01248) +
                                                            (((-0.01972) *
                                                              x111)) +
                                                            (((-5382.39999992775) *
                                                              x114)) +
                                                            (((0.01972) *
                                                              x110)) +
                                                            (((3015.9999999184) *
                                                              gconst1)) +
                                                            (((0.033322) *
                                                              gconst2)) +
                                                            (((14999.0) *
                                                              pz)))))) -
                                                      1) <=
                                                      IKFAST_SINCOS_THRESH)
                                                continue;
                                              j1array[0] = IKatan2(
                                                  ((x115.value) *
                                                   (((224970000.990784) +
                                                     (((-3944.0) * x114)) +
                                                     (((4454.4) * gconst1)) +
                                                     (((0.01632) * gconst2)) +
                                                     (((-0.007225) * x111)) +
                                                     (((-538240000.0) *
                                                       x110))))),
                                                  ((x116.value) *
                                                   (((-0.01248) +
                                                     (((-0.01972) * x111)) +
                                                     (((-5382.39999992775) *
                                                       x114)) +
                                                     (((0.01972) * x110)) +
                                                     (((3015.9999999184) *
                                                       gconst1)) +
                                                     (((0.033322) * gconst2)) +
                                                     (((14999.0) * pz))))));
                                              sj1array[0] = IKsin(j1array[0]);
                                              cj1array[0] = IKcos(j1array[0]);
                                              if (j1array[0] > IKPI) {
                                                j1array[0] -= IK2PI;
                                              } else if (j1array[0] < -IKPI) {
                                                j1array[0] += IK2PI;
                                              }
                                              j1valid[0] = true;
                                              for (int ij1 = 0; ij1 < 1;
                                                   ++ij1) {
                                                if (!j1valid[ij1]) {
                                                  continue;
                                                }
                                                _ij1[0] = ij1;
                                                _ij1[1] = -1;
                                                for (int iij1 = ij1 + 1;
                                                     iij1 < 1; ++iij1) {
                                                  if (j1valid[iij1] &&
                                                      IKabs(cj1array[ij1] -
                                                            cj1array[iij1]) <
                                                          IKFAST_SOLUTION_THRESH &&
                                                      IKabs(sj1array[ij1] -
                                                            sj1array[iij1]) <
                                                          IKFAST_SOLUTION_THRESH) {
                                                    j1valid[iij1] = false;
                                                    _ij1[1] = iij1;
                                                    break;
                                                  }
                                                }
                                                j1 = j1array[ij1];
                                                cj1 = cj1array[ij1];
                                                sj1 = sj1array[ij1];
                                                {
                                                  IkReal evalcond[5];
                                                  IkReal x117 = IKcos(j1);
                                                  IkReal x118 = IKsin(j1);
                                                  IkReal x119 = ((1.0) * pz);
                                                  IkReal x120 =
                                                      ((0.232) * gconst1);
                                                  IkReal x121 =
                                                      ((8.5e-7) * gconst2);
                                                  IkReal x122 =
                                                      ((8.5e-7) * gconst1);
                                                  IkReal x123 =
                                                      ((0.232) * gconst2);
                                                  evalcond[0] =
                                                      ((-0.13) +
                                                       (((-1.0) * x117 *
                                                         x119)) +
                                                       (((-1.0) * x122)) +
                                                       x123 +
                                                       (((-0.14999) * x118)));
                                                  evalcond[1] =
                                                      ((-9.6e-7) +
                                                       (((-1.0) * x118 *
                                                         x119)) +
                                                       (((0.14999) * x117)) +
                                                       x120 + x121);
                                                  evalcond[2] =
                                                      ((0.0172892498998009) +
                                                       (((-0.26) * pz * x117)) +
                                                       (((2.879808e-7) *
                                                         x117)) +
                                                       (((-0.0389974) * x118)) +
                                                       (((-1.92e-6) * pz *
                                                         x118)) +
                                                       (((-1.0) * pz * x119)));
                                                  evalcond[3] =
                                                      ((-0.14999) +
                                                       (((-1.0) * x118 *
                                                         x122)) +
                                                       (((9.6e-7) * x117)) +
                                                       (((-0.13) * x118)) +
                                                       (((-1.0) * x117 *
                                                         x120)) +
                                                       (((-1.0) * x117 *
                                                         x121)) +
                                                       ((x118 * x123)));
                                                  evalcond[4] =
                                                      (((x117 * x123)) +
                                                       (((-9.6e-7) * x118)) +
                                                       (((-0.13) * x117)) +
                                                       (((-1.0) * x119)) +
                                                       (((-1.0) * x117 *
                                                         x122)) +
                                                       ((x118 * x121)) +
                                                       ((x118 * x120)));
                                                  if (IKabs(evalcond[0]) >
                                                          IKFAST_EVALCOND_THRESH ||
                                                      IKabs(evalcond[1]) >
                                                          IKFAST_EVALCOND_THRESH ||
                                                      IKabs(evalcond[2]) >
                                                          IKFAST_EVALCOND_THRESH ||
                                                      IKabs(evalcond[3]) >
                                                          IKFAST_EVALCOND_THRESH ||
                                                      IKabs(evalcond[4]) >
                                                          IKFAST_EVALCOND_THRESH) {
                                                    continue;
                                                  }
                                                }

                                                {
                                                  IkReal j0array[1],
                                                      cj0array[1], sj0array[1];
                                                  bool j0valid[1] = {false};
                                                  _nj0 = 1;
                                                  j0array[0] = 0;
                                                  sj0array[0] =
                                                      IKsin(j0array[0]);
                                                  cj0array[0] =
                                                      IKcos(j0array[0]);
                                                  if (j0array[0] > IKPI) {
                                                    j0array[0] -= IK2PI;
                                                  } else if (j0array[0] <
                                                             -IKPI) {
                                                    j0array[0] += IK2PI;
                                                  }
                                                  j0valid[0] = true;
                                                  for (int ij0 = 0; ij0 < 1;
                                                       ++ij0) {
                                                    if (!j0valid[ij0]) {
                                                      continue;
                                                    }
                                                    _ij0[0] = ij0;
                                                    _ij0[1] = -1;
                                                    for (int iij0 = ij0 + 1;
                                                         iij0 < 1; ++iij0) {
                                                      if (j0valid[iij0] &&
                                                          IKabs(
                                                              cj0array[ij0] -
                                                              cj0array[iij0]) <
                                                              IKFAST_SOLUTION_THRESH &&
                                                          IKabs(
                                                              sj0array[ij0] -
                                                              sj0array[iij0]) <
                                                              IKFAST_SOLUTION_THRESH) {
                                                        j0valid[iij0] = false;
                                                        _ij0[1] = iij0;
                                                        break;
                                                      }
                                                    }
                                                    j0 = j0array[ij0];
                                                    cj0 = cj0array[ij0];
                                                    sj0 = sj0array[ij0];

                                                    {
                                                      std::vector<
                                                          IkSingleDOFSolutionBase<
                                                              IkReal>>
                                                          vinfos(3);
                                                      vinfos[0].jointtype = 1;
                                                      vinfos[0].foffset = j0;
                                                      vinfos[0].indices[0] =
                                                          _ij0[0];
                                                      vinfos[0].indices[1] =
                                                          _ij0[1];
                                                      vinfos[0].maxsolutions =
                                                          _nj0;
                                                      vinfos[1].jointtype = 1;
                                                      vinfos[1].foffset = j1;
                                                      vinfos[1].indices[0] =
                                                          _ij1[0];
                                                      vinfos[1].indices[1] =
                                                          _ij1[1];
                                                      vinfos[1].maxsolutions =
                                                          _nj1;
                                                      vinfos[2].jointtype = 1;
                                                      vinfos[2].foffset = j2;
                                                      vinfos[2].indices[0] =
                                                          _ij2[0];
                                                      vinfos[2].indices[1] =
                                                          _ij2[1];
                                                      vinfos[2].maxsolutions =
                                                          _nj2;
                                                      std::vector<int> vfree(0);
                                                      solutions.AddSolution(
                                                          vinfos, vfree);
                                                    }
                                                  }
                                                }
                                              }
                                            }

                                          } else {
                                            {
                                              IkReal j1array[1], cj1array[1],
                                                  sj1array[1];
                                              bool j1valid[1] = {false};
                                              _nj1 = 1;
                                              IkReal x124 = (gconst2 * pz);
                                              IkReal x125 = (gconst1 * pz);
                                              CheckValue<IkReal> x126 =
                                                  IKatan2WithCheck(
                                                      IkReal(
                                                          ((-194987000.0) +
                                                           (((347976800.0) *
                                                             gconst2)) +
                                                           (((2320000000.0) *
                                                             x125)) +
                                                           (((-9600.0) * pz)) +
                                                           (((8500.0) * x124)) +
                                                           (((-1274.915) *
                                                             gconst1)))),
                                                      IkReal((
                                                          (1439.904) +
                                                          (((-8500.0) * x125)) +
                                                          (((2320000000.0) *
                                                            x124)) +
                                                          (((-1300000000.0) *
                                                            pz)) +
                                                          (((-347976800.0) *
                                                            gconst1)) +
                                                          (((-1274.915) *
                                                            gconst2)))),
                                                      IKFAST_ATAN2_MAGTHRESH);
                                              if (!x126.valid) {
                                                continue;
                                              }
                                              CheckValue<IkReal> x127 =
                                                  IKPowWithIntegerCheck(
                                                      IKsign(
                                                          ((224970001.0) +
                                                           (((10000000000.0) *
                                                             (pz * pz))))),
                                                      -1);
                                              if (!x127.valid) {
                                                continue;
                                              }
                                              j1array[0] =
                                                  ((-1.5707963267949) +
                                                   (x126.value) +
                                                   (((1.5707963267949) *
                                                     (x127.value))));
                                              sj1array[0] = IKsin(j1array[0]);
                                              cj1array[0] = IKcos(j1array[0]);
                                              if (j1array[0] > IKPI) {
                                                j1array[0] -= IK2PI;
                                              } else if (j1array[0] < -IKPI) {
                                                j1array[0] += IK2PI;
                                              }
                                              j1valid[0] = true;
                                              for (int ij1 = 0; ij1 < 1;
                                                   ++ij1) {
                                                if (!j1valid[ij1]) {
                                                  continue;
                                                }
                                                _ij1[0] = ij1;
                                                _ij1[1] = -1;
                                                for (int iij1 = ij1 + 1;
                                                     iij1 < 1; ++iij1) {
                                                  if (j1valid[iij1] &&
                                                      IKabs(cj1array[ij1] -
                                                            cj1array[iij1]) <
                                                          IKFAST_SOLUTION_THRESH &&
                                                      IKabs(sj1array[ij1] -
                                                            sj1array[iij1]) <
                                                          IKFAST_SOLUTION_THRESH) {
                                                    j1valid[iij1] = false;
                                                    _ij1[1] = iij1;
                                                    break;
                                                  }
                                                }
                                                j1 = j1array[ij1];
                                                cj1 = cj1array[ij1];
                                                sj1 = sj1array[ij1];
                                                {
                                                  IkReal evalcond[5];
                                                  IkReal x128 = IKcos(j1);
                                                  IkReal x129 = IKsin(j1);
                                                  IkReal x130 = ((1.0) * pz);
                                                  IkReal x131 =
                                                      ((0.232) * gconst1);
                                                  IkReal x132 =
                                                      ((8.5e-7) * gconst2);
                                                  IkReal x133 =
                                                      ((8.5e-7) * gconst1);
                                                  IkReal x134 =
                                                      ((0.232) * gconst2);
                                                  evalcond[0] =
                                                      ((-0.13) + x134 +
                                                       (((-1.0) * x128 *
                                                         x130)) +
                                                       (((-1.0) * x133)) +
                                                       (((-0.14999) * x129)));
                                                  evalcond[1] =
                                                      ((-9.6e-7) + x131 + x132 +
                                                       (((0.14999) * x128)) +
                                                       (((-1.0) * x129 *
                                                         x130)));
                                                  evalcond[2] =
                                                      ((0.0172892498998009) +
                                                       (((-1.0) * pz * x130)) +
                                                       (((-0.26) * pz * x128)) +
                                                       (((-0.0389974) * x129)) +
                                                       (((2.879808e-7) *
                                                         x128)) +
                                                       (((-1.92e-6) * pz *
                                                         x129)));
                                                  evalcond[3] =
                                                      ((-0.14999) +
                                                       ((x129 * x134)) +
                                                       (((9.6e-7) * x128)) +
                                                       (((-0.13) * x129)) +
                                                       (((-1.0) * x128 *
                                                         x132)) +
                                                       (((-1.0) * x128 *
                                                         x131)) +
                                                       (((-1.0) * x129 *
                                                         x133)));
                                                  evalcond[4] =
                                                      (((x129 * x131)) +
                                                       ((x129 * x132)) +
                                                       (((-0.13) * x128)) +
                                                       (((-9.6e-7) * x129)) +
                                                       (((-1.0) * x128 *
                                                         x133)) +
                                                       (((-1.0) * x130)) +
                                                       ((x128 * x134)));
                                                  if (IKabs(evalcond[0]) >
                                                          IKFAST_EVALCOND_THRESH ||
                                                      IKabs(evalcond[1]) >
                                                          IKFAST_EVALCOND_THRESH ||
                                                      IKabs(evalcond[2]) >
                                                          IKFAST_EVALCOND_THRESH ||
                                                      IKabs(evalcond[3]) >
                                                          IKFAST_EVALCOND_THRESH ||
                                                      IKabs(evalcond[4]) >
                                                          IKFAST_EVALCOND_THRESH) {
                                                    continue;
                                                  }
                                                }

                                                {
                                                  IkReal j0array[1],
                                                      cj0array[1], sj0array[1];
                                                  bool j0valid[1] = {false};
                                                  _nj0 = 1;
                                                  j0array[0] = 0;
                                                  sj0array[0] =
                                                      IKsin(j0array[0]);
                                                  cj0array[0] =
                                                      IKcos(j0array[0]);
                                                  if (j0array[0] > IKPI) {
                                                    j0array[0] -= IK2PI;
                                                  } else if (j0array[0] <
                                                             -IKPI) {
                                                    j0array[0] += IK2PI;
                                                  }
                                                  j0valid[0] = true;
                                                  for (int ij0 = 0; ij0 < 1;
                                                       ++ij0) {
                                                    if (!j0valid[ij0]) {
                                                      continue;
                                                    }
                                                    _ij0[0] = ij0;
                                                    _ij0[1] = -1;
                                                    for (int iij0 = ij0 + 1;
                                                         iij0 < 1; ++iij0) {
                                                      if (j0valid[iij0] &&
                                                          IKabs(
                                                              cj0array[ij0] -
                                                              cj0array[iij0]) <
                                                              IKFAST_SOLUTION_THRESH &&
                                                          IKabs(
                                                              sj0array[ij0] -
                                                              sj0array[iij0]) <
                                                              IKFAST_SOLUTION_THRESH) {
                                                        j0valid[iij0] = false;
                                                        _ij0[1] = iij0;
                                                        break;
                                                      }
                                                    }
                                                    j0 = j0array[ij0];
                                                    cj0 = cj0array[ij0];
                                                    sj0 = sj0array[ij0];

                                                    {
                                                      std::vector<
                                                          IkSingleDOFSolutionBase<
                                                              IkReal>>
                                                          vinfos(3);
                                                      vinfos[0].jointtype = 1;
                                                      vinfos[0].foffset = j0;
                                                      vinfos[0].indices[0] =
                                                          _ij0[0];
                                                      vinfos[0].indices[1] =
                                                          _ij0[1];
                                                      vinfos[0].maxsolutions =
                                                          _nj0;
                                                      vinfos[1].jointtype = 1;
                                                      vinfos[1].foffset = j1;
                                                      vinfos[1].indices[0] =
                                                          _ij1[0];
                                                      vinfos[1].indices[1] =
                                                          _ij1[1];
                                                      vinfos[1].maxsolutions =
                                                          _nj1;
                                                      vinfos[2].jointtype = 1;
                                                      vinfos[2].foffset = j2;
                                                      vinfos[2].indices[0] =
                                                          _ij2[0];
                                                      vinfos[2].indices[1] =
                                                          _ij2[1];
                                                      vinfos[2].maxsolutions =
                                                          _nj2;
                                                      std::vector<int> vfree(0);
                                                      solutions.AddSolution(
                                                          vinfos, vfree);
                                                    }
                                                  }
                                                }
                                              }
                                            }
                                          }
                                        }
                                      }
                                    } while (0);
                                    if (bgotonextstatement) {
                                      bool bgotonextstatement = true;
                                      do {
                                        IkReal x135 =
                                            ((-1.0) +
                                             (((-1819729.15841448) * pz)));
                                        IkReal x136 =
                                            ((-272941.176470588) +
                                             (((6.66711114074272) * pz)));
                                        IkReal x137 =
                                            ((1.12941176470588) +
                                             (((1019675.82152536) * pz)));
                                        IkReal x138 =
                                            ((x135 * x135) + (x136 * x136));
                                        if ((((74496885814.1488) +
                                              (((3311414210028.33) *
                                                (pz * pz))))) < -0.00001)
                                          continue;
                                        IkReal x139 =
                                            IKabs(IKsqrt(((74496885814.1488) +
                                                          (((3311414210028.33) *
                                                            (pz * pz))))));
                                        CheckValue<IkReal> x146 =
                                            IKPowWithIntegerCheck(x139, -1);
                                        if (!x146.valid) {
                                          continue;
                                        }
                                        IkReal x140 = x146.value;
                                        IkReal x147 = x138;
                                        if (IKabs(x147) == 0) {
                                          continue;
                                        }
                                        IkReal x141 = pow(x147, -0.5);
                                        IkReal x142 = ((1.0) * x141);
                                        IkReal x143 = (x137 * x140);
                                        IkReal x144 = (x135 * x142);
                                        if ((((1.0) +
                                              (((-1.0) * (x143 * x143))))) <
                                            -0.00001)
                                          continue;
                                        IkReal x145 = IKsqrt(
                                            ((1.0) +
                                             (((-1.0) * (x143 * x143)))));
                                        CheckValue<IkReal> x148 =
                                            IKatan2WithCheck(
                                                IkReal(x135), IkReal(x136),
                                                IKFAST_ATAN2_MAGTHRESH);
                                        if (!x148.valid) {
                                          continue;
                                        }
                                        if ((x138) < -0.00001)
                                          continue;
                                        CheckValue<IkReal> x149 =
                                            IKPowWithIntegerCheck(
                                                IKabs(IKsqrt(x138)), -1);
                                        if (!x149.valid) {
                                          continue;
                                        }
                                        if (((x137 * (x149.value))) <
                                                -1 - IKFAST_SINCOS_THRESH ||
                                            ((x137 * (x149.value))) >
                                                1 + IKFAST_SINCOS_THRESH)
                                          continue;
                                        IkReal gconst3 =
                                            ((3.14159265358979) +
                                             (((-1.0) * (x148.value))) +
                                             (IKasin((x137 * (x149.value)))));
                                        IkReal gconst4 =
                                            (((x144 * x145)) +
                                             (((-1.0) * x136 * x142 * x143)));
                                        IkReal gconst5 =
                                            ((((-1.0) * x143 * x144)) +
                                             (((-1.0) * x136 * x142 * x145)));
                                        IkReal x150 = x135;
                                        IkReal x151 =
                                            ((-272941.176470588) +
                                             (((6.66711114074272) * pz)));
                                        CheckValue<IkReal> x152 =
                                            IKatan2WithCheck(
                                                IkReal(x150), IkReal(x151),
                                                IKFAST_ATAN2_MAGTHRESH);
                                        if (!x152.valid) {
                                          continue;
                                        }
                                        if ((((x150 * x150) + (x151 * x151))) <
                                            -0.00001)
                                          continue;
                                        CheckValue<IkReal> x153 =
                                            IKPowWithIntegerCheck(
                                                IKabs(IKsqrt(((x150 * x150) +
                                                              (x151 * x151)))),
                                                -1);
                                        if (!x153.valid) {
                                          continue;
                                        }
                                        if ((((x153.value) *
                                              (((1.12941176470588) +
                                                (((1019675.82152536) *
                                                  pz)))))) <
                                                -1 - IKFAST_SINCOS_THRESH ||
                                            (((x153.value) *
                                              (((1.12941176470588) +
                                                (((1019675.82152536) *
                                                  pz)))))) >
                                                1 + IKFAST_SINCOS_THRESH)
                                          continue;
                                        evalcond[0] =
                                            ((-3.14159265358979) +
                                             (IKfmod(
                                                 ((3.14159265358979) +
                                                  (IKabs((
                                                      (-3.14159265358979) +
                                                      (x152.value) +
                                                      (((-1.0) *
                                                        (IKasin((
                                                            (x153.value) *
                                                            (((1.12941176470588) +
                                                              (((1019675.82152536) *
                                                                pz))))))))) +
                                                      j2)))),
                                                 6.28318530717959)));
                                        if (IKabs(evalcond[0]) <
                                            0.0000050000000000) {
                                          bgotonextstatement = false;
                                          {
                                            IkReal j1eval[2];
                                            IkReal x154 =
                                                ((6.66711114074272) * pz);
                                            IkReal x155 = x135;
                                            IkReal x156 = pz * pz;
                                            IkReal x157 =
                                                ((-272941.176470588) + x154);
                                            IkReal x158 =
                                                ((1.12941176470588) +
                                                 (((1019675.82152536) * pz)));
                                            IkReal x159 =
                                                ((x157 * x157) + (x155 * x155));
                                            IkReal x160 = x139;
                                            CheckValue<IkReal> x168 =
                                                IKPowWithIntegerCheck(x160, -1);
                                            if (!x168.valid) {
                                              continue;
                                            }
                                            IkReal x161 = x168.value;
                                            IkReal x169 = x159;
                                            if (IKabs(x169) == 0) {
                                              continue;
                                            }
                                            IkReal x162 = pow(x169, -0.5);
                                            IkReal x163 = ((1.0) * x162);
                                            IkReal x164 = (x158 * x161);
                                            if ((x159) < -0.00001)
                                              continue;
                                            CheckValue<IkReal> x170 =
                                                IKPowWithIntegerCheck(
                                                    IKabs(IKsqrt(x159)), -1);
                                            if (!x170.valid) {
                                              continue;
                                            }
                                            IkReal x165 = x170.value;
                                            IkReal x166 = (x155 * x163);
                                            if ((((1.0) +
                                                  (((-1.0) * (x164 * x164))))) <
                                                -0.00001)
                                              continue;
                                            IkReal x167 = IKsqrt(
                                                ((1.0) +
                                                 (((-1.0) * (x164 * x164)))));
                                            px = 0;
                                            py = 0;
                                            pp = x156;
                                            sj2 = gconst4;
                                            cj2 = gconst5;
                                            if (((x165 *
                                                  (((1.12941176) +
                                                    (((1020408.16326531) *
                                                      pz)))))) <
                                                    -1 - IKFAST_SINCOS_THRESH ||
                                                ((x165 *
                                                  (((1.12941176) +
                                                    (((1020408.16326531) *
                                                      pz)))))) >
                                                    1 + IKFAST_SINCOS_THRESH)
                                              continue;
                                            CheckValue<IkReal> x171 =
                                                IKatan2WithCheck(
                                                    IkReal(
                                                        ((-1.0) +
                                                         (((-1818181.81818182) *
                                                           pz)))),
                                                    IkReal(
                                                        ((-273224.043715847) +
                                                         x154)),
                                                    IKFAST_ATAN2_MAGTHRESH);
                                            if (!x171.valid) {
                                              continue;
                                            }
                                            j2 = ((3.14159265) +
                                                  (IKasin(
                                                      (x165 *
                                                       (((1.12941176) +
                                                         (((1020408.16326531) *
                                                           pz))))))) +
                                                  (((-1.0) * (x171.value))));
                                            CheckValue<IkReal> x172 =
                                                IKatan2WithCheck(
                                                    IkReal(x155), IkReal(x157),
                                                    IKFAST_ATAN2_MAGTHRESH);
                                            if (!x172.valid) {
                                              continue;
                                            }
                                            if (((x158 * x165)) <
                                                    -1 - IKFAST_SINCOS_THRESH ||
                                                ((x158 * x165)) >
                                                    1 + IKFAST_SINCOS_THRESH)
                                              continue;
                                            IkReal gconst3 =
                                                ((3.14159265358979) +
                                                 (((-1.0) * (x172.value))) +
                                                 (IKasin((x158 * x165))));
                                            IkReal gconst4 = (((x166 * x167)) +
                                                              (((-1.0) * x157 *
                                                                x163 * x164)));
                                            IkReal gconst5 =
                                                ((((-1.0) * x157 * x163 *
                                                   x167)) +
                                                 (((-1.0) * x164 * x166)));
                                            IkReal x173 = pz * pz;
                                            j1eval[0] =
                                                ((1.0) +
                                                 (((44.4503709630156) * x173)));
                                            j1eval[1] = IKsign(
                                                ((224970001.0) +
                                                 (((10000000000.0) * x173))));
                                            if (IKabs(j1eval[0]) <
                                                    0.0000010000000000 ||
                                                IKabs(j1eval[1]) <
                                                    0.0000010000000000) {
                                              {
                                                IkReal j1array[1], cj1array[1],
                                                    sj1array[1];
                                                bool j1valid[1] = {false};
                                                _nj1 = 1;
                                                IkReal x174 = gconst4 * gconst4;
                                                IkReal x175 = gconst5 * gconst5;
                                                IkReal x176 =
                                                    (gconst4 * gconst5);
                                                IkReal x177 = (gconst4 * pz);
                                                IkReal x178 = (gconst5 * pz);
                                                CheckValue<IkReal> x179 =
                                                    IKPowWithIntegerCheck(
                                                        ((0.01439904) +
                                                         (((13000.0) * pz)) +
                                                         (((-3479.768) *
                                                           gconst4)) +
                                                         (((-0.01274915) *
                                                           gconst5)) +
                                                         (((0.085) * x177)) +
                                                         (((-23200.0) * x178))),
                                                        -1);
                                                if (!x179.valid) {
                                                  continue;
                                                }
                                                CheckValue<IkReal> x180 =
                                                    IKPowWithIntegerCheck(
                                                        ((1439.904) +
                                                         (((-2320000000.0) *
                                                           x178)) +
                                                         (((8500.0) * x177)) +
                                                         (((-347976800.0) *
                                                           gconst4)) +
                                                         (((1300000000.0) *
                                                           pz)) +
                                                         (((-1274.915) *
                                                           gconst5))),
                                                        -1);
                                                if (!x180.valid) {
                                                  continue;
                                                }
                                                if (IKabs((
                                                        (x179.value) *
                                                        (((-0.01248) +
                                                          (((-5382.39999992775) *
                                                            x176)) +
                                                          (((0.01972) * x174)) +
                                                          (((-0.01972) *
                                                            x175)) +
                                                          (((3015.9999999184) *
                                                            gconst4)) +
                                                          (((-14999.0) * pz)) +
                                                          (((0.033322) *
                                                            gconst5)))))) <
                                                        IKFAST_ATAN2_MAGTHRESH &&
                                                    IKabs(
                                                        ((x180.value) *
                                                         (((55970001.0) +
                                                           (((-2210.0) *
                                                             gconst4)) +
                                                           (((-538240000.0) *
                                                             x175)) +
                                                           (((603200000.0) *
                                                             gconst5)) +
                                                           (((3944.0) * x176)) +
                                                           (((-0.007225) *
                                                             x174)))))) <
                                                        IKFAST_ATAN2_MAGTHRESH &&
                                                    IKabs(
                                                        IKsqr((
                                                            (x179.value) *
                                                            (((-0.01248) +
                                                              (((-5382.39999992775) *
                                                                x176)) +
                                                              (((0.01972) *
                                                                x174)) +
                                                              (((-0.01972) *
                                                                x175)) +
                                                              (((3015.9999999184) *
                                                                gconst4)) +
                                                              (((-14999.0) *
                                                                pz)) +
                                                              (((0.033322) *
                                                                gconst5)))))) +
                                                        IKsqr((
                                                            (x180.value) *
                                                            (((55970001.0) +
                                                              (((-2210.0) *
                                                                gconst4)) +
                                                              (((-538240000.0) *
                                                                x175)) +
                                                              (((603200000.0) *
                                                                gconst5)) +
                                                              (((3944.0) *
                                                                x176)) +
                                                              (((-0.007225) *
                                                                x174)))))) -
                                                        1) <=
                                                        IKFAST_SINCOS_THRESH)
                                                  continue;
                                                j1array[0] = IKatan2(
                                                    ((x179.value) *
                                                     (((-0.01248) +
                                                       (((-5382.39999992775) *
                                                         x176)) +
                                                       (((0.01972) * x174)) +
                                                       (((-0.01972) * x175)) +
                                                       (((3015.9999999184) *
                                                         gconst4)) +
                                                       (((-14999.0) * pz)) +
                                                       (((0.033322) *
                                                         gconst5))))),
                                                    ((x180.value) *
                                                     (((55970001.0) +
                                                       (((-2210.0) * gconst4)) +
                                                       (((-538240000.0) *
                                                         x175)) +
                                                       (((603200000.0) *
                                                         gconst5)) +
                                                       (((3944.0) * x176)) +
                                                       (((-0.007225) *
                                                         x174))))));
                                                sj1array[0] = IKsin(j1array[0]);
                                                cj1array[0] = IKcos(j1array[0]);
                                                if (j1array[0] > IKPI) {
                                                  j1array[0] -= IK2PI;
                                                } else if (j1array[0] < -IKPI) {
                                                  j1array[0] += IK2PI;
                                                }
                                                j1valid[0] = true;
                                                for (int ij1 = 0; ij1 < 1;
                                                     ++ij1) {
                                                  if (!j1valid[ij1]) {
                                                    continue;
                                                  }
                                                  _ij1[0] = ij1;
                                                  _ij1[1] = -1;
                                                  for (int iij1 = ij1 + 1;
                                                       iij1 < 1; ++iij1) {
                                                    if (j1valid[iij1] &&
                                                        IKabs(cj1array[ij1] -
                                                              cj1array[iij1]) <
                                                            IKFAST_SOLUTION_THRESH &&
                                                        IKabs(sj1array[ij1] -
                                                              sj1array[iij1]) <
                                                            IKFAST_SOLUTION_THRESH) {
                                                      j1valid[iij1] = false;
                                                      _ij1[1] = iij1;
                                                      break;
                                                    }
                                                  }
                                                  j1 = j1array[ij1];
                                                  cj1 = cj1array[ij1];
                                                  sj1 = sj1array[ij1];
                                                  {
                                                    IkReal evalcond[5];
                                                    IkReal x181 = IKcos(j1);
                                                    IkReal x182 = IKsin(j1);
                                                    IkReal x183 =
                                                        ((8.5e-7) * gconst5);
                                                    IkReal x184 =
                                                        ((0.232) * gconst4);
                                                    IkReal x185 = ((1.0) * pz);
                                                    IkReal x186 =
                                                        ((8.5e-7) * gconst4);
                                                    IkReal x187 =
                                                        ((0.232) * gconst5);
                                                    evalcond[0] =
                                                        ((-0.13) +
                                                         (((-0.14999) * x182)) +
                                                         (((-1.0) * x186)) +
                                                         x187 +
                                                         (((-1.0) * x181 *
                                                           x185)));
                                                    evalcond[1] =
                                                        ((-9.6e-7) +
                                                         (((0.14999) * x181)) +
                                                         x184 + x183 +
                                                         (((-1.0) * x182 *
                                                           x185)));
                                                    evalcond[2] =
                                                        ((0.0172892498998009) +
                                                         (((-1.92e-6) * pz *
                                                           x182)) +
                                                         (((-1.0) * pz *
                                                           x185)) +
                                                         (((-0.0389974) *
                                                           x182)) +
                                                         (((-0.26) * pz *
                                                           x181)) +
                                                         (((2.879808e-7) *
                                                           x181)));
                                                    evalcond[3] =
                                                        ((-0.14999) +
                                                         (((9.6e-7) * x181)) +
                                                         (((-0.13) * x182)) +
                                                         (((-1.0) * x182 *
                                                           x186)) +
                                                         ((x182 * x187)) +
                                                         (((-1.0) * x181 *
                                                           x184)) +
                                                         (((-1.0) * x181 *
                                                           x183)));
                                                    evalcond[4] =
                                                        ((((-1.0) * x185)) +
                                                         (((-0.13) * x181)) +
                                                         (((-9.6e-7) * x182)) +
                                                         ((x182 * x183)) +
                                                         ((x182 * x184)) +
                                                         ((x181 * x187)) +
                                                         (((-1.0) * x181 *
                                                           x186)));
                                                    if (IKabs(evalcond[0]) >
                                                            IKFAST_EVALCOND_THRESH ||
                                                        IKabs(evalcond[1]) >
                                                            IKFAST_EVALCOND_THRESH ||
                                                        IKabs(evalcond[2]) >
                                                            IKFAST_EVALCOND_THRESH ||
                                                        IKabs(evalcond[3]) >
                                                            IKFAST_EVALCOND_THRESH ||
                                                        IKabs(evalcond[4]) >
                                                            IKFAST_EVALCOND_THRESH) {
                                                      continue;
                                                    }
                                                  }

                                                  {
                                                    IkReal j0array[1],
                                                        cj0array[1],
                                                        sj0array[1];
                                                    bool j0valid[1] = {false};
                                                    _nj0 = 1;
                                                    j0array[0] = 0;
                                                    sj0array[0] =
                                                        IKsin(j0array[0]);
                                                    cj0array[0] =
                                                        IKcos(j0array[0]);
                                                    if (j0array[0] > IKPI) {
                                                      j0array[0] -= IK2PI;
                                                    } else if (j0array[0] <
                                                               -IKPI) {
                                                      j0array[0] += IK2PI;
                                                    }
                                                    j0valid[0] = true;
                                                    for (int ij0 = 0; ij0 < 1;
                                                         ++ij0) {
                                                      if (!j0valid[ij0]) {
                                                        continue;
                                                      }
                                                      _ij0[0] = ij0;
                                                      _ij0[1] = -1;
                                                      for (int iij0 = ij0 + 1;
                                                           iij0 < 1; ++iij0) {
                                                        if (j0valid[iij0] &&
                                                            IKabs(
                                                                cj0array[ij0] -
                                                                cj0array
                                                                    [iij0]) <
                                                                IKFAST_SOLUTION_THRESH &&
                                                            IKabs(
                                                                sj0array[ij0] -
                                                                sj0array
                                                                    [iij0]) <
                                                                IKFAST_SOLUTION_THRESH) {
                                                          j0valid[iij0] = false;
                                                          _ij0[1] = iij0;
                                                          break;
                                                        }
                                                      }
                                                      j0 = j0array[ij0];
                                                      cj0 = cj0array[ij0];
                                                      sj0 = sj0array[ij0];

                                                      {
                                                        std::vector<
                                                            IkSingleDOFSolutionBase<
                                                                IkReal>>
                                                            vinfos(3);
                                                        vinfos[0].jointtype = 1;
                                                        vinfos[0].foffset = j0;
                                                        vinfos[0].indices[0] =
                                                            _ij0[0];
                                                        vinfos[0].indices[1] =
                                                            _ij0[1];
                                                        vinfos[0].maxsolutions =
                                                            _nj0;
                                                        vinfos[1].jointtype = 1;
                                                        vinfos[1].foffset = j1;
                                                        vinfos[1].indices[0] =
                                                            _ij1[0];
                                                        vinfos[1].indices[1] =
                                                            _ij1[1];
                                                        vinfos[1].maxsolutions =
                                                            _nj1;
                                                        vinfos[2].jointtype = 1;
                                                        vinfos[2].foffset = j2;
                                                        vinfos[2].indices[0] =
                                                            _ij2[0];
                                                        vinfos[2].indices[1] =
                                                            _ij2[1];
                                                        vinfos[2].maxsolutions =
                                                            _nj2;
                                                        std::vector<int> vfree(
                                                            0);
                                                        solutions.AddSolution(
                                                            vinfos, vfree);
                                                      }
                                                    }
                                                  }
                                                }
                                              }

                                            } else {
                                              {
                                                IkReal j1array[1], cj1array[1],
                                                    sj1array[1];
                                                bool j1valid[1] = {false};
                                                _nj1 = 1;
                                                IkReal x188 = (gconst4 * pz);
                                                IkReal x189 = (gconst5 * pz);
                                                CheckValue<IkReal> x190 =
                                                    IKPowWithIntegerCheck(
                                                        IKsign(
                                                            ((224970001.0) +
                                                             (((10000000000.0) *
                                                               (pz * pz))))),
                                                        -1);
                                                if (!x190.valid) {
                                                  continue;
                                                }
                                                CheckValue<IkReal> x191 =
                                                    IKatan2WithCheck(
                                                        IkReal((
                                                            (-194987000.0) +
                                                            (((347976800.0) *
                                                              gconst5)) +
                                                            (((8500.0) *
                                                              x189)) +
                                                            (((2320000000.0) *
                                                              x188)) +
                                                            (((-9600.0) * pz)) +
                                                            (((-1274.915) *
                                                              gconst4)))),
                                                        IkReal(
                                                            ((1439.904) +
                                                             (((-8500.0) *
                                                               x188)) +
                                                             (((2320000000.0) *
                                                               x189)) +
                                                             (((-1300000000.0) *
                                                               pz)) +
                                                             (((-347976800.0) *
                                                               gconst4)) +
                                                             (((-1274.915) *
                                                               gconst5)))),
                                                        IKFAST_ATAN2_MAGTHRESH);
                                                if (!x191.valid) {
                                                  continue;
                                                }
                                                j1array[0] =
                                                    ((-1.5707963267949) +
                                                     (((1.5707963267949) *
                                                       (x190.value))) +
                                                     (x191.value));
                                                sj1array[0] = IKsin(j1array[0]);
                                                cj1array[0] = IKcos(j1array[0]);
                                                if (j1array[0] > IKPI) {
                                                  j1array[0] -= IK2PI;
                                                } else if (j1array[0] < -IKPI) {
                                                  j1array[0] += IK2PI;
                                                }
                                                j1valid[0] = true;
                                                for (int ij1 = 0; ij1 < 1;
                                                     ++ij1) {
                                                  if (!j1valid[ij1]) {
                                                    continue;
                                                  }
                                                  _ij1[0] = ij1;
                                                  _ij1[1] = -1;
                                                  for (int iij1 = ij1 + 1;
                                                       iij1 < 1; ++iij1) {
                                                    if (j1valid[iij1] &&
                                                        IKabs(cj1array[ij1] -
                                                              cj1array[iij1]) <
                                                            IKFAST_SOLUTION_THRESH &&
                                                        IKabs(sj1array[ij1] -
                                                              sj1array[iij1]) <
                                                            IKFAST_SOLUTION_THRESH) {
                                                      j1valid[iij1] = false;
                                                      _ij1[1] = iij1;
                                                      break;
                                                    }
                                                  }
                                                  j1 = j1array[ij1];
                                                  cj1 = cj1array[ij1];
                                                  sj1 = sj1array[ij1];
                                                  {
                                                    IkReal evalcond[5];
                                                    IkReal x192 = IKcos(j1);
                                                    IkReal x193 = IKsin(j1);
                                                    IkReal x194 =
                                                        ((8.5e-7) * gconst5);
                                                    IkReal x195 =
                                                        ((0.232) * gconst4);
                                                    IkReal x196 = ((1.0) * pz);
                                                    IkReal x197 =
                                                        ((8.5e-7) * gconst4);
                                                    IkReal x198 =
                                                        ((0.232) * gconst5);
                                                    evalcond[0] =
                                                        ((-0.13) +
                                                         (((-1.0) * x192 *
                                                           x196)) +
                                                         (((-0.14999) * x193)) +
                                                         (((-1.0) * x197)) +
                                                         x198);
                                                    evalcond[1] =
                                                        ((-9.6e-7) +
                                                         (((0.14999) * x192)) +
                                                         x195 + x194 +
                                                         (((-1.0) * x193 *
                                                           x196)));
                                                    evalcond[2] =
                                                        ((0.0172892498998009) +
                                                         (((-1.0) * pz *
                                                           x196)) +
                                                         (((-1.92e-6) * pz *
                                                           x193)) +
                                                         (((-0.0389974) *
                                                           x193)) +
                                                         (((-0.26) * pz *
                                                           x192)) +
                                                         (((2.879808e-7) *
                                                           x192)));
                                                    evalcond[3] =
                                                        ((-0.14999) +
                                                         (((9.6e-7) * x192)) +
                                                         ((x193 * x198)) +
                                                         (((-1.0) * x192 *
                                                           x194)) +
                                                         (((-1.0) * x192 *
                                                           x195)) +
                                                         (((-0.13) * x193)) +
                                                         (((-1.0) * x193 *
                                                           x197)));
                                                    evalcond[4] =
                                                        (((x193 * x194)) +
                                                         ((x193 * x195)) +
                                                         (((-1.0) * x192 *
                                                           x197)) +
                                                         (((-1.0) * x196)) +
                                                         ((x192 * x198)) +
                                                         (((-0.13) * x192)) +
                                                         (((-9.6e-7) * x193)));
                                                    if (IKabs(evalcond[0]) >
                                                            IKFAST_EVALCOND_THRESH ||
                                                        IKabs(evalcond[1]) >
                                                            IKFAST_EVALCOND_THRESH ||
                                                        IKabs(evalcond[2]) >
                                                            IKFAST_EVALCOND_THRESH ||
                                                        IKabs(evalcond[3]) >
                                                            IKFAST_EVALCOND_THRESH ||
                                                        IKabs(evalcond[4]) >
                                                            IKFAST_EVALCOND_THRESH) {
                                                      continue;
                                                    }
                                                  }

                                                  {
                                                    IkReal j0array[1],
                                                        cj0array[1],
                                                        sj0array[1];
                                                    bool j0valid[1] = {false};
                                                    _nj0 = 1;
                                                    j0array[0] = 0;
                                                    sj0array[0] =
                                                        IKsin(j0array[0]);
                                                    cj0array[0] =
                                                        IKcos(j0array[0]);
                                                    if (j0array[0] > IKPI) {
                                                      j0array[0] -= IK2PI;
                                                    } else if (j0array[0] <
                                                               -IKPI) {
                                                      j0array[0] += IK2PI;
                                                    }
                                                    j0valid[0] = true;
                                                    for (int ij0 = 0; ij0 < 1;
                                                         ++ij0) {
                                                      if (!j0valid[ij0]) {
                                                        continue;
                                                      }
                                                      _ij0[0] = ij0;
                                                      _ij0[1] = -1;
                                                      for (int iij0 = ij0 + 1;
                                                           iij0 < 1; ++iij0) {
                                                        if (j0valid[iij0] &&
                                                            IKabs(
                                                                cj0array[ij0] -
                                                                cj0array
                                                                    [iij0]) <
                                                                IKFAST_SOLUTION_THRESH &&
                                                            IKabs(
                                                                sj0array[ij0] -
                                                                sj0array
                                                                    [iij0]) <
                                                                IKFAST_SOLUTION_THRESH) {
                                                          j0valid[iij0] = false;
                                                          _ij0[1] = iij0;
                                                          break;
                                                        }
                                                      }
                                                      j0 = j0array[ij0];
                                                      cj0 = cj0array[ij0];
                                                      sj0 = sj0array[ij0];

                                                      {
                                                        std::vector<
                                                            IkSingleDOFSolutionBase<
                                                                IkReal>>
                                                            vinfos(3);
                                                        vinfos[0].jointtype = 1;
                                                        vinfos[0].foffset = j0;
                                                        vinfos[0].indices[0] =
                                                            _ij0[0];
                                                        vinfos[0].indices[1] =
                                                            _ij0[1];
                                                        vinfos[0].maxsolutions =
                                                            _nj0;
                                                        vinfos[1].jointtype = 1;
                                                        vinfos[1].foffset = j1;
                                                        vinfos[1].indices[0] =
                                                            _ij1[0];
                                                        vinfos[1].indices[1] =
                                                            _ij1[1];
                                                        vinfos[1].maxsolutions =
                                                            _nj1;
                                                        vinfos[2].jointtype = 1;
                                                        vinfos[2].foffset = j2;
                                                        vinfos[2].indices[0] =
                                                            _ij2[0];
                                                        vinfos[2].indices[1] =
                                                            _ij2[1];
                                                        vinfos[2].maxsolutions =
                                                            _nj2;
                                                        std::vector<int> vfree(
                                                            0);
                                                        solutions.AddSolution(
                                                            vinfos, vfree);
                                                      }
                                                    }
                                                  }
                                                }
                                              }
                                            }
                                          }
                                        }
                                      } while (0);
                                      if (bgotonextstatement) {
                                        bool bgotonextstatement = true;
                                        do {
                                          if (1) {
                                            bgotonextstatement = false;
                                            continue; // branch miss [j0, j1]
                                          }
                                        } while (0);
                                        if (bgotonextstatement) {
                                        }
                                      }
                                    }
                                  }

                                } else {
                                  {
                                    IkReal j1array[1], cj1array[1], sj1array[1];
                                    bool j1valid[1] = {false};
                                    _nj1 = 1;
                                    IkReal x199 = cj2 * cj2;
                                    IkReal x200 = (cj2 * pz);
                                    IkReal x201 = (pz * sj2);
                                    IkReal x202 = (cj2 * sj2);
                                    CheckValue<IkReal> x203 =
                                        IKPowWithIntegerCheck(
                                            ((0.01439904) +
                                             (((-3479.768) * sj2)) +
                                             (((13000.0) * pz)) +
                                             (((-0.01274915) * cj2)) +
                                             (((-23200.0) * x200)) +
                                             (((0.085) * x201))),
                                            -1);
                                    if (!x203.valid) {
                                      continue;
                                    }
                                    CheckValue<IkReal> x204 =
                                        IKPowWithIntegerCheck(
                                            ((1439.904) +
                                             (((-1274.915) * cj2)) +
                                             (((8500.0) * x201)) +
                                             (((-347976800.0) * sj2)) +
                                             (((-2320000000.0) * x200)) +
                                             (((1300000000.0) * pz))),
                                            -1);
                                    if (!x204.valid) {
                                      continue;
                                    }
                                    if (IKabs(
                                            ((x203.value) *
                                             (((0.00724) +
                                               (((-0.03944) * x199)) +
                                               (((3015.9999999184) * sj2)) +
                                               (((-14999.0) * pz)) +
                                               (((-5382.39999992775) * x202)) +
                                               (((0.033322) * cj2)))))) <
                                            IKFAST_ATAN2_MAGTHRESH &&
                                        IKabs(
                                            ((x204.value) *
                                             (((55970000.992775) +
                                               (((-2210.0) * sj2)) +
                                               (((603200000.0) * cj2)) +
                                               (((-538239999.992775) * x199)) +
                                               (((3944.0) * x202)))))) <
                                            IKFAST_ATAN2_MAGTHRESH &&
                                        IKabs(
                                            IKsqr(
                                                ((x203.value) *
                                                 (((0.00724) +
                                                   (((-0.03944) * x199)) +
                                                   (((3015.9999999184) * sj2)) +
                                                   (((-14999.0) * pz)) +
                                                   (((-5382.39999992775) *
                                                     x202)) +
                                                   (((0.033322) * cj2)))))) +
                                            IKsqr(((x204.value) *
                                                   (((55970000.992775) +
                                                     (((-2210.0) * sj2)) +
                                                     (((603200000.0) * cj2)) +
                                                     (((-538239999.992775) *
                                                       x199)) +
                                                     (((3944.0) * x202)))))) -
                                            1) <= IKFAST_SINCOS_THRESH)
                                      continue;
                                    j1array[0] = IKatan2(
                                        ((x203.value) *
                                         (((0.00724) + (((-0.03944) * x199)) +
                                           (((3015.9999999184) * sj2)) +
                                           (((-14999.0) * pz)) +
                                           (((-5382.39999992775) * x202)) +
                                           (((0.033322) * cj2))))),
                                        ((x204.value) *
                                         (((55970000.992775) +
                                           (((-2210.0) * sj2)) +
                                           (((603200000.0) * cj2)) +
                                           (((-538239999.992775) * x199)) +
                                           (((3944.0) * x202))))));
                                    sj1array[0] = IKsin(j1array[0]);
                                    cj1array[0] = IKcos(j1array[0]);
                                    if (j1array[0] > IKPI) {
                                      j1array[0] -= IK2PI;
                                    } else if (j1array[0] < -IKPI) {
                                      j1array[0] += IK2PI;
                                    }
                                    j1valid[0] = true;
                                    for (int ij1 = 0; ij1 < 1; ++ij1) {
                                      if (!j1valid[ij1]) {
                                        continue;
                                      }
                                      _ij1[0] = ij1;
                                      _ij1[1] = -1;
                                      for (int iij1 = ij1 + 1; iij1 < 1;
                                           ++iij1) {
                                        if (j1valid[iij1] &&
                                            IKabs(cj1array[ij1] -
                                                  cj1array[iij1]) <
                                                IKFAST_SOLUTION_THRESH &&
                                            IKabs(sj1array[ij1] -
                                                  sj1array[iij1]) <
                                                IKFAST_SOLUTION_THRESH) {
                                          j1valid[iij1] = false;
                                          _ij1[1] = iij1;
                                          break;
                                        }
                                      }
                                      j1 = j1array[ij1];
                                      cj1 = cj1array[ij1];
                                      sj1 = sj1array[ij1];
                                      {
                                        IkReal evalcond[5];
                                        IkReal x205 = IKcos(j1);
                                        IkReal x206 = IKsin(j1);
                                        IkReal x207 = ((0.232) * cj2);
                                        IkReal x208 = ((1.0) * pz);
                                        IkReal x209 = ((0.232) * sj2);
                                        IkReal x210 = ((8.5e-7) * cj2);
                                        IkReal x211 = ((8.5e-7) * sj2);
                                        IkReal x212 = ((8.5e-7) * x205);
                                        evalcond[0] =
                                            ((-0.13) + x207 +
                                             (((-0.14999) * x206)) +
                                             (((-1.0) * x205 * x208)) +
                                             (((-1.0) * x211)));
                                        evalcond[1] =
                                            ((-9.6e-7) + (((0.14999) * x205)) +
                                             (((-1.0) * x206 * x208)) + x210 +
                                             x209);
                                        evalcond[2] =
                                            ((0.0172892498998009) +
                                             (((-0.26) * pz * x205)) +
                                             (((-1.0) * pz * x208)) +
                                             (((2.879808e-7) * x205)) +
                                             (((-1.92e-6) * pz * x206)) +
                                             (((-0.0389974) * x206)));
                                        evalcond[3] =
                                            ((-0.14999) +
                                             (((-1.0) * x206 * x211)) +
                                             (((-0.13) * x206)) +
                                             (((-1.0) * x205 * x210)) +
                                             (((9.6e-7) * x205)) +
                                             (((-1.0) * x205 * x209)) +
                                             ((x206 * x207)));
                                        evalcond[4] =
                                            ((((-9.6e-7) * x206)) +
                                             ((x205 * x207)) + ((x206 * x210)) +
                                             (((-0.13) * x205)) +
                                             (((-1.0) * x205 * x211)) +
                                             (((-1.0) * x208)) +
                                             ((x206 * x209)));
                                        if (IKabs(evalcond[0]) >
                                                IKFAST_EVALCOND_THRESH ||
                                            IKabs(evalcond[1]) >
                                                IKFAST_EVALCOND_THRESH ||
                                            IKabs(evalcond[2]) >
                                                IKFAST_EVALCOND_THRESH ||
                                            IKabs(evalcond[3]) >
                                                IKFAST_EVALCOND_THRESH ||
                                            IKabs(evalcond[4]) >
                                                IKFAST_EVALCOND_THRESH) {
                                          continue;
                                        }
                                      }

                                      {
                                        IkReal j0array[1], cj0array[1],
                                            sj0array[1];
                                        bool j0valid[1] = {false};
                                        _nj0 = 1;
                                        j0array[0] = 0;
                                        sj0array[0] = IKsin(j0array[0]);
                                        cj0array[0] = IKcos(j0array[0]);
                                        if (j0array[0] > IKPI) {
                                          j0array[0] -= IK2PI;
                                        } else if (j0array[0] < -IKPI) {
                                          j0array[0] += IK2PI;
                                        }
                                        j0valid[0] = true;
                                        for (int ij0 = 0; ij0 < 1; ++ij0) {
                                          if (!j0valid[ij0]) {
                                            continue;
                                          }
                                          _ij0[0] = ij0;
                                          _ij0[1] = -1;
                                          for (int iij0 = ij0 + 1; iij0 < 1;
                                               ++iij0) {
                                            if (j0valid[iij0] &&
                                                IKabs(cj0array[ij0] -
                                                      cj0array[iij0]) <
                                                    IKFAST_SOLUTION_THRESH &&
                                                IKabs(sj0array[ij0] -
                                                      sj0array[iij0]) <
                                                    IKFAST_SOLUTION_THRESH) {
                                              j0valid[iij0] = false;
                                              _ij0[1] = iij0;
                                              break;
                                            }
                                          }
                                          j0 = j0array[ij0];
                                          cj0 = cj0array[ij0];
                                          sj0 = sj0array[ij0];

                                          {
                                            std::vector<
                                                IkSingleDOFSolutionBase<IkReal>>
                                                vinfos(3);
                                            vinfos[0].jointtype = 1;
                                            vinfos[0].foffset = j0;
                                            vinfos[0].indices[0] = _ij0[0];
                                            vinfos[0].indices[1] = _ij0[1];
                                            vinfos[0].maxsolutions = _nj0;
                                            vinfos[1].jointtype = 1;
                                            vinfos[1].foffset = j1;
                                            vinfos[1].indices[0] = _ij1[0];
                                            vinfos[1].indices[1] = _ij1[1];
                                            vinfos[1].maxsolutions = _nj1;
                                            vinfos[2].jointtype = 1;
                                            vinfos[2].foffset = j2;
                                            vinfos[2].indices[0] = _ij2[0];
                                            vinfos[2].indices[1] = _ij2[1];
                                            vinfos[2].maxsolutions = _nj2;
                                            std::vector<int> vfree(0);
                                            solutions.AddSolution(vinfos,
                                                                  vfree);
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }

                            } else {
                              {
                                IkReal j1array[1], cj1array[1], sj1array[1];
                                bool j1valid[1] = {false};
                                _nj1 = 1;
                                IkReal x213 = pz * pz;
                                IkReal x214 = (cj2 * pz);
                                IkReal x215 = (pz * sj2);
                                CheckValue<IkReal> x216 = IKPowWithIntegerCheck(
                                    IKsign(((0.0043194240192) +
                                            (((0.192) * x213)))),
                                    -1);
                                if (!x216.valid) {
                                  continue;
                                }
                                CheckValue<IkReal> x217 = IKatan2WithCheck(
                                    IkReal(((-0.0037437504) +
                                            (((0.0221) * x215)) +
                                            (((5108.92498998009) * pz)) +
                                            (((-6032.0) * x214)) +
                                            (((-100000.0) * (pz * pz * pz))) +
                                            (((0.00668115456) * cj2)) +
                                            (((-2.4478368e-8) * sj2)))),
                                    IkReal(((-766.287659247114) +
                                            (((-0.02496) * pz)) +
                                            (((-0.003314779) * sj2)) +
                                            (((14999.0) * x213)) +
                                            (((904.73968) * cj2)) +
                                            (((-1.632e-7) * x215)) +
                                            (((0.044544) * x214)))),
                                    IKFAST_ATAN2_MAGTHRESH);
                                if (!x217.valid) {
                                  continue;
                                }
                                j1array[0] =
                                    ((-1.5707963267949) +
                                     (((1.5707963267949) * (x216.value))) +
                                     (x217.value));
                                sj1array[0] = IKsin(j1array[0]);
                                cj1array[0] = IKcos(j1array[0]);
                                if (j1array[0] > IKPI) {
                                  j1array[0] -= IK2PI;
                                } else if (j1array[0] < -IKPI) {
                                  j1array[0] += IK2PI;
                                }
                                j1valid[0] = true;
                                for (int ij1 = 0; ij1 < 1; ++ij1) {
                                  if (!j1valid[ij1]) {
                                    continue;
                                  }
                                  _ij1[0] = ij1;
                                  _ij1[1] = -1;
                                  for (int iij1 = ij1 + 1; iij1 < 1; ++iij1) {
                                    if (j1valid[iij1] &&
                                        IKabs(cj1array[ij1] - cj1array[iij1]) <
                                            IKFAST_SOLUTION_THRESH &&
                                        IKabs(sj1array[ij1] - sj1array[iij1]) <
                                            IKFAST_SOLUTION_THRESH) {
                                      j1valid[iij1] = false;
                                      _ij1[1] = iij1;
                                      break;
                                    }
                                  }
                                  j1 = j1array[ij1];
                                  cj1 = cj1array[ij1];
                                  sj1 = sj1array[ij1];
                                  {
                                    IkReal evalcond[5];
                                    IkReal x218 = IKcos(j1);
                                    IkReal x219 = IKsin(j1);
                                    IkReal x220 = ((0.232) * cj2);
                                    IkReal x221 = ((1.0) * pz);
                                    IkReal x222 = ((0.232) * sj2);
                                    IkReal x223 = ((8.5e-7) * cj2);
                                    IkReal x224 = ((8.5e-7) * sj2);
                                    IkReal x225 = ((8.5e-7) * x218);
                                    evalcond[0] = ((-0.13) + x220 +
                                                   (((-1.0) * x218 * x221)) +
                                                   (((-1.0) * x224)) +
                                                   (((-0.14999) * x219)));
                                    evalcond[1] =
                                        ((-9.6e-7) + (((-1.0) * x219 * x221)) +
                                         x223 + x222 + (((0.14999) * x218)));
                                    evalcond[2] = ((0.0172892498998009) +
                                                   (((-0.0389974) * x219)) +
                                                   (((-0.26) * pz * x218)) +
                                                   (((-1.92e-6) * pz * x219)) +
                                                   (((-1.0) * pz * x221)) +
                                                   (((2.879808e-7) * x218)));
                                    evalcond[3] =
                                        ((-0.14999) + (((-1.0) * x219 * x224)) +
                                         ((x219 * x220)) +
                                         (((-1.0) * x218 * x223)) +
                                         (((-1.0) * x218 * x222)) +
                                         (((-0.13) * x219)) +
                                         (((9.6e-7) * x218)));
                                    evalcond[4] =
                                        (((x219 * x223)) + ((x219 * x222)) +
                                         (((-9.6e-7) * x219)) +
                                         (((-1.0) * x218 * x224)) +
                                         ((x218 * x220)) + (((-0.13) * x218)) +
                                         (((-1.0) * x221)));
                                    if (IKabs(evalcond[0]) >
                                            IKFAST_EVALCOND_THRESH ||
                                        IKabs(evalcond[1]) >
                                            IKFAST_EVALCOND_THRESH ||
                                        IKabs(evalcond[2]) >
                                            IKFAST_EVALCOND_THRESH ||
                                        IKabs(evalcond[3]) >
                                            IKFAST_EVALCOND_THRESH ||
                                        IKabs(evalcond[4]) >
                                            IKFAST_EVALCOND_THRESH) {
                                      continue;
                                    }
                                  }

                                  {
                                    IkReal j0array[1], cj0array[1], sj0array[1];
                                    bool j0valid[1] = {false};
                                    _nj0 = 1;
                                    j0array[0] = 0;
                                    sj0array[0] = IKsin(j0array[0]);
                                    cj0array[0] = IKcos(j0array[0]);
                                    if (j0array[0] > IKPI) {
                                      j0array[0] -= IK2PI;
                                    } else if (j0array[0] < -IKPI) {
                                      j0array[0] += IK2PI;
                                    }
                                    j0valid[0] = true;
                                    for (int ij0 = 0; ij0 < 1; ++ij0) {
                                      if (!j0valid[ij0]) {
                                        continue;
                                      }
                                      _ij0[0] = ij0;
                                      _ij0[1] = -1;
                                      for (int iij0 = ij0 + 1; iij0 < 1;
                                           ++iij0) {
                                        if (j0valid[iij0] &&
                                            IKabs(cj0array[ij0] -
                                                  cj0array[iij0]) <
                                                IKFAST_SOLUTION_THRESH &&
                                            IKabs(sj0array[ij0] -
                                                  sj0array[iij0]) <
                                                IKFAST_SOLUTION_THRESH) {
                                          j0valid[iij0] = false;
                                          _ij0[1] = iij0;
                                          break;
                                        }
                                      }
                                      j0 = j0array[ij0];
                                      cj0 = cj0array[ij0];
                                      sj0 = sj0array[ij0];

                                      {
                                        std::vector<
                                            IkSingleDOFSolutionBase<IkReal>>
                                            vinfos(3);
                                        vinfos[0].jointtype = 1;
                                        vinfos[0].foffset = j0;
                                        vinfos[0].indices[0] = _ij0[0];
                                        vinfos[0].indices[1] = _ij0[1];
                                        vinfos[0].maxsolutions = _nj0;
                                        vinfos[1].jointtype = 1;
                                        vinfos[1].foffset = j1;
                                        vinfos[1].indices[0] = _ij1[0];
                                        vinfos[1].indices[1] = _ij1[1];
                                        vinfos[1].maxsolutions = _nj1;
                                        vinfos[2].jointtype = 1;
                                        vinfos[2].foffset = j2;
                                        vinfos[2].indices[0] = _ij2[0];
                                        vinfos[2].indices[1] = _ij2[1];
                                        vinfos[2].maxsolutions = _nj2;
                                        std::vector<int> vfree(0);
                                        solutions.AddSolution(vinfos, vfree);
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }

                        } else {
                          {
                            IkReal j1array[1], cj1array[1], sj1array[1];
                            bool j1valid[1] = {false};
                            _nj1 = 1;
                            IkReal x226 = (pz * sj2);
                            IkReal x227 = (cj2 * pz);
                            CheckValue<IkReal> x228 = IKPowWithIntegerCheck(
                                IKsign(((224970001.0) +
                                        (((10000000000.0) * (pz * pz))))),
                                -1);
                            if (!x228.valid) {
                              continue;
                            }
                            CheckValue<IkReal> x229 = IKatan2WithCheck(
                                IkReal(((-194987000.0) + (((-1274.915) * sj2)) +
                                        (((2320000000.0) * x226)) +
                                        (((8500.0) * x227)) +
                                        (((-9600.0) * pz)) +
                                        (((347976800.0) * cj2)))),
                                IkReal(((1439.904) + (((-1274.915) * cj2)) +
                                        (((2320000000.0) * x227)) +
                                        (((-347976800.0) * sj2)) +
                                        (((-8500.0) * x226)) +
                                        (((-1300000000.0) * pz)))),
                                IKFAST_ATAN2_MAGTHRESH);
                            if (!x229.valid) {
                              continue;
                            }
                            j1array[0] = ((-1.5707963267949) +
                                          (((1.5707963267949) * (x228.value))) +
                                          (x229.value));
                            sj1array[0] = IKsin(j1array[0]);
                            cj1array[0] = IKcos(j1array[0]);
                            if (j1array[0] > IKPI) {
                              j1array[0] -= IK2PI;
                            } else if (j1array[0] < -IKPI) {
                              j1array[0] += IK2PI;
                            }
                            j1valid[0] = true;
                            for (int ij1 = 0; ij1 < 1; ++ij1) {
                              if (!j1valid[ij1]) {
                                continue;
                              }
                              _ij1[0] = ij1;
                              _ij1[1] = -1;
                              for (int iij1 = ij1 + 1; iij1 < 1; ++iij1) {
                                if (j1valid[iij1] &&
                                    IKabs(cj1array[ij1] - cj1array[iij1]) <
                                        IKFAST_SOLUTION_THRESH &&
                                    IKabs(sj1array[ij1] - sj1array[iij1]) <
                                        IKFAST_SOLUTION_THRESH) {
                                  j1valid[iij1] = false;
                                  _ij1[1] = iij1;
                                  break;
                                }
                              }
                              j1 = j1array[ij1];
                              cj1 = cj1array[ij1];
                              sj1 = sj1array[ij1];
                              {
                                IkReal evalcond[5];
                                IkReal x230 = IKcos(j1);
                                IkReal x231 = IKsin(j1);
                                IkReal x232 = ((0.232) * cj2);
                                IkReal x233 = ((1.0) * pz);
                                IkReal x234 = ((0.232) * sj2);
                                IkReal x235 = ((8.5e-7) * cj2);
                                IkReal x236 = ((8.5e-7) * sj2);
                                IkReal x237 = ((8.5e-7) * x230);
                                evalcond[0] =
                                    ((-0.13) + (((-1.0) * x230 * x233)) +
                                     (((-1.0) * x236)) + x232 +
                                     (((-0.14999) * x231)));
                                evalcond[1] =
                                    ((-9.6e-7) + (((0.14999) * x230)) + x234 +
                                     x235 + (((-1.0) * x231 * x233)));
                                evalcond[2] = ((0.0172892498998009) +
                                               (((-0.26) * pz * x230)) +
                                               (((-1.0) * pz * x233)) +
                                               (((2.879808e-7) * x230)) +
                                               (((-1.92e-6) * pz * x231)) +
                                               (((-0.0389974) * x231)));
                                evalcond[3] =
                                    ((-0.14999) + (((-1.0) * x230 * x235)) +
                                     (((-1.0) * x230 * x234)) +
                                     (((-0.13) * x231)) +
                                     (((-1.0) * x231 * x236)) +
                                     (((9.6e-7) * x230)) + ((x231 * x232)));
                                evalcond[4] =
                                    ((((-9.6e-7) * x231)) +
                                     (((-1.0) * x230 * x236)) +
                                     (((-0.13) * x230)) + ((x230 * x232)) +
                                     (((-1.0) * x233)) + ((x231 * x234)) +
                                     ((x231 * x235)));
                                if (IKabs(evalcond[0]) >
                                        IKFAST_EVALCOND_THRESH ||
                                    IKabs(evalcond[1]) >
                                        IKFAST_EVALCOND_THRESH ||
                                    IKabs(evalcond[2]) >
                                        IKFAST_EVALCOND_THRESH ||
                                    IKabs(evalcond[3]) >
                                        IKFAST_EVALCOND_THRESH ||
                                    IKabs(evalcond[4]) >
                                        IKFAST_EVALCOND_THRESH) {
                                  continue;
                                }
                              }

                              {
                                IkReal j0array[1], cj0array[1], sj0array[1];
                                bool j0valid[1] = {false};
                                _nj0 = 1;
                                j0array[0] = 0;
                                sj0array[0] = IKsin(j0array[0]);
                                cj0array[0] = IKcos(j0array[0]);
                                if (j0array[0] > IKPI) {
                                  j0array[0] -= IK2PI;
                                } else if (j0array[0] < -IKPI) {
                                  j0array[0] += IK2PI;
                                }
                                j0valid[0] = true;
                                for (int ij0 = 0; ij0 < 1; ++ij0) {
                                  if (!j0valid[ij0]) {
                                    continue;
                                  }
                                  _ij0[0] = ij0;
                                  _ij0[1] = -1;
                                  for (int iij0 = ij0 + 1; iij0 < 1; ++iij0) {
                                    if (j0valid[iij0] &&
                                        IKabs(cj0array[ij0] - cj0array[iij0]) <
                                            IKFAST_SOLUTION_THRESH &&
                                        IKabs(sj0array[ij0] - sj0array[iij0]) <
                                            IKFAST_SOLUTION_THRESH) {
                                      j0valid[iij0] = false;
                                      _ij0[1] = iij0;
                                      break;
                                    }
                                  }
                                  j0 = j0array[ij0];
                                  cj0 = cj0array[ij0];
                                  sj0 = sj0array[ij0];

                                  {
                                    std::vector<IkSingleDOFSolutionBase<IkReal>>
                                        vinfos(3);
                                    vinfos[0].jointtype = 1;
                                    vinfos[0].foffset = j0;
                                    vinfos[0].indices[0] = _ij0[0];
                                    vinfos[0].indices[1] = _ij0[1];
                                    vinfos[0].maxsolutions = _nj0;
                                    vinfos[1].jointtype = 1;
                                    vinfos[1].foffset = j1;
                                    vinfos[1].indices[0] = _ij1[0];
                                    vinfos[1].indices[1] = _ij1[1];
                                    vinfos[1].maxsolutions = _nj1;
                                    vinfos[2].jointtype = 1;
                                    vinfos[2].foffset = j2;
                                    vinfos[2].indices[0] = _ij2[0];
                                    vinfos[2].indices[1] = _ij2[1];
                                    vinfos[2].maxsolutions = _nj2;
                                    std::vector<int> vfree(0);
                                    solutions.AddSolution(vinfos, vfree);
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              } while (0);
              if (bgotonextstatement) {
                bool bgotonextstatement = true;
                do {
                  if (1) {
                    bgotonextstatement = false;
                    continue; // branch miss [j0, j1, j2]
                  }
                } while (0);
                if (bgotonextstatement) {
                }
              }
            }
          }

        } else {
          {
            IkReal j0array[2], cj0array[2], sj0array[2];
            bool j0valid[2] = {false};
            _nj0 = 2;
            CheckValue<IkReal> x240 =
                IKatan2WithCheck(IkReal(((-1.0) * px)), IkReal(((-1.0) * py)),
                                 IKFAST_ATAN2_MAGTHRESH);
            if (!x240.valid) {
              continue;
            }
            IkReal x238 = ((1.0) * (x240.value));
            if ((((px * px) + (py * py))) < -0.00001)
              continue;
            CheckValue<IkReal> x241 = IKPowWithIntegerCheck(
                IKabs(IKsqrt(((px * px) + (py * py)))), -1);
            if (!x241.valid) {
              continue;
            }
            if ((((0.0535) * (x241.value))) < -1 - IKFAST_SINCOS_THRESH ||
                (((0.0535) * (x241.value))) > 1 + IKFAST_SINCOS_THRESH)
              continue;
            IkReal x239 = IKasin(((0.0535) * (x241.value)));
            j0array[0] = ((((-1.0) * x238)) + x239);
            sj0array[0] = IKsin(j0array[0]);
            cj0array[0] = IKcos(j0array[0]);
            j0array[1] =
                ((3.14159265358979) + (((-1.0) * x238)) + (((-1.0) * x239)));
            sj0array[1] = IKsin(j0array[1]);
            cj0array[1] = IKcos(j0array[1]);
            if (j0array[0] > IKPI) {
              j0array[0] -= IK2PI;
            } else if (j0array[0] < -IKPI) {
              j0array[0] += IK2PI;
            }
            j0valid[0] = true;
            if (j0array[1] > IKPI) {
              j0array[1] -= IK2PI;
            } else if (j0array[1] < -IKPI) {
              j0array[1] += IK2PI;
            }
            j0valid[1] = true;
            for (int ij0 = 0; ij0 < 2; ++ij0) {
              if (!j0valid[ij0]) {
                continue;
              }
              _ij0[0] = ij0;
              _ij0[1] = -1;
              for (int iij0 = ij0 + 1; iij0 < 2; ++iij0) {
                if (j0valid[iij0] &&
                    IKabs(cj0array[ij0] - cj0array[iij0]) <
                        IKFAST_SOLUTION_THRESH &&
                    IKabs(sj0array[ij0] - sj0array[iij0]) <
                        IKFAST_SOLUTION_THRESH) {
                  j0valid[iij0] = false;
                  _ij0[1] = iij0;
                  break;
                }
              }
              j0 = j0array[ij0];
              cj0 = cj0array[ij0];
              sj0 = sj0array[ij0];

              {
                IkReal j2array[2], cj2array[2], sj2array[2];
                bool j2valid[2] = {false};
                _nj2 = 2;
                if ((((-0.84697032327434) + (((16.5782493363067) * (px * px))) +
                      (((4.97314323590529) * cj0 * py)) +
                      (((-4.97314323590529) * px * sj0)) +
                      (((16.5782493363067) * (pz * pz))) +
                      (((16.5782493363067) * (py * py))))) <
                        -1 - IKFAST_SINCOS_THRESH ||
                    (((-0.84697032327434) + (((16.5782493363067) * (px * px))) +
                      (((4.97314323590529) * cj0 * py)) +
                      (((-4.97314323590529) * px * sj0)) +
                      (((16.5782493363067) * (pz * pz))) +
                      (((16.5782493363067) * (py * py))))) >
                        1 + IKFAST_SINCOS_THRESH)
                  continue;
                IkReal x242 = IKasin(((-0.84697032327434) +
                                      (((16.5782493363067) * (px * px))) +
                                      (((4.97314323590529) * cj0 * py)) +
                                      (((-4.97314323590529) * px * sj0)) +
                                      (((16.5782493363067) * (pz * pz))) +
                                      (((16.5782493363067) * (py * py)))));
                j2array[0] = ((1.57080004761718) + (((1.0) * x242)));
                sj2array[0] = IKsin(j2array[0]);
                cj2array[0] = IKcos(j2array[0]);
                j2array[1] = ((4.71239270120697) + (((-1.0) * x242)));
                sj2array[1] = IKsin(j2array[1]);
                cj2array[1] = IKcos(j2array[1]);
                if (j2array[0] > IKPI) {
                  j2array[0] -= IK2PI;
                } else if (j2array[0] < -IKPI) {
                  j2array[0] += IK2PI;
                }
                j2valid[0] = true;
                if (j2array[1] > IKPI) {
                  j2array[1] -= IK2PI;
                } else if (j2array[1] < -IKPI) {
                  j2array[1] += IK2PI;
                }
                j2valid[1] = true;
                for (int ij2 = 0; ij2 < 2; ++ij2) {
                  if (!j2valid[ij2]) {
                    continue;
                  }
                  _ij2[0] = ij2;
                  _ij2[1] = -1;
                  for (int iij2 = ij2 + 1; iij2 < 2; ++iij2) {
                    if (j2valid[iij2] &&
                        IKabs(cj2array[ij2] - cj2array[iij2]) <
                            IKFAST_SOLUTION_THRESH &&
                        IKabs(sj2array[ij2] - sj2array[iij2]) <
                            IKFAST_SOLUTION_THRESH) {
                      j2valid[iij2] = false;
                      _ij2[1] = iij2;
                      break;
                    }
                  }
                  j2 = j2array[ij2];
                  cj2 = cj2array[ij2];
                  sj2 = sj2array[ij2];

                  {
                    IkReal j1eval[3];
                    IkReal x243 = cj2 * cj2;
                    IkReal x244 = (px * sj0);
                    IkReal x245 = ((6.66711114074272) * sj2);
                    IkReal x246 = (pz * sj2);
                    IkReal x247 = ((0.232) * cj2);
                    IkReal x248 = (cj0 * py);
                    IkReal x249 = ((8.5e-7) * sj2);
                    IkReal x250 = ((1819729.15841448) * cj2);
                    IkReal x251 = (cj2 * pz);
                    IkReal x252 = (cj2 * sj2);
                    j1eval[0] =
                        ((152941.176470588) + (((-1019675.82152536) * x244)) +
                         sj2 + (((-7.52991375895648) * pz)) +
                         (((-1.0) * x248 * x250)) + ((x245 * x248)) +
                         ((x244 * x250)) + (((6.66711114074272) * x251)) +
                         (((1019675.82152536) * x248)) +
                         (((-1.0) * x244 * x245)) +
                         (((1819729.15841448) * x246)) +
                         (((-272941.176470588) * cj2)));
                    j1eval[1] = IKsign(
                        ((0.0194987) + (((-0.13) * x244)) + ((x248 * x249)) +
                         (((-9.6e-7) * pz)) + (((8.5e-7) * x251)) +
                         (((1.274915e-7) * sj2)) + ((x244 * x247)) +
                         (((0.13) * x248)) + (((0.232) * x246)) +
                         (((-0.03479768) * cj2)) + (((-1.0) * x247 * x248)) +
                         (((-1.0) * x244 * x249))));
                    j1eval[2] =
                        ((IKabs(((-0.0169000000007225) + (pz * pz) +
                                 (((-2.21e-7) * sj2)) +
                                 (((-0.0538239999992775) * x243)) +
                                 (((3.944e-7) * x252)) +
                                 (((0.06032) * cj2))))) +
                         (IKabs(((-7.24e-8) + (((0.0538239999992775) * x252)) +
                                 (((-1.0) * pz * x248)) +
                                 (((-0.030159999999184) * sj2)) +
                                 (((-0.14999) * pz)) + (((-3.3322e-7) * cj2)) +
                                 ((pz * x244)) + (((3.944e-7) * x243))))));
                    if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                        IKabs(j1eval[1]) < 0.0000010000000000 ||
                        IKabs(j1eval[2]) < 0.0000010000000000) {
                      {
                        IkReal j1eval[3];
                        IkReal x253 = cj2 * cj2;
                        IkReal x254 = (cj0 * py);
                        IkReal x255 = ((0.232) * sj2);
                        IkReal x256 = (cj2 * pz);
                        IkReal x257 = ((6.66711114074272) * cj2);
                        IkReal x258 = (pz * sj2);
                        IkReal x259 = ((8.5e-7) * cj2);
                        IkReal x260 = ((1819729.15841448) * sj2);
                        IkReal x261 = (px * sj0);
                        IkReal x262 = (cj2 * sj2);
                        j1eval[0] =
                            ((-1.12941176470588) +
                             (((-1019675.82152536) * pz)) + cj2 +
                             (((-7.52991375895648) * x254)) +
                             (((1819729.15841448) * x256)) + ((x254 * x260)) +
                             (((-6.66711114074272) * x258)) +
                             (((7.52991375895648) * x261)) + ((x254 * x257)) +
                             (((-1.0) * x257 * x261)) +
                             (((272941.176470588) * sj2)) +
                             (((-1.0) * x260 * x261)));
                        j1eval[1] =
                            ((IKabs(((-7.24e-8) + (((-1.0) * pz * x261)) +
                                     ((pz * x254)) + (((0.14999) * pz)) +
                                     (((-0.030159999999184) * sj2)) +
                                     (((0.0538239999992775) * x262)) +
                                     (((3.944e-7) * x253)) +
                                     (((-3.3322e-7) * cj2))))) +
                             (IKabs(((-0.0538240000009216) +
                                     (((4.4544e-7) * sj2)) +
                                     (((0.0538239999992775) * x253)) +
                                     (pz * pz) + (((-3.944e-7) * x262)) +
                                     (((1.632e-12) * cj2))))));
                        j1eval[2] = IKsign(
                            ((-1.439904e-7) + (((-1.0) * x255 * x261)) +
                             (((0.232) * x256)) + (((-1.0) * x259 * x261)) +
                             (((1.274915e-7) * cj2)) + (((-8.5e-7) * x258)) +
                             (((-9.6e-7) * x254)) + ((x254 * x259)) +
                             ((x254 * x255)) + (((9.6e-7) * x261)) +
                             (((0.03479768) * sj2)) + (((-0.13) * pz))));
                        if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                            IKabs(j1eval[1]) < 0.0000010000000000 ||
                            IKabs(j1eval[2]) < 0.0000010000000000) {
                          {
                            IkReal j1eval[2];
                            IkReal x263 = (cj0 * py);
                            IkReal x264 = ((0.232) * sj2);
                            IkReal x265 = ((1819729.15841448) * sj2);
                            IkReal x266 = (px * sj0);
                            IkReal x267 = (cj2 * pz);
                            IkReal x268 = ((6.66711114074272) * cj2);
                            IkReal x269 = (pz * sj2);
                            IkReal x270 = ((8.5e-7) * cj2);
                            j1eval[0] =
                                ((1.12941176470588) +
                                 (((-272941.176470588) * sj2)) +
                                 (((-1.0) * x263 * x265)) +
                                 (((-1.0) * x263 * x268)) + ((x266 * x268)) +
                                 (((6.66711114074272) * x269)) +
                                 (((-1819729.15841448) * x267)) +
                                 (((7.52991375895648) * x263)) +
                                 (((-7.52991375895648) * x266)) +
                                 ((x265 * x266)) + (((1019675.82152536) * pz)) +
                                 (((-1.0) * cj2)));
                            j1eval[1] = IKsign((
                                (1.439904e-7) + (((-1.0) * x263 * x264)) +
                                (((-1.0) * x263 * x270)) + ((x266 * x270)) +
                                (((-0.232) * x267)) + (((-1.274915e-7) * cj2)) +
                                (((-9.6e-7) * x266)) + (((-0.03479768) * sj2)) +
                                ((x264 * x266)) + (((0.13) * pz)) +
                                (((9.6e-7) * x263)) + (((8.5e-7) * x269))));
                            if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                                IKabs(j1eval[1]) < 0.0000010000000000) {
                              continue; // no branches [j1]

                            } else {
                              {
                                IkReal j1array[1], cj1array[1], sj1array[1];
                                bool j1valid[1] = {false};
                                _nj1 = 1;
                                IkReal x271 = cj2 * cj2;
                                IkReal x272 = cj0 * cj0;
                                IkReal x273 = px * px;
                                IkReal x274 = (cj0 * py);
                                IkReal x275 = ((0.232) * sj2);
                                IkReal x276 = (px * sj0);
                                IkReal x277 = ((8.5e-7) * cj2);
                                IkReal x278 = (cj2 * sj2);
                                CheckValue<IkReal> x279 = IKPowWithIntegerCheck(
                                    IKsign(((1.439904e-7) +
                                            (((-1.0) * x274 * x275)) +
                                            (((-1.0) * x274 * x277)) +
                                            ((x275 * x276)) +
                                            (((-1.274915e-7) * cj2)) +
                                            (((-9.6e-7) * x276)) +
                                            (((-0.03479768) * sj2)) +
                                            (((-0.232) * cj2 * pz)) +
                                            ((x276 * x277)) + (((0.13) * pz)) +
                                            (((9.6e-7) * x274)) +
                                            (((8.5e-7) * pz * sj2)))),
                                    -1);
                                if (!x279.valid) {
                                  continue;
                                }
                                CheckValue<IkReal> x280 = IKatan2WithCheck(
                                    IkReal(((7.24e-8) +
                                            (((0.030159999999184) * sj2)) +
                                            (((-1.0) * pz * x274)) +
                                            (((-0.14999) * pz)) +
                                            (((-3.944e-7) * x271)) +
                                            (((-0.0538239999992775) * x278)) +
                                            (((3.3322e-7) * cj2)) +
                                            ((pz * x276)))),
                                    IkReal(((0.0055970000992775) +
                                            (((-2.0) * x274 * x276)) +
                                            ((x272 * (py * py))) +
                                            (((0.29998) * x274)) +
                                            (((-0.29998) * x276)) + x273 +
                                            (((-2.21e-7) * sj2)) +
                                            (((-0.0538239999992775) * x271)) +
                                            (((3.944e-7) * x278)) +
                                            (((-1.0) * x272 * x273)) +
                                            (((0.06032) * cj2)))),
                                    IKFAST_ATAN2_MAGTHRESH);
                                if (!x280.valid) {
                                  continue;
                                }
                                j1array[0] =
                                    ((-1.5707963267949) +
                                     (((1.5707963267949) * (x279.value))) +
                                     (x280.value));
                                sj1array[0] = IKsin(j1array[0]);
                                cj1array[0] = IKcos(j1array[0]);
                                if (j1array[0] > IKPI) {
                                  j1array[0] -= IK2PI;
                                } else if (j1array[0] < -IKPI) {
                                  j1array[0] += IK2PI;
                                }
                                j1valid[0] = true;
                                for (int ij1 = 0; ij1 < 1; ++ij1) {
                                  if (!j1valid[ij1]) {
                                    continue;
                                  }
                                  _ij1[0] = ij1;
                                  _ij1[1] = -1;
                                  for (int iij1 = ij1 + 1; iij1 < 1; ++iij1) {
                                    if (j1valid[iij1] &&
                                        IKabs(cj1array[ij1] - cj1array[iij1]) <
                                            IKFAST_SOLUTION_THRESH &&
                                        IKabs(sj1array[ij1] - sj1array[iij1]) <
                                            IKFAST_SOLUTION_THRESH) {
                                      j1valid[iij1] = false;
                                      _ij1[1] = iij1;
                                      break;
                                    }
                                  }
                                  j1 = j1array[ij1];
                                  cj1 = cj1array[ij1];
                                  sj1 = sj1array[ij1];
                                  {
                                    IkReal evalcond[5];
                                    IkReal x281 = IKsin(j1);
                                    IkReal x282 = IKcos(j1);
                                    IkReal x283 = ((0.232) * cj2);
                                    IkReal x284 = (px * sj0);
                                    IkReal x285 = ((1.0) * pz);
                                    IkReal x286 = ((0.232) * sj2);
                                    IkReal x287 = (cj0 * py);
                                    IkReal x288 = ((8.5e-7) * cj2);
                                    IkReal x289 = ((8.5e-7) * sj2);
                                    IkReal x290 = ((8.5e-7) * x282);
                                    IkReal x291 = ((0.26) * x281);
                                    IkReal x292 = ((1.92e-6) * x282);
                                    evalcond[0] =
                                        ((-0.13) + (((-1.0) * x281 * x287)) +
                                         x283 + (((-1.0) * x282 * x285)) +
                                         (((-1.0) * x289)) + ((x281 * x284)) +
                                         (((-0.14999) * x281)));
                                    evalcond[1] =
                                        ((-9.6e-7) + (((-1.0) * x282 * x284)) +
                                         (((-1.0) * x281 * x285)) + x288 +
                                         x286 + (((0.14999) * x282)) +
                                         ((x282 * x287)));
                                    evalcond[2] =
                                        ((((-9.6e-7) * x281)) +
                                         (((-1.0) * x282 * x289)) +
                                         (((-1.0) * x285)) + ((x281 * x288)) +
                                         ((x281 * x286)) + ((x282 * x283)) +
                                         (((-0.13) * x282)));
                                    evalcond[3] =
                                        ((-0.14999) + (((-1.0) * x281 * x289)) +
                                         x284 + (((-1.0) * x282 * x286)) +
                                         (((-1.0) * x282 * x288)) +
                                         (((-1.0) * x287)) +
                                         (((9.6e-7) * x282)) + ((x281 * x283)) +
                                         (((-0.13) * x281)));
                                    evalcond[4] = ((0.0172892498998009) +
                                                   (((2.879808e-7) * x282)) +
                                                   (((-1.0) * x284 * x292)) +
                                                   (((-0.29998) * x287)) +
                                                   (((-1.0) * (px * px))) +
                                                   (((0.29998) * x284)) +
                                                   ((x287 * x292)) +
                                                   (((-1.92e-6) * pz * x281)) +
                                                   ((x284 * x291)) +
                                                   (((-0.0389974) * x281)) +
                                                   (((-0.26) * pz * x282)) +
                                                   (((-1.0) * x287 * x291)) +
                                                   (((-1.0) * pz * x285)) +
                                                   (((-1.0) * (py * py))));
                                    if (IKabs(evalcond[0]) >
                                            IKFAST_EVALCOND_THRESH ||
                                        IKabs(evalcond[1]) >
                                            IKFAST_EVALCOND_THRESH ||
                                        IKabs(evalcond[2]) >
                                            IKFAST_EVALCOND_THRESH ||
                                        IKabs(evalcond[3]) >
                                            IKFAST_EVALCOND_THRESH ||
                                        IKabs(evalcond[4]) >
                                            IKFAST_EVALCOND_THRESH) {
                                      continue;
                                    }
                                  }

                                  {
                                    std::vector<IkSingleDOFSolutionBase<IkReal>>
                                        vinfos(3);
                                    vinfos[0].jointtype = 1;
                                    vinfos[0].foffset = j0;
                                    vinfos[0].indices[0] = _ij0[0];
                                    vinfos[0].indices[1] = _ij0[1];
                                    vinfos[0].maxsolutions = _nj0;
                                    vinfos[1].jointtype = 1;
                                    vinfos[1].foffset = j1;
                                    vinfos[1].indices[0] = _ij1[0];
                                    vinfos[1].indices[1] = _ij1[1];
                                    vinfos[1].maxsolutions = _nj1;
                                    vinfos[2].jointtype = 1;
                                    vinfos[2].foffset = j2;
                                    vinfos[2].indices[0] = _ij2[0];
                                    vinfos[2].indices[1] = _ij2[1];
                                    vinfos[2].maxsolutions = _nj2;
                                    std::vector<int> vfree(0);
                                    solutions.AddSolution(vinfos, vfree);
                                  }
                                }
                              }
                            }
                          }

                        } else {
                          {
                            IkReal j1array[1], cj1array[1], sj1array[1];
                            bool j1valid[1] = {false};
                            _nj1 = 1;
                            IkReal x293 = cj2 * cj2;
                            IkReal x294 = (cj0 * py);
                            IkReal x295 = ((0.232) * sj2);
                            IkReal x296 = ((8.5e-7) * cj2);
                            IkReal x297 = (px * sj0);
                            IkReal x298 = (cj2 * sj2);
                            CheckValue<IkReal> x299 = IKatan2WithCheck(
                                IkReal(((-7.24e-8) + (((-1.0) * pz * x297)) +
                                        (((0.14999) * pz)) +
                                        (((-0.030159999999184) * sj2)) +
                                        ((pz * x294)) + (((3.944e-7) * x293)) +
                                        (((-3.3322e-7) * cj2)) +
                                        (((0.0538239999992775) * x298)))),
                                IkReal(((-0.0538240000009216) +
                                        (((4.4544e-7) * sj2)) + (pz * pz) +
                                        (((1.632e-12) * cj2)) +
                                        (((0.0538239999992775) * x293)) +
                                        (((-3.944e-7) * x298)))),
                                IKFAST_ATAN2_MAGTHRESH);
                            if (!x299.valid) {
                              continue;
                            }
                            CheckValue<IkReal> x300 = IKPowWithIntegerCheck(
                                IKsign(
                                    ((-1.439904e-7) + (((-1.0) * x296 * x297)) +
                                     (((1.274915e-7) * cj2)) +
                                     (((9.6e-7) * x297)) +
                                     (((-8.5e-7) * pz * sj2)) +
                                     (((-1.0) * x295 * x297)) +
                                     (((0.232) * cj2 * pz)) + ((x294 * x295)) +
                                     ((x294 * x296)) + (((0.03479768) * sj2)) +
                                     (((-0.13) * pz)) + (((-9.6e-7) * x294)))),
                                -1);
                            if (!x300.valid) {
                              continue;
                            }
                            j1array[0] = ((-1.5707963267949) + (x299.value) +
                                          (((1.5707963267949) * (x300.value))));
                            sj1array[0] = IKsin(j1array[0]);
                            cj1array[0] = IKcos(j1array[0]);
                            if (j1array[0] > IKPI) {
                              j1array[0] -= IK2PI;
                            } else if (j1array[0] < -IKPI) {
                              j1array[0] += IK2PI;
                            }
                            j1valid[0] = true;
                            for (int ij1 = 0; ij1 < 1; ++ij1) {
                              if (!j1valid[ij1]) {
                                continue;
                              }
                              _ij1[0] = ij1;
                              _ij1[1] = -1;
                              for (int iij1 = ij1 + 1; iij1 < 1; ++iij1) {
                                if (j1valid[iij1] &&
                                    IKabs(cj1array[ij1] - cj1array[iij1]) <
                                        IKFAST_SOLUTION_THRESH &&
                                    IKabs(sj1array[ij1] - sj1array[iij1]) <
                                        IKFAST_SOLUTION_THRESH) {
                                  j1valid[iij1] = false;
                                  _ij1[1] = iij1;
                                  break;
                                }
                              }
                              j1 = j1array[ij1];
                              cj1 = cj1array[ij1];
                              sj1 = sj1array[ij1];
                              {
                                IkReal evalcond[5];
                                IkReal x301 = IKsin(j1);
                                IkReal x302 = IKcos(j1);
                                IkReal x303 = ((0.232) * cj2);
                                IkReal x304 = (px * sj0);
                                IkReal x305 = ((1.0) * pz);
                                IkReal x306 = ((0.232) * sj2);
                                IkReal x307 = (cj0 * py);
                                IkReal x308 = ((8.5e-7) * cj2);
                                IkReal x309 = ((8.5e-7) * sj2);
                                IkReal x310 = ((8.5e-7) * x302);
                                IkReal x311 = ((0.26) * x301);
                                IkReal x312 = ((1.92e-6) * x302);
                                evalcond[0] =
                                    ((-0.13) + ((x301 * x304)) +
                                     (((-1.0) * x302 * x305)) + x303 +
                                     (((-1.0) * x301 * x307)) +
                                     (((-0.14999) * x301)) + (((-1.0) * x309)));
                                evalcond[1] =
                                    ((-9.6e-7) + (((-1.0) * x301 * x305)) +
                                     (((-1.0) * x302 * x304)) +
                                     (((0.14999) * x302)) + x308 + x306 +
                                     ((x302 * x307)));
                                evalcond[2] =
                                    (((x301 * x306)) + ((x301 * x308)) +
                                     (((-0.13) * x302)) +
                                     (((-1.0) * x302 * x309)) +
                                     (((-9.6e-7) * x301)) + ((x302 * x303)) +
                                     (((-1.0) * x305)));
                                evalcond[3] = ((-0.14999) + ((x301 * x303)) +
                                               (((-1.0) * x301 * x309)) +
                                               (((-0.13) * x301)) +
                                               (((-1.0) * x302 * x308)) +
                                               (((-1.0) * x302 * x306)) +
                                               (((-1.0) * x307)) + x304 +
                                               (((9.6e-7) * x302)));
                                evalcond[4] =
                                    ((0.0172892498998009) +
                                     (((-1.0) * (px * px))) + ((x304 * x311)) +
                                     (((-0.26) * pz * x302)) +
                                     (((2.879808e-7) * x302)) +
                                     (((-0.29998) * x307)) + ((x307 * x312)) +
                                     (((-1.0) * x304 * x312)) +
                                     (((-1.0) * pz * x305)) +
                                     (((-1.0) * x307 * x311)) +
                                     (((-1.0) * (py * py))) +
                                     (((0.29998) * x304)) +
                                     (((-1.92e-6) * pz * x301)) +
                                     (((-0.0389974) * x301)));
                                if (IKabs(evalcond[0]) >
                                        IKFAST_EVALCOND_THRESH ||
                                    IKabs(evalcond[1]) >
                                        IKFAST_EVALCOND_THRESH ||
                                    IKabs(evalcond[2]) >
                                        IKFAST_EVALCOND_THRESH ||
                                    IKabs(evalcond[3]) >
                                        IKFAST_EVALCOND_THRESH ||
                                    IKabs(evalcond[4]) >
                                        IKFAST_EVALCOND_THRESH) {
                                  continue;
                                }
                              }

                              {
                                std::vector<IkSingleDOFSolutionBase<IkReal>>
                                    vinfos(3);
                                vinfos[0].jointtype = 1;
                                vinfos[0].foffset = j0;
                                vinfos[0].indices[0] = _ij0[0];
                                vinfos[0].indices[1] = _ij0[1];
                                vinfos[0].maxsolutions = _nj0;
                                vinfos[1].jointtype = 1;
                                vinfos[1].foffset = j1;
                                vinfos[1].indices[0] = _ij1[0];
                                vinfos[1].indices[1] = _ij1[1];
                                vinfos[1].maxsolutions = _nj1;
                                vinfos[2].jointtype = 1;
                                vinfos[2].foffset = j2;
                                vinfos[2].indices[0] = _ij2[0];
                                vinfos[2].indices[1] = _ij2[1];
                                vinfos[2].maxsolutions = _nj2;
                                std::vector<int> vfree(0);
                                solutions.AddSolution(vinfos, vfree);
                              }
                            }
                          }
                        }
                      }

                    } else {
                      {
                        IkReal j1array[1], cj1array[1], sj1array[1];
                        bool j1valid[1] = {false};
                        _nj1 = 1;
                        IkReal x313 = cj2 * cj2;
                        IkReal x314 = (px * sj0);
                        IkReal x315 = ((0.232) * cj2);
                        IkReal x316 = ((8.5e-7) * sj2);
                        IkReal x317 = (cj0 * py);
                        IkReal x318 = (cj2 * sj2);
                        CheckValue<IkReal> x319 = IKatan2WithCheck(
                            IkReal(((-0.0169000000007225) +
                                    (((3.944e-7) * x318)) + (pz * pz) +
                                    (((-2.21e-7) * sj2)) +
                                    (((-0.0538239999992775) * x313)) +
                                    (((0.06032) * cj2)))),
                            IkReal(
                                ((-7.24e-8) + (((0.0538239999992775) * x318)) +
                                 (((3.944e-7) * x313)) +
                                 (((-0.030159999999184) * sj2)) +
                                 (((-0.14999) * pz)) + (((-1.0) * pz * x317)) +
                                 ((pz * x314)) + (((-3.3322e-7) * cj2)))),
                            IKFAST_ATAN2_MAGTHRESH);
                        if (!x319.valid) {
                          continue;
                        }
                        CheckValue<IkReal> x320 = IKPowWithIntegerCheck(
                            IKsign(((0.0194987) + (((8.5e-7) * cj2 * pz)) +
                                    (((-9.6e-7) * pz)) + (((-0.13) * x314)) +
                                    (((1.274915e-7) * sj2)) +
                                    (((-1.0) * x315 * x317)) +
                                    (((0.13) * x317)) + (((0.232) * pz * sj2)) +
                                    (((-1.0) * x314 * x316)) +
                                    (((-0.03479768) * cj2)) + ((x314 * x315)) +
                                    ((x316 * x317)))),
                            -1);
                        if (!x320.valid) {
                          continue;
                        }
                        j1array[0] = ((-1.5707963267949) + (x319.value) +
                                      (((1.5707963267949) * (x320.value))));
                        sj1array[0] = IKsin(j1array[0]);
                        cj1array[0] = IKcos(j1array[0]);
                        if (j1array[0] > IKPI) {
                          j1array[0] -= IK2PI;
                        } else if (j1array[0] < -IKPI) {
                          j1array[0] += IK2PI;
                        }
                        j1valid[0] = true;
                        for (int ij1 = 0; ij1 < 1; ++ij1) {
                          if (!j1valid[ij1]) {
                            continue;
                          }
                          _ij1[0] = ij1;
                          _ij1[1] = -1;
                          for (int iij1 = ij1 + 1; iij1 < 1; ++iij1) {
                            if (j1valid[iij1] &&
                                IKabs(cj1array[ij1] - cj1array[iij1]) <
                                    IKFAST_SOLUTION_THRESH &&
                                IKabs(sj1array[ij1] - sj1array[iij1]) <
                                    IKFAST_SOLUTION_THRESH) {
                              j1valid[iij1] = false;
                              _ij1[1] = iij1;
                              break;
                            }
                          }
                          j1 = j1array[ij1];
                          cj1 = cj1array[ij1];
                          sj1 = sj1array[ij1];
                          {
                            IkReal evalcond[5];
                            IkReal x321 = IKsin(j1);
                            IkReal x322 = IKcos(j1);
                            IkReal x323 = ((0.232) * cj2);
                            IkReal x324 = (px * sj0);
                            IkReal x325 = ((1.0) * pz);
                            IkReal x326 = ((0.232) * sj2);
                            IkReal x327 = (cj0 * py);
                            IkReal x328 = ((8.5e-7) * cj2);
                            IkReal x329 = ((8.5e-7) * sj2);
                            IkReal x330 = ((8.5e-7) * x322);
                            IkReal x331 = ((0.26) * x321);
                            IkReal x332 = ((1.92e-6) * x322);
                            evalcond[0] =
                                ((-0.13) + (((-1.0) * x321 * x327)) +
                                 (((-1.0) * x329)) + (((-1.0) * x322 * x325)) +
                                 x323 + (((-0.14999) * x321)) +
                                 ((x321 * x324)));
                            evalcond[1] =
                                ((-9.6e-7) + (((-1.0) * x322 * x324)) +
                                 ((x322 * x327)) + (((-1.0) * x321 * x325)) +
                                 x326 + x328 + (((0.14999) * x322)));
                            evalcond[2] =
                                (((x322 * x323)) + (((-1.0) * x325)) +
                                 (((-1.0) * x322 * x329)) + (((-0.13) * x322)) +
                                 (((-9.6e-7) * x321)) + ((x321 * x328)) +
                                 ((x321 * x326)));
                            evalcond[3] =
                                ((-0.14999) + (((-1.0) * x321 * x329)) +
                                 (((-1.0) * x327)) + (((-1.0) * x322 * x326)) +
                                 (((-1.0) * x322 * x328)) + (((-0.13) * x321)) +
                                 x324 + (((9.6e-7) * x322)) + ((x321 * x323)));
                            evalcond[4] =
                                ((0.0172892498998009) + (((-1.0) * (px * px))) +
                                 (((-0.29998) * x327)) + ((x324 * x331)) +
                                 (((0.29998) * x324)) +
                                 (((-0.26) * pz * x322)) +
                                 (((-1.0) * x324 * x332)) +
                                 (((-1.92e-6) * pz * x321)) + ((x327 * x332)) +
                                 (((2.879808e-7) * x322)) +
                                 (((-1.0) * x327 * x331)) +
                                 (((-1.0) * (py * py))) +
                                 (((-1.0) * pz * x325)) +
                                 (((-0.0389974) * x321)));
                            if (IKabs(evalcond[0]) > IKFAST_EVALCOND_THRESH ||
                                IKabs(evalcond[1]) > IKFAST_EVALCOND_THRESH ||
                                IKabs(evalcond[2]) > IKFAST_EVALCOND_THRESH ||
                                IKabs(evalcond[3]) > IKFAST_EVALCOND_THRESH ||
                                IKabs(evalcond[4]) > IKFAST_EVALCOND_THRESH) {
                              continue;
                            }
                          }

                          {
                            std::vector<IkSingleDOFSolutionBase<IkReal>> vinfos(
                                3);
                            vinfos[0].jointtype = 1;
                            vinfos[0].foffset = j0;
                            vinfos[0].indices[0] = _ij0[0];
                            vinfos[0].indices[1] = _ij0[1];
                            vinfos[0].maxsolutions = _nj0;
                            vinfos[1].jointtype = 1;
                            vinfos[1].foffset = j1;
                            vinfos[1].indices[0] = _ij1[0];
                            vinfos[1].indices[1] = _ij1[1];
                            vinfos[1].maxsolutions = _nj1;
                            vinfos[2].jointtype = 1;
                            vinfos[2].foffset = j2;
                            vinfos[2].indices[0] = _ij2[0];
                            vinfos[2].indices[1] = _ij2[1];
                            vinfos[2].maxsolutions = _nj2;
                            std::vector<int> vfree(0);
                            solutions.AddSolution(vinfos, vfree);
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    return solutions.GetNumSolutions() > 0;
  }
};

/// solves the inverse kinematics equations.
/// \param pfree is an array specifying the free joints of the chain.
IKFAST_API bool el_mini_ComputeIk_lb(const IkReal *eetrans, const IkReal *eerot,
                                     const IkReal *pfree,
                                     IkSolutionListBase<IkReal> &solutions) {
  IKSolver_lb solver;
  return solver.ComputeIk(eetrans, eerot, pfree, solutions);
}

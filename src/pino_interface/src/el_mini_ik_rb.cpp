
#include "pino_interface/ik_fast_comm.h"
class IKSolver_rb {
public:
  IkReal j0, cj0, sj0, htj0, j0mul, j1, cj1, sj1, htj1, j1mul, j2, cj2, sj2,
      htj2, j2mul, new_px, px, npx, new_py, py, npy, new_pz, pz, npz, pp;
  unsigned char _ij0[2], _nj0, _ij1[2], _nj1, _ij2[2], _nj2;

  IkReal j100, cj100, sj100;
  unsigned char _ij100[2], _nj100;
  bool ComputeIk(const IkReal *eetrans, const IkReal *eerot,
                 const IkReal *pfree, IkSolutionListBase<IkReal> &solutions) {

    // printf("IKSolver_rb");

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
      new_py = py;
      new_pz = ((0.045056) + pz);
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
                  IkReal x16 = IKasin(((0.84697032327434) +
                                       (((-16.5782493363067) * (pz * pz)))));
                  j2array[0] = ((1.57078527838641) + (((-1.0) * x16)));
                  sj2array[0] = IKsin(j2array[0]);
                  cj2array[0] = IKcos(j2array[0]);
                  j2array[1] = ((4.7123779319762) + x16);
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
                      IkReal x17 = pz * pz;
                      IkReal x18 = (pz * sj2);
                      IkReal x19 = (cj2 * pz);
                      j1eval[0] = ((1.0) + (((44.4503709630156) * x17)));
                      j1eval[1] =
                          IKsign(((224970001.0) + (((10000000000.0) * x17))));
                      j1eval[2] =
                          ((IKabs(((-1439.904) + (((-1274.915) * cj2)) +
                                   (((-347976800.0) * sj2)) +
                                   (((8500.0) * x18)) +
                                   (((1300000000.0) * pz)) +
                                   (((-2320000000.0) * x19))))) +
                           (IKabs(((-194987000.0) + (((-1274.915) * sj2)) +
                                   (((-8500.0) * x19)) + (((-9600.0) * pz)) +
                                   (((347976800.0) * cj2)) +
                                   (((-2320000000.0) * x18))))));
                      if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                          IKabs(j1eval[1]) < 0.0000010000000000 ||
                          IKabs(j1eval[2]) < 0.0000010000000000) {
                        {
                          IkReal j1eval[3];
                          px = 0;
                          py = 0;
                          pp = pz * pz;
                          IkReal x20 = cj2 * cj2;
                          IkReal x21 = (cj2 * sj2);
                          IkReal x22 = (pz * sj2);
                          IkReal x23 = (cj2 * pz);
                          j1eval[0] = ((-152941.176470588) + (((-1.0) * sj2)) +
                                       (((6.66711114074272) * x23)) +
                                       (((1819729.15841448) * x22)) +
                                       (((272941.176470588) * cj2)) +
                                       (((7.52991375895648) * pz)));
                          j1eval[1] =
                              ((IKabs(((0.0322) + (((-0.011222) * cj2)) +
                                       (((-0.03944) * x20)) +
                                       (((-5382.39999992775) * x21)) +
                                       (((-14999.0) * pz)) +
                                       (((3016.0000000816) * sj2))))) +
                               (IKabs(((1690.00000007225) + (((0.0221) * sj2)) +
                                       (((5382.39999992775) * x20)) +
                                       (((-0.03944) * x21)) +
                                       (((-6032.0) * cj2)) +
                                       (((-100000.0) * (pz * pz)))))));
                          j1eval[2] = IKsign(
                              ((-1949.87) + (((3479.768) * cj2)) +
                               (((23200.0) * x22)) + (((0.085) * x23)) +
                               (((-0.01274915) * sj2)) + (((0.096) * pz))));
                          if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                              IKabs(j1eval[1]) < 0.0000010000000000 ||
                              IKabs(j1eval[2]) < 0.0000010000000000) {
                            {
                              IkReal j1eval[3];
                              px = 0;
                              py = 0;
                              pp = pz * pz;
                              IkReal x24 = cj2 * cj2;
                              IkReal x25 = (pz * sj2);
                              IkReal x26 = (cj2 * sj2);
                              IkReal x27 = (cj2 * pz);
                              j1eval[0] = ((-1.12941176470588) +
                                           (((-272941.176470588) * sj2)) +
                                           (((-1019675.82152536) * pz)) +
                                           (((1819729.15841448) * x27)) +
                                           (((-1.0) * cj2)) +
                                           (((-6.66711114074272) * x25)));
                              j1eval[1] = IKsign(
                                  ((-0.01439904) + (((23200.0) * x27)) +
                                   (((-3479.768) * sj2)) +
                                   (((-0.01274915) * cj2)) +
                                   (((-13000.0) * pz)) + (((-0.085) * x25))));
                              j1eval[2] =
                                  ((IKabs(((5382.40000009216) +
                                           (((1.632e-7) * cj2)) +
                                           (((0.044544) * sj2)) +
                                           (((-100000.0) * (pz * pz))) +
                                           (((-5382.39999992775) * x24)) +
                                           (((0.03944) * x26))))) +
                                   (IKabs(((0.0322) + (((-0.011222) * cj2)) +
                                           (((-0.03944) * x24)) +
                                           (((-5382.39999992775) * x26)) +
                                           (((3016.0000000816) * sj2)) +
                                           (((14999.0) * pz))))));
                              if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                                  IKabs(j1eval[1]) < 0.0000010000000000 ||
                                  IKabs(j1eval[2]) < 0.0000010000000000) {
                                continue; // no branches [j0, j1]

                              } else {
                                {
                                  IkReal j1array[1], cj1array[1], sj1array[1];
                                  bool j1valid[1] = {false};
                                  _nj1 = 1;
                                  IkReal x28 = cj2 * cj2;
                                  IkReal x29 = (cj2 * sj2);
                                  CheckValue<IkReal> x30 =
                                      IKPowWithIntegerCheck(
                                          IKsign(((-0.01439904) +
                                                  (((23200.0) * cj2 * pz)) +
                                                  (((-3479.768) * sj2)) +
                                                  (((-0.085) * pz * sj2)) +
                                                  (((-0.01274915) * cj2)) +
                                                  (((-13000.0) * pz)))),
                                          -1);
                                  if (!x30.valid) {
                                    continue;
                                  }
                                  CheckValue<IkReal> x31 = IKatan2WithCheck(
                                      IkReal(((0.0322) + (((-0.011222) * cj2)) +
                                              (((-0.03944) * x28)) +
                                              (((-5382.39999992775) * x29)) +
                                              (((3016.0000000816) * sj2)) +
                                              (((14999.0) * pz)))),
                                      IkReal(((5382.40000009216) +
                                              (((1.632e-7) * cj2)) +
                                              (((0.044544) * sj2)) +
                                              (((-100000.0) * (pz * pz))) +
                                              (((-5382.39999992775) * x28)) +
                                              (((0.03944) * x29)))),
                                      IKFAST_ATAN2_MAGTHRESH);
                                  if (!x31.valid) {
                                    continue;
                                  }
                                  j1array[0] =
                                      ((-1.5707963267949) +
                                       (((1.5707963267949) * (x30.value))) +
                                       (x31.value));
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
                                      IkReal x32 = IKsin(j1);
                                      IkReal x33 = IKcos(j1);
                                      IkReal x34 = ((0.232) * cj2);
                                      IkReal x35 = ((0.232) * sj2);
                                      IkReal x36 = ((8.5e-7) * cj2);
                                      IkReal x37 = ((8.5e-7) * sj2);
                                      IkReal x38 = ((8.5e-7) * x33);
                                      IkReal x39 = (pz * x33);
                                      IkReal x40 = (pz * x32);
                                      evalcond[0] =
                                          ((-0.13) + (((-1.0) * x37)) + x39 +
                                           x34 + (((-0.14999) * x32)));
                                      evalcond[1] =
                                          ((9.6e-7) + (((0.14999) * x33)) +
                                           x36 + x35 + x40);
                                      evalcond[2] = ((0.0172892498998009) +
                                                     (((-2.879808e-7) * x33)) +
                                                     (((0.26) * x39)) +
                                                     (((-1.0) * (pz * pz))) +
                                                     (((-1.92e-6) * x40)) +
                                                     (((-0.0389974) * x32)));
                                      evalcond[3] =
                                          ((((9.6e-7) * x32)) +
                                           (((-0.13) * x33)) +
                                           (((-1.0) * x33 * x37)) + pz +
                                           ((x32 * x35)) + ((x32 * x36)) +
                                           ((x33 * x34)));
                                      evalcond[4] =
                                          ((-0.14999) + (((-1.0) * x32 * x37)) +
                                           (((-9.6e-7) * x33)) +
                                           (((-0.13) * x32)) +
                                           (((-1.0) * x33 * x35)) +
                                           (((-1.0) * x33 * x36)) +
                                           ((x32 * x34)));
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
                              IkReal x41 = cj2 * cj2;
                              IkReal x42 = (cj2 * sj2);
                              CheckValue<IkReal> x43 = IKatan2WithCheck(
                                  IkReal(((1690.00000007225) +
                                          (((0.0221) * sj2)) +
                                          (((-6032.0) * cj2)) +
                                          (((-100000.0) * (pz * pz))) +
                                          (((5382.39999992775) * x41)) +
                                          (((-0.03944) * x42)))),
                                  IkReal(((0.0322) + (((-0.011222) * cj2)) +
                                          (((-5382.39999992775) * x42)) +
                                          (((-0.03944) * x41)) +
                                          (((-14999.0) * pz)) +
                                          (((3016.0000000816) * sj2)))),
                                  IKFAST_ATAN2_MAGTHRESH);
                              if (!x43.valid) {
                                continue;
                              }
                              CheckValue<IkReal> x44 = IKPowWithIntegerCheck(
                                  IKsign(((-1949.87) + (((3479.768) * cj2)) +
                                          (((23200.0) * pz * sj2)) +
                                          (((-0.01274915) * sj2)) +
                                          (((0.096) * pz)) +
                                          (((0.085) * cj2 * pz)))),
                                  -1);
                              if (!x44.valid) {
                                continue;
                              }
                              j1array[0] =
                                  ((-1.5707963267949) + (x43.value) +
                                   (((1.5707963267949) * (x44.value))));
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
                                  IkReal x45 = IKsin(j1);
                                  IkReal x46 = IKcos(j1);
                                  IkReal x47 = ((0.232) * cj2);
                                  IkReal x48 = ((0.232) * sj2);
                                  IkReal x49 = ((8.5e-7) * cj2);
                                  IkReal x50 = ((8.5e-7) * sj2);
                                  IkReal x51 = ((8.5e-7) * x46);
                                  IkReal x52 = (pz * x46);
                                  IkReal x53 = (pz * x45);
                                  evalcond[0] =
                                      ((-0.13) + (((-0.14999) * x45)) + x47 +
                                       x52 + (((-1.0) * x50)));
                                  evalcond[1] = ((9.6e-7) + x48 + x49 + x53 +
                                                 (((0.14999) * x46)));
                                  evalcond[2] = ((0.0172892498998009) +
                                                 (((-2.879808e-7) * x46)) +
                                                 (((-1.92e-6) * x53)) +
                                                 (((-1.0) * (pz * pz))) +
                                                 (((-0.0389974) * x45)) +
                                                 (((0.26) * x52)));
                                  evalcond[3] =
                                      (((x46 * x47)) + (((9.6e-7) * x45)) +
                                       (((-1.0) * x46 * x50)) +
                                       (((-0.13) * x46)) + pz + ((x45 * x49)) +
                                       ((x45 * x48)));
                                  evalcond[4] =
                                      ((-0.14999) + (((-1.0) * x45 * x50)) +
                                       (((-9.6e-7) * x46)) + (((-0.13) * x45)) +
                                       ((x45 * x47)) + (((-1.0) * x46 * x49)) +
                                       (((-1.0) * x46 * x48)));
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
                          IkReal x54 = (pz * sj2);
                          IkReal x55 = (cj2 * pz);
                          CheckValue<IkReal> x56 = IKPowWithIntegerCheck(
                              IKsign(((224970001.0) +
                                      (((10000000000.0) * (pz * pz))))),
                              -1);
                          if (!x56.valid) {
                            continue;
                          }
                          CheckValue<IkReal> x57 = IKatan2WithCheck(
                              IkReal(((-194987000.0) + (((-1274.915) * sj2)) +
                                      (((-2320000000.0) * x54)) +
                                      (((-9600.0) * pz)) + (((-8500.0) * x55)) +
                                      (((347976800.0) * cj2)))),
                              IkReal(((-1439.904) + (((-2320000000.0) * x55)) +
                                      (((-1274.915) * cj2)) +
                                      (((-347976800.0) * sj2)) +
                                      (((8500.0) * x54)) +
                                      (((1300000000.0) * pz)))),
                              IKFAST_ATAN2_MAGTHRESH);
                          if (!x57.valid) {
                            continue;
                          }
                          j1array[0] = ((-1.5707963267949) +
                                        (((1.5707963267949) * (x56.value))) +
                                        (x57.value));
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
                              IkReal x58 = IKsin(j1);
                              IkReal x59 = IKcos(j1);
                              IkReal x60 = ((0.232) * cj2);
                              IkReal x61 = ((0.232) * sj2);
                              IkReal x62 = ((8.5e-7) * cj2);
                              IkReal x63 = ((8.5e-7) * sj2);
                              IkReal x64 = ((8.5e-7) * x59);
                              IkReal x65 = (pz * x59);
                              IkReal x66 = (pz * x58);
                              evalcond[0] =
                                  ((-0.13) + x60 + x65 + (((-0.14999) * x58)) +
                                   (((-1.0) * x63)));
                              evalcond[1] = ((9.6e-7) + (((0.14999) * x59)) +
                                             x61 + x62 + x66);
                              evalcond[2] =
                                  ((0.0172892498998009) +
                                   (((-2.879808e-7) * x59)) +
                                   (((-1.92e-6) * x66)) +
                                   (((-0.0389974) * x58)) +
                                   (((-1.0) * (pz * pz))) + (((0.26) * x65)));
                              evalcond[3] =
                                  (((x59 * x60)) + ((x58 * x62)) +
                                   ((x58 * x61)) + (((-1.0) * x59 * x63)) +
                                   (((9.6e-7) * x58)) + pz + (((-0.13) * x59)));
                              evalcond[4] =
                                  ((-0.14999) + ((x58 * x60)) +
                                   (((-1.0) * x59 * x62)) +
                                   (((-1.0) * x59 * x61)) +
                                   (((-1.0) * x58 * x63)) +
                                   (((-9.6e-7) * x59)) + (((-0.13) * x58)));
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
                    IkReal x67 = IKasin(((0.84697032327434) +
                                         (((-16.5782493363067) * (pz * pz)))));
                    j2array[0] = ((1.57078527838641) + (((-1.0) * x67)));
                    sj2array[0] = IKsin(j2array[0]);
                    cj2array[0] = IKcos(j2array[0]);
                    j2array[1] = ((4.7123779319762) + x67);
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
                        IkReal x68 = pz * pz;
                        IkReal x69 = (pz * sj2);
                        IkReal x70 = (cj2 * pz);
                        j1eval[0] = ((1.0) + (((44.4503709630156) * x68)));
                        j1eval[1] =
                            IKsign(((224970001.0) + (((10000000000.0) * x68))));
                        j1eval[2] =
                            ((IKabs(((-1439.904) + (((-1274.915) * cj2)) +
                                     (((-2320000000.0) * x70)) +
                                     (((-347976800.0) * sj2)) +
                                     (((8500.0) * x69)) +
                                     (((1300000000.0) * pz))))) +
                             (IKabs(((-194987000.0) + (((-1274.915) * sj2)) +
                                     (((-2320000000.0) * x69)) +
                                     (((-8500.0) * x70)) + (((-9600.0) * pz)) +
                                     (((347976800.0) * cj2))))));
                        if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                            IKabs(j1eval[1]) < 0.0000010000000000 ||
                            IKabs(j1eval[2]) < 0.0000010000000000) {
                          {
                            IkReal j1eval[3];
                            px = 0;
                            py = 0;
                            pp = pz * pz;
                            IkReal x71 = cj2 * cj2;
                            IkReal x72 = (cj2 * sj2);
                            IkReal x73 = (pz * sj2);
                            IkReal x74 = (cj2 * pz);
                            j1eval[0] = ((-152941.176470588) +
                                         (((1819729.15841448) * x73)) +
                                         (((-1.0) * sj2)) +
                                         (((272941.176470588) * cj2)) +
                                         (((7.52991375895648) * pz)) +
                                         (((6.66711114074272) * x74)));
                            j1eval[1] =
                                ((IKabs(((1690.00000007225) +
                                         (((0.0221) * sj2)) +
                                         (((5382.39999992775) * x71)) +
                                         (((-6032.0) * cj2)) +
                                         (((-100000.0) * (pz * pz))) +
                                         (((-0.03944) * x72))))) +
                                 (IKabs(((0.0322) + (((-0.011222) * cj2)) +
                                         (((-5382.39999992775) * x72)) +
                                         (((-0.03944) * x71)) +
                                         (((-14999.0) * pz)) +
                                         (((3016.0000000816) * sj2))))));
                            j1eval[2] = IKsign((
                                (-1949.87) + (((23200.0) * x73)) +
                                (((3479.768) * cj2)) + (((-0.01274915) * sj2)) +
                                (((0.096) * pz)) + (((0.085) * x74))));
                            if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                                IKabs(j1eval[1]) < 0.0000010000000000 ||
                                IKabs(j1eval[2]) < 0.0000010000000000) {
                              {
                                IkReal j1eval[3];
                                px = 0;
                                py = 0;
                                pp = pz * pz;
                                IkReal x75 = cj2 * cj2;
                                IkReal x76 = (pz * sj2);
                                IkReal x77 = (cj2 * sj2);
                                IkReal x78 = (cj2 * pz);
                                j1eval[0] = ((-1.12941176470588) +
                                             (((-272941.176470588) * sj2)) +
                                             (((1819729.15841448) * x78)) +
                                             (((-1019675.82152536) * pz)) +
                                             (((-6.66711114074272) * x76)) +
                                             (((-1.0) * cj2)));
                                j1eval[1] = IKsign(
                                    ((-0.01439904) + (((23200.0) * x78)) +
                                     (((-3479.768) * sj2)) +
                                     (((-0.01274915) * cj2)) +
                                     (((-13000.0) * pz)) + (((-0.085) * x76))));
                                j1eval[2] =
                                    ((IKabs(((5382.40000009216) +
                                             (((1.632e-7) * cj2)) +
                                             (((0.044544) * sj2)) +
                                             (((-5382.39999992775) * x75)) +
                                             (((0.03944) * x77)) +
                                             (((-100000.0) * (pz * pz)))))) +
                                     (IKabs(((0.0322) + (((-0.011222) * cj2)) +
                                             (((-5382.39999992775) * x77)) +
                                             (((-0.03944) * x75)) +
                                             (((3016.0000000816) * sj2)) +
                                             (((14999.0) * pz))))));
                                if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                                    IKabs(j1eval[1]) < 0.0000010000000000 ||
                                    IKabs(j1eval[2]) < 0.0000010000000000) {
                                  {
                                    IkReal evalcond[1];
                                    bool bgotonextstatement = true;
                                    do {
                                      IkReal x79 =
                                          ((-272941.176470588) +
                                           (((-6.66711114074272) * pz)));
                                      IkReal x80 =
                                          ((-1.12941176470588) +
                                           (((-1019675.82152536) * pz)));
                                      IkReal x81 =
                                          ((-1.0) +
                                           (((1819729.15841448) * pz)));
                                      IkReal x82 = x80 * x80;
                                      IkReal x83 = ((x81 * x81) + (x79 * x79));
                                      if ((((74496885814.1488) +
                                            (((3311414210028.33) *
                                              (pz * pz))))) < -0.00001)
                                        continue;
                                      IkReal x84 = IKabs(IKsqrt((
                                          (74496885814.1488) +
                                          (((3311414210028.33) * (pz * pz))))));
                                      CheckValue<IkReal> x92 =
                                          IKPowWithIntegerCheck(x84, -1);
                                      if (!x92.valid) {
                                        continue;
                                      }
                                      IkReal x85 = x92.value;
                                      CheckValue<IkReal> x93 =
                                          IKPowWithIntegerCheck(x84, -2);
                                      if (!x93.valid) {
                                        continue;
                                      }
                                      IkReal x86 = x93.value;
                                      IkReal x87 = (x80 * x85);
                                      IkReal x94 = x83;
                                      if (IKabs(x94) == 0) {
                                        continue;
                                      }
                                      IkReal x88 = pow(x94, -0.5);
                                      IkReal x89 = ((1.0) * x88);
                                      IkReal x90 = (x79 * x88);
                                      IkReal x91 = (x82 * x86);
                                      CheckValue<IkReal> x95 = IKatan2WithCheck(
                                          IkReal(x81), IkReal(x79),
                                          IKFAST_ATAN2_MAGTHRESH);
                                      if (!x95.valid) {
                                        continue;
                                      }
                                      if ((x83) < -0.00001)
                                        continue;
                                      CheckValue<IkReal> x96 =
                                          IKPowWithIntegerCheck(
                                              IKabs(IKsqrt(x83)), -1);
                                      if (!x96.valid) {
                                        continue;
                                      }
                                      if ((((1.0) * x80 * (x96.value))) <
                                              -1 - IKFAST_SINCOS_THRESH ||
                                          (((1.0) * x80 * (x96.value))) >
                                              1 + IKFAST_SINCOS_THRESH)
                                        continue;
                                      IkReal gconst0 =
                                          ((((-1.0) * (x95.value))) +
                                           (((-1.0) *
                                             (IKasin(((1.0) * x80 *
                                                      (x96.value)))))));
                                      if ((((1.0) + (((-1.0) * x91)))) <
                                          -0.00001)
                                        continue;
                                      IkReal gconst1 =
                                          ((((-1.0) * x79 * x87 * x89)) +
                                           (((-1.0) * x81 * x89 *
                                             (IKsqrt(((1.0) +
                                                      (((-1.0) * x91))))))));
                                      if ((((1.0) + (((-1.0) * x91)))) <
                                          -0.00001)
                                        continue;
                                      IkReal gconst2 =
                                          (((x90 *
                                             (IKsqrt(((1.0) +
                                                      (((-1.0) * x91))))))) +
                                           (((-1.0) * x81 * x87 * x89)));
                                      IkReal x97 = x79;
                                      IkReal x98 =
                                          ((-1.0) +
                                           (((1819729.15841448) * pz)));
                                      if ((((x97 * x97) + (x98 * x98))) <
                                          -0.00001)
                                        continue;
                                      CheckValue<IkReal> x99 =
                                          IKPowWithIntegerCheck(
                                              IKabs(IKsqrt(
                                                  ((x97 * x97) + (x98 * x98)))),
                                              -1);
                                      if (!x99.valid) {
                                        continue;
                                      }
                                      if ((((-1.0) * (x99.value) *
                                            (((-1.12941176470588) +
                                              (((-1019675.82152536) * pz)))))) <
                                              -1 - IKFAST_SINCOS_THRESH ||
                                          (((-1.0) * (x99.value) *
                                            (((-1.12941176470588) +
                                              (((-1019675.82152536) * pz)))))) >
                                              1 + IKFAST_SINCOS_THRESH)
                                        continue;
                                      CheckValue<IkReal> x100 =
                                          IKatan2WithCheck(
                                              IkReal(x98), IkReal(x97),
                                              IKFAST_ATAN2_MAGTHRESH);
                                      if (!x100.valid) {
                                        continue;
                                      }
                                      if ((((x97 * x97) + (x98 * x98))) <
                                          -0.00001)
                                        continue;
                                      CheckValue<IkReal> x101 =
                                          IKPowWithIntegerCheck(
                                              IKabs(IKsqrt(
                                                  ((x97 * x97) + (x98 * x98)))),
                                              -1);
                                      if (!x101.valid) {
                                        continue;
                                      }
                                      if ((((-1.0) * (x101.value) *
                                            (((-1.12941176470588) +
                                              (((-1019675.82152536) * pz)))))) <
                                              -1 - IKFAST_SINCOS_THRESH ||
                                          (((-1.0) * (x101.value) *
                                            (((-1.12941176470588) +
                                              (((-1019675.82152536) * pz)))))) >
                                              1 + IKFAST_SINCOS_THRESH)
                                        continue;
                                      CheckValue<IkReal> x102 =
                                          IKatan2WithCheck(
                                              IkReal(x98), IkReal(x97),
                                              IKFAST_ATAN2_MAGTHRESH);
                                      if (!x102.valid) {
                                        continue;
                                      }
                                      CheckValue<IkReal> x103 =
                                          IKatan2WithCheck(
                                              IkReal(x98), IkReal(x97),
                                              IKFAST_ATAN2_MAGTHRESH);
                                      if (!x103.valid) {
                                        continue;
                                      }
                                      CheckValue<IkReal> x104 =
                                          IKatan2WithCheck(
                                              IkReal(x98), IkReal(x97),
                                              IKFAST_ATAN2_MAGTHRESH);
                                      if (!x104.valid) {
                                        continue;
                                      }
                                      if ((((x97 * x97) + (x98 * x98))) <
                                          -0.00001)
                                        continue;
                                      CheckValue<IkReal> x105 =
                                          IKPowWithIntegerCheck(
                                              IKabs(IKsqrt(
                                                  ((x97 * x97) + (x98 * x98)))),
                                              -1);
                                      if (!x105.valid) {
                                        continue;
                                      }
                                      if ((((-1.0) * (x105.value) *
                                            (((-1.12941176470588) +
                                              (((-1019675.82152536) * pz)))))) <
                                              -1 - IKFAST_SINCOS_THRESH ||
                                          (((-1.0) * (x105.value) *
                                            (((-1.12941176470588) +
                                              (((-1019675.82152536) * pz)))))) >
                                              1 + IKFAST_SINCOS_THRESH)
                                        continue;
                                      CheckValue<IkReal> x106 =
                                          IKatan2WithCheck(
                                              IkReal(x98), IkReal(x97),
                                              IKFAST_ATAN2_MAGTHRESH);
                                      if (!x106.valid) {
                                        continue;
                                      }
                                      if ((((x97 * x97) + (x98 * x98))) <
                                          -0.00001)
                                        continue;
                                      CheckValue<IkReal> x107 =
                                          IKPowWithIntegerCheck(
                                              IKabs(IKsqrt(
                                                  ((x97 * x97) + (x98 * x98)))),
                                              -1);
                                      if (!x107.valid) {
                                        continue;
                                      }
                                      if ((((-1.0) * (x107.value) *
                                            (((-1.12941176470588) +
                                              (((-1019675.82152536) * pz)))))) <
                                              -1 - IKFAST_SINCOS_THRESH ||
                                          (((-1.0) * (x107.value) *
                                            (((-1.12941176470588) +
                                              (((-1019675.82152536) * pz)))))) >
                                              1 + IKFAST_SINCOS_THRESH)
                                        continue;
                                      CheckValue<IkReal> x108 =
                                          IKatan2WithCheck(
                                              IkReal(x98), IkReal(x97),
                                              IKFAST_ATAN2_MAGTHRESH);
                                      if (!x108.valid) {
                                        continue;
                                      }
                                      if ((((x97 * x97) + (x98 * x98))) <
                                          -0.00001)
                                        continue;
                                      CheckValue<IkReal> x109 =
                                          IKPowWithIntegerCheck(
                                              IKabs(IKsqrt(
                                                  ((x97 * x97) + (x98 * x98)))),
                                              -1);
                                      if (!x109.valid) {
                                        continue;
                                      }
                                      if ((((-1.0) * (x109.value) *
                                            (((-1.12941176470588) +
                                              (((-1019675.82152536) * pz)))))) <
                                              -1 - IKFAST_SINCOS_THRESH ||
                                          (((-1.0) * (x109.value) *
                                            (((-1.12941176470588) +
                                              (((-1019675.82152536) * pz)))))) >
                                              1 + IKFAST_SINCOS_THRESH)
                                        continue;
                                      if ((((x97 * x97) + (x98 * x98))) <
                                          -0.00001)
                                        continue;
                                      CheckValue<IkReal> x110 =
                                          IKPowWithIntegerCheck(
                                              IKabs(IKsqrt(
                                                  ((x97 * x97) + (x98 * x98)))),
                                              -1);
                                      if (!x110.valid) {
                                        continue;
                                      }
                                      if ((((-1.0) * (x110.value) *
                                            (((-1.12941176470588) +
                                              (((-1019675.82152536) * pz)))))) <
                                              -1 - IKFAST_SINCOS_THRESH ||
                                          (((-1.0) * (x110.value) *
                                            (((-1.12941176470588) +
                                              (((-1019675.82152536) * pz)))))) >
                                              1 + IKFAST_SINCOS_THRESH)
                                        continue;
                                      if ((((((-1.0) *
                                              (IKasin(((-1.0) * (x99.value) *
                                                       (((-1.12941176470588) +
                                                         (((-1019675.82152536) *
                                                           pz))))))) *
                                              (x100.value))) +
                                            ((j2 * (j2))) +
                                            (((-1.0) *
                                              (IKasin(((-1.0) * (x101.value) *
                                                       (((-1.12941176470588) +
                                                         (((-1019675.82152536) *
                                                           pz))))))) *
                                              (j2))) +
                                            (((x102.value) * (x103.value))) +
                                            (((-1.0) * (x104.value) *
                                              (IKasin(((-1.0) * (x105.value) *
                                                       (((-1.12941176470588) +
                                                         (((-1019675.82152536) *
                                                           pz))))))))) +
                                            (((x106.value) * (j2))) +
                                            (((-1.0) * j2 *
                                              (IKasin(((-1.0) * (x107.value) *
                                                       (((-1.12941176470588) +
                                                         (((-1019675.82152536) *
                                                           pz))))))))) +
                                            ((j2 * (x108.value))) +
                                            (((1.0) *
                                              (IKasin(((-1.0) * (x109.value) *
                                                       (((-1.12941176470588) +
                                                         (((-1019675.82152536) *
                                                           pz))))))) *
                                              (IKasin(((-1.0) * (x110.value) *
                                                       (((-1.12941176470588) +
                                                         (((-1019675.82152536) *
                                                           pz))))))))))) <
                                          -0.00001)
                                        continue;
                                      evalcond[0] =
                                          ((-3.14159265358979) +
                                           (IKfmod(
                                               ((3.14159265358979) +
                                                (IKsqrt((
                                                    (((-1.0) *
                                                      (IKasin((
                                                          (-1.0) * (x99.value) *
                                                          (((-1.12941176470588) +
                                                            (((-1019675.82152536) *
                                                              pz))))))) *
                                                      (x100.value))) +
                                                    ((j2 * (j2))) +
                                                    (((-1.0) *
                                                      (IKasin((
                                                          (-1.0) *
                                                          (x101.value) *
                                                          (((-1.12941176470588) +
                                                            (((-1019675.82152536) *
                                                              pz))))))) *
                                                      (j2))) +
                                                    (((x102.value) *
                                                      (x103.value))) +
                                                    (((-1.0) * (x104.value) *
                                                      (IKasin((
                                                          (-1.0) *
                                                          (x105.value) *
                                                          (((-1.12941176470588) +
                                                            (((-1019675.82152536) *
                                                              pz))))))))) +
                                                    (((x106.value) * (j2))) +
                                                    (((-1.0) * j2 *
                                                      (IKasin((
                                                          (-1.0) *
                                                          (x107.value) *
                                                          (((-1.12941176470588) +
                                                            (((-1019675.82152536) *
                                                              pz))))))))) +
                                                    ((j2 * (x108.value))) +
                                                    (((1.0) *
                                                      (IKasin((
                                                          (-1.0) *
                                                          (x109.value) *
                                                          (((-1.12941176470588) +
                                                            (((-1019675.82152536) *
                                                              pz))))))) *
                                                      (IKasin((
                                                          (-1.0) *
                                                          (x110.value) *
                                                          (((-1.12941176470588) +
                                                            (((-1019675.82152536) *
                                                              pz))))))))))))),
                                               6.28318530717959)));
                                      if (IKabs(evalcond[0]) <
                                          0.0000050000000000) {
                                        bgotonextstatement = false;
                                        {
                                          IkReal j1eval[2];
                                          IkReal x111 =
                                              ((6.66711114074272) * pz);
                                          IkReal x112 = ((-272941.176470588) +
                                                         (((-1.0) * x111)));
                                          IkReal x113 = pz * pz;
                                          IkReal x114 = x80;
                                          IkReal x115 =
                                              ((-1.0) +
                                               (((1819729.15841448) * pz)));
                                          IkReal x116 = x114 * x114;
                                          IkReal x117 =
                                              ((x112 * x112) + (x115 * x115));
                                          IkReal x118 = x84;
                                          CheckValue<IkReal> x128 =
                                              IKPowWithIntegerCheck(x118, -1);
                                          if (!x128.valid) {
                                            continue;
                                          }
                                          IkReal x119 = x128.value;
                                          CheckValue<IkReal> x129 =
                                              IKPowWithIntegerCheck(x118, -2);
                                          if (!x129.valid) {
                                            continue;
                                          }
                                          IkReal x120 = x129.value;
                                          IkReal x121 = (x114 * x119);
                                          IkReal x130 = x117;
                                          if (IKabs(x130) == 0) {
                                            continue;
                                          }
                                          IkReal x122 = pow(x130, -0.5);
                                          IkReal x123 = ((1.0) * x122);
                                          if ((x117) < -0.00001)
                                            continue;
                                          CheckValue<IkReal> x131 =
                                              IKPowWithIntegerCheck(
                                                  IKabs(IKsqrt(x117)), -1);
                                          if (!x131.valid) {
                                            continue;
                                          }
                                          IkReal x124 = x131.value;
                                          IkReal x125 = (x116 * x120);
                                          IkReal x126 = (x112 * x122);
                                          IkReal x127 = ((1.0) * x124);
                                          px = 0;
                                          py = 0;
                                          pp = x113;
                                          sj2 = gconst1;
                                          cj2 = gconst2;
                                          if ((((-1.0) * x127 *
                                                (((-1.12941176) +
                                                  (((-1020408.16326531) *
                                                    pz)))))) <
                                                  -1 - IKFAST_SINCOS_THRESH ||
                                              (((-1.0) * x127 *
                                                (((-1.12941176) +
                                                  (((-1020408.16326531) *
                                                    pz)))))) >
                                                  1 + IKFAST_SINCOS_THRESH)
                                            continue;
                                          CheckValue<IkReal> x132 =
                                              IKatan2WithCheck(
                                                  IkReal(((-1.0) +
                                                          (((1818181.81818182) *
                                                            pz)))),
                                                  IkReal(((-273224.043715847) +
                                                          (((-1.0) * x111)))),
                                                  IKFAST_ATAN2_MAGTHRESH);
                                          if (!x132.valid) {
                                            continue;
                                          }
                                          j2 = ((IKasin(
                                                    ((-1.0) * x127 *
                                                     (((-1.12941176) +
                                                       (((-1020408.16326531) *
                                                         pz))))))) +
                                                (((-1.0) * (x132.value))));
                                          CheckValue<IkReal> x133 =
                                              IKatan2WithCheck(
                                                  IkReal(x115), IkReal(x112),
                                                  IKFAST_ATAN2_MAGTHRESH);
                                          if (!x133.valid) {
                                            continue;
                                          }
                                          if (((x114 * x127)) <
                                                  -1 - IKFAST_SINCOS_THRESH ||
                                              ((x114 * x127)) >
                                                  1 + IKFAST_SINCOS_THRESH)
                                            continue;
                                          IkReal gconst0 =
                                              ((((-1.0) * (x133.value))) +
                                               (((-1.0) *
                                                 (IKasin((x114 * x127))))));
                                          if ((((1.0) + (((-1.0) * x125)))) <
                                              -0.00001)
                                            continue;
                                          IkReal gconst1 =
                                              ((((-1.0) * x115 * x123 *
                                                 (IKsqrt(
                                                     ((1.0) +
                                                      (((-1.0) * x125))))))) +
                                               (((-1.0) * x112 * x121 * x123)));
                                          if ((((1.0) + (((-1.0) * x125)))) <
                                              -0.00001)
                                            continue;
                                          IkReal gconst2 =
                                              (((x126 * (IKsqrt(((1.0) +
                                                                 (((-1.0) *
                                                                   x125))))))) +
                                               (((-1.0) * x115 * x121 * x123)));
                                          IkReal x134 = pz * pz;
                                          j1eval[0] =
                                              ((1.0) +
                                               (((44.4503709630156) * x134)));
                                          j1eval[1] = IKsign(
                                              ((224970001.0) +
                                               (((10000000000.0) * x134))));
                                          if (IKabs(j1eval[0]) <
                                                  0.0000010000000000 ||
                                              IKabs(j1eval[1]) <
                                                  0.0000010000000000) {
                                            {
                                              IkReal j1array[1], cj1array[1],
                                                  sj1array[1];
                                              bool j1valid[1] = {false};
                                              _nj1 = 1;
                                              IkReal x135 = gconst1 * gconst1;
                                              IkReal x136 = gconst2 * gconst2;
                                              IkReal x137 = (gconst1 * gconst2);
                                              CheckValue<IkReal> x138 =
                                                  IKPowWithIntegerCheck(
                                                      IKsign(
                                                          ((-1949.87) +
                                                           (((3479.768) *
                                                             gconst2)) +
                                                           (((-0.01274915) *
                                                             gconst1)) +
                                                           (((0.096) * pz)) +
                                                           (((23200.0) *
                                                             gconst1 * pz)) +
                                                           (((0.085) * gconst2 *
                                                             pz)))),
                                                      -1);
                                              if (!x138.valid) {
                                                continue;
                                              }
                                              CheckValue<IkReal> x139 =
                                                  IKatan2WithCheck(
                                                      IkReal(
                                                          ((1690.0) +
                                                           (((0.0221) *
                                                             gconst1)) +
                                                           (((-6032.0) *
                                                             gconst2)) +
                                                           (((-0.03944) *
                                                             x137)) +
                                                           (((5382.4) * x136)) +
                                                           (((7.225e-8) *
                                                             x135)) +
                                                           (((-100000.0) *
                                                             (pz * pz))))),
                                                      IkReal((
                                                          (0.01248) +
                                                          (((3016.0000000816) *
                                                            gconst1)) +
                                                          (((-0.011222) *
                                                            gconst2)) +
                                                          (((-0.01972) *
                                                            x136)) +
                                                          (((-5382.39999992775) *
                                                            x137)) +
                                                          (((0.01972) * x135)) +
                                                          (((-14999.0) * pz)))),
                                                      IKFAST_ATAN2_MAGTHRESH);
                                              if (!x139.valid) {
                                                continue;
                                              }
                                              j1array[0] =
                                                  ((-1.5707963267949) +
                                                   (((1.5707963267949) *
                                                     (x138.value))) +
                                                   (x139.value));
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
                                                  IkReal x140 = IKsin(j1);
                                                  IkReal x141 = IKcos(j1);
                                                  IkReal x142 =
                                                      ((0.232) * gconst1);
                                                  IkReal x143 =
                                                      ((8.5e-7) * gconst2);
                                                  IkReal x144 =
                                                      ((8.5e-7) * gconst1);
                                                  IkReal x145 =
                                                      ((0.232) * gconst2);
                                                  IkReal x146 = (pz * x141);
                                                  IkReal x147 = (pz * x140);
                                                  evalcond[0] =
                                                      ((-0.13) +
                                                       (((-1.0) * x144)) +
                                                       x145 + x146 +
                                                       (((-0.14999) * x140)));
                                                  evalcond[1] =
                                                      ((9.6e-7) +
                                                       (((0.14999) * x141)) +
                                                       x142 + x143 + x147);
                                                  evalcond[2] =
                                                      ((0.0172892498998009) +
                                                       (((0.26) * x146)) +
                                                       (((-0.0389974) * x140)) +
                                                       (((-1.92e-6) * x147)) +
                                                       (((-2.879808e-7) *
                                                         x141)) +
                                                       (((-1.0) * (pz * pz))));
                                                  evalcond[3] =
                                                      (((x140 * x143)) +
                                                       ((x140 * x142)) +
                                                       ((x141 * x145)) + pz +
                                                       (((-0.13) * x141)) +
                                                       (((-1.0) * x141 *
                                                         x144)) +
                                                       (((9.6e-7) * x140)));
                                                  evalcond[4] =
                                                      ((-0.14999) +
                                                       ((x140 * x145)) +
                                                       (((-1.0) * x140 *
                                                         x144)) +
                                                       (((-0.13) * x140)) +
                                                       (((-1.0) * x141 *
                                                         x143)) +
                                                       (((-1.0) * x141 *
                                                         x142)) +
                                                       (((-9.6e-7) * x141)));
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
                                              IkReal x148 = (gconst2 * pz);
                                              IkReal x149 = (gconst1 * pz);
                                              CheckValue<IkReal> x150 =
                                                  IKPowWithIntegerCheck(
                                                      IKsign(
                                                          ((224970001.0) +
                                                           (((10000000000.0) *
                                                             (pz * pz))))),
                                                      -1);
                                              if (!x150.valid) {
                                                continue;
                                              }
                                              CheckValue<IkReal> x151 =
                                                  IKatan2WithCheck(
                                                      IkReal((
                                                          (-194987000.0) +
                                                          (((347976800.0) *
                                                            gconst2)) +
                                                          (((-8500.0) * x148)) +
                                                          (((-9600.0) * pz)) +
                                                          (((-2320000000.0) *
                                                            x149)) +
                                                          (((-1274.915) *
                                                            gconst1)))),
                                                      IkReal(
                                                          ((-1439.904) +
                                                           (((8500.0) * x149)) +
                                                           (((-347976800.0) *
                                                             gconst1)) +
                                                           (((-2320000000.0) *
                                                             x148)) +
                                                           (((1300000000.0) *
                                                             pz)) +
                                                           (((-1274.915) *
                                                             gconst2)))),
                                                      IKFAST_ATAN2_MAGTHRESH);
                                              if (!x151.valid) {
                                                continue;
                                              }
                                              j1array[0] =
                                                  ((-1.5707963267949) +
                                                   (((1.5707963267949) *
                                                     (x150.value))) +
                                                   (x151.value));
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
                                                  IkReal x152 = IKsin(j1);
                                                  IkReal x153 = IKcos(j1);
                                                  IkReal x154 =
                                                      ((0.232) * gconst1);
                                                  IkReal x155 =
                                                      ((8.5e-7) * gconst2);
                                                  IkReal x156 =
                                                      ((8.5e-7) * gconst1);
                                                  IkReal x157 =
                                                      ((0.232) * gconst2);
                                                  IkReal x158 = (pz * x153);
                                                  IkReal x159 = (pz * x152);
                                                  evalcond[0] =
                                                      ((-0.13) +
                                                       (((-0.14999) * x152)) +
                                                       (((-1.0) * x156)) +
                                                       x157 + x158);
                                                  evalcond[1] =
                                                      ((9.6e-7) +
                                                       (((0.14999) * x153)) +
                                                       x155 + x154 + x159);
                                                  evalcond[2] =
                                                      ((0.0172892498998009) +
                                                       (((0.26) * x158)) +
                                                       (((-0.0389974) * x152)) +
                                                       (((-2.879808e-7) *
                                                         x153)) +
                                                       (((-1.92e-6) * x159)) +
                                                       (((-1.0) * (pz * pz))));
                                                  evalcond[3] =
                                                      (((x153 * x157)) +
                                                       (((-1.0) * x153 *
                                                         x156)) +
                                                       ((x152 * x154)) +
                                                       ((x152 * x155)) +
                                                       (((-0.13) * x153)) + pz +
                                                       (((9.6e-7) * x152)));
                                                  evalcond[4] =
                                                      ((-0.14999) +
                                                       (((-1.0) * x153 *
                                                         x154)) +
                                                       (((-1.0) * x153 *
                                                         x155)) +
                                                       ((x152 * x157)) +
                                                       (((-1.0) * x152 *
                                                         x156)) +
                                                       (((-0.13) * x152)) +
                                                       (((-9.6e-7) * x153)));
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
                                        IkReal x160 =
                                            ((-272941.176470588) +
                                             (((-6.66711114074272) * pz)));
                                        IkReal x161 =
                                            ((-1.12941176470588) +
                                             (((-1019675.82152536) * pz)));
                                        IkReal x162 =
                                            ((-1.0) +
                                             (((1819729.15841448) * pz)));
                                        IkReal x163 = x161 * x161;
                                        IkReal x164 = ((1.0) * x161);
                                        IkReal x165 =
                                            ((x160 * x160) + (x162 * x162));
                                        if ((((74496885814.1488) +
                                              (((3311414210028.33) *
                                                (pz * pz))))) < -0.00001)
                                          continue;
                                        IkReal x166 =
                                            IKabs(IKsqrt(((74496885814.1488) +
                                                          (((3311414210028.33) *
                                                            (pz * pz))))));
                                        CheckValue<IkReal> x173 =
                                            IKPowWithIntegerCheck(x166, -1);
                                        if (!x173.valid) {
                                          continue;
                                        }
                                        IkReal x167 = x173.value;
                                        CheckValue<IkReal> x174 =
                                            IKPowWithIntegerCheck(x166, -2);
                                        if (!x174.valid) {
                                          continue;
                                        }
                                        IkReal x168 = x174.value;
                                        IkReal x169 = (x161 * x167);
                                        IkReal x175 = x165;
                                        if (IKabs(x175) == 0) {
                                          continue;
                                        }
                                        IkReal x170 = pow(x175, -0.5);
                                        IkReal x171 = ((1.0) * x170);
                                        IkReal x172 = (x163 * x168);
                                        if ((x165) < -0.00001)
                                          continue;
                                        CheckValue<IkReal> x176 =
                                            IKPowWithIntegerCheck(
                                                IKabs(IKsqrt(x165)), -1);
                                        if (!x176.valid) {
                                          continue;
                                        }
                                        if (((x164 * (x176.value))) <
                                                -1 - IKFAST_SINCOS_THRESH ||
                                            ((x164 * (x176.value))) >
                                                1 + IKFAST_SINCOS_THRESH)
                                          continue;
                                        CheckValue<IkReal> x177 =
                                            IKatan2WithCheck(
                                                IkReal(x162), IkReal(x160),
                                                IKFAST_ATAN2_MAGTHRESH);
                                        if (!x177.valid) {
                                          continue;
                                        }
                                        IkReal gconst3 =
                                            ((3.14159265358979) +
                                             (((1.0) *
                                               (IKasin(
                                                   (x164 * (x176.value)))))) +
                                             (((-1.0) * (x177.value))));
                                        if ((((1.0) + (((-1.0) * x172)))) <
                                            -0.00001)
                                          continue;
                                        IkReal gconst4 =
                                            (((x162 * x171 *
                                               (IKsqrt(((1.0) +
                                                        (((-1.0) * x172))))))) +
                                             (((-1.0) * x160 * x164 * x167 *
                                               x170)));
                                        if ((((1.0) + (((-1.0) * x172)))) <
                                            -0.00001)
                                          continue;
                                        IkReal gconst5 =
                                            ((((-1.0) * x162 * x164 * x167 *
                                               x170)) +
                                             (((-1.0) * x160 * x171 *
                                               (IKsqrt(((1.0) +
                                                        (((-1.0) * x172))))))));
                                        IkReal x178 = j2;
                                        CheckValue<IkReal> x184 =
                                            IKatan2WithCheck(
                                                IkReal(((-1.0) +
                                                        (((1819729.15841448) *
                                                          pz)))),
                                                IkReal(((-272941.176470588) +
                                                        (((-6.66711114074272) *
                                                          pz)))),
                                                IKFAST_ATAN2_MAGTHRESH);
                                        if (!x184.valid) {
                                          continue;
                                        }
                                        IkReal x179 = x184.value;
                                        IkReal x180 = x179;
                                        if ((((74496885814.1488) +
                                              (((3311414210028.33) *
                                                (pz * pz))))) < -0.00001)
                                          continue;
                                        CheckValue<IkReal> x185 =
                                            IKPowWithIntegerCheck(
                                                IKabs(IKsqrt(
                                                    ((74496885814.1488) +
                                                     (((3311414210028.33) *
                                                       (pz * pz)))))),
                                                -1);
                                        if (!x185.valid) {
                                          continue;
                                        }
                                        IkReal x181 = x185.value;
                                        if ((((((1019675.82152536) * pz *
                                                x181)) +
                                              (((1.12941176470588) * x181)))) <
                                                -1 - IKFAST_SINCOS_THRESH ||
                                            (((((1019675.82152536) * pz *
                                                x181)) +
                                              (((1.12941176470588) * x181)))) >
                                                1 + IKFAST_SINCOS_THRESH)
                                          continue;
                                        IkReal x182 = IKasin((
                                            (((1019675.82152536) * pz * x181)) +
                                            (((1.12941176470588) * x181))));
                                        IkReal x183 = x182;
                                        if ((((9.86960440108936) +
                                              ((x180 * x182)) +
                                              ((x178 * x182)) +
                                              ((x178 * x179)) +
                                              ((x179 * x183)) +
                                              ((x179 * x180)) +
                                              (((-3.14159265358979) * j2)) +
                                              ((j2 * x178)) + ((j2 * x183)) +
                                              ((j2 * x180)) +
                                              (((-3.14159265358979) * x179)) +
                                              (((-3.14159265358979) * x178)) +
                                              ((x182 * x183)) +
                                              (((-3.14159265358979) * x182)) +
                                              (((-3.14159265358979) * x183)) +
                                              (((-3.14159265358979) * x180)))) <
                                            -0.00001)
                                          continue;
                                        evalcond[0] =
                                            ((-3.14159265358979) +
                                             (IKfmod(
                                                 ((3.14159265358979) +
                                                  (IKsqrt(
                                                      ((9.86960440108936) +
                                                       ((x180 * x182)) +
                                                       ((x178 * x182)) +
                                                       ((x178 * x179)) +
                                                       ((x179 * x183)) +
                                                       ((x179 * x180)) +
                                                       (((-3.14159265358979) *
                                                         j2)) +
                                                       ((j2 * x178)) +
                                                       ((j2 * x183)) +
                                                       ((j2 * x180)) +
                                                       (((-3.14159265358979) *
                                                         x179)) +
                                                       (((-3.14159265358979) *
                                                         x178)) +
                                                       ((x182 * x183)) +
                                                       (((-3.14159265358979) *
                                                         x182)) +
                                                       (((-3.14159265358979) *
                                                         x183)) +
                                                       (((-3.14159265358979) *
                                                         x180)))))),
                                                 6.28318530717959)));
                                        if (IKabs(evalcond[0]) <
                                            0.0000050000000000) {
                                          bgotonextstatement = false;
                                          {
                                            IkReal j1eval[2];
                                            IkReal x186 =
                                                ((6.66711114074272) * pz);
                                            IkReal x187 = ((-272941.176470588) +
                                                           (((-1.0) * x186)));
                                            IkReal x188 = pz * pz;
                                            IkReal x189 = x161;
                                            IkReal x190 =
                                                ((-1.0) +
                                                 (((1819729.15841448) * pz)));
                                            IkReal x191 = x189 * x189;
                                            IkReal x192 = ((1.0) * x189);
                                            IkReal x193 =
                                                ((x187 * x187) + (x190 * x190));
                                            IkReal x194 = x166;
                                            CheckValue<IkReal> x203 =
                                                IKPowWithIntegerCheck(x194, -1);
                                            if (!x203.valid) {
                                              continue;
                                            }
                                            IkReal x195 = x203.value;
                                            CheckValue<IkReal> x204 =
                                                IKPowWithIntegerCheck(x194, -2);
                                            if (!x204.valid) {
                                              continue;
                                            }
                                            IkReal x196 = x204.value;
                                            IkReal x197 = (x189 * x195);
                                            IkReal x205 = x193;
                                            if (IKabs(x205) == 0) {
                                              continue;
                                            }
                                            IkReal x198 = pow(x205, -0.5);
                                            IkReal x199 = ((1.0) * x198);
                                            if ((x193) < -0.00001)
                                              continue;
                                            CheckValue<IkReal> x206 =
                                                IKPowWithIntegerCheck(
                                                    IKabs(IKsqrt(x193)), -1);
                                            if (!x206.valid) {
                                              continue;
                                            }
                                            IkReal x200 = x206.value;
                                            IkReal x201 = ((1.0) * x200);
                                            IkReal x202 = (x191 * x196);
                                            px = 0;
                                            py = 0;
                                            pp = x188;
                                            sj2 = gconst4;
                                            cj2 = gconst5;
                                            CheckValue<IkReal> x207 =
                                                IKatan2WithCheck(
                                                    IkReal(
                                                        ((-1.0) +
                                                         (((1818181.81818182) *
                                                           pz)))),
                                                    IkReal(
                                                        ((-273224.043715847) +
                                                         (((-1.0) * x186)))),
                                                    IKFAST_ATAN2_MAGTHRESH);
                                            if (!x207.valid) {
                                              continue;
                                            }
                                            if ((((-1.0) * x201 *
                                                  (((-1.12941176) +
                                                    (((-1020408.16326531) *
                                                      pz)))))) <
                                                    -1 - IKFAST_SINCOS_THRESH ||
                                                (((-1.0) * x201 *
                                                  (((-1.12941176) +
                                                    (((-1020408.16326531) *
                                                      pz)))))) >
                                                    1 + IKFAST_SINCOS_THRESH)
                                              continue;
                                            j2 =
                                                ((3.14159265) +
                                                 (((-1.0) * (x207.value))) +
                                                 (((-1.0) *
                                                   (IKasin((
                                                       (-1.0) * x201 *
                                                       (((-1.12941176) +
                                                         (((-1020408.16326531) *
                                                           pz))))))))));
                                            if (((x192 * x200)) <
                                                    -1 - IKFAST_SINCOS_THRESH ||
                                                ((x192 * x200)) >
                                                    1 + IKFAST_SINCOS_THRESH)
                                              continue;
                                            CheckValue<IkReal> x208 =
                                                IKatan2WithCheck(
                                                    IkReal(x190), IkReal(x187),
                                                    IKFAST_ATAN2_MAGTHRESH);
                                            if (!x208.valid) {
                                              continue;
                                            }
                                            IkReal gconst3 =
                                                ((3.14159265358979) +
                                                 (((1.0) *
                                                   (IKasin((x192 * x200))))) +
                                                 (((-1.0) * (x208.value))));
                                            if ((((1.0) + (((-1.0) * x202)))) <
                                                -0.00001)
                                              continue;
                                            IkReal gconst4 =
                                                (((x190 * x199 *
                                                   (IKsqrt(
                                                       ((1.0) +
                                                        (((-1.0) * x202))))))) +
                                                 (((-1.0) * x187 * x192 * x195 *
                                                   x198)));
                                            if ((((1.0) + (((-1.0) * x202)))) <
                                                -0.00001)
                                              continue;
                                            IkReal gconst5 =
                                                ((((-1.0) * x190 * x192 * x195 *
                                                   x198)) +
                                                 (((-1.0) * x187 * x199 *
                                                   (IKsqrt(
                                                       ((1.0) +
                                                        (((-1.0) * x202))))))));
                                            IkReal x209 = pz * pz;
                                            j1eval[0] =
                                                ((1.0) +
                                                 (((44.4503709630156) * x209)));
                                            j1eval[1] = IKsign(
                                                ((224970001.0) +
                                                 (((10000000000.0) * x209))));
                                            if (IKabs(j1eval[0]) <
                                                    0.0000010000000000 ||
                                                IKabs(j1eval[1]) <
                                                    0.0000010000000000) {
                                              {
                                                IkReal j1array[1], cj1array[1],
                                                    sj1array[1];
                                                bool j1valid[1] = {false};
                                                _nj1 = 1;
                                                IkReal x210 = gconst4 * gconst4;
                                                IkReal x211 = gconst5 * gconst5;
                                                IkReal x212 =
                                                    (gconst4 * gconst5);
                                                CheckValue<IkReal> x213 =
                                                    IKatan2WithCheck(
                                                        IkReal(((1690.0) +
                                                                (((0.0221) *
                                                                  gconst4)) +
                                                                (((-6032.0) *
                                                                  gconst5)) +
                                                                (((7.225e-8) *
                                                                  x210)) +
                                                                (((-100000.0) *
                                                                  (pz * pz))) +
                                                                (((-0.03944) *
                                                                  x212)) +
                                                                (((5382.4) *
                                                                  x211)))),
                                                        IkReal((
                                                            (0.01248) +
                                                            (((3016.0000000816) *
                                                              gconst4)) +
                                                            (((-0.011222) *
                                                              gconst5)) +
                                                            (((-5382.39999992775) *
                                                              x212)) +
                                                            (((0.01972) *
                                                              x210)) +
                                                            (((-14999.0) *
                                                              pz)) +
                                                            (((-0.01972) *
                                                              x211)))),
                                                        IKFAST_ATAN2_MAGTHRESH);
                                                if (!x213.valid) {
                                                  continue;
                                                }
                                                CheckValue<IkReal> x214 =
                                                    IKPowWithIntegerCheck(
                                                        IKsign(
                                                            ((-1949.87) +
                                                             (((23200.0) *
                                                               gconst4 * pz)) +
                                                             (((3479.768) *
                                                               gconst5)) +
                                                             (((-0.01274915) *
                                                               gconst4)) +
                                                             (((0.085) *
                                                               gconst5 * pz)) +
                                                             (((0.096) * pz)))),
                                                        -1);
                                                if (!x214.valid) {
                                                  continue;
                                                }
                                                j1array[0] =
                                                    ((-1.5707963267949) +
                                                     (x213.value) +
                                                     (((1.5707963267949) *
                                                       (x214.value))));
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
                                                    IkReal x215 = IKsin(j1);
                                                    IkReal x216 = IKcos(j1);
                                                    IkReal x217 =
                                                        ((8.5e-7) * gconst5);
                                                    IkReal x218 =
                                                        ((0.232) * gconst4);
                                                    IkReal x219 =
                                                        ((8.5e-7) * gconst4);
                                                    IkReal x220 =
                                                        ((0.232) * gconst5);
                                                    IkReal x221 = (pz * x216);
                                                    IkReal x222 = (pz * x215);
                                                    evalcond[0] =
                                                        ((-0.13) + x221 + x220 +
                                                         (((-1.0) * x219)) +
                                                         (((-0.14999) * x215)));
                                                    evalcond[1] =
                                                        ((9.6e-7) + x217 +
                                                         x218 + x222 +
                                                         (((0.14999) * x216)));
                                                    evalcond[2] =
                                                        ((0.0172892498998009) +
                                                         (((-0.0389974) *
                                                           x215)) +
                                                         (((-1.0) *
                                                           (pz * pz))) +
                                                         (((0.26) * x221)) +
                                                         (((-1.92e-6) * x222)) +
                                                         (((-2.879808e-7) *
                                                           x216)));
                                                    evalcond[3] =
                                                        (((x216 * x220)) +
                                                         ((x215 * x217)) +
                                                         ((x215 * x218)) + pz +
                                                         (((-0.13) * x216)) +
                                                         (((9.6e-7) * x215)) +
                                                         (((-1.0) * x216 *
                                                           x219)));
                                                    evalcond[4] =
                                                        ((-0.14999) +
                                                         ((x215 * x220)) +
                                                         (((-9.6e-7) * x216)) +
                                                         (((-0.13) * x215)) +
                                                         (((-1.0) * x216 *
                                                           x217)) +
                                                         (((-1.0) * x216 *
                                                           x218)) +
                                                         (((-1.0) * x215 *
                                                           x219)));
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
                                                IkReal x223 = (gconst4 * pz);
                                                IkReal x224 = (gconst5 * pz);
                                                CheckValue<IkReal> x225 =
                                                    IKatan2WithCheck(
                                                        IkReal((
                                                            (-194987000.0) +
                                                            (((-2320000000.0) *
                                                              x223)) +
                                                            (((347976800.0) *
                                                              gconst5)) +
                                                            (((-8500.0) *
                                                              x224)) +
                                                            (((-9600.0) * pz)) +
                                                            (((-1274.915) *
                                                              gconst4)))),
                                                        IkReal(
                                                            ((-1439.904) +
                                                             (((-2320000000.0) *
                                                               x224)) +
                                                             (((8500.0) *
                                                               x223)) +
                                                             (((-347976800.0) *
                                                               gconst4)) +
                                                             (((1300000000.0) *
                                                               pz)) +
                                                             (((-1274.915) *
                                                               gconst5)))),
                                                        IKFAST_ATAN2_MAGTHRESH);
                                                if (!x225.valid) {
                                                  continue;
                                                }
                                                CheckValue<IkReal> x226 =
                                                    IKPowWithIntegerCheck(
                                                        IKsign(
                                                            ((224970001.0) +
                                                             (((10000000000.0) *
                                                               (pz * pz))))),
                                                        -1);
                                                if (!x226.valid) {
                                                  continue;
                                                }
                                                j1array[0] =
                                                    ((-1.5707963267949) +
                                                     (x225.value) +
                                                     (((1.5707963267949) *
                                                       (x226.value))));
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
                                                    IkReal x227 = IKsin(j1);
                                                    IkReal x228 = IKcos(j1);
                                                    IkReal x229 =
                                                        ((8.5e-7) * gconst5);
                                                    IkReal x230 =
                                                        ((0.232) * gconst4);
                                                    IkReal x231 =
                                                        ((8.5e-7) * gconst4);
                                                    IkReal x232 =
                                                        ((0.232) * gconst5);
                                                    IkReal x233 = (pz * x228);
                                                    IkReal x234 = (pz * x227);
                                                    evalcond[0] =
                                                        ((-0.13) +
                                                         (((-0.14999) * x227)) +
                                                         (((-1.0) * x231)) +
                                                         x232 + x233);
                                                    evalcond[1] =
                                                        ((9.6e-7) + x229 +
                                                         x230 + x234 +
                                                         (((0.14999) * x228)));
                                                    evalcond[2] =
                                                        ((0.0172892498998009) +
                                                         (((-2.879808e-7) *
                                                           x228)) +
                                                         (((-1.0) *
                                                           (pz * pz))) +
                                                         (((0.26) * x233)) +
                                                         (((-0.0389974) *
                                                           x227)) +
                                                         (((-1.92e-6) * x234)));
                                                    evalcond[3] =
                                                        ((((9.6e-7) * x227)) +
                                                         ((x227 * x229)) +
                                                         ((x228 * x232)) +
                                                         (((-0.13) * x228)) +
                                                         pz + ((x227 * x230)) +
                                                         (((-1.0) * x228 *
                                                           x231)));
                                                    evalcond[4] =
                                                        ((-0.14999) +
                                                         (((-1.0) * x227 *
                                                           x231)) +
                                                         (((-1.0) * x228 *
                                                           x229)) +
                                                         (((-0.13) * x227)) +
                                                         ((x227 * x232)) +
                                                         (((-9.6e-7) * x228)) +
                                                         (((-1.0) * x228 *
                                                           x230)));
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
                                          IkReal x235 =
                                              ((272941.176470588) +
                                               (((6.66711114074272) * pz)));
                                          IkReal x236 =
                                              ((-1.0) +
                                               (((1819729.15841448) * pz)));
                                          IkReal x237 =
                                              ((-152941.176470588) +
                                               (((7.52991375895648) * pz)));
                                          if ((((74496885814.1488) +
                                                (((3311414210028.33) *
                                                  (pz * pz))))) < -0.00001)
                                            continue;
                                          IkReal x238 = IKabs(
                                              IKsqrt(((74496885814.1488) +
                                                      (((3311414210028.33) *
                                                        (pz * pz))))));
                                          IkReal x239 =
                                              ((x235 * x235) + (x236 * x236));
                                          CheckValue<IkReal> x245 =
                                              IKPowWithIntegerCheck(x238, -1);
                                          if (!x245.valid) {
                                            continue;
                                          }
                                          IkReal x240 = x245.value;
                                          IkReal x246 = x239;
                                          if (IKabs(x246) == 0) {
                                            continue;
                                          }
                                          IkReal x241 = pow(x246, -0.5);
                                          if ((((1.0) +
                                                (((-1.0) * (x237 * x237) *
                                                  (x240 * x240))))) < -0.00001)
                                            continue;
                                          IkReal x242 = IKsqrt((
                                              (1.0) + (((-1.0) * (x237 * x237) *
                                                        (x240 * x240)))));
                                          IkReal x243 = (x237 * x240 * x241);
                                          IkReal x244 = (x241 * x242);
                                          CheckValue<IkReal> x247 =
                                              IKatan2WithCheck(
                                                  IkReal(x235), IkReal(x236),
                                                  IKFAST_ATAN2_MAGTHRESH);
                                          if (!x247.valid) {
                                            continue;
                                          }
                                          if ((x239) < -0.00001)
                                            continue;
                                          CheckValue<IkReal> x248 =
                                              IKPowWithIntegerCheck(
                                                  IKabs(IKsqrt(x239)), -1);
                                          if (!x248.valid) {
                                            continue;
                                          }
                                          if (((x237 * (x248.value))) <
                                                  -1 - IKFAST_SINCOS_THRESH ||
                                              ((x237 * (x248.value))) >
                                                  1 + IKFAST_SINCOS_THRESH)
                                            continue;
                                          IkReal gconst6 =
                                              ((((-1.0) * (x247.value))) +
                                               (((-1.0) *
                                                 (IKasin(
                                                     (x237 * (x248.value)))))));
                                          IkReal gconst7 =
                                              ((((-1.0) * x236 * x243)) +
                                               (((-1.0) * x235 * x244)));
                                          IkReal gconst8 =
                                              (((x236 * x244)) +
                                               (((-1.0) * x235 * x243)));
                                          if ((((74496885814.1488) +
                                                (((3311414210028.33) *
                                                  (pz * pz))))) < -0.00001)
                                            continue;
                                          CheckValue<IkReal> x250 =
                                              IKPowWithIntegerCheck(
                                                  IKabs(IKsqrt(
                                                      ((74496885814.1488) +
                                                       (((3311414210028.33) *
                                                         (pz * pz)))))),
                                                  -1);
                                          if (!x250.valid) {
                                            continue;
                                          }
                                          IkReal x249 = x250.value;
                                          CheckValue<IkReal> x251 =
                                              IKatan2WithCheck(
                                                  IkReal(((272941.176470588) +
                                                          (((6.66711114074272) *
                                                            pz)))),
                                                  IkReal(((-1.0) +
                                                          (((1819729.15841448) *
                                                            pz)))),
                                                  IKFAST_ATAN2_MAGTHRESH);
                                          if (!x251.valid) {
                                            continue;
                                          }
                                          if ((((((7.52991375895648) * pz *
                                                  x249)) +
                                                (((-152941.176470588) *
                                                  x249)))) <
                                                  -1 - IKFAST_SINCOS_THRESH ||
                                              (((((7.52991375895648) * pz *
                                                  x249)) +
                                                (((-152941.176470588) *
                                                  x249)))) >
                                                  1 + IKFAST_SINCOS_THRESH)
                                            continue;
                                          CheckValue<IkReal> x252 =
                                              IKatan2WithCheck(
                                                  IkReal(((272941.176470588) +
                                                          (((6.66711114074272) *
                                                            pz)))),
                                                  IkReal(((-1.0) +
                                                          (((1819729.15841448) *
                                                            pz)))),
                                                  IKFAST_ATAN2_MAGTHRESH);
                                          if (!x252.valid) {
                                            continue;
                                          }
                                          if ((((((7.52991375895648) * pz *
                                                  x249)) +
                                                (((-152941.176470588) *
                                                  x249)))) <
                                                  -1 - IKFAST_SINCOS_THRESH ||
                                              (((((7.52991375895648) * pz *
                                                  x249)) +
                                                (((-152941.176470588) *
                                                  x249)))) >
                                                  1 + IKFAST_SINCOS_THRESH)
                                            continue;
                                          if ((((((7.52991375895648) * pz *
                                                  x249)) +
                                                (((-152941.176470588) *
                                                  x249)))) <
                                                  -1 - IKFAST_SINCOS_THRESH ||
                                              (((((7.52991375895648) * pz *
                                                  x249)) +
                                                (((-152941.176470588) *
                                                  x249)))) >
                                                  1 + IKFAST_SINCOS_THRESH)
                                            continue;
                                          CheckValue<IkReal> x253 =
                                              IKatan2WithCheck(
                                                  IkReal(((272941.176470588) +
                                                          (((6.66711114074272) *
                                                            pz)))),
                                                  IkReal(((-1.0) +
                                                          (((1819729.15841448) *
                                                            pz)))),
                                                  IKFAST_ATAN2_MAGTHRESH);
                                          if (!x253.valid) {
                                            continue;
                                          }
                                          if ((((((7.52991375895648) * pz *
                                                  x249)) +
                                                (((-152941.176470588) *
                                                  x249)))) <
                                                  -1 - IKFAST_SINCOS_THRESH ||
                                              (((((7.52991375895648) * pz *
                                                  x249)) +
                                                (((-152941.176470588) *
                                                  x249)))) >
                                                  1 + IKFAST_SINCOS_THRESH)
                                            continue;
                                          if ((((((7.52991375895648) * pz *
                                                  x249)) +
                                                (((-152941.176470588) *
                                                  x249)))) <
                                                  -1 - IKFAST_SINCOS_THRESH ||
                                              (((((7.52991375895648) * pz *
                                                  x249)) +
                                                (((-152941.176470588) *
                                                  x249)))) >
                                                  1 + IKFAST_SINCOS_THRESH)
                                            continue;
                                          if ((((((7.52991375895648) * pz *
                                                  x249)) +
                                                (((-152941.176470588) *
                                                  x249)))) <
                                                  -1 - IKFAST_SINCOS_THRESH ||
                                              (((((7.52991375895648) * pz *
                                                  x249)) +
                                                (((-152941.176470588) *
                                                  x249)))) >
                                                  1 + IKFAST_SINCOS_THRESH)
                                            continue;
                                          CheckValue<IkReal> x254 =
                                              IKatan2WithCheck(
                                                  IkReal(((272941.176470588) +
                                                          (((6.66711114074272) *
                                                            pz)))),
                                                  IkReal(((-1.0) +
                                                          (((1819729.15841448) *
                                                            pz)))),
                                                  IKFAST_ATAN2_MAGTHRESH);
                                          if (!x254.valid) {
                                            continue;
                                          }
                                          CheckValue<IkReal> x255 =
                                              IKatan2WithCheck(
                                                  IkReal(((272941.176470588) +
                                                          (((6.66711114074272) *
                                                            pz)))),
                                                  IkReal(((-1.0) +
                                                          (((1819729.15841448) *
                                                            pz)))),
                                                  IKFAST_ATAN2_MAGTHRESH);
                                          if (!x255.valid) {
                                            continue;
                                          }
                                          CheckValue<IkReal> x256 =
                                              IKatan2WithCheck(
                                                  IkReal(((272941.176470588) +
                                                          (((6.66711114074272) *
                                                            pz)))),
                                                  IkReal(((-1.0) +
                                                          (((1819729.15841448) *
                                                            pz)))),
                                                  IKFAST_ATAN2_MAGTHRESH);
                                          if (!x256.valid) {
                                            continue;
                                          }
                                          if ((((((x251.value) *
                                                  (IKasin(
                                                      ((((7.52991375895648) *
                                                         pz * x249)) +
                                                       (((-152941.176470588) *
                                                         x249))))))) +
                                                (((x252.value) * (j2))) +
                                                ((j2 * (j2))) +
                                                (((IKasin(
                                                      ((((7.52991375895648) *
                                                         pz * x249)) +
                                                       (((-152941.176470588) *
                                                         x249))))) *
                                                  (IKasin(
                                                      ((((7.52991375895648) *
                                                         pz * x249)) +
                                                       (((-152941.176470588) *
                                                         x249))))))) +
                                                ((j2 * (x253.value))) +
                                                (((IKasin(
                                                      ((((7.52991375895648) *
                                                         pz * x249)) +
                                                       (((-152941.176470588) *
                                                         x249))))) *
                                                  (j2))) +
                                                ((j2 *
                                                  (IKasin(
                                                      ((((7.52991375895648) *
                                                         pz * x249)) +
                                                       (((-152941.176470588) *
                                                         x249))))))) +
                                                (((IKasin(
                                                      ((((7.52991375895648) *
                                                         pz * x249)) +
                                                       (((-152941.176470588) *
                                                         x249))))) *
                                                  (x254.value))) +
                                                (((x255.value) *
                                                  (x256.value))))) < -0.00001)
                                            continue;
                                          evalcond[0] =
                                              ((-3.14159265358979) +
                                               (IKfmod(
                                                   ((3.14159265358979) +
                                                    (IKsqrt((
                                                        (((x251.value) *
                                                          (IKasin((
                                                              (((7.52991375895648) *
                                                                pz * x249)) +
                                                              (((-152941.176470588) *
                                                                x249))))))) +
                                                        (((x252.value) *
                                                          (j2))) +
                                                        ((j2 * (j2))) +
                                                        (((IKasin((
                                                              (((7.52991375895648) *
                                                                pz * x249)) +
                                                              (((-152941.176470588) *
                                                                x249))))) *
                                                          (IKasin((
                                                              (((7.52991375895648) *
                                                                pz * x249)) +
                                                              (((-152941.176470588) *
                                                                x249))))))) +
                                                        ((j2 * (x253.value))) +
                                                        (((IKasin((
                                                              (((7.52991375895648) *
                                                                pz * x249)) +
                                                              (((-152941.176470588) *
                                                                x249))))) *
                                                          (j2))) +
                                                        ((j2 *
                                                          (IKasin((
                                                              (((7.52991375895648) *
                                                                pz * x249)) +
                                                              (((-152941.176470588) *
                                                                x249))))))) +
                                                        (((IKasin((
                                                              (((7.52991375895648) *
                                                                pz * x249)) +
                                                              (((-152941.176470588) *
                                                                x249))))) *
                                                          (x254.value))) +
                                                        (((x255.value) *
                                                          (x256.value))))))),
                                                   6.28318530717959)));
                                          if (IKabs(evalcond[0]) <
                                              0.0000050000000000) {
                                            bgotonextstatement = false;
                                            {
                                              IkReal j1eval[2];
                                              IkReal x257 =
                                                  ((6.66711114074272) * pz);
                                              IkReal x258 = pz * pz;
                                              IkReal x259 =
                                                  ((272941.176470588) + x257);
                                              IkReal x260 =
                                                  ((-1.0) +
                                                   (((1819729.15841448) * pz)));
                                              IkReal x261 =
                                                  ((-152941.176470588) +
                                                   (((7.52991375895648) * pz)));
                                              IkReal x262 = x238;
                                              IkReal x263 = ((x259 * x259) +
                                                             (x260 * x260));
                                              CheckValue<IkReal> x270 =
                                                  IKPowWithIntegerCheck(x262,
                                                                        -1);
                                              if (!x270.valid) {
                                                continue;
                                              }
                                              IkReal x264 = x270.value;
                                              IkReal x271 = x263;
                                              if (IKabs(x271) == 0) {
                                                continue;
                                              }
                                              IkReal x265 = pow(x271, -0.5);
                                              if ((x263) < -0.00001)
                                                continue;
                                              CheckValue<IkReal> x272 =
                                                  IKPowWithIntegerCheck(
                                                      IKabs(IKsqrt(x263)), -1);
                                              if (!x272.valid) {
                                                continue;
                                              }
                                              IkReal x266 = x272.value;
                                              if ((((1.0) +
                                                    (((-1.0) * (x261 * x261) *
                                                      (x264 * x264))))) <
                                                  -0.00001)
                                                continue;
                                              IkReal x267 = IKsqrt(
                                                  ((1.0) +
                                                   (((-1.0) * (x261 * x261) *
                                                     (x264 * x264)))));
                                              IkReal x268 =
                                                  (x261 * x264 * x265);
                                              IkReal x269 = (x265 * x267);
                                              px = 0;
                                              py = 0;
                                              pp = x258;
                                              sj2 = gconst7;
                                              cj2 = gconst8;
                                              CheckValue<IkReal> x273 =
                                                  IKatan2WithCheck(
                                                      IkReal(
                                                          ((273224.043715847) +
                                                           x257)),
                                                      IkReal((
                                                          (-1.0) +
                                                          (((1818181.81818182) *
                                                            pz)))),
                                                      IKFAST_ATAN2_MAGTHRESH);
                                              if (!x273.valid) {
                                                continue;
                                              }
                                              if (((x266 *
                                                    (((-152905.198776758) +
                                                      (((7.52991352270815) *
                                                        pz)))))) <
                                                      -1 -
                                                          IKFAST_SINCOS_THRESH ||
                                                  ((x266 *
                                                    (((-152905.198776758) +
                                                      (((7.52991352270815) *
                                                        pz)))))) >
                                                      1 + IKFAST_SINCOS_THRESH)
                                                continue;
                                              j2 =
                                                  ((((-1.0) * (x273.value))) +
                                                   (((-1.0) *
                                                     (IKasin((
                                                         x266 *
                                                         (((-152905.198776758) +
                                                           (((7.52991352270815) *
                                                             pz))))))))));
                                              if (((x261 * x266)) <
                                                      -1 -
                                                          IKFAST_SINCOS_THRESH ||
                                                  ((x261 * x266)) >
                                                      1 + IKFAST_SINCOS_THRESH)
                                                continue;
                                              CheckValue<IkReal> x274 =
                                                  IKatan2WithCheck(
                                                      IkReal(x259),
                                                      IkReal(x260),
                                                      IKFAST_ATAN2_MAGTHRESH);
                                              if (!x274.valid) {
                                                continue;
                                              }
                                              IkReal gconst6 =
                                                  ((((-1.0) *
                                                     (IKasin((x261 * x266))))) +
                                                   (((-1.0) * (x274.value))));
                                              IkReal gconst7 =
                                                  ((((-1.0) * x259 * x269)) +
                                                   (((-1.0) * x260 * x268)));
                                              IkReal gconst8 =
                                                  (((x260 * x269)) +
                                                   (((-1.0) * x259 * x268)));
                                              IkReal x275 = pz * pz;
                                              j1eval[0] =
                                                  ((1.0) +
                                                   (((44.4503709630156) *
                                                     x275)));
                                              j1eval[1] = IKsign(
                                                  ((224970001.0) +
                                                   (((10000000000.0) * x275))));
                                              if (IKabs(j1eval[0]) <
                                                      0.0000010000000000 ||
                                                  IKabs(j1eval[1]) <
                                                      0.0000010000000000) {
                                                {
                                                  IkReal j1array[1],
                                                      cj1array[1], sj1array[1];
                                                  bool j1valid[1] = {false};
                                                  _nj1 = 1;
                                                  IkReal x276 =
                                                      gconst7 * gconst7;
                                                  IkReal x277 =
                                                      gconst8 * gconst8;
                                                  IkReal x278 =
                                                      (gconst7 * gconst8);
                                                  CheckValue<IkReal> x279 =
                                                      IKPowWithIntegerCheck(
                                                          IKsign((
                                                              (-0.01439904) +
                                                              (((-3479.768) *
                                                                gconst7)) +
                                                              (((-0.01274915) *
                                                                gconst8)) +
                                                              (((-13000.0) *
                                                                pz)) +
                                                              (((23200.0) *
                                                                gconst8 * pz)) +
                                                              (((-0.085) *
                                                                gconst7 *
                                                                pz)))),
                                                          -1);
                                                  if (!x279.valid) {
                                                    continue;
                                                  }
                                                  CheckValue<IkReal> x280 =
                                                      IKatan2WithCheck(
                                                          IkReal((
                                                              (0.01248) +
                                                              (((3016.0000000816) *
                                                                gconst7)) +
                                                              (((0.01972) *
                                                                x276)) +
                                                              (((-0.011222) *
                                                                gconst8)) +
                                                              (((-0.01972) *
                                                                x277)) +
                                                              (((14999.0) *
                                                                pz)) +
                                                              (((-5382.39999992775) *
                                                                x278)))),
                                                          IkReal(
                                                              ((9.216e-8) +
                                                               (((5382.4) *
                                                                 x276)) +
                                                               (((7.225e-8) *
                                                                 x277)) +
                                                               (((1.632e-7) *
                                                                 gconst8)) +
                                                               (((-100000.0) *
                                                                 (pz * pz))) +
                                                               (((0.03944) *
                                                                 x278)) +
                                                               (((0.044544) *
                                                                 gconst7)))),
                                                          IKFAST_ATAN2_MAGTHRESH);
                                                  if (!x280.valid) {
                                                    continue;
                                                  }
                                                  j1array[0] =
                                                      ((-1.5707963267949) +
                                                       (((1.5707963267949) *
                                                         (x279.value))) +
                                                       (x280.value));
                                                  sj1array[0] =
                                                      IKsin(j1array[0]);
                                                  cj1array[0] =
                                                      IKcos(j1array[0]);
                                                  if (j1array[0] > IKPI) {
                                                    j1array[0] -= IK2PI;
                                                  } else if (j1array[0] <
                                                             -IKPI) {
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
                                                          IKabs(
                                                              cj1array[ij1] -
                                                              cj1array[iij1]) <
                                                              IKFAST_SOLUTION_THRESH &&
                                                          IKabs(
                                                              sj1array[ij1] -
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
                                                      IkReal x281 = IKsin(j1);
                                                      IkReal x282 = IKcos(j1);
                                                      IkReal x283 =
                                                          ((8.5e-7) * gconst7);
                                                      IkReal x284 =
                                                          ((8.5e-7) * gconst8);
                                                      IkReal x285 =
                                                          ((0.232) * gconst7);
                                                      IkReal x286 =
                                                          ((0.232) * gconst8);
                                                      IkReal x287 =
                                                          ((0.232) * x281);
                                                      IkReal x288 = (pz * x282);
                                                      IkReal x289 = (pz * x281);
                                                      evalcond[0] =
                                                          ((-0.13) + x288 +
                                                           x286 +
                                                           (((-1.0) * x283)) +
                                                           (((-0.14999) *
                                                             x281)));
                                                      evalcond[1] =
                                                          ((9.6e-7) + x289 +
                                                           x285 + x284 +
                                                           (((0.14999) *
                                                             x282)));
                                                      evalcond[2] =
                                                          ((0.0172892498998009) +
                                                           (((0.26) * x288)) +
                                                           (((-1.92e-6) *
                                                             x289)) +
                                                           (((-0.0389974) *
                                                             x281)) +
                                                           (((-1.0) *
                                                             (pz * pz))) +
                                                           (((-2.879808e-7) *
                                                             x282)));
                                                      evalcond[3] =
                                                          ((((-1.0) * x282 *
                                                             x283)) +
                                                           pz +
                                                           (((9.6e-7) * x281)) +
                                                           ((x281 * x284)) +
                                                           ((x281 * x285)) +
                                                           ((x282 * x286)) +
                                                           (((-0.13) * x282)));
                                                      evalcond[4] =
                                                          ((-0.14999) +
                                                           (((-9.6e-7) *
                                                             x282)) +
                                                           (((-1.0) * x281 *
                                                             x283)) +
                                                           (((-1.0) * x282 *
                                                             x285)) +
                                                           (((-1.0) * x282 *
                                                             x284)) +
                                                           ((x281 * x286)) +
                                                           (((-0.13) * x281)));
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
                                                                  cj0array
                                                                      [ij0] -
                                                                  cj0array
                                                                      [iij0]) <
                                                                  IKFAST_SOLUTION_THRESH &&
                                                              IKabs(
                                                                  sj0array
                                                                      [ij0] -
                                                                  sj0array
                                                                      [iij0]) <
                                                                  IKFAST_SOLUTION_THRESH) {
                                                            j0valid[iij0] =
                                                                false;
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
                                                          vinfos[0].jointtype =
                                                              1;
                                                          vinfos[0].foffset =
                                                              j0;
                                                          vinfos[0].indices[0] =
                                                              _ij0[0];
                                                          vinfos[0].indices[1] =
                                                              _ij0[1];
                                                          vinfos[0]
                                                              .maxsolutions =
                                                              _nj0;
                                                          vinfos[1].jointtype =
                                                              1;
                                                          vinfos[1].foffset =
                                                              j1;
                                                          vinfos[1].indices[0] =
                                                              _ij1[0];
                                                          vinfos[1].indices[1] =
                                                              _ij1[1];
                                                          vinfos[1]
                                                              .maxsolutions =
                                                              _nj1;
                                                          vinfos[2].jointtype =
                                                              1;
                                                          vinfos[2].foffset =
                                                              j2;
                                                          vinfos[2].indices[0] =
                                                              _ij2[0];
                                                          vinfos[2].indices[1] =
                                                              _ij2[1];
                                                          vinfos[2]
                                                              .maxsolutions =
                                                              _nj2;
                                                          std::vector<int>
                                                              vfree(0);
                                                          solutions.AddSolution(
                                                              vinfos, vfree);
                                                        }
                                                      }
                                                    }
                                                  }
                                                }

                                              } else {
                                                {
                                                  IkReal j1array[1],
                                                      cj1array[1], sj1array[1];
                                                  bool j1valid[1] = {false};
                                                  _nj1 = 1;
                                                  IkReal x290 = ((8500.0) * pz);
                                                  IkReal x291 =
                                                      ((2320000000.0) * pz);
                                                  CheckValue<IkReal> x292 =
                                                      IKPowWithIntegerCheck(
                                                          IKsign((
                                                              (224970001.0) +
                                                              (((10000000000.0) *
                                                                (pz * pz))))),
                                                          -1);
                                                  if (!x292.valid) {
                                                    continue;
                                                  }
                                                  CheckValue<IkReal> x293 =
                                                      IKatan2WithCheck(
                                                          IkReal(
                                                              ((-194987000.0) +
                                                               (((347976800.0) *
                                                                 gconst8)) +
                                                               (((-1.0) *
                                                                 gconst8 *
                                                                 x290)) +
                                                               (((-9600.0) *
                                                                 pz)) +
                                                               (((-1.0) *
                                                                 gconst7 *
                                                                 x291)) +
                                                               (((-1274.915) *
                                                                 gconst7)))),
                                                          IkReal((
                                                              (-1439.904) +
                                                              ((gconst7 *
                                                                x290)) +
                                                              (((-1.0) *
                                                                gconst8 *
                                                                x291)) +
                                                              (((-347976800.0) *
                                                                gconst7)) +
                                                              (((1300000000.0) *
                                                                pz)) +
                                                              (((-1274.915) *
                                                                gconst8)))),
                                                          IKFAST_ATAN2_MAGTHRESH);
                                                  if (!x293.valid) {
                                                    continue;
                                                  }
                                                  j1array[0] =
                                                      ((-1.5707963267949) +
                                                       (((1.5707963267949) *
                                                         (x292.value))) +
                                                       (x293.value));
                                                  sj1array[0] =
                                                      IKsin(j1array[0]);
                                                  cj1array[0] =
                                                      IKcos(j1array[0]);
                                                  if (j1array[0] > IKPI) {
                                                    j1array[0] -= IK2PI;
                                                  } else if (j1array[0] <
                                                             -IKPI) {
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
                                                          IKabs(
                                                              cj1array[ij1] -
                                                              cj1array[iij1]) <
                                                              IKFAST_SOLUTION_THRESH &&
                                                          IKabs(
                                                              sj1array[ij1] -
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
                                                      IkReal x294 = IKsin(j1);
                                                      IkReal x295 = IKcos(j1);
                                                      IkReal x296 =
                                                          ((8.5e-7) * gconst7);
                                                      IkReal x297 =
                                                          ((8.5e-7) * gconst8);
                                                      IkReal x298 =
                                                          ((0.232) * gconst7);
                                                      IkReal x299 =
                                                          ((0.232) * gconst8);
                                                      IkReal x300 =
                                                          ((0.232) * x294);
                                                      IkReal x301 = (pz * x295);
                                                      IkReal x302 = (pz * x294);
                                                      evalcond[0] =
                                                          ((-0.13) +
                                                           (((-1.0) * x296)) +
                                                           (((-0.14999) *
                                                             x294)) +
                                                           x301 + x299);
                                                      evalcond[1] =
                                                          ((9.6e-7) +
                                                           (((0.14999) *
                                                             x295)) +
                                                           x302 + x298 + x297);
                                                      evalcond[2] =
                                                          ((0.0172892498998009) +
                                                           (((-1.92e-6) *
                                                             x302)) +
                                                           (((0.26) * x301)) +
                                                           (((-2.879808e-7) *
                                                             x295)) +
                                                           (((-1.0) *
                                                             (pz * pz))) +
                                                           (((-0.0389974) *
                                                             x294)));
                                                      evalcond[3] =
                                                          ((((9.6e-7) * x294)) +
                                                           (((-0.13) * x295)) +
                                                           pz +
                                                           ((x295 * x299)) +
                                                           (((-1.0) * x295 *
                                                             x296)) +
                                                           ((x294 * x298)) +
                                                           ((x294 * x297)));
                                                      evalcond[4] =
                                                          ((-0.14999) +
                                                           (((-0.13) * x294)) +
                                                           (((-1.0) * x295 *
                                                             x298)) +
                                                           (((-1.0) * x295 *
                                                             x297)) +
                                                           ((x294 * x299)) +
                                                           (((-1.0) * x294 *
                                                             x296)) +
                                                           (((-9.6e-7) *
                                                             x295)));
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
                                                                  cj0array
                                                                      [ij0] -
                                                                  cj0array
                                                                      [iij0]) <
                                                                  IKFAST_SOLUTION_THRESH &&
                                                              IKabs(
                                                                  sj0array
                                                                      [ij0] -
                                                                  sj0array
                                                                      [iij0]) <
                                                                  IKFAST_SOLUTION_THRESH) {
                                                            j0valid[iij0] =
                                                                false;
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
                                                          vinfos[0].jointtype =
                                                              1;
                                                          vinfos[0].foffset =
                                                              j0;
                                                          vinfos[0].indices[0] =
                                                              _ij0[0];
                                                          vinfos[0].indices[1] =
                                                              _ij0[1];
                                                          vinfos[0]
                                                              .maxsolutions =
                                                              _nj0;
                                                          vinfos[1].jointtype =
                                                              1;
                                                          vinfos[1].foffset =
                                                              j1;
                                                          vinfos[1].indices[0] =
                                                              _ij1[0];
                                                          vinfos[1].indices[1] =
                                                              _ij1[1];
                                                          vinfos[1]
                                                              .maxsolutions =
                                                              _nj1;
                                                          vinfos[2].jointtype =
                                                              1;
                                                          vinfos[2].foffset =
                                                              j2;
                                                          vinfos[2].indices[0] =
                                                              _ij2[0];
                                                          vinfos[2].indices[1] =
                                                              _ij2[1];
                                                          vinfos[2]
                                                              .maxsolutions =
                                                              _nj2;
                                                          std::vector<int>
                                                              vfree(0);
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
                                            IkReal x303 =
                                                ((272941.176470588) +
                                                 (((6.66711114074272) * pz)));
                                            IkReal x304 =
                                                ((-1.0) +
                                                 (((1819729.15841448) * pz)));
                                            IkReal x305 =
                                                ((-152941.176470588) +
                                                 (((7.52991375895648) * pz)));
                                            if ((((74496885814.1488) +
                                                  (((3311414210028.33) *
                                                    (pz * pz))))) < -0.00001)
                                              continue;
                                            IkReal x306 = IKabs(
                                                IKsqrt(((74496885814.1488) +
                                                        (((3311414210028.33) *
                                                          (pz * pz))))));
                                            IkReal x307 =
                                                ((x303 * x303) + (x304 * x304));
                                            CheckValue<IkReal> x314 =
                                                IKPowWithIntegerCheck(x306, -1);
                                            if (!x314.valid) {
                                              continue;
                                            }
                                            IkReal x308 = x314.value;
                                            IkReal x309 = (x305 * x308);
                                            IkReal x315 = x307;
                                            if (IKabs(x315) == 0) {
                                              continue;
                                            }
                                            IkReal x310 = pow(x315, -0.5);
                                            IkReal x311 = ((1.0) * x303 * x310);
                                            IkReal x312 = ((1.0) * x304 * x310);
                                            if ((((1.0) +
                                                  (((-1.0) * (x309 * x309))))) <
                                                -0.00001)
                                              continue;
                                            IkReal x313 = IKsqrt(
                                                ((1.0) +
                                                 (((-1.0) * (x309 * x309)))));
                                            if ((x307) < -0.00001)
                                              continue;
                                            CheckValue<IkReal> x316 =
                                                IKPowWithIntegerCheck(
                                                    IKabs(IKsqrt(x307)), -1);
                                            if (!x316.valid) {
                                              continue;
                                            }
                                            if (((x305 * (x316.value))) <
                                                    -1 - IKFAST_SINCOS_THRESH ||
                                                ((x305 * (x316.value))) >
                                                    1 + IKFAST_SINCOS_THRESH)
                                              continue;
                                            CheckValue<IkReal> x317 =
                                                IKatan2WithCheck(
                                                    IkReal(x303), IkReal(x304),
                                                    IKFAST_ATAN2_MAGTHRESH);
                                            if (!x317.valid) {
                                              continue;
                                            }
                                            IkReal gconst9 =
                                                ((3.14159265358979) +
                                                 (IKasin(
                                                     (x305 * (x316.value)))) +
                                                 (((-1.0) * (x317.value))));
                                            IkReal gconst10 =
                                                (((x311 * x313)) +
                                                 (((-1.0) * x309 * x312)));
                                            IkReal gconst11 =
                                                ((((-1.0) * x312 * x313)) +
                                                 (((-1.0) * x309 * x311)));
                                            IkReal x318 =
                                                ((272941.176470588) +
                                                 (((6.66711114074272) * pz)));
                                            IkReal x319 =
                                                ((-1.0) +
                                                 (((1819729.15841448) * pz)));
                                            if ((((x319 * x319) +
                                                  (x318 * x318))) < -0.00001)
                                              continue;
                                            CheckValue<IkReal> x320 =
                                                IKPowWithIntegerCheck(
                                                    IKabs(IKsqrt(
                                                        ((x319 * x319) +
                                                         (x318 * x318)))),
                                                    -1);
                                            if (!x320.valid) {
                                              continue;
                                            }
                                            if ((((x320.value) *
                                                  (((-152941.176470588) +
                                                    (((7.52991375895648) *
                                                      pz)))))) <
                                                    -1 - IKFAST_SINCOS_THRESH ||
                                                (((x320.value) *
                                                  (((-152941.176470588) +
                                                    (((7.52991375895648) *
                                                      pz)))))) >
                                                    1 + IKFAST_SINCOS_THRESH)
                                              continue;
                                            CheckValue<IkReal> x321 =
                                                IKatan2WithCheck(
                                                    IkReal(x318), IkReal(x319),
                                                    IKFAST_ATAN2_MAGTHRESH);
                                            if (!x321.valid) {
                                              continue;
                                            }
                                            evalcond[0] =
                                                ((-3.14159265358979) +
                                                 (IKfmod(
                                                     ((3.14159265358979) +
                                                      (IKabs((
                                                          (-3.14159265358979) +
                                                          (((-1.0) *
                                                            (IKasin((
                                                                (x320.value) *
                                                                (((-152941.176470588) +
                                                                  (((7.52991375895648) *
                                                                    pz))))))))) +
                                                          (x321.value) + j2)))),
                                                     6.28318530717959)));
                                            if (IKabs(evalcond[0]) <
                                                0.0000050000000000) {
                                              bgotonextstatement = false;
                                              {
                                                IkReal j1eval[2];
                                                IkReal x322 =
                                                    ((6.66711114074272) * pz);
                                                IkReal x323 = pz * pz;
                                                IkReal x324 =
                                                    ((272941.176470588) + x322);
                                                IkReal x325 =
                                                    ((-1.0) +
                                                     (((1819729.15841448) *
                                                       pz)));
                                                IkReal x326 =
                                                    ((-152941.176470588) +
                                                     (((7.52991375895648) *
                                                       pz)));
                                                IkReal x327 = x306;
                                                IkReal x328 = ((x325 * x325) +
                                                               (x324 * x324));
                                                CheckValue<IkReal> x336 =
                                                    IKPowWithIntegerCheck(x327,
                                                                          -1);
                                                if (!x336.valid) {
                                                  continue;
                                                }
                                                IkReal x329 = x336.value;
                                                IkReal x330 = (x326 * x329);
                                                IkReal x337 = x328;
                                                if (IKabs(x337) == 0) {
                                                  continue;
                                                }
                                                IkReal x331 = pow(x337, -0.5);
                                                if ((x328) < -0.00001)
                                                  continue;
                                                CheckValue<IkReal> x338 =
                                                    IKPowWithIntegerCheck(
                                                        IKabs(IKsqrt(x328)),
                                                        -1);
                                                if (!x338.valid) {
                                                  continue;
                                                }
                                                IkReal x332 = x338.value;
                                                IkReal x333 =
                                                    ((1.0) * x324 * x331);
                                                IkReal x334 =
                                                    ((1.0) * x325 * x331);
                                                if ((((1.0) +
                                                      (((-1.0) *
                                                        (x330 * x330))))) <
                                                    -0.00001)
                                                  continue;
                                                IkReal x335 = IKsqrt((
                                                    (1.0) + (((-1.0) *
                                                              (x330 * x330)))));
                                                px = 0;
                                                py = 0;
                                                pp = x323;
                                                sj2 = gconst10;
                                                cj2 = gconst11;
                                                if (((x332 *
                                                      (((-152905.198776758) +
                                                        (((7.52991352270815) *
                                                          pz)))))) <
                                                        -1 -
                                                            IKFAST_SINCOS_THRESH ||
                                                    ((x332 *
                                                      (((-152905.198776758) +
                                                        (((7.52991352270815) *
                                                          pz)))))) >
                                                        1 + IKFAST_SINCOS_THRESH)
                                                  continue;
                                                CheckValue<IkReal> x339 =
                                                    IKatan2WithCheck(
                                                        IkReal((
                                                            (273224.043715847) +
                                                            x322)),
                                                        IkReal((
                                                            (-1.0) +
                                                            (((1818181.81818182) *
                                                              pz)))),
                                                        IKFAST_ATAN2_MAGTHRESH);
                                                if (!x339.valid) {
                                                  continue;
                                                }
                                                j2 =
                                                    ((3.14159265) +
                                                     (IKasin((
                                                         x332 *
                                                         (((-152905.198776758) +
                                                           (((7.52991352270815) *
                                                             pz))))))) +
                                                     (((-1.0) * (x339.value))));
                                                CheckValue<IkReal> x340 =
                                                    IKatan2WithCheck(
                                                        IkReal(x324),
                                                        IkReal(x325),
                                                        IKFAST_ATAN2_MAGTHRESH);
                                                if (!x340.valid) {
                                                  continue;
                                                }
                                                if (((x326 * x332)) <
                                                        -1 -
                                                            IKFAST_SINCOS_THRESH ||
                                                    ((x326 * x332)) >
                                                        1 + IKFAST_SINCOS_THRESH)
                                                  continue;
                                                IkReal gconst9 =
                                                    ((3.14159265358979) +
                                                     (((-1.0) * (x340.value))) +
                                                     (IKasin((x326 * x332))));
                                                IkReal gconst10 =
                                                    ((((-1.0) * x330 * x334)) +
                                                     ((x333 * x335)));
                                                IkReal gconst11 =
                                                    ((((-1.0) * x330 * x333)) +
                                                     (((-1.0) * x334 * x335)));
                                                IkReal x341 = pz * pz;
                                                j1eval[0] =
                                                    ((1.0) +
                                                     (((44.4503709630156) *
                                                       x341)));
                                                j1eval[1] =
                                                    IKsign(((224970001.0) +
                                                            (((10000000000.0) *
                                                              x341))));
                                                if (IKabs(j1eval[0]) <
                                                        0.0000010000000000 ||
                                                    IKabs(j1eval[1]) <
                                                        0.0000010000000000) {
                                                  {
                                                    IkReal j1array[1],
                                                        cj1array[1],
                                                        sj1array[1];
                                                    bool j1valid[1] = {false};
                                                    _nj1 = 1;
                                                    IkReal x342 =
                                                        gconst10 * gconst10;
                                                    IkReal x343 =
                                                        gconst11 * gconst11;
                                                    IkReal x344 =
                                                        (gconst10 * gconst11);
                                                    CheckValue<IkReal> x345 =
                                                        IKPowWithIntegerCheck(
                                                            IKsign((
                                                                (-1949.87) +
                                                                (((3479.768) *
                                                                  gconst11)) +
                                                                (((23200.0) *
                                                                  gconst10 *
                                                                  pz)) +
                                                                (((0.096) *
                                                                  pz)) +
                                                                (((-0.01274915) *
                                                                  gconst10)) +
                                                                (((0.085) *
                                                                  gconst11 *
                                                                  pz)))),
                                                            -1);
                                                    if (!x345.valid) {
                                                      continue;
                                                    }
                                                    CheckValue<IkReal> x346 = IKatan2WithCheck(
                                                        IkReal(((1690.0) +
                                                                (((7.225e-8) *
                                                                  x342)) +
                                                                (((0.0221) *
                                                                  gconst10)) +
                                                                (((-100000.0) *
                                                                  (pz * pz))) +
                                                                (((-6032.0) *
                                                                  gconst11)) +
                                                                (((-0.03944) *
                                                                  x344)) +
                                                                (((5382.4) *
                                                                  x343)))),
                                                        IkReal((
                                                            (0.01248) +
                                                            (((-0.01972) *
                                                              x343)) +
                                                            (((3016.0000000816) *
                                                              gconst10)) +
                                                            (((-0.011222) *
                                                              gconst11)) +
                                                            (((-5382.39999992775) *
                                                              x344)) +
                                                            (((-14999.0) *
                                                              pz)) +
                                                            (((0.01972) *
                                                              x342)))),
                                                        IKFAST_ATAN2_MAGTHRESH);
                                                    if (!x346.valid) {
                                                      continue;
                                                    }
                                                    j1array[0] =
                                                        ((-1.5707963267949) +
                                                         (((1.5707963267949) *
                                                           (x345.value))) +
                                                         (x346.value));
                                                    sj1array[0] =
                                                        IKsin(j1array[0]);
                                                    cj1array[0] =
                                                        IKcos(j1array[0]);
                                                    if (j1array[0] > IKPI) {
                                                      j1array[0] -= IK2PI;
                                                    } else if (j1array[0] <
                                                               -IKPI) {
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
                                                            IKabs(
                                                                cj1array[ij1] -
                                                                cj1array
                                                                    [iij1]) <
                                                                IKFAST_SOLUTION_THRESH &&
                                                            IKabs(
                                                                sj1array[ij1] -
                                                                sj1array
                                                                    [iij1]) <
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
                                                        IkReal x347 = IKsin(j1);
                                                        IkReal x348 = IKcos(j1);
                                                        IkReal x349 =
                                                            ((8.5e-7) *
                                                             gconst11);
                                                        IkReal x350 =
                                                            ((0.232) *
                                                             gconst10);
                                                        IkReal x351 =
                                                            ((0.232) *
                                                             gconst11);
                                                        IkReal x352 =
                                                            ((8.5e-7) *
                                                             gconst10);
                                                        IkReal x353 =
                                                            ((0.232) * x347);
                                                        IkReal x354 =
                                                            ((8.5e-7) * x347);
                                                        IkReal x355 =
                                                            (pz * x348);
                                                        IkReal x356 =
                                                            (pz * x347);
                                                        evalcond[0] =
                                                            ((-0.13) +
                                                             (((-0.14999) *
                                                               x347)) +
                                                             (((-1.0) * x352)) +
                                                             x351 + x355);
                                                        evalcond[1] =
                                                            ((9.6e-7) +
                                                             (((0.14999) *
                                                               x348)) +
                                                             x350 + x356 +
                                                             x349);
                                                        evalcond[2] =
                                                            ((0.0172892498998009) +
                                                             (((-0.0389974) *
                                                               x347)) +
                                                             (((-1.0) *
                                                               (pz * pz))) +
                                                             (((-2.879808e-7) *
                                                               x348)) +
                                                             (((0.26) * x355)) +
                                                             (((-1.92e-6) *
                                                               x356)));
                                                        evalcond[3] =
                                                            ((((-1.0) * x348 *
                                                               x352)) +
                                                             (((9.6e-7) *
                                                               x347)) +
                                                             ((x347 * x350)) +
                                                             pz +
                                                             ((x347 * x349)) +
                                                             ((x348 * x351)) +
                                                             (((-0.13) *
                                                               x348)));
                                                        evalcond[4] =
                                                            ((-0.14999) +
                                                             (((-1.0) * x348 *
                                                               x350)) +
                                                             (((-1.0) * x348 *
                                                               x349)) +
                                                             ((x347 * x351)) +
                                                             (((-9.6e-7) *
                                                               x348)) +
                                                             (((-1.0) * x347 *
                                                               x352)) +
                                                             (((-0.13) *
                                                               x347)));
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
                                                        bool j0valid[1] = {
                                                            false};
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
                                                        for (int ij0 = 0;
                                                             ij0 < 1; ++ij0) {
                                                          if (!j0valid[ij0]) {
                                                            continue;
                                                          }
                                                          _ij0[0] = ij0;
                                                          _ij0[1] = -1;
                                                          for (int iij0 =
                                                                   ij0 + 1;
                                                               iij0 < 1;
                                                               ++iij0) {
                                                            if (j0valid[iij0] &&
                                                                IKabs(
                                                                    cj0array
                                                                        [ij0] -
                                                                    cj0array
                                                                        [iij0]) <
                                                                    IKFAST_SOLUTION_THRESH &&
                                                                IKabs(
                                                                    sj0array
                                                                        [ij0] -
                                                                    sj0array
                                                                        [iij0]) <
                                                                    IKFAST_SOLUTION_THRESH) {
                                                              j0valid[iij0] =
                                                                  false;
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
                                                            vinfos[0]
                                                                .jointtype = 1;
                                                            vinfos[0].foffset =
                                                                j0;
                                                            vinfos[0]
                                                                .indices[0] =
                                                                _ij0[0];
                                                            vinfos[0]
                                                                .indices[1] =
                                                                _ij0[1];
                                                            vinfos[0]
                                                                .maxsolutions =
                                                                _nj0;
                                                            vinfos[1]
                                                                .jointtype = 1;
                                                            vinfos[1].foffset =
                                                                j1;
                                                            vinfos[1]
                                                                .indices[0] =
                                                                _ij1[0];
                                                            vinfos[1]
                                                                .indices[1] =
                                                                _ij1[1];
                                                            vinfos[1]
                                                                .maxsolutions =
                                                                _nj1;
                                                            vinfos[2]
                                                                .jointtype = 1;
                                                            vinfos[2].foffset =
                                                                j2;
                                                            vinfos[2]
                                                                .indices[0] =
                                                                _ij2[0];
                                                            vinfos[2]
                                                                .indices[1] =
                                                                _ij2[1];
                                                            vinfos[2]
                                                                .maxsolutions =
                                                                _nj2;
                                                            std::vector<int>
                                                                vfree(0);
                                                            solutions
                                                                .AddSolution(
                                                                    vinfos,
                                                                    vfree);
                                                          }
                                                        }
                                                      }
                                                    }
                                                  }

                                                } else {
                                                  {
                                                    IkReal j1array[1],
                                                        cj1array[1],
                                                        sj1array[1];
                                                    bool j1valid[1] = {false};
                                                    _nj1 = 1;
                                                    IkReal x357 =
                                                        (gconst10 * pz);
                                                    IkReal x358 =
                                                        (gconst11 * pz);
                                                    CheckValue<IkReal> x359 =
                                                        IKatan2WithCheck(
                                                            IkReal((
                                                                (-194987000.0) +
                                                                (((-1274.915) *
                                                                  gconst10)) +
                                                                (((-8500.0) *
                                                                  x358)) +
                                                                (((-9600.0) *
                                                                  pz)) +
                                                                (((-2320000000.0) *
                                                                  x357)) +
                                                                (((347976800.0) *
                                                                  gconst11)))),
                                                            IkReal((
                                                                (-1439.904) +
                                                                (((-1274.915) *
                                                                  gconst11)) +
                                                                (((8500.0) *
                                                                  x357)) +
                                                                (((-2320000000.0) *
                                                                  x358)) +
                                                                (((1300000000.0) *
                                                                  pz)) +
                                                                (((-347976800.0) *
                                                                  gconst10)))),
                                                            IKFAST_ATAN2_MAGTHRESH);
                                                    if (!x359.valid) {
                                                      continue;
                                                    }
                                                    CheckValue<IkReal> x360 =
                                                        IKPowWithIntegerCheck(
                                                            IKsign((
                                                                (224970001.0) +
                                                                (((10000000000.0) *
                                                                  (pz * pz))))),
                                                            -1);
                                                    if (!x360.valid) {
                                                      continue;
                                                    }
                                                    j1array[0] =
                                                        ((-1.5707963267949) +
                                                         (x359.value) +
                                                         (((1.5707963267949) *
                                                           (x360.value))));
                                                    sj1array[0] =
                                                        IKsin(j1array[0]);
                                                    cj1array[0] =
                                                        IKcos(j1array[0]);
                                                    if (j1array[0] > IKPI) {
                                                      j1array[0] -= IK2PI;
                                                    } else if (j1array[0] <
                                                               -IKPI) {
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
                                                            IKabs(
                                                                cj1array[ij1] -
                                                                cj1array
                                                                    [iij1]) <
                                                                IKFAST_SOLUTION_THRESH &&
                                                            IKabs(
                                                                sj1array[ij1] -
                                                                sj1array
                                                                    [iij1]) <
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
                                                        IkReal x361 = IKsin(j1);
                                                        IkReal x362 = IKcos(j1);
                                                        IkReal x363 =
                                                            ((8.5e-7) *
                                                             gconst11);
                                                        IkReal x364 =
                                                            ((0.232) *
                                                             gconst10);
                                                        IkReal x365 =
                                                            ((0.232) *
                                                             gconst11);
                                                        IkReal x366 =
                                                            ((8.5e-7) *
                                                             gconst10);
                                                        IkReal x367 =
                                                            ((0.232) * x361);
                                                        IkReal x368 =
                                                            ((8.5e-7) * x361);
                                                        IkReal x369 =
                                                            (pz * x362);
                                                        IkReal x370 =
                                                            (pz * x361);
                                                        evalcond[0] =
                                                            ((-0.13) +
                                                             (((-1.0) * x366)) +
                                                             x369 + x365 +
                                                             (((-0.14999) *
                                                               x361)));
                                                        evalcond[1] =
                                                            ((9.6e-7) + x370 +
                                                             x363 + x364 +
                                                             (((0.14999) *
                                                               x362)));
                                                        evalcond[2] =
                                                            ((0.0172892498998009) +
                                                             (((-0.0389974) *
                                                               x361)) +
                                                             (((-1.92e-6) *
                                                               x370)) +
                                                             (((-1.0) *
                                                               (pz * pz))) +
                                                             (((0.26) * x369)) +
                                                             (((-2.879808e-7) *
                                                               x362)));
                                                        evalcond[3] =
                                                            (((x361 * x363)) +
                                                             ((x361 * x364)) +
                                                             (((-1.0) * x362 *
                                                               x366)) +
                                                             (((-0.13) *
                                                               x362)) +
                                                             pz +
                                                             ((x362 * x365)) +
                                                             (((9.6e-7) *
                                                               x361)));
                                                        evalcond[4] =
                                                            ((-0.14999) +
                                                             ((x361 * x365)) +
                                                             (((-9.6e-7) *
                                                               x362)) +
                                                             (((-1.0) * x362 *
                                                               x363)) +
                                                             (((-1.0) * x362 *
                                                               x364)) +
                                                             (((-1.0) * x361 *
                                                               x366)) +
                                                             (((-0.13) *
                                                               x361)));
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
                                                        bool j0valid[1] = {
                                                            false};
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
                                                        for (int ij0 = 0;
                                                             ij0 < 1; ++ij0) {
                                                          if (!j0valid[ij0]) {
                                                            continue;
                                                          }
                                                          _ij0[0] = ij0;
                                                          _ij0[1] = -1;
                                                          for (int iij0 =
                                                                   ij0 + 1;
                                                               iij0 < 1;
                                                               ++iij0) {
                                                            if (j0valid[iij0] &&
                                                                IKabs(
                                                                    cj0array
                                                                        [ij0] -
                                                                    cj0array
                                                                        [iij0]) <
                                                                    IKFAST_SOLUTION_THRESH &&
                                                                IKabs(
                                                                    sj0array
                                                                        [ij0] -
                                                                    sj0array
                                                                        [iij0]) <
                                                                    IKFAST_SOLUTION_THRESH) {
                                                              j0valid[iij0] =
                                                                  false;
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
                                                            vinfos[0]
                                                                .jointtype = 1;
                                                            vinfos[0].foffset =
                                                                j0;
                                                            vinfos[0]
                                                                .indices[0] =
                                                                _ij0[0];
                                                            vinfos[0]
                                                                .indices[1] =
                                                                _ij0[1];
                                                            vinfos[0]
                                                                .maxsolutions =
                                                                _nj0;
                                                            vinfos[1]
                                                                .jointtype = 1;
                                                            vinfos[1].foffset =
                                                                j1;
                                                            vinfos[1]
                                                                .indices[0] =
                                                                _ij1[0];
                                                            vinfos[1]
                                                                .indices[1] =
                                                                _ij1[1];
                                                            vinfos[1]
                                                                .maxsolutions =
                                                                _nj1;
                                                            vinfos[2]
                                                                .jointtype = 1;
                                                            vinfos[2].foffset =
                                                                j2;
                                                            vinfos[2]
                                                                .indices[0] =
                                                                _ij2[0];
                                                            vinfos[2]
                                                                .indices[1] =
                                                                _ij2[1];
                                                            vinfos[2]
                                                                .maxsolutions =
                                                                _nj2;
                                                            std::vector<int>
                                                                vfree(0);
                                                            solutions
                                                                .AddSolution(
                                                                    vinfos,
                                                                    vfree);
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
                                                continue; // branch miss [j0,
                                                          // j1]
                                              }
                                            } while (0);
                                            if (bgotonextstatement) {
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
                                    IkReal x371 = cj2 * cj2;
                                    IkReal x372 = (cj2 * sj2);
                                    CheckValue<IkReal> x373 = IKatan2WithCheck(
                                        IkReal(((0.0322) +
                                                (((-0.03944) * x371)) +
                                                (((-0.011222) * cj2)) +
                                                (((-5382.39999992775) * x372)) +
                                                (((3016.0000000816) * sj2)) +
                                                (((14999.0) * pz)))),
                                        IkReal(
                                            ((5382.40000009216) +
                                             (((1.632e-7) * cj2)) +
                                             (((0.044544) * sj2)) +
                                             (((-100000.0) * (pz * pz))) +
                                             (((0.03944) * x372)) +
                                             (((-5382.39999992775) * x371)))),
                                        IKFAST_ATAN2_MAGTHRESH);
                                    if (!x373.valid) {
                                      continue;
                                    }
                                    CheckValue<IkReal> x374 =
                                        IKPowWithIntegerCheck(
                                            IKsign(((-0.01439904) +
                                                    (((23200.0) * cj2 * pz)) +
                                                    (((-3479.768) * sj2)) +
                                                    (((-0.085) * pz * sj2)) +
                                                    (((-0.01274915) * cj2)) +
                                                    (((-13000.0) * pz)))),
                                            -1);
                                    if (!x374.valid) {
                                      continue;
                                    }
                                    j1array[0] =
                                        ((-1.5707963267949) + (x373.value) +
                                         (((1.5707963267949) * (x374.value))));
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
                                        IkReal x375 = IKsin(j1);
                                        IkReal x376 = IKcos(j1);
                                        IkReal x377 = ((0.232) * cj2);
                                        IkReal x378 = ((0.232) * sj2);
                                        IkReal x379 = ((8.5e-7) * cj2);
                                        IkReal x380 = ((8.5e-7) * sj2);
                                        IkReal x381 = ((8.5e-7) * x376);
                                        IkReal x382 = (pz * x376);
                                        IkReal x383 = (pz * x375);
                                        evalcond[0] =
                                            ((-0.13) + (((-0.14999) * x375)) +
                                             x382 + x377 + (((-1.0) * x380)));
                                        evalcond[1] =
                                            ((9.6e-7) + x383 + x379 + x378 +
                                             (((0.14999) * x376)));
                                        evalcond[2] =
                                            ((0.0172892498998009) +
                                             (((-1.92e-6) * x383)) +
                                             (((0.26) * x382)) +
                                             (((-0.0389974) * x375)) +
                                             (((-1.0) * (pz * pz))) +
                                             (((-2.879808e-7) * x376)));
                                        evalcond[3] =
                                            ((((9.6e-7) * x375)) +
                                             ((x375 * x379)) + ((x375 * x378)) +
                                             (((-1.0) * x376 * x380)) +
                                             (((-0.13) * x376)) + pz +
                                             ((x376 * x377)));
                                        evalcond[4] =
                                            ((-0.14999) + ((x375 * x377)) +
                                             (((-1.0) * x375 * x380)) +
                                             (((-0.13) * x375)) +
                                             (((-1.0) * x376 * x378)) +
                                             (((-1.0) * x376 * x379)) +
                                             (((-9.6e-7) * x376)));
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
                                IkReal x384 = cj2 * cj2;
                                IkReal x385 = (cj2 * sj2);
                                CheckValue<IkReal> x386 = IKatan2WithCheck(
                                    IkReal(((1690.00000007225) +
                                            (((0.0221) * sj2)) +
                                            (((5382.39999992775) * x384)) +
                                            (((-0.03944) * x385)) +
                                            (((-6032.0) * cj2)) +
                                            (((-100000.0) * (pz * pz))))),
                                    IkReal(((0.0322) + (((-0.03944) * x384)) +
                                            (((-5382.39999992775) * x385)) +
                                            (((-0.011222) * cj2)) +
                                            (((-14999.0) * pz)) +
                                            (((3016.0000000816) * sj2)))),
                                    IKFAST_ATAN2_MAGTHRESH);
                                if (!x386.valid) {
                                  continue;
                                }
                                CheckValue<IkReal> x387 = IKPowWithIntegerCheck(
                                    IKsign(((-1949.87) + (((3479.768) * cj2)) +
                                            (((23200.0) * pz * sj2)) +
                                            (((-0.01274915) * sj2)) +
                                            (((0.096) * pz)) +
                                            (((0.085) * cj2 * pz)))),
                                    -1);
                                if (!x387.valid) {
                                  continue;
                                }
                                j1array[0] =
                                    ((-1.5707963267949) + (x386.value) +
                                     (((1.5707963267949) * (x387.value))));
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
                                    IkReal x388 = IKsin(j1);
                                    IkReal x389 = IKcos(j1);
                                    IkReal x390 = ((0.232) * cj2);
                                    IkReal x391 = ((0.232) * sj2);
                                    IkReal x392 = ((8.5e-7) * cj2);
                                    IkReal x393 = ((8.5e-7) * sj2);
                                    IkReal x394 = ((8.5e-7) * x389);
                                    IkReal x395 = (pz * x389);
                                    IkReal x396 = (pz * x388);
                                    evalcond[0] =
                                        ((-0.13) + (((-1.0) * x393)) +
                                         (((-0.14999) * x388)) + x395 + x390);
                                    evalcond[1] =
                                        ((9.6e-7) + (((0.14999) * x389)) +
                                         x396 + x391 + x392);
                                    evalcond[2] = ((0.0172892498998009) +
                                                   (((-0.0389974) * x388)) +
                                                   (((-2.879808e-7) * x389)) +
                                                   (((-1.92e-6) * x396)) +
                                                   (((-1.0) * (pz * pz))) +
                                                   (((0.26) * x395)));
                                    evalcond[3] =
                                        ((((-0.13) * x389)) +
                                         (((-1.0) * x389 * x393)) +
                                         (((9.6e-7) * x388)) + pz +
                                         ((x389 * x390)) + ((x388 * x392)) +
                                         ((x388 * x391)));
                                    evalcond[4] =
                                        ((-0.14999) + (((-0.13) * x388)) +
                                         (((-1.0) * x389 * x391)) +
                                         (((-1.0) * x389 * x392)) +
                                         (((-9.6e-7) * x389)) +
                                         (((-1.0) * x388 * x393)) +
                                         ((x388 * x390)));
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
                            IkReal x397 = (pz * sj2);
                            IkReal x398 = (cj2 * pz);
                            CheckValue<IkReal> x399 = IKPowWithIntegerCheck(
                                IKsign(((224970001.0) +
                                        (((10000000000.0) * (pz * pz))))),
                                -1);
                            if (!x399.valid) {
                              continue;
                            }
                            CheckValue<IkReal> x400 = IKatan2WithCheck(
                                IkReal(((-194987000.0) + (((-1274.915) * sj2)) +
                                        (((-2320000000.0) * x397)) +
                                        (((-8500.0) * x398)) +
                                        (((-9600.0) * pz)) +
                                        (((347976800.0) * cj2)))),
                                IkReal(((-1439.904) +
                                        (((-2320000000.0) * x398)) +
                                        (((-1274.915) * cj2)) +
                                        (((8500.0) * x397)) +
                                        (((-347976800.0) * sj2)) +
                                        (((1300000000.0) * pz)))),
                                IKFAST_ATAN2_MAGTHRESH);
                            if (!x400.valid) {
                              continue;
                            }
                            j1array[0] = ((-1.5707963267949) +
                                          (((1.5707963267949) * (x399.value))) +
                                          (x400.value));
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
                                IkReal x401 = IKsin(j1);
                                IkReal x402 = IKcos(j1);
                                IkReal x403 = ((0.232) * cj2);
                                IkReal x404 = ((0.232) * sj2);
                                IkReal x405 = ((8.5e-7) * cj2);
                                IkReal x406 = ((8.5e-7) * sj2);
                                IkReal x407 = ((8.5e-7) * x402);
                                IkReal x408 = (pz * x402);
                                IkReal x409 = (pz * x401);
                                evalcond[0] =
                                    ((-0.13) + x403 + x408 +
                                     (((-0.14999) * x401)) + (((-1.0) * x406)));
                                evalcond[1] = ((9.6e-7) + x405 + x404 + x409 +
                                               (((0.14999) * x402)));
                                evalcond[2] = ((0.0172892498998009) +
                                               (((-0.0389974) * x401)) +
                                               (((0.26) * x408)) +
                                               (((-2.879808e-7) * x402)) +
                                               (((-1.0) * (pz * pz))) +
                                               (((-1.92e-6) * x409)));
                                evalcond[3] =
                                    ((((-1.0) * x402 * x406)) + pz +
                                     (((9.6e-7) * x401)) + (((-0.13) * x402)) +
                                     ((x402 * x403)) + ((x401 * x405)) +
                                     ((x401 * x404)));
                                evalcond[4] =
                                    ((-0.14999) + (((-9.6e-7) * x402)) +
                                     (((-1.0) * x402 * x404)) +
                                     (((-1.0) * x402 * x405)) +
                                     (((-1.0) * x401 * x406)) +
                                     (((-0.13) * x401)) + ((x401 * x403)));
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
            CheckValue<IkReal> x412 = IKatan2WithCheck(IkReal(px), IkReal(py),
                                                       IKFAST_ATAN2_MAGTHRESH);
            if (!x412.valid) {
              continue;
            }
            IkReal x410 = ((1.0) * (x412.value));
            if ((((px * px) + (py * py))) < -0.00001)
              continue;
            CheckValue<IkReal> x413 = IKPowWithIntegerCheck(
                IKabs(IKsqrt(((px * px) + (py * py)))), -1);
            if (!x413.valid) {
              continue;
            }
            if ((((0.0535) * (x413.value))) < -1 - IKFAST_SINCOS_THRESH ||
                (((0.0535) * (x413.value))) > 1 + IKFAST_SINCOS_THRESH)
              continue;
            IkReal x411 = IKasin(((0.0535) * (x413.value)));
            j0array[0] = ((((-1.0) * x410)) + (((-1.0) * x411)));
            sj0array[0] = IKsin(j0array[0]);
            cj0array[0] = IKcos(j0array[0]);
            j0array[1] = ((3.14159265358979) + x411 + (((-1.0) * x410)));
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
                IkReal x414 = IKasin(((-0.84697032327434) +
                                      (((16.5782493363067) * (px * px))) +
                                      (((4.97314323590529) * cj0 * py)) +
                                      (((-4.97314323590529) * px * sj0)) +
                                      (((16.5782493363067) * (pz * pz))) +
                                      (((16.5782493363067) * (py * py)))));
                j2array[0] = ((1.57078527838641) + (((1.0) * x414)));
                sj2array[0] = IKsin(j2array[0]);
                cj2array[0] = IKcos(j2array[0]);
                j2array[1] = ((4.7123779319762) + (((-1.0) * x414)));
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
                    IkReal x415 = cj2 * cj2;
                    IkReal x416 = (px * sj0);
                    IkReal x417 = ((6.66711114074272) * sj2);
                    IkReal x418 = (pz * sj2);
                    IkReal x419 = (cj0 * py);
                    IkReal x420 = ((0.232) * cj2);
                    IkReal x421 = ((8.5e-7) * sj2);
                    IkReal x422 = ((1819729.15841448) * cj2);
                    IkReal x423 = (cj2 * pz);
                    IkReal x424 = (cj2 * sj2);
                    j1eval[0] =
                        ((152941.176470588) + sj2 +
                         (((-1019675.82152536) * x416)) +
                         (((-7.52991375895648) * pz)) + ((x416 * x422)) +
                         ((x417 * x419)) + (((-6.66711114074272) * x423)) +
                         (((-1.0) * x416 * x417)) + (((-1.0) * x419 * x422)) +
                         (((-1819729.15841448) * x418)) +
                         (((1019675.82152536) * x419)) +
                         (((-272941.176470588) * cj2)));
                    j1eval[1] =
                        ((IKabs(((-0.0169000000007225) +
                                 (((-0.0538239999992775) * x415)) + (pz * pz) +
                                 (((-2.21e-7) * sj2)) + (((3.944e-7) * x424)) +
                                 (((0.06032) * cj2))))) +
                         (IKabs(((-3.22e-7) + ((pz * x419)) +
                                 (((-1.0) * pz * x416)) + (((0.14999) * pz)) +
                                 (((-0.030160000000816) * sj2)) +
                                 (((3.944e-7) * x415)) + (((1.1222e-7) * cj2)) +
                                 (((0.0538239999992775) * x424))))));
                    j1eval[2] = IKsign(
                        ((0.0194987) + (((-1.0) * x416 * x421)) +
                         (((-8.5e-7) * x423)) + (((-9.6e-7) * pz)) +
                         (((1.274915e-7) * sj2)) + (((-0.232) * x418)) +
                         ((x419 * x421)) + ((x416 * x420)) +
                         (((-1.0) * x419 * x420)) + (((-0.03479768) * cj2)) +
                         (((-0.13) * x416)) + (((0.13) * x419))));
                    if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                        IKabs(j1eval[1]) < 0.0000010000000000 ||
                        IKabs(j1eval[2]) < 0.0000010000000000) {
                      {
                        IkReal j1eval[3];
                        IkReal x425 = cj2 * cj2;
                        IkReal x426 = (cj0 * py);
                        IkReal x427 = ((0.232) * sj2);
                        IkReal x428 = (px * sj0);
                        IkReal x429 = (cj2 * pz);
                        IkReal x430 = ((6.66711114074272) * cj2);
                        IkReal x431 = (pz * sj2);
                        IkReal x432 = ((8.5e-7) * cj2);
                        IkReal x433 = ((1819729.15841448) * sj2);
                        IkReal x434 = (cj2 * sj2);
                        j1eval[0] =
                            ((1.12941176470588) + cj2 + ((x426 * x433)) +
                             ((x426 * x430)) + (((7.52991375895648) * x426)) +
                             (((-7.52991375895648) * x428)) +
                             (((-1819729.15841448) * x429)) +
                             (((6.66711114074272) * x431)) +
                             (((-1.0) * x428 * x433)) +
                             (((-1.0) * x428 * x430)) +
                             (((1019675.82152536) * pz)) +
                             (((272941.176470588) * sj2)));
                        j1eval[1] =
                            ((IKabs(((-3.22e-7) + (((-1.0) * pz * x426)) +
                                     (((-0.14999) * pz)) +
                                     (((-0.030160000000816) * sj2)) +
                                     ((pz * x428)) + (((1.1222e-7) * cj2)) +
                                     (((3.944e-7) * x425)) +
                                     (((0.0538239999992775) * x434))))) +
                             (IKabs(((-0.0538240000009216) + (pz * pz) +
                                     (((-4.4544e-7) * sj2)) +
                                     (((-1.632e-12) * cj2)) +
                                     (((0.0538239999992775) * x425)) +
                                     (((-3.944e-7) * x434))))));
                        j1eval[2] = IKsign(
                            ((1.439904e-7) + ((x426 * x427)) + ((x426 * x432)) +
                             (((-1.0) * x427 * x428)) +
                             (((1.274915e-7) * cj2)) + (((-0.232) * x429)) +
                             (((-9.6e-7) * x428)) + (((8.5e-7) * x431)) +
                             (((-1.0) * x428 * x432)) + (((0.13) * pz)) +
                             (((0.03479768) * sj2)) + (((9.6e-7) * x426))));
                        if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                            IKabs(j1eval[1]) < 0.0000010000000000 ||
                            IKabs(j1eval[2]) < 0.0000010000000000) {
                          {
                            IkReal j1eval[2];
                            j1eval[0] = ((106122.08151018) + sj2 +
                                         (((-90510.7736605966) * cj2)));
                            j1eval[1] = IKsign(((0.0707240000016441) +
                                                (((-0.060319999998368) * cj2)) +
                                                (((6.6644e-7) * sj2))));
                            if (IKabs(j1eval[0]) < 0.0000010000000000 ||
                                IKabs(j1eval[1]) < 0.0000010000000000) {
                              continue; // no branches [j1]

                            } else {
                              {
                                IkReal j1array[1], cj1array[1], sj1array[1];
                                bool j1valid[1] = {false};
                                _nj1 = 1;
                                IkReal x435 = (cj0 * py);
                                IkReal x436 = ((0.232) * sj2);
                                IkReal x437 = ((0.232) * cj2);
                                IkReal x438 = (px * sj0);
                                IkReal x439 = ((8.5e-7) * cj2);
                                IkReal x440 = ((8.5e-7) * sj2);
                                CheckValue<IkReal> x441 = IKPowWithIntegerCheck(
                                    IKsign(((0.0707240000016441) +
                                            (((-0.060319999998368) * cj2)) +
                                            (((6.6644e-7) * sj2)))),
                                    -1);
                                if (!x441.valid) {
                                  continue;
                                }
                                CheckValue<IkReal> x442 = IKatan2WithCheck(
                                    IkReal(((-0.0194987) +
                                            (((-1.0) * x437 * x438)) +
                                            (((-1.274915e-7) * sj2)) +
                                            (((-1.0) * x435 * x440)) +
                                            (((0.13) * x438)) +
                                            (((-9.6e-7) * pz)) +
                                            (((-0.13) * x435)) +
                                            (((-1.0) * pz * x439)) +
                                            (((-1.0) * pz * x436)) +
                                            ((x438 * x440)) +
                                            (((0.03479768) * cj2)) +
                                            ((x435 * x437)))),
                                    IkReal(
                                        ((-1.439904e-7) + ((x436 * x438)) +
                                         (((-1.0) * x435 * x439)) +
                                         (((-1.0) * x435 * x436)) +
                                         (((-1.274915e-7) * cj2)) +
                                         (((-0.03479768) * sj2)) +
                                         ((pz * x440)) + (((-9.6e-7) * x435)) +
                                         (((-1.0) * pz * x437)) +
                                         (((0.13) * pz)) + (((9.6e-7) * x438)) +
                                         ((x438 * x439)))),
                                    IKFAST_ATAN2_MAGTHRESH);
                                if (!x442.valid) {
                                  continue;
                                }
                                j1array[0] =
                                    ((-1.5707963267949) +
                                     (((1.5707963267949) * (x441.value))) +
                                     (x442.value));
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
                                    IkReal x443 = IKsin(j1);
                                    IkReal x444 = IKcos(j1);
                                    IkReal x445 = ((0.232) * cj2);
                                    IkReal x446 = (px * sj0);
                                    IkReal x447 = ((0.232) * sj2);
                                    IkReal x448 = (cj0 * py);
                                    IkReal x449 = ((8.5e-7) * cj2);
                                    IkReal x450 = ((8.5e-7) * sj2);
                                    IkReal x451 = ((8.5e-7) * x444);
                                    IkReal x452 = (pz * x444);
                                    IkReal x453 = (pz * x443);
                                    IkReal x454 = ((0.26) * x443);
                                    IkReal x455 = ((1.92e-6) * x444);
                                    evalcond[0] =
                                        ((-0.13) + (((-1.0) * x450)) + x445 +
                                         x452 + (((-0.14999) * x443)) +
                                         (((-1.0) * x443 * x448)) +
                                         ((x443 * x446)));
                                    evalcond[1] =
                                        ((9.6e-7) + (((-1.0) * x444 * x446)) +
                                         x449 + x447 + x453 + ((x444 * x448)) +
                                         (((0.14999) * x444)));
                                    evalcond[2] =
                                        ((((-0.13) * x444)) + ((x444 * x445)) +
                                         pz + (((9.6e-7) * x443)) +
                                         ((x443 * x449)) + ((x443 * x447)) +
                                         (((-1.0) * x444 * x450)));
                                    evalcond[3] =
                                        ((-0.14999) + (((-1.0) * x443 * x450)) +
                                         (((-1.0) * x444 * x449)) +
                                         (((-1.0) * x444 * x447)) +
                                         (((-1.0) * x448)) +
                                         (((-9.6e-7) * x444)) +
                                         (((-0.13) * x443)) + x446 +
                                         ((x443 * x445)));
                                    evalcond[4] =
                                        ((0.0172892498998009) +
                                         (((0.29998) * x446)) +
                                         (((-1.0) * (px * px))) +
                                         (((0.26) * x452)) +
                                         (((-1.92e-6) * x453)) +
                                         ((x446 * x455)) + ((x446 * x454)) +
                                         (((-0.0389974) * x443)) +
                                         (((-1.0) * (pz * pz))) +
                                         (((-0.29998) * x448)) +
                                         (((-2.879808e-7) * x444)) +
                                         (((-1.0) * (py * py))) +
                                         (((-1.0) * x448 * x454)) +
                                         (((-1.0) * x448 * x455)));
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
                            IkReal x456 = cj2 * cj2;
                            IkReal x457 = (cj0 * py);
                            IkReal x458 = ((0.232) * sj2);
                            IkReal x459 = (px * sj0);
                            IkReal x460 = ((8.5e-7) * cj2);
                            IkReal x461 = (cj2 * sj2);
                            CheckValue<IkReal> x462 = IKPowWithIntegerCheck(
                                IKsign(
                                    ((1.439904e-7) + (((9.6e-7) * x457)) +
                                     (((1.274915e-7) * cj2)) +
                                     (((-1.0) * x459 * x460)) +
                                     (((-0.232) * cj2 * pz)) + (((0.13) * pz)) +
                                     (((-1.0) * x458 * x459)) +
                                     (((-9.6e-7) * x459)) +
                                     (((8.5e-7) * pz * sj2)) + ((x457 * x460)) +
                                     (((0.03479768) * sj2)) + ((x457 * x458)))),
                                -1);
                            if (!x462.valid) {
                              continue;
                            }
                            CheckValue<IkReal> x463 = IKatan2WithCheck(
                                IkReal(((-3.22e-7) + (((-1.0) * pz * x457)) +
                                        (((-0.14999) * pz)) +
                                        (((-0.030160000000816) * sj2)) +
                                        (((0.0538239999992775) * x461)) +
                                        (((1.1222e-7) * cj2)) + ((pz * x459)) +
                                        (((3.944e-7) * x456)))),
                                IkReal(((-0.0538240000009216) +
                                        (((-3.944e-7) * x461)) + (pz * pz) +
                                        (((0.0538239999992775) * x456)) +
                                        (((-4.4544e-7) * sj2)) +
                                        (((-1.632e-12) * cj2)))),
                                IKFAST_ATAN2_MAGTHRESH);
                            if (!x463.valid) {
                              continue;
                            }
                            j1array[0] = ((-1.5707963267949) +
                                          (((1.5707963267949) * (x462.value))) +
                                          (x463.value));
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
                                IkReal x464 = IKsin(j1);
                                IkReal x465 = IKcos(j1);
                                IkReal x466 = ((0.232) * cj2);
                                IkReal x467 = (px * sj0);
                                IkReal x468 = ((0.232) * sj2);
                                IkReal x469 = (cj0 * py);
                                IkReal x470 = ((8.5e-7) * cj2);
                                IkReal x471 = ((8.5e-7) * sj2);
                                IkReal x472 = ((8.5e-7) * x465);
                                IkReal x473 = (pz * x465);
                                IkReal x474 = (pz * x464);
                                IkReal x475 = ((0.26) * x464);
                                IkReal x476 = ((1.92e-6) * x465);
                                evalcond[0] = ((-0.13) + (((-1.0) * x471)) +
                                               (((-0.14999) * x464)) +
                                               (((-1.0) * x464 * x469)) + x466 +
                                               x473 + ((x464 * x467)));
                                evalcond[1] = ((9.6e-7) + (((0.14999) * x465)) +
                                               ((x465 * x469)) + x468 + x474 +
                                               x470 + (((-1.0) * x465 * x467)));
                                evalcond[2] =
                                    ((((-0.13) * x465)) + (((9.6e-7) * x464)) +
                                     ((x465 * x466)) + pz + ((x464 * x470)) +
                                     ((x464 * x468)) +
                                     (((-1.0) * x465 * x471)));
                                evalcond[3] =
                                    ((-0.14999) + (((-1.0) * x465 * x468)) +
                                     (((-0.13) * x464)) + (((-1.0) * x469)) +
                                     (((-1.0) * x464 * x471)) + x467 +
                                     ((x464 * x466)) +
                                     (((-1.0) * x465 * x470)) +
                                     (((-9.6e-7) * x465)));
                                evalcond[4] =
                                    ((0.0172892498998009) +
                                     (((-1.0) * (px * px))) +
                                     (((0.26) * x473)) +
                                     (((-1.0) * x469 * x475)) +
                                     (((-1.0) * x469 * x476)) +
                                     (((0.29998) * x467)) + ((x467 * x475)) +
                                     ((x467 * x476)) + (((-1.92e-6) * x474)) +
                                     (((-0.29998) * x469)) +
                                     (((-1.0) * (pz * pz))) +
                                     (((-0.0389974) * x464)) +
                                     (((-1.0) * (py * py))) +
                                     (((-2.879808e-7) * x465)));
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
                        IkReal x477 = cj2 * cj2;
                        IkReal x478 = (cj0 * py);
                        IkReal x479 = (px * sj0);
                        IkReal x480 = ((0.232) * cj2);
                        IkReal x481 = ((8.5e-7) * sj2);
                        IkReal x482 = (cj2 * sj2);
                        CheckValue<IkReal> x483 = IKatan2WithCheck(
                            IkReal(((-0.0169000000007225) +
                                    (((3.944e-7) * x482)) + (pz * pz) +
                                    (((-2.21e-7) * sj2)) +
                                    (((-0.0538239999992775) * x477)) +
                                    (((0.06032) * cj2)))),
                            IkReal(((-3.22e-7) + (((3.944e-7) * x477)) +
                                    (((0.0538239999992775) * x482)) +
                                    (((0.14999) * pz)) +
                                    (((-0.030160000000816) * sj2)) +
                                    (((-1.0) * pz * x479)) +
                                    (((1.1222e-7) * cj2)) + ((pz * x478)))),
                            IKFAST_ATAN2_MAGTHRESH);
                        if (!x483.valid) {
                          continue;
                        }
                        CheckValue<IkReal> x484 = IKPowWithIntegerCheck(
                            IKsign(((0.0194987) + (((-0.13) * x479)) +
                                    (((-9.6e-7) * pz)) +
                                    (((1.274915e-7) * sj2)) +
                                    (((0.13) * x478)) + ((x478 * x481)) +
                                    (((-8.5e-7) * cj2 * pz)) +
                                    (((-1.0) * x479 * x481)) +
                                    (((-0.232) * pz * sj2)) + ((x479 * x480)) +
                                    (((-0.03479768) * cj2)) +
                                    (((-1.0) * x478 * x480)))),
                            -1);
                        if (!x484.valid) {
                          continue;
                        }
                        j1array[0] = ((-1.5707963267949) + (x483.value) +
                                      (((1.5707963267949) * (x484.value))));
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
                            IkReal x485 = IKsin(j1);
                            IkReal x486 = IKcos(j1);
                            IkReal x487 = ((0.232) * cj2);
                            IkReal x488 = (px * sj0);
                            IkReal x489 = ((0.232) * sj2);
                            IkReal x490 = (cj0 * py);
                            IkReal x491 = ((8.5e-7) * cj2);
                            IkReal x492 = ((8.5e-7) * sj2);
                            IkReal x493 = ((8.5e-7) * x486);
                            IkReal x494 = (pz * x486);
                            IkReal x495 = (pz * x485);
                            IkReal x496 = ((0.26) * x485);
                            IkReal x497 = ((1.92e-6) * x486);
                            evalcond[0] =
                                ((-0.13) + ((x485 * x488)) + (((-1.0) * x492)) +
                                 (((-1.0) * x485 * x490)) +
                                 (((-0.14999) * x485)) + x487 + x494);
                            evalcond[1] =
                                ((9.6e-7) + (((-1.0) * x486 * x488)) + x489 +
                                 x491 + x495 + ((x486 * x490)) +
                                 (((0.14999) * x486)));
                            evalcond[2] =
                                (((x485 * x489)) + (((-1.0) * x486 * x492)) +
                                 (((-0.13) * x486)) + pz + (((9.6e-7) * x485)) +
                                 ((x486 * x487)) + ((x485 * x491)));
                            evalcond[3] =
                                ((-0.14999) + (((-1.0) * x486 * x489)) +
                                 ((x485 * x487)) + (((-9.6e-7) * x486)) +
                                 (((-1.0) * x486 * x491)) + (((-0.13) * x485)) +
                                 (((-1.0) * x490)) + x488 +
                                 (((-1.0) * x485 * x492)));
                            evalcond[4] =
                                ((0.0172892498998009) +
                                 (((-2.879808e-7) * x486)) + ((x488 * x496)) +
                                 ((x488 * x497)) + (((-1.0) * (px * px))) +
                                 (((0.26) * x494)) + (((-0.0389974) * x485)) +
                                 (((-1.92e-6) * x495)) +
                                 (((-1.0) * x490 * x496)) +
                                 (((-1.0) * x490 * x497)) +
                                 (((-1.0) * (pz * pz))) +
                                 (((-0.29998) * x490)) +
                                 (((-1.0) * (py * py))) + (((0.29998) * x488)));
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
IKFAST_API bool el_mini_ComputeIK_rb(const IkReal *eetrans, const IkReal *eerot,
                                     const IkReal *pfree,
                                     IkSolutionListBase<IkReal> &solutions) {
  IKSolver_rb solver;
  return solver.ComputeIk(eetrans, eerot, pfree, solutions);
}

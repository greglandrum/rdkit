//
//  Copyright (c) 2016, Guillaume GODIN.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Institue of Cancer Research.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Guillaume GODIN access the AutoCorrelation 3D descriptors names in Dragon TDB

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

#include "AUTOCORR3D.h"
#include "MolData3Ddescriptors.h"

#include "GraphMol/PartialCharges/GasteigerCharges.h"
#include "GraphMol/PartialCharges/GasteigerParams.h"
#include <Numerics/EigenSolvers/PowerEigenSolver.h>

#include <Numerics/Matrix.h>
#include <Numerics/SquareMatrix.h>
#include <Numerics/SymmMatrix.h>
#include <boost/foreach.hpp>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/QR>

using namespace Eigen;
namespace RDKit {
    namespace Descriptors {

        namespace {

            MolData3Ddescriptors moldata3D;

            VectorXd getEigenVect(std::vector<double> v) {
                double* varray_ptr = &v[0];
                Map<VectorXd> V(varray_ptr, v.size());
                return V;
            }

            double* GetGeodesicMatrix(double* dist, int lag, int numAtoms) {
                int sizeArray = numAtoms * numAtoms;
                double* Geodesic = new double[sizeArray];
                for (int i = 0; i < sizeArray; i++) {
                    if (dist[i] == lag)
                        Geodesic[i] = 1.0;
                    else
                        Geodesic[i] = 0.0;
                }
                return Geodesic;
            }

            // matrix prod that mimic the v2 version
            void get3DautocorrelationDesc(double* dist3D, double* topologicaldistance, int numAtoms,
                                          const ROMol& mol, std::vector<double>& res) {
                Map<MatrixXd> dm(dist3D, numAtoms, numAtoms);
                Map<MatrixXd> di(topologicaldistance, numAtoms, numAtoms);

                std::vector<double> wp = moldata3D.GetRelativePol(mol);
                VectorXd Wp = getEigenVect(wp);

                std::vector<double> wm = moldata3D.GetRelativeMW(mol);
                VectorXd Wm = getEigenVect(wm);

                std::vector<double> wi = moldata3D.GetRelativeIonPol(mol);
                VectorXd Wi = getEigenVect(wi);

                std::vector<double> wv = moldata3D.GetRelativeVdW(mol);
                VectorXd Wv = getEigenVect(wv);

                std::vector<double> we = moldata3D.GetRelativeENeg(mol);
                VectorXd We = getEigenVect(we);

                std::vector<double> wu = moldata3D.GetUn(numAtoms);

                //std::cout << "Topo Distance Matrix : " << numAtoms << "\n";
                //std::cout << di << "\n";
                //std::cout << "\n";

                VectorXd Wu = getEigenVect(wu);

                std::vector<double> ws = moldata3D.GetIState(mol);
                VectorXd Ws = getEigenVect(ws);

                std::vector<double> wr = moldata3D.GetRelativeRcov(mol);
                VectorXd Wr = getEigenVect(wr);

                MatrixXd Bi;
                MatrixXd tmp;
                double TDBmat[8][10];
                double dtmp;

                for (int i = 0; i < 10; i++) {
                    double* Bimat = GetGeodesicMatrix(topologicaldistance, i + 1, numAtoms);
                    Map<MatrixXd> Bi(Bimat, numAtoms, numAtoms);
                    MatrixXd RBi = Bi.cwiseProduct(dm);
                    double Bicount = (double)Bi.sum();

                    tmp = Wu.transpose() * RBi * Wu / Bicount;
                    dtmp = (double)tmp(0);
                    if (std::isnan(dtmp)) dtmp = 0.0;
                    TDBmat[0][i] = dtmp;

                    tmp = Wm.transpose() * RBi * Wm / Bicount;
                    dtmp = (double)tmp(0);
                    if (std::isnan(dtmp)) dtmp = 0.0;
                    TDBmat[1][i] = dtmp;

                    tmp = Wv.transpose() * RBi * Wv / Bicount;
                    dtmp = (double)tmp(0);
                    if (std::isnan(dtmp)) dtmp = 0.0;
                    TDBmat[2][i] = dtmp;

                    tmp = We.transpose() * RBi * We / Bicount;
                    dtmp = (double)tmp(0);
                    if (std::isnan(dtmp)) dtmp = 0.0;
                    TDBmat[3][i] = dtmp;

                    tmp = Wp.transpose() * RBi * Wp / Bicount;
                    dtmp = (double)tmp(0);
                    if (std::isnan(dtmp)) dtmp = 0.0;
                    TDBmat[4][i] = dtmp;

                    tmp = Wi.transpose() * RBi * Wi / Bicount;
                    dtmp = (double)tmp(0);
                    if (std::isnan(dtmp)) dtmp = 0.0;
                    TDBmat[5][i] = dtmp;

                    tmp = Ws.transpose() * RBi * Ws / Bicount;
                    dtmp = (double)tmp(0);
                    if (std::isnan(dtmp)) dtmp = 0.0;
                    TDBmat[6][i] = dtmp;

                    tmp = Wr.transpose() * RBi * Wr / Bicount;
                    dtmp = (double)tmp(0);
                    if (std::isnan(dtmp)) dtmp = 0.0;
                    TDBmat[7][i] = dtmp;
                }

                // update the Output vector!
                for (unsigned int j = 0; j < 8; ++j) {
                    for (unsigned int i = 0; i < 10; ++i) {
                        res[j * 10 + i] = TDBmat[j][i];
                    }
                }
            }

            // same as in the publication
            // J. Chem. Inf. Compunt. Sci. 2004, 44, 200-209 _--(Klein, Kaiser, Ecker)
            // identical results with Padel on the same dataset
            void get3DautocorrelationDesc2(double* dist3D, double* dist, int numAtoms,
                                          const ROMol& mol, std::vector<double>& res) {


                std::vector<double> wp = moldata3D.GetRelativePol(mol);
                std::vector<double> wm = moldata3D.GetRelativeMW(mol);
                std::vector<double> wv = moldata3D.GetRelativeVdW(mol);
                std::vector<double> wi = moldata3D.GetRelativeIonPol(mol);
                std::vector<double> we = moldata3D.GetRelativeENeg(mol);
                std::vector<double> wu = moldata3D.GetUn(numAtoms);
                std::vector<double> ws = moldata3D.GetIState(mol);
                std::vector<double> wr = moldata3D.GetRelativeRcov(mol);
                double w[8][numAtoms];

                for (unsigned int i = 0; i < numAtoms; i++) {
                    w[0][i] = wu[i];
                    w[1][i] = wm[i];
                    w[2][i] = wv[i];
                    w[3][i] = we[i];
                    w[4][i] = wp[i];
                    w[5][i] = wi[i];
                    w[6][i] = ws[i];
                    w[7][i] = wr[i];
                }

                double TDBmat[8][10];
                for (unsigned int k = 0; k < 10; k++) {
                    int maxkVertexPairs = 0;
                    for (unsigned int i = 0; i < numAtoms - 1 ; i++)
                    {
                        for (unsigned int j = i + 1; j < numAtoms; j++)
                        {
                            if (dist[j * numAtoms + i]==k + 1)
                            {
                                for (unsigned int t = 0; t < 8; t++)
                                {
                                    TDBmat[t][k] += w[t][i] * dist3D[j * numAtoms + i] * w[t][j];
                                    //TDBmat[t][k] += dist3D[j * numAtoms + i];

                                }
                                maxkVertexPairs += 1;
                            }
                        }
                    }

                    for (unsigned int t = 0; t < 8; t++)
                    {
                        if (maxkVertexPairs>0) {
                            TDBmat[t][k] = TDBmat[t][k] / maxkVertexPairs;
                        }
                        else {
                            TDBmat[t][k] =  0;
                        }
                    }
                }

                // update the Output vector!
                for (unsigned int j = 0; j < 8; ++j) {
                    for (unsigned int i = 0; i < 10; ++i) {
                        //std::cout << TDBmat[j][i] << ",";
                        res[j * 10 + i] = TDBmat[j][i];
                    }
                }
            }

            // same as in the article Descriptors from Molecular Geometry (Todeschini & Consonni) equation 50
            // there is no link between to two
            void get3DautocorrelationDesc3(double* dist3D, double* dist, int numAtoms,
                                           const ROMol& mol, std::vector<double>& res) {


                std::vector<double> wp = moldata3D.GetRelativePol(mol);
                std::vector<double> wm = moldata3D.GetRelativeMW(mol);
                std::vector<double> wv = moldata3D.GetRelativeVdW(mol);
                std::vector<double> wi = moldata3D.GetRelativeIonPol(mol);
                std::vector<double> we = moldata3D.GetRelativeENeg(mol);
                std::vector<double> wu = moldata3D.GetUn(numAtoms);
                std::vector<double> ws = moldata3D.GetIState(mol);
                std::vector<double> wr = moldata3D.GetRelativeRcov(mol);
                double w[8][numAtoms];

                for (unsigned int i = 0; i < numAtoms; i++) {
                    w[0][i] = wu[i];
                    w[1][i] = wm[i];
                    w[2][i] = wv[i];
                    w[3][i] = we[i];
                    w[4][i] = wp[i];
                    w[5][i] = wi[i];
                    w[6][i] = ws[i];
                    w[7][i] = wr[i];
                }

                double TDBmat[8][10];
                for (unsigned int k = 0; k < 10; k++) {
                    int maxkVertexPairs = 0;
                    for (unsigned int i = 0; i < numAtoms - 1 ; ++i)
                    {
                        for (unsigned int j = i + 1; j < numAtoms; ++j)
                        {
                            //std::cout << "test distance matrix" << dist[j * numAtoms + i] << std::endl;
                            if (dist3D[j * numAtoms + i]>=k + 1 and dist3D[j * numAtoms + i]>k)
                            {
                                for (unsigned int t = 0; t < 8; ++t)
                                {
                                    //TDBmat[t][k] += w[t][i] * dist3D[j * numAtoms + i] * w[t][j];
                                    TDBmat[t][k] += w[t][i] * w[t][j];

                                    //TDBmat[t][k] += dist3D[j * numAtoms + i];

                                }
                                maxkVertexPairs += 1;
                            }
                        }
                    }

                    for (unsigned int t = 0; t < 8; ++t)
                    {
                        if (maxkVertexPairs>0) {
                            TDBmat[t][k] = TDBmat[t][k] / maxkVertexPairs;
                        }
                        else {
                            TDBmat[t][k] =  0.0;
                        }
                    }
                }

                // update the Output vector!
                for (unsigned int j = 0; j < 8; ++j) {
                    for (unsigned int i = 0; i < 10; ++i) {
                        res[j * 10 + i] = TDBmat[j][i];
                    }
                }
            }





            void Get3Dauto(double* dist3D, double* topologicaldistance, int numAtoms, const ROMol& mol,
                           std::vector<double>& res) {
                // std::vector<std::string>
                // AUTOCORRNAMES={"TDB01u","TDB02u","TDB03u","TDB04u","TDB05u","TDB06u","TDB07u","TDB08u","TDB09u","TDB10u","TDB01m","TDB02m","TDB03m","TDB04m","TDB05m","TDB06m","TDB07m","TDB08m","TDB09m","TDB10m","TDB01v","TDB02v","TDB03v","TDB04v","TDB05v","TDB06v","TDB07v","TDB08v","TDB09v","TDB10v","TDB01e","TDB02e","TDB03e","TDB04e","TDB05e","TDB06e","TDB07e","TDB08e","TDB09e","TDB10e","TDB01p","TDB02p","TDB03p","TDB04p","TDB05p","TDB06p","TDB07p","TDB08p","TDB09p","TDB10p","TDB01i","TDB02i","TDB03i","TDB04i","TDB05i","TDB06i","TDB07i","TDB08i","TDB09i","TDB10i","TDB01s","TDB02s","TDB03s","TDB04s","TDB05s","TDB06s","TDB07s","TDB08s","TDB09s","TDB10s","TDB01r","TDB02r","TDB03r","TDB04r","TDB05r","TDB06r","TDB07r","TDB08r","TDB09r","TDB10r"};
                get3DautocorrelationDesc(dist3D, topologicaldistance, numAtoms, mol, res);
            }

        }  // end of anonymous namespace

        void AUTOCORR3D(const ROMol& mol, std::vector<double>& res, int confId) {
            PRECONDITION(mol.getNumConformers() >= 1, "molecule has no conformers")
            int numAtoms = mol.getNumAtoms();

            const Conformer& conf = mol.getConformer(confId);

            double* topologicaldistance = MolOps::getDistanceMat(mol, false);  // topological matrix
            double* dist3D = MolOps::get3DDistanceMat(mol, confId, false, true);  // 3D distance matrix

            res.clear();
            res.resize(80);
            Get3Dauto(dist3D, topologicaldistance, numAtoms, mol, res);
        }
    }  // end of Descriptors namespace
}  // end of RDKit namespace

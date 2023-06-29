// ================================================================
//
// This software was developed by employees of the National Institute of
// Standards and Technology (NIST), an agency of the Federal Government.
// Pursuant to title 17 United States Code Section 105, works of NIST employees
// are not subject to copyright protection in the United States and are
// considered to be in the public domain. Permission to freely use, copy,
// modify, and distribute this software and its documentation without fee is
// hereby granted, provided that this notice and disclaimer of warranty appears
// in all copies.
//
// THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
// EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY
// WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM
// FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO
// THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO
// EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO,
// DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF,
// RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT
// BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS
// SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS
// SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR
// SERVICES PROVIDED HEREUNDER.
//
// Distributions of NIST software should also include copyright and licensing
// statements of any third-party software that are legally bundled with the
// code in compliance with the conditions of those licenses.
// 
// ================================================================

// ================================================================
// 
// Authors: Derek Juba <derek.juba@nist.gov>
// Created: Wed Apr 22 11:11:48 2015 EDT
//
// ================================================================

#include "ResultsZeno.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace zeno;

/// Constructs the class to collect results with the given bounding sphere 
/// from the given number of threads,
/// and optionally save the hit point locations.
///
ResultsZeno::
ResultsZeno(Sphere<double> const & boundingSphere,
	    int numThreads,
	    bool saveHitPoints) 
  : boundingSphere(boundingSphere),
    numThreads(numThreads),
    saveHitPoints(saveHitPoints),
    numWalks(NULL),
    hitMissMean(NULL),
    hitMissM2(NULL),
    KPlus(NULL),
    KMinus(NULL),
    KPlusMean(NULL),
    KMinusMean(NULL),
    KPlusM2(NULL),
    KMinusM2(NULL),
    VPlus(NULL),
    VMinus(NULL),
    VPlusMean(NULL),
    VMinusMean(NULL),
    VPlusM2(NULL),
    VMinusM2(NULL),
    numWalksReduced(0),
    numHitsReduced(0),
    numHitsVarianceReduced(0),
    KPlusReduced(0, 0, 0),
    KMinusReduced(0, 0, 0),
    VPlusReduced(0, 0, 0, 
		 0, 0, 0, 
		 0, 0, 0),
    VMinusReduced(0, 0, 0, 
		  0, 0, 0, 
		  0, 0, 0),
    points(NULL),
    charges(NULL),
    gatheredPoints(),
    gatheredCharges(),
    reduced(true),
    hitPointsGathered(true),
    totalSteps(NULL), //Added by mvk1-nist
    hitSteps(NULL), // Added by mvk1-nist
    missSteps(NULL), // Added by mvk1-nist
    totalStepsMean(NULL), // Added by mvk1-nist
    totalStepsM2(NULL), // Added by mvk1-nist
    hitStepsMean(NULL), // Added by mvk1-nist
    hitStepsM2(NULL), // Added by mvk1-nist
    missStepsMean(NULL), // Added by mvk1-nist
    missStepsM2(NULL), // Added by mvk1-nist
    totalStepsReduced(0), // Added by mvk1-nist
    totalStepsVarianceReduced(0), // Added by mvk1-nist
    hitStepsReduced(0), // Added by mvk1-nist
    hitStepsVarianceReduced(0), // Added by mvk1-nist
    missStepsReduced(0), // Added by mvk1-nist
    missStepsVarianceReduced(0) { // Added by mvk1-nist

  totalSteps = new double[numThreads]; // Added by mvk1-nist
  hitSteps = new double[numThreads]; // Added by mvk1-nist
  missSteps = new double[numThreads]; // Added by mvk1-nist

  totalStepsMean = new double[numThreads]; // Added by mvk1-nist
  totalStepsM2 = new double[numThreads]; // Added by mvk1-nist
  hitStepsMean = new double[numThreads]; // Added by mvk1-nist
  hitStepsM2 = new double[numThreads]; // Added by mvk1-nist
  missStepsMean = new double[numThreads]; // Added by mvk1-nist
  missStepsM2 = new double[numThreads]; // Added by mvk1-nist


  numWalks = new double[numThreads];

  hitMissMean = new double[numThreads];

  hitMissM2 = new double[numThreads];

  KPlus  = new Vector3<double>[numThreads];
  KMinus = new Vector3<double>[numThreads];

  KPlusMean  = new Vector3<double>[numThreads];
  KMinusMean = new Vector3<double>[numThreads];  

  KPlusM2  = new Vector3<double>[numThreads];
  KMinusM2 = new Vector3<double>[numThreads];

  VPlus  = new Matrix3x3<double>[numThreads];
  VMinus = new Matrix3x3<double>[numThreads];

  VPlusMean  = new Matrix3x3<double>[numThreads];
  VMinusMean = new Matrix3x3<double>[numThreads];

  VPlusM2  = new Matrix3x3<double>[numThreads];
  VMinusM2 = new Matrix3x3<double>[numThreads];

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    totalSteps[threadNum] = 0; // Added by mvk1-nist
    hitSteps[threadNum] = 0; // Added by mvk1-nist
    missSteps[threadNum] = 0; // Added by mvk1-nist

    totalStepsMean[threadNum] = 0; // Added by mvk1-nist
    totalStepsM2[threadNum] = 0; // Added by mvk1-nist

    hitStepsMean[threadNum] = 0; // Added by mvk1-nist
    hitStepsM2[threadNum] = 0; // Added by mvk1-nist

    missStepsMean[threadNum] = 0; // Added by mvk1-nist
    missStepsM2[threadNum] = 0; // Added by mvk1-nist

    numWalks[threadNum] = 0;

    hitMissMean[threadNum] = 0;

    hitMissM2[threadNum] = 0;

    KPlus[threadNum].setXYZ(0, 0, 0);
    KMinus[threadNum].setXYZ(0, 0, 0);

    KPlusMean[threadNum].setXYZ(0, 0, 0);
    KMinusMean[threadNum].setXYZ(0, 0, 0);

    KPlusM2[threadNum].setXYZ(0, 0, 0);
    KMinusM2[threadNum].setXYZ(0, 0, 0);

    for (int component = 0; component < 3*3; component++) {
      VPlus[threadNum].set(component, 0);
      VMinus[threadNum].set(component, 0);

      VPlusMean[threadNum].set(component, 0);
      VMinusMean[threadNum].set(component, 0);

      VPlusM2[threadNum].set(component, 0);
      VMinusM2[threadNum].set(component, 0);
    }
  }

  points  = new std::vector<Vector3<double> >[numThreads];
  charges = new std::vector<Vector3<char> >[numThreads];
}

ResultsZeno::
~ResultsZeno() {
  delete [] totalSteps; // Added by mvk1-nist
  delete [] hitSteps; // Added by mvk1-nist
  delete [] missSteps; // Added by mvk1-nist

  delete [] numWalks;

  delete [] hitMissMean;

  delete [] hitMissM2;

  delete [] KPlus;
  delete [] KMinus;

  delete [] KPlusMean;
  delete [] KMinusMean;

  delete [] KPlusM2;
  delete [] KMinusM2;

  delete [] VPlus;
  delete [] VMinus;

  delete [] VPlusMean;
  delete [] VMinusMean;

  delete [] VPlusM2;
  delete [] VMinusM2;

  delete [] points;
  delete [] charges;
}

/// Record a miss from the given thread number.
///
void 
ResultsZeno::
recordMiss(int threadNum, double missStepCount) { // Modified by mvk1-nist
  assert(threadNum >= 0 && threadNum < numThreads);

  reduced = false;

  double hitMissData = 0;

  Vector3<double> KPlusData(0, 0, 0);
  Vector3<double> KMinusData(0, 0, 0);

  Matrix3x3<double> VPlusData(0, 0, 0, 0, 0, 0, 0, 0, 0);
  Matrix3x3<double> VMinusData(0, 0, 0, 0, 0, 0, 0, 0, 0);

  totalSteps[threadNum] += missStepCount; // Added by mvk1-nist
  missSteps[threadNum] += missStepCount; // Added by mvk1-nist

  numWalks[threadNum]++;

  updateVariance(threadNum,
		 hitMissData,
		 KPlusData, 
		 KMinusData,
		 VPlusData, 
		 VMinusData,
		 missStepCount, // Modified by mvk1-nist
		 false); // Modified by mvk1-nist
}

/// Perform a parallel reduction on the hit counts and other statistics and 
/// corresponding variances across threads and MPI nodes.
///
void 
ResultsZeno::
reduce() {
  if (reduced) {
    return;
  }

  numWalksReduced = 0;

  numHitsReduced = 0;

  numHitsVarianceReduced = 0;

  KPlusReduced.setXYZ(0, 0, 0);
  KMinusReduced.setXYZ(0, 0, 0);

  KPlusVarianceReduced.setXYZ(0, 0, 0);
  KMinusVarianceReduced.setXYZ(0, 0, 0);

  for (int component = 0; component < 3*3; component++) {
    VPlusReduced.set(component, 0);
    VMinusReduced.set(component, 0);

    VPlusVarianceReduced.set(component, 0);
    VMinusVarianceReduced.set(component, 0);
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    const double nn1 = (double)numWalks[threadNum] / (numWalks[threadNum] - 1);

    numWalksReduced += numWalks[threadNum];

    numHitsReduced += hitMissMean[threadNum] * numWalks[threadNum];

    numHitsVarianceReduced += hitMissM2[threadNum] * nn1;

    KPlusReduced  += KPlus[threadNum];
    KMinusReduced += KMinus[threadNum];

    KPlusVarianceReduced  += KPlusM2[threadNum] * nn1;
    KMinusVarianceReduced += KMinusM2[threadNum] * nn1;

    VPlusReduced  += VPlus[threadNum];
    VMinusReduced += VMinus[threadNum];

    VPlusVarianceReduced  += VPlusM2[threadNum] * nn1;
    VMinusVarianceReduced += VMinusM2[threadNum] * nn1;

    totalStepsReduced += totalSteps[threadNum]; // Added by mvk1-nist
    totalStepsVarianceReduced += totalStepsM2[threadNum] * nn1; // Added by mvk1-nist

    hitStepsReduced += hitSteps[threadNum]; // Added by mvk1-nist
    hitStepsVarianceReduced += hitStepsM2[threadNum] * nn1; // Added by mvk1-nist

    missStepsReduced += missSteps[threadNum]; // Added by mvk1-nist
    missStepsVarianceReduced += missStepsM2[threadNum] * nn1; // Added by mvk1-nist
  }

#ifdef USE_MPI
  const int mpiBufferSize = 51;

  double sendbuf[mpiBufferSize];

  int offset = 0;

  sendbuf[offset++] = numWalksReduced;

  sendbuf[offset++] = numHitsReduced;
  sendbuf[offset++] = numHitsVarianceReduced;

  for (int i = 0; i < 3; i++) {
    sendbuf[offset++] = KPlusReduced.get(i);
    sendbuf[offset++] = KPlusVarianceReduced.get(i);
  }

  for (int i = 0; i < 3; i++) {
    sendbuf[offset++] = KMinusReduced.get(i);
    sendbuf[offset++] = KMinusVarianceReduced.get(i);
  }

  for (int i = 0; i < 9; i++) {
    sendbuf[offset++] = VPlusReduced.get(i);
    sendbuf[offset++] = VPlusVarianceReduced.get(i);
  }

  for (int i = 0; i < 9; i++) {
    sendbuf[offset++] = VMinusReduced.get(i);
    sendbuf[offset++] = VMinusVarianceReduced.get(i);
  }

  double recvbuf[mpiBufferSize];

  for (int i = 0; i < mpiBufferSize; i++) {
    recvbuf[i] = 0;
  }

  // MPI_Reduce(sendbuf, recvbuf, mpiBufferSize, MPI_DOUBLE,
  // 	     MPI_SUM, 0, MPI_COMM_WORLD);

  MPI_Allreduce(sendbuf, recvbuf, mpiBufferSize, MPI_DOUBLE,
		MPI_SUM, MPI_COMM_WORLD);

  offset = 0;

  numWalksReduced = recvbuf[offset++];

  numHitsReduced         = recvbuf[offset++];
  numHitsVarianceReduced = recvbuf[offset++];

  for (int i = 0; i < 3; i++) {
    KPlusReduced.set(i, recvbuf[offset++]);
    KPlusVarianceReduced.set(i, recvbuf[offset++]);
  }

  for (int i = 0; i < 3; i++) {
    KMinusReduced.set(i, recvbuf[offset++]);
    KMinusVarianceReduced.set(i, recvbuf[offset++]);
  }

  for (int i = 0; i < 9; i++) {
    VPlusReduced.set(i, recvbuf[offset++]);
    VPlusVarianceReduced.set(i, recvbuf[offset++]);
  }

  for (int i = 0; i < 9; i++) {
    VMinusReduced.set(i, recvbuf[offset++]);
    VMinusVarianceReduced.set(i, recvbuf[offset++]);
  }
#endif

  reduced = true;
}

/// Gather the hit locations from all threads and MPI nodes.
///
void 
ResultsZeno::
gatherHitPoints() {
  if (hitPointsGathered) {
    return;
  }

  for (int threadNum = 0; threadNum < numThreads; threadNum++) {
    gatheredPoints.insert(gatheredPoints.end(), 
			 points[threadNum].begin(), 
			 points[threadNum].end());

    gatheredCharges.insert(gatheredCharges.end(),
			  charges[threadNum].begin(),
			  charges[threadNum].end());
  }

#ifdef USE_MPI
  int mpiSize = 0;
  int mpiRank = 0;

  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

  int * arrayLengths = new int[mpiSize];

  int arrayLength = gatheredPoints.size() * 3;

  MPI_Allgather(&arrayLength, 1, MPI_INT,
		arrayLengths, 1, MPI_INT,
		MPI_COMM_WORLD);

  int combinedArrayLength = 0;

  for (int i = 0; i < mpiSize; i++) {
    combinedArrayLength += arrayLengths[i];
  }

  double * combinedPointArray  = new double[combinedArrayLength];
  char *   combinedChargeArray = new char[combinedArrayLength];

  int * combinedArrayOffsets = new int[mpiSize];

  combinedArrayOffsets[0] = 0;

  for (int i = 1; i < mpiSize; i++) {
    combinedArrayOffsets[i] = combinedArrayOffsets[i - 1] + arrayLengths[i - 1];
  }

  double * pointArray  = new double[arrayLengths[mpiRank]];
  char *   chargeArray = new char[arrayLengths[mpiRank]];

  for (int index = 0; index < arrayLengths[mpiRank]; index++) {
    pointArray[index]  = gatheredPoints[index / 3].get(index % 3);
    chargeArray[index] = gatheredCharges[index / 3].get(index % 3);
  }

  MPI_Allgatherv(pointArray, 
		 arrayLengths[mpiRank], MPI_DOUBLE,
		 combinedPointArray, 
		 arrayLengths, combinedArrayOffsets, MPI_DOUBLE,
		 MPI_COMM_WORLD);

  MPI_Allgatherv(chargeArray,
		 arrayLengths[mpiRank], MPI_BYTE,
		 combinedChargeArray,
		 arrayLengths, combinedArrayOffsets, MPI_BYTE,
		 MPI_COMM_WORLD);

  gatheredPoints.clear();
  gatheredCharges.clear();

  gatheredPoints.reserve(combinedArrayLength / 3);
  gatheredCharges.reserve(combinedArrayLength / 3);

  for (int index = 0; index < combinedArrayLength; index += 3) {
    gatheredPoints.emplace_back(combinedPointArray[index + 0],
				combinedPointArray[index + 1],
				combinedPointArray[index + 2]);

    gatheredCharges.emplace_back(combinedChargeArray[index + 0],
				 combinedChargeArray[index + 1],
				 combinedChargeArray[index + 2]);
  }

  delete [] arrayLengths;

  delete [] combinedPointArray;
  delete [] combinedChargeArray;

  delete [] combinedArrayOffsets;

  delete [] pointArray;
  delete [] chargeArray;
#endif

  hitPointsGathered = true;
}

void 
ResultsZeno::
updateVariance(int threadNum,
	       double hitMissData,
	       Vector3<double> const & KPlusData, 
	       Vector3<double> const & KMinusData,
	       Matrix3x3<double> const & VPlusData, 
	       Matrix3x3<double> const & VMinusData,
	       double stepData, // Modified by mvk1-nist
	       bool hitNotMiss) { // Modified by mvk1-nist

  updateItemVariance(hitMissData,
		     numWalks[threadNum],
		     &(hitMissMean[threadNum]),
		     &(hitMissM2[threadNum]));

  updateItemVariance(KPlusData,
		     numWalks[threadNum],
		     &(KPlusMean[threadNum]),
		     &(KPlusM2[threadNum]));

  updateItemVariance(KMinusData,
		     numWalks[threadNum],
		     &(KMinusMean[threadNum]),
		     &(KMinusM2[threadNum]));

  updateItemVariance(VPlusData,
		     numWalks[threadNum],
		     &(VPlusMean[threadNum]),
		     &(VPlusM2[threadNum]));

  updateItemVariance(VMinusData,
		     numWalks[threadNum],
		     &(VMinusMean[threadNum]),
		     &(VMinusM2[threadNum]));

  updateItemVariance(stepData,
             numWalks[threadNum],
             &(totalStepsMean[threadNum]),
             &(totalStepsM2[threadNum])); // Added by mvk1-nist

  if(hitNotMiss) {
      updateItemVariance(stepData,
                 numWalks[threadNum],
                 &(hitStepsMean[threadNum]),
                 &(hitStepsM2[threadNum])); // Added by mvk1-nist
  } else {
      updateItemVariance(stepData,
                 numWalks[threadNum],
                 &(missStepsMean[threadNum]),
                 &(missStepsM2[threadNum])); // Added by mvk1-nist
  } // Added by mvk1-nist



}

Uncertain<double> 
ResultsZeno::
getTotalSteps() const {
  assert(reduced);

  return Uncertain<double>(totalStepsReduced, totalStepsVarianceReduced);
} // Added by mvk1-nist

Uncertain<double> 
ResultsZeno::
getHitSteps() const {
  assert(reduced);

  return Uncertain<double>(hitStepsReduced, hitStepsVarianceReduced);
} // Added by mvk1-nist

Uncertain<double> 
ResultsZeno::
getMissSteps() const {
  assert(reduced);

  return Uncertain<double>(missStepsReduced, missStepsVarianceReduced);
} // Added by mvk1-nist

double 
ResultsZeno::
getNumWalks() const {
  assert(reduced);

  return numWalksReduced;
}

Uncertain<double> 
ResultsZeno::
getNumHits() const {
  assert(reduced);

  return Uncertain<double>(numHitsReduced, numHitsVarianceReduced);
}

Vector3<Uncertain<double> > 
ResultsZeno::
getKPlus() const {
  assert(reduced);

  return Uncertain<double>::zip(KPlusReduced, KPlusVarianceReduced);
}

Vector3<Uncertain<double> > 
ResultsZeno::
getKMinus() const {
  assert(reduced);

  return Uncertain<double>::zip(KMinusReduced, KMinusVarianceReduced);
}

Matrix3x3<Uncertain<double> > 
ResultsZeno::
getVPlus() const {
  assert(reduced);

  return Uncertain<double>::zip(VPlusReduced, VPlusVarianceReduced);
}

Matrix3x3<Uncertain<double> > 
ResultsZeno::
getVMinus() const {
  assert(reduced);

  return Uncertain<double>::zip(VMinusReduced, VMinusVarianceReduced);
}

bool
ResultsZeno::
getSaveHitPoints() const {

  return saveHitPoints;
}

Sphere<double>
ResultsZeno::
getBoundingSphere() const {
  return boundingSphere;
}

std::vector<Vector3<double> > const * 
ResultsZeno::
getPoints() const {
  assert(saveHitPoints);
  assert(hitPointsGathered);

  return &gatheredPoints;
}

std::vector<Vector3<char> > const * 
ResultsZeno:: 
getCharges() const {
  assert(saveHitPoints);
  assert(hitPointsGathered);

  return &gatheredCharges;
}


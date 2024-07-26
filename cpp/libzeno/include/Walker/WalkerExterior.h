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
// Created: Fri Feb 13 13:31:22 2015 EDT
//
// ================================================================

#ifndef WALKER_EXTERIOR_H
#define WALKER_EXTERIOR_H

#include <vector>
#include <cmath>

#include "../Geometry/Sphere.h"
#include "../Geometry/Vector3.h"

namespace zeno {

/// Performs random walks starting on a bounding sphere and determines whether 
/// they hit an object, allowing for a given relative error in distance.
///
template <class T, 
  class RandomNumberGenerator,
  class NearestSurfacePointFinder,
  class RandomSpherePointGenerator,
  class BiasedSpherePointGenerator>
class WalkerExterior {
 public:
  WalkerExterior(RandomNumberGenerator * randomNumberGenerator, 
		 Sphere<T> const & boundingSphere, 
		 NearestSurfacePointFinder const & nearestSurfacePointFinder,
		 T shellThickness);

  ~WalkerExterior();

  void walk(bool * hitObject, int * numSteps,
	    Vector3<T> * startPoint, Vector3<T> * endPoint,
        int expansion); // Modified by mvk1-nist

 private:
  RandomNumberGenerator * randomNumberGenerator;
  Sphere<T> const * boundingSphere;
  NearestSurfacePointFinder const * nearestSurfacePointFinder;
  T shellThickness;
};

template <class T, 
  class RandomNumberGenerator,
  class NearestSurfacePointFinder,
  class RandomSpherePointGenerator,
  class BiasedSpherePointGenerator>
WalkerExterior<T, 
               RandomNumberGenerator,
               NearestSurfacePointFinder,
               RandomSpherePointGenerator,
               BiasedSpherePointGenerator>::
  WalkerExterior(RandomNumberGenerator * randomNumberGenerator, 
		 Sphere<T> const & boundingSphere, 
		 NearestSurfacePointFinder const & nearestSurfacePointFinder,
		 T shellThickness) :
  randomNumberGenerator(randomNumberGenerator), 
  boundingSphere(&boundingSphere),
  nearestSurfacePointFinder(&nearestSurfacePointFinder),
  shellThickness(shellThickness) {

}

template <class T, 
  class RandomNumberGenerator,
  class NearestSurfacePointFinder,
  class RandomSpherePointGenerator,
  class BiasedSpherePointGenerator>
WalkerExterior<T, 
               RandomNumberGenerator,
               NearestSurfacePointFinder,
               RandomSpherePointGenerator,
               BiasedSpherePointGenerator>::
  ~WalkerExterior() {

}

/// Perform a random walk and determine whether it hits the object, the number
/// of steps it took, its start and end points, and the surface normal of its
/// hit point.
///
template <class T, 
  class RandomNumberGenerator,
  class NearestSurfacePointFinder,
  class RandomSpherePointGenerator,
  class BiasedSpherePointGenerator>
void 
WalkerExterior<T, 
               RandomNumberGenerator,
               NearestSurfacePointFinder,
               RandomSpherePointGenerator,
               BiasedSpherePointGenerator>::
  walk(bool * hitObject, int * numSteps,
       Vector3<T> * startPoint, Vector3<T> * endPoint,
       int expansion) { // Modified by mvk1-nist

  *hitObject = false;
  *numSteps  = 0;

  Vector3<T> position = 
    RandomSpherePointGenerator::generate(randomNumberGenerator, 
					 *boundingSphere);

  *startPoint = position;

  for (;;) {

    T minDistanceSqr = 0;

    nearestSurfacePointFinder->findNearestPointSigned(position,
						&minDistanceSqr);

    T minDistance = std::sqrt(minDistanceSqr);
    
    if (minDistance < shellThickness) {
      //walker is absorbed
      
      *endPoint  = position;
      *hitObject = true;
      
      return;
    }

    //Expanding Variable
    T delta; // Added by mvk1-nist
    switch(expansion) {
        case 0: // Original case
            delta = 0;
            break;
        case 1: // Constant case 0.1e
            delta = shellThickness * 0.1;
            break;
        case 2: // Constant case 0.25e
            delta = shellThickness * 0.25;
            break;
        case 3: // Constant case 0.5e
            delta = shellThickness * 0.5;
            break;
        case 4: // Constant case 1e
            delta = shellThickness * 1;
            break;
        case 5: // Constant case 2e
            delta = shellThickness * 2;
            break;
        case 6: // Constant case 3e
            delta = shellThickness * 3;
            break;
        case 7: // Constant case 4e
            delta = shellThickness * 4;
            break;
        case 8: // Constant case 5e
            delta = shellThickness * 5;
            break;
        case 9: // Proportional case 10%
            delta = 0.1 * minDistance;
            break;
        case 10: // Proportional case 5%
            delta = 0.05 * minDistance;
            break;
        case 11: // Proportional case 4%
            delta = 0.04 * minDistance;
            break;
        case 12: // Proportional case 3%
            delta = 0.03 * minDistance;
            break;
        case 13: // Proportional case 2%
            delta = 0.02 * minDistance;
            break;
        case 14: // Proportional case 1%
            delta = 0.01 * minDistance;
            break;
        case 15: // Proportional case 0.5%
            delta = 0.005 * minDistance;
            break;
        case 16: // Proportional case 0.25%
            delta = 0.0025 * minDistance;
            break;
        case 17: // Proportional case 0.1%
            delta = 0.001 * minDistance;
            break;
        case 18: // Random case uniform between -0.5e and +0.5e
            delta = randomNumberGenerator->getRandInRange(-0.5 * shellThickness, 0.5 * shellThickness);
            break;
        case 19: // Random case uniform between -e and +e
            delta = randomNumberGenerator->getRandInRange(-1 * shellThickness, 1 * shellThickness);
            break;
        case 20: // Random case uniform between -2e and +2e
            delta = randomNumberGenerator->getRandInRange(-2 * shellThickness, 2 * shellThickness);
            break;
        case 21: // Random case uniform between 0 and +0.5e
            delta = randomNumberGenerator->getRandInRange(0, 0.5 * shellThickness);
            break;
        case 22: // Random case uniform between 0 and +e
            delta = randomNumberGenerator->getRandInRange(0, 1 * shellThickness);
            break;
        case 23: // Random case uniform between 0 and +2e
            delta = randomNumberGenerator->getRandInRange(0, 2 * shellThickness);
            break;
        case 24: // Random case normal with mean 0 and sd 0.5e
            delta = 0.5 * shellThickness * std::sqrt(-2 * std::log(randomNumberGenerator->getRandInRange(0, 1))) * std::cos(2 * M_PI * randomNumberGenerator->getRandInRange(0, 1));
            break;
        case 25: // Random case normal with mean 0 and sd 1e
            delta = 1 * shellThickness * std::sqrt(-2 * std::log(randomNumberGenerator->getRandInRange(0, 1))) * std::cos(2 * M_PI * randomNumberGenerator->getRandInRange(0, 1));
            break;
        case 26: // Random case normal with mean 0 and sd 2e
            delta = 2 * shellThickness * std::sqrt(-2 * std::log(randomNumberGenerator->getRandInRange(0, 1))) * std::cos(2 * M_PI * randomNumberGenerator->getRandInRange(0, 1));
            break;
        case 27: // Piecewise case (Far: Proportional case 2%; Close: Constant case 1e; Boundary: 4e)
            if (minDistance > 4 * shellThickness) {
                delta = 0.02 * minDistance;
            } else {
                delta = shellThickness * 1;
            }
            break;
        case 28: // Piecewise case (Far: Proportional case 2%; Close: Constant case 1e; Boundary: 16e)
            if (minDistance > 16 * shellThickness) {
                delta = 0.02 * minDistance;
            } else {
                delta = shellThickness * 1;
            }
            break;
        case 29: // Piecewise case (Far: Proportional case 2%; Close: Constant case 1e; Boundary: 64e)
            if (minDistance > 64 * shellThickness) {
                delta = 0.02 * minDistance;
            } else {
                delta = shellThickness * 1;
            }
            break;
        default: // No option
            exit(1);
    } // Added by mvk1-nist

    minDistance += delta;

    (*numSteps)++;

    Sphere<T> stepSphere(position, minDistance);

    position = RandomSpherePointGenerator::generate(randomNumberGenerator, 
						    stepSphere);

    T centerDistSqr = 
      (position - boundingSphere->getCenter()).getMagnitudeSqr();

    if (centerDistSqr > boundingSphere->getRadiusSqr()) {
      //walker left bounding sphere

      T alpha = boundingSphere->getRadius() / std::sqrt(centerDistSqr);

      if (randomNumberGenerator->getRandInRange(0, 1) > (1 - alpha)) {
	//walker is replaced

	position = BiasedSpherePointGenerator::generate(randomNumberGenerator, 
							*boundingSphere,
							position,
							alpha);
      }
      else {
	//walker escapes

	*hitObject = false;
	return;
      }
    }
  }
}

}

#endif


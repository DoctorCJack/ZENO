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
// Created: 2019-08-02
//
// ================================================================

// Generated by Bisonc++ V6.01.00 on Fri, 02 Aug 2019 12:47:24 -0400

#ifndef bod_parserBodParser_h_included
#define bod_parserBodParser_h_included

// $insert baseclass
#include "BodParserbase.h"
// $insert scanner.h
#include "BodScanner.h"

#undef BodParser
    // CAVEAT: between the baseclass-include directive and the 
    // #undef directive in the previous line references to BodParser 
    // are read as BodParserBase.
    // If you need to include additional headers in this file 
    // you should do so after these comment-lines.

#include <vector>
#include <string>

#include "../ParametersWalkOnSpheres.h"
#include "../ParametersInteriorSampling.h"
#include "../ParametersResults.h"
#include "../ParametersLocal.h"
#include "../Geometry/MixedModel.h"

// $insert namespace-open
namespace bod_parser
{

class BodParser: public BodParserBase
{
    // $insert scannerobject
    BodScanner d_scanner;

    zeno::ParametersWalkOnSpheres * parametersWalkOnSpheres;
    zeno::ParametersInteriorSampling * parametersInteriorSampling;
    zeno::ParametersResults * parametersResults;
    ParametersLocal const * parametersLocal;
    zeno::MixedModel<double> * model;

    public:
        BodParser(ParametersLocal const & parametersLocal,
	       std::istream &in,
	       zeno::ParametersWalkOnSpheres * parametersWalkOnSpheres,
	       zeno::ParametersInteriorSampling * parametersInteriorSampling,
	       zeno::ParametersResults * parametersResults,
	       zeno::MixedModel<double> * model);

        int parse();

    private:
        void error();                   // called on (syntax) errors
        int lex();                      // returns the next token from the
                                        // lexical scanner. 
        void print();                   // use, e.g., d_token, d_loc
        void exceptionHandler(std::exception const &exc);

        void addSphere(double x, double y, double z, double r);
	void addCube(double x, double y, double z, double s);
	void addCuboid(double x1, double y1, double z1,
		       double x2, double y2, double z2);
	void addVoxels(std::string voxelsFileName);
	void setST(double skinThickness);
	void setRLAUNCH(double launchRadius);
	void setHUNITS(double number, std::string unitString);
	void setUNITS(std::string unitString);
	void setTEMP(double number, std::string unitString);
	void setMASS(double number, std::string unitString);
	void setVISCOSITY(double number, std::string unitString);
	void setBF(double buoyancyFactor);

    // support functions for parse():
        void executeAction__(int ruleNr);
        void errorRecovery__();
        void nextCycle__();
        void nextToken__();
        void print__();
};

// $insert namespace-close
}

#endif

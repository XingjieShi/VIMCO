//
//  plinkfun.hpp
//  PlinkRead
//
//  Created by DaiMingwei on 16/10/6.
//  Copyright © 2016年 daviddai. All rights reserved.
//




#ifndef plinkfun_hpp
#define plinkfun_hpp



#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <bitset>
//#include <boost/algorithm/string.hpp>
//// [[Rcpp::depends(BH)]] 

using namespace std;
using namespace Rcpp;
using namespace arma;


void getFourGentype(int* geno, std::bitset<8> bits);
vector<string> read_snpnames(string filename, int P);
void readPlink(string stringname, int N, int P, unsigned* X);
int getLineNum(string filename);
void ReadPlinkFamFile3(std::string stringname, CharacterVector FID, CharacterVector IID,
	NumericMatrix pheno, int nrows, int nPheno);


#endif /* plinkfun_hpp */

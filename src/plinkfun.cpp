//
//  plinkfun.cpp
//  PlinkRead
//
//  Created by DaiMingwei on 16/10/6.
//  Copyright © 2016年 daviddai. All rights reserved.
//


// We can now use the BH package
// [[Rcpp::depends(BH, RcppArmadillo)]] 

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "plinkfun.hpp"
#include <stdio.h>
#include <math.h>
#include <bitset>
#include <boost/algorithm/string.hpp>

#define MAX_LEN 20

using namespace std;


int getLineNum(string filename){
    FILE *pf = fopen(filename.c_str(), "r"); // 打开文件
    char buf[10000];
    int lineCnt = 0;
    if (!pf) // 判断是否打开成功
        return -1;
    while (fgets(buf, 10000, pf)) // fgets循环读取，直到文件最后，才会返回NULL
        lineCnt++; // 累计行数
    fclose(pf);
    return lineCnt;
}

void getFourGentype(int* geno, std::bitset<8> bits){
    int idx = 0;
    for (int j=0; j < 8; j = j + 2) {
        if(bits[j] && bits[j+1]){
            geno[idx] = 0;
        }else if(!bits[j] && !bits[j+1]){
            geno[idx] = 2;
        }else if(!bits[j] && bits[j+1]){
            geno[idx] = 1;
        }else if(bits[j] && !bits[j+1]){
            geno[idx] = 3;
        }
        idx++;
    }
}


vector<string> read_snpnames(string filename, int P){
    vector<string> snpnames;
    std::ifstream ifs(filename.c_str());
    std::string line;
    int chromsome;
    string snpname;
    vector <string> fields;
    while(std::getline(ifs, line)) // read one line from ifs
    {
        std::istringstream iss(line); // access line as a stream
        boost::split( fields, line, boost::is_any_of(" \t *"));
        chromsome = (int)atoi(fields[0].c_str());
        snpname = fields[1];
        snpnames.push_back(snpname);
    }
    ifs.close();
    return snpnames;
}


void readPlink(string stringname,int N, int P, unsigned* X){

   // string stringname = dir + dataname;
    FILE *fp;
    unsigned char buff[3];
    string bedfile = stringname +".bed";
    fp = fopen(bedfile.c_str(), "rb");
    if (!fp) return;
    fread(buff, sizeof(char), 3, fp);

    std::bitset<8> magic1(buff[0]);
    std::bitset<8> magic2(buff[1]);
    std::bitset<8> mode0(buff[2]);

    if(magic1.to_ulong() != 108 || magic2.to_ulong() != 27){
     //   cout <<"Error Identifier of plink binary file" << endl;
    }

    unsigned long mode =  mode0.to_ulong();
    if(mode == 0){
        printf ("individual-Major Order:improper type of plink file");
        exit (EXIT_FAILURE);
    }
    //     cout << "SNP-Major Order" << endl;
    // }else if(mode == 0){
    //    cout << "individual-Major Order" << endl;
    // }
// X = new int[N*P];
    int n = 0;
    long charNum = ceil(N*1.0/4)*10000;
    int leftGenoNum = ceil(N*1.0/4)*P;
    int nblock = ceil(N*1.0/4);
    int nSNP = 0;
    while (!feof(fp)) {
        if(leftGenoNum <= 0)
            break;
        if(leftGenoNum <= charNum){
            charNum  = leftGenoNum;
        }
        char* genotype = new char[charNum];
        fread(genotype, sizeof(char), charNum, fp);
        int* geno = new int[4];
        int nSNPc = int(charNum / nblock); //number of SNPs of this iteration
        int idx = 0;
        for (int i=0; i < nSNPc; i++) {
            for(int j=0; j < nblock - 1; j++){
                std::bitset<8> bits(genotype[idx]);
                getFourGentype(geno,bits);
                memcpy(X + nSNP * N + j*4, geno, 4*sizeof(int));
                idx++;
                leftGenoNum -= 1;
            }
            int left = N - (nblock - 1)*4;
            std::bitset<8> bits(genotype[idx]);
            getFourGentype(geno,bits);
            memcpy(X + nSNP * N + (nblock - 1)*4, geno, left*sizeof(int));
            idx++;
            leftGenoNum -= 1;
            nSNP ++;
        }
        delete[] geno;
        delete[] genotype;
        n++;
    //    cout <<n << " processing"<<endl;
    }


}


// [[Rcpp::export]]
void ReadPlinkFamFile3(std::string stringname, CharacterVector FID, CharacterVector IID, 
	NumericMatrix pheno, int nrows, int nPheno){

	// check the number of columns in pheno is the same to nPheno;
	if ( pheno.nrow()!= nPheno){
		
	}

	std::ifstream myfile(stringname.c_str());
	std::string line;

	clock_t t1 = clock();

	int nrow_ind = 0;
	vector <string> tmp;

	if (myfile.is_open()){
		while (nrow_ind < nrows){
			if (nrow_ind % 1000 == 0 && nrow_ind != 0){
				cout << nrow_ind << "-th individual" << ",";
				cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
			}

			getline(myfile, line);
			boost::split(tmp, line, boost::is_any_of(" \t *"));

			FID(nrow_ind) = tmp[0];
			IID(nrow_ind) = tmp[1];
			// sex(nrow_ind) = atoi(tmp[4].c_str());
			for (int j = 0; j < nPheno; j ++){
				pheno(nrow_ind, j) = atof(tmp[5 + j].c_str());
			}
			
			// pheno(nrow_ind) = atof(tmp[4 + whCol].c_str());

			//cout << "value: " << tmp[0] << ";" << tmp[1] << ";" << tmp[2] << ";" << tmp[3] << ";" << tmp[4] 
			//	<< ";" << tmp[5] << ";" << tmp[6] << endl;
			nrow_ind++;
		}
	}
}



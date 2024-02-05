#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame INTERNAL_varTableFromAlns(DataFrame aln){
  CharacterVector s1s = aln["sAlnSeq"];
  CharacterVector s2s = aln["qAlnSeq"];
  IntegerVector s1starts = aln["sStart"];
  IntegerVector s2starts = aln["qStart"];
  IntegerVector inc1s = aln["qInc"];
  IntegerVector inc2s = aln["qInc"];

  std::vector<int> outPos1;
  std::vector<int> outPos2;
  CharacterVector outState1;
  CharacterVector outState2;

  for ( int i=0; i<aln.nrows(); i++ ){

    long int s1start = s1starts[i];
    long int s2start = s2starts[i];

    std::string s1{s1s[i]};
    std::string s2{s2s[i]};

    int inc1 = inc1s[i];
    int inc2 = inc2s[i];

    long int s1pos = s1start;
    long int s2pos = s2start;

    long int s1len = s1.length();
    long int s2len = s2.length();

    // Rcpp::Rcout << s1start << "::";
    // Rcpp::Rcout << s2start << "\n";
    // printf("%s::",s1);
    // printf("%s\n",s2);
    // std::cout << inc1 << "::";
    // std::cout << inc2 << "\n";
    // std::cout << s1pos << "::";
    // std::cout << s2pos << "\n";
    // std::cout << s1len << "::";
    // std::cout << s2len << "\n";
    // std::cout << "----------------------------------------------------\n\n";

    if(s1len!=s2len){
      stop("Sequences must be of equal length");
    }
    int gap1d = 0;
    int gap2i = 0;

    for(int i=0; i<s1len; i++ ){
      int nogap = 0;


      if(s1[i]!='-'){
        //printf("Increment 1\t");
        ++nogap;
        s1pos = s1pos + inc1;
      } else {
        ++gap1d;
      }

      if(s2[i]!='-'){
        //printf("Increment 2\t");
        ++nogap;
        s2pos = s2pos + inc2;
      }  else {
        ++gap2i;
      }

      // fprintf(stderr,"i is %d\n",i);
      // fprintf(stderr,"chars are %c:%c\n",s1[i],s2[i]);
      // fprintf(stderr,"positions are %d:%d\n",s1pos,s2pos);
      // fprintf(stderr,"gaps i/d are are %d:%d\n",gap1d,gap2i);
      // fprintf(stderr,"nogap is %d\n",nogap);

      if(nogap==2){ // no gaps at pos
        // SNP, not next to gap
        if ((s1[i]!=s2[i]) && ((i==0 || i==(s1len-1) ) || (s1[i-1]!='-' && s2[i-1]!='-' && s1[i+1]!='-' && s2[i+1]!='-' ))){
          // printf("%d\t%d\t%c\t%c\n",s1pos,s2pos,s1[i],s2[i]);
          outPos1.push_back(s1pos);
          outPos2.push_back(s2pos);
          outState1.push_back(std::string{s1[i]});
          outState2.push_back(std::string{s2[i]});
        }
        // In ... (s1 is ref; always reported on ref)
        if(gap1d>0){
          // printf("%d\t%s\t%s:%d\t%s\n",s1pos-(1*inc1),"NA","d",gap1d,"NA");
          outPos1.push_back(s1pos-(1*inc1));
          outPos2.push_back(NA_INTEGER);
          outState1.push_back("d:"+std::to_string(gap1d));
          outState2.push_back(NA_STRING);

        } else if (gap2i>0){
          // printf("%d\t%s\t%s:%d\t%s\n",s1pos-((gap2i+1)*inc1),"NA","i",gap2i,"NA");
          outPos1.push_back(s1pos-((gap2i+1)*inc1));
          outPos2.push_back(NA_INTEGER);
          outState1.push_back("i:"+std::to_string(gap2i));
          outState2.push_back(NA_STRING);
        }
        gap1d = 0;
        gap2i = 0;
      }

      //fprintf(stderr,"\n");
    }
  }

  DataFrame outDF = DataFrame::create(
    Named("sPos")     = outPos1,
    Named("qPos")    = outPos2,
    Named("sState")     = outState1,
    Named("qState")     = outState2
  );
  return outDF;
}

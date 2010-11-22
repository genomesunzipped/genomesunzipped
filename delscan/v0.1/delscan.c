// Delscan 0.1: Program to detect deletions and output various useful information regrading
// them. Input is a set of genotype data and a pedigree file.
// Only works for trios currently.
// Author: Don Conrad
// Date: October 27, 2010

// INPUT ARGS: 
// ARGV[1] Name of file containing genotypes

// ARGV[2] Name of pedigree file

// ARGV[3] Output format
// Output format types 
// 0- output summary stats about each family
// 1- output trio genotype configs for putative deletion "stretches"
// 2 -output location of stretches involving at least two interesting
// genotypes-from set of gte or mie 
// 3- 
// 4-
// 5-output type II MIs

// ARGV[4] Genotype file format 
// Genotype format types
// 0 HapMap format
// 1 PedCNV
// 2 BEAGLE format

// ARGV[5] Chromosome number


#define MAXNAME 30
#define MAXLINE 2000

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <cstring>
//#include "delscan.h"
//#include "delpow.h"

//vector<fxn*> dbsnp;

using namespace std;
typedef struct {string famid; string mom; string dad; string kid; int sex; int pop;} Trio;
struct hapmap_locus {string rs; string SNPalleles;  string chrom; string snppos; string strand; string genome_build; string center; string protLSID; string assayLSID; string panelLSID; string QC_code;};
typedef struct hapmap_locus HapMap_Locus;
typedef struct {string mom; string dad; string kid;} family;
struct seqinfo {vector< vector<int> >seqs; vector< vector<int> > seq2; vector< vector<string> > gt; vector<string> name;};
typedef struct seqinfo Seqinfo;
struct newseqs {int **mat; int **pat; int **mat2; int **pat2;};
typedef struct newseqs Newseq;
 

int read_fasta(Seqinfo *,FILE *, int, int, char **);
int **imatrix(int,int,int,int);
char **cmatrix(int,int,int,int);
double **dmatrix(int,int,int,int);
void nrerror(char *);

void allele_count(Seqinfo * s,int,int,int **,int,int);
void get_freq(int,int **,double **,int,int,int);

int bsim(int);

void read_hapmap(ifstream &in, Seqinfo &s, vector<string>& seqnames, vector<HapMap_Locus>& loci,string *header);
void read_beagle(ifstream &in, Seqinfo &s, vector<string> & seqnames, vector<HapMap_Locus>& locs,string* header);
void readped(ifstream &pedin,map <string, Trio>& tm);
void read_family(ifstream &famin, vector<family>& tt, vector<string>& pednames);
void read_pedcnv(ifstream &in,Seqinfo &s,vector<string>&seqnames,vector<string> pednames, map<string,Trio> ped);


int main (int argc, char *argv[])

{
  FILE *fptr, *ofp;
  int **seqs, **seq2, **snpmk,**nall,**inf_mat,**inf3_mat;
  int nseq=0, lseq=0,gte=0,hd,hz,it=0,sk=0,i,j,z=1,miss=0,temp=0, test,test2=0,test3=0,test4,test5,temp1=0,temp2=0;
  int mid,did,dst=0,cid,k,dc=0,mp=0,u=0,u2=0,v=0,w=0,y=0,x=0,p1=0,p2=0,p12=0,p22=0,q=0,t=0,vt=0,tv=0,start_mi=0,latest_mi=0,seen=0;
  int *match_d,*match_c,*match_m,*pa_rs,pa1,pa2,end_flag,flag,match=0,seen_deletion_flag,dist;
  int middel=0,del=0, d=0, length=0,del_dst[30000],dc_violation=0,del_family[30000],del_dc[30000];  
  long int score[11], score2[11], read[11];

  double **delhz,**freq_vec;
  vector<HapMap_Locus> locs;
  char line[MAXLINE], fname[MAXNAME], pid[MAXNAME];
  map<string, Trio> tm;
  string header;
  Newseq news;
  Seqinfo si;

  vector<string> pednames;
  vector<string> seqnames;
  vector<family> tt;



  for (i=0;i<12;i++){score[i]=0;score2[i]=0;read[i]=0;}

ifstream pedin(argv[2]);  
if(!pedin) {cerr<<" cannot open pedigree file.\n"; return 1;}
readped(pedin,tm);

ifstream famin(argv[2]);
if(!famin) {cerr<<" cannot open family file.\n"; return 1;}
read_family(famin,tt,pednames);



 ifstream gtfile(argv[1]); 

 if(!gtfile){cerr<<" Can't open gt file\n";}
 if (atoi(argv[4])==0){ 
     read_hapmap(gtfile, si, seqnames, locs,&header);}
 else if (atoi(argv[4])==1){ 
    read_pedcnv(gtfile, si, seqnames,pednames,tm);}
 else if (atoi(argv[4])==2){
   read_beagle(gtfile, si,seqnames,locs,&header); }
 else {cerr<<" Invalid file format specified "<<endl; return 1;}
 

   gtfile.close();

 nseq=seqnames.size();
 lseq=si.seqs[0].size();



vector<int> start;
vector<int> end;
vector<int> twohit(lseq+1);
vector<float> mom_hzarray(lseq+1);
vector<float> dad_hzarray(lseq+1);

 vector<int> inf(lseq,0);
 vector<int> inf2(lseq,0);
match_c  = (int *) calloc(91,sizeof(int));
match_m  = (int *) calloc(91,sizeof(int));
match_d  = (int *) calloc(91,sizeof(int));
pa_rs   = (int *) calloc(4990,sizeof(int));

snpmk     = (int **) imatrix(0,30000,0,500);
delhz     = (double **) dmatrix(0,30000,0,500);
nall      = (int **) imatrix(0,lseq,0,5);
freq_vec  = (double **) dmatrix(0,30000,0,500);
inf_mat  = (int **) imatrix(0,3,0,lseq);
inf3_mat = (int **) imatrix(0,11,0,9);



 start.push_back(0);
 end.push_back(0);


 for (i=0;i<tt.size();i++){
   inf.clear();
   inf2.clear();
       mom_hzarray.clear();
   dad_hzarray.clear();

     //did,cid,mid defines individuals position in the input data set (vs. position in pedinfo file)
     
     j=0; while(seqnames[j].compare(tt[i].dad)!=0) {j++;} did=j;
     j=0; while(seqnames[j].compare(tt[i].mom)!=0) {j++;} mid=j;
     j=0; while(seqnames[j].compare(tt[i].kid)!=0) {j++;} cid=j;





   for (k=0;k<lseq;k++){
     end_flag=0; match=0;

     //skip NRX  
     //if ((atoi(argv[4])==23 ) && (locs[k].snppos > 2290000) && (locs[k].snppos < 153430143)){continue;}

    
     //Skip Ig and TCR genes
     //          if ((atoi(argv[4])==2 ) && (locs[k].snppos > 89000000) && (locs[k].snppos < 89500000)){continue;}
     //if ((atoi(argv[4])==7 ) && (locs[k].snppos > 37920000) && (locs[k].snppos < 38100000)){continue;}
     //if ((atoi(argv[4])==7 ) && (locs[k].snppos > 141600000) && (locs[k].snppos < 141960000)){continue;}
     //if ((atoi(argv[4])==14 ) && (locs[k].snppos > 20550000) && (locs[k].snppos < 2101000)){continue;} 
     //if ((atoi(argv[4])==14 ) && (locs[k].snppos > 104000000) && (locs[k].snppos < 105000000)){continue;}
     //if ((atoi(argv[4])==22 ) && (locs[k].snppos > 20700000) && (locs[k].snppos < 22000000)){continue;}  
     

     //skip duplicated loci
     //      if(k>1){if(locs[k].snppos==locs[k-1].snppos && k < (lseq-1)) {cout<<"dup "<<locs[k].snppos<<endl;u++;continue;}} 
     //    if(k>1){if(locs[k].rs==locs[k-1].rs && k <(lseq-1)) {cout<<"DUP RS"<<endl;u2++;continue;}} 




     //inf[k]=1 informative
     //inf[k]=2 uninformative


     inf[k]=0;
     inf2[k]=0;
       
	 miss=0;        //miss=1 indicates a parent is missing data; miss=2, child 


    //missing loci 
    // mp = 1 if parent who is missing/homozygous is dad; 2=mom;3=both
    // hz = 1 means non-missing parent is homozygous
     
     // 1. Child is missing data or homozygous
	 if (si.seqs[cid][k]==si.seq2[cid][k]){
        
     // 1b.mom and dad homozygous or missing data; class D
     if (si.seqs[did][k] == si.seq2[did][k] && si.seqs[mid][k] == si.seq2[mid][k]){inf[k]=4;}
     
     // 1c. mom homozygous or missing data, dad heterozygous; class E
     else if (si.seqs[did][k]!=si.seq2[did][k] && si.seqs[mid][k]==si.seq2[mid][k]){inf[k]=5;}
 
      // 1d. dad homozygous or missing data, mom heterozygous; class F
     else if (si.seqs[mid][k]!=si.seq2[mid][k] && si.seqs[did][k]==si.seq2[did][k]){inf[k]=6;}
  
     //  both parents het   
         else if (si.seqs[mid][k]!=si.seq2[mid][k] && si.seqs[did][k]!=si.seq2[did][k]){inf[k]=7;}
    }
     
	 else {inf[k]=7;}
   
     

     
     if (inf[k]==0){cerr<<"problems with assigning trio config in main "<<si.seqs[cid][k]<<" "<<si.seq2[cid][k]<<" "<<si.seqs[mid][k]<<" "<<si.seq2[mid][k]<<" "<<si.seqs[did][k]<<" "<<si.seq2[did][k]<<endl; 
                               exit(1);}


     //adding probabilities for transmission for MCMC output
     if (si.seqs[did][k]==si.seq2[did][k]){
     if (si.seqs[did][k]==5){dad_hzarray[k]=.75;}
     else {dad_hzarray[k]=1;}}
     else {dad_hzarray[k]=.5;}
      

     if (si.seqs[mid][k]==si.seq2[mid][k]){
     if (si.seqs[mid][k]==5){mom_hzarray[k]=.75;}
     else {mom_hzarray[k]=1;}}
     else {mom_hzarray[k]=.5;}



     if (si.seqs[did][k]==5 || si.seqs[mid][k]==5){miss=1;}
     if (si.seqs[did][k]==5 && si.seqs[mid][k]==5){miss=2;}
     if (si.seqs[cid][k]==5){miss=3;}
      
         if (miss==1){       
         //reduced MIE search for case of missing data   
        
         if (si.seqs[did][k]==5 && si.seqs[cid][k]==si.seq2[cid][k]){ 
     
         if (si.seqs[cid][k]!=si.seqs[mid][k] && si.seqs[cid][k]!=si.seq2[mid][k]){
	 z++;d=1;dc=2;inf[k]=1;}}
     
         else if (si.seqs[mid][k]==5 && si.seqs[cid][k]==si.seq2[cid][k]){
     
         if (si.seqs[cid][k]!=si.seqs[did][k]&&si.seqs[cid][k]!=si.seq2[did][k]){
	 z++;d=1;dc=1;inf[k]=2;}}}

         else if (miss==0){
        
	 p1=0;p2=0;p12=0;p22=0;

         if ((si.seqs[cid][k] == si.seq2[mid][k]) || (si.seqs[cid][k] == si.seqs[mid][k])) {p1++;}
         if ((si.seqs[cid][k] == si.seq2[did][k]) || (si.seqs[cid][k] == si.seqs[did][k])) {p2++;}
         if ((p1==0) && (p2==0)) {test2=1;} else {test2=0;}    
     
         if (si.seq2[cid][k] == si.seqs[mid][k]){p12++;}
         if (si.seq2[cid][k] == si.seq2[mid][k]){p12++;}
         if (si.seq2[cid][k] == si.seqs[did][k]){p22++;}
         if (si.seq2[cid][k] == si.seq2[did][k]){p22++;}
         if((p12==0) && (p22==0)){test3=1;} else {test3=0;}
         p1=p1+p12;
         p2=p2+p22;      
         
         if((test2==0)||(test3==0)){
	 if (p1==0){z++;d=1;dc=2;inf[k]=1;} //dc=2 indicates mom is deletion carrier
         else if (p2==0){z++;d=1;dc=1;inf[k]=2;}
	 
         }
	 if (test2==1 || test3==1) {gte=1;w++;inf[k]=3;}
         }
   

         inf_mat[i][k]=inf[k];
            
         if (inf[k]==0){cout<<"Error in assigning inf"<<endl<<miss<<" "<<inf[k]<<" "<<locs[k].snppos<<" "<<si.seqs[cid][k]<<" "<<si.seq2[cid][k]<<" "<<si.seqs[mid][k]<<" "<<si.seq2[mid][k]<<" "<<si.seqs[did][k]<<" "<<si.seq2[did][k]<<endl; }


     if (d==0) {
     if (si.seqs[mid][k]==si.seq2[mid][k]){hz=1;}
     else {hz=0;}     

     //print location of GTEs if switch is thrown
     if (gte==1 && atoi(argv[3])==5){ cout<<atoi(argv[4])<<" "<<locs[k].snppos<<" "<<locs[k].rs<<" "<<i<<" "<<tt[i].kid<<" "<<hz<<endl;}
     else if (gte==1 && atoi(argv[3])==3){
       cout<<tt[i].kid<<" "<<tt[i].dad<<" "<<tt[i].mom<<" "<<inf[k]<<" "<<locs[k].snppos<<" "<<si.seqs[cid][k]<<" "<<si.seq2[cid][k]<<" "<<si.seqs[did][k]<<" "<< si.seq2[did][k]<<" "<<si.seqs[mid][k]<<" "<< si.seq2[mid][k]<<endl;}
}  


     //If locus is MIE that fits our model, print data and keep track of deletion tract
     if (d==1){ 
     switch(atoi(argv[3])){
     case 1:cout<<tt[i].kid<<" "<<tt[i].dad<<" "<<tt[i].mom<<" del "<<inf[k]<<" "<<locs[k].snppos<<" "<<si.seqs[cid][k]<<" "<<si.seq2[cid][k]<<" "<<si.seqs[did][k]<<" "<< si.seq2[did][k]<<" "<<si.seqs[mid][k]<<" "<< si.seq2[mid][k]<<endl; break;
     case 4: cout <<endl;break;      
     case 5: cout<<atoi(argv[4])<<" "<<locs[k].snppos<<" "<<locs[k].rs<<" "<<i<<" "<<tt[i].kid<<" "<<hz<<" 1"<<endl; break;
     }
    
    
    //end flags
    //1-normal termination
    //2-disrupted by wrong MI class?
    //3-disrupted by MI not compatible with transmission?
  
     if (inf[k-1]==1 && inf[k]==2){

     if (length>1){end_flag=2;}  
     if (length==1 &&  atoi(argv[3])==2){ 
     cout<<argv[4]<<" "<<locs[k-1].snppos<<" "<<locs[k].snppos<<" "<<"2 FP 2 FP "<<tt[i].kid<<" 2"<<endl;}
     middel=0;length=0;}

     if (inf[k-1]==2 && inf[k]==1){
     if (length>1){end_flag=2;}
     if (length==1 &&  atoi(argv[3])==2){ 
       cout<<argv[4]<<" "<<locs[k-1].snppos<<" "<<locs[k].snppos<<" "<<"2 FP 2 FP "<<tt[i].kid<<" 2"<<endl;}
 middel=0;length=0;}  
 
     switch(middel){
     case 1: if (length==1){
          dst=1;
          del++; 
          del_family[del]=i;
          del_dc[del]=dc;

     if (atoi(argv[3])==1){	 printf("Start reverse search\n|\n|\n|\n|\n+\n");}
 
       temp=k; 

     while (inf[temp]!=7 && inf[temp]!=3 && dc_violation==0 && temp >1 ){

       //CENT     if (atoi(locs[k].snppos.c_str()) > cenout[atoi(argv[4])] && atoi(locs[temp-1].snppos.c_str()) < cenin[atoi(argv[4])]){dc_violation=2;}
       
      temp--;
      switch(dc){
      case 1:if (si.seqs[did][temp]!=si.seq2[did][temp]){dc_violation=1;
	  if (atoi(argv[3])==1){printf("dc viol %d %d\n",dc_violation,locs[temp].snppos.c_str());}} break;
      case 2:if (si.seqs[mid][temp]!=si.seq2[mid][temp]){dc_violation=1;
	  if (atoi(argv[3])==1){printf("dc viol %d %d\n",dc_violation,locs[temp].snppos.c_str());}} break;
      default:break;} 
      

      if(atoi(argv[3])==1){ cout<<" m-u "<<locs[temp].snppos<<" "<<locs[temp].rs<<" "<<inf[temp]<<" "<<si.seqs[cid][temp]<<" "<<si.seq2[cid][temp]<<" "<<si.seqs[did][temp]<<" "<< si.seq2[did][temp]<<" "<<si.seqs[mid][temp]<<" "<< si.seq2[mid][temp]<<endl;}

      }

      if (dc_violation==2){temp++;}
 
      //CENT      if (dc_violation==1 && atoi(locs[temp+1].snppos.c_str()) > cenout[atoi(argv[4])] && atoi(locs[temp].snppos.c_str()) < cenin[atoi(argv[4])]){temp++;}

      if (atoi(argv[3])==1){        cout<<"End of reverse search "<<dc_violation <<" push back "<<locs[temp].snppos<<endl;}
         start.push_back(temp); 
         temp=0;

       }
      if (length>1){

      if (del_dc[del]==1 && inf[k]==1){end_flag=2;}
      if (del_dc[del]==2 && inf[k]==2){end_flag=2;}
   }

      dst++; length++; break;

      case 0: middel=1;length++; break;
      default: printf("error in deletion length estimation\n"); break; }}
  
  
      //If locus is not MIE, but we are in the middle of a deletion stretch
      else if ((middel) && length >= 1){

         //data is uninformative
	 if (inf[k]!=7 && inf[k]!=3){
       
            if (inf[k-1]==1){
	       if (inf[k]==4 || inf[k]==5){match=1;}}

             if (inf[k-1]==2){
	       if (inf[k]==4 || inf[k]==6){match=1;}}
		     
              
	   if (length==1 && match==1){
                        
	     length++; 	del++;    
             dst=1;del_family[del]=i; del_dc[del]=dc;
 
	   

	 if (atoi(argv[3])==1){	 printf("Start reverse search\n|\n|\n|\n|\n+\n");}

         temp=k-1; 

         while (inf[temp]!=7 && inf[temp]!=3 && dc_violation==0 && temp >1){

  temp--;
  //CENT   if (atoi(locs[k].snppos.c_str()) > cenout[atoi(argv[4])] && atoi(locs[temp-1].snppos.c_str()) < cenin[atoi(argv[4])]){dc_violation=2;}

      switch(dc){
      case 1:if (si.seqs[did][temp]!=si.seq2[did][temp]){dc_violation=1;
	  if (atoi(argv[3])==1){printf("dc viol %d %d\n",dc_violation,locs[temp].snppos.c_str());}} break;
      case 2:if (si.seqs[mid][temp]!=si.seq2[mid][temp]){dc_violation=1;
	  if (atoi(argv[3])==1){printf("dc viol %d %d\n",dc_violation,locs[temp].snppos.c_str());}} break;
      default:break;} 
    
         
      if(atoi(argv[3])==1){ 	 cout<<"m-u"<< locs[temp].snppos<<" "<<locs[temp].rs<<" "<<inf[temp]<<" "<<si.seqs[cid][temp]<<" "<<si.seq2[cid][temp]<<" "<<si.seqs[did][temp]<<" "<< si.seq2[did][temp]<<" "<<si.seqs[mid][temp]<<" "<< si.seq2[mid][temp]<<endl;}

      }

 if (dc_violation==2){temp++;}

 //CENT  if (dc_violation==1 && atoi(locs[temp+1].snppos.c_str()) > cenout[atoi(argv[4])] && atoi(locs[temp].snppos.c_str()) < cenin[atoi(argv[4])]){temp++;}

 if (atoi(argv[3])==1){         printf("End Reverse Search %d push back %d\n|\n|\n",dc_violation,locs[temp].snppos.c_str());}
         start.push_back(temp); 
         temp=0;
	   
	   }
	

         else if (length==1 && match==0){
	   if (atoi(argv[3])==1){printf("\n");}
         middel=0;length=0;dst=0;dc=0;}

	else if (length > 1){
     
     if (del_dc[del]==1 && inf[k]==5){end_flag=1;}
     if (del_dc[del]==2 && inf[k]==6){end_flag=1;}
     
         length++;

       switch(atoi(argv[3])){
       case 1:  cout<<"m-u"<< locs[temp].snppos<<" "<<locs[temp].rs<<" "<<inf[temp]<<" "<<si.seqs[cid][temp]<<" "<<si.seq2[cid][temp]<<" "<<si.seqs[did][temp]<<" "<< si.seq2[did][temp]<<" "<<si.seqs[mid][temp]<<" "<< si.seq2[mid][temp]<<endl; break;
    
       }}

	 
	 
	 }

       //data is informative
       else if (inf[k]==7){
	 
    	   
	 if (length > 1){end_flag=1;}

       
         else if (length==1){ if (atoi(argv[3])==1){printf("\n");}}

	 middel=0;length=0;mp=0;

       }
       
      // data is a genotype error not compatible with deletion transmission
	 else if (inf[k]==3){ 
      
	   if (length==1){if (atoi(argv[3])==1){printf("\n");}}

          

	  else if (length >1){end_flag=3;} middel=0;length=0;	 
	  }

       
	 // single MIE 
	 else if (length==1){ if (atoi(argv[3])==1){printf("\n");} dst=0;middel=0;length=0;mp=0;}}
   
 
     //special case where deletion affects last SNP on chromosome
     if (k == (lseq-1) && middel) 
       {
     if (length>1){end_flag=1;}
    
    length=0;middel=0;}

     //special case where deletion overlaps centromere
     //     if (atoi(locs[k].snppos.c_str()) < cenin[atoi(argv[5])] && atoi(locs[k+1].snppos.c_str()) > cenout[atoi(argv[5])] && length > 1){
     //end_flag=1;length=0;middel=0;}

     //special case where deletion is in PAR1 of X chromosome
     //     if ((atoi(locs[k+1].snppos.c_str()) > 2290000) && middel && (atoi(argv[4])==23))
     // {
     //if (length>1){
     //end_flag=1;}
     //length=0;middel=0; 
     //}


     //Ending Deletion Process

     if (end_flag){


          end.push_back(k);

	  del_dst[del]=dst; //was equal to "it" 
          
          if (end_flag==3 && dst >1) {del_dst[del]=-1;}    
	  else if (end_flag >1 ) { del_dst[del]=1;}
          
         for (sk=start[del],it=0;sk<=end[del];sk++){

           snpmk[del][it]=inf[sk];

           if(dc==2){     delhz[del][it]=dad_hzarray[sk];

            get_freq(si.seqs[mid][sk],nall,freq_vec,sk,it,del);}

           if(dc==1){     delhz[del][it]=mom_hzarray[sk];

           get_freq(si.seqs[did][sk],nall,freq_vec,sk,it,del);}

      it++;
}   
 
          length=end[del]-start[del]+1;

       switch (atoi(argv[3])){
       case 1:
       cout<<"end "<<locs[k].snppos<<" "<<locs[k].rs<<" "<<inf[k]<<" "<<si.seqs[cid][k]<<" "<<si.seq2[cid][k]<<" "<< si.seqs[did][k]<<" "<< si.seq2[did][k]<<" "<<si.seqs[mid][k]<<" "<<si.seq2[mid][k]<<endl; break;
       case 2: if (dst>1){cout<<atoi(argv[5])<<"\t"<<locs[start[del]].snppos<<"\t"<<locs[end[del]].snppos<<"\t"<<length<<"\t"<<atoi(locs[end[del]].snppos.c_str())-atoi(locs[start[del]].snppos.c_str())<<"\t"<<dst<<"\t"<<dc<<"\t"<<tt[i].kid<<endl;}break;
       }

     
     length=0;dst=0;middel=0;mp=0;d=0;
   
     }

     if (miss>0){vt++;}
     d=0;mp=0; p1=0;p2=0;t=0;p12=0;p22=0;test2=0;test3=0;miss=0;hz=0;gte=0; dc_violation=0;
    
 
   }
   


   dc=0;length=0;middel=0;

if(atoi(argv[3])==0){
  printf ("CID %5d hom sites %5d repeats %5d ? %5d MI %5d %5d num MI sites %4d gte %3d lseq %6d\n",cid,y-x,u,u2,z,v-vt,vt,w,lseq); vt=0;u=0;v=0;w=0;y=0;x=0;z=0;u2=0;}   /*u=repeated genotype,v = missing data, w=nondeletion gte,y is total number of homozygous sites, x is number that are 5/5, z is mendelian errors, vt is number of loci per family where at least one member has deletion*/
 }

 //     cout<<" last k "<<k<<endl;


}

/*----------------------------------------------------------------------------------------------------------------------*/

void read_pedcnv(ifstream &in,Seqinfo &s,vector<string>&seqnames,vector<string> pednames, map<string,Trio> ped){


  string line,sblah;
  int iblah,pos=0,i=0,mpos=0;
 vector<int> tivec;
 vector<string> tvec;
 vector<string>::iterator it;
 vector<string> master_seqnames;   
  std::stringstream out;

 getline(in,line);
istringstream iss(line);       

 iss>>sblah; //first entry should be "CNVE"

while(iss>>sblah){

     master_seqnames.push_back(sblah);
     if (ped[sblah].famid.compare("999")==0){continue;}
     seqnames.push_back(sblah);
 }

//initialize space required for gts
 for (i=0;i<seqnames.size();i++){
    s.seqs.push_back(tivec);
   s.seq2.push_back(tivec);
   s.gt.push_back(tvec);
 }

 while(   getline(in,line)){
   istringstream iss(line);       
   pos=0;mpos=0;
   iss>>sblah;
   s.name.push_back(sblah);
   while(iss>>iblah){

     if (ped[master_seqnames[mpos]].famid.compare("999")==0){mpos++;continue;} //skip singletons
     switch(iblah){        
     case 0: { s.seqs[pos].push_back(0); s.seq2[pos].push_back(0); break;}
     case 1: { s.seqs[pos].push_back(0); s.seq2[pos].push_back(1); break;}
     case 2: { s.seqs[pos].push_back(1); s.seq2[pos].push_back(1); break;}
     case 3: { s.seqs[pos].push_back(1); s.seq2[pos].push_back(2); break;}
     case 4: { s.seqs[pos].push_back(2); s.seq2[pos].push_back(2); break;}
     case -9: { s.seqs[pos].push_back(5); s.seq2[pos].push_back(5); break;}
     case -99: { s.seqs[pos].push_back(5); s.seq2[pos].push_back(5); break;}
     default: { s.seqs[pos].push_back(5); s.seq2[pos].push_back(5); break;}}

      out << iblah;
          sblah = out.str();
	  out.str(""); //empty the string
         s.gt[pos].push_back(sblah);   
     
	 pos++;  mpos++;     
   }


 }


  cerr<<"Read "<<seqnames.size()<<" samples from cnv file"<<endl;
  cerr<<"Read "<<s.seqs[0].size()<<" loci from cnv file"<<endl;
//cerr<<"last intensity datapoint "<<ints[ints.size()-1][ints[0].size()-1]<<endl;  
 cerr<<"First Sample Name "<<seqnames[0]<<" Last Sample Name "<<seqnames[seqnames.size()-1]<<endl;

}



void read_hapmap(ifstream &in, Seqinfo &s, vector<string> & seqnames, vector<HapMap_Locus>& loci,string* header){

  string line;
  string sblah;
  char *c;
  int i=0,j=0,k=0;
 HapMap_Locus loc;
 map<char,int> decode;
 decode['A']=0;
 decode['C']=1;
 decode['G']=2;
 decode['T']=3;
 decode['N']=5;

   //header
   getline(in,line);
   istringstream iss(line);       
   *header=line;
   for (j=0;j<11;j++){iss>>sblah;}
       
       while (iss>>sblah){
	 seqnames.push_back(sblah);
         
} 

 vector<int> tivec;
 for (i=0;i<seqnames.size();i++){s.seqs.push_back(tivec);}
 for (i=0;i<seqnames.size();i++){s.seq2.push_back(tivec);}

 vector<string> tvec;
 for (i=0;i<seqnames.size();i++){s.gt.push_back(tvec);}

 while(getline(in,line)){
    istringstream iss(line);       

    iss>>loc.rs;
    iss>>loc.SNPalleles;
    iss>>loc.chrom;
    iss>>loc.snppos;
    iss>>loc.strand;    
    iss>>loc.genome_build;
    iss>>loc.center;
    iss>>loc.protLSID;
    iss>>loc.assayLSID;
    iss>>loc.panelLSID;
    iss>>loc.QC_code;
    loci.push_back(loc);

    for (k=0;k<seqnames.size();k++){

    iss>>sblah;
    s.gt[k].push_back(sblah);

    strncpy(c,sblah.substr(0,1).c_str(),1);
    s.seqs[k].push_back(decode[*c]);
    //    cout<<s.gt[k][s.gt[k].size()-1]<<" c1 "<<*c<<" ";
    strncpy(c,sblah.substr(1,1).c_str(),1);
    //cout<<"c2 "<<*c<<endl;
    s.seq2[k].push_back(decode[*c]); 

    }

 }

  cerr<<"Read "<<seqnames.size()<<" samples from snp file"<<endl;
  cerr<<"Read "<<loci.size()<<" loci from snp file"<<endl;
//cerr<<"last intensity datapoint "<<ints[ints.size()-1][ints[0].size()-1]<<endl;  
 cerr<<"First seqname "<<seqnames[0]<<" Last Findiv "<<seqnames[seqnames.size()-1]<<endl;
cerr<<"First Locus "<<loci[0].rs<<" Last Locus "<<loci[loci.size()-1].rs<<endl;  




}

void read_beagle(ifstream &in, Seqinfo &s, vector<string> & seqnames, vector<HapMap_Locus>& loci,string* header){
 string line;
  string sblah;
  char c[1];
  int i=0,j=0,k=0,skip=1;
 HapMap_Locus loc;
 map<char,int> decode;
 decode['A']=0;
 decode['C']=1;
 decode['G']=2;
 decode['T']=3;
 decode['N']=5;

   //header
   getline(in,line);
   istringstream iss(line);       
   *header=line;
   for (j=0;j<2;j++){iss>>sblah;}
       while (iss>>sblah){
	 if (skip==1){
	   seqnames.push_back(sblah);}
	 skip*=-1;
       } 

 vector<int> tivec;
 for (i=0;i<seqnames.size();i++){s.seqs.push_back(tivec);}
 for (i=0;i<seqnames.size();i++){s.seq2.push_back(tivec);}

 vector<string> tvec;
 for (i=0;i<seqnames.size();i++){s.gt.push_back(tvec);}

 while(getline(in,line)){
    istringstream iss(line);       

    iss>>loc.assayLSID;
    iss>>loc.snppos;

    loci.push_back(loc);

    for (k=0;k<seqnames.size();k++){

    iss>>c;
    s.seqs[k].push_back(decode[*c]);
    iss>>c;
    s.seq2[k].push_back(decode[*c]); 
    //    s.gt[k].push_back(sblah);
    }

 }

 cerr<<"Read "<<seqnames.size()<<" samples from snp file"<<endl;
 cerr<<"Read "<<loci.size()<<" loci from snp file"<<endl;
//cerr<<"last intensity datapoint "<<ints[ints.size()-1][ints[0].size()-1]<<endl;  
 cerr<<"First seqname "<<seqnames[0]<<" Last Findiv "<<seqnames[seqnames.size()-1]<<endl;
 cerr<<"First Locus "<<loci[0].snppos<<" Last Locus "<<loci[loci.size()-1].snppos<<endl;  




}

int read_fasta(Seqinfo * s, FILE *ifp,int lseq,int nseq,char **seqnames) 

{

	int h=0,i,z=0,site=1, seq=1, cts[5];
	char line[MAXLINE], *c, bases[6]="-TCAG";
	/*was bases[5]="TCAG-"*/
	/*	printf("\n\nReading sequences in fasta format\n\n");*/
	for (i=0;i<5;i++) cts[i]=0;


	while (!feof(ifp) && (seq<=nseq)) {
		fgets(line, MAXLINE, ifp);
		if ((c = (char *) strchr(line, '>')) != NULL) {
		    if  (z > 1 && z % 2 == 0) {seq++;h=0;}

		    /*	    printf("Sequence :%3i ", seq);*/
			strncpy(seqnames[seq], (c+1), MAXNAME);
			for (i=1;i<=MAXNAME;i++) if(seqnames[seq][i]=='\n') seqnames[seq][i]='\0';
			/*			printf("%s: ", seqnames[seq]);*/

			site=1;
			for (i=0;i<5;i++) cts[i]=0;
			while (site<lseq) {
				fgets(line, MAXLINE, ifp);
				for (c=line; (*c)!='\0'; c++) {
				  /*printf("%s %d %d %d\n",c,site,z,h);*/
				   switch(*c) {
					case 'T': case 't': case '0': {
					  if (h==0) {s->seqs[seq][site] = 3;}
                                          else {s->seq2[seq][site]=3;}
						site++;
						cts[0]++;
						break;
					}
					case 'C': case 'c': case '1': {
					  if (h==0){s->seqs[seq][site] = 1;}
                                          else {s->seq2[seq][site]=1;}        
						site++;
						cts[1]++;
						break;
					}
					case 'A': case 'a': case '2': {
                                          if (h==0) {
					    s->seqs[seq][site] = 0;}
                                          else {s->seq2[seq][site]=0;}
						site++;
						cts[2]++;
						break;
					}
					case 'G': case 'g': case '3' :{
					  if (h==0){s->seqs[seq][site] = 2;}
					  else {s->seq2[seq][site] = 2;}
					        site++;
						cts[3]++;
						break;
					}
					case'-': case'N': case'n': case'?': case'R': case'Y': case'M': case'K': case'S': case'W': case'H': case'B': case'V': case'D' :{
					  if (h==0) {s->seqs[seq][site]=5;}
				          else {s->seq2[seq][site]=5;}		site++;
						cts[4]++;
						break;
					}
					case '>': {printf("\n Error in sequence file(%i of %i bases read)\n",site,lseq); exit(1);}
					default: {
						break;
					}
				   }
				}
			}
			if (site != lseq+1) {printf("\nSequences incorrect length (%i %i)\n\n",site,lseq); exit(1);}
			/*for (i=0;i<5;i++) printf("%c:%5i ",bases[i], cts[i]);
			  printf("\n");*/  h++;z++;
		}
	}
	if (seq!=nseq) {printf("\n\nDid not read %i sequences \n\n",nseq); exit(1);}
	return 1;
}            
  


/*----------------------------------------------------------------------------------------------------------------------------------------------*/
/* functions from the tools.c group courtesy of Paul Fearnhead*/
int **imatrix(int nrl,int nrh,int ncl,int nch)
     /*int nrl,nrh,ncl,nch;*/
{
	int i,**m;

	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

char **cmatrix(int nrl,int nrh,int ncl,int nch)

{
        int i;
	char **m;

        m=(char **)malloc((unsigned) (nrh-nrl+1)*sizeof(char*));
        if (!m) nrerror("allocation failure 1 in cmatrix()");
        m -= nrl;

        for(i=nrl;i<=nrh;i++) {
                m[i]=(char *)malloc((unsigned) (nch-ncl+1)*sizeof(char));
                if (!m[i]) nrerror("allocation failure 2 in cmatrix()");
                m[i] -= ncl;
        }
        return m; 
}

void nrerror(char error_text[])
     /*char error_text[];*/
{
  

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


/*---------------------------------------------------------------------------------------------------------*/

void allele_count(Seqinfo * s,int nseq,int lseq,int **nall,int fl,int hd) 
{
        
  int seq=0, site=0, i=0;
        
        for (site=1; site<=lseq; site++) {
	  for(i=0; i<=5; i++) {nall[site][i]=0;}
	      for (seq=1; seq<=nseq; seq++) {
		nall[site][s->seqs[seq][site]]++;
                nall[site][s->seq2[seq][site]]++;                                   
                        }
	      //  cout<<"a's:  "<<nall[site][0]<<" c's: "<<nall[site][1]<<endl;
   }   
              
	
	 	 
	  }


double **dmatrix(int nrl,int nrh,int ncl,int nch)

{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

void get_freq(int base,int **nall,double **freq_vec,int k,int it,int del){
  //dc=1 means father is carrier, dc=2 mom.
  int denom=0;

  denom=180-nall[k][5];

  //missing data case
  if (base==5){freq_vec[del][it]=1;}
  else{
    freq_vec[del][it]=(double)nall[k][base]/(double)denom;
  }
}


void readped (ifstream &pedin,map <string, Trio>& tm){

string line;
Trio temp_ped;
int temp=0,i=0;


 while(getline(pedin,line)){
    istringstream iss(line);       
   
  iss >>  temp_ped.famid;
  iss >> temp_ped.kid;
  iss>>temp_ped.dad;
  iss>>temp_ped.mom;
 
  iss>>temp_ped.sex;
  iss>>temp_ped.pop;


    tm[temp_ped.kid]=temp_ped;
  // tm[i]=temp_ped;
  i++;
 }
     cerr<<"Recorded "<<tm.size()<<" pedigree entries"<<endl;
 pedin.close();

}



void read_family(ifstream &famin, vector<family>& tt, vector<string>& pednames){

  int i=0,j=0;
 string line,pid;
istringstream iss(line);       
 int temp;
 family tfam;
 vector<string> pvec;

 while(getline(famin,line)){
    istringstream iss(line);       
    iss >> temp;
    if (temp==999){ 
      //for (i=0;i<5;i++){iss>>temp;}
  continue;}
  
  else { 
      pvec.clear();
      for (i=0;i<5;i++){     iss >>pid;
	pvec.push_back(pid);}
      if (pvec[1].compare("0")!=0 && pvec[2].compare("0")!=0){
	tfam.kid=pvec[0];
	tfam.dad=pvec[1];
	tfam.mom=pvec[2];
	tt.push_back(tfam);
      }
    }
 }
 famin.close();
 cerr<<"Read "<<tt.size()<<" Families "<<endl;
 }

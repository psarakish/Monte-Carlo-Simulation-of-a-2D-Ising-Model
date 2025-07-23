#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <fstream>
#include <cstdlib>
using namespace std;

double f (int n, vector <int> s, int i, int j){
    double F=0;
    F = ( s[(i+1)+n*j] + s[(i-1)+n*j] + s[i+(j+1)*n] + s[i+(j-1)*n] );
    return F;
}



double energy(int n, vector <int> s, int J, int B){
    double sumj=0;
    double sumb=0;
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            sumj+= s[i+n*j]* ( s[(i+1)+j*n] + s[(i-1)+j*n] + s[i+(j+1)*n] +s[i+(j-1)*n]);                                
        }
    }
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            sumb+= s[i+j*n];
        }
    }
    return -J*sumj - B*sumb;
}



double magnitikh(int n, vector <int> s){
    double summ=0;
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            summ+= s[i+n*j];
        }
    }
    return abs(summ); //κανονικοποιημενη μαγνητικη ροπη
}



void Spin (int n, vector <int> & s){
    for (int i{0}; i<n; i++){
        for (int j{0}; j<n; j++){
            s[i+n*j]=1; // αρχικοποιηση πινακα σπιν με ολα τα σπιν πανω
        }
    }

    return;
}



int main () {
    int loops {10000}; //πληθος προσομοιωσεων
    int n{10};
    double J{1.}; //σταθερα συζευξης
    double B{0.}; //μαγνητικο πεδιο
    vector <int> s(n*n); //πλεγμα
    Spin(n,s); //συναρτηση για αρχικοποιηση πινακα σπιν
    // περιοδικες συνοριακες
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
             s[(n-1)+j*n]=s[0+j*n];
            s[0+j*n]=s[(n-1)+j*n];
            s[i+(n-1)*j]=s[i+0*j];
            s[i+0*j]=s[i+(n-1)*j];
          }
    }
    ofstream out1{"montecarlo1(10)"};
    ofstream out2{"montecarlo2(10)"};
    double h,r,meshE,meshm,EP;

    for (double T=0.1; T<=4; T+=0.1){
        double b_const = 1/T ; // σταθερα Boltzmann εστω 1
        double C=0;
        double sumE{0.} , summm{0.} , sumE2{0.} , summm2{0.} ; 
        int population{0}; //αριθμος μετρησεων


        for (int ii=1       ; ii<loops; ii++){ //ξεκιναμε τις σαρωσεις

            
            
            
            // metropolis algorithm
            for (int i=0; i<n; i++){
                for (int j=0; j<n; j++) {
                    r = exp(-2.0*b_const*s[i+j*n]*(J*f(n,s,i,j)+B));
                    if (r>1){
                        s[i+j*n]*=-1;
                    }
                    else {
                        h=1.0*rand()/RAND_MAX; // τυχαιο αριθμο απο 0-1
                        if (r>h){   
                            s[i+j*n]*=-1;
                        }
                    }            
                }
                
            }
           
            
            if (ii>3000){

                double E,m;
                E=energy(n,s,J,B); // ενεργεια καθε σαρωσης
                m=magnitikh(n,s); //μαγνητικη ροπη καθε σαρωσης

                sumE += E ;
                sumE2 += pow(E,2);
                summm += m ; 
                summm2 += pow(m,2);
                population++ ;
            }
        } //τελος σαρωσεων

        meshE = sumE/(population*n*n) ; // (n*n) για να κανονικοποιηθει
        meshm = summm/population;
        out1<< T <<"       "<< meshE/(n*n) <<"         "<< meshm/(n*n) << endl; // στελνει τη μαγνητιση

        
        C = 1/(T*T) * ( (sumE2/population) - pow((sumE/population),2) ); // ειδικη θερμοτητα
        EP = b_const * ( (summm2/population)  - pow((summm/population),2) ); // μαγνητικη επιδεκτικοτητα

        out2<< T <<"         "<< C/(n*n) <<"           " << EP/(n*n) << endl; 
    }
    

    return 0;
}
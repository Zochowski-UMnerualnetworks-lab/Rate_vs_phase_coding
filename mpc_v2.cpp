//
//  mpc_v1.cpp
//  
//
//  Created by roachjp on 6/4/18.
//

#include <iostream>
#include <vector>
#include <iomanip>  // for setw if I use it for output
#include <fstream>  //for ifstream and ofstream
#include <string>  // for the string filename
#include <sstream> // for ostreamstring
#include <stdlib.h> // for Qk sort
#include <cstdlib> // for abs()
#include <math.h> //for sqrt
#include <time.h>

using namespace std;

//#include "mpc_v1.hpp"
#define pi 3.14159

int argc_global;
const char** argv_global;

int main(int argc, char *argv[])
{
    argc_global = argc;
    argv_global = (const char**) argv;
    
    //  --------------------------------------------------------------
    //  1) the input portion
    //  --------------------------------------------------------------
    //  Should read in a user selected file into vectors
    //  one vector will be an int, the second a dbl vector
    //  will also create a count of how many elements in the vector first is 0.
    //  the actual number of elements will be numel+1
    //  numel refers to the positioning in the vector
    
    string filename =argv[1];
    ifstream input(filename.c_str());
    string train, time;
    
    clock_t start,end;
    start = clock();
    
    //initalize vectors s(int) for train number and t(double) for time
    vector <int> S(0);
    S.reserve(40000);
    vector <double> T(0);
    T.reserve(40000);
    int numel=0; // this will give us the spikes in the file
    while (input >> train >> time)
    {   S.push_back(atoi(train.c_str()));
        T.push_back(atof(time.c_str()));
        numel ++;//numel = T.size()
    }
    
    cout << "Spikes: " << S.size() << endl;
    cout << "" << endl;
    
    
    
    //Close ifstream
    input.close(); //check this to see that it works
    input.clear();// this part isn't looped, so I don't need to do this
    
    bool test=true; // this is the test for when you want to run the FCA, FCA ends when this is false
    
    
    // --------------------------------------------------------
    // Calculate MPC
    // --------------------------------------------------------
    
    vector <int> spikes(1); // has the zeroth element equal to 0, no elements before the first train
    spikes.reserve(100);
    int y;
    y=S[0];
    for(int z=0; z<numel; z++)
    {
        if (S[z] != y)
        {
            spikes.push_back(z);
            y=S[z];
        }
        if (z == numel-1) //fence post for the last element in spikes
        {
            spikes.push_back(z+1);
            // the +1 is needed to account for the last line never getting
            // the z++ one last time
        }
    }
    // double delta_t;
    //spikes.size() is one more than the number of trains
    
    double **mpc; //allocate memory, create matrix
    double **phase; //allocate memory, create matrix
    mpc = new double*[spikes.size()-1];
    phase = new double*[spikes.size()-1];
    for (int i=0;i<spikes.size()-1;i++)
    {
        mpc[i]=new double[spikes.size()-1];
        phase[i]=new double[spikes.size()-1];
    }
    double min_mpc = 1.1;
    double max_mpc = 0.0;
    int    max_spks = 0;
    int    min_spks = 0;
    for (int ii=0; ii<spikes.size()-1; ii++)
    {
        for (int jj=0; jj<spikes.size()-1; jj++)
        {
            if (ii==jj)
            {
                mpc[ii][jj] = -1.0;
            }
            else
            {
                double cosSum   = 0.0;
                double sinSum   = 0.0;
                double phaseSum = 0.0;
                int norm_iter = 0;
                int z = 0;
                for (int k=spikes[ii];k<spikes[ii+1];k++)
                {
                    while (T[k]>T[spikes[jj]+z])
                    {
                        if (z==spikes[jj+1]-spikes[jj]-1)
                        {
                            break;
                        }
                        z++;
                    }
                    if (z!=0) // we haven't landed on the first spike of reference
                    {
                        double a = T[k];
                        double b = T[spikes[jj]+z-1]; // reference spike right before T[k]
                        double c = T[spikes[jj]+z];   // reference spike right after T[k]
                        
                        double phi = 2.0*pi*(a-b)/(c-b);
                        phaseSum += phil;
                        cosSum += cos(phi);
                        sinSum += sin(phi);
                        norm_iter++;
                    }
                }
                if (norm_iter > 0)
                {
                    mpc[ii][jj] = sqrt((cosSum/norm_iter)*(cosSum/norm_iter)+(sinSum/norm_iter)*(sinSum/norm_iter));
                    phase[ii][jj] = phaseSum/norm_iter;
                    //cout << mpc[ii][jj] << " " << norm_iter << " " << z << endl;
                    if (mpc[ii][jj] > max_mpc)
                    {
                        max_mpc = mpc[ii][jj];
                        max_spks = norm_iter;
                    }
                    if (mpc[ii][jj] < min_mpc)
                    {
                        min_mpc = mpc[ii][jj];
                        min_spks = norm_iter;
                    }
                }
            }
        }
    }
    cout << "Max MPC is " << max_mpc << " based on " << max_spks << " spikes" << endl;
    cout << "Min MPC is " << min_mpc << " based on " << min_spks << " spikes"  << endl;
    
    string file_name = argv_global[1];
    file_name.resize(file_name.size()-4);
    ofstream output;
    string fileout1 = file_name + "_mpc.dat";
    string fileout2 = file_name + "_phase.dat";
    output1.open(fileout1.c_str());
    output2.open(fileout2.c_str());
    
    for(int k=0;k<spikes.size()-1;k++)
    {
        output1 << setprecision(8) << S[spikes[k]] << "\t"; //cell ID
        output2 << setprecision(8) << S[spikes[k]] << "\t"; //cell ID
    }
    output << "\n";
    
    for(int k=0;k<spikes.size()-1;k++)
    {
        output1 << setprecision(5) << spikes[k+1]-spikes[k]<<"\t"; //spikes per cell
        output2 << setprecision(5) << spikes[k+1]-spikes[k]<<"\t"; //spikes per cell
    }
    output << "\n";
    
    for(int i=0;i<spikes.size()-1;i++)
    {
        output << "\n";//add a newline before each new row
        for(int j=0;j<spikes.size()-1;j++)
        {
            output1<< mpc[i][j] << "\t" ;
            output2<< phase[i][j] << "\t" ;
        }
    }
    
    output1.close();
    output2.close();
    
    
    end = clock();
    cout << "Total run time = " ;
    cout << double(start-end)/(CLOCKS_PER_SEC)<<endl;
    
    return 0;
}

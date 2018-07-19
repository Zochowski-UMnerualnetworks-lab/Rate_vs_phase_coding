//
//  resonance_singleCELL.cpp
//  
//
//  Created by james roach on 5/10/16.
//
//

#include "network_dyn.hpp"



int main(int argc, char* argv[])
{

    if (argc != 9)
    {
        cout << "wrong number of input arguments" << endl;
        return -1;
    }

    
    double gks_e   = atof(argv[1]);
    double freq_in = atof(argv[2]);
    int xl=1*time(NULL);
    //int xn=25;
    int xn=atoi(argv[3]);
    double dc_spread = atof(argv[4]);
    double Pe        = atof(argv[5]);
    double Pi        = atof(argv[6]);
    double we        = atof(argv[7]);
    double wi        = atof(argv[8]);
    
    
    Ran ICRand(xl);
    Ran NetRan(xn);
    
    
    double gks   = gks_e;
    
    ofstream prams("params.txt");
    prams << gks       << endl;
    prams << freq_in       << endl;
    prams << xl << endl;
    prams << xn << endl;
    prams << dc_spread << endl;
    prams << Pe << endl;
    prams << Pi << endl;
    prams << we << endl;
    prams << wi << endl;
    
    
    prams.close();
    
    double tfinal = 24000.0;
    
    
    double tstep     =   0.05;
    int no_steps     = 1*(tfinal/tstep);
    double time_now;
    // network variables are global
    const int num_ex = 300;
    const int num_ih = 1*0.25*num_ex;
    const int n      = num_ex + num_ih;
    int this_clust;
    double driving_curr = 0.0;
    
    double noise_freq = 0./1000.0; //  Hz/ms !!!
    double noise_STEP_prob = noise_freq*tstep;
    
//    double Pe         = 0.2;
//    double Pi         = 1.0;
    double EI_coorRAT = 4.0;
    double out_conP   = 0.5;
    double out_RAD[2] = {0.03*num_ex,0.06*num_ex}; // E->E,E->I,I->E,I->I (2X disired connections)
    double wint[4] = {we,we,wi,wi}; // E->E,E->I,I->E,I->I

    
    double Ee   = 0.0;
    double Ei   = -70.0;
    
    // need 3 synapses ampa, nmda, gaba.
    
    double tauF = 0.3;
    double tauS = 3.0;
    
    
    int    syn_pl = 1*7*tauS/tstep;
    double syn_pulse[syn_pl];
    double max_syn = 0.0;
    
    

    for (int ii=0; ii<syn_pl; ii++)
    {
        syn_pulse[ii] = exp((-1.0*ii*tstep)/tauS) - exp((-1.0*ii*tstep)/tauF);
        if (syn_pulse[ii] < 0.0)
        {
            syn_pulse[ii] = 0.0;
        }
        if (syn_pulse[ii]>max_syn)
        {
            max_syn = syn_pulse[ii];
        }
    }
    for (int ii=0;ii<syn_pl;ii++)
    {
        syn_pulse[ii] = syn_pulse[ii]/max_syn;
    }

    
    
    double cells_bio[n][6]; // E ring coordinate, I ring coordinate, Idrive,gL, cluster, E/I
    int**    connect = new int*[n];
    double** w_mat = new double*[n];
    int** out_listE = new int*[n];
    int** out_listI = new int*[n];
    int** inn_listE = new int*[n];
    int** inn_listI = new int*[n];
    
    int** out_NOTlistE = new int*[n];
    int** out_NOTlistI = new int*[n];
    
    int out_degE[n];
    int out_degI[n];
    int inn_degE[n];
    int inn_degI[n];
    int out_NOTdegE[n];
    int out_NOTdegI[n];
    
    
    solutions out;
    double V[n];
    double h[n];
    double n_v[n];
    double s[n];
    double Vprev[n];
    double Isyn[n];
    double gL[n];
    double Idrive[n];
    double gsynE[n];
    double gsynI[n];
    int numspikes = 0;
    
    int    last_spike_s[n];
    double last_spike_t[n];
    
    
    string linef;
    string linei;
    string lineg;
    int colf = 0;
    int rowf = 0;
    int coli = 0;
    int rowi= 0;
    int colg = 0;
    int rowg = 0;
    double x_inp;
    double freq[400][400]={{0}};
    double curr[400][400]={{0}};
    double mcond[400][400]={{0}};
    string filename;
    ifstream fileINf;
    ifstream fileINi;
    ifstream fileINg;
    
    //filename = "ifg_mat.txt";
    
    fileINf.open("ifg_mat.txt");
    
    if (fileINf.fail())
    {
        cerr << " * file cannot be found or opened" << endl;
        exit(1);
    }
    
    while (fileINf.good())
    {
        while (getline(fileINf,linef))
        {
            stringstream streamf(linef);
            colf = 0;
            while (streamf >> x_inp)
            {
                freq[rowf][colf] = x_inp;
                colf++;
            }
            rowf++;
        }
    }
    cout << "# of rows [freq] " << rowf << endl;
    cout << "# of columns [freq] " << colf << endl;
    cout << endl;
    
    //----------------------------------------------
    //filename = "curr.txt";
    
    fileINi.open("curr.txt");
    
    if (fileINi.fail())
    {
        cerr << " * file cannot be found or opened" << endl;
        exit(1);
    }
    
    while (fileINi.good())
    {
        while (getline(fileINi,linei))
        {
            stringstream streami(linei);
            coli = 0;
            while (streami >> x_inp)
            {
                curr[rowi][coli] = x_inp;
                coli++;
            }
            rowi++;
        }
    }
    cout << "# of rows [Idrive] " << rowi << endl;
    cout << "# of columns [Idrive] " << coli << endl;
    cout << endl;
    
    //----------------------------------------------
    //filename = "mcond.txt";
    
    fileINg.open("mcond.txt");
    
    if (fileINg.fail())
    {
        cerr << " * file cannot be found or opened" << endl;
        exit(1);
    }
    
    while (fileINg.good())
    {
        while (getline(fileINg,lineg))
        {
            stringstream streamg(lineg);
            colg = 0;
            while (streamg >> x_inp)
            {
                mcond[rowg][colg] = x_inp;
                colg++;
            }
            rowg++;
        }
    }
    cout << "# of rows [gKs] " << rowg << endl;
    cout << "# of columns [gKs] " << colg << endl;
    cout << endl;
    
    
    
    // find I values for chosen frequency
    
    double Istim[colf];
    double I_select[n][colf];
    double f_diff = 1000.0; // very big to start
    int    f_diff_i = 0;
    double f_diff_t = 1001.0;
    
    for (int kk=0; kk<n; kk++)
    {
        for (int ii=0; ii < colf; ii++)
        {
            f_diff   = 1000.0;
            f_diff_i = 0;
            f_diff_t = 1001.0;
            for (int jj=0; jj <rowf; jj++)
            {
                f_diff_t = fabs(freq[jj][ii] - freq_in);
                if (f_diff_t <= f_diff)
                {
                    f_diff   = f_diff_t;
                    f_diff_i = jj;
                }
            }
            I_select[kk][ii] = curr[0][f_diff_i];
        }
    }
    
    int gks_i = 0;
    
    for (int ii=0; ii<colg; ii++)
    {
        if (mcond[0][ii] == gks)
        {
            gks_i = ii;
        }
    }
//    for (int ii=0; ii<colf; ii++)
//    {
//        f_diff   = 1000.0;
//        f_diff_i = 0;
//        f_diff_t = 1001.0;
//        for (int jj=0; jj<rowf; jj++)
//        {
//            f_diff_t = fabs(freq[jj][ii] - freq_stim);
//            if (f_diff_t <= f_diff)
//            {
//                f_diff   = f_diff_t;
//                f_diff_i = jj;
//            }
//        }
//        Istim[ii] = curr[0][f_diff_i];
//        //cout << Istim[ii] << endl;
//    }
//    double I_drive[n];
//    for (int ii=0; ii<n; ii++)
//    {
//        I_drive[ii] = I_select[ii][gks_i];
//    }
    
    
    for (int i=0; i<n; i++)
    {
        connect[i]     = new int[n];
        w_mat[i]       = new double[n];
        out_listE[i]    = new int[n];
        out_listI[i]    = new int[n];
        inn_listE[i]    = new int[n];
        inn_listI[i]    = new int[n];
        out_NOTlistE[i] = new int[n];
        out_NOTlistI[i] = new int[n];
        out_degE[i]     = 0;
        out_degI[i]     = 0;
        inn_degE[i]     = 0;
        inn_degI[i]     = 0;
        out_NOTdegE[i]  = 0;
        out_NOTdegI[i]  = 0;
        
        //Idrive[i] = 0.*NetRan.doub();
        //Idrive[i] = (amp_dc - 0.1*amp_dc) + (0.1*amp_dc)*NetRan.doub();
        //Idrive[i] = (amp_dc - 0.5*amp_dc) + (0.5*amp_dc)*NetRan.doub();
        Idrive[i] = I_select[i][gks_i] + (dc_spread)*NetRan.doub();
        cells_bio[i][2] = Idrive[i];
        gL[i]    = 0.02;
        //gL[i]    = 0.02 + NetRan.doub()*0.005;
        cells_bio[i][3] = gL[i];
        V[i]     = -70.0 + 10.0*ICRand.doub();
        h[i]     =  0.5*ICRand.doub();
        n_v[i]   =  0.5*ICRand.doub();
        s[i]     =  0.5*ICRand.doub();
        Vprev[i] = V[i];
        last_spike_s[i] = -1;
        last_spike_t[i] = -1.0;
    }
    
    // place cells on rings
    int e_placed=0;
    for (int ii=0;ii<num_ex;ii++)
    {
        cells_bio[ii][5]=0.0;
        cells_bio[ii][0]=1.0*e_placed;
        cells_bio[ii][1]=(1.0/EI_coorRAT)*e_placed;
        cells_bio[ii][4]=-1.0;
        e_placed++;
    }
    int i_placed=0;
    for (int ii=num_ex;ii<n;ii++)
    {
        cells_bio[ii][5]=1.0;
        cells_bio[ii][0]=1.0*EI_coorRAT*i_placed;
        cells_bio[ii][1]=1.0*i_placed;
        cells_bio[ii][4]=-1.0;
        i_placed++;
    }
    
    for (int ii=0;ii<num_ex;ii++)
    {
        connect[ii][ii] = 0;
        for (int jj=0;jj<num_ex;jj++)
        {
            if (ii!=jj)
            {
                double dx = fabs(cells_bio[ii][0]-cells_bio[jj][0]);
                if (num_ex-dx < dx)
                {
                    dx = num_ex-dx;
                }
                if ((dx<=out_RAD[0])&&(NetRan.doub()<=out_conP))
                {
                    connect[ii][jj] = 1; // pre --> post
                    out_listE[ii][out_degE[ii]] = jj;
                    out_degE[ii]++;
                }
                else
                {
                    connect[ii][jj] = 0;
                    out_NOTlistE[ii][out_NOTdegE[ii]] = jj;
                    out_NOTdegE[ii]++;
                }
            }
        }
        for (int jj=num_ex;jj<n;jj++)
        {
            double dx = fabs(cells_bio[ii][0]-cells_bio[jj][0]);
            if (num_ex-dx < dx)
            {
                dx = num_ex-dx;
            }
            if ((dx<=out_RAD[0])&&(NetRan.doub()<=out_conP))
            {
                connect[ii][jj] = 1; // pre --> post
                out_listI[ii][out_degI[ii]] = jj;
                out_degI[ii]++;
            }
            else
            {
                connect[ii][jj] = 0;
                out_NOTlistI[ii][out_NOTdegI[ii]] = jj;
                out_NOTdegI[ii]++;
            }
        }
    }
    for (int ii=num_ex;ii<n;ii++)
    {
        for (int jj=0;jj<num_ex;jj++)
        {
            double dx = fabs(cells_bio[ii][0]-cells_bio[jj][0]);
            if (num_ex-dx < dx)
            {
                dx = num_ex-dx;
            }
            if ((dx<=out_RAD[1])&&(NetRan.doub()<=out_conP))
            {
                connect[ii][jj] = 1; // pre --> post
                out_listE[ii][out_degE[ii]] = jj;
                out_degE[ii]++;
            }
            else
            {
                connect[ii][jj] = 0;
                out_NOTlistE[ii][out_NOTdegE[ii]] = jj;
                out_NOTdegE[ii]++;
            }
        }
        for (int jj=num_ex;jj<n;jj++)
        {
            if (ii!=jj)
            {
                double dx = fabs(cells_bio[ii][0]-cells_bio[jj][0]);
                if (num_ex-dx < dx)
                {
                    dx = num_ex-dx;
                }
                if ((dx<=out_RAD[1])&&(NetRan.doub()<=out_conP))
                {
                    connect[ii][jj] = 1; // pre --> post
                    out_listI[ii][out_degI[ii]] = jj;
                    out_degI[ii]++;
                }
                else
                {
                    connect[ii][jj] = 0;
                    out_NOTlistI[ii][out_NOTdegI[ii]] = jj;
                    out_NOTdegI[ii]++;
                }
            }
        }
    }
    
    
    if (Pe>0.0)
    {
        for (int ii=0;ii<num_ex;ii++)
        {
            // randomise vector of nonoutputs
            int eNouts[out_NOTdegE[ii]];
            int iNouts[out_NOTdegI[ii]];
            for (int oo=0;oo<out_NOTdegE[ii];oo++)
            {
                eNouts[oo]=out_NOTlistE[ii][oo];
            }
            for (int oo=0;oo<out_NOTdegI[ii];oo++)
            {
                iNouts[oo]=out_NOTlistI[ii][oo];
            }
            randomize(eNouts,out_NOTdegE[ii],NetRan);
            randomize(iNouts,out_NOTdegI[ii],NetRan);
            // rewire connections
            int erw_cnt =0;
            int irw_cnt =0;
            //int erw_cnt =out_NOTdegE[ii];
            //int irw_cnt =out_NOTdegI[ii];
            for (int oo=0;oo<out_degE[ii];oo++)
            {
                if (NetRan.doub()<Pe)
                {
                    connect[ii][out_listE[ii][oo]] =0;
                    connect[ii][eNouts[erw_cnt]]   =1;
                    erw_cnt++;
                    //erw_cnt--;
                }
            }
            for (int oo=0;oo<out_degI[ii];oo++)
            {
                if (NetRan.doub()<Pe)
                {
                    connect[ii][out_listI[ii][oo]] =0;
                    connect[ii][iNouts[irw_cnt]]   =1;
                    irw_cnt++;
                    //irw_cnt--;
                }
            }
        }
    }
    if (Pi>0.0)
    {
        for (int ii=num_ex;ii<n;ii++)
        {
            // randomise vector of nonoutputs
            int eNouts[out_NOTdegE[ii]];
            int iNouts[out_NOTdegI[ii]];
            for (int oo=0;oo<out_NOTdegE[ii];oo++)
            {
                eNouts[oo]=out_NOTlistE[ii][oo];
            }
            for (int oo=0;oo<out_NOTdegI[ii];oo++)
            {
                iNouts[oo]=out_NOTlistI[ii][oo];
            }
            randomize(eNouts,out_NOTdegE[ii],NetRan);
            randomize(iNouts,out_NOTdegI[ii],NetRan);
            // rewire connections
            int erw_cnt =0;
            int irw_cnt =0;
            //int erw_cnt =out_NOTdegE[ii];
            //int irw_cnt =out_NOTdegI[ii];
            for (int oo=0;oo<out_degE[ii];oo++)
            {
                if (NetRan.doub()<Pi)
                {
                    connect[ii][out_listE[ii][oo]] =0;
                    connect[ii][eNouts[erw_cnt]]   =1;
                    erw_cnt++;
                    //erw_cnt--;
                }
            }
            for (int oo=0;oo<out_degI[ii];oo++)
            {
                if (NetRan.doub()<Pi)
                {
                    connect[ii][out_listI[ii][oo]] =0;
                    connect[ii][iNouts[irw_cnt]]   =1;
                    irw_cnt++;
                    //irw_cnt--;
                }
            }
        }
    }
    
    // remake connection list, write network files
    ofstream cfile("connect.txt");
    ofstream cbfile("cells_bio.txt");
    ofstream wfile("w_initial.txt");
    for (int ii=0;ii<n;ii++)
    {
        out_degE[ii] = 0;
        out_degE[ii] = 0;
        inn_degE[ii] = 0;
        inn_degI[ii] = 0;
        for (int jj=0;jj<num_ex;jj++)
        {
            if (ii!=jj)
            {
                if (connect[ii][jj]==1)
                {
                    if (ii<num_ex)
                    {
                        w_mat[ii][jj] = wint[0];
                    }
                    else
                    {
                        w_mat[ii][jj] = wint[2];
                    }
                    out_listE[ii][out_degE[ii]] = jj;
                    out_degE[ii]++;
                }
                if (connect[jj][ii]==1)
                {
                    inn_listE[ii][inn_degE[ii]] = jj;
                    inn_degE[ii]++;
                }
            }
            cfile << connect[ii][jj] << " ";
            wfile << w_mat[ii][jj] << " ";
        }
        for (int jj=num_ex;jj<n;jj++)
        {
            if (ii!=jj)
            {
                if (connect[ii][jj]==1)
                {
                    if (ii<num_ex)
                    {
                        w_mat[ii][jj] = wint[1];
                    }
                    else
                    {
                        w_mat[ii][jj] = wint[3];
                    }
                    out_listI[ii][out_degI[ii]] = jj;
                    out_degI[ii]++;
                }
                if (connect[jj][ii]==1)
                {
                    inn_listI[ii][inn_degI[ii]] = jj;
                    inn_degI[ii]++;
                }
            }
            cfile << connect[ii][jj] << " ";
            wfile << w_mat[ii][jj] << " ";
        }
        cfile << endl;
        wfile << endl;
        cbfile << cells_bio[ii][0] << " ";
        cbfile << cells_bio[ii][1] << " ";
        cbfile << cells_bio[ii][2] << " ";
        cbfile << cells_bio[ii][3] << " ";
        cbfile << cells_bio[ii][4] << " ";
        cbfile << cells_bio[ii][5] << " ";
        cbfile << endl;
        
    }
    cfile.close();
    cbfile.close();
    wfile.close();
    
    cout << "bloop" << endl;
    //ofstream vfile("V.txt");
    ofstream rast_file("raster_dat.txt");
    //vfile.precision(9);
    rast_file.precision(9);
    vector<double> raster[n];
    for (int i_t=0; i_t<no_steps; i_t++)
    {
        time_now = i_t*tstep;
       
        for (int ii=0;ii<n;ii++)
        {
            if ((Vprev[ii]<0.0)&&(V[ii]>=0.0)&&(time_now>=100.0))
            {
                //rast_file << time_now << " " << ii << endl;
                raster[ii].push_back(time_now);
                last_spike_t[ii] = time_now;
                last_spike_s[ii] = i_t;
                ++ numspikes;
            }
        }
        // calculate synaptic conductances
        for (int ii=0;ii<n;ii++)
        {
            gsynE[ii] = 0.0;
            gsynI[ii] = 0.0;
            for (int jj=0;jj<inn_degE[ii];jj++)
            {
                if (i_t-last_spike_s[inn_listE[ii][jj]]<syn_pl) // just nmda
                {
                    gsynE[ii] = gsynE[ii] + w_mat[inn_listE[ii][jj]][ii]*syn_pulse[i_t-last_spike_s[inn_listE[ii][jj]]];
                }
            }
            for (int jj=0;jj<inn_degI[ii];jj++)
            {
                if (i_t-last_spike_s[inn_listI[ii][jj]]<syn_pl)
                {
                    gsynI[ii] = gsynI[ii] + w_mat[inn_listI[ii][jj]][ii]*syn_pulse[i_t-last_spike_s[inn_listI[ii][jj]]];
                }
            }
            Isyn[ii] = gsynE[ii]*(Ee-V[ii]) + gsynI[ii]*(Ei-V[ii]);
        }
        //vfile << time_now << " "<< I_osc[0] << " ";
        for (int ii=0;ii<n;ii++)
        {
//            out = v_integrate (gks,gL[ii],V[ii],h[ii],n_v[ii],s[ii],I_osc+Idrive[ii]+Isyn[ii],tstep);
            
            driving_curr = Idrive[ii];
            
//            cout << this_clust << " " << driving_curr << endl;
            out = v_integrate (gks,gL[ii],V[ii],h[ii],n_v[ii],s[ii],driving_curr+Isyn[ii],tstep);
            Vprev[ii] = V[ii];
            V[ii]   = out.Vnew;
            h[ii]   = out.hnew;
            n_v[ii] = out.nnew;
            s[ii]   = out.snew;
            if ((Vprev[ii]<0.0)&&(V[ii]<0.0)&&(time_now>=100.0)&&(noise_STEP_prob>ICRand.doub()))
            {
                V[ii] = 20.0;
            }
            //vfile << V[ii] << " ";
        }
        //vfile << endl;
        //vfile << time_now << " " << I_osc[0]+amp_dc+dc_spread << " " << V[0] << endl;
        if (fmod(time_now,500.0)==0.0)
        {
            cout << time_now << " " << numspikes << endl;
        }
    }
    cout << numspikes << endl;
    for (int ii=0;ii<n;ii++)
    {
        for (int jj=0;jj<raster[ii].size();jj++)
        {
            rast_file << ii << " " << raster[ii][jj] << endl;
        }
    }
    //vfile.close();
    rast_file.close();
}

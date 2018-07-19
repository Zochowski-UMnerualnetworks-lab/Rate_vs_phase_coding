//
//  resonance_singleCELL.hpp
//  
//
//  Created by james roach on 5/10/16.
//
//

#ifndef resonance_singleCELL_hpp
#define resonance_singleCELL_hpp

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <string>
#include <cstdio>
#include <sstream>
//#include <array>
#include <iomanip>
#include <vector>
#include <stdlib.h>
#include <algorithm>

#include <stdio.h>

using namespace std;
typedef unsigned long long int Ullong;
typedef double Doub;
typedef unsigned int Uint;
#define pi 3.14159

double r4_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R4_MIN returns the minimum of two R4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R4_MIN, the minimum of X and Y.
//
{
    double value;
    
    if ( y < x )
    {
        value = y;
    }
    else
    {
        value = x;
    }
    return value;
}


struct Ran {
    
    Ullong u,v,w;
    Ran(Ullong j) : v(4101842887655102017LL), w(1) {
        u = j ^ v; int64();
        v = u; int64();
        w = v; int64();
    }
    inline unsigned long long int64() {
        u = u * 2862933555777941757LL + 7046029254386353087LL;
        v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
        w = 429457665U * (w & 0xffffffff) + (w >> 32);
        Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
        return (x+v)^w;
    }
    inline Doub doub() { return 5.42101086242752217E-20 * int64();}
    
    inline Uint int32() { return (Uint)int64(); }
    
    
};

int random_int( Ran &uni_gen, int low, int high )
{
    double ld  = 1.0*low;
    double hd = 1.0*high;
    double y = (ld + (uni_gen.doub()*(hd-ld))) ;//+ 0.5;
    int x = (int) y;
    return(x);
}

int binornd ( int n, double pp, Ran &uni_gen )

//****************************************************************************80
//
//  Purpose:
//
//    IGNBIN generates a binomial random deviate.
//
//  Discussion:
//
//    This procedure generates a single random deviate from a binomial
//    distribution whose number of trials is N and whose
//    probability of an event in each trial is P.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    Original FORTRAN77 version by Barry Brown, James Lovato.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Voratas Kachitvichyanukul, Bruce Schmeiser,
//    Binomial Random Variate Generation,
//    Communications of the ACM,
//    Volume 31, Number 2, February 1988, pages 216-222.
//
//  Parameters:
//
//    Input, int N, the number of binomial trials, from which a
//    random deviate will be generated.
//    0 < N.
//
//    Input, double PP, the probability of an event in each trial of
//    the binomial distribution from which a random deviate is to be generated.
//    0.0 < PP < 1.0.
//
//    Output, int IGNBIN, a random deviate from the
//    distribution.
//
{
    double al;
    double alv;
    double amaxp;
    double c;
    double f;
    double f1;
    double f2;
    double ffm;
    double fm;
    double g;
    int i;
    int ix;
    int ix1;
    int k;
    int m;
    int mp;
    double p;
    double p1;
    double p2;
    double p3;
    double p4;
    double q;
    double qn;
    double r;
    double t;
    double u;
    double v;
    int value;
    double w;
    double w2;
    double x;
    double x1;
    double x2;
    double xl;
    double xll;
    double xlr;
    double xm;
    double xnp;
    double xnpq;
    double xr;
    double ynorm;
    double z;
    double z2;
    
    if ( pp <= 0.0 || 1.0 <= pp )
    {
        cerr << "\n";
        cerr << "IGNBIN - Fatal error!\n";
        cerr << "  PP is out of range.\n";
        exit ( 1 );
    }
    
    p = r4_min ( pp, 1.0 - pp );
    q = 1.0 - p;
    xnp = ( double ) ( n ) * p;
    
    if ( xnp < 30.0 )
    {
        qn = pow ( q, n );
        r = p / q;
        g = r * ( double ) ( n + 1 );
        
        for ( ; ; )
        {
            ix = 0;
            f = qn;
            u = uni_gen.doub();;
            
            for ( ; ; )
            {
                if ( u < f )
                {
                    if ( 0.5 < pp )
                    {
                        ix = n - ix;
                    }
                    value = ix;
                    return value;
                }
                
                if ( 110 < ix )
                {
                    break;
                }
                u = u - f;
                ix = ix + 1;
                f = f * ( g / ( double ) ( ix ) - r );
            }
        }
    }
    ffm = xnp + p;
    m = ffm;
    fm = m;
    xnpq = xnp * q;
    p1 = ( int ) ( 2.195 * sqrt ( xnpq ) - 4.6 * q ) + 0.5;
    xm = fm + 0.5;
    xl = xm - p1;
    xr = xm + p1;
    c = 0.134 + 20.5 / ( 15.3 + fm );
    al = ( ffm - xl ) / ( ffm - xl * p );
    xll = al * ( 1.0 + 0.5 * al );
    al = ( xr - ffm ) / ( xr * q );
    xlr = al * ( 1.0 + 0.5 * al );
    p2 = p1 * ( 1.0 + c + c );
    p3 = p2 + c / xll;
    p4 = p3 + c / xlr;
    //
    //  Generate a variate.
    //
    for ( ; ; )
    {
        u = uni_gen.doub() * p4;
        v = uni_gen.doub();
        //
        //  Triangle
        //
        if ( u < p1 )
        {
            ix = xm - p1 * v + u;
            if ( 0.5 < pp )
            {
                ix = n - ix;
            }
            value = ix;
            return value;
        }
        //
        //  Parallelogram
        //
        if ( u <= p2 )
        {
            x = xl + ( u - p1 ) / c;
            v = v * c + 1.0 - fabs ( xm - x ) / p1;
            
            if ( v <= 0.0 || 1.0 < v )
            {
                continue;
            }
            ix = x;
        }
        else if ( u <= p3 )
        {
            ix = xl + log ( v ) / xll;
            if ( ix < 0 )
            {
                continue;
            }
            v = v * ( u - p2 ) * xll;
        }
        else
        {
            ix = xr - log ( v ) / xlr;
            if ( n < ix )
            {
                continue;
            }
            v = v * ( u - p3 ) * xlr;
        }
        k = abs ( ix - m );
        
        if ( k <= 20 || xnpq / 2.0 - 1.0 <= k )
        {
            f = 1.0;
            r = p / q;
            g = ( n + 1 ) * r;
            
            if ( m < ix )
            {
                mp = m + 1;
                for ( i = mp; i <= ix; i++ )
                {
                    f = f * ( g / i - r );
                }
            }
            else if ( ix < m )
            {
                ix1 = ix + 1;
                for ( i = ix1; i <= m; i++ )
                {
                    f = f / ( g / i - r );
                }
            }
            
            if ( v <= f )
            {
                if ( 0.5 < pp )
                {
                    ix = n - ix;
                }
                value = ix;
                return value;
            }
        }
        else
        {
            amaxp = ( k / xnpq ) * ( ( k * ( k / 3.0
                                            + 0.625 ) + 0.1666666666666 ) / xnpq + 0.5 );
            ynorm = - ( double ) ( k * k ) / ( 2.0 * xnpq );
            alv = log ( v );
            
            if ( alv < ynorm - amaxp )
            {
                if ( 0.5 < pp )
                {
                    ix = n - ix;
                }
                value = ix;
                return value;
            }
            
            if ( ynorm + amaxp < alv )
            {
                continue;
            }
            
            x1 = ( double ) ( ix + 1 );
            f1 = fm + 1.0;
            z = ( double ) ( n + 1 ) - fm;
            w = ( double ) ( n - ix + 1 );
            z2 = z * z;
            x2 = x1 * x1;
            f2 = f1 * f1;
            w2 = w * w;
            
            t = xm * log ( f1 / x1 ) + ( n - m + 0.5 ) * log ( z / w )
            + ( double ) ( ix - m ) * log ( w * p / ( x1 * q ))
            + ( 13860.0 - ( 462.0 - ( 132.0 - ( 99.0 - 140.0
                                               / f2 ) / f2 ) / f2 ) / f2 ) / f1 / 166320.0
            + ( 13860.0 - ( 462.0 - ( 132.0 - ( 99.0 - 140.0
                                               / z2 ) / z2 ) / z2 ) / z2 ) / z / 166320.0
            + ( 13860.0 - ( 462.0 - ( 132.0 - ( 99.0 - 140.0
                                               / x2 ) / x2 ) / x2 ) / x2 ) / x1 / 166320.0
            + ( 13860.0 - ( 462.0 - ( 132.0 - ( 99.0 - 140.0
                                               / w2 ) / w2 ) / w2 ) / w2 ) / w / 166320.0;
            
            if ( alv <= t )
            {
                if ( 0.5 < pp )
                {
                    ix = n - ix;
                }
                value = ix;
                return value;
            }
        }
    }
    return value;
}
void swap (int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

void randomize ( int arr[], int n, Ran &uni_gen )
{
    // Use a different seed value so that we don't get same
    // result each time we run this program
    // Start from the last element and swap one by one. We don't
    // need to run for the first element that's why i > 0
    for (int i = n-1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        int j = random_int(uni_gen, 0, i );
        //cout << i << " " << j << endl;
        // Swap arr[i] with the element at random index
        swap(&arr[i], &arr[j]);
    }
}

struct solutions
{
    double Vnew;
    double hnew;
    double nnew;
    double snew;
};

double minf (double V)
{
    double x;
    x=1.0/(1.0+exp( (-30.0-V)/9.5 ));
    return(x);
}


double Vdot (double V, double h, double n, double s, double Iext, double C, double gNa, double gKdr, double gKs, double gL, double VNa, double VK, double VL)
{
    double x;
    x=(1/C)*( -gNa*pow((minf(V)),3)*h*(V-VNa)-gKdr *pow(n,4)*(V-VK)-gKs*s*(V-VK) -gL*(V-VL)+Iext);
    return(x);
}

double hinf_ach (double V)
{
    double x;
    x=1.0/(1.0+exp( (53.0+V)/7.0 ));
    return(x);
}

double tau_h (double V)
{
    double x;
    x=0.37+2.78/(1.0+exp( (40.5+V)/6.0 ));
    return(x);
}

double hdot (double V, double h)
{
    double x;
    x=(hinf_ach(V)-h)/tau_h(V);
    return(x);
}

double tau_n (double V)
{
    double x;
    x=0.37+1.85/(1.0+exp( (27.0+V)/15.0 ));
    return(x);
}

double ninf_ach (double V)
{
    double x;
    x=1.0/(1.0+exp( (-30.0-V)/10.0 ));
    return(x);
}

double ndot (double V, double n)
{
    double x;
    x=(ninf_ach(V)-n)/tau_n(V);
    return(x);
}

double sinf_ach (double V)
{
    double x;
    x=1.0/(1.0+exp((-39.0-V)/5.0 ));
    return(x);
}

double sdot(double V, double s)
{
    double x;
    x=(sinf_ach(V)-s)/75.0;
    return(x);
}


solutions v_integrate (double gKs_temp, double gL, double V, double h, double n, double s, double Iext, double tstep)
{
    double C     = 1.0; //microF/cm^2
    double gNa   = 24.0; //mS/cm^2
    double gKdr  = 3.0; //mS/cm^2
    //double gKs   = 0.0; // 0 for type 1; 1.5 for type 2 (in mS/cm^2)
    double gKs = gKs_temp;
//    double gL    = 0.02; //mS/cm^2
    double VNa   = 55.0; //mV
    double VK    = -90.0; //mV
    double VL    = -60.0; //mV
    double Vk1, hk1, nk1, sk1, tstepovertwo;
    double Vk2, hk2, nk2, sk2;
    double Vk3, hk3, nk3, sk3;
    double Vk4, hk4, nk4, sk4;
    double V_new, h_new, n_new, s_new, tstepover6;
    solutions out;
    
    //use Runge-Kutta algorithm to integrate forward in time
    tstepovertwo=tstep/2.0;
    Vk1=Vdot(V,h,n,s,Iext,C,gNa,gKdr,gKs,gL,VNa,VK,VL);
    hk1=hdot(V,h);
    nk1=ndot(V,n);
    sk1=sdot(V,s);

    Vk2=Vdot(V+tstepovertwo*Vk1,h+tstepovertwo*hk1,n+tstepovertwo*nk1,s+tstepovertwo*sk1,Iext,C,gNa,gKdr,gKs,gL,VNa,VK,VL);
    hk2=hdot(V+tstepovertwo*Vk1,h+tstepovertwo*hk1);
    nk2=ndot(V+tstepovertwo*Vk1,n+tstepovertwo*nk1);
    sk2=sdot(V+tstepovertwo*Vk1,s+tstepovertwo*sk1);

    Vk3=Vdot(V+tstepovertwo*Vk2,h+tstepovertwo*hk2,n+tstepovertwo*nk2,s+tstepovertwo*sk2,Iext,C,gNa,gKdr,gKs,gL,VNa,VK,VL);
    hk3=hdot(V+tstepovertwo*Vk2,h+tstepovertwo*hk2);
    nk3=ndot(V+tstepovertwo*Vk2,n+tstepovertwo*nk2);
    sk3=sdot(V+tstepovertwo*Vk2,s+tstepovertwo*sk2);

    Vk4=Vdot(V+tstep*Vk3,h+tstep*hk3,n+tstep*nk3,s+tstep*sk3,Iext,C,gNa,gKdr,gKs,gL,VNa,VK,VL);
    hk4=hdot(V+tstep*Vk3,h+tstep*hk3);
    nk4=ndot(V+tstep*Vk3,n+tstep*nk3);
    sk4=sdot(V+tstep*Vk3,s+tstep*sk3);

    //step the variables forward in time
    tstepover6=tstep/6.0;
    out.Vnew= V+tstepover6*(Vk1+2.0*Vk2+2.0*Vk3+Vk4);
    out.hnew= h+tstepover6*(hk1+2.0*hk2+2.0*hk3+hk4);
    out.nnew= n+tstepover6*(nk1+2.0*nk2+2.0*nk3+nk4);
    out.snew= s+tstepover6*(sk1+2.0*sk2+2.0*sk3+sk4);
    return(out);
}


#endif /* resonance_singleCELL_hpp */


/*
    By J. Ludwig, Copyright 2023. 
*/

/*
    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

*/

/*
    This program calculates the p-value for the Fisher exact test, Chi squared, and G test for a 2x2 contingency table; the Chi square and G test are approximated. This code has not been extensively tested and the approximations could be more accurate. 
*/

#include <iostream>
#include <cmath>
#include <iomanip>      // std::setprecision

long double fisher(long double a, long double b, long double c, long double d);
long double factorial(long double n);
long double logfact(long double n);
long double chi(long double a, long double b, long double c, long double d, long double e, long double f, long double g, long double h, long double n);
long double gtest(long double a, long double b, long double c, long double d, long double e, long double f, long double g, long double h, long double n);



int main(int argc, char* argv[]){
    int a{std::stoi(argv[1])};
    int b{std::stoi(argv[2])};
    int c{std::stoi(argv[3])};
    int d{std::stoi(argv[4])};

    long double e = a + b; // these can be converted to int and other code verified that math involving these numbers doesn't need to be cast for the decimal
    long double f = c + d;
    long double g = a + c;
    long double h = b + d;
    long double n = a + b + c + d;

    long double fish_ans;
    fish_ans = fisher(a, b, c, d);
    std::cout << "Fisher exact test p-value" << "\t" << fish_ans << '\n';

    long double chi_p = chi(a, b, c, d, e, f, g, h, n);
    std::cout << "Chi square p-value" << "\t" << chi_p << '\n';

    long double gt = gtest(a, b, c, d, e, f, g, h, n);
    std::cout << "G test p-value" << "\t" << gt << '\n';
}

long double fisher(long double a, long double b, long double c, long double d){

    int e = a + b; // repeat from main, can be passed from main and removed from here but it matters none
    int f = c + d;
    int g = a + c;
    int h = b + d;
    int n = a + b + c + d + 1;

    int gm; //group min
    int sm; // single min
    if(a < b && a < c && a < d){ // find min and find which combination has the min
        sm = a;
        if(e < g){
            gm = 1;  // a and e low
        }else{
            gm = 2;  // a and g low
        }
    }else if(b < c && b < d){
        sm = b;
        if(e < h){
            gm = 3;  // b and e low
        }else{
            gm = 4;  // b and h low
        }
    }else if(c < d){
        sm = c;
        if(f < g){
            gm = 5;  // c and f low
        }else{
            gm = 6;  // c and g low
        }
    }else{
        sm = d;
        if(f < h){
            gm = 7;  // d and f low
        }else{
            gm = 8;  // d and h low
        }
    }

    int na;  // new a
    int nb;
    int nc;
    int nd;

    long double lnfc[n];  // Natural log factorial
    for(int i = 0; i < n; ++i){
        long double temp = i + 1;
        long double temp2 = lgamma(temp); 
        lnfc[i] = temp2;
        //cout << i << "\t" << lnfc[i] << endl;
    }
    /*
    I found lgamma to be more accuate than the other methods I tried.
    Tried different ways to find factorial for Fisher exact test. (C++ lgamma, writing own function that sums log(n) from 1 to n, writing own Ramanujan approx, writing own Sterling approx).
    Using Wolframalpha for n=2000 as the most correct answer (with precision and accuracy set to 20): lgamma was the most correct (16 decimal places), and my own sum log(n) function was second (14 decimal places).

    long double logfact(long double n){
        long double ans{0};
        do {
            //cout <<std::setprecision(50)<< ans << endl;
            ans += log(n);
            --n;
        } while(n > 0);
        return ans;
    }

    */

    if(gm == 1){
        na = 0;
        nb = e;
        nc = g;
        nd = h - e;
    }else if(gm == 2){
        na = 0;
        nb = e;
        nc = g;
        nd = f - g;
    }else if(gm == 3){
        na = e;
        nb = 0;
        nc = g - e;
        nd = h;
    }else if(gm == 4){
        na = e;
        nb = 0;
        nc = f - h;
        nd = h;
    }else if(gm == 5){
        na = g;
        nb = h - f;
        nc = 0;
        nd = f;
    }else if(gm == 6){
        na = g;
        nb = e - g;
        nc = 0;
        nd = f;
    }else if(gm == 7){
        na = g - f;
        nb = h;
        nc = f;
        nd = 0;
    }else if(gm == 8){
        na = e - h;
        nb = h;
        nc = f;
        nd = 0;
    }

    long double comp;
    long double n_lnfc = lnfc[(int)a + (int)b + (int)c + (int)d];  // N (toal) Natural log factorial
    //std::cout << n_lnfc << endl;
    long double start_p = lnfc[(int)a + (int)b] + lnfc[(int)a + (int)c] + lnfc[(int)c + (int)d] + lnfc[(int)b + (int)d] - lnfc[(int)a] - lnfc[(int)b] - lnfc[(int)c] - lnfc[(int)d] - n_lnfc;
    long double ans = exp(start_p); // Answer

    do{ // loop through more extreme
        //std::cout << cout << na << "\t" << nb << endl;
        //std::cout << nc << "\t" << nd << endl;
        comp = lnfc[na + nb] + lnfc[na + nc] + lnfc[nc + nd] + lnfc[nb + nd] - lnfc[na] - lnfc[nb] - lnfc[nc] - lnfc[nd] - n_lnfc;
        //std::cout << lnfc[na + nb] << "\t" << lnfc[na + nc] << "\t" << lnfc[nc + nd]<< "\t" << lnfc[nb + nd] << "\t" << lnfc[na] << "\t" << lnfc[nb] << "\t" << lnfc[nc] << "\t" << lnfc[nd] << "\t" << n_lnfc <<endl;
        //std::cout << "TOP" << lnfc[na + nb] + lnfc[na + nc] + lnfc[nc + nd] + lnfc[nb + nd] << endl;
        //std::cout << "BOT" << -lnfc[na] - lnfc[nb] - lnfc[nc] - lnfc[nd] - n_lnfc << endl;
        //std::cout << comp << "\t" << setprecision(20) << exp(comp) << endl;
        if(comp < start_p){
            ans += exp(comp);
        }
        if(gm == 3 || gm == 4 || gm == 5 || gm == 6){
            --na;
            ++nb;
            ++nc;
            --nd;
        }else{
            ++na;
            --nb;
            --nc;
            ++nd;
        }
    }while(comp < start_p);

    if(gm == 1){ // set up for two sided other more extreme
        na = e;
        nb = 0;
        nc = g - e;
        nd = h;
    }else if(gm == 2){
        na = g;
        nb = e - g;
        nc = 0;
        nd = f;
    }else if(gm == 3){
        na = 0;
        nb = e;
        nc = g;
        nd = h - e;
    }else if(gm == 4){
        na = e - h;
        nb = h;
        nc = f;
        nd = 0;
    }else if(gm == 5){
        na = g - f;
        nb = h;
        nc = f;
        nd = 0;
    }else if(gm == 6){
        na = 0;
        nb = e;
        nc = g;
        nd = f - g;
    }else if(gm == 7){
        na = g;
        nb = h - f;
        nc = 0;
        nd = f;
    }else if(gm == 8){
        na = e;
        nb = 0;
        nc = f- h;
        nd = h;
    }

    //std::cout << "REAL BOT" << "\t" << n_lnfc << endl;
    //cout << "START P " << "\t" << start_p <<  "\t" << exp(start_p) << endl;
    do{ // loop through more extreme on two tail
        //std::cout << na << "\t" << nb << endl;
        //std::cout << nc << "\t" << nd << endl;
        comp = lnfc[na + nb] + lnfc[na + nc] + lnfc[nc + nd] + lnfc[nb + nd] - lnfc[na] - lnfc[nb] - lnfc[nc] - lnfc[nd] - n_lnfc;
        //std::cout << lnfc[na + nb] << "\t" << lnfc[na + nc] << "\t" << lnfc[nc + nd]<< "\t" << lnfc[nb + nd] << "\t" << lnfc[na] << "\t" << lnfc[nb] << "\t" << lnfc[nc] << "\t" << lnfc[nd] << "\t" << n_lnfc <<endl;
        //std::cout << "TOP" << lnfc[na + nb] + lnfc[na + nc] + lnfc[nc + nd] + lnfc[nb + nd] << endl;
        //std::cout << "BOT" << -lnfc[na] - lnfc[nb] - lnfc[nc] - lnfc[nd] - n_lnfc << endl;
        //std::cout << comp << "\t" << setprecision(20) << exp(comp) << endl;
        if(comp < start_p){
            ans += exp(comp);
        }
        if(gm == 3 || gm == 4 || gm == 5 || gm == 6){
            ++na;
            --nb;
            --nc;
            ++nd;
        }else{
            --na;
            ++nb;
            ++nc;
            --nd;
        }
    }while(comp < start_p);

  return ans;

}

// Unused, found lgamma to be more accurate 
long double logfact(long double n){
    long double ans{0};
    do {
        //std::cout <<std::setprecision(50)<< ans << endl;
        ans += log(n);
        --n;
    }while(n > 0);
  return ans;
}

//long double factorial(long double n){
//  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
//}

long double chi(long double a, long double b, long double c, long double d, long double e, long double f, long double g, long double h, long double n){
    std::cout << std::setprecision(10);
    //std::cout << a << "\t" << b << "\t" << c << "\t" << d << "\t" << e << "\t" << f << "\t" << g << "\t" << h << "\t" << n << endl;

    long double chi_sum = (long double)n * ((a * d) - (b * c)) * ((a * d) - (b * c)) / (e * f * g * h);
    long double chi_p;
    long double Ha{-1.37266}; // Hoaglin 1977 aprox
    long double Hb{1.06807};
    long double Hc{2.13161};
    long double Hd{-0.04589};

    // Approx from Beh 2018. Exploring How to Simply Approximate the P-value of a Chi-Squared Statistic. 
    // https://doi.org/10.17713/ajs.v47i3.757 
    // Equation 7 has typo, should not be 0.22 but should be 0.47. 1/(2.13161 - 0.04589) = 0.47
    std::cout << "Chi square value" << '\t' << chi_sum << '\n';
    if(chi_sum < (-1.37266 + 1.06807) * (-1.37266 + 1.06807)){
        chi_p = 1;
    }else{
        chi_p = pow((1/10.), (((sqrt(chi_sum) - (-1.37266 + 1.06807)) / (2.13161 - 0.04589)) * ((sqrt(chi_sum) - (-1.37266 + 1.06807)) / (2.13161 - 0.04589))));
    }
    return chi_p;
}

long double gtest(long double a, long double b, long double c, long double d, long double e, long double f, long double g, long double h, long double n){
    long double na = ((long double)e/n) * g;
    long double nb = ((long double)e/n) * h;
    long double nc = ((long double)f/n) * g;
    long double nd = ((long double)f/n) * h;
    long double Gtest = -2 * ((a * log(na / a)) + (b * log(nb / b)) + (c * log(nc /c)) + (d * log(nd / d)));
    std::cout << "G test value" << '\t' << Gtest << '\n';

    // From Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables. 1964.
    // Gauss error function approximation. The function is the prob that a value falls between two intervals
    long double sum_e = sqrt(Gtest);
    long double y = sqrt(0.5 * Gtest);
    long double S1 = 0.0706230784 * y;
    long double S2 = 0.0422820123 * y * y;
    long double S3 = 0.0092705272 * pow(y, 3);
    long double S4 = 0.0001520143 * pow(y, 4);
    long double S5 = 0.0002765672 * pow(y, 5);
    long double S6 = 0.0000430638 * pow(y, 6);
    long double TEG = 1 / pow((1 + S1 + S2 + S3 + S4 + S5 + S6), 16);
    //std::cout << TEG << endl;
    //std::cout << 1 - TEG - sqrt(2 / M_PI) * (1 / exp(0.5 * Gtest)) * sum_e << endl;
    //std::cout << (sqrt(0.5 * Gtest) * (1 / exp(0.5 * Gtest))) << endl;
    //std::cout << (sqrt(0.5 * Gtest) * (1 / exp(0.5 * Gtest))) / 0.886225 << endl;
    long double P = TEG - sqrt(2 / M_PI) * (1 / exp(0.5 * Gtest)) * sum_e + (sqrt(0.5 * Gtest) * (1 / exp(0.5 * Gtest))) / 0.886225;
    //long double P = 1 - TEG - sqrt(2 / M_PI) * (1 / exp(0.5 * Gtest)) * sum_e + (sqrt(0.5 * Gtest) * (1 / exp(0.5 * Gtest))) / 0.886225;
    return P;
}


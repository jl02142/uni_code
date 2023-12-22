/*
    By J. Ludwig, Copyright 2018. 
*/

/*
    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

*/

/*
    This program was written for the alg class and clustered data in a file based on GC content and temperature.
*/

#include <iostream>
#include <list>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <math.h>    // for sqrt
#include <time.h>    // for random seed
#include <iomanip>   // for setprecision

using namespace std;

void euc_dist(int i, int j, int col_counter, double* euc_dist_arr[], double* centroid_arr[], string* data_arr[]);
int checker(int cluster_1_arr[], int cluster_2_arr[], int row_counter);

int main(int argc, char* argv[]){
    cout << setprecision(5);  //number of digits to print in output
    ifstream file (argv[1]);
    //ifstream file ("Bacteria+Archaea.txt"); // hard coded file name

    //Start of column count
    string temp;
    int col_counter{0};
    getline(file, temp, '\n');
    for(int i = 0; i < temp.length(); ++i){
        if(temp.at(i) == '\t'){
            ++col_counter;
        }
    }
    // end of column count

    //start of row count
    int row_counter{0};
    string line_length;
    while (getline(file, line_length)){
        row_counter++;
    }
    //end of row count

    file.close();

    ifstream file2 (argv[1]);
    //ifstream file2 ("Bacteria+Archaea.txt");

    //dynamically allocate a 2D array
    string** data_arr = new string*[row_counter];
    for(int i = 0; i < row_counter; ++i){
        data_arr[i] = new string[col_counter];
    }
    // end array allocation

    //array to hold header line
    string * header = new string[col_counter+1];
    //end array to hold header line
    //array to hold sample name
    string * sample = new string[row_counter];
    // end array to hold sample name

    //store header line
    //start reading in data
    string temp2;
    getline(file2,temp2,'\n');
    stringstream ss(temp2);
    int temp_i{0};
    while(getline(ss, temp2, '\t')){
        header[temp_i] = temp2;
        ++temp_i;
    }
    //end store header line

    //begin reading data
    int row_counter2{0};
    while(getline(file2, temp2, '\n')){
        int col_counter2{0};
        string flag = "first";
        stringstream ss(temp2);
        while(getline(ss, temp2, '\t')){
            if(flag == "first"){
                sample[row_counter2] = temp2;
                flag = "second";
            } else {
                data_arr[row_counter2][col_counter2] = temp2;
                ++col_counter2;
            }
        }
        ++row_counter2;
    }
    // end reading data

    // start normalization
    for(int i = 0; i < col_counter; ++i){
        double sum{0};
        double count_of{0};
        double mean{0};
        double var_sum{0};
        double variance{0};
        double stand_dev{0};
        for(int j = 0; j < row_counter; ++j){
            sum = sum + stod(data_arr[j][i]);
            ++count_of;
        }
        mean = sum/count_of;
        for(int m = 0; m < row_counter; ++m){
            var_sum = var_sum + ((stod(data_arr[m][i])-mean)*(stod(data_arr[m][i])-mean));
        }
        variance = var_sum/count_of;
        stand_dev = sqrt(variance);
        for(int h = 0; h < row_counter; ++h){
            data_arr[h][i] = to_string(((stod(data_arr[h][i])-mean)/stand_dev));
        }
    }
    // end normalization

    list <int> WCSS_list;
    list <int> AIC_list;
    int min_k{strtol(argv[2], NULL, 10)};
    int max_k{strtol(argv[3], NULL, 10)};
    srand(time(NULL));
    for(int k = min_k; k < max_k+1; ++k){
        string outfile = string(argv[4]) + to_string(k);
        ofstream myfile (outfile);
        //start dynamically allocate centroid array
        double ** centroid_arr = new double*[k];
        for(int i = 0; i < k; ++i){
            centroid_arr[i] = new double[col_counter];
        }
        //end centroid array

        //start euclidean distance array
        double ** euc_dist_arr = new double*[row_counter];
        for(int i = 0; i < row_counter; ++i){
            euc_dist_arr[i] = new double[k];
        }
        //end ecu distance array

        //start centroid k-mean++
        int c1 = rand()%((row_counter + 1) - 0) + 0;//centroid 1
        int c2 = rand()%((row_counter + 1) - 0) + 0;
        for(int i = 0; i < col_counter; ++i){
            centroid_arr[0][i] = stod(data_arr[c1][i]);
        }
        // array to store number line
        double** number_line_arr = new double*[row_counter];
        for(int t = 0; t < row_counter; ++t){
            number_line_arr[t] = new double[k];
        }
        //end number line array

        for(int i = 0; i < k; ++i){
            for(int j = 0; j < row_counter; ++j){
                euc_dist(i, j, col_counter, euc_dist_arr, centroid_arr, data_arr);
            }
            if(i < k-1){
                double weight_max{0};
                for(int m = 0; m < row_counter; ++m){
                    weight_max = weight_max + euc_dist_arr[m][i];
                }
                for(int n = 0; n < row_counter; ++n){
                    euc_dist_arr[n][i] = euc_dist_arr[n][i] / weight_max;
                    //cout << "euc_dist_arr " << euc_dist_arr[n][i] << "\n";
                }
                double number_line{0};
                for(int y = 0; y < row_counter; ++y){
                    number_line = number_line + euc_dist_arr[y][i];
                    //cout << "numer line " << number_line << "\n";
                    number_line_arr[y][i] = euc_dist_arr[y][i] + number_line;
                    //cout << "euc_dist_arr" << euc_dist_arr[y][i] << "\n";
                }
                double rand_weight_num = ((double) rand() / (RAND_MAX)); // see https://stackoverflow.com/questions/9878965/rand-between-0-and-1
                //cout << "random " << rand_weight_num << "\n";
                for (int r = 0; r < row_counter; ++r){
                    if(number_line_arr[r][i] > rand_weight_num){
                        for(int w = 0; w < col_counter; ++w){
                            centroid_arr[i+1][w] = stod(data_arr[r][w]);
                        }
                        break;
                    }
                }
            }

        }
        // end initial centroids

        //start create array to store what cluster each sample will belong to
        int* cluster_1_arr = new int[row_counter];
        int* cluster_2_arr = new int[row_counter];
        //end
        string flag = "first";
        int tracker{0};

        //start clustering
        string final;
        do{
            if(flag == "first"){
                flag = "second";
                for(int i = 0; i < row_counter; ++i){
                    double min_val{0};
                    int min_index{-1};
                    for(int e = 0; e < k; ++e){
                        if(min_index == -1){
                            //cout << "!!!!! " << euc_dist_arr[i][e] << "\n";
                            min_val = euc_dist_arr[i][e];
                            min_index = e;
                            cluster_1_arr[i] = min_index;
                        } else {
                            //cout << "@@@@" << euc_dist_arr[i][e] << "\n";
                            if(min_val > euc_dist_arr[i][e]){
                                min_val = euc_dist_arr[i][e];
                                min_index = e;
                                cluster_1_arr[i] = min_index;
                            }
                        }
                    }

                }
            }else{
                flag = "first";
                for(int i = 0; i < row_counter; ++i){
                    double min_val{0};
                    int min_index{-1};
                    for(int e = 0; e < k; ++e){
                        if(min_index == -1){
                            //cout << "!!!!! " << euc_dist_arr[i][e] << "\n";
                            min_val = euc_dist_arr[i][e];
                            min_index = e;
                            cluster_2_arr[i] = min_index;
                        } else {
                            //cout << "@@@@@" << euc_dist_arr[i][e] << "\n";
                            if(min_val > euc_dist_arr[i][e]){
                                min_val = euc_dist_arr[i][e];
                                min_index = e;
                                cluster_2_arr[i] = min_index;
                            }
                        }
                    }
                }
            }
            ++tracker; // how many times the clustering loop executes

            if(flag == "second"){
                final = "1";
                for(int m = 0; m < k; ++m){
                    list <int> g1;
                    for(int f = 0; f < row_counter; ++f){
                        if (cluster_1_arr[f] == m){
                            g1.push_back(f);
                        }
                    }
                    for(int w = 0; w < col_counter; ++w){
                        double cent_col_sum{0};
                        double cent_col_count{0};
                        for(auto v : g1){
                            cent_col_sum = cent_col_sum + stod(data_arr[v][w]);
                            ++cent_col_count;
                        }
                        centroid_arr[m][w] = cent_col_sum / cent_col_count;
                    }
                    if(g1.size() == 0){
                        int c2 = rand()%((row_counter + 1) - 0) + 0;//centroid 1
                        for(int i = 0; i < col_counter; ++i){
                            centroid_arr[m][i] = stod(data_arr[c2][i]);
                        }
                    }
                    for(int n = 0; n < row_counter; ++n){
                        euc_dist(m, n, col_counter, euc_dist_arr, centroid_arr, data_arr);
                    }
                }
                } else {
                    final = "2";
                    for(int m = 0; m < k; ++m){
                        list <int> g1;
                        for(int f = 0; f < row_counter; ++f){
                            if (cluster_2_arr[f] == m){
                                g1.push_back(f);
                            }
                        }
                        for(int w = 0; w < col_counter; ++w){
                            double cent_col_sum{0};
                            double cent_col_count{0};
                            for(auto v : g1){
                                cent_col_sum = cent_col_sum + stod(data_arr[v][w]);
                                ++cent_col_count;
                            }
                            centroid_arr[m][w] = cent_col_sum / cent_col_count;
                        }
                        if(g1.size() == 0){
                            int c2 = rand()%((row_counter + 1) - 0) + 0;//centroid 1
                            for(int i = 0; i < col_counter; ++i){
                                centroid_arr[m][i] = stod(data_arr[c2][i]);
                            }
                        }
                        for(int n = 0; n < row_counter; ++n){
                            euc_dist(m, n, col_counter, euc_dist_arr, centroid_arr, data_arr);
                        }
                    }
            }

        }while(checker(cluster_1_arr, cluster_2_arr, row_counter) == 0);

        // end clustering
        double WCSS{0};
        if(final == "1"){
            for(int m = 0; m < k; ++m){
                list <int> g2;
                double within_sum{0};
                for(int f = 0; f < row_counter; ++f){
                    if (cluster_1_arr[f] == m){
                        g2.push_back(f);
                    }
                }
                //myfile << "!!!!!!!!!!!!!!!!!!!!!!!!" << "\n";
                for(auto v : g2){
                    //cout << "v1 " << v << "\n";
                    within_sum = within_sum + pow(euc_dist_arr[v][m], 2);
                    //cout << "euc dist " << euc_dist_arr[v][m] << "\n";
                    myfile << euc_dist_arr[v][m] << "\t" << sample[v] << "\n";
                }
            //cout << "ws pre " << within_sum << "\n";
            WCSS = WCSS + within_sum;
            }
        }else{
            for(int m = 0; m < k; ++m){
                list <int> g2;
                double within_sum{0};
                for(int f = 0; f < row_counter; ++f){
                    if (cluster_2_arr[f] == m){
                        g2.push_back(f);
                    }
                }
                //myfile << "!!!!!!!!!!!!!!!!!!!!!!!!" << "\n";
                for(auto v : g2){
                   //cout << "v2 " << v << "\n";
                   within_sum = within_sum + pow(euc_dist_arr[v][m], 2);
                   //cout << "euc dist " << euc_dist_arr[v][m] << "\n";
                   myfile << euc_dist_arr[v][m] << "\t" << sample[v] << "\n";
                }
            //cout << "ws pre " << within_sum << "\n";
            WCSS = WCSS + within_sum;
            }
        }
        double AIC{0};
        AIC = (2*k*col_counter) + WCSS;
        WCSS_list.push_back(WCSS);
        WCSS_list.push_back(AIC);
    }

    // print table to stdout
    cout << "k" << "\t" << "WCSS" << "\t" << "AIC" << "\n";
    int pc{0};
    int min_k2{strtol(argv[2], NULL, 10)};
    for(auto v : WCSS_list){
        if(pc % 2 == 0){
            cout << min_k2 << "\t" << v << "\t";
            ++min_k2;
        } else {
            cout << v << "\n";
        }
        ++pc;
    }

}

int checker(int cluster_1_arr[], int cluster_2_arr[], int row_counter){
    for(int u = 0; u < row_counter; ++u){
        //cout << "u " << u << "\n";
        //cout << "clust1 " << cluster_1_arr[u] << "\t" << "clust2 " << cluster_2_arr[u] << "\n";
        //cout << "clust2 " << cluster_2_arr[u] << "\n";
        //cout << "new " << "\n";
        if(cluster_1_arr[u] != cluster_2_arr[u]){
            //cout << u << "\t" << cluster_1_arr[u] << "\t" << cluster_2_arr[u] << "\n";
            return 0;
        }
    }
    for(int u = 0; u < row_counter; ++u){
        //cout << cluster_1_arr[u] << "\t" << cluster_2_arr[u] << "\n";
    }
    return 1;
}

void euc_dist(int i, int j, int col_counter, double* euc_dist_arr[], double* centroid_arr[], string* data_arr[]){
    double euc_dist_sum{0};
    for(int p = 0; p < col_counter; ++p){
        //cout << "check " << centroid_arr[prev_centroid][i] << "\n";
        //cout<< "data here " << stod(data_arr[j][p]) << "\n";
        //cout << "cent here " << centroid_arr[i][p] << "\n";
        euc_dist_sum = euc_dist_sum + pow(stod(data_arr[j][p]) - centroid_arr[i][p] ,2);
        //cout << "sum here " << euc_dist_sum << "\n";
    }
    euc_dist_arr[j][i] = sqrt(euc_dist_sum);
}


// J ludwig 2018

// This program takes in two sequences and aligns them using the Needleman-Wunsch algorithm. The program takes in the sequences from a fasta file and outputs the aligned sequences to a file. 

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
// Have to edit fasta files as follows to get program to work. Add a newline character after the sequence name and remove all newline characters from the sequence.
/*
cat $infile_1.fa | sed 's/^\(>.*\)/\1 HERE/' | tr -d '\n' | sed 's/HERE/\n/' | sed -e '$a\' > temp

mv temp $infile_1.fa

 

cat $infile_2.fa | sed 's/^\(>.*\)/\1 HERE/' | tr -d '\n' | sed 's/HERE/\n/' | sed -e '$a\' > temp

mv temp $infile_2.fa

 

cat $infile_1.fa > $infile.fa

cat $infile_2.fa >> $infile.fa

*/

double max_num(double diag, double left, double top);
int dir(double diag, double left, double top);

int main(int argc, char* argv[]){

    //Take match, mismatch and gap from arguments
    int match = stoi(argv[2]);
    int missmatch = stoi(argv[3]);
    int gap = stoi(argv[4]);
    //int match = 1;
    //int missmatch = -1;
    //int gap = -1;
    //

    //Take infile from arguments
    ifstream file (argv[1]);
    //ifstream file ("");
    //

    //Count length of both sequences. Set longer length sequence to be first (s1_data_arr[]).
    int temp_max_seq_len_1{0};
    int temp_max_seq_len_2{0};
    int max_seq_len;
    int min_seq_len;
    string temp;
    string temp2;
    char c;
    getline(file, temp);
    getline(file, temp2);
    stringstream ss (temp2);
    while(ss >> c){
        ++temp_max_seq_len_1;
    }
    getline(file, temp);
    getline(file, temp2);
    stringstream ss2 (temp2);
    while(ss2 >> c){
        ++temp_max_seq_len_2;
    }
    if(temp_max_seq_len_1 > temp_max_seq_len_2){
        max_seq_len = temp_max_seq_len_1;
        min_seq_len = temp_max_seq_len_2;
    } else {
        min_seq_len = temp_max_seq_len_1;
        max_seq_len = temp_max_seq_len_2;
    }
    //cout << max_seq_len << "\n";
    //cout << min_seq_len << "\n";
    file.close();

    ifstream file2 (argv[1]);

    //Read in first sequence
    string* s1_name_arr = new string[1];
    string* s2_name_arr = new string[1];
    string** s1_data_arr = new string*[1];
    string** s2_data_arr = new string*[1];
    s1_data_arr[1] = new string [max_seq_len];
    s2_data_arr[1] = new string [min_seq_len];

    getline(file2, s1_name_arr[0]);
    getline(file2, temp2);
    stringstream ss3 (temp2);
    int temp_counter{0};
    if(temp_max_seq_len_1 > temp_max_seq_len_2){
        while(ss3 >> c){
            s1_data_arr[1][temp_counter] = c;
            ++temp_counter;
        }
    } else {
        while(ss3 >> c){
            s2_data_arr[1][temp_counter] = c;
            ++temp_counter;
        }
    }
    //Read in second sequence
    getline(file2, s2_name_arr[0]);
    getline(file2, temp2);
    stringstream ss4 (temp2);
    int temp_counter2{0};
    if(temp_max_seq_len_1 > temp_max_seq_len_2){
        while(ss4 >> c){
            s2_data_arr[1][temp_counter2] = c;
            ++temp_counter2;
        }
    } else {
        while(ss4 >> c){
            s1_data_arr[1][temp_counter2] = c;
            ++temp_counter2;
        }
    }
    //Allocate score and directional matrix
    ++max_seq_len;
    ++min_seq_len;
    double** score_matrix = new double* [max_seq_len];
    for(int i = 0; i < max_seq_len; ++i){
        score_matrix[i] = new double[min_seq_len];
    }

    double** dir_matrix = new double* [max_seq_len];
    for(int i = 0; i < max_seq_len; ++i){
        dir_matrix[i] = new double[min_seq_len];
    }
    //Fill in first row and column of matrix
    int temp_neg{0 + gap};
    score_matrix[0][0] = 0;
    dir_matrix[0][0] = 0;
    for(int i = 1; i < max_seq_len; ++i){
        score_matrix[i][0] = temp_neg;
        dir_matrix[i][0] = 0;
        temp_neg = temp_neg + gap;
    }
    temp_neg = 0 + gap;
    for(int i = 1; i < min_seq_len; ++i){
        score_matrix[0][i] = temp_neg;
        dir_matrix[0][i] = 0;
        temp_neg = temp_neg + gap;
    }

    //Fill in rest of score and direction matrix
    --max_seq_len;
    --min_seq_len;
    double diag{0};
    double left{0};
    double top{0};
    for(int i = 0; i < max_seq_len; ++i){
        for(int j = 0; j < min_seq_len; ++j){
            if(s1_data_arr[1][i] == s2_data_arr[1][j]){ //match
                diag = score_matrix[i][j] + match;
                left = score_matrix[i+1][j] + gap;
                top = score_matrix[i][j+1] + gap;
                score_matrix[i+1][j+1] = max_num(diag, left, top);
                dir_matrix[i+1][j+1] = dir(diag, left, top);
            } else {                                   //mismatch
                diag = score_matrix[i][j] + missmatch;
                left = score_matrix[i+1][j] + gap;
                top = score_matrix[i][j+1] + gap;
                score_matrix[i+1][j+1] = max_num(diag, left, top);
                dir_matrix[i+1][j+1] = dir(diag, left, top);
            }
        }
    }

    // Count number of times top or left movement occurs in direction matrix
    int tracker{0};
    int t1{max_seq_len};
    int t2{min_seq_len};
    while(t1 != 0 || t2 != 0){
        if(dir_matrix[t1][t2] == 1){
            t1 = t1 - 1;
            t2 = t2 - 1;
        } else if(dir_matrix[t1][t2] == 2){
            ++tracker;
            t2 = t2 - 1;
        } else {
            ++tracker;
            t1 = t1 - 1;
        }
    }


    //Allocate output matrix
    string** output_matrix = new string* [2];
    for(int i = 0; i < 2; ++i){
        output_matrix[i] = new string[max_seq_len + tracker];
    }

    //Retrace through direction matrix to find output
    int seq_length_checker{max_seq_len + tracker - 1};
    int seq_length_checker2{max_seq_len-1};
    int seq_length_checker3{min_seq_len-1};
    t1 = max_seq_len;
    t2 = min_seq_len;
    while(t1 != 0 || t2 != 0){
        if(dir_matrix[t1][t2] == 1){
            output_matrix[0][seq_length_checker] = s1_data_arr[1][seq_length_checker2];
            output_matrix[1][seq_length_checker] = s2_data_arr[1][seq_length_checker3];
            --seq_length_checker;
            --seq_length_checker2;
            --seq_length_checker3;
            t1 = t1 - 1;
            t2 = t2 - 1;
        } else if(dir_matrix[t1][t2] == 2){
            output_matrix[0][seq_length_checker] = "-";
            output_matrix[1][seq_length_checker] = s2_data_arr[1][seq_length_checker3];
            --seq_length_checker;
            output_matrix[0][seq_length_checker] = s1_data_arr[1][seq_length_checker2];
            --seq_length_checker3;
            t2 = t2 - 1;
        } else if(dir_matrix[t1][t2] == 3){
            output_matrix[0][seq_length_checker] =  s1_data_arr[1][seq_length_checker2];
            output_matrix[1][seq_length_checker] = "-";
            --seq_length_checker;
            output_matrix[1][seq_length_checker] = s2_data_arr[1][seq_length_checker3];
            --seq_length_checker2;
            t1 = t1 - 1;
        }
    }

    //Get sequence name for printing
    string seq_name_1;
    string seq_name_2;
    if(temp_max_seq_len_1 > temp_max_seq_len_2){
       seq_name_1 = s1_name_arr[0];
       seq_name_2 = s2_name_arr[0];
    } else {
        seq_name_1 = s2_name_arr[0];
        seq_name_2 = s1_name_arr[0];
    }
    //Print aligned sequences to outfile
    string outfile = argv[5];  //Outfile name from argument
    ofstream myfile (outfile);
    myfile << seq_name_1 <<"\n";
    for(int i = 0; i < max_seq_len + tracker; ++i){
        myfile << output_matrix[0][i];
    }
    myfile << "\n";

    myfile << seq_name_2 <<"\n";
    for(int i = 0; i < max_seq_len + tracker; ++i){
        myfile << output_matrix[1][i];
    }
    myfile << "\n";

    cout << "Score " << score_matrix[max_seq_len][min_seq_len] << "\n";

}

//Function to return max number, returns the max number
double max_num(double diag, double left, double top){
    if(diag >= left && diag >= top){
        return diag;
    } else if(left >= diag && left >= top){
        return left;
    } else {
        return top;
    }
}
//Function to find max number, returns the direction as a number (1 = diag, 2 = left, 3 = top).
int dir(double diag, double left, double top){
    if(diag >= left && diag >= top){
        return 1;
    } else if(left >= diag && left >= top){
        return 2;
    } else {
        return 3;
    }
}

// J Ludwig 2018

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>  // log
#include <iomanip>   // for setprecision

using namespace std;
void printer(int start, int finish, string sample_data[]);
void printer2(int i, int max_seq_length, string* motif_data[]);

int main(int argc, char* argv[])
{

    cout << setprecision(3);  //number of digits to print in output
    cout << showpoint;
    cout << fixed;


    double score = stod((argv[3]));

    ifstream file(argv[1]);                              //Read infile and get size of file
    //ifstream file ("ecoK12-MG1655.fasta");

    string sample_name;
    getline(file, sample_name, '\n');
    string temp;
    char c;
    int temp_char_counter{0};
    while(getline(file, temp, '\n')){
        stringstream ss (temp);
        while(ss >> c){
            ++temp_char_counter;
        }
    }

    file.close();
    long int Acount{0};
    long int Tcount{0};
    long int Gcount{0};
    long int Ccount{0};
    double Totalcount{0};
    ifstream file2(argv[1]);                       //Read query infile and store data
    //ifstream file2 ("ecoK12-MG1655.fasta");
    string* sample_data = new string[temp_char_counter];
    string temp2;
    int temp_counter{0};
    getline(file2, temp2);
    while(getline(file2, temp2, '\n')){
        stringstream ss (temp2);
        while(ss >> c){
            sample_data[temp_counter] = c;
            ++temp_counter;
            ++Totalcount;
            if(c == 'A' || c == 'a'){
                ++Acount;
            } else if(c == 'T' || c == 't'){
                ++Tcount;
            } else if(c == 'G' || c == 'g'){
                ++Gcount;
            } else if(c == 'C' || c == 'c'){
                ++Ccount;
            }
        }
    }
    file2.close();
    double Apercent = Acount / Totalcount;
    double Tpercent = Tcount / Totalcount;
    double Gpercent = Gcount / Totalcount;
    double Cpercent = Ccount / Totalcount;

    //cout << "A count " << Acount << "\n";
    //cout << "Tot count " << Totalcount << "\n";
    //cout << "A% " << Apercent << "\n";

    string temp3;
    ifstream file3(argv[2]);                           //Read motif (database) infile and get size. Can read both fasta and txt files
    //ifstream file3 ("sigma54.fasta");
    getline(file3, temp3, '\n');
    stringstream ss (temp3);
    ss >> c;
    int row_counter{0};
    int max_seq_length{0};
    bool flag_fasta = 0;
    if(c == '>'){
        flag_fasta = 1;
        while(getline(file3, temp3, '\n')){
            int temp_seq_length{0};
            stringstream ss (temp3);
            while(ss >> c){
                ++temp_seq_length;
            }
            getline(file3, temp, '\n');
            ++row_counter;
            if(temp_seq_length > max_seq_length){
                max_seq_length = temp_seq_length;
            }
        }
    } else {
        ++row_counter;
        int temp_seq_length{0};
        while(ss >> c){
            ++temp_seq_length;
        }
        if(temp_seq_length > max_seq_length){
            max_seq_length = temp_seq_length;
        }
        while(getline(file3, temp3, '\n')){
            int temp_seq_length{0};
            stringstream ss (temp3);
            while(ss >> c){
                ++temp_seq_length;
            }
            ++row_counter;
            if(temp_seq_length > max_seq_length){
                max_seq_length = temp_seq_length;
            }
        }
    }

    //cout << "row counter " << row_counter << "\n";
    //cout << "max seq length " << max_seq_length << "\n";

    file3.close();

    ifstream file4(argv[2]);                                //Read motif infile and store data
    //ifstream file4 ("sigma54.fasta");
    string* motif_name = new string[row_counter];
    string** motif_data = new string*[row_counter];
    for(int i = 0; i < row_counter; ++i){
        motif_data[i] = new string[max_seq_length];
    }
    double** freq_matrix = new double*[max_seq_length];
    for(int i = 0; i < max_seq_length; ++i){
        freq_matrix[i] = new double[4]{Apercent, Tpercent, Gpercent, Cpercent};
    }
    temp_counter = 0;
    if(flag_fasta == 1){                                                                           //Two similar loops. One for fasta format and one for other format
        bool alt = 1;
        while(getline(file4, temp3, '\n')){
            if(alt == 1){
                motif_name[temp_counter] = temp3;
                alt = 0;
            } else {
                int temp_counter2{0};
                stringstream ss (temp3);
                while(ss >> c){
                    motif_data[temp_counter][temp_counter2] = c;                                    //Count number of each nucleotide
                    if(c == 'A' || c == 'a'){
                        freq_matrix[temp_counter2][0] = freq_matrix[temp_counter2][0] + 1;
                    } else if(c == 'T' || c == 't'){
                        freq_matrix[temp_counter2][1] = freq_matrix[temp_counter2][1] + 1;
                    } else if(c == 'G' || c == 'g'){
                        freq_matrix[temp_counter2][2] = freq_matrix[temp_counter2][2] + 1;
                    } else if (c == 'C' || c == 'c'){
                        freq_matrix[temp_counter2][3] = freq_matrix[temp_counter2][3] + 1;
                    }
                    ++temp_counter2;
                }
                alt = 1;
                ++temp_counter;
            }
        }
    } else {
        while(getline(file4, temp3, '\n')){
            int temp_counter2{0};
            stringstream ss (temp3);
            while(ss >> c){
                motif_data[temp_counter][temp_counter2] = c;
                if(c == 'A' || c == 'a'){
                    freq_matrix[temp_counter2][0] = freq_matrix[temp_counter2][0] + 1;
                } else if(c == 'T' || c == 't'){
                    freq_matrix[temp_counter2][1] = freq_matrix[temp_counter2][1] + 1;
                } else if(c == 'G' || c == 'g'){
                    freq_matrix[temp_counter2][2] = freq_matrix[temp_counter2][2] + 1;
                } else if (c == 'C' || c == 'c'){
                    freq_matrix[temp_counter2][3] = freq_matrix[temp_counter2][3] + 1;
                }
                ++temp_counter2;
            }
        ++temp_counter;
        }
    }

    /* Uncomment to see motif data printed to stdout
    for(int i = 0; i < row_counter; ++i){
        cout << motif_name[i] << "\n";
        for(int j = 0; j < max_seq_length; ++j){
            cout << motif_data[i][j];
        }
        cout << "\n";
    }
    */

    cout << "\n";                                                           //Print frequency matrix
    cout << "Freq matrix" << "\n";
    cout << "\t\t" << "A \t" << "T \t" << "G \t" << "C \t" << "\n";
    for(int i = 0; i < max_seq_length; ++i){
        cout << "Position " << i << "\t";
        for(int j = 0; j < 4; ++j){
            cout << setw(7) << freq_matrix[i][j] << "\t";
        }
        cout << "\n";
    }
    cout << "\n";


    //cout << "Sample A count " << Acount << "\n";
    //cout << "Sample T count " << Tcount << "\n";
    //cout << "Sample G count " << Gcount << "\n";
    //cout << "Sample C count " << Ccount << "\n";
    //cout << "Sample Total count " << Totalcount << "\n";

    for(int i = 0; i < max_seq_length; ++i){
        double row_sum{0};
        for(int j = 0; j < 4; ++j){
            row_sum = row_sum + freq_matrix[i][j];
        }
        for(int k = 0; k < 4; ++k){
            freq_matrix[i][k] = freq_matrix[i][k] / row_sum;
        }
    }

    //for(int i = 0; i < max_seq_length; ++i){
    //    cout << "\n";
    //    for(int j = 0; j < 4; ++j){
    //        cout << freq_matrix[i][j];
    //    }
    //}

    double** freq_matrix_inverse = new double*[max_seq_length];
    for(int i = 0; i < max_seq_length; ++i){
        freq_matrix_inverse[i] = new double[4];
    }

    int g = 0;
    int f = max_seq_length-1;                                                //Transpose the matrix instead of reverse comp
    while(f > -1){
        double background = freq_matrix[g][0] + freq_matrix[g][1] + freq_matrix[g][2] + freq_matrix[g][3];
        freq_matrix_inverse[f][0] = freq_matrix[g][0] / background;
        freq_matrix_inverse[f][1] = freq_matrix[g][1] / background;
        freq_matrix_inverse[f][2] = freq_matrix[g][2] / background;
        freq_matrix_inverse[f][3] = freq_matrix[g][3] / background;
    }

    cout << "PSSM matrix" << "\n";                                               //Print PSSM
    cout << "\t\t" << "A \t" << "T \t" << "G \t" << "C \t" << "\n";
    for(int i = 0; i < max_seq_length; ++i){
        cout << "Position " << i << "\t";
        cout << setw(7) << log2(freq_matrix[i][0] / Apercent) << "\t";
        cout << setw(7) << log2(freq_matrix[i][1] / Tpercent) << "\t";
        cout << setw(7) << log2(freq_matrix[i][2] / Gpercent) << "\t";
        cout << setw(7) << log2(freq_matrix[i][3] / Cpercent) << "\n";
    }
    cout << "\n";

    for(int i = 0; i < row_counter; ++i){                                   //Print motif score
        int m{0};
        double x{0};
        while(m < max_seq_length){
            string thing = motif_data[i][m];
            if(thing == "A" || thing == "a"){
                x  = x  + log2(freq_matrix[m][0] / Apercent);
            } else if(thing == "T" || thing == "t"){
                x  = x  + log2(freq_matrix[m][1] / Tpercent);
            } else if(thing == "G" || thing == "g"){
                x  = x  + log2(freq_matrix[m][2] / Gpercent);
            } else if(thing == "C" || thing == "c"){
                x  = x  + log2(freq_matrix[m][3] / Cpercent);
            }
            ++m;
        }
        if(x > score){
            cout << "Motif score " << x << " ";

            printer2(i, max_seq_length, motif_data);
            cout << "\n";
        }
    }

    cout << "\n";
    int start = 0;
    int finish = max_seq_length;
    
}


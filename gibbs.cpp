/*
Gibbs sampler
Copyright 2018 J Ludwig
*/

//All seq must be of the same length

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>    // log
#include <iomanip>   // for setprecision

using namespace std;

int main(int argc, char* argv[]){

  //seed counter control to end seeding on line 182, loop counter control to end loop on line 495

  //Uncomment 3 below lines to better control (format) printing matrix. Default has the matrix printing commented out.
  //cout << setprecision(3);  //number of digits to print in output
  //cout << showpoint;
  //cout << fixed;

  int      run_motif = stoi(argv[3]);                       //Initial motif length estimate from user.
  int  ori_motif_len = stoi(argv[3]);                       //Initial motif length estimate from user. Reseed loop will always revert starting motif length back to this number.
  int  low_motif_len = stoi(argv[2]);                       //Low estimate and lower bounds; will not test below this number.
  int high_motif_len = stoi(argv[4]);                       //High estimate and higher bounds; will not test above this number.
  int motif_diff = high_motif_len+1 - low_motif_len;

  int motif_len_range[motif_diff];                          //This loop creates an array that contains each number between the low and high estimate.  Will use later to adjust motif length.
  int temp_motif_counter = low_motif_len;
  for(int i = 0; i < motif_diff; ++i){
      motif_len_range[i] = temp_motif_counter;
      ++temp_motif_counter;
  }

  srand(time(NULL));

  string thing;
  ifstream file(argv[1]);
  //ifstream file("H.pyloriRpoN-sequences-10-300nt.fasta"); //Hardcoded file name
  //ifstream file("E.coliRpoN-sequences-16-100nt.fasta");
  file.seekg(0, ios::end);
  thing.resize(file.tellg());
  file.seekg(0, ios::beg);
  file.read(&thing[0], thing.size());                       //Store the entire contents of the input file into a string called thing.
  file.close();


  stringstream ss (thing);
  char c;
  int seq_count{0};
  int max_seq_len{0};
  int temp_seq_length{0};
  while(ss >> noskipws >> c){
    if(c == '>'){
      ++seq_count;                                          //Count number of sequences by checking for ">".
      if (temp_seq_length > max_seq_len){                   //Checks for longest seq length.  All sequences must be of the same length. Wanted to allow diff length but this was not implemented.
        max_seq_len = temp_seq_length;
      }
      temp_seq_length = 0;
      while(ss >> noskipws >> c && c != '\n'){
        continue;
      }
    }else if(c != '\n'){
      ++temp_seq_length;                                    //Count sequence length.
    }
  }

  string* seq_name_data = new string[seq_count];            //1-D array for sequence names.
  int** seq_nuc_data = new int*[seq_count];                 //2-D array for sequence nucleotides.
  for(int i = 0; i < seq_count; ++i){
    seq_nuc_data[i] = new int[max_seq_len];
  }


  double** background_prob = new double*[seq_count];        //2-D array for background probs. Pseudo counts of 0.25 are used. Very rarely but without pseudo counts - if one nucleotide (e.g. "A") never appears - program would attempt to divide by 0.
  for(int i = 0; i < seq_count; ++i){
    background_prob[i] = new double[4]{0.02, 0.02, 0.02, 0.02};
  }

  ss.str((thing));
  ss.clear();                                               //Reset string "thing"; avoid closing and reopening file.
  int n{0};                                                 //Counter for sequence number; used to store sequence name.
  int e{-1};                                                //Counter for sequence number; used to store nucleotide data. -1 because I couldn't figure out different way to set up loop to keep track of which sequence reading from.
  int m{0};                                                 //Counter for sequence position.
  string temp_name;
  while(ss >> noskipws >> c){

    if(c == '>'){                                           //Store sequence name
      temp_name = c;
      while(ss >> noskipws >> c && c != '\n'){
        temp_name = temp_name + c;
      }
      seq_name_data[n] = temp_name;
      m = 0;
      ++e;
      ++n;
    }else if(c != '\n'){                                    //Store nucleotide. "A" is stored as 0, "T" is 1, "G" is 2, "C" is 3. Compared to string storage, this should consume less memory and also run faster when comparing.
      if(c == 'A'){
        seq_nuc_data[e][m] = 0;
        background_prob[e][0] += 1;
      }else if(c == 'T'){
        seq_nuc_data[e][m] = 1;
        background_prob[e][1] += 1;
      }else if(c == 'G'){
        seq_nuc_data[e][m] = 2;
        background_prob[e][2] += 1;
      }else if(c == 'C'){
        seq_nuc_data[e][m] = 3;
        background_prob[e][3] += 1;
      }
      ++m;
    }

    if(ss.peek() == '>' || ss.peek() == -1){                //Calculate background prob if next character is ">" or NULL (eof).
      for(int k = 0; k < 4; ++k){
        background_prob[e][k] = background_prob[e][k] / m;
      }
    }

  }

  /*                                                        //Uncomment to print background probs.
  cout << "\n";
  cout << "Background Prob" << "\n";
  for(int i = 0; i < seq_count; ++i){
    cout << "Sample " << i << "\t";
    for(int k = 0; k < 4; ++k){
      cout << background_prob[i][k] << "\t";
    }
    cout << "\n";
  }
  cout << "\n";
  */

  int** motif_data_1 = new int*[seq_count];                 //2-D array for motif_1 (current run).
  int** motif_data_2 = new int*[seq_count];                 //2-D array for motif_2 (previous run).
  int** motif_data_3 = new int*[seq_count];                 //2-D array for motif_3 (best run).
  for(int i = 0; i < seq_count; ++i){
    motif_data_1[i] = new int[high_motif_len];
    motif_data_2[i] = new int[high_motif_len];
    motif_data_3[i] = new int[high_motif_len];
  }

  double** freq_matrix = new double*[high_motif_len];       //2-D array for frequency of nucleotides of motifs.
  double** PSSM_matrix = new double*[high_motif_len];       //2-D array for PSSM of motifs.
  double** in_PSSM_mat = new double*[high_motif_len];       //2-D array for inverse of PSSM (instead of reverse comp).
  for(int i = 0; i < high_motif_len; ++i){
    freq_matrix[i] = new double[4];
    PSSM_matrix[i] = new double[4];
    in_PSSM_mat[i] = new double[4];
  }

  int* pos_storage = new int[seq_count];                    //1-D array for position (current run).
  int* pos_storage2 = new int[seq_count];                   //1-D array for position (best run).
  double* score_to_line = new double[2*max_seq_len];        //1-D array for storing PSSM scores.  All scores will be normalized to 1; create a weighted number line from 0 to 1.
  double* score_storage_1 = new double[seq_count];          //1-D array for storing scores (current run).
  double* score_storage_2 = new double[seq_count];          //1-D array for storing scores (previous run).
  double* score_storage_3 = new double[seq_count];          //1-D array for storing scores (best run).

  int print_motif{0};                                       //Store final motif into this.
  int seed_counter{0};
  while(seed_counter < 100){                                //Main data processing / analysis loop.
    int perm_motif = ori_motif_len;                         //Set perm_motif to the original motif provided by user as argument (above).

    for(int j = 0; j < seq_count; ++j){                     //Randomly initialize motif_data from sequence data.  Performed once per seed loop. Max position allowed is no greater than sequence_max_length - motif_length.
      double min_pos = rand() % (max_seq_len - ori_motif_len);
      double max_pos = min_pos + ori_motif_len;
      int a{0};
      for(int i = min_pos; i < max_pos; ++i){
        motif_data_1[j][a] = seq_nuc_data[j][i];
        ++a;
      }
    }

    int loop_counter{0};
    bool loop_flag{0};                                      //Flag to end loop.
    bool adj_mot_len_flag{0};                               //Flag if motif length was adjusted.
    bool adj_mot_pos_flag{0};                               //Flag if motif position was walked.
    while(loop_flag == 0){                                  //Second main data processing / analysis loop.
      ++loop_counter;
                                                            //First run is normal (no motif length adjustment or walking). Second is adjusting motif length. Third is motif position walking. Repeat.
      if(loop_counter % 2 == 0 && loop_counter % 4 != 0){   //Adjust motif length.
        int rand_num = rand() % (motif_diff);               //Random number between 0 and difference between user provided low and high estimate of motif length.
        run_motif = motif_len_range[rand_num];
        for(int j = 0; j < seq_count; ++j){
          int h = pos_storage[j];
          if(h >= max_seq_len - run_motif){                 //If position is greater than max_seq_len - run_motif (no partial motif allowed; no overhang past known nucleotides), set position to last position that still fits run_motif length.
            h = max_seq_len - run_motif;
          }
          for(int a = 0; a < run_motif; ++a){               //Loop to store nucleotide data into motif array with new motif length (run_motif).
            motif_data_1[j][a] = seq_nuc_data[j][h];
            ++h;
          }
        }
        adj_mot_len_flag = 1;
      }else{
        run_motif = perm_motif;
        adj_mot_len_flag = 0;
      }

      if(adj_mot_len_flag == 0 && loop_counter % 4 == 0){   //Walk position along sequence where motif nucleotides were pulled from.
        int rand_num = rand() % 1;                          //Random number 0 or 1.
        if(rand_num == 0){                                  //If random number is 0, walk to the left (-position).
          for(int j = 0; j < seq_count; ++j){
            int pos = pos_storage[j];
            if(pos == 0){                                   //Can't walk if position is already 0 index.
              continue;
            }else if(pos >= max_seq_len - run_motif){       //If position is greater than max_seq_len - run_motif (no partial motif allowed; no overhang past known nucleotides), set position to last position that still fits run_motif length.
              pos = max_seq_len - (run_motif + 1);
            }else if(pos >= (2 * max_seq_len) - run_motif){ //If position is greater than 2*max_seq_len - run_motif (reverse comp)(no partial motif allowed; no overhang past known nucleotides), set position to last position that still fits run_motif length.
              pos = (2 * max_seq_len) - (run_motif + 1);
            }else{
              --pos;                                        //Walk to the left one
            }
            for(int a = 0; a < run_motif; ++a){
              motif_data_1[j][a] = seq_nuc_data[j][pos + a];//Loop to store nucleotide data in motif array with new position.
            }
          }
        }else{                                              //Same as above but walk right (+position).
          for(int j = 0; j < seq_count; ++j){
            int pos = pos_storage[j];
            if(pos >= max_seq_len - run_motif){
              pos = max_seq_len - run_motif;
            }else if(pos >= (2 * max_seq_len) - run_motif){
              pos = (2 * max_seq_len) - run_motif;
            }else{
              ++pos;
            }
            for(int a = 0; a < run_motif; ++a){
              motif_data_1[j][a] = seq_nuc_data[j][pos + a];
            }
          }
        }
        adj_mot_pos_flag = 1;
      }else{
        adj_mot_pos_flag = 0;
      }

      for(int k = 0; k < seq_count; ++k){                   //Third main loop for data processing / analyze.
        for(int a = 0; a < run_motif; ++a){
          for(int m = 0; m < 4; ++m){                       //Set frequency matrix elements to pseudo counts of 0.25. This also clears the frequency matrix each loop of k.
            freq_matrix[a][m] = 0.02;
          }
        }
        for(int j = 0; j < seq_count; ++j){                 //j and k are the sequences.  Checking if j and k  are equal, if they are "skip". Can't use own sequence motif when scoring same sequence.
          if(j == k){
            ++j;
          }
          if(j >= seq_count){                               //If last sequence element, skip. Can't use own sequence motif when scoring same sequence.
            goto skipping;
          }
          for(int a = 0; a < run_motif; ++a){
            if(motif_data_1[j][a] == 0){
              freq_matrix[a][0] += 1;
            }else if(motif_data_1[j][a] == 1){
              freq_matrix[a][1] += 1;
            }else if(motif_data_1[j][a] == 2){
              freq_matrix[a][2] += 1;
            }else if(motif_data_1[j][a] == 3){
              freq_matrix[a][3] += 1;
            }
          }
        }
        skipping:;

        /*                                                  //Uncomment to print frequency matrix.
        cout << "\n";
        cout << "Freq matrix" << "\n";
        cout << "\t\t" << "A \t" << "T \t" << "G \t" << "C \t" << "\n";
        for(int a = 0; a < run_motif; ++a){
          cout << "Position " << a << "\t";
          for(int m = 0; m < 4; ++m){
            cout << freq_matrix[a][m] << "\t";
          }
          cout << "\n";
        }
        cout << "\n";
        */

        for(int a = 0; a < run_motif; ++a){                 //Convert frequency to prob.
          double row_sum{0};
          for(int m = 0; m < 4; ++m){
            row_sum += freq_matrix[a][m];
          }
          for(int m = 0; m < 4; ++m){
            freq_matrix[a][m] = freq_matrix[a][m] / row_sum;
          }
        }

        /*                                                  //Uncomment to print prob matrix.
        cout << "\n";
        cout << "Prob matrix" << "\n";
        cout << "\t\t" << "A \t" << "T \t" << "G \t" << "C \t" << "\n";
        for(int a = 0; a < run_motif; ++a){
          cout << "Position " << a << "\t";
          for(int m = 0; m < 4; ++m){
            cout << freq_matrix[a][m] << "\t";
          }
          cout << "\n";
        }
        cout << "\n";
        */

        for(int a = 0; a < run_motif; ++a){                 //Calculate PSSM matrix.
          for(int m = 0; m < 4; ++m){
            PSSM_matrix[a][m] = log2(freq_matrix[a][m] / background_prob[k][m]);
          }
        }

        /*                                                  //Uncomment to print PSSM.
        cout << "\n";
        cout << "PSSM matrix" << "\n";
        cout << "\t\t" << "A \t" << "T \t" << "G \t" << "C \t" << "\n";
        for(int a = 0; a < run_motif; ++a){
          cout << "Position " << a << "\t";
          for(int m = 0; m < 4; ++m){
            cout << PSSM_matrix[a][m] << "\t";
          }
          cout << "\n";
        }
        cout << "\n";
        */

        int g{0};                                           //Inverse PSSM instead of reverse comp.
        int f = run_motif - 1;
        while( f > -1){
          for(int m = 0; m < 4; ++m){
            in_PSSM_mat[f][m] = PSSM_matrix[g][m];
          }
          ++g;
          --f;
        }

        int pos_start{0};                                   //Score a sequence (k) using the other sequence's motifs. Score is sum of each nucleotide in the motif length (window). Will not calculate score when motif length is longer than remaining nucleotides.
        int finish = run_motif;
        while(finish < max_seq_len - run_motif){
          int temp_nuc;
          int a{0};
          double for_score{0};
          double rev_score{0};
          while(a < run_motif){
            temp_nuc = seq_nuc_data[k][pos_start + a];
            if(temp_nuc == 0){
              for_score += PSSM_matrix[a][0];
              rev_score += in_PSSM_mat[a][1];
            }else if(temp_nuc == 1){
              for_score += PSSM_matrix[a][1];
              rev_score += in_PSSM_mat[a][0];
            }else if(temp_nuc == 2){
              for_score += PSSM_matrix[a][2];
              rev_score += in_PSSM_mat[a][3];
            }else if (temp_nuc == 3){
              for_score += PSSM_matrix[a][3];
              rev_score += in_PSSM_mat[a][2];
            }
            ++a;
          }
          score_to_line[pos_start] = for_score;             //Forward score stored in index 0 to index length of sequences. Reverse score (instead of reverse comp) stored in index max_seq_length to 2*max_seq_length-run_motif.
          score_to_line[pos_start + max_seq_len] = rev_score;
          ++pos_start;
          ++finish;
        }

        double score_sum{0};                                //For number line, exp2 the scores to remove negative numbers.
        for(int h = 0; h < 2 * (max_seq_len - run_motif); ++h){
          score_to_line[h] = exp2(score_to_line[h]);
          score_sum += score_to_line[h];
        }

        for(int h = 0; h < 2 * (max_seq_len - run_motif); ++h){
          score_to_line[h] = score_to_line[h] / score_sum;  //Normalize so sum of all elements is 1.
        }

        double number_line_sum{0};                          //Create number line between 0 and 1. The amount of space an element (start position of window motif (e.g. length) as scored by PSSM) occupies is proportional to its score.
        for(int h = 0; h < 2 * (max_seq_len - run_motif); ++h){
          if(h == 0){
            number_line_sum = score_to_line[h];
          }else{
            score_to_line[h] += number_line_sum;
            number_line_sum = score_to_line[h];
          }
        }
                                                            //Generate random number between 0 and 1.

        for(int h = 0; h < 2 * (max_seq_len - run_motif); ++h){
          if(score_to_line[h] > rand_weight_num){           //If element on number line (must range from 0 to 1) is greater than the random number.  Number line is weighted with more weight given to higher PSSM score.
            if(h < max_seq_len - run_motif){                //If position less seq len - motif length (no partial motif allowed; no overhang past known nucleotides)
              pos_storage[k] = h;                           //Store position
              for(int a = 0; a < run_motif; ++a){           //Store new motif data for next loop
                motif_data_1[k][a] = seq_nuc_data[k][h];
                ++h;
              }
            }else if(h < max_seq_len){                      //End of forward sequence, accounting for overhang in for loop by going backwards.
              pos_storage[k] = h;
              h += run_motif;
              for(int a = run_motif - 1; a > 0; --a){
                motif_data_1[k][a] = seq_nuc_data[k][h];
                --h;
              }
            }else{                                          //If element is greater than seq len - motif len. (Mostly covers reverse comp but also end of forward sequence).
              h = h / 2;
              if(h <  max_seq_len - run_motif){
                pos_storage[k] = h;
              }else{
                pos_storage[k] = h - run_motif;             //Set position as max position where (sequence position + motif length) does not go out of bounds.
              }
              for(int a = run_motif - 1; a > 0; --a){
                motif_data_1[k][a] = seq_nuc_data[k][h];
                --h;
              }
            }
          break;
          }
        }
      }

      for(int j = 0; j < seq_count; ++j){                   //Score the newly selected motifs
        double score{0};
        for(int a = 0; a < run_motif; ++a){
          if(motif_data_1[j][a] == 0){
            score += PSSM_matrix[a][0];
          }else if(motif_data_1[j][a] == 1){
            score += PSSM_matrix[a][1];
          }else if(motif_data_1[j][a] == 2){
            score += PSSM_matrix[a][2];
          }else if(motif_data_1[j][a] == 3){
            score += PSSM_matrix[a][3];
          }
        }
      score_storage_1[j] = score;                           //Store score in storage for current loop.
      }


      if(adj_mot_pos_flag == 1){                            //If motif position was adjusted.
        double sum1{0};
        double sum2{0};
        for(int j = 0; j < seq_count; ++j){
          sum1 += score_storage_1[j];                       //Sum of scores from current run.
          sum2 += score_storage_2[j];                       //Sum of scores from previous run.
        }
        if(sum1 < sum2){                                    //If previous run has higher score,
          for(int j = 0; j < seq_count; ++j){
            for(int a = 0; a < run_motif; ++a){
              motif_data_1[j][a] = motif_data_2[j][a];      //Store the previous run motif data (2) back into the current run data storage (1).
            }
          }
        }

      }

      if(adj_mot_len_flag == 1){                            //If motif length was adjusted.
        double sum1{0};
        double sum2{0};
        double avg1{0};
        double avg2{0};
        for(int j = 0; j < seq_count; ++j){
          sum1 += score_storage_1[j];                       //Sum of scores from current run.
          sum2 += score_storage_2[j];                       //Sum of score from previous run.
        }
        avg1 = sum1 / run_motif;                            //Average to attempt to control for the length of the motif (run v perm) influencing the score. Want nucleotide composition and not length to influence score.
        avg2 = sum2 / perm_motif;
        if(avg1 > avg2){
          perm_motif = run_motif;                           //If current run with different motif length was higher score, set perm_motif length (previous) to current run.
        }

      }

      for(int j = 0; j < seq_count; ++j){
        bool end_flag = 1;                                  //Flag to check if second main loop should end (below in elseif)
        for(int a = 0; a < run_motif; ++a){
          if(motif_data_1[j][a] != motif_data_2[j][a]){     //If current run motif data is not equal to previous run motif data
            motif_data_2[j][a] = motif_data_1[j][a];        //Store current into previous.
            score_storage_2[j] = score_storage_1[j];
            end_flag = 0;                                   //Do not allow second main loop to end. Line below: requires current motif data to match previous run (cycle all sequences and motif position {j & a}) and the second main loop must have run at least 20 times (allow motif length and walking to happen at least 5 times).
          }else if(j == seq_count - 1 && a == run_motif - 1 && end_flag == 1 && loop_counter > 20){
            //cout << "counter " << loop_counter << "\n";   //Uncomment to see how many times second main loop runs.
            loop_flag = 1;
          }
        }
      }
      skipping2:;                                           //Do not check to end the second main loop if the motif position or length was adjusted. Run the loop again.

    }
                                                            //After end of second main loop, make comparisons between seed runs.
    double sum1{0};
    double sum2{0};
    for(int j = 0; j < seq_count; ++j){
      sum1 += score_storage_1[j];                           //Scores of last run of second main loop.
      sum2 += score_storage_3[j];                           //Best scores from best run
    }
    if(sum1 > sum2){                                        //If current run was better than previous seeds
      print_motif = run_motif;                              //print motif used when printing and is equal to the len of motif from the best run
      for(int j = 0; j < seq_count; ++j){                   //Store the best scores and position
        score_storage_3[j] = score_storage_1[j];
        pos_storage2[j] = pos_storage[j];
        for(int a = 0; a < run_motif; ++a){
          motif_data_3[j][a] = motif_data_1[j][a];
        }
      }
    }


  }


}


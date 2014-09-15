#include <iostream>
#include <iomanip>
#include <thread>
#include <vector>
#include <cstring>
#include <fstream>
#include "kernels.h"

static const int work_size = 256;
using namespace std;

template <typename T>
void print_array(const T * const __restrict array, const int len, bool pc)
{
    if (!pc) {
        cout << "[ ";

    }
    for (int i = 0; i < len; i++) {
        if (!pc) {
            cout << setfill(' ') << setw(2) << array[i];
            if (i != len - 1) {
                cout << ',';
            }
            cout << ' ';
        } else {
            cout << i << ' ' << setfill(' ') << setw(2) << array[i] << endl;
        }
    }

    if (!pc) {
        cout << "] ";
    }
    cout << endl;
}

template<typename T>
void read_vector3(const char *path, int *length, T **h0, T **h1, T **h2)
{
    fstream f(path);
    int id;
    f >> *length;
    *h0 = new T[*length];
    *h1 = new T[*length];
    *h2 = new T[*length];
    for (int i = 0; i < *length; i++) {
        f >> id >> (*h0)[i] >> (*h1)[i] >> (*h2)[i];
    }
    f.close();
}

template <typename ... Args>
void p_call(bool serialize, int totalWork, Args&&... argv)
{
    vector<thread> threads;
    int workDone = 0;

    while (totalWork > workDone) {
        for (int i = 0; i < work_size; ++i) {
            threads.push_back(thread(argv..., workDone + i));
            if (serialize) {
                threads.back().join();
            }
        }

        if (!serialize) {
            for (uint i = 0; i < work_size; ++i) {
                threads[i].join();
            }
        }

        workDone += work_size;
//        cout << "work done " << workDone << endl;
        threads.clear();
    }
}

int main()
{
    //int layerOffsets[] = {0, 4, 8, 12, 16};
    //const int nSeedingSets = 5;
    //int hitLayer[] = {0};
    int nHits;
    double *hitsX;
    double *hitsY;
    double *hitsZ;

    int nTriplets;
    int *h0;
    int *h1;
    int *h2;
    read_vector3("hits", &nHits, &hitsX, &hitsY, &hitsZ);
    read_vector3("triplets", &nTriplets, &h0, &h1, &h2);

    cout << "number of hits " << nHits << endl;
    cout << "number of triplets " << nTriplets << endl;
    //print_array(hitsX, nHits, true);
    //print_array(triplets, nTriplets, true);
    int *connectivity_count = new int[nTriplets];
    int *connectivity_oracle = new int[nTriplets];
    p_call(false, nTriplets, vector_init<int>, connectivity_count, 0, nTriplets);
    p_call(false, nTriplets, vector_init<int>, connectivity_oracle, 0, nTriplets);

    //forward

    // connectivity_count holds the number of right neigours of the triplet i
    // connectivity_oracle[i] == true <==> triplet[i] is connectable whith another
    p_call(true, nTriplets, connectivity_count_tight,
           h1, h2, h0, h1,
           //hitLayer, layerOffsets, nSeedingSets,
           nTriplets, connectivity_count, connectivity_oracle);
    //backward
    //p_call(true, nTriplets, connectivity_count_tight, h0, h1, h1, h2,
    //       hitLayer, layerOffsets, nSeedingSets, nTriplets, connectivity_count);

    //print_array(connectivity_count, nTriplets, false);
    //print_array(connectivity_oracle, nTriplets, false);
    int * const connectivity_count_oracle_prefix_sum = new int [nTriplets];
    int total_connectable_triplets = 0;

    prefix_sum(connectivity_oracle,
               connectivity_count_oracle_prefix_sum,
               nTriplets,
               &total_connectable_triplets);

    int * const connectable_triplet_indices = new int[total_connectable_triplets];
    p_call(false, nTriplets,
           predicate_to_valid_index,
           connectivity_oracle,
           connectivity_count_oracle_prefix_sum,
           connectable_triplet_indices,
           nTriplets);

    //print_array(connectable_triplet_indices, total_connectable_triplets, true);
    cout << "connectable triplets " << total_connectable_triplets << endl;
    
    int *connectivity_prefixsum = new int[nTriplets];
    int total_combinations;
    prefix_sum(connectivity_count, connectivity_prefixsum, nTriplets, &total_combinations);
    //print_array(connectivity_prefixsum, nTriplets, false);
    cout << "total combinations " << total_combinations << endl;
    /*
     cout << "stream filter" << endl;
     print_array(connectivity_count, nTriplets, false);
     int *connectable_triplet_indices = new int[total_combinations];
     stream_compaction(connectivity_count, nTriplets, connectable_triplet_indices, total_combinations);
     //print_array(connectable_triplet_indices, total_combinations);
    */

    int *triplet_pairs_base = new int[total_combinations];
    p_call(false, total_combinations,
           where_do_I_belong_function_name_needed,
           connectivity_prefixsum, triplet_pairs_base, nTriplets, total_combinations);

    int *triplet_pairs_follower = new int[total_combinations];
    //p_call(false, total_combinations, vector_init<int>, triplet_pairs_follower, 0, total_combinations);
    p_call(false, total_connectable_triplets,
           make_triplet_pairs_tigh_connectivity,
           connectable_triplet_indices, connectivity_prefixsum, h1, h2, h0, h1,
           total_connectable_triplets, nTriplets, triplet_pairs_follower);
//#define p
#ifdef p
    for (int i = 0; i < total_combinations; i++) {
        std::cout << i << ' ' << triplet_pairs_base[i] << ' ' << triplet_pairs_follower[i] << endl;
    }
#endif
    double * const triplets_eta = new double[nTriplets];
    double * const triplets_pt = new double[nTriplets];

    //init to -1 just for printing purposes
    p_call(false, nTriplets, vector_init<double>, triplets_eta, -1, nTriplets);
    p_call(false, nTriplets, vector_init<double>, triplets_pt, -1, nTriplets);

    p_call(true, total_connectable_triplets,
           circle_fit,
           hitsX, hitsY, hitsZ,
           h0, h1, h2,
           connectable_triplet_indices, total_connectable_triplets,
           triplets_pt, triplets_eta);

    cout << "triplets_pt:" << endl;
    //print_array(triplets_pt, nTriplets, true);

    cout << "triplets_eta:" << endl;
    //print_array(triplets_eta, nTriplets, true);

    // count the final amount of compatible triplets
    int * const compatible_triplet_oracle = new int[total_combinations];
    p_call(true, total_combinations,
           compatible_triplet_filter_eta,
           triplet_pairs_base, triplet_pairs_follower,
           total_combinations,
           triplets_eta,
           (double)0.0256,
           compatible_triplet_oracle);
    //print_array(compatible_triplet_oracle, total_combinations, true);
    int * const compatible_triplet_prefix_sum = new int[total_combinations];
    int compatible_triplet_total = 0;
    prefix_sum(compatible_triplet_oracle,
               compatible_triplet_prefix_sum,
               total_combinations,
               &compatible_triplet_total);

    // compact the stream leaving only those compatible
    int * const compatible_triplets_base = new int [compatible_triplet_total];
    int * const compatible_triplets_follower = new int [compatible_triplet_total];
    p_call(false, total_combinations,
           stream_compaction,
           triplet_pairs_base,
           compatible_triplet_oracle,
           compatible_triplet_prefix_sum,
           total_combinations,
           compatible_triplets_base);

    p_call(false, total_combinations,
           stream_compaction,
           triplet_pairs_follower,
           compatible_triplet_oracle,
           compatible_triplet_prefix_sum,
           total_combinations,
           compatible_triplets_follower);
#define p
#ifdef p
    for (int i = 0; i < compatible_triplet_total; i++) {
        std::cout << i << ' ' << compatible_triplets_base[i] << ' ' << compatible_triplets_follower[i] << endl;
    }
#endif

#ifdef not_implemented


    // start the CA
    p_call(false, total_connectable_triplets,
           ca - automaton
           connectable_triplet_indices, total_connectable_triplets
#endif
           return 0;
}












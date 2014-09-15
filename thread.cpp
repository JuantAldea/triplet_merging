
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
    cout << "[ ";
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

    cout << "] ";
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
    int layerOffsets[] = {0, 4, 8, 12, 16};
    const int nSeedingSets = 5;
//#define mock
#ifdef mock
    const int nTriplets = 20;
    int h0[] = {0, 1,  2,  3,  4,  4,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    int h1[] = {4, 4,  6,  7,  8,  8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};
    int h2[] = {8, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27};
    int hitLayer[] = {
        0, 0, 0, 0,
        1, 1, 1, 1,
        2, 2, 2, 2,
        3, 3, 3, 3,
        4, 4, 4, 4,
        5, 5, 5, 5,
        8, 8, 8, 8
    };

    for (int i = 0; i < nTriplets; i++) {
        cout << setfill(' ') << setw(2) << i << ' ';
    }
    cout << endl;
    cout << "-----" << endl;
    print_array(h0, nTriplets, false);
    print_array(h1, nTriplets, false);
    print_array(h2, nTriplets, false);
    cout << "-----" << endl;
#else
    int nHits;
    float *hitsX;
    float *hitsY;
    float *hitsZ;

    int nTriplets;
    int *h0;
    int *h1;
    int *h2;
    read_vector3("hits", &nHits, &hitsX, &hitsY, &hitsZ);
    read_vector3("log_7977", &nTriplets, &h0, &h1, &h2);

    //TODO think of filtering the list of triplets leaving only those
    // with connectivity
    //double * const tripletsEta = new double[nTriplets];
    //double * const tripletsPt = new double[nTriplets];

    int hitLayer[] = {0};
#endif
    cout << "number of hits " << nHits << endl;
    cout << "number of triplets " << nTriplets << endl;

    int *connectivity_count = new int[nTriplets];
    int *connectivity_oracle = new int[nTriplets];
    p_call(false, nTriplets, vector_init, connectivity_count, 0, nTriplets);
    p_call(false, nTriplets, vector_init, connectivity_oracle, 0, nTriplets);
    //forward

    p_call(true, nTriplets, connectivity_count_tight,
           h1, h2, h0, h1,
           hitLayer, layerOffsets, nSeedingSets, nTriplets, connectivity_count, connectivity_oracle);
    //backward
    //p_call(true, nTriplets, connectivity_count_tight, h0, h1, h1, h2,
    //       hitLayer, layerOffsets, nSeedingSets, nTriplets, connectivity_count);

    //print_array(connectivity_count, nTriplets, false);
    //print_array(connectivity_oracle, nTriplets, false);
    int *connectable_triplet_indices = new int[nTriplets];
    int nConnectables = 0;
    stream_compaction(connectivity_oracle, nTriplets, connectable_triplet_indices, nConnectables);
    //print_array(connectable_triplet_indices, nConnectables, true);
    cout << "connectable triplets " << nConnectables << endl;
    int *connectivity_prefixsum = new int[nTriplets];
    int total_combinations;
    prefix_sum(connectivity_count, nTriplets, connectivity_prefixsum, &total_combinations);
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
    //p_call(false, total_combinations, vector_init, triplet_pairs_follower, 0, total_combinations);
    p_call(false, nConnectables,
           make_triplet_pairs_tigh_connectivity,
           connectable_triplet_indices, connectivity_prefixsum, h1, h2, h0, h1,
           nConnectables, nTriplets, triplet_pairs_follower);
#define p
#ifdef p
    for (int i = 0; i < total_combinations; i++) {
        std::cout << i << ' ' << triplet_pairs_base[i] << ' ' << triplet_pairs_follower[i] << endl;
    }
#endif

#ifdef not_implemented
    p_call(false, nConnectables,
           circle_fit,
           connectable_triplet_indices, nConnectables, h0, h1, h2, triplet_pt, triplet_eta);
    // count the final amount of compatible triplets
    p_call(false, total_combinations,
           compatible_triplet_test_count,
           triplet_pair_base, triplet_pair_follower,
           triplet_pt, triplet_eta,
           total_combinations,
           compatible_triplet_oracle);

    prefix_sum(compatible_triplet_oracle, compatible_triplet_prefix_sum, compatible_triplet_total);

    // compact the stream leaving only those compatible
    int *compatible_triplet_base = new int[compatible_triplet_total];
    int *compatible_triplet_followers = new int[compatible_triplet_total];
    stream_compaction(false, total_combinations,
                      compatible_triplet_oracle,
                      compatible_triplet_prefix_sum,
                      triplet_pairs_base, triplet_pairs_follower,
                      compatible_triplet_base, compatible_triplet_followers);


    // start the CA
    p_call(false, nConnectables,
           ca - automaton
           connectable_triplet_indices, nConnectables
#endif
           /*
               int kf_result[total_combinations] = {0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1};
               int kf_result_prefix_sum[total_combinations];
               int total_combinations_kf;
               prefix_sum(kf_result, total_combinations, kf_result_prefix_sum, &total_combinations_kf);
               print_array(connectivity_prefixsum, total_combinations_kf, false);
               cout << "total_combinations_kf " << total_combinations_kf << endl;

               int *after_kf_filter = new int [total_combinations_kf];

               stream_compaction(kf_result, total_combinations, after_kf_filter, total_combinations_kf);
               print_array(after_kf_filter, total_combinations_kf, false);
           */
           return 0;
}









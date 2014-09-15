#pragma once

#include <iostream>
#include <mutex>
using namespace std;

void vector_init(int * const __restrict vector,
                 const int value,
                 const int len,
                 const int gid)
{
    if (gid >= len) {
        return;
    }

    vector[gid] = value;
}

void where_do_I_belong_function_name_needed(
    const int * const __restrict prefixsum,
    int * const __restrict expanded_name_needed,
    const int len, const int total,
    const int gid)
{
    //gid < total....
    if (gid >= total) {
        return;
    }

    expanded_name_needed[gid] = 0;
    int where_do_I_belong = 0;
    // NEW!
    // 100% Branch divergence free on GPU!
    // Advertised on TV!
    for (int i = 0; i < len; i++) {
        const int test1 = (prefixsum[i] <= gid);
        int test2 = 1;
        if (i < len - 1) {
            test2 = (prefixsum[i + 1] > gid);
        } else {
            test2 = total > gid;
        }
        where_do_I_belong += i * test1 * test2 ;
    }

    expanded_name_needed[gid] = where_do_I_belong;
}

/*
void connectivity_count_wide(const int gid,
                             const int * const __restrict hitsBasis,
                             const int * const __restrict hitsFollowers,
                             const int * const __restrict hitLayer,
                             const int * const __restrict layerOffsets,
                             const int nSeedingLayers, const int nTriplets,
                             int * __restrict connectivityCount)
{
    if (gid >= nTriplets) {
        return;
    }

    const int tripletIndex = gid;
    const int hitIndex = hitsBasis[tripletIndex];
    const int thisHitLayer = hitLayer[hitIndex];

    int compatibleLayerOffsetIndex = 0;
    for (int i = 1; i < nSeedingLayers; i++) {
        compatibleLayerOffsetIndex += i * (thisHitLayer == hitLayer[layerOffsets[i]]);// * (myLayerIndex < i);
    }
    const int searchOffsetBegin = layerOffsets[compatibleLayerOffsetIndex];
    const bool compatibleLayerIsTheLastLayer = (compatibleLayerOffsetIndex + 1) == nSeedingLayers;
    const int searchOffsetEnd =
        !compatibleLayerIsTheLastLayer *
        layerOffsets[compatibleLayerOffsetIndex + !compatibleLayerIsTheLastLayer]
        + compatibleLayerIsTheLastLayer * nTriplets;

    int connectivityCountLocal = 0;
    for (int i = searchOffsetBegin; i < searchOffsetEnd; i++) {
        connectivityCountLocal += (hitIndex == hitsFollowers[i]);
    }

    connectivityCount[gid] = connectivityCountLocal;
}
*/

void connectivity_count_tight(const int * const __restrict hitsBasisH1,
                              const int * const __restrict hitsBasisH2,
                              const int * const __restrict hitsFollowersH0,
                              const int * const __restrict hitsFollowersH1,
                              const int * const __restrict hitLayer,
                              const int * const __restrict layerOffsets,
                              const int nSeedingLayers,
                              const int nTriplets,
                              int * const __restrict connectivityCount,
                              int * const __restrict connectivityOracle,
                              const int gid)
{
    if (gid >= nTriplets) {
        return;
    }

    const int tripletIndex = gid;
    const int hitIndexInner = hitsBasisH1[tripletIndex];
    const int hitIndexOuter = hitsBasisH2[tripletIndex];
    const int tripletInnerHitLayer = hitLayer[hitIndexInner];

    /*
    int myLayerIndex = 0;
    for (int i = 0; i < nSeedingLayers - 1; i++) {
        myLayerIndex += i * (layerOffsets[i] <= hitIndex) * (hitIndex < layerOffsets[i + 1]);
    }

    myLayerIndex += (nSeedingLayers - 1) * (layerOffsets[nSeedingLayers- 1] <= hitIndex);
    */

    int compatibleLayerOffsetIndex = 0;
    for (int i = 1; i < nSeedingLayers; i++) {
        compatibleLayerOffsetIndex +=
            i * (tripletInnerHitLayer == hitLayer[layerOffsets[i]]);// * (myLayerIndex < i);
    }

    /*
    const int searchOffsetBegin = layerOffsets[compatibleLayerOffsetIndex];
    const bool compatibleLayerIsLastLayer = (compatibleLayerOffsetIndex + 1) == nSeedingLayers;
    const int searchOffsetEnd =
        !compatibleLayerIsLastLayer * layerOffsets[compatibleLayerOffsetIndex + !compatibleLayerIsLastLayer]
        + compatibleLayerIsLastLayer * nTriplets;
    */
    int connectivityCountLocal = 0;
    //connectivityCount[gid] = 0;
    //for (int i = searchOffsetBegin; i < searchOffsetEnd; i++) {
    for (int i = 0; i < nTriplets; i++) {
        const bool test = (hitIndexInner == hitsFollowersH0[i]) * (hitIndexOuter == hitsFollowersH1[i]);
        if (gid == nTriplets - 1) {

            //    cout << i <<  " ultimo " << test <<  ' ' << connectivityCountLocal << endl;
        }
        connectivityCountLocal += test;
        connectivityOracle[i] |= test;
    }

    connectivityOracle[gid] |= connectivityCountLocal > 0;
    connectivityCount[gid] += connectivityCountLocal;
}

void make_triplet_pairs_tigh_connectivity(
    const int * const __restrict connectable_triplets,
    const int * const __restrict connectiviy_prefix_sum,
    const int * const __restrict hitsBasisH1,
    const int * const __restrict hitsBasisH2,
    const int * const __restrict hitsFollowersH0,
    const int * const __restrict hitsFollowersH1,
    const int nConnectableTriplets,
    const int nTriplets,
    int * const __restrict follower_triplet_indices,
    const int gid)
{
    if (gid >= nConnectableTriplets) {
        return;
    }

    const int tripletIndex = connectable_triplets[gid];
    const int hitIndexInner = hitsBasisH1[tripletIndex];
    const int hitIndexOuter = hitsBasisH2[tripletIndex];
    const int offset = connectiviy_prefix_sum[tripletIndex];
    /* cout << gid
          << " index " << tripletIndex
          << " inner " << hitIndexInner
          << " outer " << hitIndexOuter
          << " offset " << offset
          << endl;
          */
    //todo use layer information
    int connectivityCountLocal = 0;
    //TODO add layer information
    int searchOffsetBegin = tripletIndex + 1;
    int searchOffsetEnd = nTriplets;
    for (int i = searchOffsetBegin; i < searchOffsetEnd; i++) {
        const bool test = (hitIndexInner == hitsFollowersH0[i]) * (hitIndexOuter == hitsFollowersH1[i]);
        connectivityCountLocal += test;
        if (test) {
            follower_triplet_indices[offset + connectivityCountLocal - 1] = i;
        }
    }
}


void prefix_sum(const int * const __restrict array,
                const int array_len,
                int * const __restrict output,
                int * const total)
{
    output[0] = 0;
    for (int i = 0 ; i < array_len - 1; i++) {
        output[i + 1] = output[i] + array[i];
    }
    *total = output[array_len - 1] + array[array_len - 1];
}

void stream_compaction(const int * const __restrict predicate,
                       const int  input_length,
                       int * __restrict output,
                       int &output_size)
{
    int output_index = 0;
    for (int i = 0; i < input_length; i++) {
        if (predicate[i] > 0) {
            output[output_index] = i;
            output_index++;
        }
    }
    output_size = output_index;
}


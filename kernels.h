#pragma once

#include <iostream>
#include <cmath>
using namespace std;
template<typename T>
void vector_init(T * const __restrict vector,
                 const T value,
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
 //                             const int * const __restrict hitLayer,
 //                             const int * const __restrict layerOffsets,
 //                             const int nSeedingLayers,
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
//    const int tripletInnerHitLayer = hitLayer[hitIndexInner];

    /*
    int myLayerIndex = 0;
    for (int i = 0; i < nSeedingLayers - 1; i++) {
        myLayerIndex += i * (layerOffsets[i] <= hitIndex) * (hitIndex < layerOffsets[i + 1]);
    }

    myLayerIndex += (nSeedingLayers - 1) * (layerOffsets[nSeedingLayers- 1] <= hitIndex);
    */
/*
    int compatibleLayerOffsetIndex = 0;
    for (int i = 1; i < nSeedingLayers; i++) {
        compatibleLayerOffsetIndex +=
            i * (tripletInnerHitLayer == hitLayer[layerOffsets[i]]);// * (myLayerIndex < i);
    }
*/
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
                int * const __restrict total)
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

template<typename T>
struct vector3 {
    vector3(const T ox, const T oy, const T oz)
    {
        this->x = ox;
        this->y = oy;
        this->z = oz;
    }
    T x;
    T y;
    T z;
};

template<typename T>
struct vector2 {
    vector2(const T ox, const T oy)
    {
        this->x = ox;
        this->y = oy;
    }

    T x;
    T y;
};

void circle_fit(const double * const __restrict hitX,
                const double * const __restrict hitY,
                const double * const __restrict hitZ,
                const int * const __restrict h0,
                const int * const __restrict h1,
                const int * const __restrict h2,
                const int * const __restrict connectable_triplets_indices,
                const int n_connectable_triplets,
                double * const __restrict triplet_pt,
                double * const __restrict triplet_eta,
                const int gid)
{
    if (gid >= n_connectable_triplets) {
        return;
    }

    const int triplet_index = connectable_triplets_indices[gid];
    const int innerHit = h0[triplet_index];
    const int middleHit = h1[triplet_index];
    const int outerHit = h2[triplet_index];
    //some room for opencl vector operations: dot, substraction...
    const vector3<double> hit0(hitX[innerHit], hitY[innerHit], hitZ[innerHit]);
    const vector3<double> hit1(hitX[middleHit], hitY[middleHit], hitZ[middleHit]);
    const vector3<double> hit2(hitX[outerHit], hitY[outerHit], hitZ[outerHit]);

    //calculate Pt
    const vector3<double> pP0(hit0.x, hit0.y, hit0.x * hit0.x + hit0.y * hit0.y);
    const vector3<double> pP1(hit1.x, hit1.y, hit1.x * hit1.x + hit1.y * hit1.y);
    const vector3<double> pP2(hit2.x, hit2.y, hit2.x * hit2.x + hit2.y * hit2.y);

    const vector3<double> a(pP1.x - pP0.x, pP1.y - pP0.y, pP1.z - pP0.z);
    const vector3<double> b(pP2.x - pP0.x, pP2.y - pP0.y, pP2.z - pP0.z);

    const vector3<double> n(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    const double n_module = sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
    const vector3<double> unit_n(n.x / n_module, n.y / n_module, n.z / n_module);
    const vector2<double> circle_center(-unit_n.x / (2 * unit_n.z), -unit_n.y / (2 * unit_n.z));
    const double c = -(unit_n.x * pP0.x + unit_n.y * pP0.y + unit_n.z * pP0.z);
    const double circle_radius = sqrt((1 - unit_n.z * unit_n.z - 4 * c * unit_n.z) / (4 * unit_n.z * unit_n.z));

    // Ideal Magnetic Field [T] (-0,-0, 3.8112)
    const double BZ = 3.8112;
    // e = 1.602177×10^-19 C (coulombs)
    const double Q = 1.602177E-19;
    // c = 2.998×10^8 m/s (meters per second)
    //const double C = 2.998E8;
    // 1 GeV/c = 5.344286×10^-19 J s/m (joule seconds per meter)
    const double GEV_C = 5.344286E-19;

    triplet_pt[triplet_index] = Q * BZ * (circle_radius * 1E-2) / GEV_C;

    // calculate eta
    const vector3<double> p(hit2.x - hit0.x, hit2.y - hit0.y, hit2.z - hit0.z);
    const double t = p.z / sqrt(p.x * p.x + p.y * p.y);
    triplet_eta[triplet_index] = asinh(t);
}


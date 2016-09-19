#pragma once
#ifndef HASH_EXTEND_H
#define HASH_EXTEND_H

#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <algorithm>

namespace std {
    template<typename T> struct hash<vector<T>> {
        inline size_t operator() (const vector<T>& vec) const {
            hash<T> hasher;
            size_t seed = 0;
            for (auto& i : vec) {
                seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    template<typename T> struct hash<unordered_multiset<T>> {
        inline size_t operator() (const unordered_multiset<T>& uset) const {
            hash<T> hasher;
            vector<size_t> hashes;
            for (auto& i : uset) {
                hashes.push_back(hasher(i));
            }
            sort(hashes.begin(), hashes.end());
            hash<vector<size_t>> h;

            return h(hashes);
        }
    };
}

#endif
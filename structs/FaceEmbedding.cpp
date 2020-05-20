//
// Created by pavel on 22.09.2019.
//

#include "FaceEmbedding.h"


embedding::FaceEmbedding::FaceEmbedding(unsigned long id, std::string name,
                                        std::vector<double> &arr) : id(id), name(std::move(name)) {
    dim = arr.size();
    embedding = new double[arr.size()];
    for (size_t i = 0; i < arr.size(); i++) {
        embedding[i] = arr[i];
    }
}



embedding::FaceEmbedding::FaceEmbedding(unsigned long id, std::string name,
                                        const double arr[], size_t n) : id(id), name(std::move(name)) {
    dim = n;
    embedding = new double[n];
    for (size_t i = 0; i < n; i++) {
        embedding[i] = arr[i];
    }
}

embedding::FaceEmbedding::~FaceEmbedding() {
    delete [] embedding;
}

embedding::FaceEmbedding::FaceEmbedding(const embedding::FaceEmbedding &object) {
    dim = object.length();
    embedding = new double[object.length()];
    for (size_t i = 0; i < object.length(); i++) {
        embedding[i] = object[i];
    }
}

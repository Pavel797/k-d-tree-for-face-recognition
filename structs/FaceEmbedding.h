//
// Created by pavel on 22.09.2019.
//

#ifndef FACESDATABASE_FACEEMBEDDING_H
#define FACESDATABASE_FACEEMBEDDING_H

#include <array>
#include <cassert>
#include <vector>

namespace embedding {
    class FaceEmbedding  {
        std::string name;
        unsigned long id;
        double *embedding;
        size_t dim;
    public:
        FaceEmbedding(unsigned long id, std::string name, std::vector<double> &arr);

        FaceEmbedding(unsigned long id, std::string name, const double arr[], size_t n);

        FaceEmbedding(const FaceEmbedding & object);

        ~FaceEmbedding();

        double operator[](size_t i) const {
            return embedding[i];
        }

        [[nodiscard]] size_t length() const {
            return dim;
        }
    };
}
#endif //FACESDATABASE_FACEEMBEDDING_H

#ifndef METRIC_H
#define METRIC_H

#include <vector>

template<class T>
class Metric {
    public:
        virtual T factor(long long i, long long j) = 0;
        T factorByMask (auto i, auto j) {
            return 1;
        }

};

template<class T>
class NonOrthogonalMetric : public Metric<T> {

    private:
        std::vector<std::vector<T>> Q;

    public:
        NonOrthogonalMetric() {}
        NonOrthogonalMetric(std::vector<std::vector<T>> Q){
            this->Q = Q;
        }
        T factor(long long i, long long j){
            return Q[i][j];
        }

        T factorByMask (auto i, auto j) {
            return 1;
        }

};


template<class T>
class OrthogonalMetric : public Metric<T> {

    private:
        std::vector<T> d;

    public:
        OrthogonalMetric() {}
        OrthogonalMetric(std::vector<T> d) {
            this->d = d;
        }

        T factor (long long i, long long j) {
            if (i == j) {
                return d[i];
            }
            return 0;
        }

        T factorByMask (auto i, auto j) {
            return i == j;
        }


};

template<class T>
class OrthonormalMetric : public OrthogonalMetric<T> {

    public :
        OrthonormalMetric() {}

        T factor (long long i, long long j) {
            return i == j;
        }

        T factorByMask (auto i, auto j) {
            return i == j;
        }
};


#endif // METRIC_H


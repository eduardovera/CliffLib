#ifndef METRIC_H
#define METRIC_H

#include <vector>

template<class T>
class Metric {
    public:
        Metric(){}
        virtual T factor(int i, int j){}
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

        T factor (int i, int j) {
            if (i == j) {
                return d[i];
            }
            return 0;
        }
};

template<class T>
class OrthonormalMetric : public OrthogonalMetric<T> {

    public :
        OrthonormalMetric() {}

        T factor (int i, int j) {
            return i == j;
        }
};


#endif // METRIC_H


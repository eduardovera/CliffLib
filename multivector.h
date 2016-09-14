#ifndef MULTIVECTOR_H
#define MULTIVECTOR_H

#include <iostream>
#include <map>
#include <bitset>
#include <algorithm>

using namespace std;

template<class coeff_type>
class multivector {

    private:
        std::map<uint, coeff_type> M;

        string get_base(int index) const {
            string ret = "";
            std::string binary = std::bitset<32>(index).to_string();
            reverse(binary.begin(), binary.end());
            string s = "1";

            size_t pos = binary.find(s);
            if (pos == string::npos) {
            }
            while(pos != string::npos) {
                ret.append("e");
                ret.append(to_string(pos + 1));
                pos = binary.find(s, pos+1);
                if (pos != string::npos) {
                    ret.append("^");
                }
            }
            return ret;
        }

    public:

        multivector() {
        }

        multivector(const std::initializer_list<coeff_type> &l) {
            int index = 0;
            for (auto it = l.begin(); it != l.end(); ++it) {
                if ((*it) != 0.0) {
                    M[index] = (*it);
                }
                index++;
            }
        }

    template<class T> friend std::ostream& operator << (std::ostream &, const multivector<T> &);
    template<class T> friend multivector<T> e(int);
};


template<class coeff_type>
std::ostream& operator << (std::ostream &os, const multivector<coeff_type> &m) {
    for (auto it = m.M.begin(); it != m.M.end(); ++it) {
        os << ((*it).second > 0.0 ? "+" : "") << (*it).second << "(" << m.get_base((*it).first) << ")";
    }
    return os;
}

template<class coeff_type = double>
multivector<coeff_type> scalar(coeff_type value) {
    multivector<coeff_type> r;
    r.M[0] = value;
    return r;
}


template<class coeff_type = double>
multivector<coeff_type> e(int index) {
    multivector<coeff_type> r;
    r.M[1 << index] = 1;
    return r;
}

//template<class coeff_type1, class coeff_type2>
//multivector<std::common_type<coeff_type1, coeff_type2>::type> operator + (const multivector<coeff_type1> &m1, const multivector<coeff_type1> &m2) {

//}

//template<class coeff_type1, class coeff_type2>
//multivector<std::common_type<coeff_type1, coeff_type2>::type> operator + (const coeff_type1 m1, const multivector<coeff_type1> &m2) {
//    return scalar(m1) + m2;
//}


//template<class coeff_type1, class coeff_type2>
//multivector<std::common_type<coeff_type1, coeff_type2>::type> operator + (const multivector<coeff_type1> &m1, const coeff_type1 m2) {
//    return m1 + scalar(m2);
//}

#endif // MULTIVECTOR_H
